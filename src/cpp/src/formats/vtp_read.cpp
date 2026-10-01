//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "../detail/region_field_data.hpp"
#include "../detail/vtk_xml_read.hpp"

namespace meshioplusplus {

namespace {

using detail::vtu_to_int64;

// The codec is resolved once from the root's compressor= attribute.
NDArray vtp_read_data_array(const pugi::xml_node& rDa, const detail::VtuContext& rCtx,
                            int& rNumComponents) {
    return detail::vtu_read_data_array(rDa, rCtx, rNumComponents);
}

// One PolyData section's connectivity + VTK end-offsets.
struct VtpPiece {
    std::vector<std::int64_t> mConn;
    std::vector<std::int64_t> mOffsets;
    bool mPresent = false;
};

/**
 * @param offsets_only Skip the `connectivity` array. PolyData has no `types`
 *        array -- cell types are synthesized from each section's per-cell size
 *        -- so `offsets` alone (one value per cell) is enough to summarize a
 *        section, while `connectivity` is the bulk of the section's bytes.
 */
VtpPiece vtp_read_section(const pugi::xml_node& rSection, const detail::VtuContext& rCtx,
                          bool offsets_only = false) {
    VtpPiece out;
    if (!rSection)
        return out;
    out.mPresent = true;
    for (pugi::xml_node da : rSection.children("DataArray")) {
        std::string name = da.attribute("Name").as_string();
        if (offsets_only && name != "offsets")
            continue;
        int nc = 0;
        NDArray arr = vtp_read_data_array(da, rCtx, nc);
        if (name == "connectivity")
            out.mConn = vtu_to_int64(arr);
        else if (name == "offsets")
            out.mOffsets = vtu_to_int64(arr);
    }
    return out;
}

/** @brief `<Piece>` plus the framing attributes; mirrors `vtu_parse_header`. */
struct vtp_header {
    pugi::xml_node mGrid;
    pugi::xml_node mPiece;
    detail::VtuContext mCtx;
    std::size_t mNumPoints = 0;
};

vtp_header vtp_parse_header(const detail::VtuSource& rSource) {
    const auto& rDoc = rSource.mDoc;
    pugi::xml_node root = rDoc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile'");
    if (std::string(root.attribute("type").as_string()) != "PolyData")
        throw ReadError("Expected type PolyData");

    vtp_header h;
    pugi::xml_node grid = root.child("PolyData");
    if (!grid)
        throw ReadError("No PolyData found");

    h.mCtx = detail::vtk_xml_read_context(rSource, "VTP");

    h.mGrid = grid;
    h.mPiece = grid.child("Piece");
    if (!h.mPiece)
        throw ReadError("No Piece found");
    // A single piece is supported; multiple pieces -> Python reader.
    if (h.mPiece.next_sibling("Piece"))
        throw ReadError("multi-piece VTP not supported by the C++ reader");

    h.mNumPoints = static_cast<std::size_t>(h.mPiece.attribute("NumberOfPoints").as_ullong());
    return h;
}

/** @brief `<DataArray>` `Name` attributes under @p pSection, sorted. */
std::vector<std::string> vtp_array_names(const pugi::xml_node& rPiece, const char* pSection) {
    std::vector<std::string> names;
    for (pugi::xml_node da : rPiece.child(pSection).children("DataArray"))
        names.emplace_back(da.attribute("Name").as_string());
    std::sort(names.begin(), names.end());
    return names;
}

/** @brief Whether a `<DataArray type=>` is one of the ten numeric types meshio++ holds. */
bool vtp_is_numeric_type(const std::string& rType) {
    return rType == "Float32" || rType == "Float64" || rType == "Int8" || rType == "Int16" ||
           rType == "Int32" || rType == "Int64" || rType == "UInt8" || rType == "UInt16" ||
           rType == "UInt32" || rType == "UInt64";
}

/**
 * @brief Read the `<FieldData>` arrays under @p rNode into `mesh.field_data`.
 *
 * Field data belongs to the dataset, not to a piece: VTK writes it on the
 * `<PolyData>` element, before the `<Piece>`, and also accepts it inside one, so the
 * reader looks at both (the piece's overriding the grid's, since `AddFieldData`
 * is insert-or-assign). A non-numeric array (`type="String"`, `"Bit"`) has no
 * meshio++ dtype: it is skipped with a warning rather than failing a read that
 * used to succeed by ignoring the whole section.
 */
void vtp_read_field_data(const pugi::xml_node& rNode, const detail::VtuContext& rCtx,
                         const ReadOptions& rOpts, bool WantData,
                         std::vector<std::pair<std::string, NDArray>>& rOut) {
    for (pugi::xml_node da : rNode.child("FieldData").children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        // Region arrays are topology, not data (detail/region_field_data.hpp).
        if (!detail::is_region_field_name(name) && (!WantData || !rOpts.WantsArray(name)))
            continue;
        if (!vtp_is_numeric_type(da.attribute("type").as_string())) {
            log::warn(
                "meshio++: VTP: skipping <FieldData> array '{}' of type '{}' (only numeric "
                "arrays are read)",
                name, da.attribute("type").as_string());
            continue;
        }
        int nc = 0;
        NDArray arr = vtp_read_data_array(da, rCtx, nc);
        if (nc > 1)
            arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
        rOut.emplace_back(name, std::move(arr));
    }
}

/**
 * @brief The field-data names a real read would return: the numeric arrays of the
 * grid's and the piece's `<FieldData>`, sorted and unique.
 *
 * Numeric only, like `vtp_read_field_data`, so a summary never names an array the
 * read skips.
 */
std::vector<std::string> vtp_field_data_names(const vtp_header& rHeader) {
    std::vector<std::string> names;
    for (const pugi::xml_node& rNode : {rHeader.mGrid, rHeader.mPiece})
        for (pugi::xml_node da : rNode.child("FieldData").children("DataArray"))
            if (vtp_is_numeric_type(da.attribute("type").as_string()) &&
                !detail::is_region_field_name(da.attribute("Name").as_string()))
                names.emplace_back(da.attribute("Name").as_string());
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
}

/**
 * @brief Synthesize `types`/`offsets` for the three PolyData sections.
 *
 * Factored out of `read_vtp` so the mesh and metadata paths derive cell types
 * from section sizes in exactly one place, and therefore cannot disagree.
 */
void vtp_build_types(const VtpPiece& rSec, int kind, std::vector<std::int64_t>& rConn,
                     std::vector<std::int64_t>& rOffsets, std::vector<std::int64_t>& rTypes) {
    const std::int64_t conn_base = static_cast<std::int64_t>(rConn.size());
    std::int64_t prev = 0;
    for (std::int64_t end : rSec.mOffsets) {
        const std::int64_t sz = end - prev;
        prev = end;
        std::int64_t vtk_type = 0;
        if (kind == 0) {
            if (sz != 1)
                throw ReadError("poly-vertex VTP cells not supported by the C++ reader");
            vtk_type = 1;  // VTK_VERTEX
        } else if (kind == 1) {
            if (sz != 2)
                throw ReadError("poly-line VTP cells not supported by the C++ reader");
            vtk_type = 3;  // VTK_LINE
        } else {
            vtk_type = sz == 3 ? 5 : sz == 4 ? 9 : 7;  // triangle / quad / polygon
        }
        rTypes.push_back(vtk_type);
        rOffsets.push_back(conn_base + end);
    }
    rConn.insert(rConn.end(), rSec.mConn.begin(), rSec.mConn.end());
}

}  // namespace

Mesh read_vtp(const std::string& rPath, const ReadOptions& rOpts) {
    detail::VtuSource source;
    detail::vtu_load(rPath, pugi::parse_default, source, "PolyData", "VTP");
    const vtp_header h = vtp_parse_header(source);
    const pugi::xml_node piece = h.mPiece;
    const detail::VtuContext& ctx = h.mCtx;
    const std::size_t num_points = h.mNumPoints;
    const bool want_data = rOpts.WantsAnyData();

    Mesh mesh;
    std::unordered_map<std::string, NDArray> cell_data_raw;

    for (pugi::xml_node child : piece.children()) {
        std::string tag = child.name();
        if (tag == "Points") {
            pugi::xml_node da = child.child("DataArray");
            int nc = 0;
            NDArray pts = vtp_read_data_array(da, ctx, nc);
            if (nc <= 0)
                nc = 3;
            if (pts.Size() / static_cast<std::size_t>(nc) != num_points ||
                pts.Size() % static_cast<std::size_t>(nc) != 0)
                throw ReadError("VTP Points length differs from NumberOfPoints");
            pts.Reshape({num_points, static_cast<std::size_t>(nc)});
            mesh.AssignPoints(std::move(pts));
        } else if (tag == "PointData") {
            if (!want_data)
                continue;
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                // Name is readable before the payload -- skipping is free.
                if (!rOpts.WantsArray(name))
                    continue;
                NDArray arr = vtp_read_data_array(da, ctx, nc);
                if (arr.Size() / static_cast<std::size_t>(nc > 0 ? nc : 1) != num_points)
                    throw ReadError("VTP PointData length differs from NumberOfPoints");
                if (nc > 1)
                    arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                mesh.AddPointData(name, std::move(arr));
            }
        } else if (tag == "CellData") {
            if (!want_data)
                continue;
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                if (!rOpts.WantsArray(name))
                    continue;
                NDArray arr = vtp_read_data_array(da, ctx, nc);
                if (nc > 1)
                    arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                cell_data_raw.emplace(name, std::move(arr));
            }
        }
    }

    std::vector<std::pair<std::string, NDArray>> field_arrays;
    vtp_read_field_data(h.mGrid, ctx, rOpts, want_data, field_arrays);
    vtp_read_field_data(piece, ctx, rOpts, want_data, field_arrays);

    VtpPiece verts = vtp_read_section(piece.child("Verts"), ctx);
    VtpPiece lines = vtp_read_section(piece.child("Lines"), ctx);
    VtpPiece polys = vtp_read_section(piece.child("Polys"), ctx);
    VtpPiece strips = vtp_read_section(piece.child("Strips"), ctx);
    if (!strips.mOffsets.empty())
        throw ReadError("triangle-strip VTP cells not supported by the C++ reader");

    // Concatenate sections in VTK's canonical PolyData cell order (Verts,
    // Lines, Polys), synthesizing a VTK type id per row so the shared
    // reconstruction (detail/vtk_cells.hpp) can build the blocks and split
    // cell_data.
    std::vector<std::int64_t> conn, offsets, types;
    vtp_build_types(verts, 0, conn, offsets, types);
    vtp_build_types(lines, 1, conn, offsets, types);
    vtp_build_types(polys, 2, conn, offsets, types);

    detail::check_vtk_cell_arrays(conn.size(), offsets, types, cell_data_raw);
    static const std::vector<std::int64_t> kNoFaceOffsets;
    std::vector<std::int64_t> file_to_global;
    detail::reconstruct_cells(conn.data(), offsets, types, cell_data_raw, nullptr, kNoFaceOffsets,
                              mesh, &file_to_global);
    detail::regions_from_field_arrays(mesh, field_arrays, &file_to_global, "vtp");
    return mesh;
}

MeshMetadata read_vtp_metadata(const std::string& rPath, const ReadOptions&) {
    detail::VtuSource source;
    // See read_vtu_metadata: parse_minimal trims text conversions, but the
    // saving that matters is skipping the array bodies below.
    detail::vtu_load(rPath, pugi::parse_minimal, source, "PolyData", "VTP");
    const vtp_header h = vtp_parse_header(source);

    MeshMetadata meta;
    meta.mNumPoints = h.mNumPoints;  // an attribute -- free

    pugi::xml_node points_da = h.mPiece.child("Points").child("DataArray");
    const int point_nc = points_da ? points_da.attribute("NumberOfComponents").as_int(0) : 0;
    meta.mPointDim = point_nc > 0 ? static_cast<std::size_t>(point_nc) : 3;

    // PolyData carries no `types` array; cell types follow from each section's
    // per-cell size, so reading `offsets` alone suffices and the connectivity --
    // the bulk of the bytes -- is never decoded.
    const VtpPiece verts = vtp_read_section(h.mPiece.child("Verts"), h.mCtx,
                                            /*offsets_only=*/true);
    const VtpPiece lines = vtp_read_section(h.mPiece.child("Lines"), h.mCtx,
                                            /*offsets_only=*/true);
    const VtpPiece polys = vtp_read_section(h.mPiece.child("Polys"), h.mCtx,
                                            /*offsets_only=*/true);
    if (h.mPiece.child("Strips") && !vtp_read_section(h.mPiece.child("Strips"), h.mCtx,
                                                      /*offsets_only=*/true)
                                         .mOffsets.empty())
        throw ReadError("triangle-strip VTP cells not supported by the C++ reader");

    std::vector<std::int64_t> conn, offsets, types;
    vtp_build_types(verts, 0, conn, offsets, types);
    vtp_build_types(lines, 1, conn, offsets, types);
    vtp_build_types(polys, 2, conn, offsets, types);
    meta.mCellBlocks = detail::summarize_cells(offsets, types);

    meta.mPointDataNames = vtp_array_names(h.mPiece, "PointData");
    meta.mCellDataNames = vtp_array_names(h.mPiece, "CellData");
    meta.mFieldDataNames = vtp_field_data_names(h);
    {
        std::vector<std::pair<std::string, std::size_t>> arrays;
        for (const pugi::xml_node& rNode : {h.mGrid, h.mPiece})
            for (pugi::xml_node da : rNode.child("FieldData").children("DataArray"))
                arrays.emplace_back(da.attribute("Name").as_string(),
                                    da.attribute("NumberOfTuples").as_ullong(0));
        meta.mRegions = detail::region_summaries_from_field_names(arrays);
    }

    // No bbox: it would mean decoding the point coordinates. See read_options.hpp.
    meta.mHasBBox = false;
    return meta;
}

}  // namespace meshioplusplus
