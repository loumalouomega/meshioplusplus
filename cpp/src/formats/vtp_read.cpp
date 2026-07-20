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
#include "meshioplusplus/formats/vtp.hpp"

namespace meshioplusplus {

namespace {

using detail::vtu_to_int64;

// The codec is resolved once from the root's compressor= attribute.
NDArray vtp_read_data_array(const pugi::xml_node& rDa, detail::VtkCodec codec, std::size_t hsz,
                            int& rNumComponents) {
    std::string fmt = rDa.attribute("format").as_string("ascii");
    DType dt = detail::dtype_from_vtu(rDa.attribute("type").as_string());
    rNumComponents = rDa.attribute("NumberOfComponents").as_int(0);

    if (fmt == "ascii")
        return detail::vtu_parse_ascii(rDa.text().get(), dt);
    if (fmt == "binary")
        return detail::vtu_parse_binary(detail::vtu_strip(rDa.text().get()), dt, codec, hsz);
    throw ReadError("VTP '" + fmt + "' data is not supported by the C++ reader");
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
VtpPiece vtp_read_section(const pugi::xml_node& rSection, detail::VtkCodec codec, std::size_t hsz,
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
        NDArray arr = vtp_read_data_array(da, codec, hsz, nc);
        if (name == "connectivity")
            out.mConn = vtu_to_int64(arr);
        else if (name == "offsets")
            out.mOffsets = vtu_to_int64(arr);
    }
    return out;
}

/** @brief `<Piece>` plus the framing attributes; mirrors `vtu_parse_header`. */
struct vtp_header {
    pugi::xml_node mPiece;
    detail::VtkCodec mCodec = detail::VtkCodec::None;
    std::size_t mHeaderSize = 4;
    std::size_t mNumPoints = 0;
};

vtp_header vtp_parse_header(const pugi::xml_document& rDoc) {
    pugi::xml_node root = rDoc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile'");
    if (std::string(root.attribute("type").as_string()) != "PolyData")
        throw ReadError("Expected type PolyData");

    vtp_header h;
    const std::string compressor = root.attribute("compressor").as_string("");
    if (compressor.empty())
        h.mCodec = detail::VtkCodec::None;
    else if (compressor == detail::vtk_codec_compressor(detail::VtkCodec::Zlib))
        h.mCodec = detail::VtkCodec::Zlib;
    else if (compressor == detail::vtk_codec_compressor(detail::VtkCodec::LZ4))
        h.mCodec = detail::VtkCodec::LZ4;
    else if (compressor == detail::vtk_codec_compressor(detail::VtkCodec::ZSTD))
        h.mCodec = detail::VtkCodec::ZSTD;
    else if (compressor == detail::vtk_codec_compressor(detail::VtkCodec::LZMA))
        throw ReadError("lzma-compressed VTP not supported by the C++ reader");
    else
        throw ReadError("Unknown VTP compressor '" + compressor + "'");
    // Fail early and actionably when the file needs a codec this build lacks,
    // rather than at the first array body.
    detail::vtk_codec_require_read(h.mCodec);

    std::string header_type = root.attribute("header_type").as_string("UInt32");
    h.mHeaderSize = (header_type == "UInt64") ? 8 : 4;

    pugi::xml_node grid = root.child("PolyData");
    if (!grid)
        throw ReadError("No PolyData found");

    // Appended data is not handled here -> let the Python reader take over.
    if (grid.parent().child("AppendedData") || root.child("AppendedData"))
        throw ReadError("appended VTP data not supported by the C++ reader");

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
    pugi::xml_document doc;
    pugi::xml_parse_result res = doc.load_file(rPath.c_str());
    if (!res)
        throw ReadError(std::string("VTP XML parse failed: ") + res.description());

    const vtp_header h = vtp_parse_header(doc);
    const pugi::xml_node piece = h.mPiece;
    const detail::VtkCodec codec = h.mCodec;
    const std::size_t hsz = h.mHeaderSize;
    const std::size_t num_points = h.mNumPoints;
    const bool want_data = rOpts.WantsAnyData();

    Mesh mesh;
    std::unordered_map<std::string, NDArray> cell_data_raw;

    for (pugi::xml_node child : piece.children()) {
        std::string tag = child.name();
        if (tag == "Points") {
            pugi::xml_node da = child.child("DataArray");
            int nc = 0;
            NDArray pts = vtp_read_data_array(da, codec, hsz, nc);
            if (nc <= 0)
                nc = 3;
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
                NDArray arr = vtp_read_data_array(da, codec, hsz, nc);
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
                NDArray arr = vtp_read_data_array(da, codec, hsz, nc);
                if (nc > 1)
                    arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                cell_data_raw.emplace(name, std::move(arr));
            }
        }
    }

    VtpPiece verts = vtp_read_section(piece.child("Verts"), codec, hsz);
    VtpPiece lines = vtp_read_section(piece.child("Lines"), codec, hsz);
    VtpPiece polys = vtp_read_section(piece.child("Polys"), codec, hsz);
    VtpPiece strips = vtp_read_section(piece.child("Strips"), codec, hsz);
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

    detail::reconstruct_cells(conn.data(), offsets, types, cell_data_raw, mesh);
    return mesh;
}

MeshMetadata read_vtp_metadata(const std::string& rPath, const ReadOptions&) {
    pugi::xml_document doc;
    // See read_vtu_metadata: parse_minimal trims text conversions, but the
    // saving that matters is skipping the array bodies below.
    pugi::xml_parse_result res = doc.load_file(rPath.c_str(), pugi::parse_minimal);
    if (!res)
        throw ReadError(std::string("VTP XML parse failed: ") + res.description());

    const vtp_header h = vtp_parse_header(doc);

    MeshMetadata meta;
    meta.mNumPoints = h.mNumPoints;  // an attribute -- free

    pugi::xml_node points_da = h.mPiece.child("Points").child("DataArray");
    const int point_nc = points_da ? points_da.attribute("NumberOfComponents").as_int(0) : 0;
    meta.mPointDim = point_nc > 0 ? static_cast<std::size_t>(point_nc) : 3;

    // PolyData carries no `types` array; cell types follow from each section's
    // per-cell size, so reading `offsets` alone suffices and the connectivity --
    // the bulk of the bytes -- is never decoded.
    const VtpPiece verts = vtp_read_section(h.mPiece.child("Verts"), h.mCodec, h.mHeaderSize,
                                            /*offsets_only=*/true);
    const VtpPiece lines = vtp_read_section(h.mPiece.child("Lines"), h.mCodec, h.mHeaderSize,
                                            /*offsets_only=*/true);
    const VtpPiece polys = vtp_read_section(h.mPiece.child("Polys"), h.mCodec, h.mHeaderSize,
                                            /*offsets_only=*/true);
    if (h.mPiece.child("Strips") &&
        !vtp_read_section(h.mPiece.child("Strips"), h.mCodec, h.mHeaderSize,
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

    // No bbox: it would mean decoding the point coordinates. See read_options.hpp.
    meta.mHasBBox = false;
    return meta;
}

}  // namespace meshioplusplus
