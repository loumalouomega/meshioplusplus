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
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtu.hpp"

namespace meshioplusplus {

namespace {

using detail::vtu_to_int64;

// The codec is resolved once from the root's compressor= attribute.
NDArray vtu_read_data_array(const pugi::xml_node& rDa, detail::VtkCodec codec, std::size_t hsz,
                            int& rNumComponents) {
    std::string fmt = rDa.attribute("format").as_string("ascii");
    DType dt = detail::dtype_from_vtu(rDa.attribute("type").as_string());
    rNumComponents = rDa.attribute("NumberOfComponents").as_int(0);

    if (fmt == "ascii")
        return detail::vtu_parse_ascii(rDa.text().get(), dt);
    if (fmt == "binary")
        return detail::vtu_parse_binary(detail::vtu_strip(rDa.text().get()), dt, codec, hsz);
    throw ReadError("VTU '" + fmt + "' data is not supported by the C++ reader");
}

/**
 * @brief The `<Piece>` node plus the framing attributes every path needs.
 *
 * Shared by the mesh and metadata readers so the two cannot disagree about
 * which files they accept -- a metadata summary must never succeed on a file
 * `read_vtu` would reject.
 */
struct vtu_header {
    pugi::xml_node mPiece;
    detail::VtkCodec mCodec = detail::VtkCodec::None;
    std::size_t mHeaderSize = 4;
    std::size_t mNumPoints = 0;
};

vtu_header vtu_parse_header(const pugi::xml_document& rDoc) {
    pugi::xml_node root = rDoc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile'");
    if (std::string(root.attribute("type").as_string()) != "UnstructuredGrid")
        throw ReadError("Expected type UnstructuredGrid");

    vtu_header h;
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
        throw ReadError("lzma-compressed VTU not supported by the C++ reader");
    else
        throw ReadError("Unknown VTU compressor '" + compressor + "'");
    // Fail early and actionably when the file needs a codec this build lacks,
    // rather than at the first array body.
    detail::vtk_codec_require_read(h.mCodec);

    std::string header_type = root.attribute("header_type").as_string("UInt32");
    h.mHeaderSize = (header_type == "UInt64") ? 8 : 4;

    pugi::xml_node grid = root.child("UnstructuredGrid");
    if (!grid)
        throw ReadError("No UnstructuredGrid found");

    // Appended data is not handled here -> let the Python reader take over.
    if (grid.parent().child("AppendedData") || root.child("AppendedData"))
        throw ReadError("appended VTU data not supported by the C++ reader");

    h.mPiece = grid.child("Piece");
    if (!h.mPiece)
        throw ReadError("No Piece found");
    // A single piece is supported; multiple pieces -> Python reader.
    if (h.mPiece.next_sibling("Piece"))
        throw ReadError("multi-piece VTU not supported by the C++ reader");

    h.mNumPoints = static_cast<std::size_t>(h.mPiece.attribute("NumberOfPoints").as_ullong());
    return h;
}

/** @brief `<DataArray>` `Name` attributes under @p rSection, in document order. */
std::vector<std::string> vtu_array_names(const pugi::xml_node& rPiece, const char* pSection) {
    std::vector<std::string> names;
    for (pugi::xml_node da : rPiece.child(pSection).children("DataArray"))
        names.emplace_back(da.attribute("Name").as_string());
    // The uniform mesh API hands back sorted names; match it so a summary and a
    // real read report data arrays in the same order.
    std::sort(names.begin(), names.end());
    return names;
}

}  // namespace

Mesh read_vtu(const std::string& rPath, const ReadOptions& rOpts) {
    pugi::xml_document doc;
    pugi::xml_parse_result res = doc.load_file(rPath.c_str());
    if (!res)
        throw ReadError(std::string("VTU XML parse failed: ") + res.description());

    const vtu_header h = vtu_parse_header(doc);
    const pugi::xml_node piece = h.mPiece;
    const detail::VtkCodec codec = h.mCodec;
    const std::size_t hsz = h.mHeaderSize;
    const std::size_t num_points = h.mNumPoints;
    const bool want_data = rOpts.WantsAnyData();

    Mesh mesh;
    std::vector<std::int64_t> conn, offsets, types;
    // VTU's polyhedral stream; empty when the file has none.
    std::vector<std::int64_t> faces, face_offsets;
    std::unordered_map<std::string, NDArray> cell_data_raw;

    for (pugi::xml_node child : piece.children()) {
        std::string tag = child.name();
        if (tag == "Points") {
            pugi::xml_node da = child.child("DataArray");
            int nc = 0;
            NDArray pts = vtu_read_data_array(da, codec, hsz, nc);
            if (nc <= 0)
                nc = 3;
            pts.Reshape({num_points, static_cast<std::size_t>(nc)});
            mesh.AssignPoints(std::move(pts));
        } else if (tag == "Cells") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = vtu_read_data_array(da, codec, hsz, nc);
                if (name == "connectivity")
                    conn = vtu_to_int64(arr);
                else if (name == "offsets")
                    offsets = vtu_to_int64(arr);
                else if (name == "types")
                    types = vtu_to_int64(arr);
                else if (name == "faces")
                    faces = vtu_to_int64(arr);
                else if (name == "faceoffsets")
                    face_offsets = vtu_to_int64(arr);
            }
        } else if (tag == "PointData") {
            if (!want_data)
                continue;
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                // The Name attribute is available before the body is touched,
                // so an unwanted array costs nothing but the attribute read --
                // no base64 decode, no inflate, no allocation.
                if (!rOpts.WantsArray(name))
                    continue;
                NDArray arr = vtu_read_data_array(da, codec, hsz, nc);
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
                NDArray arr = vtu_read_data_array(da, codec, hsz, nc);
                if (nc > 1)
                    arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                cell_data_raw.emplace(name, std::move(arr));
            }
        }
    }

    detail::reconstruct_cells(conn.data(), offsets, types, cell_data_raw,
                              faces.empty() ? nullptr : &faces, face_offsets, mesh);
    return mesh;
}

MeshMetadata read_vtu_metadata(const std::string& rPath, const ReadOptions&) {
    pugi::xml_document doc;
    // parse_minimal skips escape expansion and EOL normalization over the
    // base64 bodies. It does NOT avoid reading the file: pugixml always
    // materializes PCDATA. The real saving below is skipping base64 decode,
    // decompression, allocation and byte-swapping for every array we don't
    // touch -- a solid multiple, not an asymptotic change. See read_options.hpp.
    pugi::xml_parse_result res = doc.load_file(rPath.c_str(), pugi::parse_minimal);
    if (!res)
        throw ReadError(std::string("VTU XML parse failed: ") + res.description());

    const vtu_header h = vtu_parse_header(doc);

    MeshMetadata meta;
    meta.mNumPoints = h.mNumPoints;  // an attribute -- free

    // Point dimension comes from the Points DataArray's NumberOfComponents
    // attribute, so the coordinates themselves are never decoded.
    pugi::xml_node points_da = h.mPiece.child("Points").child("DataArray");
    const int point_nc = points_da ? points_da.attribute("NumberOfComponents").as_int(0) : 0;
    meta.mPointDim = point_nc > 0 ? static_cast<std::size_t>(point_nc) : 3;

    // 'types' is one value per cell and is the only array a summary must decode;
    // 'offsets' is additionally needed only when a variable-node-count type
    // (polygon, VTK_LAGRANGE_*) is present.
    std::vector<std::int64_t> types, offsets;
    pugi::xml_node cells = h.mPiece.child("Cells");
    for (pugi::xml_node da : cells.children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        if (name != "types")
            continue;
        int nc = 0;
        types = vtu_to_int64(vtu_read_data_array(da, h.mCodec, h.mHeaderSize, nc));
        break;
    }
    if (detail::cells_need_offsets(types)) {
        for (pugi::xml_node da : cells.children("DataArray")) {
            const std::string name = da.attribute("Name").as_string();
            if (name != "offsets")
                continue;
            int nc = 0;
            offsets = vtu_to_int64(vtu_read_data_array(da, h.mCodec, h.mHeaderSize, nc));
            break;
        }
    }
    meta.mCellBlocks = detail::summarize_cells(offsets, types);

    meta.mPointDataNames = vtu_array_names(h.mPiece, "PointData");
    meta.mCellDataNames = vtu_array_names(h.mPiece, "CellData");
    // VTU has no field-data section the C++ reader consumes, so the list stays
    // empty rather than claiming an unknown.

    // No bounding box: it would require decoding the point coordinates, which
    // are usually the largest array in the file -- exactly what this path exists
    // to avoid. metadata_from_mesh fills it in on the full-read fallback.
    meta.mHasBBox = false;
    return meta;
}

}  // namespace meshioplusplus
