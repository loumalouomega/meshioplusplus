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
// VTK XML ImageData. Reader, writer and metadata reader in one file, unlike VTU
// (whose two halves are large enough to want separating): the geometry here is
// three attributes, so all three entry points fit comfortably together.
//
// Anonymous-namespace helpers are prefixed `vti_` -- the amalgamation
// concatenates every translation unit, and `vtu_` is taken.

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/formats/vti.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::vtu_ascii_ndarray;
using detail::vtu_type_str;

// Parse a whitespace-separated run of N numbers from an XML attribute. VTK's own
// files use plain spaces; being liberal here costs nothing and a mis-parse would
// silently relocate the whole grid.
template <class T>
bool vti_parse_n(const char* pText, T* pOut, std::size_t Count) {
    if (pText == nullptr)
        return false;
    std::istringstream is(pText);
    for (std::size_t i = 0; i < Count; ++i)
        if (!(is >> pOut[i]))
            return false;
    return true;
}

// The framing every ImageData path needs, resolved once so the mesh reader and
// the metadata reader cannot disagree about which files they accept.
struct vti_header {
    pugi::xml_node mPiece;
    detail::VtkCodec mCodec = detail::VtkCodec::None;
    std::size_t mHeaderSize = 4;
    detail::LatticeSpec mSpec;
    std::size_t mNumPoints = 0;
    std::size_t mNumCells = 0;
};

vti_header vti_parse_header(const pugi::xml_document& rDoc) {
    pugi::xml_node root = rDoc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile'");
    if (std::string(root.attribute("type").as_string()) != "ImageData")
        throw ReadError("Expected type ImageData");

    vti_header h;
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
        throw ReadError("lzma-compressed VTI not supported by the C++ reader");
    else
        throw ReadError("Unknown VTI compressor '" + compressor + "'");
    // Fail early and actionably when the file needs a codec this build lacks.
    detail::vtk_codec_require_read(h.mCodec);

    const std::string header_type = root.attribute("header_type").as_string("UInt32");
    h.mHeaderSize = (header_type == "UInt64") ? 8 : 4;

    if (root.child("AppendedData"))
        throw ReadError("appended VTI data not supported by the C++ reader");

    pugi::xml_node grid = root.child("ImageData");
    if (!grid)
        throw ReadError("No ImageData found");

    std::int64_t whole[6] = {0, 0, 0, 0, 0, 0};
    if (!vti_parse_n(grid.attribute("WholeExtent").as_string(nullptr), whole, 6))
        throw ReadError("ImageData has no readable WholeExtent");
    double origin[3] = {0.0, 0.0, 0.0};
    double spacing[3] = {1.0, 1.0, 1.0};
    if (grid.attribute("Origin"))
        vti_parse_n(grid.attribute("Origin").as_string(), origin, 3);
    if (grid.attribute("Spacing"))
        vti_parse_n(grid.attribute("Spacing").as_string(), spacing, 3);
    if (grid.attribute("Direction")) {
        // A non-identity direction matrix rotates the lattice, which an
        // axis-aligned hexahedron grid cannot express without baking the
        // rotation into the coordinates -- a different mesh from the one the
        // file describes. Refuse rather than silently drop the rotation.
        double dir[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        if (vti_parse_n(grid.attribute("Direction").as_string(), dir, 9)) {
            const double id[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
            for (std::size_t i = 0; i < 9; ++i)
                if (dir[i] != id[i])
                    throw ReadError(
                        "VTI with a non-identity Direction is not supported by the "
                        "C++ reader");
        }
    }

    h.mPiece = grid.child("Piece");
    if (!h.mPiece)
        throw ReadError("No Piece found");
    if (h.mPiece.next_sibling("Piece"))
        throw ReadError("multi-piece VTI not supported by the C++ reader");
    // The piece's own extent may legally be a sub-box of the whole extent; the
    // arrays are then sized to the PIECE, and reading them against the whole
    // extent would be silently misaligned. One piece covering everything is what
    // every writer emits, so decline the rest by name.
    if (h.mPiece.attribute("Extent")) {
        std::int64_t piece[6] = {0, 0, 0, 0, 0, 0};
        if (vti_parse_n(h.mPiece.attribute("Extent").as_string(), piece, 6))
            for (std::size_t i = 0; i < 6; ++i)
                if (piece[i] != whole[i])
                    throw ReadError(
                        "VTI Piece Extent differs from WholeExtent; a partial piece "
                        "is not supported by the C++ reader");
    }

    for (std::size_t k = 0; k < 3; ++k) {
        const std::int64_t n = whole[2 * k + 1] - whole[2 * k];
        if (n < 0)
            throw ReadError("VTI WholeExtent is inverted on axis " + std::to_string(k));
        h.mSpec.mDims[k] = n;
        // The extent may start away from zero; the point at extent index i sits
        // at Origin + i * Spacing, so the mesh's own lo corner is offset by the
        // extent's start. Dropping that offset would translate the whole grid.
        h.mSpec.mOrigin[k] = origin[k] + static_cast<double>(whole[2 * k]) * spacing[k];
        h.mSpec.mSpacing[k] = spacing[k];
    }
    h.mNumPoints = static_cast<std::size_t>((h.mSpec.mDims[0] + 1) * (h.mSpec.mDims[1] + 1) *
                                            (h.mSpec.mDims[2] + 1));
    h.mNumCells = static_cast<std::size_t>(detail::lattice_num_cells(h.mSpec));
    return h;
}

NDArray vti_read_data_array(const pugi::xml_node& rDa, detail::VtkCodec codec, std::size_t hsz,
                            int& rNumComponents) {
    const std::string fmt = rDa.attribute("format").as_string("ascii");
    const DType dt = detail::dtype_from_vtu(rDa.attribute("type").as_string());
    rNumComponents = rDa.attribute("NumberOfComponents").as_int(0);
    if (fmt == "ascii")
        return detail::vtu_parse_ascii(rDa.text().get(), dt);
    if (fmt == "binary")
        return detail::vtu_parse_binary(detail::vtu_strip(rDa.text().get()), dt, codec, hsz);
    throw ReadError("VTI '" + fmt + "' data is not supported by the C++ reader");
}

// One geometry attribute value. `%.17g` rather than the stream's default six
// significant digits, which would lose ~10 digits of a real origin -- a grid
// placed 1e-7 off its own points, which nothing downstream would flag. 17 is the
// round-trip width for a double, and the identical spelling is what the numpy
// twin uses, so the two writers' attributes agree character for character.
std::string vti_num(double Value) {
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.17g", Value);
    return buf;
}

std::vector<std::string> vti_array_names(const pugi::xml_node& rPiece, const char* pSection) {
    std::vector<std::string> names;
    for (pugi::xml_node da : rPiece.child(pSection).children("DataArray"))
        names.emplace_back(da.attribute("Name").as_string());
    // The uniform mesh API hands back sorted names; match it so a summary and a
    // real read report data arrays in the same order.
    std::sort(names.begin(), names.end());
    return names;
}

}  // namespace

void write_vti(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    write_vti_codec(rPath, rMesh, binary, zlib ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
}

void write_vti_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                     detail::VtkCodec codec) {
    detail::LatticeSpec spec;
    if (!detail::lattice_from_mesh(rMesh, spec))
        throw WriteError(
            "ImageData is a regular lattice, and this mesh is not one: it needs exactly one "
            "hexahedron block whose points tile an axis-aligned box with uniform spacing. A "
            "partial grid (voxelize's 'surface'/'inside' fill, or an octree) cannot be written "
            "as .vti either -- write it as .vtu, which stores the cells explicitly.");
    if (binary && codec != detail::VtkCodec::None)
        detail::vtk_codec_require_write(codec);

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const char* fmt = binary ? "binary" : "ascii";
    auto da_header = [&](const char* type, const std::string& name, int ncomp) {
        os << "<DataArray type=\"" << type << "\" Name=\"" << name << "\"";
        if (ncomp > 0)
            os << " NumberOfComponents=\"" << ncomp << "\"";
        os << " format=\"" << fmt << "\">\n";
    };
    auto emit_bin = [&](const unsigned char* d, std::size_t n) {
        os << detail::vtu_encode_binary(d, n, binary ? codec : detail::VtkCodec::None) << "\n";
    };

    std::ostringstream ext;
    ext << "0 " << spec.mDims[0] << " 0 " << spec.mDims[1] << " 0 " << spec.mDims[2];

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\"";
    if (binary && codec != detail::VtkCodec::None)
        os << " compressor=\"" << detail::vtk_codec_compressor(codec) << "\"";
    os << ">\n";
    os << detail::provenance_render_xml_comment(detail::SlotTier::Block) << "\n";
    // Origin/Spacing/WholeExtent ARE the geometry: no Points section exists, and
    // that is the whole reason this format is worth having for a grid.
    os << "<ImageData WholeExtent=\"" << ext.str() << "\" Origin=\"" << vti_num(spec.mOrigin[0])
       << " " << vti_num(spec.mOrigin[1]) << " " << vti_num(spec.mOrigin[2]) << "\" Spacing=\""
       << vti_num(spec.mSpacing[0]) << " " << vti_num(spec.mSpacing[1]) << " "
       << vti_num(spec.mSpacing[2]) << "\">\n";
    os << "<Piece Extent=\"" << ext.str() << "\">\n";

    if (rMesh.NumPointData() != 0) {
        os << "<PointData>\n";
        for (const auto& name : rMesh.PointDataNames()) {
            const NDArray& d = rMesh.PointData(name);
            const int ncomp = (d.Shape().size() == 2) ? static_cast<int>(cols(d)) : 0;
            da_header(vtu_type_str(d.Dtype()), name, ncomp);
            if (binary)
                emit_bin(reinterpret_cast<const unsigned char*>(d.Data()), d.Nbytes());
            else
                vtu_ascii_ndarray(os, d);
            os << "</DataArray>\n";
        }
        os << "</PointData>\n";
    }

    if (rMesh.NumCellData() != 0) {
        os << "<CellData>\n";
        for (const auto& name : rMesh.CellDataNames()) {
            const std::size_t nblocks = rMesh.CellDataNumBlocks(name);
            if (nblocks == 0)
                continue;
            // A lattice has exactly one block, so there is nothing to
            // concatenate -- but iterate anyway rather than assume, since an
            // array that does not cover the block is a caller error worth not
            // writing silently truncated.
            const NDArray& first = rMesh.CellData(name, 0);
            const int ncomp = (first.Shape().size() == 2) ? static_cast<int>(cols(first)) : 0;
            da_header(vtu_type_str(first.Dtype()), name, ncomp);
            if (binary) {
                std::vector<unsigned char> buf;
                for (std::size_t bi = 0; bi < nblocks; ++bi) {
                    const NDArray& blk = rMesh.CellData(name, bi);
                    const auto* p = reinterpret_cast<const unsigned char*>(blk.Data());
                    buf.insert(buf.end(), p, p + blk.Nbytes());
                }
                emit_bin(buf.data(), buf.size());
            } else {
                for (std::size_t bi = 0; bi < nblocks; ++bi)
                    vtu_ascii_ndarray(os, rMesh.CellData(name, bi));
            }
            os << "</DataArray>\n";
        }
        os << "</CellData>\n";
    }

    os << "</Piece>\n</ImageData>\n</VTKFile>\n";
}

Mesh read_vti(const std::string& rPath, const ReadOptions& rOpts) {
    pugi::xml_document doc;
    const pugi::xml_parse_result res = doc.load_file(rPath.c_str());
    if (!res)
        throw ReadError(std::string("VTI XML parse failed: ") + res.description());

    const vti_header h = vti_parse_header(doc);

    // The extent is expanded into explicit points and hexahedra through the same
    // helper `grid()` and `voxelize()` use, so a .vti read and a grid() call of
    // the same shape produce byte-identical meshes.
    Mesh mesh = detail::lattice_build_mesh(h.mSpec);
    if (!rOpts.WantsAnyData())
        return mesh;

    for (pugi::xml_node da : h.mPiece.child("PointData").children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        if (!rOpts.WantsArray(name))
            continue;
        int nc = 0;
        NDArray arr = vti_read_data_array(da, h.mCodec, h.mHeaderSize, nc);
        if (nc > 1)
            arr.Reshape({arr.Size() / static_cast<std::size_t>(nc), static_cast<std::size_t>(nc)});
        if (arr.Size() != 0 && detail::rows(arr) != h.mNumPoints)
            throw ReadError("VTI point array '" + name + "' has " +
                            std::to_string(detail::rows(arr)) + " rows, but the extent has " +
                            std::to_string(h.mNumPoints) + " points");
        mesh.AddPointData(name, std::move(arr));
    }
    for (pugi::xml_node da : h.mPiece.child("CellData").children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        if (!rOpts.WantsArray(name))
            continue;
        int nc = 0;
        NDArray arr = vti_read_data_array(da, h.mCodec, h.mHeaderSize, nc);
        if (nc > 1)
            arr.Reshape({arr.Size() / static_cast<std::size_t>(nc), static_cast<std::size_t>(nc)});
        if (arr.Size() != 0 && detail::rows(arr) != h.mNumCells)
            throw ReadError("VTI cell array '" + name + "' has " +
                            std::to_string(detail::rows(arr)) + " rows, but the extent has " +
                            std::to_string(h.mNumCells) + " cells");
        if (h.mNumCells == 0)
            continue;  // no cell block to attach it to
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(arr));
        mesh.AddCellData(name, std::move(blocks));
    }
    return mesh;
}

MeshMetadata read_vti_metadata(const std::string& rPath, const ReadOptions&) {
    pugi::xml_document doc;
    // parse_minimal skips escape expansion over the base64 bodies. Unlike VTU's
    // metadata path this decodes NOTHING at all: the extent attribute alone
    // gives both counts.
    const pugi::xml_parse_result res = doc.load_file(rPath.c_str(), pugi::parse_minimal);
    if (!res)
        throw ReadError(std::string("VTI XML parse failed: ") + res.description());

    const vti_header h = vti_parse_header(doc);

    MeshMetadata meta;
    meta.mNumPoints = h.mNumPoints;
    meta.mPointDim = 3;
    if (h.mNumCells != 0) {
        CellBlockInfo info;
        info.mType = "hexahedron";
        info.mNumCells = h.mNumCells;
        info.mNodesPerCell = 8;
        info.mRagged = false;
        meta.mCellBlocks.push_back(std::move(info));
    }
    meta.mPointDataNames = vti_array_names(h.mPiece, "PointData");
    meta.mCellDataNames = vti_array_names(h.mPiece, "CellData");

    // The bounding box IS the extent here, so unlike every other native metadata
    // path this one can report it for free rather than declining.
    meta.mHasBBox = true;
    for (std::size_t k = 0; k < 3; ++k) {
        meta.mBBoxMin[k] = h.mSpec.mOrigin[k];
        meta.mBBoxMax[k] =
            h.mSpec.mOrigin[k] + static_cast<double>(h.mSpec.mDims[k]) * h.mSpec.mSpacing[k];
    }
    return meta;
}

}  // namespace meshioplusplus
