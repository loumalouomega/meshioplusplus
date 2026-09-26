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
// VTK XML StructuredGrid. Reader, writer and metadata reader in one file, as
// .vti's: the geometry logic here is small enough that splitting would only
// add link noise.
//
// Anonymous-namespace helpers are prefixed `vts_` -- the amalgamation
// concatenates every translation unit, and `vtu_`/`vti_` are taken.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/formats/vts.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "vtk_preflight.hpp"
#include "../detail/vtu_decode.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::vtu_ascii_ndarray;
using detail::vtu_type_str;

// Parse a whitespace-separated run of N numbers from an XML attribute.
template <class T>
bool vts_parse_n(const char* pText, T* pOut, std::size_t Count) {
    if (pText == nullptr)
        return false;
    auto is = detail::make_classic_istringstream(pText);
    for (std::size_t i = 0; i < Count; ++i)
        if (!(is >> pOut[i]))
            return false;
    return true;
}

// The framing every StructuredGrid path needs, resolved once so the mesh
// reader and the metadata reader cannot disagree about which files they
// accept. `mDims` are CELL counts per axis (as `LatticeSpec` uses them);
// `mOrigin`/`mSpacing` are unused here (points are explicit) but the type is
// shared with `.vti` for the writer's `lattice_from_mesh` call.
struct vts_header {
    pugi::xml_node mPiece;
    detail::VtkCodec mCodec = detail::VtkCodec::None;
    std::size_t mHeaderSize = 4;
    std::array<std::int64_t, 3> mDims{{0, 0, 0}};
    std::size_t mNumPoints = 0;
    std::size_t mNumCells = 0;
};

vts_header vts_parse_header(const pugi::xml_document& rDoc) {
    pugi::xml_node root = rDoc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile'");
    if (std::string(root.attribute("type").as_string()) != "StructuredGrid")
        throw ReadError("Expected type StructuredGrid");

    vts_header h;
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
        throw ReadError("lzma-compressed VTS not supported by the C++ reader");
    else
        throw ReadError("Unknown VTS compressor '" + compressor + "'");
    detail::vtk_codec_require_read(h.mCodec);

    const std::string header_type = root.attribute("header_type").as_string("UInt32");
    h.mHeaderSize = (header_type == "UInt64") ? 8 : 4;

    if (root.child("AppendedData"))
        throw ReadError("appended VTS data not supported by the C++ reader");

    pugi::xml_node grid = root.child("StructuredGrid");
    if (!grid)
        throw ReadError("No StructuredGrid found");

    std::int64_t whole[6] = {0, 0, 0, 0, 0, 0};
    if (!vts_parse_n(grid.attribute("WholeExtent").as_string(nullptr), whole, 6))
        throw ReadError("StructuredGrid has no readable WholeExtent");

    h.mPiece = grid.child("Piece");
    if (!h.mPiece)
        throw ReadError("No Piece found");
    if (h.mPiece.next_sibling("Piece"))
        throw ReadError("multi-piece VTS not supported by the C++ reader");
    if (h.mPiece.attribute("Extent")) {
        std::int64_t piece[6] = {0, 0, 0, 0, 0, 0};
        if (vts_parse_n(h.mPiece.attribute("Extent").as_string(), piece, 6))
            for (std::size_t i = 0; i < 6; ++i)
                if (piece[i] != whole[i])
                    throw ReadError(
                        "VTS Piece Extent differs from WholeExtent; a partial piece "
                        "is not supported by the C++ reader");
    }

    for (std::size_t k = 0; k < 3; ++k) {
        const std::int64_t n = whole[2 * k + 1] - whole[2 * k];
        if (n < 0)
            throw ReadError("VTS WholeExtent is inverted on axis " + std::to_string(k));
        h.mDims[k] = n;
    }
    h.mNumPoints = static_cast<std::size_t>((h.mDims[0] + 1) * (h.mDims[1] + 1) * (h.mDims[2] + 1));
    h.mNumCells = static_cast<std::size_t>(h.mDims[0] * h.mDims[1] * h.mDims[2]);
    return h;
}

NDArray vts_read_data_array(const pugi::xml_node& rDa, detail::VtkCodec codec, std::size_t hsz,
                            int& rNumComponents) {
    const std::string fmt = rDa.attribute("format").as_string("ascii");
    const DType dt = detail::dtype_from_vtu(rDa.attribute("type").as_string());
    rNumComponents = rDa.attribute("NumberOfComponents").as_int(0);
    if (fmt == "ascii")
        return detail::vtu_parse_ascii(rDa.text().get(), dt);
    if (fmt == "binary")
        return detail::vtu_decode_bin_view(detail::vtu_strip_view(rDa.text().get()), dt, codec,
                                           hsz);
    throw ReadError("VTS '" + fmt + "' data is not supported by the C++ reader");
}

std::vector<std::string> vts_array_names(const pugi::xml_node& rPiece, const char* pSection) {
    std::vector<std::string> names;
    for (pugi::xml_node da : rPiece.child(pSection).children("DataArray"))
        names.emplace_back(da.attribute("Name").as_string());
    std::sort(names.begin(), names.end());
    return names;
}

// hexahedron connectivity for cell (i, j, k), the same index formula
// `detail/grid_lattice.hpp` uses -- points come from the file, not from a
// recomputed origin/spacing, but the CONNECTIVITY formula is identical
// regardless of where the points came from.
void vts_hex_conn(std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t px, std::int64_t py,
                  std::int64_t* pOut) {
    const std::int64_t base = (k * py + j) * px + i;
    const std::int64_t top = base + px * py;
    pOut[0] = base;
    pOut[1] = base + 1;
    pOut[2] = base + px + 1;
    pOut[3] = base + px;
    pOut[4] = top;
    pOut[5] = top + 1;
    pOut[6] = top + px + 1;
    pOut[7] = top + px;
}

}  // namespace

void write_vts(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    write_vts_codec(rPath, rMesh, binary, zlib ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
}

void write_vts_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                     detail::VtkCodec codec) {
    detail::LatticeSpec spec;
    if (!detail::lattice_from_mesh(rMesh, spec))
        throw WriteError(
            "StructuredGrid is a lattice, and this mesh is not one: it needs exactly one "
            "hexahedron block whose points tile an axis-aligned box with uniform spacing "
            "(the writer recovers WholeExtent from it, then writes the mesh's own points "
            "unchanged). A partial grid (voxelize's 'surface'/'inside' fill, or an octree) "
            "cannot be written as .vts either -- write it as .vtu, which stores the cells "
            "explicitly.");
    if (binary && codec != detail::VtkCodec::None)
        detail::vtk_codec_require_write(codec);

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const char* fmt = binary ? "binary" : "ascii";
    // UInt64 size headers only where an uncompressed array could pass 4 GiB
    // (compressed ones count 32 KiB blocks); everything else keeps its bytes.
    const std::size_t hsz =
        (binary && codec == detail::VtkCodec::None) ? detail::vtk_xml_header_bytes(rMesh) : 4;
    auto da_header = [&](const char* type, const std::string& name, int ncomp) {
        os << "<DataArray type=\"" << type << "\" Name=\"" << name << "\"";
        if (ncomp > 0)
            os << " NumberOfComponents=\"" << ncomp << "\"";
        os << " format=\"" << fmt << "\">\n";
    };
    auto emit_bin = [&](const unsigned char* d, std::size_t n) {
        os << detail::vtu_encode_binary(d, n, binary ? codec : detail::VtkCodec::None, hsz) << "\n";
    };

    auto ext = detail::make_classic_ostringstream();
    ext << "0 " << spec.mDims[0] << " 0 " << spec.mDims[1] << " 0 " << spec.mDims[2];

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"";
    if (binary && codec != detail::VtkCodec::None)
        os << " compressor=\"" << detail::vtk_codec_compressor(codec) << "\"";
    if (hsz == 8)
        os << " header_type=\"UInt64\"";
    os << ">\n";
    os << detail::provenance_render_xml_comment(detail::SlotTier::Block) << "\n";
    os << "<StructuredGrid WholeExtent=\"" << ext.str() << "\">\n";
    os << "<Piece Extent=\"" << ext.str() << "\">\n";

    const NDArray& pts = rMesh.Points();
    const std::size_t np = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    os << "<Points>\n";
    da_header(vtu_type_str(DType::Float64), "Points", 3);
    if (binary) {
        std::vector<double> buf(np * 3, 0.0);
        detail::dispatch_dtype(pts.Dtype(), [&]<class T>() {
            const T* src = pts.As<T>();
            parallel_for_bw(np, [&](std::size_t r) {
                for (std::size_t c = 0; c < pdim && c < 3; ++c)
                    buf[r * 3 + c] = static_cast<double>(src[r * pdim + c]);
            });
        });
        emit_bin(reinterpret_cast<const unsigned char*>(buf.data()), buf.size() * sizeof(double));
    } else {
        for (std::size_t r = 0; r < np; ++r) {
            for (std::size_t c = 0; c < 3; ++c)
                os << (c ? " " : "") << (c < pdim ? detail::read_double(pts, r * pdim + c) : 0.0);
            os << "\n";
        }
    }
    os << "</DataArray>\n</Points>\n";

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

    os << "</Piece>\n</StructuredGrid>\n</VTKFile>\n";
}

Mesh read_vts(const std::string& rPath, const ReadOptions& rOpts) {
    pugi::xml_document doc;
    detail::vtk_preflight(rPath, "StructuredGrid",
                          "lzma-compressed VTS not supported by the C++ reader");
    const pugi::xml_parse_result res = doc.load_file(rPath.c_str());
    if (!res)
        throw ReadError(std::string("VTS XML parse failed: ") + res.description());

    const vts_header h = vts_parse_header(doc);

    pugi::xml_node points_da = h.mPiece.child("Points").child("DataArray");
    if (!points_da)
        throw ReadError("VTS Piece has no Points/DataArray");
    int pnc = 0;
    NDArray pts = vts_read_data_array(points_da, h.mCodec, h.mHeaderSize, pnc);
    if (pnc > 1)
        pts.Reshape({pts.Size() / static_cast<std::size_t>(pnc), static_cast<std::size_t>(pnc)});
    if (detail::rows(pts) != h.mNumPoints)
        throw ReadError("VTS Points has " + std::to_string(detail::rows(pts)) +
                        " rows, but WholeExtent has " + std::to_string(h.mNumPoints) + " points");

    Mesh mesh;
    mesh.AssignPoints(std::move(pts));

    if (h.mNumCells != 0) {
        const std::int64_t px = h.mDims[0] + 1;
        const std::int64_t py = h.mDims[1] + 1;
        NDArray conn = NDArray::Uninit(DType::Int64, {h.mNumCells, std::size_t{8}});
        std::int64_t* dst = conn.As<std::int64_t>();
        std::size_t c = 0;
        for (std::int64_t k = 0; k < h.mDims[2]; ++k)
            for (std::int64_t j = 0; j < h.mDims[1]; ++j)
                for (std::int64_t i = 0; i < h.mDims[0]; ++i, ++c)
                    vts_hex_conn(i, j, k, px, py, dst + c * 8);
        mesh.AddCellBlock("hexahedron", std::move(conn));
    }

    if (!rOpts.WantsAnyData())
        return mesh;

    for (pugi::xml_node da : h.mPiece.child("PointData").children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        if (!rOpts.WantsArray(name))
            continue;
        int nc = 0;
        NDArray arr = vts_read_data_array(da, h.mCodec, h.mHeaderSize, nc);
        if (nc > 1)
            arr.Reshape({arr.Size() / static_cast<std::size_t>(nc), static_cast<std::size_t>(nc)});
        if (arr.Size() != 0 && detail::rows(arr) != h.mNumPoints)
            throw ReadError("VTS point array '" + name + "' has " +
                            std::to_string(detail::rows(arr)) + " rows, but the extent has " +
                            std::to_string(h.mNumPoints) + " points");
        mesh.AddPointData(name, std::move(arr));
    }
    for (pugi::xml_node da : h.mPiece.child("CellData").children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        if (!rOpts.WantsArray(name))
            continue;
        int nc = 0;
        NDArray arr = vts_read_data_array(da, h.mCodec, h.mHeaderSize, nc);
        if (nc > 1)
            arr.Reshape({arr.Size() / static_cast<std::size_t>(nc), static_cast<std::size_t>(nc)});
        if (arr.Size() != 0 && detail::rows(arr) != h.mNumCells)
            throw ReadError("VTS cell array '" + name + "' has " +
                            std::to_string(detail::rows(arr)) + " rows, but the extent has " +
                            std::to_string(h.mNumCells) + " cells");
        if (h.mNumCells == 0)
            continue;
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(arr));
        mesh.AddCellData(name, std::move(blocks));
    }
    return mesh;
}

MeshMetadata read_vts_metadata(const std::string& rPath, const ReadOptions&) {
    pugi::xml_document doc;
    detail::vtk_preflight(rPath, "StructuredGrid",
                          "lzma-compressed VTS not supported by the C++ reader");
    const pugi::xml_parse_result res = doc.load_file(rPath.c_str(), pugi::parse_minimal);
    if (!res)
        throw ReadError(std::string("VTS XML parse failed: ") + res.description());

    const vts_header h = vts_parse_header(doc);

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
    meta.mPointDataNames = vts_array_names(h.mPiece, "PointData");
    meta.mCellDataNames = vts_array_names(h.mPiece, "CellData");
    return meta;
}

}  // namespace meshioplusplus
