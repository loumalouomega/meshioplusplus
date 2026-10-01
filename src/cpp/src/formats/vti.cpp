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
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/vtk_xml_read.hpp"
#include "../detail/vtk_xml_pieces.hpp"
#include "../detail/text_cursor.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::vtu_ascii_ndarray;
using detail::vtu_type_str;

// One geometry attribute value. `%.17g` rather than the stream's default six
// significant digits, which would lose ~10 digits of a real origin -- a grid
// placed 1e-7 off its own points, which nothing downstream would flag. 17 is the
// round-trip width for a double, and the identical spelling is what the numpy
// twin uses, so the two writers' attributes agree character for character.
std::string vti_num(double Value) {
    char buf[40];
    detail::snprintf_c(buf, sizeof(buf), "%.17g", Value);
    return buf;
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
    os << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\"";
    if (binary && codec != detail::VtkCodec::None)
        os << " compressor=\"" << detail::vtk_codec_compressor(codec) << "\"";
    if (hsz == 8)
        os << " header_type=\"UInt64\"";
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
    return detail::vtk_xml_read_pieces(rPath, rOpts, "ImageData", "VTI");
}


MeshMetadata read_vti_metadata(const std::string& rPath, const ReadOptions&) {
    return detail::vtk_xml_pieces_metadata(rPath, "ImageData", "VTI");
}

}  // namespace meshioplusplus
