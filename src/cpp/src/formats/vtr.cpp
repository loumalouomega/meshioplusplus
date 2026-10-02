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
// VTK XML RectilinearGrid.
//
// Anonymous-namespace helpers are prefixed `vtr_` -- the amalgamation
// concatenates every translation unit, and `vtu_`/`vti_`/`vts_` are taken.

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/vtr.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/vtk_xml_pieces.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::vtu_ascii_ndarray;
using detail::vtu_type_str;

}  // namespace

void write_vtr(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    write_vtr_codec(rPath, rMesh, binary, zlib ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
}

void write_vtr_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                     detail::VtkCodec codec) {
    detail::LatticeSpec spec;
    if (!detail::lattice_from_mesh(rMesh, spec))
        throw WriteError(
            "RectilinearGrid needs a UNIFORM dense lattice today: exactly one hexahedron "
            "block whose points tile an axis-aligned box with uniform per-axis spacing. A "
            "genuinely graded (non-uniform) mesh cannot be written as .vtr yet -- a "
            "documented follow-up, see doc/roadmap.md -- and a partial grid (voxelize's "
            "'surface'/'inside' fill, or an octree) cannot be written as .vtr either -- "
            "write it as .vtu, which stores the cells explicitly.");
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
    os << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\"";
    if (binary && codec != detail::VtkCodec::None)
        os << " compressor=\"" << detail::vtk_codec_compressor(codec) << "\"";
    if (hsz == 8)
        os << " header_type=\"UInt64\"";
    os << ">\n";
    os << detail::provenance_render_xml_comment(detail::SlotTier::Block) << "\n";
    os << "<RectilinearGrid WholeExtent=\"" << ext.str() << "\">\n";
    os << "<Piece Extent=\"" << ext.str() << "\">\n";

    os << "<Coordinates>\n";
    const char* axis_names[3] = {"x_coordinates", "y_coordinates", "z_coordinates"};
    for (std::size_t k = 0; k < 3; ++k) {
        const std::int64_t n = spec.mDims[k] + 1;
        std::vector<double> axis(static_cast<std::size_t>(n));
        for (std::int64_t i = 0; i < n; ++i)
            axis[static_cast<std::size_t>(i)] =
                spec.mOrigin[k] + static_cast<double>(i) * spec.mSpacing[k];
        da_header(vtu_type_str(DType::Float64), axis_names[k], 0);
        if (binary)
            emit_bin(reinterpret_cast<const unsigned char*>(axis.data()),
                     axis.size() * sizeof(double));
        else
            for (double v : axis)
                os << v << "\n";
        os << "</DataArray>\n";
    }
    os << "</Coordinates>\n";

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

    os << "</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
}

Mesh read_vtr(const std::string& rPath, const ReadOptions& rOpts) {
    return detail::vtk_xml_read_pieces(rPath, rOpts, "RectilinearGrid", "VTR");
}

MeshMetadata read_vtr_metadata(const std::string& rPath, const ReadOptions&) {
    return detail::vtk_xml_pieces_metadata(rPath, "RectilinearGrid", "VTR");
}

}  // namespace meshioplusplus
