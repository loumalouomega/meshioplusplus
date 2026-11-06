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
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/vtk_xml_pieces.hpp"
#include "../detail/vtk_xml_write.hpp"

namespace meshioplusplus {

namespace {

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

    const std::size_t hsz =
        detail::vtk_xml_lattice_preamble(os, rMesh, spec, "RectilinearGrid", "vtr", binary, codec);
    const detail::VtkXmlEmitter emit{os, binary, codec, hsz};

    os << "<Coordinates>\n";
    const char* axis_names[3] = {"x_coordinates", "y_coordinates", "z_coordinates"};
    for (std::size_t k = 0; k < 3; ++k) {
        const std::int64_t n = spec.mDims[k] + 1;
        std::vector<double> axis(static_cast<std::size_t>(n));
        for (std::int64_t i = 0; i < n; ++i)
            axis[static_cast<std::size_t>(i)] =
                spec.mOrigin[k] + static_cast<double>(i) * spec.mSpacing[k];
        emit.Header(vtu_type_str(DType::Float64), axis_names[k], 0);
        if (binary)
            emit.Bin(reinterpret_cast<const unsigned char*>(axis.data()),
                     axis.size() * sizeof(double));
        else
            for (double v : axis)
                os << v << "\n";
        os << "</DataArray>\n";
    }
    os << "</Coordinates>\n";

    detail::vtk_xml_lattice_data(emit, rMesh);

    os << "</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
}

Mesh read_vtr(const std::string& rPath, const ReadOptions& rOpts) {
    return detail::vtk_xml_read_pieces(rPath, rOpts, "RectilinearGrid", "VTR");
}

MeshMetadata read_vtr_metadata(const std::string& rPath, const ReadOptions&) {
    return detail::vtk_xml_pieces_metadata(rPath, "RectilinearGrid", "VTR");
}

}  // namespace meshioplusplus
