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
//                   SPDX-License-Identifier: MIT
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
#include "../detail/structured_grid.hpp"
#include "../detail/vtk_xml_pieces.hpp"
#include "../detail/vtk_xml_write.hpp"

namespace meshioplusplus {

namespace {

using detail::vtu_type_str;

}  // namespace

void write_vtr(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    write_vtr_codec(rPath, rMesh, binary, detail::vtk_codec_from_zlib(zlib));
}

void write_vtr_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                     detail::VtkCodec codec) {
    std::array<std::int64_t, 3> dims{};
    std::vector<double> xs, ys, zs;
    if (!detail::rectilinear_from_mesh(rMesh, dims, xs, ys, zs))
        throw WriteError(
            "RectilinearGrid needs a dense rectilinear grid: one hexahedron block whose "
            "points lie on strictly increasing axis-aligned planes (either annotated with "
            "field_data['meshio:structured:0'] or geometrically uniform). Write unstructured "
            "or deformed meshes as .vtu.");
    const std::vector<double>* axes[3] = {&xs, &ys, &zs};
    detail::vtk_xml_write_structured_dims(
        rPath, rMesh, dims, binary, codec, "RectilinearGrid", "vtr",
        [&](const detail::VtkXmlEmitter& rEmit, const detail::LatticeSpec&) {
            std::ostream& os = rEmit.mOs;
            os << "<Coordinates>\n";
            const char* axis_names[3] = {"x_coordinates", "y_coordinates", "z_coordinates"};
            for (std::size_t k = 0; k < 3; ++k) {
                const std::vector<double>& axis = *axes[k];
                rEmit.Header(vtu_type_str(DType::Float64), axis_names[k], 0);
                if (binary)
                    rEmit.Bin(reinterpret_cast<const unsigned char*>(axis.data()),
                              axis.size() * sizeof(double));
                else
                    for (double v : axis)
                        os << v << "\n";
                os << "</DataArray>\n";
            }
            os << "</Coordinates>\n";
        });
}

Mesh read_vtr(const std::string& rPath, const ReadOptions& rOpts) {
    return detail::vtk_xml_read_pieces(rPath, rOpts, "RectilinearGrid", "VTR");
}

MeshMetadata read_vtr_metadata(const std::string& rPath, const ReadOptions&) {
    return detail::vtk_xml_pieces_metadata(rPath, "RectilinearGrid", "VTR");
}

}  // namespace meshioplusplus
