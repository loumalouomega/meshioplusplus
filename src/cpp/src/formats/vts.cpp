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
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/vts.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/parallel.hpp"
#include "../detail/vtk_xml_pieces.hpp"
#include "../detail/typed_view.hpp"
#include "../detail/vtk_xml_write.hpp"

namespace meshioplusplus {

namespace {

using detail::vtu_type_str;

}  // namespace

void write_vts(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    write_vts_codec(rPath, rMesh, binary, detail::vtk_codec_from_zlib(zlib));
}

void write_vts_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                     detail::VtkCodec codec) {
    detail::vtk_xml_write_lattice(
        rPath, rMesh, binary, codec, "StructuredGrid", "vts",
        "StructuredGrid is a lattice, and this mesh is not one: it needs exactly one "
        "hexahedron block whose points tile an axis-aligned box with uniform spacing "
        "(the writer recovers WholeExtent from it, then writes the mesh's own points "
        "unchanged). A partial grid (voxelize's 'surface'/'inside' fill, or an octree) "
        "cannot be written as .vts either -- write it as .vtu, which stores the cells "
        "explicitly.",
        [&](const detail::VtkXmlEmitter& rEmit, const detail::LatticeSpec&) {
            std::ostream& os = rEmit.mOs;
            const NDArray& pts = rMesh.Points();
            const std::size_t np = rMesh.NumPoints();
            const std::size_t pdim = rMesh.PointDim();
            os << "<Points>\n";
            rEmit.Header(vtu_type_str(DType::Float64), "Points", 3);
            if (binary) {
                std::vector<double> buf(np * 3, 0.0);
                detail::dispatch_dtype(pts.Dtype(), [&]<class T>() {
                    const T* src = pts.As<T>();
                    parallel_for_bw(np, [&](std::size_t r) {
                        for (std::size_t c = 0; c < pdim && c < 3; ++c)
                            buf[r * 3 + c] = static_cast<double>(src[r * pdim + c]);
                    });
                });
                rEmit.Bin(reinterpret_cast<const unsigned char*>(buf.data()),
                          buf.size() * sizeof(double));
            } else {
                const detail::DoubleView pts_values(pts);
                for (std::size_t r = 0; r < np; ++r) {
                    for (std::size_t c = 0; c < 3; ++c)
                        os << (c ? " " : "") << (c < pdim ? pts_values[r * pdim + c] : 0.0);
                    os << "\n";
                }
            }
            os << "</DataArray>\n</Points>\n";
        });
}

Mesh read_vts(const std::string& rPath, const ReadOptions& rOpts) {
    return detail::vtk_xml_read_pieces(rPath, rOpts, "StructuredGrid", "VTS");
}

MeshMetadata read_vts_metadata(const std::string& rPath, const ReadOptions&) {
    return detail::vtk_xml_pieces_metadata(rPath, "StructuredGrid", "VTS");
}

}  // namespace meshioplusplus
