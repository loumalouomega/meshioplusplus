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
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/region_field_data.hpp"
#include "../detail/vtk_xml_pieces.hpp"
#include "../detail/typed_view.hpp"
#include "../detail/vtk_xml_write.hpp"
#include "../detail/vtu_decode.hpp"

namespace meshioplusplus {

namespace {

using detail::vtu_type_str;

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

    const std::size_t hsz = detail::vtk_xml_size_header_bytes(rMesh, binary, codec);
    const detail::VtkXmlEmitter emit{os, binary, codec, hsz};
    const std::string ext = detail::vtk_xml_lattice_extent(spec);

    detail::vtk_xml_prolog(os, "StructuredGrid", "0.1", binary, codec, hsz);
    os << "<StructuredGrid WholeExtent=\"" << ext << "\">\n";
    detail::write_field_data_section(os, rMesh, "vts", binary, codec, hsz);
    os << "<Piece Extent=\"" << ext << "\">\n";

    const NDArray& pts = rMesh.Points();
    const std::size_t np = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    os << "<Points>\n";
    emit.Header(vtu_type_str(DType::Float64), "Points", 3);
    if (binary) {
        std::vector<double> buf(np * 3, 0.0);
        detail::dispatch_dtype(pts.Dtype(), [&]<class T>() {
            const T* src = pts.As<T>();
            parallel_for_bw(np, [&](std::size_t r) {
                for (std::size_t c = 0; c < pdim && c < 3; ++c)
                    buf[r * 3 + c] = static_cast<double>(src[r * pdim + c]);
            });
        });
        emit.Bin(reinterpret_cast<const unsigned char*>(buf.data()), buf.size() * sizeof(double));
    } else {
        const detail::DoubleView pts_values(pts);
        for (std::size_t r = 0; r < np; ++r) {
            for (std::size_t c = 0; c < 3; ++c)
                os << (c ? " " : "") << (c < pdim ? pts_values[r * pdim + c] : 0.0);
            os << "\n";
        }
    }
    os << "</DataArray>\n</Points>\n";

    detail::vtk_xml_lattice_data(emit, rMesh);

    os << "</Piece>\n</StructuredGrid>\n</VTKFile>\n";
}

Mesh read_vts(const std::string& rPath, const ReadOptions& rOpts) {
    return detail::vtk_xml_read_pieces(rPath, rOpts, "StructuredGrid", "VTS");
}

MeshMetadata read_vts_metadata(const std::string& rPath, const ReadOptions&) {
    return detail::vtk_xml_pieces_metadata(rPath, "StructuredGrid", "VTS");
}

}  // namespace meshioplusplus
