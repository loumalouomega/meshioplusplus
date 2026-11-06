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
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/vti.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/region_field_data.hpp"
#include "../detail/vtk_xml_pieces.hpp"
#include "../detail/vtk_xml_write.hpp"
#include "../detail/vtu_decode.hpp"

namespace meshioplusplus {

namespace {

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

    const std::size_t hsz = detail::vtk_xml_size_header_bytes(rMesh, binary, codec);
    const detail::VtkXmlEmitter emit{os, binary, codec, hsz};
    const std::string ext = detail::vtk_xml_lattice_extent(spec);

    detail::vtk_xml_prolog(os, "ImageData", "0.1", binary, codec, hsz);
    // Origin/Spacing/WholeExtent ARE the geometry: no Points section exists, and
    // that is the whole reason this format is worth having for a grid.
    os << "<ImageData WholeExtent=\"" << ext << "\" Origin=\"" << vti_num(spec.mOrigin[0]) << " "
       << vti_num(spec.mOrigin[1]) << " " << vti_num(spec.mOrigin[2]) << "\" Spacing=\""
       << vti_num(spec.mSpacing[0]) << " " << vti_num(spec.mSpacing[1]) << " "
       << vti_num(spec.mSpacing[2]) << "\">\n";
    detail::write_field_data_section(os, rMesh, "vti", binary, codec, hsz);
    os << "<Piece Extent=\"" << ext << "\">\n";

    detail::vtk_xml_lattice_data(emit, rMesh);

    os << "</Piece>\n</ImageData>\n</VTKFile>\n";
}

Mesh read_vti(const std::string& rPath, const ReadOptions& rOpts) {
    return detail::vtk_xml_read_pieces(rPath, rOpts, "ImageData", "VTI");
}

MeshMetadata read_vti_metadata(const std::string& rPath, const ReadOptions&) {
    return detail::vtk_xml_pieces_metadata(rPath, "ImageData", "VTI");
}

}  // namespace meshioplusplus
