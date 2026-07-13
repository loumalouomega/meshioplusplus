#pragma once
//
// CGNS (.cgns) reader/writer, mirroring meshio's h5py-based dialect
// (Base/Zone1/GridCoordinates + GridElements, datasets named " data",
// tetrahedra only). Available only when built with MESHIO_HAS_HDF5.

#ifdef MESHIO_HAS_HDF5

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_cgns(const std::string& path, const Mesh& mesh, int gzip_level);
Mesh read_cgns(const std::string& path);

}  // namespace meshio

#endif  // MESHIO_HAS_HDF5
