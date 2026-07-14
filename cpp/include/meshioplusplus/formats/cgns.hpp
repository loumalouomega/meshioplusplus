#pragma once
//
// CGNS (.cgns) reader/writer, mirroring meshio's h5py-based dialect
// (Base/Zone1/GridCoordinates + GridElements, datasets named " data",
// tetrahedra only). Available only when built with MESHIOPLUSPLUS_HAS_HDF5.

#ifdef MESHIOPLUSPLUS_HAS_HDF5

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_cgns(const std::string& path, const Mesh& mesh, int gzip_level);
Mesh read_cgns(const std::string& path);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
