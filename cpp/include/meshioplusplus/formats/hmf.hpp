#pragma once
//
// HMF (.hmf) reader/writer — meshio's experimental HDF5 container using the
// XDMF topology names: /domain/grid with Topology{k} datasets, a Geometry
// dataset, and NodeAttributes/CellAttributes groups. Available only when
// built with MESHIOPLUSPLUS_HAS_HDF5.

#ifdef MESHIOPLUSPLUS_HAS_HDF5

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_hmf(const std::string& path, const Mesh& mesh, int gzip_level);
Mesh read_hmf(const std::string& path);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
