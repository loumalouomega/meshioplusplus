#pragma once
//
// HMF (.hmf) reader/writer — meshio's experimental HDF5 container using the
// XDMF topology names: /domain/grid with Topology{k} datasets, a Geometry
// dataset, and NodeAttributes/CellAttributes groups. Available only when
// built with MESHIO_HAS_HDF5.

#ifdef MESHIO_HAS_HDF5

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_hmf(const std::string& path, const Mesh& mesh, int gzip_level);
Mesh read_hmf(const std::string& path);

}  // namespace meshio

#endif  // MESHIO_HAS_HDF5
