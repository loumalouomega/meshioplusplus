#pragma once
//
// MOAB h5m (.h5m) reader/writer, mirroring meshio's h5py dialect: the `tstt`
// root group with nodes/coordinates (+ tags), elements/<Type>/connectivity
// (1-based, enum-typed element_type attribute), committed tag datatypes, a
// history dataset and max_id. Write supports line/triangle/tetra (the meshio
// limitation). Available only when built with MESHIOPLUSPLUS_HAS_HDF5.

#ifdef MESHIOPLUSPLUS_HAS_HDF5

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_h5m(const std::string& path, const Mesh& mesh, bool add_global_ids,
               int gzip_level);
Mesh read_h5m(const std::string& path);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
