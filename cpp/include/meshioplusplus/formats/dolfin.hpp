#pragma once
//
// DOLFIN legacy XML (.xml) reader/writer. The mesh lives in one file
// (<dolfin><mesh><vertices/><cells/></mesh></dolfin>); each cell_data array is
// stored in a sibling file "<stem>_<name>.xml" holding a <mesh_function>.
// Only triangle and tetrahedron cells are supported (as in meshio).

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_dolfin(const std::string& path, const Mesh& mesh);
Mesh read_dolfin(const std::string& path);

}  // namespace meshioplusplus
