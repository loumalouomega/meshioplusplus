#pragma once
//
// FreeFem++ .msh reader/writer. ASCII: header `nver n_el1 n_el2`, then vertices
// (coords + ref), volume elements (triangle in 2D / tetra in 3D) and boundary
// elements (line in 2D / triangle in 3D), each 1-based and carrying an integer
// reference exposed as point_data/cell_data "freefem:ref".

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_freefem(const std::string& path, const Mesh& mesh);
Mesh read_freefem(const std::string& path);

}  // namespace meshioplusplus
