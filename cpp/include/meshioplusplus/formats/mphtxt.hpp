#pragma once
//
// COMSOL .mphtxt reader/writer. Token-stream ASCII: version, tag/type tables,
// then a mesh object (sdim, node coordinates, and hybrid element-type blocks
// with per-element geometric entity indices -> cell_data["mphtxt:geom"]).
// Comments run from '#' to end of line.

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_mphtxt(const std::string& path, const Mesh& mesh);
Mesh read_mphtxt(const std::string& path);

}  // namespace meshioplusplus
