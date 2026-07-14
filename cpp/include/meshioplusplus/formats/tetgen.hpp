#pragma once
//
// TetGen (.node / .ele) reader/writer. TetGen stores a mesh as a pair of
// sibling files sharing a stem: <stem>.node (points + attributes + boundary
// markers) and <stem>.ele (tetrahedra + region attributes). Either path
// selects the pair.

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_tetgen(const std::string& path, const Mesh& mesh);
Mesh read_tetgen(const std::string& path);

}  // namespace meshioplusplus
