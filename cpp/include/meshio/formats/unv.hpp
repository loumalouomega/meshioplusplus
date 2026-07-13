#pragma once
//
// I-DEAS Universal (.unv) reader/writer — datasets 2411 (nodes) and 2412
// (elements) with the FE-descriptor -> meshio type map and the parabolic
// "sandwich" node reordering. Property ids ride in cell_data["unv:pid"].
// Permanent groups (2467) and fields need Python: the reader throws when it
// meets a 2467/2414 dataset and the shim routes meshes with point_sets/
// cell_sets to the Python writer.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_unv(const std::string& path, const Mesh& mesh);
Mesh read_unv(const std::string& path);

}  // namespace meshio
