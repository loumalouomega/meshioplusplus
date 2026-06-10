#pragma once
//
// AVS-UCD (.avs) ascii reader/writer. Handles nodes, cells (with per-cell
// material -> avsucd:material), node data and cell data sections, and the
// avsucd<->meshio node orderings.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_avsucd(const std::string& path, const Mesh& mesh);
Mesh read_avsucd(const std::string& path);

}  // namespace meshio
