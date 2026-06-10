#pragma once
//
// Tecplot ASCII finite-element reader/writer. Handles a single FE zone
// (FETRIANGLE / FEQUADRILATERAL / FETETRAHEDRON / FEBRICK), BLOCK and POINT
// data packing, point_data and cell-centered cell_data (VARLOCATION). Meshes
// with multiple cell types fall back to the Python implementation.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_tecplot(const std::string& path, const Mesh& mesh);
Mesh read_tecplot(const std::string& path);

}  // namespace meshio
