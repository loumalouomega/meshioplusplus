#pragma once
//
// Abaqus .inp reader/writer (ascii). Covers *NODE and *ELEMENT (with the abaqus
// element-type names). Files using *NSET / *ELSET / *INCLUDE, and meshes
// carrying point_sets / cell_sets, are handled by the Python fallback.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_abaqus(const std::string& path, const Mesh& mesh);
Mesh read_abaqus(const std::string& path);

}  // namespace meshio
