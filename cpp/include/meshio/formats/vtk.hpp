#pragma once
//
// VTK legacy (.vtk) writer. First cut: ASCII, format version 5.1,
// UNSTRUCTURED_GRID, for non-polyhedron meshes. Binary, version 4.2, and the
// 2-component vector padding behaviour are handled by the Python fallback.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

// Write `mesh` to `path` as an ASCII VTK 5.1 file. Throws meshio::WriteError on
// unsupported input.
void write_vtk_ascii_51(const std::string& path, const Mesh& mesh);

}  // namespace meshio
