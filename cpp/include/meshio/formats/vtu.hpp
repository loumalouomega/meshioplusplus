#pragma once
//
// VTU (VTK XML UnstructuredGrid) writer. This first cut implements the ASCII
// format for non-polyhedron meshes (points, standard/polygon/Lagrange cells,
// point_data and cell_data). Binary, compression, and polyhedron support are
// handled by the Python fallback until ported.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

// Write `mesh` to `path` as an ASCII .vtu file. Throws meshio::WriteError on
// unsupported input (e.g. polyhedron blocks).
void write_vtu_ascii(const std::string& path, const Mesh& mesh);

}  // namespace meshio
