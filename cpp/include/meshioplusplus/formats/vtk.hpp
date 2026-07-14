#pragma once
//
// VTK legacy (.vtk) reader/writer. Writer: format version 5.1,
// UNSTRUCTURED_GRID, ASCII and big-endian binary, for non-polyhedron meshes.
// Reader: version 5.1 UNSTRUCTURED_GRID, ASCII and big-endian binary. The
// Python implementation remains the fallback (version 4.2, structured grids,
// 2-component vector padding, polyhedron).

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

// Write `mesh` to `path` as a VTK legacy UNSTRUCTURED_GRID file (ascii or
// big-endian binary). `v51` selects version 5.1 (OFFSETS/CONNECTIVITY) vs 4.2
// (interleaved CELLS). Throws meshioplusplus::WriteError on unsupported input.
void write_vtk(const std::string& path, const Mesh& mesh, bool binary, bool v51);

// Read a VTK legacy UNSTRUCTURED_GRID file (versions 4.2 and 5.1, ascii or
// big-endian binary). Throws meshioplusplus::ReadError on anything it doesn't handle so
// the Python reader can take over.
Mesh read_vtk(const std::string& path);

}  // namespace meshioplusplus
