#pragma once
//
// Nastran bulk-data (.bdf/.fem/.nas) writer + reader. The writer emits meshio's
// default layout (fixed-large GRID* points, fixed-small element cards) plus a
// sentinel comment; the reader only parses files carrying that sentinel.
// Arbitrary Nastran files (reference meshes, buffers) use the Python fallback.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_nastran(const std::string& path, const Mesh& mesh);
Mesh read_nastran(const std::string& path);

}  // namespace meshio
