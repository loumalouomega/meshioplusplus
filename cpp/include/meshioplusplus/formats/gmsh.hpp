#pragma once
//
// Gmsh .msh reader/writer. This implementation covers format version 2.2
// (ascii and binary): nodes, elements with physical/geometrical tags,
// $PhysicalNames (field_data), and $NodeData / $ElementData. Periodic meshes
// and format versions 4.0 / 4.1 are handled by the Python fallback.

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

// Write `mesh` to `path` as a Gmsh 2.2 .msh file (ascii or binary).
void write_gmsh22(const std::string& path, const Mesh& mesh, bool binary);

// Write `mesh` to `path` as a Gmsh 4.1 .msh file (ascii or binary). Intended for
// meshes without entity information (no gmsh:dim_tags); $Entities is not emitted.
void write_gmsh41(const std::string& path, const Mesh& mesh, bool binary);

// Read a Gmsh .msh file. Throws meshioplusplus::ReadError for anything not handled
// (version != 2.2, $Periodic) so the Python reader can take over.
Mesh read_gmsh(const std::string& path);

}  // namespace meshioplusplus
