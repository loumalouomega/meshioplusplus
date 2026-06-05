#pragma once
//
// STL reader/writer (ascii and binary). STL stores raw triangles; on read the
// vertices are de-duplicated in first-occurrence order (matching the net result
// of meshio's numpy implementation).

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_stl(const std::string& path, const Mesh& mesh, bool binary);
Mesh read_stl(const std::string& path);

}  // namespace meshio
