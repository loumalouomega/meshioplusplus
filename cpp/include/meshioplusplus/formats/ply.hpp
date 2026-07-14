#pragma once
//
// PLY (Polygon File Format) reader/writer, ascii and binary (little/big
// endian). Vertex properties beyond x/y/z become point_data; faces (a list
// property) become cell blocks grouped by vertex count. Faces with extra
// (non-index) properties fall back to the Python implementation.

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_ply(const std::string& path, const Mesh& mesh, bool binary);
Mesh read_ply(const std::string& path);

}  // namespace meshioplusplus
