#pragma once
//
// Ansys/Fluent .msh reader/writer. Hexadecimal, parenthesis-nested sections:
// (1 header), (2 dim), (10/2010/3010 nodes), (12/2012/3012 cells). ASCII and
// binary (float32/float64 nodes, int32/int64 cells). Face sections (13) with a
// data body are deferred to the Python fallback (the writer never emits them).

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_ansys(const std::string& path, const Mesh& mesh, bool binary);
Mesh read_ansys(const std::string& path);

}  // namespace meshioplusplus
