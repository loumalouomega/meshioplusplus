#pragma once
//
// SU2 (.su2) ascii mesh reader/writer. Handles NDIME / NPOIN / NELEM (volume
// cells) and NMARK boundary markers (-> su2:tag cell_data).

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_su2(const std::string& path, const Mesh& mesh);
Mesh read_su2(const std::string& path);

}  // namespace meshioplusplus
