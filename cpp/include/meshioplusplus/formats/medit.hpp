#pragma once
//
// Medit / GMF (.mesh) ascii reader/writer. The binary .meshb variant (libMeshb
// GMF with position fields and version-typed records) is handled by the Python
// fallback.

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_medit_ascii(const std::string& path, const Mesh& mesh);
Mesh read_medit_ascii(const std::string& path);

}  // namespace meshioplusplus
