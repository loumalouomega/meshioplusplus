#pragma once
//
// WKT (Well-Known Text) TIN reader/writer. Parses
// `TIN (((x y z, x y z, x y z, x y z)), ...)` triangulated irregular networks:
// points are de-duplicated in first-occurrence order and each triangle drops
// its repeated closing vertex.

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_wkt(const std::string& path, const Mesh& mesh);
Mesh read_wkt(const std::string& path);

}  // namespace meshioplusplus
