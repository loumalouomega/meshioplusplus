#pragma once
//
// Modulef Formatted Mesh (.mfm) reader/writer — ASCII, single linear element
// type. Header of 8 ints, then vertex connectivity, reference arrays (emitted
// as zeros), vertex coordinates and a per-element subdomain array
// (cell_data["mfm:ref"]). Higher-order/multi-type meshes are rejected.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_mfm(const std::string& path, const Mesh& mesh, const std::string& float_fmt);
Mesh read_mfm(const std::string& path);

}  // namespace meshio
