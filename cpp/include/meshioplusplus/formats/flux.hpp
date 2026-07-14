#pragma once
//
// FLUX .pf3 reader/writer. ASCII with French keyword headers; each element is a
// 12-int record (7th field = type descriptor, 8th = node count) followed by its
// 1-based connectivity, and nodes are listed under "COORDONNEES DES NOEUDS".
// Per-element region references ride in cell_data["pf3:ref"].

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_flux(const std::string& path, const Mesh& mesh);
Mesh read_flux(const std::string& path);

}  // namespace meshioplusplus
