#pragma once
//
// FLAC3D (.f3grid) reader/writer — common path. Handles the ASCII and binary
// gridpoint/zone/face sections with the meshio<->FLAC3D node orderings (incl.
// the determinant-based right-handed zone reorder). Cell groups (ZGROUP/
// FGROUP -> cell_sets) are left to the Python fallback: the reader throws when
// it sees a group and the shim routes meshes carrying cell_sets to Python.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_flac3d(const std::string& path, const Mesh& mesh,
                  const std::string& float_fmt, bool binary);
Mesh read_flac3d(const std::string& path);

}  // namespace meshio
