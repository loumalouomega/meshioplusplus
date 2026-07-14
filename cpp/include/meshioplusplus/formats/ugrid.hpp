#pragma once
//
// AFLR UGRID (.ugrid) reader/writer. Handles the ascii variant and every
// binary variant encoded in the file's penultimate suffix (b8l/b8/b4/lb8l/
// lb8/lb4 C-type, r8/r4/lr8/lr4 Fortran-record type) — i.e. {big,little}
// endian, {4,8}-byte ints, {4,8}-byte floats, with optional Fortran record
// length markers.

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_ugrid(const std::string& path, const Mesh& mesh);
Mesh read_ugrid(const std::string& path);

}  // namespace meshioplusplus
