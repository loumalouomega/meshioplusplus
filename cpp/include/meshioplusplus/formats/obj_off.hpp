#pragma once
//
// Wavefront OBJ and Geomview OFF readers/writers (ascii text formats).

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_off(const std::string& path, const Mesh& mesh);
Mesh read_off(const std::string& path);

void write_obj(const std::string& path, const Mesh& mesh);
Mesh read_obj(const std::string& path);

}  // namespace meshioplusplus
