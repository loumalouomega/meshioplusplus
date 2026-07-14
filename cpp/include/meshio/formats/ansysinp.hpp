#pragma once
//
// Ansys MAPDL "coded database" (.cdb / .inp) reader/writer, mirroring
// src/meshio/ansysInp/_ansysInp.py: ET/ETBLOCK element-type declarations,
// NBLOCK nodes, EBLOCK elements (with continuation lines for >8-node cells),
// and CMBLOCK named components.
//
// CMBLOCK components map to point_sets / cell_sets, which are custom
// attributes on the Python Mesh (not carried by the Mesh conversion layer), so
// they travel through the AnsysInfo side-channel (the MedInfo pattern).

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "meshio/mesh.hpp"

namespace meshio {

struct AnsysInfo {
    // Component name -> node indices (0-based), from `CMBLOCK ...,NODE`.
    std::map<std::string, std::vector<std::int64_t>> point_sets;
    // Component name -> per-cell-block lists of local cell indices (0-based),
    // one inner list per mesh cell block in block order, from
    // `CMBLOCK ...,ELEM`.
    std::map<std::string, std::vector<std::vector<std::int64_t>>> cell_sets;
};

Mesh read_ansysinp(const std::string& path, AnsysInfo& info);
void write_ansysinp(const std::string& path, const Mesh& mesh, const AnsysInfo& info);

}  // namespace meshio
