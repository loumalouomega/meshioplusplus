#pragma once
//
// OpenFOAM polyMesh reader (read-only), mirroring
// src/meshio/openfoam/_openfoam.py. A polyMesh is a directory of sibling files
// (points, faces, owner, neighbour, boundary) in ASCII or binary
// (little-endian, label=32/64, scalar=32/64). Cells are reconstructed from the
// face-based representation into tetra/pyramid/wedge/hexahedron (rectangular)
// or general polyhedra (ragged), with boundary faces grouped into
// triangle/quad/polygon blocks.
//
// Boundary patch names map to mesh.cell_tags (a custom Python Mesh attribute),
// carried through the OpenFoamInfo side-channel.

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "meshio/mesh.hpp"

namespace meshio {

struct OpenFoamInfo {
    // MED-style negative family id -> {patch name}.
    std::map<std::int64_t, std::vector<std::string>> cell_tags;
};

// `path` may be a `.foam` marker file, a case directory, or a polyMesh
// directory (resolved like the Python reader).
Mesh read_openfoam(const std::string& path, OpenFoamInfo& info);

}  // namespace meshio
