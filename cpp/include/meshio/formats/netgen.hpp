#pragma once
//
// Netgen neutral mesh (.vol) reader/writer — common path only. Handles the
// mesh3d/dimension/geomtype/points/{point,edge,surface,volume}elements
// sections with the netgen<->meshio node permutations and the single integer
// cell index (`netgen:index`). Advanced content (identifications, codim
// material/bc names, the two-line edgesegmentsgi2 variant, and the gzip
// `.vol.gz` container) is left to the Python fallback: the reader throws on
// those tokens and the writer is gated off by the shim.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_netgen(const std::string& path, const Mesh& mesh, const std::string& float_fmt);
Mesh read_netgen(const std::string& path);

}  // namespace meshio
