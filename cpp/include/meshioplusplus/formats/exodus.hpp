#pragma once
//
// Exodus II (.e/.exo/.ex2) reader/writer via netCDF, mirroring meshio's
// netCDF4-python dialect: dims num_nodes/num_dim/..., `coord` (transposed),
// `connect{k}` blocks with an elem_type attribute (1-based), point data via
// name_nod_var/vals_nod_var{k} (first time step) with the X/Y/Z–R/Z name
// recombination, and cell data via vals_elem_var{i}eb{k}. Node sets
// (point_sets) and info/qa records are not representable in the conversion
// layer: the reader throws on them and the shim defers to Python. Available
// only when built with MESHIOPLUSPLUS_HAS_NETCDF.

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

void write_exodus(const std::string& path, const Mesh& mesh);
Mesh read_exodus(const std::string& path);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
