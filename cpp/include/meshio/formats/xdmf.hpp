#pragma once
//
// XDMF3 reader/writer — XML and Binary data paths only (no HDF5 linkage).
// DataItem Format="XML" (inline) and Format="Binary" (raw external .bin) are
// handled in C++ via pugixml; Format="HDF" is left to the Python/h5py
// fallback (the reader throws on it and the shim routes data_format="HDF" to
// Python). Supports single-type and Mixed topology, Node/Cell attributes, and
// XY/XYZ geometry.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_xdmf(const std::string& path, const Mesh& mesh, const std::string& data_format);
Mesh read_xdmf(const std::string& path);

}  // namespace meshio
