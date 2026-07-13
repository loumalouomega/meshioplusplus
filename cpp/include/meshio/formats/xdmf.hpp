#pragma once
//
// XDMF3 reader/writer. DataItem Format="XML" (inline), Format="Binary" (raw
// external .bin) and — when built with MESHIO_HAS_HDF5 — Format="HDF"
// (sibling .h5 file, optional gzip) are handled in C++; without HDF5 the HDF
// paths throw and the Python/h5py fallback takes over. Supports single-type
// and Mixed topology, Node/Cell attributes, and XY/XYZ geometry.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

void write_xdmf(const std::string& path, const Mesh& mesh,
                const std::string& data_format, int gzip_level = -1);
Mesh read_xdmf(const std::string& path);

}  // namespace meshio
