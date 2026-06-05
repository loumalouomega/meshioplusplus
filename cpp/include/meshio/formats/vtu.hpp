#pragma once
//
// VTU (VTK XML UnstructuredGrid) writer. Supports ASCII and binary (raw and
// zlib-compressed) output for non-polyhedron meshes (points, standard/polygon/
// Lagrange cells, point_data and cell_data). lzma compression and polyhedron
// cells are handled by the Python fallback.

#include <string>

#include "meshio/mesh.hpp"

namespace meshio {

// Write `mesh` to `path` as a .vtu file. `binary` selects base64 binary
// DataArrays; `zlib` additionally zlib-compresses them (ignored when ascii).
// Throws meshio::WriteError on unsupported input (e.g. polyhedron blocks).
void write_vtu(const std::string& path, const Mesh& mesh, bool binary, bool zlib);

// Read a .vtu file. Supports ascii and inline binary (uncompressed or zlib)
// DataArrays for non-polyhedron meshes. Throws meshio::ReadError on anything
// it doesn't handle (lzma, appended/raw data, polyhedron) so the Python reader
// can take over.
Mesh read_vtu(const std::string& path);

}  // namespace meshio
