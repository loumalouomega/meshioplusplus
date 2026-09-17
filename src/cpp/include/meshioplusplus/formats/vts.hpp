//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file formats/vts.hpp
 * @brief VTK XML StructuredGrid (`.vts`): a lattice with explicit points but
 * implicit (index-formula) connectivity (v11.6.0, roadmap §1 tier B4).
 *
 * `.vti` (ImageData) states a lattice's geometry as three attributes
 * (`Origin`/`Spacing`/`WholeExtent`); `.vts` states the SAME topology --
 * `nx * ny * nz` hexahedra over a `WholeExtent` corner grid, points ordered
 * x-fastest (`detail/grid_lattice.hpp`'s own numbering) -- but with an
 * explicit `<Points>` array instead, exactly as `.vtu` writes one. That is
 * the entire difference: connectivity is still implicit (the index formula,
 * never written), so a StructuredGrid is a lattice whose points are allowed
 * to be curved or non-uniformly spaced while ImageData's cannot.
 *
 * ### The mesh side of the deal
 *
 * - **`read_vts` expands** `WholeExtent` into hexahedron connectivity via the
 *   same index formula `detail/grid_lattice.hpp` uses, but reads the point
 *   *positions* verbatim from the file's own `<Points>` array rather than
 *   recomputing them -- unlike `.vti`, a `.vts` file's points are not
 *   required to be an even grid at all (VTK itself allows a curved
 *   structured mesh; this reader accepts that geometry, it just never
 *   needs to verify it, since implicit connectivity does not depend on it).
 * - **`write_vts` requires a lattice**, exactly as `.vti` does (via
 *   `detail::lattice_from_mesh`): the writer has no way to recover which
 *   `(nx, ny, nz)` a mesh's points were meant to tile without one. Once
 *   confirmed, the mesh's own points are written unchanged -- there is
 *   nothing to recompute, unlike `.vti`'s Origin/Spacing attributes.
 *
 * Data arrays reuse the same `detail/vtk_xml.hpp`/`detail/vtu_binary.hpp`
 * codec machinery `.vti`/`.vtu` already use.
 *
 * ### Deliberately not supported (both raise, so a shim falls back to Python)
 *
 * Identical to `.vti`'s list: `<AppendedData>`, more than one `<Piece>` or a
 * piece whose `Extent` is not the `WholeExtent`, lzma and any codec this
 * build lacks. `header_type="UInt64"` is honoured on read; the writer always
 * emits the default `UInt32`.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as VTK XML StructuredGrid.
 * @param rPath the output path.
 * @param rMesh the mesh; must be a dense lattice (see the file docs).
 * @param binary base64-encode the arrays instead of writing them as text.
 * @param zlib compress the binary blocks. Ignored when @p binary is false.
 * @throws WriteError when @p rMesh is not a dense lattice, or when zlib was
 *         requested and this build has none.
 */
MESHIOPLUSPLUS_API void write_vts(const std::string& rPath, const Mesh& rMesh, bool binary = true,
                                  bool zlib = true);

/**
 * @brief Write a mesh as VTK XML StructuredGrid with an explicit block codec.
 * @param rPath the output path.
 * @param rMesh the mesh; must be a dense lattice.
 * @param binary base64-encode the arrays instead of writing them as text.
 * @param codec the block compressor; `None` writes uncompressed base64.
 * @throws WriteError as `write_vts`, and when @p codec is not in this build.
 */
MESHIOPLUSPLUS_API void write_vts_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                                        detail::VtkCodec codec);

/**
 * @brief Read a VTK XML StructuredGrid file.
 * @param rPath the input path.
 * @param rOpts selective-read options; `mPointsOnly` and `mDataArrays` apply.
 * @return a mesh with one `hexahedron` block in `detail/grid_lattice.hpp`'s
 *         index-formula numbering, or a point-only mesh when the extent has
 *         no cells.
 * @throws ReadError on a construct the C++ reader declines (see the file docs).
 */
MESHIOPLUSPLUS_API Mesh read_vts(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief Summarize a VTK XML StructuredGrid file without decoding its arrays.
 *
 * `WholeExtent` gives both the point and the cell count without decoding
 * `<Points>` or any data array; unlike `.vti`, the bounding box is NOT free
 * (points are explicit and would have to be decoded), so it is not reported.
 */
MESHIOPLUSPLUS_API MeshMetadata read_vts_metadata(const std::string& rPath,
                                                  const ReadOptions& rOpts = {});

}  // namespace meshioplusplus
