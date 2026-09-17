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
 * @file formats/vtr.hpp
 * @brief VTK XML RectilinearGrid (`.vtr`): a lattice whose per-axis point
 * coordinates are three 1-D arrays rather than a uniform `Origin`/`Spacing`
 * pair (v11.6.0, roadmap §1 tier B4).
 *
 * A RectilinearGrid states the same `nx * ny * nz` hexahedron topology
 * `.vti`/`.vts` do -- points ordered x-fastest, `detail/grid_lattice.hpp`'s
 * own numbering -- but its `<Coordinates>` are three independent, only
 * *monotonic* 1-D arrays (`x_coordinates`/`y_coordinates`/`z_coordinates`),
 * evaluated as a tensor product: point `(i, j, k)` sits at
 * `(xs[i], ys[j], zs[k])`. That is genuinely more general than `.vti`'s
 * uniform spacing -- a graded mesh (finer near a wall, coarser far from it)
 * is a RectilinearGrid, never an ImageData.
 *
 * ### The mesh side of the deal
 *
 * - **`read_vtr` is fully general**: it builds points from the tensor
 *   product of the file's own three coordinate arrays, with no uniformity
 *   check at all -- a genuinely graded grid reads correctly.
 * - **`write_vtr` requires a *uniform* lattice**, exactly as `.vti`'s writer
 *   does (via `detail::lattice_from_mesh`): recovering three arbitrary
 *   per-axis coordinate arrays from an unstructured point set, rather than
 *   one `Origin`/`Spacing` pair, needs the same lattice detection this
 *   writer does not re-derive. **A genuinely non-uniform (graded) mesh
 *   cannot be written as `.vtr` today** -- a documented follow-up, not a
 *   silent gap: `lattice_from_mesh` is the single owner of "is this mesh a
 *   dense lattice" and extending it to recover ungraded per-axis arrays is
 *   future work, tracked in `doc/roadmap.md`.
 *
 * Data arrays reuse the same `detail/vtk_xml.hpp`/`detail/vtu_binary.hpp`
 * codec machinery `.vti`/`.vts`/`.vtu` already use.
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
 * @brief Write a mesh as VTK XML RectilinearGrid.
 * @param rPath the output path.
 * @param rMesh the mesh; must be a dense, UNIFORM lattice (see the file docs
 *        -- a graded mesh cannot be written today).
 * @param binary base64-encode the arrays instead of writing them as text.
 * @param zlib compress the binary blocks. Ignored when @p binary is false.
 * @throws WriteError when @p rMesh is not a dense uniform lattice, or when
 *         zlib was requested and this build has none.
 */
MESHIOPLUSPLUS_API void write_vtr(const std::string& rPath, const Mesh& rMesh, bool binary = true,
                                  bool zlib = true);

/**
 * @brief Write a mesh as VTK XML RectilinearGrid with an explicit block codec.
 * @param rPath the output path.
 * @param rMesh the mesh; must be a dense, uniform lattice.
 * @param binary base64-encode the arrays instead of writing them as text.
 * @param codec the block compressor; `None` writes uncompressed base64.
 * @throws WriteError as `write_vtr`, and when @p codec is not in this build.
 */
MESHIOPLUSPLUS_API void write_vtr_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                                        detail::VtkCodec codec);

/**
 * @brief Read a VTK XML RectilinearGrid file.
 * @param rPath the input path.
 * @param rOpts selective-read options; `mPointsOnly` and `mDataArrays` apply.
 * @return a mesh with one `hexahedron` block in `detail/grid_lattice.hpp`'s
 *         index-formula numbering, or a point-only mesh when the extent has
 *         no cells. Points are the tensor product of the file's own
 *         per-axis coordinate arrays -- no uniformity is assumed or checked.
 * @throws ReadError on a construct the C++ reader declines (see the file docs).
 */
MESHIOPLUSPLUS_API Mesh read_vtr(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief Summarize a VTK XML RectilinearGrid file without decoding its arrays.
 *
 * `WholeExtent` gives both the point and the cell count without decoding
 * `<Coordinates>` or any data array.
 */
MESHIOPLUSPLUS_API MeshMetadata read_vtr_metadata(const std::string& rPath,
                                                  const ReadOptions& rOpts = {});

}  // namespace meshioplusplus
