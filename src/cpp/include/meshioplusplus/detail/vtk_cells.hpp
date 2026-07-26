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
 * @file vtk_cells.hpp
 * @brief Shared reconstruction of meshio cell blocks from the VTK/VTU
 * connectivity + end-offsets + types representation.
 *
 * Both the VTU reader and the VTK 5.1 legacy reader store cells in the same
 * layout — a flat `connectivity` array of node indices, an `offsets` array
 * giving each cell's end position within it, and a `types` array giving each
 * cell's VTK type id — so this header's `detail::reconstruct_cells` (ported
 * from `vtk_cells_from_data` in `_vtk_common.py`) is the single place that
 * turns that layout back into meshio's per-type cell-block representation
 * (appended straight onto the output `Mesh`),
 * grouping consecutive same-type runs and further splitting runs of
 * variable-node-count types (polygon, VTK_LAGRANGE_*) by per-cell size.
 * It leans heavily on `parallel_for_bw`/`parallel_copy_i64` (memory-gather
 * and memory-fault-bound work) since reconstructing connectivity is pure
 * data movement, not compute.
 *
 * Every function here is called once per file read (or once per contiguous
 * run within it), never once per element, so bodies live in
 * `src/cpp/src/detail/vtk_cells.cpp` rather than inline here — moving them does
 * not change the granularity of their own internal `parallel_for_bw` loops,
 * which are unaffected by which translation unit compiles the enclosing
 * function.
 */

// System includes
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Copies `n` `int64_t` elements from `src` to `dst`, splitting the
 * copy into large contiguous chunks run across `parallel_for_bw`'s
 * bandwidth-capped threads.
 *
 * `dst` is assumed to be a fresh allocation, so most of the wall-clock cost
 * is first-touch page faults (the OS zeroing/mapping pages on first write)
 * rather than the memcpy itself — servicing those faults concurrently across
 * a few threads beats one thread doing a single serial `memcpy`. Falls back
 * to a single sequential `memcpy` when `n` doesn't even fill one 4 MiB chunk
 * (`nchunks <= 1`). Uses `grain=1` so every chunk (already coarse at 512Ki
 * elements) dispatches individually rather than being batched further by the
 * default grain.
 * @param pDst Destination buffer, at least `n` elements, ideally freshly
 *            allocated (unfaulted) memory.
 * @param pSrc Source buffer, at least `n` elements.
 * @param n Number of `int64_t` elements to copy.
 */
MESHIOPLUSPLUS_API void parallel_copy_i64(std::int64_t* pDst, const std::int64_t* pSrc, std::size_t n);

/**
 * @brief Extracts rows `[r0, r1)` of a 2-D (or column-vector) `NDArray` into
 * a new, freshly-allocated `NDArray`.
 *
 * The output buffer is allocated via `NDArray::Uninit` (skipping the
 * zero-fill) since the single `memcpy` below fully overwrites it.
 * @param rA Source array; row size is `rA.Shape()[1]` if 2-D, else 1.
 * @param r0 First row to include (inclusive).
 * @param r1 One past the last row to include (exclusive).
 * @return A new owning `NDArray` with `r1 - r0` rows, same dtype/row-width as `rA`.
 */
MESHIOPLUSPLUS_API NDArray slice_rows(const NDArray& rA, std::size_t r0, std::size_t r1);

/**
 * @brief The cell-block *shape* `reconstruct_cells` would produce, without
 *        touching connectivity.
 *
 * Backs the VTU/VTP metadata path: `types` is one byte per cell and `offsets`
 * is only consulted for the variable-node-count types, so a summary costs a
 * fraction of a full reconstruction.
 *
 * The run-grouping here **duplicates** `reconstruct_cells`' -- consecutive
 * same-type runs, further split by per-cell size for `is_special_cell` types.
 * That duplication is deliberate (the two have very different inner loops and
 * merging them would slow the hot one down) but it can drift, so it is pinned
 * by a test asserting this function agrees with `metadata_from_mesh` applied to
 * a real reconstruction, rather than by comment alone.
 *
 * @param rOffsets End offsets, one per cell; may be **empty** when no special
 *                 cell type is present, since only those consult it.
 * @param rTypes VTK cell type id for each cell.
 * @return One entry per output block, in the same order as the blocks
 *         `reconstruct_cells` would append.
 * @throws ReadError on the same unsupported types `reconstruct_cells` rejects,
 *         so a summary never claims a file is readable when it is not.
 */
MESHIOPLUSPLUS_API std::vector<CellBlockInfo> summarize_cells(const std::vector<std::int64_t>& rOffsets,
                                           const std::vector<std::int64_t>& rTypes);

/** @brief Whether any type in @p rTypes needs `offsets` to be summarized. */
MESHIOPLUSPLUS_API bool cells_need_offsets(const std::vector<std::int64_t>& rTypes);

/**
 * @brief Reconstructs meshio cell blocks (and the matching per-block
 * `cell_data`) from the VTK/VTU flat connectivity + end-offsets + types
 * representation, appending them to @p rMesh.
 *
 * Ported from `vtk_cells_from_data` in `_vtk_common.py`; shared by the VTU
 * reader and the VTK 5.1 legacy reader, which store cells identically.
 * Walks `types` and groups consecutive cells of the same VTK type into a
 * run; a run of a *fixed*-node-count type becomes one rectangular cell
 * block (data gathered per-row via `vtk_to_meshio_order`, or
 * block-copied via `parallel_copy_i64` when the run is contiguous in
 * `conn` with no reordering needed); a run of a *variable*-node-count type
 * (`is_special_cell`, e.g. polygon or VTK_LAGRANGE_*) is further split into
 * sub-runs of a single common node count each, since a rectangular cell
 * block still requires a `(num_cells, n)` layout — each such sub-run
 * is emitted as its own separate block sharing the same meshio type
 * name. Matching slices of every array in `cell_data_raw` are appended to
 * @p rMesh's cell data in lockstep with the appended blocks, via
 * `slice_rows`.
 *
 * @param pConn Flat node-index connectivity buffer. Passed as a raw
 *             `int64_t*` (rather than an `NDArray`) so callers can hand in
 *             an `NDArray`'s buffer directly — VTK 5.1 connectivity is
 *             already `vtktypeint64` — without an intermediate
 *             to-int64 copy.
 * @param rOffsets End offsets, one per cell: `rOffsets[i]` is the index in
 *                `pConn` just past cell `i`'s last node (so cell `i`'s nodes
 *                are `pConn[rOffsets[i-1] .. rOffsets[i])`, with `rOffsets[-1]`
 *                treated as 0).
 * @param rTypes VTK cell type id for each cell, same length as `rOffsets`.
 * @param rCellDataRaw Per-name cell-data arrays covering the whole mesh
 *                      (all cells concatenated), to be re-sliced per output
 *                      block.
 * @param rMesh Mesh appended to: one rectangular cell block per contiguous
 *              same-type (and, for special types, same-size) run, plus — in
 *              lockstep — for each name in `rCellDataRaw` one sliced
 *              `NDArray` per new block.
 * @throws ReadError if a cell's VTK type id is 42 (polyhedron — unsupported
 *         by the C++ reader) or is otherwise not in `vtk_to_meshio_type()`,
 *         or if a resolved meshio type has no entry in `num_nodes_per_cell()`.
 */
MESHIOPLUSPLUS_API void reconstruct_cells(const std::int64_t* pConn, const std::vector<std::int64_t>& rOffsets,
                       const std::vector<std::int64_t>& rTypes,
                       const std::unordered_map<std::string, NDArray>& rCellDataRaw, Mesh& rMesh);

}  // namespace detail
}  // namespace meshioplusplus
