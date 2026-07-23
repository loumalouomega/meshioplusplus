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
 * @file detail/cell_index.hpp
 * @brief The single owner of the **global (block-major) cell index**: block 0's
 * cells first, then block 1's, and so on.
 *
 * Several parts of meshio++ already number cells this way — `partition_labels`
 * fills its flat buffer in this order, `split`'s connected-component pass
 * union-finds over it, and `Region`s of kind `Cell` and `Side` index into it.
 * Each of those used to re-derive the prefix sum locally. They no longer do:
 * the browser viewer's `_dimension_rows` counter bug (a running index that
 * failed to advance over *skipped* blocks, so every per-cell quantity silently
 * shifted) is the cautionary tale for what a second, subtly different copy of
 * this arithmetic costs.
 *
 * Header-only inline helpers — they are called once per operation, not per
 * element, but they are small enough that a translation unit of their own would
 * only add link noise. Built on the uniform mesh API, so they compile under
 * every mesh backend.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Prefix sums of the mesh's per-block cell counts.
 *
 * Entry `b` is the global index of block `b`'s first cell; the last entry is
 * the mesh's total cell count. Every block contributes, including empty ones
 * and ragged ones, so the numbering never skips.
 * @param rMesh The mesh to measure.
 * @return `NumCellBlocks() + 1` ascending offsets, starting at 0.
 */
inline std::vector<std::int64_t> block_bases(const Mesh& rMesh) {
    std::vector<std::int64_t> bases;
    bases.reserve(rMesh.NumCellBlocks() + 1);
    bases.push_back(0);
    for (const auto cb : rMesh.CellRange())
        bases.push_back(bases.back() + static_cast<std::int64_t>(cb.NumCells()));
    return bases;
}

/// Total number of cells described by a `block_bases` table.
inline std::int64_t total_cells(const std::vector<std::int64_t>& rBases) {
    return rBases.empty() ? 0 : rBases.back();
}

/**
 * @brief Global cell index of cell @p row of block @p block.
 * @param rBases A `block_bases` table.
 * @param block Block index.
 * @param row Cell index within that block.
 * @return The global (block-major) cell index.
 */
inline std::int64_t block_row_to_global(const std::vector<std::int64_t>& rBases, std::size_t block,
                                        std::int64_t row) {
    return rBases[block] + row;
}

/**
 * @brief Split a global cell index back into `(block, row)`.
 *
 * Uses a linear scan rather than a binary search: meshes have a handful of
 * blocks, and callers that convert many indices should hoist the scan (or sort
 * their indices, which regions already are) instead.
 * @param rBases A `block_bases` table.
 * @param global The global cell index; must be in `[0, total_cells)`.
 * @return `(block, row)`, or `(npos, 0)` when @p global is out of range.
 */
inline std::pair<std::size_t, std::int64_t> global_to_block_row(
    const std::vector<std::int64_t>& rBases, std::int64_t global) {
    constexpr std::size_t npos = static_cast<std::size_t>(-1);
    if (rBases.size() < 2 || global < 0 || global >= rBases.back())
        return {npos, 0};
    for (std::size_t b = 0; b + 1 < rBases.size(); ++b)
        if (global < rBases[b + 1])
            return {b, global - rBases[b]};
    return {npos, 0};
}

}  // namespace detail
}  // namespace meshioplusplus
