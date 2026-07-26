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
 * @file detail/subset.hpp
 * @brief Shared "build a submesh from a chosen set of cells" helper used by the
 * `crop` and `split` operations: prune the points to only those referenced by
 * the kept cells, compact + remap connectivity, and gather all data arrays.
 *
 * Free functions in `meshioplusplus::detail`, each called once per crop/split
 * operation (not per element), so their bodies live in
 * `src/cpp/src/detail/subset.cpp` rather than inline here. Built on the uniform
 * mesh API only, so it works under every backend.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// The result of `build_cell_subset`: the pruned submesh plus the index maps.
struct SubsetResult {
    Mesh mMesh;
    /// Int64 shape `(num_points_in,)`, input point index -> output index (-1 if pruned).
    NDArray mPointMap;
    /// Per input block, Int64 shape `(num_cells_in_block,)`, input cell -> output
    /// index within that block (-1 if not kept).
    std::vector<NDArray> mCellMaps;
};

/// Row-preserving gather: out[j] = src[idx[j]] (whole trailing dims), dtype kept.
MESHIOPLUSPLUS_API NDArray subset_gather_rows(const NDArray& rSrc, const std::vector<std::int64_t>& rIdx);

/**
 * @brief Build a submesh containing exactly the chosen cells, pruning points to
 * only those referenced and remapping connectivity + all data.
 * @param rMesh the source mesh.
 * @param rKeptCellsPerBlock per input block, the (ascending) old local cell
 *        indices to keep.
 * @param rPointIdName if non-empty, attach an Int64 `point_data` of original
 *        point indices under this name.
 * @param rCellIdName if non-empty, attach an Int64 `cell_data` of original cell
 *        indices under this name (per block).
 * @param drop_empty_blocks omit blocks that keep no cells from the output (the
 *        `split` behaviour); the input -> output block shift this causes is
 *        handled for the named regions below.
 * @param rOpName the calling operation's name, used only in the warning
 *        `detail::remap_regions` emits when a region cannot be carried.
 * @return the pruned submesh — with its `point_data`, `cell_data`, `field_data`
 *         **and named regions** already carried over — plus the index maps.
 */
MESHIOPLUSPLUS_API SubsetResult build_cell_subset(const Mesh& rMesh,
                               const std::vector<std::vector<std::int64_t>>& rKeptCellsPerBlock,
                               const std::string& rPointIdName = "",
                               const std::string& rCellIdName = "",
                               bool drop_empty_blocks = false,
                               const std::string& rOpName = "");

}  // namespace detail
}  // namespace meshioplusplus
