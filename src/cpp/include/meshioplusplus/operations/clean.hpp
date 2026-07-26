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
 * @file operations/clean.hpp
 * @brief Dependency-free one-pass mesh cleanup: weld coincident points, drop
 * degenerate and exact-duplicate cells, and remove orphaned (unused) points.
 *
 * Each step is individually toggleable. `weld` fuses points within an absolute
 * tolerance using a spatial-hash bucket grid (never O(N^2)); the surviving point
 * keeps the value of the *first* contributing point (keep-first). Degenerate
 * cells (a repeated corner node after welding, or a near-zero measure) and
 * exact-duplicate cells (identical connectivity, keep-first per block) are
 * dropped. Orphaned points (referenced by no surviving cell) are removed and the
 * connectivity is remapped. The applied point/cell index maps are returned so a
 * caller can remap external per-point / per-cell arrays (the Python-only sets).
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles under
 * every mesh backend. This is an operation, not a file format — it is not in the
 * format registry.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Options controlling `clean` (each step is independently toggleable).
struct CleanOptions {
    bool weld = false;                 ///< Fuse coincident points within `atol`.
    double atol = 1e-8;                ///< Absolute coincidence tolerance for weld.
    bool remove_orphans = true;        ///< Drop points referenced by no surviving cell.
    bool drop_degenerate = true;       ///< Drop degenerate cells (repeated node / zero measure).
    bool drop_duplicate_cells = true;  ///< Drop exact-duplicate cells (keep-first per block).
};

/**
 * @brief The result of `clean`: the cleaned mesh, the index maps, and the counts
 * of what was welded / removed / dropped.
 */
struct CleanResult {
    /// The cleaned mesh.
    Mesh mMesh;
    /// `Int64` shape `(num_points_in,)` mapping each input point index to its
    /// output point index (or -1 if welded away / removed as an orphan).
    NDArray mPointMap;
    /// Per input cell block, an `Int64` shape `(num_cells_in_block,)` array
    /// mapping each input cell to its output index within that block (or -1 if
    /// the cell was dropped).
    std::vector<NDArray> mCellMaps;
    std::int64_t mPointsWelded = 0;            ///< Points fused away by welding.
    std::int64_t mPointsRemovedOrphan = 0;     ///< Points removed as orphans.
    std::int64_t mCellsDroppedDegenerate = 0;  ///< Cells dropped as degenerate.
    std::int64_t mCellsDroppedDuplicate = 0;   ///< Cells dropped as exact duplicates.
};

/**
 * @brief Clean a mesh in one pass (weld / drop-degenerate / drop-duplicate /
 * remove-orphans, per the options).
 * @param rMesh the mesh to clean (unmodified).
 * @param rOpts the cleanup options.
 * @return a `CleanResult` with the cleaned mesh, index maps, and removal counts.
 */
MESHIOPLUSPLUS_API CleanResult clean(const Mesh& rMesh, const CleanOptions& rOpts = {});

}  // namespace meshioplusplus
