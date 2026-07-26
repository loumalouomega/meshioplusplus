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
 * @file operations/merge.hpp
 * @brief Dependency-free mesh merge/weld: combine an ordered list of two or
 * more meshes into a single mesh.
 *
 * Two modes:
 *
 *  - **Concatenation** (default): every input's points are appended and each
 *    input's connectivity is offset by the running point count so indices stay
 *    valid; cell blocks of the same meshio type are merged (in input order).
 *    `point_data`/`cell_data` are combined per the `data_policy` (keys present
 *    in every input are concatenated; partial keys are either dropped —
 *    `Intersection` — or filled with NaN — `Fill`). A per-cell `source_mesh_id`
 *    tag records each cell's originating input (optional, on by default).
 *
 *  - **Welding** (`weld = true`): coincident points across all inputs within an
 *    absolute tolerance are detected with a spatial-hash bucket grid (never
 *    O(N^2)) and merged into one, remapping connectivity. When merged points
 *    carry `point_data`, the surviving point keeps the value of the *first*
 *    contributing point (earliest input / lowest global index) — a fixed,
 *    deterministic keep-first tie-break. `drop_duplicate_cells` additionally
 *    removes cells that become identical after welding.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format — it is not
 * in the format registry. `point_sets`/`cell_sets` are not reachable in the
 * core (they never enter the uniform mesh API); the Python shim namespaces and
 * remaps them using the `mPointMaps`/`mCellMaps` returned here.
 */

// System includes
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// How `merge` treats a `point_data`/`cell_data` key missing from some inputs.
enum class MergeDataPolicy {
    Intersection,  ///< Keep only keys present in every input (with matching cols); drop the rest.
    Fill,          ///< Keep every key; fill rows of inputs missing it with NaN (as Float64).
};

/// Options controlling `merge`.
struct MergeOptions {
    bool weld = false;       ///< Merge coincident points within `atol`.
    double atol = 1e-8;      ///< Absolute coincidence tolerance for weld.
    bool source_tag = true;  ///< Add Int64 cell_data "source_mesh_id".
    MergeDataPolicy data_policy = MergeDataPolicy::Intersection;  ///< Missing-key policy.
    bool drop_duplicate_cells = false;  ///< After welding, drop cells with identical connectivity.
};

/**
 * @brief The result of `merge`: the combined mesh plus the index maps needed to
 * remap external per-point / per-cell arrays (e.g. the Python-only sets).
 */
struct MergeResult {
    /// The merged mesh (points, connectivity, point/cell/field data).
    Mesh mMesh;
    /// Per input mesh, an `Int64` shape `(num_points_in_input,)` array mapping
    /// that input's local point index -> the output point index.
    std::vector<NDArray> mPointMaps;
    /// Per input mesh, an `Int64` shape `(num_cells_in_input,)` array mapping
    /// that input's global cell index (blocks concatenated in order) -> the
    /// output global cell index (or -1 if the cell was dropped as a duplicate).
    std::vector<NDArray> mCellMaps;
};

/**
 * @brief Merge two or more meshes into one.
 * @param rMeshes ordered, non-empty list of input meshes (pointers must be
 *                non-null; at least one mesh is required, two or more is the
 *                intended use)
 * @param rOpts merge options (concatenate vs weld, data policy, tagging)
 * @return a `MergeResult` carrying the merged mesh and the point/cell index maps
 * @throws std::invalid_argument if `rMeshes` is empty or contains a null pointer
 */
MESHIOPLUSPLUS_API MergeResult merge(const std::vector<const Mesh*>& rMeshes, const MergeOptions& rOpts = {});

}  // namespace meshioplusplus
