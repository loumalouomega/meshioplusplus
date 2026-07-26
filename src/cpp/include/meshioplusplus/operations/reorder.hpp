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
 * @file operations/reorder.hpp
 * @brief Dependency-free mesh renumbering: permute node/element indices to
 * reduce sparse-matrix bandwidth and improve cache locality.
 *
 * `reorder` computes a node permutation with one of three strategies and
 * applies it as a *pure permutation* — no nodes or cells are added or removed,
 * geometry and all data are preserved, and the applied permutations (old->new)
 * are returned so callers can remap external arrays:
 *
 *  - **RCM** (Reverse Cuthill-McKee): a bandwidth-reducing ordering derived
 *    from the node-adjacency graph (nodes sharing a cell are neighbours),
 *    started from a pseudo-peripheral node per connected component.
 *  - **Morton** (Z-order) and **Hilbert**: space-filling-curve orderings over
 *    the node coordinates (bounding box quantized to 21 bits/axis), giving good
 *    locality without needing connectivity.
 *
 * A derived element permutation orders each cell block by the minimum new node
 * index of its cells, and `point_data`/`cell_data` follow their nodes/cells.
 * Everything uses only standard C++ and the uniform mesh API, so it runs under
 * every mesh backend. This is an operation, not a file format — it is not in
 * the format registry.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Renumbering strategy selected by `reorder`.
enum class ReorderMethod {
    RCM,      ///< Reverse Cuthill-McKee on the node-adjacency graph.
    Morton,   ///< Morton (Z-order) space-filling curve on node coordinates.
    Hilbert,  ///< Hilbert space-filling curve on node coordinates.
};

/**
 * @brief The result of `reorder`: the permuted mesh plus the applied
 * permutations, so external per-node / per-cell arrays can be remapped.
 */
struct ReorderResult {
    /// The renumbered mesh (points, connectivity, point/cell/field data).
    Mesh mMesh;
    /// Node permutation, `Int64` shape `(num_points,)`, old index -> new index
    /// (i.e. `new_points[mNodePermutation[i]] == old_points[i]`).
    NDArray mNodePermutation;
    /// Per cell block, an `Int64` shape `(num_cells_in_block,)` cell permutation
    /// old index -> new index, in cell-block order.
    std::vector<NDArray> mCellPermutations;
};

/**
 * @brief Parse a method name into a `ReorderMethod`.
 * @param rName one of `"rcm"`, `"morton"`, `"hilbert"` (case-insensitive for
 *              `rcm`; `"z"`/`"zorder"` also map to Morton)
 * @return the matching enum value
 * @throws std::invalid_argument on an unknown name
 */
MESHIOPLUSPLUS_API ReorderMethod reorder_method_from_name(const std::string& rName);

/**
 * @brief Connectivity bandwidth of a mesh: the maximum over all cells of
 * `(max node index - min node index)` within a cell.
 *
 * A smaller value means node indices referenced together sit closer, which is
 * what RCM minimises. Useful as a before/after measure of a reordering.
 * @param rMesh the mesh to measure
 * @return the bandwidth (0 for an empty / node-less mesh)
 */
MESHIOPLUSPLUS_API std::int64_t compute_bandwidth(const Mesh& rMesh);

/**
 * @brief Renumber a mesh with the given strategy (pure permutation).
 * @param rMesh the mesh to renumber
 * @param method the renumbering strategy (default RCM)
 * @return a `ReorderResult` carrying the permuted mesh and the applied
 *         node/cell permutations
 */
MESHIOPLUSPLUS_API ReorderResult reorder(const Mesh& rMesh, ReorderMethod method = ReorderMethod::RCM);

}  // namespace meshioplusplus
