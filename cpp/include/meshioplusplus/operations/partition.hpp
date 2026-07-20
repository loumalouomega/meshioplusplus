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
 * @file partition.hpp
 * @brief Mesh partitioning: decompose a mesh into `nparts` balanced pieces for
 * domain decomposition — the count-driven complement to the criterion-driven
 * `split`.
 *
 * Two methods:
 *
 *  - **SFC** (always available, the fallback): cells are ordered along a
 *    Hilbert space-filling curve of their centroids (bounding box quantized to
 *    21 bits per axis, the same key transforms `reorder` uses, shared via
 *    `detail/space_filling.hpp`) and the curve is cut into `nparts` contiguous
 *    ranges. Unweighted, the cell at curve rank `r` of `n` goes to part
 *    `(r * nparts) / n`, so part sizes differ by at most one cell. With
 *    `mWeightsKey` set (a scalar numeric `cell_data` array), the cut follows
 *    the weight prefix sum: a cell goes to the ideal interval containing its
 *    weight midpoint, `min(nparts - 1, floor((S_r + w_r / 2) * nparts / W))`.
 *    Both rules are monotone along the curve, hence contiguous; a part may be
 *    empty when a single weight exceeds `W / nparts`.
 *
 *  - **KaHIP** (optional, the quality path; `-DMESHIOPLUSPLUS_WITH_KAHIP=ON`,
 *    ported from the Kratos KaHIPApplication): the mesh's **dual graph** —
 *    cells are vertices, with an edge where two cells share a face (3D, via
 *    `detail/cell_faces.hpp`) or an edge (2D, via `detail/cell_edges.hpp`) —
 *    is handed in CSR form to KaHIP's serial `kaffpa()` interface with the
 *    requested `mImbalance`/`mMode`/`mSeed`. Only the sequential interface is
 *    used — never ParHIP — so no MPI dependency exists. Cells whose block has
 *    no facet table (ragged polygons, polyhedra, lower-dimensional cells in a
 *    mixed mesh, unsupported types) become isolated graph vertices: they are
 *    still assigned a part, just without adjacency preferences. When KaHIP is
 *    not compiled in, requesting it throws an error naming the CMake option
 *    (never a silent downgrade to SFC).
 *
 * **Determinism.** The SFC path is byte-identical across mesh backends and
 * thread counts: keys are computed in `parallel_for` into disjoint slots, the
 * argsort is a serial `std::stable_sort` tie-broken by cell index, and the cut
 * is a serial integer/prefix-sum rule. The KaHIP path is deterministic for a
 * fixed KaHIP build and seed, but its assignment may differ between KaHIP
 * versions — tests must assert balance and coverage, never exact labels.
 *
 * **Pieces keep the input block structure 1:1** (`drop_empty_blocks=false`,
 * deliberately unlike `split`): every piece has exactly `NumCellBlocks()`
 * blocks in input order, so `mCellMaps` indexes by input block and
 * concatenating the pieces reproduces the input mesh (partition of unity —
 * every cell lands in exactly one piece while `mGhostLayers == 0`).
 *
 * `mGhostLayers != 0` is reserved (shared-node BFS growth tagged
 * `partition:ghost`) and throws in v1. This is an operation, not a file
 * format — it is not in the format registry.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// The partitioning backend.
enum class PartitionMethod {
    SFC,    ///< Hilbert space-filling-curve cut (always available).
    KaHIP,  ///< KaHIP `kaffpa()` on the dual graph (optional dependency).
    Auto,   ///< KaHIP when compiled in, else SFC.
};

/// KaHIP preconfiguration (ignored by SFC). Per the Kratos KaHIPApplication
/// benchmarks, `Eco`/`Strong` carry the edge-cut wins over `Fast`.
enum class PartitionMode {
    Fast,
    Eco,
    Strong,
};

/// Parse `"sfc"` / `"kahip"` / `"auto"` (throws `std::invalid_argument`).
PartitionMethod partition_method_from_name(const std::string& rName);

/// Parse `"fast"` / `"eco"` / `"strong"` (throws `std::invalid_argument`).
PartitionMode partition_mode_from_name(const std::string& rName);

/// Options for `partition` / `partition_labels`.
struct PartitionOptions {
    /// Number of parts to decompose into (>= 1).
    int mNParts = 2;
    /// Partitioning backend; `Auto` resolves to KaHIP when compiled in.
    PartitionMethod mMethod = PartitionMethod::Auto;
    /// KaHIP only: allowed imbalance fraction, in (0, 1) (0.03 = 3%).
    double mImbalance = 0.03;
    /// KaHIP only: preconfiguration (default `Eco`).
    PartitionMode mMode = PartitionMode::Eco;
    /// KaHIP only: random seed (`kaffpa` is deterministic per seed).
    int mSeed = 0;
    /// Attach Int64 `partition:original_point_id` `point_data` and
    /// `partition:original_cell_id` `cell_data` (original input indices) to
    /// every piece. Only affects `partition`, not `partition_labels`.
    bool mRecordIds = false;
    /// Reserved: grow each piece by N shared-node BFS layers tagged
    /// `partition:ghost`. Not implemented in v1 — any value != 0 throws.
    int mGhostLayers = 0;
    /// Name of a scalar numeric `cell_data` array of per-cell weights
    /// (`""` = unweighted). SFC cuts the weight prefix sum; KaHIP receives it
    /// as vertex weights (scaled to positive integers).
    std::string mWeightsKey;
};

/// One piece of a partition.
struct PartitionPiece {
    /// The part id, in [0, nparts).
    int mPartId;
    /// The piece submesh (same number of blocks as the input, in input order).
    Mesh mMesh;
    /// Int64 shape `(num_points_in,)`, input point index -> piece point index
    /// (-1 if the point is not in this piece).
    NDArray mPointMap;
    /// Per input block, Int64 shape `(num_cells_in_block,)`, input cell ->
    /// piece cell index within that block (-1 if not in this piece).
    std::vector<NDArray> mCellMaps;
};

/// The result of `partition`: exactly `nparts` pieces, part id ascending
/// (pieces may be empty meshes).
struct PartitionResult {
    std::vector<PartitionPiece> mPieces;
};

/**
 * @brief Decompose a mesh into `nparts` balanced pieces.
 * @param rMesh The mesh to partition (unchanged).
 * @param rOptions Part count, method, and method parameters.
 * @return Exactly `nparts` pieces with their point/cell index maps.
 * @throws std::invalid_argument on `nparts < 1`, `ghost_layers != 0`, a bad
 *   weights array, an out-of-range imbalance, or `method == KaHIP` in a build
 *   without KaHIP (the error names `-DMESHIOPLUSPLUS_WITH_KAHIP=ON`).
 */
PartitionResult partition(const Mesh& rMesh, const PartitionOptions& rOptions = {});

/**
 * @brief Compute the per-cell part assignment without building the pieces.
 * @param rMesh The mesh to partition (unchanged).
 * @param rOptions Same options as `partition` (`mRecordIds`/`mGhostLayers`
 *   have no effect here).
 * @return Per input block, an Int64 array of shape `(num_cells_in_block,)`
 *   with values in [0, nparts) — block-aligned, ready to attach as the
 *   `partition:part` `cell_data`.
 * @throws std::invalid_argument as `partition`.
 */
std::vector<NDArray> partition_labels(const Mesh& rMesh, const PartitionOptions& rOptions = {});

/**
 * @brief Whether this build has the KaHIP backend compiled in.
 * @return `true` iff built with `MESHIOPLUSPLUS_WITH_KAHIP=ON` and KaHIP was
 *   found.
 */
bool partition_has_kahip() noexcept;

}  // namespace meshioplusplus
