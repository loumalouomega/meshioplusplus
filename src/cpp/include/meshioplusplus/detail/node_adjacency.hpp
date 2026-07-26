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
 * @file node_adjacency.hpp
 * @brief The node-to-node graph of a mesh in compressed-sparse-row form,
 * shared by `operations/reorder.cpp` (bandwidth reduction) and
 * `operations/smooth.cpp` (the Laplacian/Taubin neighbour average).
 *
 * Two neighbour definitions are offered, because the two consumers genuinely
 * want different graphs:
 *
 * - `NodeAdjacencyKind::Clique` — every pair of nodes co-occurring in a cell is
 *   an edge. This is the sparsity pattern of an FEM stiffness matrix, which is
 *   exactly what a bandwidth-reducing renumbering must minimise, so it is what
 *   `reorder` uses.
 * - `NodeAdjacencyKind::Edge` — only nodes joined by an actual cell *edge*.
 *   A smoother must use this one: under the clique graph a hexahedron corner
 *   would be pulled by its 3 face diagonals and the body diagonal as strongly
 *   as by its 3 real edge neighbours, which collapses the element rather than
 *   regularising it.
 *
 * Edge topology comes from `detail/cell_subdivision.hpp`'s `cell_refine_edges`
 * (the same table the refine/elevate operations use, so the three cannot
 * drift). It covers all seven linear types — line, triangle, quad, tetra,
 * wedge, pyramid, hexahedron. Ragged polygon rows and polyhedron faces
 * contribute their natural ring edges, which is their true and unambiguous
 * edge topology.
 *
 * **Blocks whose edge topology is unknown — the whole higher-order family, the
 * VTK-Lagrange types, `custom` — contribute no edges at all** under
 * `Kind::Edge`; `node_edge_topology_known` reports them so a caller can pin
 * their nodes. This is deliberate, and not the same as falling back to the
 * clique: a `tetra10`'s mid-edge nodes have no meaningful clique centroid, so
 * moving them would silently distort a deliberately curved mesh, while moving
 * only the corners would leave the mid-nodes stranded off their edges. An
 * unknown neighbourhood means an undefined target, and an undefined target
 * means the node must hold still.
 *
 * Construction is the three-phase idiom shared with the operations layer:
 * a **serial** accumulation pass (so the graph is byte-identical regardless of
 * thread count), a `parallel_for` per-node sort/unique, and a prefix-sum plus
 * `parallel_for_bw` pack into CSR.
 *
 * Free functions in `meshioplusplus::detail`, each called once per
 * reorder/smooth operation (not per element), so their bodies live in
 * `src/cpp/src/detail/node_adjacency.cpp` rather than inline here — the same
 * convention `detail/subset.hpp` follows for `build_cell_subset`.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// Which nodes count as neighbours when building the graph.
enum class NodeAdjacencyKind {
    Clique,  ///< Every pair of nodes sharing a cell (FEM matrix sparsity).
    Edge     ///< Only nodes joined by a cell edge (what a smoother needs).
};

/**
 * @brief A node-to-node graph in compressed-sparse-row form.
 *
 * The neighbours of node `i` are `mAdj[mXadj[i] .. mXadj[i + 1])`, sorted
 * ascending and free of duplicates and self-loops.
 */
struct NodeAdjacency {
    std::vector<std::int64_t> mXadj;  ///< Row offsets, size `n + 1`.
    std::vector<std::int64_t> mAdj;   ///< Neighbour ids, size `mXadj[n]`.
};

/**
 * @brief Collects the node ids referenced by one cell, dropping out-of-range
 * entries so a malformed connectivity can never index a per-node array out of
 * bounds. Handles dense, ragged (polygon) and polyhedron blocks uniformly.
 * @param rBlock The cell block to read from.
 * @param Cell Index of the cell within @p rBlock.
 * @param NumPoints Node-id upper bound (ids outside `[0, NumPoints)` are dropped).
 * @param rOut Cleared and filled with the cell's node ids.
 */
MESHIOPLUSPLUS_API void cell_node_ids(const Mesh::CellView& rBlock, std::size_t Cell, std::size_t NumPoints,
                   std::vector<std::int64_t>& rOut);

/**
 * @brief Whether @p rBlock's cell edges can be enumerated exactly under
 * `NodeAdjacencyKind::Edge`.
 *
 * True for polyhedron blocks (the union of their face rings), ragged polygon
 * blocks (the row ring), and any type with a `cell_refine_edges` row (the seven
 * linear types). False for the higher-order family, the VTK-Lagrange types and
 * `custom` — those contribute no edges, and a smoother must pin their nodes
 * rather than guess a target for them.
 * @param rBlock The cell block to classify.
 * @return `true` when the block's edges are known exactly.
 */
MESHIOPLUSPLUS_API bool node_edge_topology_known(const Mesh::CellView& rBlock);

/**
 * @brief Builds the node-to-node graph of @p rMesh.
 * @param rMesh Mesh to read connectivity from (never modified).
 * @param NumPoints Number of nodes; the returned `mXadj` has `NumPoints + 1` entries.
 * @param Kind Which neighbour definition to use (see `NodeAdjacencyKind`).
 * @return The CSR graph, with each node's neighbour list sorted and deduplicated.
 */
MESHIOPLUSPLUS_API NodeAdjacency build_node_adjacency(const Mesh& rMesh, std::size_t NumPoints,
                                   NodeAdjacencyKind Kind);

}  // namespace detail
}  // namespace meshioplusplus
