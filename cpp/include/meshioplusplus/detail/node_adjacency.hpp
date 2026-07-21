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
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/value_io.hpp"

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
inline void cell_node_ids(const Mesh::CellView& rBlock, std::size_t Cell, std::size_t NumPoints,
                          std::vector<std::int64_t>& rOut) {
    rOut.clear();
    if (rBlock.IsPolyhedron()) {
        for (std::size_t f = 0; f < rBlock.NumFaces(Cell); ++f) {
            std::pair<const std::int64_t*, std::size_t> face = rBlock.Face(Cell, f);
            rOut.insert(rOut.end(), face.first, face.first + face.second);
        }
    } else if (rBlock.IsRagged()) {
        const std::int64_t* row = rBlock.Row(Cell);
        rOut.assign(row, row + rBlock.RowSize(Cell));
    } else {
        const NDArray& conn = rBlock.Conn();
        const std::size_t npc = rBlock.NodesPerCell();
        rOut.reserve(npc);
        for (std::size_t k = 0; k < npc; ++k)
            rOut.push_back(detail::read_int(conn, Cell * npc + k));
    }
    // Keep only in-range ids so out-of-bounds connectivity can never index a
    // permutation/adjacency array out of bounds.
    rOut.erase(std::remove_if(rOut.begin(), rOut.end(),
                              [NumPoints](std::int64_t v) {
                                  return v < 0 || static_cast<std::size_t>(v) >= NumPoints;
                              }),
               rOut.end());
}

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
inline bool node_edge_topology_known(const Mesh::CellView& rBlock) {
    if (rBlock.IsPolyhedron() || rBlock.IsRagged())
        return true;
    return !detail::cell_refine_edges(cell_type_from_name(rBlock.Type())).empty();
}

namespace node_adjacency_impl {

/// Adds every ordered pair of a cell's nodes (an element clique). Duplicates
/// and self-loops are removed later by the dedup pass.
inline void add_clique(std::vector<std::vector<std::int64_t>>& rAdj,
                       const std::vector<std::int64_t>& rNodes) {
    const std::size_t k = rNodes.size();
    for (std::size_t a = 0; a < k; ++a) {
        std::vector<std::int64_t>& row = rAdj[static_cast<std::size_t>(rNodes[a])];
        for (std::size_t b = 0; b < k; ++b)
            if (a != b)
                row.push_back(rNodes[b]);
    }
}

/// Adds the undirected edge `(u, v)` in both directions, skipping out-of-range
/// endpoints and self-loops.
inline void add_edge(std::vector<std::vector<std::int64_t>>& rAdj, std::size_t NumPoints,
                     std::int64_t u, std::int64_t v) {
    if (u < 0 || v < 0 || u == v)
        return;
    if (static_cast<std::size_t>(u) >= NumPoints || static_cast<std::size_t>(v) >= NumPoints)
        return;
    rAdj[static_cast<std::size_t>(u)].push_back(v);
    rAdj[static_cast<std::size_t>(v)].push_back(u);
}

/// Adds the closed ring of edges around an ordered node loop (a polygon row or
/// one face of a polyhedron).
inline void add_ring(std::vector<std::vector<std::int64_t>>& rAdj, std::size_t NumPoints,
                     const std::int64_t* pNodes, std::size_t Count) {
    if (Count < 2)
        return;
    for (std::size_t k = 0; k < Count; ++k)
        add_edge(rAdj, NumPoints, pNodes[k], pNodes[(k + 1) % Count]);
}

}  // namespace node_adjacency_impl

/**
 * @brief Builds the node-to-node graph of @p rMesh.
 * @param rMesh Mesh to read connectivity from (never modified).
 * @param NumPoints Number of nodes; the returned `mXadj` has `NumPoints + 1` entries.
 * @param Kind Which neighbour definition to use (see `NodeAdjacencyKind`).
 * @return The CSR graph, with each node's neighbour list sorted and deduplicated.
 */
inline NodeAdjacency build_node_adjacency(const Mesh& rMesh, std::size_t NumPoints,
                                          NodeAdjacencyKind Kind) {
    // Phase 1 (serial, deterministic): accumulate raw adjacency lists.
    std::vector<std::vector<std::int64_t>> adj(NumPoints);
    std::vector<std::int64_t> nodes;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t nc = cb.NumCells();
        const bool ragged = cb.IsRagged();
        const bool polyhedron = cb.IsPolyhedron();

        // Hoist the per-block edge table out of the per-cell loop: it is a
        // table lookup by cell type, which cannot change within a block.
        const std::vector<CellEdgePair>* p_edges = nullptr;
        if (Kind == NodeAdjacencyKind::Edge && !ragged && !polyhedron) {
            const std::vector<CellEdgePair>& e =
                detail::cell_refine_edges(cell_type_from_name(cb.Type()));
            if (e.empty())
                continue;  // unknown edge topology: contributes nothing (see the file docs)
            p_edges = &e;
        }

        for (std::size_t c = 0; c < nc; ++c) {
            if (Kind == NodeAdjacencyKind::Clique) {
                cell_node_ids(cb, c, NumPoints, nodes);
                node_adjacency_impl::add_clique(adj, nodes);
                continue;
            }

            // --- edge adjacency ---
            if (polyhedron) {
                // Each face contributes its own ring of edges.
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    std::pair<const std::int64_t*, std::size_t> face = cb.Face(c, f);
                    node_adjacency_impl::add_ring(adj, NumPoints, face.first, face.second);
                }
            } else if (ragged) {
                // A jagged polygon row is a closed loop of nodes.
                node_adjacency_impl::add_ring(adj, NumPoints, cb.Row(c), cb.RowSize(c));
            } else {
                // A known edge table (all seven linear types); blocks without
                // one were skipped wholesale above.
                const NDArray& conn = cb.Conn();
                const std::size_t row = c * cb.NodesPerCell();
                for (const CellEdgePair& e : *p_edges)
                    node_adjacency_impl::add_edge(adj, NumPoints,
                                                  detail::read_int(conn, row + e[0]),
                                                  detail::read_int(conn, row + e[1]));
            }
        }
    }

    // Phase 2 (parallel): sort + unique each node's list, dropping self-loops.
    parallel_for(NumPoints, [&](std::size_t i) {
        std::vector<std::int64_t>& v = adj[i];
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        v.erase(std::remove(v.begin(), v.end(), static_cast<std::int64_t>(i)), v.end());
    });

    // Phase 3: pack into CSR.
    NodeAdjacency csr;
    csr.mXadj.resize(NumPoints + 1);
    csr.mXadj[0] = 0;
    for (std::size_t i = 0; i < NumPoints; ++i)
        csr.mXadj[i + 1] = csr.mXadj[i] + static_cast<std::int64_t>(adj[i].size());
    csr.mAdj.resize(static_cast<std::size_t>(csr.mXadj[NumPoints]));
    parallel_for_bw(NumPoints, [&](std::size_t i) {
        std::copy(adj[i].begin(), adj[i].end(),
                  csr.mAdj.begin() + static_cast<std::ptrdiff_t>(csr.mXadj[i]));
    });
    return csr;
}

}  // namespace detail
}  // namespace meshioplusplus
