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

// System includes
#include <algorithm>
#include <utility>

// Project includes
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// Adds every ordered pair of a cell's nodes (an element clique). Duplicates
// and self-loops are removed later by the dedup pass.
void nodeadj_add_clique(std::vector<std::vector<std::int64_t>>& rAdj,
                        const std::vector<std::int64_t>& rNodes) {
    const std::size_t k = rNodes.size();
    for (std::size_t a = 0; a < k; ++a) {
        std::vector<std::int64_t>& row = rAdj[static_cast<std::size_t>(rNodes[a])];
        for (std::size_t b = 0; b < k; ++b)
            if (a != b)
                row.push_back(rNodes[b]);
    }
}

// Adds the undirected edge `(u, v)` in both directions, skipping out-of-range
// endpoints and self-loops.
void nodeadj_add_edge(std::vector<std::vector<std::int64_t>>& rAdj, std::size_t NumPoints,
                      std::int64_t u, std::int64_t v) {
    if (u < 0 || v < 0 || u == v)
        return;
    if (static_cast<std::size_t>(u) >= NumPoints || static_cast<std::size_t>(v) >= NumPoints)
        return;
    rAdj[static_cast<std::size_t>(u)].push_back(v);
    rAdj[static_cast<std::size_t>(v)].push_back(u);
}

// Adds the closed ring of edges around an ordered node loop (a polygon row or
// one face of a polyhedron).
void nodeadj_add_ring(std::vector<std::vector<std::int64_t>>& rAdj, std::size_t NumPoints,
                      const std::int64_t* pNodes, std::size_t Count) {
    if (Count < 2)
        return;
    for (std::size_t k = 0; k < Count; ++k)
        nodeadj_add_edge(rAdj, NumPoints, pNodes[k], pNodes[(k + 1) % Count]);
}

}  // namespace

void cell_node_ids(const Mesh::CellView& rBlock, std::size_t Cell, std::size_t NumPoints,
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

bool node_edge_topology_known(const Mesh::CellView& rBlock) {
    if (rBlock.IsPolyhedron() || rBlock.IsRagged())
        return true;
    return !detail::cell_refine_edges(cell_type_from_name(rBlock.Type())).empty();
}

NodeAdjacency build_node_adjacency(const Mesh& rMesh, std::size_t NumPoints,
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
                nodeadj_add_clique(adj, nodes);
                continue;
            }

            // --- edge adjacency ---
            if (polyhedron) {
                // Each face contributes its own ring of edges.
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    std::pair<const std::int64_t*, std::size_t> face = cb.Face(c, f);
                    nodeadj_add_ring(adj, NumPoints, face.first, face.second);
                }
            } else if (ragged) {
                // A jagged polygon row is a closed loop of nodes.
                nodeadj_add_ring(adj, NumPoints, cb.Row(c), cb.RowSize(c));
            } else {
                // A known edge table (all seven linear types); blocks without
                // one were skipped wholesale above.
                const NDArray& conn = cb.Conn();
                const std::size_t row = c * cb.NodesPerCell();
                for (const CellEdgePair& e : *p_edges)
                    nodeadj_add_edge(adj, NumPoints, detail::read_int(conn, row + e[0]),
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
