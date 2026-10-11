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
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// Calls `rEmit(u, v)` once per raw adjacency entry "v is a neighbour of u":
// every ordered pair of a cell's distinct nodes (Clique), or both directions
// of every in-range, non-loop edge (Edge). The order is fixed (blocks, cells,
// then pairs/edges), but only the set each node ends up with matters: rows
// are sorted and deduplicated afterwards.
template <class Emit>
void nodeadj_visit(const Mesh& rMesh, std::size_t NumPoints, NodeAdjacencyKind Kind,
                   const Emit& rEmit) {
    const auto edge = [&](std::int64_t u, std::int64_t v) {
        if (u < 0 || v < 0 || u == v)
            return;
        if (static_cast<std::size_t>(u) >= NumPoints || static_cast<std::size_t>(v) >= NumPoints)
            return;
        rEmit(u, v);
        rEmit(v, u);
    };
    // The closed ring of edges around an ordered node loop (a polygon row or
    // one face of a polyhedron).
    const auto ring = [&](const std::int64_t* pNodes, std::size_t Count) {
        if (Count < 2)
            return;
        for (std::size_t k = 0; k < Count; ++k)
            edge(pNodes[k], pNodes[(k + 1) % Count]);
    };
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
                // The clique over the cell's DISTINCT nodes: a polyhedron
                // lists each corner once per face, and pairing them by
                // position would add every edge k^2 times and self-loops.
                cell_node_ids(cb, c, NumPoints, nodes);
                std::sort(nodes.begin(), nodes.end());
                nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
                for (const std::int64_t u : nodes)
                    for (const std::int64_t v : nodes)
                        if (u != v)
                            rEmit(u, v);
                continue;
            }

            // --- edge adjacency ---
            if (polyhedron) {
                // Each face contributes its own ring of edges.
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    std::pair<const std::int64_t*, std::size_t> face = cb.Face(c, f);
                    ring(face.first, face.second);
                }
            } else if (ragged) {
                // A jagged polygon row is a closed loop of nodes.
                ring(cb.Row(c), cb.RowSize(c));
            } else {
                // A known edge table (all seven linear types); blocks without
                // one were skipped wholesale above.
                const NDArray& conn = cb.Conn();
                const std::size_t row = c * cb.NodesPerCell();
                for (const CellEdgePair& e : *p_edges)
                    edge(detail::read_int(conn, row + e[0]), detail::read_int(conn, row + e[1]));
            }
        }
    }
}

// The CSR build over raw entries stored as `Id` (32-bit when every node id
// fits, which halves the scratch): count each node's raw entries, fill them
// into one flat array, then sort and deduplicate each row in place and pack.
template <class Id>
NodeAdjacency nodeadj_build(const Mesh& rMesh, std::size_t NumPoints, NodeAdjacencyKind Kind) {
    // Pass 1 (serial): raw entries per node, then their start offsets.
    std::vector<std::size_t> start(NumPoints + 1, 0);
    nodeadj_visit(rMesh, NumPoints, Kind,
                  [&](std::int64_t u, std::int64_t) { ++start[static_cast<std::size_t>(u) + 1]; });
    for (std::size_t i = 0; i < NumPoints; ++i)
        start[i + 1] += start[i];

    // Pass 2 (serial): the raw entries, row by row.
    std::vector<Id> flat(start[NumPoints]);
    {
        std::vector<std::size_t> cursor(start.begin(), start.end() - 1);
        nodeadj_visit(rMesh, NumPoints, Kind, [&](std::int64_t u, std::int64_t v) {
            flat[cursor[static_cast<std::size_t>(u)]++] = static_cast<Id>(v);
        });
    }

    // Pass 3 (parallel): sort + unique each node's row in place. No self-loop
    // survives the visit, so a row is exactly the node's neighbour set.
    std::vector<std::int64_t> count(NumPoints, 0);
    parallel_for(NumPoints, [&](std::size_t i) {
        Id* lo = flat.data() + start[i];
        Id* hi = flat.data() + start[i + 1];
        std::sort(lo, hi);
        count[i] = std::unique(lo, hi) - lo;
    });

    // Pass 4: pack into CSR.
    NodeAdjacency csr;
    csr.mXadj.resize(NumPoints + 1);
    csr.mXadj[0] = 0;
    for (std::size_t i = 0; i < NumPoints; ++i)
        csr.mXadj[i + 1] = csr.mXadj[i] + count[i];
    csr.mAdj.resize(static_cast<std::size_t>(csr.mXadj[NumPoints]));
    parallel_for_bw(NumPoints, [&](std::size_t i) {
        std::copy(flat.data() + start[i], flat.data() + start[i] + count[i],
                  csr.mAdj.begin() + static_cast<std::ptrdiff_t>(csr.mXadj[i]));
    });
    return csr;
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
    // One flat array of raw entries rather than a vector per node: a
    // polyhedral mesh's per-node vectors held 8 KiB each at 10M cells
    // (roadmap §6.1). The rows, and so the CSR, are the same sets.
    if (NumPoints <= static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
        return nodeadj_build<std::uint32_t>(rMesh, NumPoints, Kind);
    return nodeadj_build<std::int64_t>(rMesh, NumPoints, Kind);
}

}  // namespace detail
}  // namespace meshioplusplus
