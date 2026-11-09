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
// Node-connected components of a mesh's cells, by union-find over the cells
// that share each node. Shared by split(by="component") and select_region. See
// cell_components.hpp for the numbering.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "cell_components.hpp"
#include "dual_graph.hpp"
#include "typed_view.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

std::int64_t ccomp_find(std::vector<std::int64_t>& rParent, std::int64_t x) {
    while (rParent[static_cast<std::size_t>(x)] != x) {
        rParent[static_cast<std::size_t>(x)] =
            rParent[static_cast<std::size_t>(rParent[static_cast<std::size_t>(x)])];
        x = rParent[static_cast<std::size_t>(x)];
    }
    return x;
}

void ccomp_union(std::vector<std::int64_t>& rParent, std::int64_t a, std::int64_t b) {
    const std::int64_t ra = ccomp_find(rParent, a), rb = ccomp_find(rParent, b);
    if (ra != rb)
        rParent[static_cast<std::size_t>(std::max(ra, rb))] = std::min(ra, rb);
}

}  // namespace

CellComponents cell_components(const Mesh& rMesh) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::size_t n = rMesh.NumPoints();

    std::vector<std::size_t> block_base(nblocks + 1, 0);
    std::size_t bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        block_base[bi + 1] = block_base[bi] + cb.NumCells();
        ++bi;
    }
    const std::size_t total = block_base[nblocks];

    std::vector<std::int64_t> parent(total);
    parallel_for_bw(total, [&](std::size_t i) { parent[i] = static_cast<std::int64_t>(i); });

    // Union cells that share a node (the first cell per node as the anchor). The
    // sweep is order-dependent and stays serial; a rectangular block's
    // connectivity is read once, not through a dtype switch per id.
    std::vector<std::int64_t> node_first(n, -1);
    bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const auto visit = [&](std::size_t c, std::int64_t node) {
            const std::int64_t gc = static_cast<std::int64_t>(block_base[bi] + c);
            std::int64_t& first = node_first[static_cast<std::size_t>(node)];
            if (first < 0)
                first = gc;
            else
                ccomp_union(parent, first, gc);
        };
        if (cb.IsPolyhedron()) {
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                for (std::size_t face = 0; face < cb.NumFaces(c); ++face) {
                    const auto fc = cb.Face(c, face);
                    for (std::size_t k = 0; k < fc.second; ++k)
                        visit(c, fc.first[k]);
                }
        } else if (cb.IsRagged()) {
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                for (std::size_t k = 0; k < cb.RowSize(c); ++k)
                    visit(c, cb.Row(c)[k]);
        } else {
            const Int64View conn(cb.Conn());
            const std::int64_t* p = conn.Data();
            const std::size_t npc = cb.NodesPerCell();
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                for (std::size_t k = 0; k < npc; ++k)
                    visit(c, p[c * npc + k]);
        }
        ++bi;
    }

    // Every root is its component's smallest cell (a union points the larger
    // root at the smaller) and parent[x] <= x, so one ascending pass resolves
    // each cell's root. Components are numbered in first-seen (ascending cell)
    // order: a component's first cell is its root, so the label of a cell is
    // the number of roots before its own -- a scan of the root flags.
    for (std::size_t gc = 0; gc < total; ++gc)
        parent[gc] = parent[static_cast<std::size_t>(parent[gc])];
    std::vector<std::uint64_t> is_root(total);
    parallel_for_bw(total, [&](std::size_t gc) {
        is_root[gc] = parent[gc] == static_cast<std::int64_t>(gc) ? 1 : 0;
    });
    std::vector<std::uint64_t> label_at(total);
    const std::uint64_t count =
        parallel_exclusive_scan(is_root.data(), total, label_at.data(), std::uint64_t{0});

    CellComponents out;
    out.mCount = static_cast<std::int64_t>(count);
    out.mLabel.resize(total);
    parallel_for_bw(total, [&](std::size_t gc) {
        out.mLabel[gc] = static_cast<std::int64_t>(label_at[static_cast<std::size_t>(parent[gc])]);
    });
    return out;
}

CellComponents cell_components_by_edge(const Mesh& rMesh) {
    std::size_t total = 0;
    for (const auto cb : rMesh.CellRange())
        total += cb.NumCells();
    const DualGraphCsr graph = build_facet_dual_graph(rMesh, total, true);
    CellComponents out;
    out.mLabel.assign(total, -1);
    std::vector<std::int64_t> stack;
    for (std::size_t seed = 0; seed < total; ++seed) {
        if (out.mLabel[seed] >= 0)
            continue;
        out.mLabel[seed] = out.mCount;
        stack.assign(1, static_cast<std::int64_t>(seed));
        while (!stack.empty()) {
            const auto c = static_cast<std::size_t>(stack.back());
            stack.pop_back();
            for (std::int64_t e = graph.mXadj[c]; e < graph.mXadj[c + 1]; ++e) {
                const auto nb =
                    static_cast<std::size_t>(graph.mAdjncy[static_cast<std::size_t>(e)]);
                if (out.mLabel[nb] < 0) {
                    out.mLabel[nb] = out.mCount;
                    stack.push_back(static_cast<std::int64_t>(nb));
                }
            }
        }
        ++out.mCount;
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
