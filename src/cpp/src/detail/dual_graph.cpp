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
// The facet-sharing dual graph (moved out of operations/partition.cpp, where
// it fed only KaHIP, so the adjacency operation can return the same graph).
// See dual_graph.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "dual_graph.hpp"
#include "slot_runs.hpp"
#include "typed_view.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The corner-node rows of a cell type's facets in the dual-graph sense: faces
// for 3D cells, edges for 2D cells, empty (isolated vertex) otherwise.
struct DualGraphFacetDef {
    std::uint8_t mNumCorners;
    std::array<std::uint8_t, 4> mNodes;
};

std::vector<DualGraphFacetDef> dualgraph_facets_for(CellType type, int dualDim) {
    std::vector<DualGraphFacetDef> defs;
    if (dualDim == 3) {
        for (const CellFaceDef& f : cell_faces(type)) {
            DualGraphFacetDef d;
            d.mNumCorners = f.mNumCorners;
            for (std::uint8_t k = 0; k < f.mNumCorners; ++k)
                d.mNodes[k] = f.mNodes[k];
            defs.push_back(d);
        }
    } else if (dualDim == 2) {
        for (const CellEdgeDef& e : cell_edges(type)) {
            DualGraphFacetDef d;
            d.mNumCorners = e.mNumCorners;
            for (std::uint8_t k = 0; k < e.mNumCorners; ++k)
                d.mNodes[k] = e.mNodes[k];
            defs.push_back(d);
        }
    }
    return defs;
}

bool dualgraph_is_polygon(const Mesh::CellView& rBlock) {
    return !rBlock.IsPolyhedron() && (rBlock.IsRagged() || rBlock.Type().rfind("polygon", 0) == 0);
}

}  // namespace

DualGraphCsr build_facet_dual_graph(const Mesh& rMesh, std::size_t total, bool PolygonRings) {
    // The dual dimension: faces when any 3D cell exists, else edges. Blocks of
    // lower dimension (or without a facet table) become isolated vertices.
    int dual_dim = 0;
    for (const auto cb : rMesh.CellRange()) {
        int d = cb.IsPolyhedron() ? 3 : cell_type_dimension(cell_type_from_name(cb.Type()));
        if (PolygonRings && dualgraph_is_polygon(cb))
            d = 2;
        dual_dim = std::max(dual_dim, d);
    }

    // Per-block facet tables + the flat facet-slot layout. A polyhedron or
    // polygon-ring block has a variable facet count per cell (mFaceStart).
    struct BlockDesc {
        const NDArray* mpConn = nullptr;
        std::size_t mNpc = 0;
        std::size_t mNumCells = 0;
        std::vector<DualGraphFacetDef> mFacets;
        std::size_t mFirstFacet = 0;
        std::int64_t mGlobalCellBase = 0;
        bool mPolyhedron = false;
        bool mRing = false;
        std::size_t mBlock = 0;
        std::vector<std::size_t> mFaceStart;  // mNumCells + 1
    };
    std::vector<BlockDesc> descs;
    std::size_t total_facets = 0;
    std::int64_t global_cell_base = 0;
    std::size_t block_index = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t this_block = block_index++;
        const std::int64_t base = global_cell_base;
        global_cell_base += static_cast<std::int64_t>(cb.NumCells());
        const bool ring = PolygonRings && dualgraph_is_polygon(cb);
        if (cb.IsPolyhedron() || ring) {
            // A polyhedron carries its own faces, and a polygon row its ring of
            // edges, so neither needs a table.
            if (dual_dim != (ring ? 2 : 3))
                continue;
            BlockDesc d;
            d.mPolyhedron = !ring;
            d.mRing = ring;
            d.mBlock = this_block;
            d.mNumCells = cb.NumCells();
            d.mFirstFacet = total_facets;
            d.mGlobalCellBase = base;
            if (ring && !cb.IsRagged()) {
                d.mpConn = &cb.Conn();
                d.mNpc = cb.NodesPerCell();
            }
            d.mFaceStart.reserve(d.mNumCells + 1);
            std::size_t at = 0;
            d.mFaceStart.push_back(0);
            for (std::size_t c = 0; c < d.mNumCells; ++c) {
                at += ring ? (cb.IsRagged() ? cb.RowSize(c) : d.mNpc) : cb.NumFaces(c);
                d.mFaceStart.push_back(at);
            }
            total_facets += at;
            descs.push_back(std::move(d));
            continue;
        }
        if (cb.IsRagged())
            continue;  // isolated vertices (no facet table)
        const CellType ct = cell_type_from_name(cb.Type());
        if (cell_type_dimension(ct) != dual_dim)
            continue;
        std::vector<DualGraphFacetDef> facets = dualgraph_facets_for(ct, dual_dim);
        if (facets.empty())
            continue;
        BlockDesc d;
        d.mpConn = &cb.Conn();
        d.mNpc = cb.NodesPerCell();
        d.mNumCells = cb.NumCells();
        d.mFacets = std::move(facets);
        d.mFirstFacet = total_facets;
        d.mGlobalCellBase = base;
        total_facets += d.mNumCells * d.mFacets.size();
        descs.push_back(std::move(d));
    }

    // Phase 1: each facet slot's owning cell, in parallel (disjoint slots).
    std::vector<std::int64_t> parent(total_facets);
    std::vector<std::size_t> first_facets(descs.size());
    using DualGraphCellView = decltype(rMesh.Cells(0));
    std::vector<std::optional<DualGraphCellView>> views(descs.size());
    std::vector<std::optional<Int64View>> conn_views(descs.size());
    for (std::size_t di = 0; di < descs.size(); ++di) {
        const BlockDesc& d = descs[di];
        first_facets[di] = d.mFirstFacet;
        if (d.mPolyhedron || d.mRing)
            views[di].emplace(rMesh.Cells(d.mBlock));
        if (d.mpConn)
            conn_views[di].emplace(*d.mpConn);
        if (d.mPolyhedron || d.mRing) {
            parallel_for(d.mNumCells, [&](std::size_t cell) {
                for (std::size_t f = d.mFaceStart[cell]; f < d.mFaceStart[cell + 1]; ++f)
                    parent[d.mFirstFacet + f] = d.mGlobalCellBase + static_cast<std::int64_t>(cell);
            });
            continue;
        }
        const std::size_t fpc = d.mFacets.size();
        parallel_for(d.mNumCells * fpc, [&](std::size_t j) {
            parent[d.mFirstFacet + j] = d.mGlobalCellBase + static_cast<std::int64_t>(j / fpc);
        });
    }

    // Phase 2: group the facet slots by key (detail::group_facets, roadmap
    // 5.4.3: the keys flattened to the width they need, no record per facet).
    // A key is the facet's corner ids, of any arity, so a hexahedron and a
    // polyhedron that meet on a face become neighbours. A run's head is its
    // first slot in stored order, the owner a first-seen map would keep; every
    // later slot whose parent differs from the owner adds an owner<->parent
    // edge, so the occurrences of a non-manifold facet beyond the second all
    // connect to the first owner.
    const auto desc_of = [&](std::size_t ri) -> std::size_t {
        if (first_facets.size() == 1)
            return 0;
        return static_cast<std::size_t>(
                   std::upper_bound(first_facets.begin(), first_facets.end(), ri) -
                   first_facets.begin()) -
               1;
    };
    const auto size_of = [&](std::size_t ri) -> std::size_t {
        const std::size_t di = desc_of(ri);
        const BlockDesc& d = descs[di];
        if (d.mRing)
            return 2;
        if (d.mPolyhedron) {
            const auto cell = static_cast<std::size_t>(parent[ri] - d.mGlobalCellBase);
            return views[di]->Face(cell, ri - d.mFirstFacet - d.mFaceStart[cell]).second;
        }
        return d.mFacets[(ri - d.mFirstFacet) % d.mFacets.size()].mNumCorners;
    };
    const auto fill = [&](std::size_t ri, std::int64_t* pOut) {
        const std::size_t di = desc_of(ri);
        const BlockDesc& d = descs[di];
        const auto cell = static_cast<std::size_t>(parent[ri] - d.mGlobalCellBase);
        if (d.mPolyhedron) {
            const auto face = views[di]->Face(cell, ri - d.mFirstFacet - d.mFaceStart[cell]);
            std::copy(face.first, face.first + face.second, pOut);
            return;
        }
        if (d.mRing) {
            const std::size_t n = d.mFaceStart[cell + 1] - d.mFaceStart[cell];
            const std::size_t k = ri - d.mFirstFacet - d.mFaceStart[cell];
            const auto node = [&](std::size_t at) {
                return conn_views[di] ? (*conn_views[di])[cell * d.mNpc + at]
                                      : views[di]->Row(cell)[at];
            };
            pOut[0] = node(k);
            pOut[1] = node((k + 1) % n);
            return;
        }
        const DualGraphFacetDef& facet = d.mFacets[(ri - d.mFirstFacet) % d.mFacets.size()];
        const Int64View& conn_v = *conn_views[di];
        for (std::uint8_t k = 0; k < facet.mNumCorners; ++k)
            pOut[k] = conn_v[cell * d.mNpc + facet.mNodes[k]];
    };
    const SlotRuns runs = group_facets(total_facets, rMesh.NumPoints(), size_of, fill);
    const std::size_t num_runs = runs.NumRuns();

    // Flat adjacency with duplicates: degree count, prefix sum, fill (both serial
    // over the runs, plain integer work, so the lists do not depend on a backend).
    std::vector<std::int64_t> dup_start(total + 1, 0);
    for (std::size_t r = 0; r < num_runs; ++r) {
        const std::int64_t owner = parent[runs.Head(r)];
        for (const std::uint64_t* p = runs.Begin(r); p != runs.End(r); ++p) {
            const std::int64_t other = parent[*p];
            if (other != owner) {
                ++dup_start[static_cast<std::size_t>(owner) + 1];
                ++dup_start[static_cast<std::size_t>(other) + 1];
            }
        }
    }
    for (std::size_t i = 0; i < total; ++i)
        dup_start[i + 1] += dup_start[i];
    std::vector<std::int64_t> dup(static_cast<std::size_t>(dup_start[total]));
    std::vector<std::int64_t> cursor(dup_start.begin(), dup_start.end() - 1);
    for (std::size_t r = 0; r < num_runs; ++r) {
        const std::int64_t owner = parent[runs.Head(r)];
        for (const std::uint64_t* p = runs.Begin(r); p != runs.End(r); ++p) {
            const std::int64_t other = parent[*p];
            if (other != owner) {
                dup[static_cast<std::size_t>(cursor[static_cast<std::size_t>(owner)]++)] = other;
                dup[static_cast<std::size_t>(cursor[static_cast<std::size_t>(other)]++)] = owner;
            }
        }
    }

    // Dedupe neighbour lists (cells can share more than one facet) in place and
    // pack CSR. Sorting makes each list independent of the order it was filled in.
    std::vector<std::int64_t> kept(total);
    parallel_for(total, [&](std::size_t i) {
        auto first = dup.begin() + dup_start[i];
        auto last = dup.begin() + dup_start[i + 1];
        std::sort(first, last);
        kept[i] = static_cast<std::int64_t>(std::unique(first, last) - first);
    });
    DualGraphCsr csr;
    csr.mXadj.resize(total + 1);
    csr.mXadj[0] = 0;
    for (std::size_t i = 0; i < total; ++i)
        csr.mXadj[i + 1] = csr.mXadj[i] + kept[i];
    csr.mAdjncy.reserve(static_cast<std::size_t>(csr.mXadj[total]));
    for (std::size_t i = 0; i < total; ++i)
        csr.mAdjncy.insert(csr.mAdjncy.end(), dup.begin() + dup_start[i],
                           dup.begin() + dup_start[i] + kept[i]);
    return csr;
}

}  // namespace detail
}  // namespace meshioplusplus
