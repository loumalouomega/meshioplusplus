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
// Sorted-corner facet lookup shared by the formats that name a facet by its
// nodes. Python twin: src/python/meshioplusplus/_facets.py.

// System includes
#include <algorithm>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "slot_runs.hpp"

namespace meshioplusplus {
namespace detail {

FacetIndex::FacetIndex(const Mesh& rMesh, const FacetIndexOptions& rOptions) {
    // Every facet occurrence, in block-major (cell, local facet) order -- the
    // order the owners are ranked in.
    std::vector<FlatFacetKey> keys;
    std::vector<FacetOwner> owners;
    std::int64_t base = 0;
    std::int64_t corners[4];
    const auto add = [&](std::size_t N, std::int64_t Cell, std::int64_t Facet) {
        keys.push_back(flat_facet_key(corners, N));
        owners.push_back(FacetOwner{Cell, Facet});
    };
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t n_cells = cb.NumCells();
        if (cb.IsRagged()) {
            base += static_cast<std::int64_t>(n_cells);
            continue;
        }
        const CellType type = cell_type_from_name(std::string(cb.Type()));
        const int dim = cell_type_dimension(type);
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        const bool self =
            rOptions.mSurfaceSelf && (type == CellType::Triangle || type == CellType::Quad);
        if (dim == 3 && rOptions.mSolidFaces) {
            const auto& faces = cell_faces(type);
            for (std::size_t r = 0; r < n_cells; ++r)
                for (std::size_t f = 0; f < faces.size(); ++f) {
                    for (std::size_t c = 0; c < faces[f].mNumCorners; ++c)
                        corners[c] = read_int(conn, r * k + faces[f].mNodes[c]);
                    add(faces[f].mNumCorners, base + static_cast<std::int64_t>(r),
                        static_cast<std::int64_t>(f));
                }
        } else if (dim == 2 && (rOptions.mSurfaceEdges || self)) {
            const auto& edges = cell_edges(type);
            for (std::size_t r = 0; r < n_cells; ++r) {
                const std::int64_t g = base + static_cast<std::int64_t>(r);
                if (self) {
                    for (std::size_t c = 0; c < k; ++c)
                        corners[c] = read_int(conn, r * k + c);
                    add(k, g, 0);
                }
                if (rOptions.mSurfaceEdges)
                    for (std::size_t e = 0; e < edges.size(); ++e) {
                        for (std::size_t c = 0; c < 2; ++c)
                            corners[c] = read_int(conn, r * k + edges[e].mNodes[c]);
                        add(2, g, static_cast<std::int64_t>(e));
                    }
            }
        }
        base += static_cast<std::int64_t>(n_cells);
    }

    // Runs of equal keys come out in ascending key order, each holding its
    // occurrences in the order above: the first two owners and the count are
    // what a serial insert kept.
    const SlotRuns runs = group_slots_sorted(keys);  // sorted: Find binary-searches it
    mKeys.resize(runs.NumRuns());
    mHits.resize(runs.NumRuns());
    parallel_for(runs.NumRuns(), [&](std::size_t r) {
        const std::uint64_t* slot = runs.Begin(r);
        mKeys[r] = keys[slot[0]];
        FacetHit& hit = mHits[r];
        hit.mFirst = owners[slot[0]];
        if (runs.Size(r) > 1)
            hit.mSecond = owners[slot[1]];
        hit.mCount = runs.Size(r);
    });
}

const FacetHit* FacetIndex::Find(const std::int64_t* pCorners, std::size_t N) const {
    if (N > FacetKey::kInline)
        return nullptr;  // no indexed facet has more than four corners
    const FlatFacetKey key = flat_facet_key(pCorners, N);
    const auto it = std::lower_bound(mKeys.begin(), mKeys.end(), key);
    if (it == mKeys.end() || *it != key)
        return nullptr;
    return &mHits[static_cast<std::size_t>(it - mKeys.begin())];
}

bool facet_nodes(const Mesh& rMesh, std::int64_t Cell, std::int64_t Facet, CellType& rType,
                 std::vector<std::int64_t>& rNodes) {
    const auto [block, row] = global_to_block_row(block_bases(rMesh), Cell);
    if (block == static_cast<std::size_t>(-1) || Facet < 0)
        return false;
    const auto cb = rMesh.Cells(block);
    if (cb.IsRagged())
        return false;
    const CellType type = cell_type_from_name(std::string(cb.Type()));
    const NDArray& conn = cb.Conn();
    const std::size_t k = cb.NodesPerCell();
    const std::size_t r = static_cast<std::size_t>(row);
    const int dim = cell_type_dimension(type);
    rNodes.clear();
    if (dim == 3) {
        const auto& faces = cell_faces(type);
        if (static_cast<std::size_t>(Facet) >= faces.size())
            return false;
        const CellFaceDef& face = faces[static_cast<std::size_t>(Facet)];
        rType = face.mFaceType;
        for (std::size_t c = 0; c < face.mNumNodes; ++c)
            rNodes.push_back(read_int(conn, r * k + face.mNodes[c]));
        return true;
    }
    if (dim == 2) {
        const auto& edges = cell_edges(type);
        if (static_cast<std::size_t>(Facet) >= edges.size())
            return false;
        const CellEdgeDef& edge = edges[static_cast<std::size_t>(Facet)];
        rType = edge.mEdgeType;
        for (std::size_t c = 0; c < edge.mNumNodes; ++c)
            rNodes.push_back(read_int(conn, r * k + edge.mNodes[c]));
        return true;
    }
    return false;
}

}  // namespace detail
}  // namespace meshioplusplus
