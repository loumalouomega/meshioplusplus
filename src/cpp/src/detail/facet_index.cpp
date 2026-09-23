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
#include <string>

// Project includes
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {
namespace detail {

FacetIndex::FacetIndex(const Mesh& rMesh, const FacetIndexOptions& rOptions) {
    std::int64_t base = 0;
    std::int64_t corners[4];
    const auto add = [&](std::size_t N, std::int64_t Cell, std::int64_t Facet) {
        FacetHit& hit = mMap[FacetKey(corners, N)];
        if (hit.mCount == 0)
            hit.mFirst = FacetOwner{Cell, Facet};
        else if (hit.mCount == 1)
            hit.mSecond = FacetOwner{Cell, Facet};
        ++hit.mCount;
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
}

const FacetHit* FacetIndex::Find(const std::int64_t* pCorners, std::size_t N) const {
    const auto it = mMap.find(FacetKey(pCorners, N));
    return it == mMap.end() ? nullptr : &it->second;
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
