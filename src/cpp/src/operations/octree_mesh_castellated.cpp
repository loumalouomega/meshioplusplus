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
// octree_mesh, castellated output: the leaves as hexahedra, the inside /
// boundary-cell / surrogate-boundary regions, hanging nodes, and the signed
// distance and closest skin point at every node. Anonymous-namespace helpers
// are prefixed `octmc_` (the amalgamation concatenates every translation unit
// into one).

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "octree_mesh_internal.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The local hexahedron facet (detail/cell_faces) on each face direction, in
// the order of `linear_octree_directions(Face)`: -z, -y, -x, +x, +y, +z.
constexpr std::int64_t kOctmcFacet[6] = {4, 2, 0, 1, 3, 5};

NDArray octmc_ids(const std::vector<std::int64_t>& rIds, std::size_t Stride) {
    NDArray a = Stride == 1 ? NDArray::Uninit(DType::Int64, {rIds.size()})
                            : NDArray::Uninit(DType::Int64, {rIds.size() / Stride, Stride});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

// Classify every leaf whose (key, level) changed since the last pass.
void octmc_classify(const OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR,
                    const LinearOctree& rTree, const std::vector<std::size_t>& rWhich,
                    std::vector<std::uint8_t>& rInside) {
    const std::vector<Vec3> all = octree_mesh_leaf_centres(rS.mFrame, rTree);
    std::vector<Vec3> pts(rWhich.size());
    for (std::size_t k = 0; k < rWhich.size(); ++k)
        pts[k] = all[rWhich[k]];
    const std::vector<double> w = octree_mesh_winding(rS, pts, rO, rR);
    for (std::size_t k = 0; k < rWhich.size(); ++k)
        rInside[rWhich[k]] = w[k] > rO.mWindingThreshold;
}

}  // namespace

void octree_mesh_castellated(OctreeMeshSetup& rS, const OctreeMeshOptions& rO,
                             OctreeMeshResult& rR) {
    AdaptiveOctree& oct = rS.mOctree;
    LinearOctree& tree = oct.mTree;
    const auto faces = linear_octree_directions(OctreeAdjacency::Face);

    // Inside by the winding number at the leaf centres.
    std::vector<std::uint8_t> inside(tree.Size(), 0);
    {
        std::vector<std::size_t> every(tree.Size());
        for (std::size_t i = 0; i < every.size(); ++i)
            every[i] = i;
        octmc_classify(rS, rO, rR, tree, every, inside);
    }
    // Split an inside leaf beside a finer leaf that is not inside until every
    // exposed face of an inside leaf is a whole face of it.
    for (;;) {
        std::vector<std::uint8_t> mark(tree.Size(), 0);
        bool any = false;
        parallel_for(tree.Size(), [&](std::size_t i) {
            if (!inside[i])
                return;
            for (const auto& d : faces)
                for (std::int64_t j : linear_octree_neighbours(tree, i, d[0], d[1], d[2]))
                    if (tree.mLevel[static_cast<std::size_t>(j)] > tree.mLevel[i] &&
                        !inside[static_cast<std::size_t>(j)]) {
                        mark[i] = 1;
                        return;
                    }
        });
        for (std::uint8_t m : mark)
            any = any || m;
        if (!any)
            break;
        const std::vector<std::uint64_t> old_key = tree.mKey;
        const std::vector<std::uint8_t> old_level = tree.mLevel;
        std::vector<std::int64_t> split_origin, balance_origin;
        linear_octree_split(tree, mark, rO.mMaxCells, &split_origin);
        linear_octree_balance(tree, OctreeAdjacency::Full, rO.mMaxCells, &balance_origin);
        std::vector<std::int64_t> origin(tree.Size());
        for (std::size_t i = 0; i < tree.Size(); ++i)
            origin[i] = split_origin[static_cast<std::size_t>(balance_origin[i])];
        octree_adaptive_refresh(oct, rS.mFrame, rS.mField.mPrims, origin, std::move(oct.mCand));
        std::vector<std::uint8_t> next(tree.Size(), 0);
        std::vector<std::size_t> changed;
        for (std::size_t i = 0; i < tree.Size(); ++i) {
            const auto o = static_cast<std::size_t>(origin[i]);
            if (old_key[o] == tree.mKey[i] && old_level[o] == tree.mLevel[i])
                next[i] = inside[o];
            else
                changed.push_back(i);
        }
        inside.swap(next);
        octmc_classify(rS, rO, rR, tree, changed, inside);
    }
    octree_mesh_report_tree(tree, rR);

    // The kept leaves, in key order.
    std::vector<std::size_t> kept;
    for (std::size_t i = 0; i < tree.Size(); ++i) {
        const bool keep = rO.mFill == OctreeMeshFill::All ||
                          (rO.mFill == OctreeMeshFill::Inside && inside[i]) ||
                          (rO.mFill == OctreeMeshFill::Cover && (inside[i] || oct.mCut[i]));
        if (keep)
            kept.push_back(i);
    }
    if (static_cast<std::int64_t>(kept.size()) > rO.mMaxTets)
        throw std::invalid_argument(std::string(kOctreeMeshPrefix) + "the output would have " +
                                    std::to_string(kept.size()) + " cells, more than max_tets");
    std::vector<std::int64_t> cell_of(tree.Size(), -1);
    for (std::size_t n = 0; n < kept.size(); ++n)
        cell_of[kept[n]] = static_cast<std::int64_t>(n);

    // Points: the kept leaves' corners on the finest lattice, in key order.
    constexpr std::int64_t kCorner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    std::vector<std::uint64_t> corner_key(kept.size() * 8);
    parallel_for(kept.size(), [&](std::size_t n) {
        const auto lo = tree.Lo(kept[n]);
        const std::int64_t s = tree.Span(kept[n]);
        for (std::size_t c = 0; c < 8; ++c)
            corner_key[n * 8 + c] = linoct_key(lo[0] + kCorner[c][0] * s, lo[1] + kCorner[c][1] * s,
                                               lo[2] + kCorner[c][2] * s);
    });
    std::vector<std::uint64_t> point_key = corner_key;
    parallel_sort(point_key.begin(), point_key.end());
    point_key.erase(std::unique(point_key.begin(), point_key.end()), point_key.end());
    const std::size_t np = point_key.size();
    std::vector<Vec3> pts(np);
    parallel_for(np, [&](std::size_t p) {
        const auto x = linoct_decode(point_key[p]);
        pts[p] = rS.mFrame.Position(x[0], x[1], x[2]);
    });
    NDArray points = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
    for (std::size_t p = 0; p < np; ++p)
        for (std::size_t k = 0; k < 3; ++k)
            points.As<double>()[p * 3 + k] = pts[p][k];
    NDArray conn = NDArray::Uninit(DType::Int64, {kept.size(), std::size_t{8}});
    std::int64_t* cp = conn.As<std::int64_t>();
    parallel_for(corner_key.size(), [&](std::size_t i) {
        cp[i] = static_cast<std::int64_t>(
            std::lower_bound(point_key.begin(), point_key.end(), corner_key[i]) -
            point_key.begin());
    });
    Mesh& m = rR.mMesh;
    m.AssignPoints(std::move(points));
    m.AddCellBlock("hexahedron", std::move(conn));

    // Cell data.
    NDArray level = NDArray::Uninit(DType::Int64, {kept.size()});
    NDArray cut = NDArray::Uninit(DType::Int64, {kept.size()});
    for (std::size_t n = 0; n < kept.size(); ++n) {
        level.As<std::int64_t>()[n] = tree.mLevel[kept[n]];
        cut.As<std::int64_t>()[n] = oct.mCut[kept[n]];
    }
    {
        std::vector<NDArray> lv, ct;
        lv.push_back(std::move(level));
        ct.push_back(std::move(cut));
        m.AddCellData(kOctreeLevelName, std::move(lv));
        m.AddCellData(kOctreeCutName, std::move(ct));
    }

    // Hanging nodes: on a kept leaf's boundary without being its corner.
    NDArray hanging = NDArray::Uninit(DType::Int64, {np});
    std::int64_t* hp = hanging.As<std::int64_t>();
    parallel_for(np, [&](std::size_t p) {
        const auto x = linoct_decode(point_key[p]);
        hp[p] = 0;
        for (std::int64_t c = 0; c < 8 && !hp[p]; ++c) {
            const std::int64_t q[3] = {x[0] - (c & 1), x[1] - ((c >> 1) & 1),
                                       x[2] - ((c >> 2) & 1)};
            const std::int64_t j = tree.Locate(q[0], q[1], q[2]);
            if (j < 0 || cell_of[static_cast<std::size_t>(j)] < 0)
                continue;
            const auto lo = tree.Lo(static_cast<std::size_t>(j));
            const std::int64_t s = tree.Span(static_cast<std::size_t>(j));
            for (std::size_t k = 0; k < 3; ++k)
                if (x[k] != lo[k] && x[k] != lo[k] + s)
                    hp[p] = 1;
        }
    });
    rR.mNumHanging = 0;
    for (std::size_t p = 0; p < np; ++p)
        rR.mNumHanging += hp[p];
    m.AddPointData(kOctreeHangingName, std::move(hanging));

    // Regions.
    std::vector<std::int64_t> in_cells, cut_cells, sides;
    for (std::size_t n = 0; n < kept.size(); ++n) {
        const std::size_t i = kept[n];
        if (oct.mCut[i])
            cut_cells.push_back(static_cast<std::int64_t>(n));
        if (!inside[i])
            continue;
        in_cells.push_back(static_cast<std::int64_t>(n));
        for (std::size_t f = 0; f < faces.size(); ++f) {
            const auto nb =
                linear_octree_neighbours(tree, i, faces[f][0], faces[f][1], faces[f][2]);
            bool exposed = true;
            for (std::int64_t j : nb)
                exposed = exposed && !inside[static_cast<std::size_t>(j)];
            if (exposed) {
                sides.push_back(static_cast<std::int64_t>(n));
                sides.push_back(kOctmcFacet[f]);
            }
        }
    }
    rR.mNumInside = static_cast<std::int64_t>(in_cells.size());
    rR.mNumBoundaryCells = static_cast<std::int64_t>(cut_cells.size());
    rR.mNumSurrogateFaces = static_cast<std::int64_t>(sides.size() / 2);
    m.AddRegion(Region(kOctreeMeshInsideRegion, RegionKind::Cell, 3, -1, octmc_ids(in_cells, 1)));
    m.AddRegion(
        Region(kOctreeMeshBoundaryCellsRegion, RegionKind::Cell, 3, -1, octmc_ids(cut_cells, 1)));
    m.AddRegion(Region(kOctreeMeshSurrogateRegion, RegionKind::Side, 2, -1, octmc_ids(sides, 2)));

    // The signed distance and the closest skin point at every node.
    if (rO.mAttachDistance) {
        std::vector<double> dist;
        std::vector<Vec3> closest;
        octree_mesh_signed_distance(rS, pts, rO, rR, dist, closest);
        NDArray d = NDArray::Uninit(DType::Float64, {np});
        NDArray c = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
        for (std::size_t p = 0; p < np; ++p) {
            d.As<double>()[p] = dist[p];
            for (std::size_t k = 0; k < 3; ++k)
                c.As<double>()[p * 3 + k] = closest[p][k];
        }
        m.AddPointData(kOctreeMeshDistanceName, std::move(d));
        m.AddPointData(kOctreeMeshClosestPointName, std::move(c));
    }
}

}  // namespace detail
}  // namespace meshioplusplus
