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

// octree_mesh, hexahedral output: the grid cells of the skin's level inside
// the skin, made manifold, pillowed, fitted to the skin and its features, and
// untangled. Anonymous-namespace helpers are prefixed `octmh_` (the
// amalgamation concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "../detail/feature_snap.hpp"
#include "../detail/hex_fit.hpp"
#include "octree_mesh_internal.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

NDArray octmh_ids(const std::vector<std::int64_t>& rIds, std::size_t Stride) {
    NDArray a = Stride == 1 ? NDArray::Uninit(DType::Int64, {rIds.size()})
                            : NDArray::Uninit(DType::Int64, {rIds.size() / Stride, Stride});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

// The nearest point to @p rP on the polyline through the skin points @p rChain.
Vec3 octmh_on_chain(const std::vector<Vec3>& rSkin, const std::vector<std::int64_t>& rChain,
                    const Vec3& rP, double& rDistance) {
    Vec3 best = rSkin[static_cast<std::size_t>(rChain.front())];
    rDistance = std::numeric_limits<double>::infinity();
    for (std::size_t s = 0; s + 1 < rChain.size(); ++s) {
        const Vec3& a = rSkin[static_cast<std::size_t>(rChain[s])];
        const Vec3& b = rSkin[static_cast<std::size_t>(rChain[s + 1])];
        const Vec3 ab = vec3_sub(b, a);
        const double l2 = vec3_norm_sq(ab);
        const double t = l2 > 0.0 ? std::clamp(vec3_dot(vec3_sub(rP, a), ab) / l2, 0.0, 1.0) : 0.0;
        const Vec3 q = vec3_add(a, vec3_scale(ab, t));
        const double d = vec3_norm(vec3_sub(rP, q));
        if (d < rDistance) {
            rDistance = d;
            best = q;
        }
    }
    return best;
}

// What the feature snap of a hexahedral boundary recovered.
struct OctmhFeatures {
    std::int64_t mCorners = 0, mCornersSnapped = 0, mChains = 0, mChainsMissed = 0, mMoved = 0;
};

// Snap the boundary nodes @p rBoundaryNodes (with the boundary quads' edges as
// their graph) onto the feature corners and chains of the skin: a corner
// takes the nearest unclaimed node within 1.5 sizes, a chain the shortest
// path between its corners' nodes in a tube of 1.25 sizes around it. Snapped
// nodes are set to 1 (corner) or 2 (chain) in @p rFeature.
OctmhFeatures octmh_snap(std::vector<Vec3>& rPoints,
                         const std::vector<std::array<std::int64_t, 4>>& rQuads,
                         const FeatureGraph& rGraph, const std::vector<Vec3>& rSkin, double H,
                         std::vector<std::uint8_t>& rFeature) {
    OctmhFeatures rep;
    rep.mCorners = static_cast<std::int64_t>(rGraph.mCorners.size());
    rep.mChains = static_cast<std::int64_t>(rGraph.mChains.size());
    std::vector<std::vector<std::int64_t>> adj(rPoints.size());
    std::vector<std::int64_t> nodes;
    for (const auto& q : rQuads)
        for (std::size_t e = 0; e < 4; ++e) {
            adj[static_cast<std::size_t>(q[e])].push_back(q[(e + 1) % 4]);
            adj[static_cast<std::size_t>(q[(e + 1) % 4])].push_back(q[e]);
            nodes.push_back(q[e]);
        }
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    for (auto& a : adj) {
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
    }
    // Corners.
    std::vector<std::int64_t> corner_node(rSkin.size(), -1);
    for (std::int64_t c : rGraph.mCorners) {
        const Vec3& x = rSkin[static_cast<std::size_t>(c)];
        std::int64_t best = -1;
        double best_d = 1.5 * H;
        for (std::int64_t v : nodes) {
            if (rFeature[static_cast<std::size_t>(v)])
                continue;
            const double d = vec3_norm(vec3_sub(rPoints[static_cast<std::size_t>(v)], x));
            if (d < best_d) {
                best_d = d;
                best = v;
            }
        }
        if (best < 0)
            continue;
        rPoints[static_cast<std::size_t>(best)] = x;
        rFeature[static_cast<std::size_t>(best)] = 1;
        corner_node[static_cast<std::size_t>(c)] = best;
        ++rep.mCornersSnapped;
        ++rep.mMoved;
    }
    // Chains: Dijkstra from one end's node to the other's inside the tube.
    for (const auto& chain : rGraph.mChains) {
        const std::int64_t s = corner_node[static_cast<std::size_t>(chain.front())];
        const std::int64_t t = corner_node[static_cast<std::size_t>(chain.back())];
        if (s < 0 || t < 0 || s == t) {
            ++rep.mChainsMissed;
            continue;
        }
        std::vector<double> dist(rPoints.size(), std::numeric_limits<double>::infinity());
        std::vector<double> tube(rPoints.size(), -1.0);
        std::vector<std::int64_t> prev(rPoints.size(), -1);
        auto off = [&](std::int64_t V) {
            double& d = tube[static_cast<std::size_t>(V)];
            if (d < 0.0)
                octmh_on_chain(rSkin, chain, rPoints[static_cast<std::size_t>(V)], d);
            return d;
        };
        using Item = std::pair<double, std::int64_t>;
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
        dist[static_cast<std::size_t>(s)] = 0.0;
        heap.emplace(0.0, s);
        while (!heap.empty()) {
            const auto [d, u] = heap.top();
            heap.pop();
            if (d > dist[static_cast<std::size_t>(u)])
                continue;
            if (u == t)
                break;
            for (std::int64_t v : adj[static_cast<std::size_t>(u)]) {
                if (v != t && rFeature[static_cast<std::size_t>(v)])
                    continue;
                if (off(v) > 1.25 * H)
                    continue;
                const double len = vec3_norm(vec3_sub(rPoints[static_cast<std::size_t>(v)],
                                                      rPoints[static_cast<std::size_t>(u)]));
                const double nd = d + len * (1.0 + 16.0 * (off(u) + off(v)) / H);
                if (nd < dist[static_cast<std::size_t>(v)]) {
                    dist[static_cast<std::size_t>(v)] = nd;
                    prev[static_cast<std::size_t>(v)] = u;
                    heap.emplace(nd, v);
                }
            }
        }
        if (prev[static_cast<std::size_t>(t)] < 0) {
            ++rep.mChainsMissed;
            continue;
        }
        for (std::int64_t v = prev[static_cast<std::size_t>(t)]; v != s;
             v = prev[static_cast<std::size_t>(v)]) {
            double d = 0.0;
            rPoints[static_cast<std::size_t>(v)] =
                octmh_on_chain(rSkin, chain, rPoints[static_cast<std::size_t>(v)], d);
            rFeature[static_cast<std::size_t>(v)] = 2;
            ++rep.mMoved;
        }
    }
    return rep;
}

}  // namespace

void octree_mesh_hexahedra(OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR) {
    const OctreeFrame& frame = rS.mFrame;
    octree_mesh_report_tree(rS.mOctree.mTree, rR);
    // The grid: every cell at the skin's level.
    const double hs = rO.mSurfaceSize > 0.0 ? rO.mSurfaceSize : rS.mField.mSkinMinSize;
    const std::int32_t level = frame.LevelFor(std::max(hs, frame.mMinSize));
    const double h = frame.mRootSize / static_cast<double>(std::int64_t{1} << level);
    std::array<std::int64_t, 3> n{};
    double cells = 1.0;
    for (std::size_t k = 0; k < 3; ++k) {
        n[k] = frame.mRoots[k] << level;
        cells *= static_cast<double>(n[k]);
    }
    if (rO.mMaxCells > 0 && cells > static_cast<double>(rO.mMaxCells))
        throw std::invalid_argument(
            std::string(kOctreeMeshPrefix) + "the hexahedral grid would have " +
            std::to_string(static_cast<std::int64_t>(cells)) + " cells, more than max_cells (" +
            std::to_string(rO.mMaxCells) + "); raise surface_size or max_cells");
    rR.mSurfaceLevel = level;
    const auto num = static_cast<std::size_t>(cells);
    std::vector<Vec3> centres(num);
    for (std::int64_t k = 0; k < n[2]; ++k)
        for (std::int64_t j = 0; j < n[1]; ++j)
            for (std::int64_t i = 0; i < n[0]; ++i)
                centres[static_cast<std::size_t>((k * n[1] + j) * n[0] + i)] = {
                    frame.mOrigin[0] + (static_cast<double>(i) + 0.5) * h,
                    frame.mOrigin[1] + (static_cast<double>(j) + 0.5) * h,
                    frame.mOrigin[2] + (static_cast<double>(k) + 0.5) * h};
    const std::vector<double> w = octree_mesh_winding(rS, centres, rO, rR);
    std::vector<std::uint8_t> kept(num, 0);
    for (std::size_t c = 0; c < num; ++c)
        kept[c] = w[c] > rO.mWindingThreshold;
    if (std::none_of(kept.begin(), kept.end(), [](std::uint8_t k) { return k != 0; }))
        throw std::invalid_argument(std::string(kOctreeMeshPrefix) + "no grid cell of edge " +
                                    std::to_string(h) +
                                    " has its centre inside the skin; lower surface_size");
    const std::int64_t added = hex_grid_make_manifold(n, kept, w);
    if (added > 0)
        log::debug("{}{} cells added for a manifold boundary", kOctreeMeshPrefix, added);

    // The hexahedra of the kept cells, on the grid vertices they use.
    const std::int64_t vx = n[0] + 1, vy = n[1] + 1;
    auto vertex = [&](std::int64_t I, std::int64_t J, std::int64_t K) {
        return (K * vy + J) * vx + I;
    };
    std::vector<std::int64_t> remap(static_cast<std::size_t>(vx * vy * (n[2] + 1)), -1);
    HexMesh hm;
    std::vector<std::int64_t> hex_cell;
    for (std::int64_t k = 0; k < n[2]; ++k)
        for (std::int64_t j = 0; j < n[1]; ++j)
            for (std::int64_t i = 0; i < n[0]; ++i)
                if (kept[static_cast<std::size_t>((k * n[1] + j) * n[0] + i)]) {
                    hex_cell.push_back((k * n[1] + j) * n[0] + i);
                    for (std::int64_t dz = 0; dz < 2; ++dz)
                        for (std::int64_t dy = 0; dy < 2; ++dy)
                            for (std::int64_t dx = 0; dx < 2; ++dx)
                                remap[static_cast<std::size_t>(vertex(i + dx, j + dy, k + dz))] = 0;
                }
    for (std::int64_t k = 0; k <= n[2]; ++k)
        for (std::int64_t j = 0; j <= n[1]; ++j)
            for (std::int64_t i = 0; i < vx; ++i) {
                std::int64_t& r = remap[static_cast<std::size_t>(vertex(i, j, k))];
                if (r < 0)
                    continue;
                r = static_cast<std::int64_t>(hm.mPoints.size());
                hm.mPoints.push_back({frame.mOrigin[0] + static_cast<double>(i) * h,
                                      frame.mOrigin[1] + static_cast<double>(j) * h,
                                      frame.mOrigin[2] + static_cast<double>(k) * h});
            }
    for (std::int64_t c : hex_cell) {
        const std::int64_t i = c % n[0], j = (c / n[0]) % n[1], k = c / (n[0] * n[1]);
        auto at = [&](std::int64_t DX, std::int64_t DY, std::int64_t DZ) {
            return remap[static_cast<std::size_t>(vertex(i + DX, j + DY, k + DZ))];
        };
        hm.mHexes.push_back({at(0, 0, 0), at(1, 0, 0), at(1, 1, 0), at(0, 1, 0), at(0, 0, 1),
                             at(1, 0, 1), at(1, 1, 1), at(0, 1, 1)});
    }
    const std::size_t num_grid = hm.mHexes.size();

    // The boundary, pillowed or not, onto the skin.
    const HexBoundary grid_boundary = hex_boundary(hm);
    std::vector<std::int64_t> boundary_nodes;
    for (const auto& q : grid_boundary.mQuads)
        boundary_nodes.insert(boundary_nodes.end(), q.begin(), q.end());
    std::sort(boundary_nodes.begin(), boundary_nodes.end());
    boundary_nodes.erase(std::unique(boundary_nodes.begin(), boundary_nodes.end()),
                         boundary_nodes.end());
    std::vector<Vec3> at_nodes;
    for (std::int64_t v : boundary_nodes)
        at_nodes.push_back(hm.mPoints[static_cast<std::size_t>(v)]);
    const std::vector<SurfaceProjection> proj = query_surface_projections(rS.mQuery, at_nodes);
    std::vector<std::int64_t> copy;
    if (rO.mPillow)
        copy = hex_pillow(hm, grid_boundary);
    std::vector<HexNodeMode> mode(hm.mPoints.size(), HexNodeMode::Free);
    for (std::size_t b = 0; b < boundary_nodes.size(); ++b) {
        const auto v = static_cast<std::size_t>(boundary_nodes[b]);
        const Vec3 q = proj[b].mPoint;
        if (rO.mPillow) {
            // The copy on the skin, the grid's node half a cell inside it.
            const Vec3 nrm = vec3_scale(vec3_normalize(proj[b].mNormal), rS.mWindingSign);
            hm.mPoints[static_cast<std::size_t>(copy[v])] = q;
            mode[static_cast<std::size_t>(copy[v])] = HexNodeMode::Slide;
            hm.mPoints[v] = vec3_sub(q, vec3_scale(nrm, 0.5 * h));
        } else {
            hm.mPoints[v] = q;
            mode[v] = HexNodeMode::Slide;
        }
    }
    const HexBoundary boundary = rO.mPillow ? hex_boundary(hm) : grid_boundary;

    // Features, then untangling with the skin nodes sliding along it.
    std::vector<std::uint8_t> feature(hm.mPoints.size(), 0);
    if (rO.mPreserveFeatures && !rS.mFeatures.empty()) {
        std::vector<std::array<std::int64_t, 2>> edges;
        for (const OctreeFeatureEdge& e : rS.mFeatures)
            edges.push_back({e.mA, e.mB});
        const FeatureGraph graph = feature_graph(rS.mSkin.mPoints, edges, rO.mCornerAngle);
        const OctmhFeatures f =
            octmh_snap(hm.mPoints, boundary.mQuads, graph, rS.mSkin.mPoints, h, feature);
        rR.mNumCorners = f.mCorners;
        rR.mNumCornersSnapped = f.mCornersSnapped;
        rR.mNumFeatureChains = f.mChains;
        rR.mNumFeatureChainsMissed = f.mChainsMissed;
        rR.mNumFeatureMoves = f.mMoved;
        for (std::size_t p = 0; p < feature.size(); ++p)
            if (feature[p])
                mode[p] = HexNodeMode::Fixed;
        if (f.mCornersSnapped < f.mCorners)
            octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) +
                                     std::to_string(f.mCorners - f.mCornersSnapped) + " of " +
                                     std::to_string(f.mCorners) + " feature corners not recovered");
        if (f.mChainsMissed > 0)
            octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(f.mChainsMissed) +
                                     " of " + std::to_string(f.mChains) +
                                     " feature chains not recovered");
    }
    HexUntangleOptions uo;
    uo.mIterations = rO.mUntangleIterations;
    const HexUntangleReport ur = hex_untangle(
        hm.mPoints, hm.mHexes, mode,
        [&](const Vec3& rP) { return query_closest_points(rS.mQuery, {rP}).front().mPoint; }, uo);
    rR.mNumOptimizedMoves = ur.mMoves;
    double worst = 1.0;
    for (const auto& hex : hm.mHexes)
        worst = std::min(worst, hex_min_scaled_jacobian(hm.mPoints, hex));
    if (worst <= 0.0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) +
                                 "untangling left a hexahedron with a non-positive corner "
                                 "Jacobian; raise untangle_iterations or lower surface_size");
    else if (worst < 0.2)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) +
                                 "the smallest corner scaled Jacobian is " + std::to_string(worst));

    // The mesh, its data and regions.
    const std::size_t np = hm.mPoints.size();
    NDArray points = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
    for (std::size_t p = 0; p < np; ++p)
        for (std::size_t k = 0; k < 3; ++k)
            points.As<double>()[p * 3 + k] = hm.mPoints[p][k];
    NDArray conn = NDArray::Uninit(DType::Int64, {hm.mHexes.size(), std::size_t{8}});
    NDArray lv = NDArray::Uninit(DType::Int64, {hm.mHexes.size()});
    NDArray pillow = NDArray::Uninit(DType::Int64, {hm.mHexes.size()});
    for (std::size_t c = 0; c < hm.mHexes.size(); ++c) {
        for (std::size_t k = 0; k < 8; ++k)
            conn.As<std::int64_t>()[c * 8 + k] = hm.mHexes[c][k];
        lv.As<std::int64_t>()[c] = level;
        pillow.As<std::int64_t>()[c] = c >= num_grid ? 1 : 0;
    }
    Mesh& m = rR.mMesh;
    m.AssignPoints(std::move(points));
    m.AddCellBlock("hexahedron", std::move(conn));
    {
        std::vector<NDArray> a, b;
        a.push_back(std::move(lv));
        b.push_back(std::move(pillow));
        m.AddCellData(kOctreeLevelName, std::move(a));
        m.AddCellData(kOctreeMeshPillowName, std::move(b));
    }
    std::vector<std::int64_t> all(hm.mHexes.size());
    for (std::size_t c = 0; c < all.size(); ++c)
        all[c] = static_cast<std::int64_t>(c);
    rR.mNumInside = static_cast<std::int64_t>(all.size());
    m.AddRegion(Region(kOctreeMeshInsideRegion, RegionKind::Cell, 3, -1, octmh_ids(all, 1)));
    std::vector<std::array<std::int64_t, 2>> sides = boundary.mSides;
    std::sort(sides.begin(), sides.end());
    std::vector<std::int64_t> flat;
    for (const auto& s : sides) {
        flat.push_back(s[0]);
        flat.push_back(s[1]);
    }
    m.AddRegion(Region(kOctreeMeshBoundaryRegion, RegionKind::Side, 2, -1, octmh_ids(flat, 2)));
    if (std::any_of(feature.begin(), feature.end(), [](std::uint8_t f) { return f != 0; })) {
        std::vector<std::int64_t> corners, chains;
        for (std::size_t p = 0; p < np; ++p)
            if (feature[p]) {
                chains.push_back(static_cast<std::int64_t>(p));
                if (feature[p] == 1)
                    corners.push_back(static_cast<std::int64_t>(p));
            }
        m.AddRegion(
            Region(kOctreeMeshCornersRegion, RegionKind::Point, 0, -1, octmh_ids(corners, 1)));
        m.AddRegion(
            Region(kOctreeMeshFeatureEdgesRegion, RegionKind::Point, 1, -1, octmh_ids(chains, 1)));
    }
    rR.mNumBoundaryEdgesOpen = boundary.mNumOpenEdges;
    rR.mNumNonManifoldBoundaryEdges = boundary.mNumNonManifoldEdges;
    std::vector<std::array<std::int64_t, 3>> tris;
    tris.reserve(boundary.mQuads.size() * 2);
    for (const auto& q : boundary.mQuads) {
        tris.push_back({q[0], q[1], q[2]});
        tris.push_back({q[0], q[2], q[3]});
    }
    if (!tris.empty())
        rR.mHausdorffToSkin = octree_mesh_hausdorff_to_skin(rS, hm.mPoints, tris);
    if (rO.mAttachDistance) {
        std::vector<double> dist;
        std::vector<Vec3> closest;
        octree_mesh_signed_distance(rS, hm.mPoints, rO, rR, dist, closest);
        NDArray d = NDArray::Uninit(DType::Float64, {np});
        std::copy(dist.begin(), dist.end(), d.As<double>());
        m.AddPointData(kOctreeMeshDistanceName, std::move(d));
    }
}

}  // namespace detail
}  // namespace meshioplusplus
