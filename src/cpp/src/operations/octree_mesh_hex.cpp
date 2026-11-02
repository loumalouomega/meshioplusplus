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
// the skin (or, graded, the dual of the octree's transition schemes), made
// manifold, pillowed, fitted to the skin and its features, and untangled. Anonymous-namespace
// helpers are prefixed `octmh_` (the amalgamation concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "../detail/feature_snap.hpp"
#include "../detail/hex_dual_schemes.hpp"
#include "../detail/hex_fit.hpp"
#include "../detail/tri_bvh.hpp"
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
    /// The boundary edges the traced chains run along, each (low, high), and
    /// the chain (an index into the feature graph's chains) of each.
    std::vector<std::array<std::int64_t, 2>> mPathEdges;
    std::vector<std::int64_t> mPathChain;
    /// Per point, the chain a snapped chain node lies on, and the skin point a
    /// snapped corner node lies on (-1: neither).
    std::vector<std::int64_t> mNodeChain, mNodeCorner;
};

// The point at arc length @p At along the polyline through the skin points
// @p rChain (clamped to its ends).
Vec3 octmh_along_chain(const std::vector<Vec3>& rSkin, const std::vector<std::int64_t>& rChain,
                       double At) {
    for (std::size_t s = 0; s + 1 < rChain.size(); ++s) {
        const Vec3& a = rSkin[static_cast<std::size_t>(rChain[s])];
        const Vec3& b = rSkin[static_cast<std::size_t>(rChain[s + 1])];
        const double l = vec3_norm(vec3_sub(b, a));
        if (At <= l || s + 2 == rChain.size())
            return l > 0.0 ? vec3_add(a, vec3_scale(vec3_sub(b, a), std::clamp(At / l, 0.0, 1.0)))
                           : a;
        At -= l;
    }
    return rSkin[static_cast<std::size_t>(rChain.front())];
}

// Snap the boundary nodes @p rBoundaryNodes (with the boundary quads' edges as
// their graph) onto the feature corners and chains of the skin: a corner
// takes the nearest unclaimed node within 1.5 of its size @p rSize, a chain the
// shortest path between its corners' nodes in a tube of 1.25 sizes around it. Snapped
// nodes are set to 1 (corner) or 2 (chain) in @p rFeature. A path's inner nodes
// go to their nearest chain points, or, with @p ArcLength, to the chain points
// at the same fraction of its length as they are of the path's, which keeps
// them in order and apart (the nearest point of a node beside a corner can be
// the corner itself).
OctmhFeatures octmh_snap(std::vector<Vec3>& rPoints,
                         const std::vector<std::array<std::int64_t, 4>>& rQuads,
                         const FeatureGraph& rGraph, const std::vector<Vec3>& rSkin,
                         const std::vector<double>& rSize, std::vector<std::uint8_t>& rFeature,
                         bool ArcLength) {
    OctmhFeatures rep;
    rep.mNodeChain.assign(rPoints.size(), -1);
    rep.mNodeCorner.assign(rPoints.size(), -1);
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
        double best_d = std::numeric_limits<double>::infinity();
        for (std::int64_t v : nodes) {
            if (rFeature[static_cast<std::size_t>(v)])
                continue;
            const double d = vec3_norm(vec3_sub(rPoints[static_cast<std::size_t>(v)], x));
            if (d < best_d && d < 1.5 * rSize[static_cast<std::size_t>(v)]) {
                best_d = d;
                best = v;
            }
        }
        if (best < 0)
            continue;
        rPoints[static_cast<std::size_t>(best)] = x;
        rFeature[static_cast<std::size_t>(best)] = 1;
        rep.mNodeCorner[static_cast<std::size_t>(best)] = c;
        corner_node[static_cast<std::size_t>(c)] = best;
        ++rep.mCornersSnapped;
        ++rep.mMoved;
    }
    // Chains: Dijkstra from one end's node to the other's inside the tube.
    for (std::size_t ci = 0; ci < rGraph.mChains.size(); ++ci) {
        const auto& chain = rGraph.mChains[ci];
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
                if (off(v) > 1.25 * rSize[static_cast<std::size_t>(v)])
                    continue;
                const double len = vec3_norm(vec3_sub(rPoints[static_cast<std::size_t>(v)],
                                                      rPoints[static_cast<std::size_t>(u)]));
                const double h =
                    0.5 * (rSize[static_cast<std::size_t>(u)] + rSize[static_cast<std::size_t>(v)]);
                const double nd = d + len * (1.0 + 16.0 * (off(u) + off(v)) / h);
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
        for (std::int64_t v = t; v != s; v = prev[static_cast<std::size_t>(v)]) {
            const std::int64_t u = prev[static_cast<std::size_t>(v)];
            rep.mPathEdges.push_back({std::min(u, v), std::max(u, v)});
            rep.mPathChain.push_back(static_cast<std::int64_t>(ci));
        }
        if (ArcLength) {
            std::vector<std::int64_t> path;
            for (std::int64_t v = t; v >= 0 && v != s; v = prev[static_cast<std::size_t>(v)])
                path.push_back(v);
            path.push_back(s);
            std::reverse(path.begin(), path.end());
            std::vector<double> at(path.size(), 0.0);
            for (std::size_t i = 1; i < path.size(); ++i)
                at[i] =
                    at[i - 1] + vec3_norm(vec3_sub(rPoints[static_cast<std::size_t>(path[i])],
                                                   rPoints[static_cast<std::size_t>(path[i - 1])]));
            double length = 0.0;
            for (std::size_t i = 0; i + 1 < chain.size(); ++i)
                length += vec3_norm(vec3_sub(rSkin[static_cast<std::size_t>(chain[i + 1])],
                                             rSkin[static_cast<std::size_t>(chain[i])]));
            for (std::size_t i = 1; i + 1 < path.size(); ++i) {
                const auto v = static_cast<std::size_t>(path[i]);
                rPoints[v] = octmh_along_chain(rSkin, chain, length * at[i] / at.back());
                rFeature[v] = 2;
                rep.mNodeChain[v] = static_cast<std::int64_t>(ci);
                ++rep.mMoved;
            }
            continue;
        }
        for (std::int64_t v = prev[static_cast<std::size_t>(t)]; v != s;
             v = prev[static_cast<std::size_t>(v)]) {
            double d = 0.0;
            rPoints[static_cast<std::size_t>(v)] =
                octmh_on_chain(rSkin, chain, rPoints[static_cast<std::size_t>(v)], d);
            rFeature[static_cast<std::size_t>(v)] = 2;
            rep.mNodeChain[static_cast<std::size_t>(v)] = static_cast<std::int64_t>(ci);
            ++rep.mMoved;
        }
    }
    return rep;
}

// The uniform grid: every cell at the skin's level, the kept ones made
// manifold. @p rLevel gets each hexahedron's level, @p rSize each node's size.
void octmh_uniform(OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR,
                   HexMesh& hm, std::vector<std::int64_t>& rLevel, std::vector<double>& rSize) {
    const OctreeFrame& frame = rS.mFrame;
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
    rLevel.assign(hm.mHexes.size(), level);
    rSize.assign(hm.mPoints.size(), h);
}

// The graded hexahedra: the all-hexahedral dual of the octree, paired and
// with its transitions replaced by Livesu, Pitzalis and Cherchi's schemes
// (detail/hex_dual_schemes.hpp). A dual hexahedron is kept when the primal
// vertex it surrounds is inside; @p rSize is each node's mean edge.
void octmh_graded(OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR,
                  HexMesh& hm, std::vector<std::int64_t>& rLevel, std::vector<double>& rSize) {
    const OctreeFrame& frame = rS.mFrame;
    LinearOctree tree = rS.mOctree.mTree;
    linear_octree_pair(tree, OctreeAdjacency::Full, rO.mMaxCells);
    const HexDualMesh dual = hex_dual_mesh(tree);
    const double unit = frame.mRootSize / static_cast<double>(std::int64_t{1} << tree.mMaxLevel) /
                        static_cast<double>(kHexDualUnit);
    auto world = [&](double X, double Y, double Z) {
        return Vec3{frame.mOrigin[0] + X * unit, frame.mOrigin[1] + Y * unit,
                    frame.mOrigin[2] + Z * unit};
    };
    std::vector<Vec3> centres(dual.mHexes.size());
    for (std::size_t c = 0; c < centres.size(); ++c)
        centres[c] =
            world(static_cast<double>(dual.mCentre[c][0]), static_cast<double>(dual.mCentre[c][1]),
                  static_cast<double>(dual.mCentre[c][2]));
    const std::vector<double> w = octree_mesh_winding(rS, centres, rO, rR);
    std::vector<std::uint8_t> kept(centres.size(), 0);
    for (std::size_t c = 0; c < kept.size(); ++c)
        kept[c] = w[c] > rO.mWindingThreshold;
    // A centre is a primal vertex, on the lattice's planes, and a skin may lie on
    // them too (a cap on a lattice plane): a hexahedron kept by a hair there
    // would be pressed flat by the projection and the pillow. Keep only centres
    // more than a tenth of their cell inside.
    {
        std::vector<std::size_t> near;
        std::vector<Vec3> at;
        for (std::size_t c = 0; c < kept.size(); ++c)
            if (kept[c]) {
                near.push_back(c);
                at.push_back(centres[c]);
            }
        const std::vector<SurfaceProjection> proj = query_surface_projections(rS.mQuery, at);
        for (std::size_t i = 0; i < near.size(); ++i) {
            const std::size_t c = near[i];
            const double cell =
                frame.mRootSize / static_cast<double>(std::int64_t{1} << dual.mLevel[c]);
            if (vec3_norm(vec3_sub(proj[i].mPoint, centres[c])) < 0.1 * cell)
                kept[c] = 0;
        }
    }
    if (std::none_of(kept.begin(), kept.end(), [](std::uint8_t k) { return k != 0; }))
        throw std::invalid_argument(std::string(kOctreeMeshPrefix) +
                                    "no graded hexahedron has its centre inside the skin; "
                                    "lower min_size");
    HexMesh all;
    all.mPoints.reserve(dual.mPoints.size());
    for (const Vec3& p : dual.mPoints)
        all.mPoints.push_back(world(p[0], p[1], p[2]));
    all.mHexes = dual.mHexes;
    const std::int64_t added = hex_make_manifold(all, kept, w);
    if (added > 0)
        log::debug("{}{} cells added for a manifold boundary", kOctreeMeshPrefix, added);
    std::int64_t level = 0;
    for (std::int32_t l : dual.mLevel)
        level = std::max<std::int64_t>(level, l);
    rR.mSurfaceLevel = level;
    std::vector<std::int64_t> remap(all.mPoints.size(), -1);
    for (std::size_t c = 0; c < kept.size(); ++c) {
        if (!kept[c])
            continue;
        std::array<std::int64_t, 8> hex{};
        for (std::size_t k = 0; k < 8; ++k) {
            std::int64_t& r = remap[static_cast<std::size_t>(all.mHexes[c][k])];
            if (r < 0) {
                r = static_cast<std::int64_t>(hm.mPoints.size());
                hm.mPoints.push_back(all.mPoints[static_cast<std::size_t>(all.mHexes[c][k])]);
            }
            hex[k] = r;
        }
        hm.mHexes.push_back(hex);
        rLevel.push_back(dual.mLevel[c]);
    }
    // Each node's size: the mean of its edges.
    rSize.assign(hm.mPoints.size(), 0.0);
    std::vector<std::int64_t> count(hm.mPoints.size(), 0);
    static constexpr int kEdges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                          {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (const auto& hex : hm.mHexes)
        for (const auto& e : kEdges) {
            const auto a = static_cast<std::size_t>(hex[static_cast<std::size_t>(e[0])]);
            const auto b = static_cast<std::size_t>(hex[static_cast<std::size_t>(e[1])]);
            const double len = vec3_norm(vec3_sub(hm.mPoints[a], hm.mPoints[b]));
            rSize[a] += len;
            rSize[b] += len;
            ++count[a];
            ++count[b];
        }
    for (std::size_t v = 0; v < rSize.size(); ++v)
        rSize[v] /= static_cast<double>(std::max<std::int64_t>(count[v], 1));
    std::int64_t schemes = 0;
    for (std::int64_t n : dual.mNumSchemes)
        schemes += n;
    log::debug("{}{} transition vertices, {} cells replaced by schemes, {} dual hexahedra",
               kOctreeMeshPrefix, dual.mNumTransitionVertices, schemes, dual.mHexes.size());
}

// The smoothing rounds of the boundary normals the robust fit projects along.
constexpr int kOctmhNormalSmoothing = 1;

// Below this smallest corner scaled Jacobian the robust fit is compared with
// the body-fitted one, and the better is kept.
constexpr double kOctmhFallbackBelow = 0.2;

// Log the quality of the hexahedra after a stage (debug level only).
void octmh_log_stage(const char* pStage, const std::vector<Vec3>& rPoints,
                     const std::vector<std::array<std::int64_t, 8>>& rHexes) {
    if (!log::enabled(log::Level::Debug))
        return;
    double worst = 1.0;
    std::int64_t poor = 0, bad = 0;
    for (const auto& hex : rHexes) {
        const double w = hex_min_scaled_jacobian(rPoints, hex);
        worst = std::min(worst, w);
        poor += w < 0.2;
        bad += w <= 0.0;
    }
    log::debug(
        "{}after {}: smallest corner scaled Jacobian {}, {} hexahedra below 0.2, {} at or "
        "below 0",
        kOctreeMeshPrefix, pStage, worst, poor, bad);
}

// How the boundary is placed on the skin.
enum class OctmhPlacement {
    /// Along each node's boundary normal, the inward offset held to the local
    /// thickness, sliding along the normal too (closest point where the
    /// normal misses).
    Robust,
    /// The closest point, as in 18.0.0.
    BodyFitted,
};

// What the fit works on: the pillowed (or bare) mesh at the grid's positions.
struct OctmhFitInput {
    const OctreeMeshSetup& rS;
    const OctreeMeshOptions& rO;
    const HexMesh& rMesh;
    const std::vector<std::int64_t>& rBoundaryNodes;  // the grid's boundary nodes
    const std::vector<std::int64_t>& rCopy;           // their pillow copies (pillow only)
    const std::vector<double>& rSize;
    const std::vector<Vec3>& rNormals;                       // the grid boundary's unit normals
    const std::vector<std::array<std::int64_t, 4>>& rQuads;  // the fitted mesh's boundary
    const FeatureGraph* pGraph;                              // null: no feature snap
};

// A fitted mesh's positions, features and quality.
struct OctmhFit {
    std::vector<Vec3> mPoints;
    std::vector<std::uint8_t> mFeature;
    OctmhFeatures mFeatures;
    std::int64_t mMoves = 0;
    double mWorst = 1.0;
    double mRobustWorst = 1.0;  // set when the body-fitted fit was kept
    bool mFellBack = false;
    std::int64_t mReleased = 0;  // feature nodes released to untangle
};

// Where the line through @p rX along the unit @p rDir first meets the skin
// from outside (a triangle facing along @p rDir), nearest @p rX within
// @p Reach either way; false when it misses.
bool octmh_ray_foot(const OctreeMeshSetup& rS, const Vec3& rX, const Vec3& rDir, double Reach,
                    Vec3& rFoot) {
    const double sign = rS.mWindingSign;
    auto facing = [&](std::int64_t, const Vec3& rN) { return sign * vec3_dot(rN, rDir) > 0.0; };
    const TriRayHit ahead = tri_bvh_ray_filtered(rS.mBvh, rX, rDir, -1e-12 * Reach, Reach, facing);
    const TriRayHit behind =
        tri_bvh_ray_filtered(rS.mBvh, rX, vec3_scale(rDir, -1.0), 0.0, Reach, facing);
    if (ahead.mTriangle < 0 && behind.mTriangle < 0)
        return false;
    const bool use_ahead =
        ahead.mTriangle >= 0 && (behind.mTriangle < 0 || std::abs(ahead.mT) <= behind.mT);
    rFoot = use_ahead ? vec3_add(rX, vec3_scale(rDir, ahead.mT))
                      : vec3_sub(rX, vec3_scale(rDir, behind.mT));
    return true;
}

// The distance from the skin point @p rQ inward along -@p rOut (unit,
// outward) to the far wall of the skin, or @p Reach when it is farther.
double octmh_thickness(const OctreeMeshSetup& rS, const Vec3& rQ, const Vec3& rOut, double Reach) {
    const Vec3 in = vec3_scale(rOut, -1.0);
    const double sign = rS.mWindingSign;
    auto far_wall = [&](std::int64_t, const Vec3& rN) { return sign * vec3_dot(rN, in) > 0.0; };
    const TriRayHit hit = tri_bvh_ray_filtered(rS.mBvh, rQ, in, 1e-9 * Reach, Reach, far_wall);
    return hit.mTriangle < 0 ? Reach : hit.mT;
}

// Put the grid node @p V of the pillow copy @p Copy half a size @p Size inside
// the skin along -@p rOut, but no more than 0.45 of the way to the far wall,
// so a part thinner than a cell does not fold. @return whether it was held.
bool octmh_inset(const OctreeMeshSetup& rS, std::vector<Vec3>& rPoints, std::size_t V,
                 std::size_t Copy, const Vec3& rOut, double Size) {
    const Vec3 q = rPoints[Copy];
    const double t = octmh_thickness(rS, q, rOut, 1.2 * Size);
    const double d = std::min(0.5 * Size, 0.45 * t);
    rPoints[V] = vec3_sub(q, vec3_scale(rOut, d));
    return d < 0.5 * Size;
}

// The robust foot of a point @p rX whose closest skin point is @p Distance
// away: where the line through it along the unit normal @p rNormal meets the
// skin, when that is no more than half a size @p Size farther than the
// closest point; false otherwise (the caller keeps the closest point).
bool octmh_robust_foot(const OctreeMeshSetup& rS, const Vec3& rX, const Vec3& rNormal,
                       double Distance, double Size, Vec3& rFoot) {
    if (!(vec3_norm(rNormal) > 0.0))
        return false;
    Vec3 f{};
    if (!octmh_ray_foot(rS, rX, rNormal, Distance + 1.5 * Size, f))
        return false;
    if (vec3_norm(vec3_sub(f, rX)) > Distance + 0.5 * Size)
        return false;
    rFoot = f;
    return true;
}

// A chain traced on a staircase turns inside a quadrilateral (along two of
// its edges at one corner), and snapping the three nodes onto the chain's
// polyline puts them on one line, which flattens the quadrilateral and leaves
// its hexahedron a corner of zero scaled Jacobian that no move of the pinned
// nodes can raise. Release the middle node of every such turn from the chain:
// it slides along the skin instead, tethered to the chain (`octmh_tether`).
// @return per point, the chain a released node was on, or -1.
std::vector<std::int64_t> octmh_release_turns(
    const std::vector<std::array<std::int64_t, 4>>& rQuads, const OctmhFeatures& rFeatures,
    std::vector<std::uint8_t>& rFeature) {
    std::vector<std::pair<std::array<std::int64_t, 2>, std::int64_t>> edges;
    for (std::size_t e = 0; e < rFeatures.mPathEdges.size(); ++e)
        edges.emplace_back(rFeatures.mPathEdges[e], rFeatures.mPathChain[e]);
    std::sort(edges.begin(), edges.end());
    // The chain of a path edge, or -1.
    auto chain_of = [&](std::int64_t A, std::int64_t B) -> std::int64_t {
        const std::array<std::int64_t, 2> e{std::min(A, B), std::max(A, B)};
        const auto it =
            std::lower_bound(edges.begin(), edges.end(), e,
                             [](const auto& rItem, const std::array<std::int64_t, 2>& rKey) {
                                 return rItem.first < rKey;
                             });
        return it != edges.end() && it->first == e ? it->second : -1;
    };
    std::vector<std::int64_t> released(rFeature.size(), -1);
    std::int64_t count = 0;
    for (const auto& q : rQuads)
        for (std::size_t k = 0; k < 4; ++k) {
            const std::int64_t o = q[k];
            if (rFeature[static_cast<std::size_t>(o)] != 2 &&
                released[static_cast<std::size_t>(o)] < 0)
                continue;
            const std::int64_t c = chain_of(q[(k + 3) % 4], o);
            if (c >= 0 && chain_of(o, q[(k + 1) % 4]) == c &&
                released[static_cast<std::size_t>(o)] < 0) {
                released[static_cast<std::size_t>(o)] = c;
                rFeature[static_cast<std::size_t>(o)] = 0;
                ++count;
            }
        }
    log::debug("{}{} chain nodes released where their chain turns inside a quadrilateral",
               kOctreeMeshPrefix, count);
    return released;
}

// How far (in sizes) a released chain node may slide from its chain.
constexpr double kOctmhTether = 0.08;

// Hold the skin point @p rP within @p Reach of the chain @p rChain: pulled
// straight back towards its nearest chain point.
Vec3 octmh_tether(const std::vector<Vec3>& rSkin, const std::vector<std::int64_t>& rChain,
                  const Vec3& rP, double Reach) {
    double d = 0.0;
    const Vec3 on = octmh_on_chain(rSkin, rChain, rP, d);
    if (d <= Reach)
        return rP;
    return vec3_add(on, vec3_scale(vec3_sub(rP, on), Reach / d));
}

// The repair after untangling: when a hexahedron's smallest corner is below
// this, the feature nodes of every such hexahedron are released (tethered
// within `kOctmhTether` sizes of their chain or corner) and the mesh is
// untangled again for a fifth of the sweeps.
constexpr double kOctmhRepairBelow = 0.05;

// Place the boundary on the skin, snap the features and untangle.
OctmhFit octmh_fit(const OctmhFitInput& rIn, OctmhPlacement Placement) {
    const OctreeMeshSetup& rs = rIn.rS;
    const OctreeMeshOptions& ro = rIn.rO;
    const bool robust = Placement == OctmhPlacement::Robust;
    OctmhFit out;
    out.mPoints = rIn.rMesh.mPoints;
    std::vector<Vec3>& pts = out.mPoints;
    const std::vector<double>& size = rIn.rSize;
    std::vector<Vec3> at_nodes;
    at_nodes.reserve(rIn.rBoundaryNodes.size());
    for (std::int64_t v : rIn.rBoundaryNodes)
        at_nodes.push_back(pts[static_cast<std::size_t>(v)]);
    const std::vector<SurfaceProjection> proj = query_surface_projections(rs.mQuery, at_nodes);
    // Per point, the direction a sliding node is projected along (robust).
    std::vector<Vec3> slide(robust ? pts.size() : 0, Vec3{0.0, 0.0, 0.0});
    std::vector<HexNodeMode> mode(pts.size(), HexNodeMode::Free);
    // Per pillow copy, its outward direction (robust).
    std::vector<Vec3> outward(robust ? pts.size() : 0, Vec3{0.0, 0.0, 0.0});
    std::int64_t along_normal = 0, held = 0;
    for (std::size_t b = 0; b < rIn.rBoundaryNodes.size(); ++b) {
        const auto v = static_cast<std::size_t>(rIn.rBoundaryNodes[b]);
        Vec3 q = proj[b].mPoint;
        if (!robust) {
            if (ro.mPillow) {
                // The copy on the skin, the grid's node half its size inside it.
                const Vec3 nrm = vec3_scale(vec3_normalize(proj[b].mNormal), rs.mWindingSign);
                const auto c = static_cast<std::size_t>(rIn.rCopy[v]);
                pts[c] = q;
                mode[c] = HexNodeMode::Slide;
                pts[v] = vec3_sub(q, vec3_scale(nrm, 0.5 * size[v]));
            } else {
                pts[v] = q;
                mode[v] = HexNodeMode::Slide;
            }
            continue;
        }
        // Robust: along the boundary's own normal, which the closest point
        // ignores where two sides of the skin are about as near (a node by an
        // edge would otherwise put its copy on the side its face does not face).
        Vec3 out_dir = vec3_scale(vec3_normalize(proj[b].mNormal), rs.mWindingSign);
        const Vec3& nb = rIn.rNormals[v];
        Vec3 f{};
        if (octmh_robust_foot(rs, pts[v], nb, proj[b].mDistance, size[v], f)) {
            q = f;
            out_dir = nb;
            ++along_normal;
        }
        const std::size_t moved = ro.mPillow ? static_cast<std::size_t>(rIn.rCopy[v]) : v;
        slide[moved] = nb;
        pts[moved] = q;
        mode[moved] = HexNodeMode::Slide;
        if (ro.mPillow) {
            outward[moved] = out_dir;
            held += octmh_inset(rs, pts, v, moved, out_dir, size[v]);
        }
    }
    if (robust)
        log::debug("{}{} of {} boundary nodes placed along their normal, {} held by the thickness",
                   kOctreeMeshPrefix, along_normal, rIn.rBoundaryNodes.size(), held);
    octmh_log_stage(robust ? "robust placement" : "body-fitted placement", pts, rIn.rMesh.mHexes);
    out.mFeature.assign(pts.size(), 0);
    // Per point, what a released feature node is tethered to: its chain, or
    // its corner's skin point (-1: free to slide).
    std::vector<std::int64_t> tether_chain, tether_corner;
    if (rIn.pGraph != nullptr) {
        out.mFeatures =
            octmh_snap(pts, rIn.rQuads, *rIn.pGraph, rs.mSkin.mPoints, size, out.mFeature, robust);
        if (robust) {
            tether_chain = octmh_release_turns(rIn.rQuads, out.mFeatures, out.mFeature);
            tether_corner.assign(pts.size(), -1);
        }
        if (robust && ro.mPillow)
            // A snapped copy can move far (up to 1.5 sizes for a corner): its
            // grid node follows it, so the pillow cell keeps its height.
            for (const std::int64_t v : rIn.rBoundaryNodes) {
                const auto c = static_cast<std::size_t>(rIn.rCopy[static_cast<std::size_t>(v)]);
                if (out.mFeature[c])
                    octmh_inset(rs, pts, static_cast<std::size_t>(v), c, outward[c],
                                size[static_cast<std::size_t>(v)]);
            }
        for (std::size_t p = 0; p < out.mFeature.size(); ++p)
            if (out.mFeature[p])
                mode[p] = HexNodeMode::Fixed;
        octmh_log_stage("feature snap", pts, rIn.rMesh.mHexes);
    }
    HexUntangleOptions uo;
    uo.mIterations = ro.mUntangleIterations;
    uo.mRobust = robust;
    HexProjector project;
    if (robust)
        project = [&](std::size_t Node, const Vec3& rP) {
            const SurfaceProjection c = query_surface_projections(rs.mQuery, {rP}).front();
            Vec3 f = c.mPoint;
            octmh_robust_foot(rs, rP, slide[Node], c.mDistance, size[Node], f);
            if (tether_chain.empty())
                return f;
            Vec3 held = f;
            const double within = kOctmhTether * size[Node];
            if (tether_chain[Node] >= 0) {
                const auto& chain =
                    rIn.pGraph->mChains[static_cast<std::size_t>(tether_chain[Node])];
                held = octmh_tether(rs.mSkin.mPoints, chain, f, within);
            } else if (tether_corner[Node] >= 0) {
                held = octmh_tether(rs.mSkin.mPoints, {tether_corner[Node], tether_corner[Node]}, f,
                                    within);
            }
            return held != f ? query_closest_points(rs.mQuery, {held}).front().mPoint : f;
        };
    else
        project = [&](std::size_t, const Vec3& rP) {
            return query_closest_points(rs.mQuery, {rP}).front().mPoint;
        };
    const HexUntangleReport ur = hex_untangle(pts, rIn.rMesh.mHexes, mode, project, uo);
    out.mMoves = ur.mMoves;
    for (const auto& hex : rIn.rMesh.mHexes)
        out.mWorst = std::min(out.mWorst, hex_min_scaled_jacobian(pts, hex));
    // The repair: a feature node pinned where the boundary cannot follow its
    // feature (a crease shorter or tighter than a cell, as a scan has) can
    // hold a corner flat or inverted. Release the feature nodes of the poorest
    // hexahedra, tethered to their feature, and untangle again.
    if (robust && rIn.pGraph != nullptr && ro.mUntangleIterations >= 5 &&
        out.mWorst < kOctmhRepairBelow) {
        std::vector<std::size_t> freeing;
        for (const auto& hex : rIn.rMesh.mHexes)
            if (hex_min_scaled_jacobian(pts, hex) < kOctmhRepairBelow)
                for (std::int64_t n : hex)
                    if (out.mFeature[static_cast<std::size_t>(n)])
                        freeing.push_back(static_cast<std::size_t>(n));
        std::sort(freeing.begin(), freeing.end());
        freeing.erase(std::unique(freeing.begin(), freeing.end()), freeing.end());
        if (!freeing.empty()) {
            for (std::size_t n : freeing) {
                out.mFeature[n] = 0;
                mode[n] = HexNodeMode::Slide;
                tether_chain[n] = out.mFeatures.mNodeChain[n];
                tether_corner[n] = out.mFeatures.mNodeCorner[n];
            }
            out.mReleased = static_cast<std::int64_t>(freeing.size());
            HexUntangleOptions again = uo;
            again.mIterations = ro.mUntangleIterations / 5;
            const HexUntangleReport more =
                hex_untangle(pts, rIn.rMesh.mHexes, mode, project, again);
            out.mMoves += more.mMoves;
            out.mWorst = 1.0;
            for (const auto& hex : rIn.rMesh.mHexes)
                out.mWorst = std::min(out.mWorst, hex_min_scaled_jacobian(pts, hex));
            log::debug("{}repair: {} feature nodes released, smallest corner {}", kOctreeMeshPrefix,
                       freeing.size(), out.mWorst);
        }
    }
    log::debug("{}{} fit: smallest corner scaled Jacobian {} after {} sweeps", kOctreeMeshPrefix,
               robust ? "robust" : "body-fitted", out.mWorst, ur.mSweeps);
    return out;
}

}  // namespace

void octree_mesh_hexahedra(OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR) {
    octree_mesh_report_tree(rS.mOctree.mTree, rR);
    HexMesh hm;
    std::vector<std::int64_t> cell_level;
    std::vector<double> size;
    if (rO.mHexGrading == OctreeMeshHexGrading::Graded)
        octmh_graded(rS, rO, rR, hm, cell_level, size);
    else
        octmh_uniform(rS, rO, rR, hm, cell_level, size);
    const std::size_t num_grid = hm.mHexes.size();

    // The boundary, pillowed or not, onto the skin.
    const HexBoundary grid_boundary = hex_boundary(hm);
    if (rO.mHexGrading == OctreeMeshHexGrading::Graded) {
        // A boundary node's size is its edges along the boundary: an edge into a
        // coarser cell inside would push its pillow past the nodes behind it.
        std::vector<double> along(hm.mPoints.size(), 0.0);
        std::vector<std::int64_t> count(hm.mPoints.size(), 0);
        for (const auto& q : grid_boundary.mQuads)
            for (std::size_t e = 0; e < 4; ++e) {
                const auto a = static_cast<std::size_t>(q[e]);
                const auto b = static_cast<std::size_t>(q[(e + 1) % 4]);
                const double len = vec3_norm(vec3_sub(hm.mPoints[a], hm.mPoints[b]));
                along[a] += len;
                along[b] += len;
                ++count[a];
                ++count[b];
            }
        for (std::size_t v = 0; v < along.size(); ++v)
            if (count[v] > 0)
                size[v] = along[v] / static_cast<double>(count[v]);
    }
    std::vector<std::int64_t> boundary_nodes;
    for (const auto& q : grid_boundary.mQuads)
        boundary_nodes.insert(boundary_nodes.end(), q.begin(), q.end());
    std::sort(boundary_nodes.begin(), boundary_nodes.end());
    boundary_nodes.erase(std::unique(boundary_nodes.begin(), boundary_nodes.end()),
                         boundary_nodes.end());
    octmh_log_stage("grid", hm.mPoints, hm.mHexes);
    // The boundary's normals, from the grid's own positions.
    const std::vector<Vec3> normals =
        hex_boundary_normals(hm.mPoints, grid_boundary.mQuads, kOctmhNormalSmoothing);
    std::vector<std::int64_t> copy;
    if (rO.mPillow) {
        copy = hex_pillow(hm, grid_boundary);
        size.resize(hm.mPoints.size());
        for (std::size_t v = 0; v < copy.size(); ++v)
            if (copy[v] >= 0)
                size[static_cast<std::size_t>(copy[v])] = size[v];
    }
    const HexBoundary boundary = rO.mPillow ? hex_boundary(hm) : grid_boundary;
    std::optional<FeatureGraph> graph;
    if (rO.mPreserveFeatures && !rS.mFeatures.empty()) {
        std::vector<std::array<std::int64_t, 2>> edges;
        for (const OctreeFeatureEdge& e : rS.mFeatures)
            edges.push_back({e.mA, e.mB});
        graph = feature_graph(rS.mSkin.mPoints, edges, rO.mCornerAngle);
    }

    // Fit the boundary to the skin and untangle; keep the body-fitted fit of
    // 18.0.0 where the robust one does no better.
    const OctmhFitInput in{rS,
                           rO,
                           hm,
                           boundary_nodes,
                           copy,
                           size,
                           normals,
                           boundary.mQuads,
                           graph ? &*graph : nullptr};
    const char* p_env = std::getenv("MESHIOPLUSPLUS_OCTREE_MESH_HEX_FIT");
    const std::string forced = p_env != nullptr ? std::string(p_env) : std::string();
    OctmhFit fit;
    if (forced == "body_fitted") {
        fit = octmh_fit(in, OctmhPlacement::BodyFitted);
    } else {
        fit = octmh_fit(in, OctmhPlacement::Robust);
        if (fit.mWorst < kOctmhFallbackBelow && forced != "robust") {
            OctmhFit body = octmh_fit(in, OctmhPlacement::BodyFitted);
            log::debug("{}robust fit: smallest corner {}, body-fitted: {}", kOctreeMeshPrefix,
                       fit.mWorst, body.mWorst);
            if (body.mWorst > fit.mWorst) {
                const double robust_worst = fit.mWorst;
                fit = std::move(body);
                fit.mRobustWorst = robust_worst;
                fit.mFellBack = true;
            }
        }
    }
    hm.mPoints = std::move(fit.mPoints);
    const std::vector<std::uint8_t>& feature = fit.mFeature;
    if (graph) {
        const OctmhFeatures& f = fit.mFeatures;
        rR.mNumCorners = f.mCorners;
        rR.mNumCornersSnapped = f.mCornersSnapped;
        rR.mNumFeatureChains = f.mChains;
        rR.mNumFeatureChainsMissed = f.mChainsMissed;
        rR.mNumFeatureMoves = f.mMoved;
        if (f.mCornersSnapped < f.mCorners)
            octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) +
                                     std::to_string(f.mCorners - f.mCornersSnapped) + " of " +
                                     std::to_string(f.mCorners) + " feature corners not recovered");
        if (f.mChainsMissed > 0)
            octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(f.mChainsMissed) +
                                     " of " + std::to_string(f.mChains) +
                                     " feature chains not recovered");
    }
    rR.mNumOptimizedMoves = fit.mMoves;
    const double worst = fit.mWorst;
    if (fit.mReleased > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(fit.mReleased) +
                                 " feature nodes were released from their feature to untangle "
                                 "the hexahedra around them");
    if (fit.mFellBack)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) +
                                 "the robust fit reached a smallest corner scaled Jacobian of " +
                                 std::to_string(fit.mRobustWorst) +
                                 "; the body-fitted fit did better and was kept");
    if (worst <= 0.0)
        octree_mesh_warn(
            rR, std::string(kOctreeMeshPrefix) +
                    "untangling left a hexahedron with a non-positive corner "
                    "Jacobian; raise untangle_iterations or lower " +
                    (rO.mHexGrading == OctreeMeshHexGrading::Graded ? "min_size" : "surface_size"));
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
        lv.As<std::int64_t>()[c] =
            c < num_grid
                ? cell_level[c]
                : cell_level[static_cast<std::size_t>(grid_boundary.mSides[c - num_grid][0])];
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
