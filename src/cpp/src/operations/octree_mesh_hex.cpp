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
                d = feature_snap_on_polyline(rSkin, chain, rPoints[static_cast<std::size_t>(V)])
                        .second;
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
            rPoints[static_cast<std::size_t>(v)] =
                feature_snap_on_polyline(rSkin, chain, rPoints[static_cast<std::size_t>(v)]).first;
            rFeature[static_cast<std::size_t>(v)] = 2;
            rep.mNodeChain[static_cast<std::size_t>(v)] = static_cast<std::int64_t>(ci);
            ++rep.mMoved;
        }
    }
    return rep;
}

// Per hexahedron of a grid, its winding number and its centre (for materials).
struct OctmhCells {
    std::vector<double> mWinding;
    std::vector<Vec3> mCentre;
};

// The uniform grid: every cell at the skin's level, the kept ones made
// manifold. @p rLevel gets each hexahedron's level, @p rSize each node's size.
void octmh_uniform(OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR,
                   HexMesh& hm, std::vector<std::int64_t>& rLevel, std::vector<double>& rSize,
                   OctmhCells& rCells) {
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
        rCells.mWinding.push_back(w[static_cast<std::size_t>(c)]);
        rCells.mCentre.push_back(centres[static_cast<std::size_t>(c)]);
    }
    rLevel.assign(hm.mHexes.size(), level);
    rSize.assign(hm.mPoints.size(), h);
}

// The graded hexahedra: the all-hexahedral dual of the octree, paired and
// with its transitions replaced by Livesu, Pitzalis and Cherchi's schemes
// (detail/hex_dual_schemes.hpp). A dual hexahedron is kept when the primal
// vertex it surrounds is inside; @p rSize is each node's mean edge.
void octmh_graded(OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR,
                  HexMesh& hm, std::vector<std::int64_t>& rLevel, std::vector<double>& rSize,
                  OctmhCells& rCells) {
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
    std::vector<double> w = octree_mesh_winding(rS, centres, rO, rR);
    // A centre is a primal vertex, on the lattice's planes, and a skin may lie on
    // them too (a cap on a lattice plane): a hexahedron kept by a hair there
    // would be pressed flat by the projection and the pillow. Keep only centres
    // more than a tenth of their cell inside the boundary. A nested shell
    // inside (cleaving's interfaces) is not the boundary, and dropping the
    // centres beside it would open a cavity there: a centre within a tenth of
    // its cell of the skin is dropped only when some point a tenth of its cell
    // away along an axis, or just past its nearest skin point, is outside (a
    // plane through the centre always leaves one axis direction outside).
    // Per close centre: six axis probes, one past its nearest skin point, and
    // eight along the cube's diagonals (off the lattice's planes, for a centre
    // on a skin lying on them, where the axis probes are as ambiguous).
    constexpr std::size_t kProbes = 15;
    std::vector<std::size_t> close;  // the centres within a tenth of a cell
    std::vector<double> close_d;
    std::vector<double> pw;  // their probes' winding numbers, kProbes each
    {
        const std::vector<SurfaceProjection> proj = query_surface_projections(rS.mQuery, centres);
        std::vector<Vec3> probes;
        for (std::size_t c = 0; c < centres.size(); ++c) {
            const double cell =
                frame.mRootSize / static_cast<double>(std::int64_t{1} << dual.mLevel[c]);
            const Vec3 to = vec3_sub(proj[c].mPoint, centres[c]);
            const double d = vec3_norm(to);
            if (!(d < 0.1 * cell))
                continue;
            close.push_back(c);
            close_d.push_back(d / cell);
            for (std::size_t k = 0; k < 6; ++k) {
                Vec3 p = centres[c];
                p[k / 2] += (k % 2 ? -0.1 : 0.1) * cell;
                probes.push_back(p);
            }
            probes.push_back(
                d > 0.0 ? vec3_add(proj[c].mPoint, vec3_scale(to, std::max(d, 0.05 * cell) / d))
                        : centres[c]);
            const double step = 0.1 * cell / std::sqrt(3.0);
            for (int k = 0; k < 8; ++k)
                probes.push_back(vec3_add(
                    centres[c],
                    Vec3{(k & 1 ? -step : step), (k & 2 ? -step : step), (k & 4 ? -step : step)}));
        }
        pw = fast_windings(rS.mBvh, probes, rO.mWindingBeta);
        for (double& v : pw)
            v *= rS.mWindingSign;
    }
    // On the skin itself the winding number is half-way between its two sides
    // only up to the sign of a zero (a point in a triangle's plane): take the
    // middle of the diagonal probes' range where it disagrees by more than a
    // quarter (a nested shell read as the boundary).
    std::vector<std::uint8_t> on_skin(close.size(), 0);
    for (std::size_t i = 0; i < close.size(); ++i) {
        if (close_d[i] > 1e-9)
            continue;
        on_skin[i] = 1;
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (std::size_t k = 7; k < kProbes; ++k) {
            lo = std::min(lo, pw[i * kProbes + k]);
            hi = std::max(hi, pw[i * kProbes + k]);
        }
        const double mid = 0.5 * (lo + hi);
        if (std::abs(mid - w[close[i]]) > 0.25)
            w[close[i]] = mid;
    }
    std::vector<std::uint8_t> kept(centres.size(), 0);
    for (std::size_t c = 0; c < kept.size(); ++c)
        kept[c] = w[c] > rO.mWindingThreshold;
    for (std::size_t i = 0; i < close.size(); ++i) {
        if (!kept[close[i]])
            continue;
        // On the skin, only the diagonal probes; off it, the axis probes and
        // the one past the nearest skin point (as before the diagonals).
        for (std::size_t k = on_skin[i] ? 7 : 0; k < (on_skin[i] ? kProbes : 7); ++k) {
            if (pw[i * kProbes + k] <= rO.mWindingThreshold) {
                kept[close[i]] = 0;
                break;
            }
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
        rCells.mWinding.push_back(w[c]);
        rCells.mCentre.push_back(centres[c]);
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

// The shells of a skin with nested parts, under cleaving: per soup triangle
// its shell (the material on its inner side: 1 for the outer skin, L for the
// interface between materials L-1 and L), per skin point the shell of its
// triangles, and a closest-point query per shell.
struct OctmhShells {
    std::vector<std::int32_t> mTriShell;
    std::vector<std::int32_t> mPointShell;
    /// Shell L's triangles at L - 1, and their query, which points into its
    /// soup: the soups are filled before any query is built and never move.
    std::vector<TriangleSoup> mSoup;
    std::vector<DistanceQuery> mQuery;
};

// The interface between materials L-1 and L, pillowed on both sides.
struct OctmhInterface {
    std::int32_t mShell = 0;
    std::vector<std::int64_t> mNodes;                  // the outer-side nodes
    std::vector<std::int64_t> mShellNode, mInnerNode;  // per node
    std::vector<Vec3> mNormals;                        // per node, unit, towards L-1
    std::vector<std::array<std::int64_t, 4>> mQuads;   // on the shell nodes, towards L-1
    std::optional<FeatureGraph> mGraph;                // the shell's features
};

// The closest-point query of shell @p Shell (the whole skin without shells).
const DistanceQuery& octmh_query(const OctreeMeshSetup& rS, const OctmhShells* pShells,
                                 std::int32_t Shell) {
    return pShells != nullptr && Shell > 0 ? pShells->mQuery[static_cast<std::size_t>(Shell - 1)]
                                           : rS.mQuery;
}

// Each hexahedron corner's three edge neighbours (VTK order).
constexpr int kOctmhCorner[8][3] = {{1, 3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2, 7},
                                    {7, 5, 0}, {4, 6, 1}, {5, 7, 2}, {6, 4, 3}};

// What the fit works on: the pillowed (or bare) mesh at the grid's positions.
struct OctmhFitInput {
    const OctreeMeshSetup& rS;
    const OctreeMeshOptions& rO;
    const HexMesh& rMesh;
    const std::vector<std::int64_t>& rBoundaryNodes;  // the grid's boundary nodes
    const std::vector<std::int64_t>& rCopy;           // their pillow copies (pillow only)
    const std::vector<double>& rSize;
    const std::vector<Vec3>& rNormals;                         // the grid boundary's unit normals
    const std::vector<std::array<std::int64_t, 4>>& rQuads;    // the fitted mesh's boundary
    const FeatureGraph* pGraph;                                // null: no feature snap
    const OctmhShells* pShells = nullptr;                      // cleaving: the shells
    const std::vector<OctmhInterface>* pInterfaces = nullptr;  // cleaving: the interfaces
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
    std::int64_t mReleased = 0;      // feature nodes released to untangle
    std::int64_t mShortChains = 0;   // chains shorter than a cell, not snapped
    std::int64_t mShortCorners = 0;  // corners only such chains end at
};

// Where the line through @p rX along the unit @p rDir first meets the skin
// from outside (a triangle facing along @p rDir), nearest @p rX within
// @p Reach either way; false when it misses.
bool octmh_ray_foot(const OctreeMeshSetup& rS, const Vec3& rX, const Vec3& rDir, double Reach,
                    Vec3& rFoot, const OctmhShells* pShells = nullptr, std::int32_t Shell = 0) {
    const double sign = rS.mWindingSign;
    auto facing = [&](std::int64_t T, const Vec3& rN) {
        return sign * vec3_dot(rN, rDir) > 0.0 &&
               (pShells == nullptr || pShells->mTriShell[static_cast<std::size_t>(T)] == Shell);
    };
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
// outward) to the far wall of the skin, or @p Reach when it is farther. With
// shells, any wall counts (a nested shell faces the other way).
double octmh_thickness(const OctreeMeshSetup& rS, const Vec3& rQ, const Vec3& rOut, double Reach,
                       const OctmhShells* pShells = nullptr) {
    const Vec3 in = vec3_scale(rOut, -1.0);
    const double sign = rS.mWindingSign;
    auto far_wall = [&](std::int64_t, const Vec3& rN) {
        return pShells != nullptr || sign * vec3_dot(rN, in) > 0.0;
    };
    const TriRayHit hit = tri_bvh_ray_filtered(rS.mBvh, rQ, in, 1e-9 * Reach, Reach, far_wall);
    return hit.mTriangle < 0 ? Reach : hit.mT;
}

// Put the grid node @p V of the pillow copy @p Copy half a size @p Size inside
// the skin along -@p rOut, but no more than 0.45 of the way to the far wall,
// so a part thinner than a cell does not fold. @return whether it was held.
bool octmh_inset(const OctreeMeshSetup& rS, std::vector<Vec3>& rPoints, std::size_t V,
                 std::size_t Copy, const Vec3& rOut, double Size,
                 const OctmhShells* pShells = nullptr) {
    const Vec3 q = rPoints[Copy];
    const double t = octmh_thickness(rS, q, rOut, 1.2 * Size, pShells);
    const double d = std::min(0.5 * Size, 0.45 * t);
    rPoints[V] = vec3_sub(q, vec3_scale(rOut, d));
    return d < 0.5 * Size;
}

// The robust foot of a point @p rX whose closest skin point is @p Distance
// away: where the line through it along the unit normal @p rNormal meets the
// skin, when that is no more than half a size @p Size farther than the
// closest point; false otherwise (the caller keeps the closest point).
bool octmh_robust_foot(const OctreeMeshSetup& rS, const Vec3& rX, const Vec3& rNormal,
                       double Distance, double Size, Vec3& rFoot,
                       const OctmhShells* pShells = nullptr, std::int32_t Shell = 0) {
    if (!(vec3_norm(rNormal) > 0.0))
        return false;
    Vec3 f{};
    if (!octmh_ray_foot(rS, rX, rNormal, Distance + 1.5 * Size, f, pShells, Shell))
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
    const auto [on, d] = feature_snap_on_polyline(rSkin, rChain, rP);
    if (d <= Reach)
        return rP;
    return vec3_add(on, vec3_scale(vec3_sub(rP, on), Reach / d));
}

// The repair after untangling: when a hexahedron's smallest corner is below
// this, the feature nodes of every such hexahedron are released (tethered
// within `kOctmhTether` sizes of their chain or corner) and the mesh is
// untangled again for a fifth of the sweeps.
constexpr double kOctmhRepairBelow = 0.05;

// The feature graph without the chains shorter than a cell: on a scanned skin
// a short crease is noise the boundary cannot follow, and pinning it holds
// the hexahedra around it flat or inverted. A corner is kept when a kept chain
// ends at it. @p Uniform is the uniform grid's edge (0: graded, where the
// cell is the leaf's at the chain's midpoint). @p rDropped counts the chains
// dropped.
FeatureGraph octmh_drop_short(const OctreeMeshSetup& rS, const FeatureGraph& rGraph, double Uniform,
                              std::int64_t& rDropped) {
    const std::vector<Vec3>& skin = rS.mSkin.mPoints;
    FeatureGraph out;
    rDropped = 0;
    std::vector<std::int64_t> ends;
    for (const auto& chain : rGraph.mChains) {
        double length = 0.0;
        for (std::size_t i = 0; i + 1 < chain.size(); ++i)
            length += vec3_norm(vec3_sub(skin[static_cast<std::size_t>(chain[i + 1])],
                                         skin[static_cast<std::size_t>(chain[i])]));
        const double cell = Uniform > 0.0 ? Uniform
                                          : octree_mesh_local_size(
                                                rS, octmh_along_chain(skin, chain, 0.5 * length));
        if (length < cell) {
            ++rDropped;
            continue;
        }
        out.mChains.push_back(chain);
        ends.push_back(chain.front());
        ends.push_back(chain.back());
    }
    std::sort(ends.begin(), ends.end());
    for (std::int64_t c : rGraph.mCorners)
        if (std::binary_search(ends.begin(), ends.end(), c))
            out.mCorners.push_back(c);
    return out;
}

// Place the boundary on the skin, snap the features and untangle.
OctmhFit octmh_fit(const OctmhFitInput& rIn, OctmhPlacement Placement) {
    const OctreeMeshSetup& rs = rIn.rS;
    const OctreeMeshOptions& ro = rIn.rO;
    const OctmhShells* p_shells = rIn.pShells;
    const bool robust = Placement == OctmhPlacement::Robust;
    OctmhFit out;
    out.mPoints = rIn.rMesh.mPoints;
    std::vector<Vec3>& pts = out.mPoints;
    const std::vector<double>& size = rIn.rSize;
    std::vector<Vec3> at_nodes;
    at_nodes.reserve(rIn.rBoundaryNodes.size());
    for (std::int64_t v : rIn.rBoundaryNodes)
        at_nodes.push_back(pts[static_cast<std::size_t>(v)]);
    const std::vector<SurfaceProjection> proj =
        query_surface_projections(octmh_query(rs, p_shells, 1), at_nodes);
    // Per point, the direction a sliding node is projected along (robust).
    std::vector<Vec3> slide(robust ? pts.size() : 0, Vec3{0.0, 0.0, 0.0});
    std::vector<HexNodeMode> mode(pts.size(), HexNodeMode::Free);
    // Per node on the skin (robust): its outward direction, and the nodes that
    // stay a pillow's depth inside and outside it.
    std::vector<Vec3> outward(robust ? pts.size() : 0, Vec3{0.0, 0.0, 0.0});
    std::vector<std::int64_t> inside_of(robust ? pts.size() : 0, -1);
    std::vector<std::int64_t> outside_of(robust ? pts.size() : 0, -1);
    // Per point, the shell a sliding node slides on (shells only).
    std::vector<std::int32_t> shell_of(p_shells != nullptr ? pts.size() : 0, 0);
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
        if (octmh_robust_foot(rs, pts[v], nb, proj[b].mDistance, size[v], f, p_shells, 1)) {
            q = f;
            out_dir = nb;
            ++along_normal;
        }
        const std::size_t moved = ro.mPillow ? static_cast<std::size_t>(rIn.rCopy[v]) : v;
        slide[moved] = nb;
        pts[moved] = q;
        mode[moved] = HexNodeMode::Slide;
        if (!shell_of.empty())
            shell_of[moved] = 1;
        if (ro.mPillow) {
            outward[moved] = out_dir;
            inside_of[moved] = static_cast<std::int64_t>(v);
            held += octmh_inset(rs, pts, v, moved, out_dir, size[v], p_shells);
        }
    }
    // The interfaces (cleaving, robust): each shell node on its shell, the
    // outer node a pillow's depth outside it, the inner copy as deep inside.
    // An interface can lie half a cell from its shell (the graded dual keeps
    // the centres more than a tenth of a cell inside), so the node a pillow's
    // depth out could land on its neighbour: it goes at most half way to the
    // nearest neighbour on its side that is not on the interface.
    std::vector<std::vector<std::int64_t>> side_nbr;
    std::vector<std::uint8_t> on_interface;
    if (robust && rIn.pInterfaces != nullptr) {
        on_interface.assign(pts.size(), 0);
        for (const OctmhInterface& in : *rIn.pInterfaces)
            for (std::size_t k = 0; k < in.mNodes.size(); ++k) {
                on_interface[static_cast<std::size_t>(in.mNodes[k])] = 1;
                on_interface[static_cast<std::size_t>(in.mShellNode[k])] = 1;
                on_interface[static_cast<std::size_t>(in.mInnerNode[k])] = 1;
            }
        side_nbr.resize(pts.size());
        for (const auto& hex : rIn.rMesh.mHexes)
            for (std::size_t c = 0; c < 8; ++c) {
                const auto u = static_cast<std::size_t>(hex[c]);
                if (!on_interface[u])
                    continue;
                for (std::size_t e = 0; e < 3; ++e) {
                    const std::int64_t n = hex[static_cast<std::size_t>(kOctmhCorner[c][e])];
                    if (!on_interface[static_cast<std::size_t>(n)])
                        side_nbr[u].push_back(n);
                }
            }
    }
    // Move @p V from the shell point @p rFoot along the unit @p rDir by
    // `octmh_inset`'s depth, held to half the gap to its side's neighbours.
    auto place_side = [&](std::size_t V, std::size_t Shell, const Vec3& rDir, double Size) {
        const Vec3 foot = pts[Shell];
        double gap = std::numeric_limits<double>::infinity();
        for (std::int64_t n : side_nbr[V]) {
            const double along = vec3_dot(vec3_sub(pts[static_cast<std::size_t>(n)], foot), rDir);
            if (along > 1e-9 * Size)
                gap = std::min(gap, along);
        }
        const double t = octmh_thickness(rs, foot, vec3_scale(rDir, -1.0), 1.2 * Size, p_shells);
        const double d = std::min({0.5 * Size, 0.45 * t, 0.5 * gap});
        pts[V] = vec3_add(foot, vec3_scale(rDir, d));
        return d < 0.5 * Size;
    };
    if (robust && rIn.pInterfaces != nullptr)
        for (const OctmhInterface& in : *rIn.pInterfaces) {
            std::vector<Vec3> at;
            for (std::int64_t v : in.mNodes)
                at.push_back(pts[static_cast<std::size_t>(v)]);
            const std::vector<SurfaceProjection> near =
                query_surface_projections(octmh_query(rs, p_shells, in.mShell), at);
            for (std::size_t k = 0; k < in.mNodes.size(); ++k) {
                const auto v = static_cast<std::size_t>(in.mNodes[k]);
                const auto sn = static_cast<std::size_t>(in.mShellNode[k]);
                const auto inner = static_cast<std::size_t>(in.mInnerNode[k]);
                Vec3 q = near[k].mPoint;
                Vec3 out_dir = vec3_scale(vec3_normalize(near[k].mNormal), rs.mWindingSign);
                const Vec3& nb = in.mNormals[k];
                Vec3 f{};
                if (octmh_robust_foot(rs, pts[v], nb, near[k].mDistance, size[v], f, p_shells,
                                      in.mShell)) {
                    q = f;
                    out_dir = nb;
                    ++along_normal;
                }
                pts[sn] = q;
                slide[sn] = nb;
                mode[sn] = HexNodeMode::Slide;
                shell_of[sn] = in.mShell;
                outward[sn] = out_dir;
                inside_of[sn] = static_cast<std::int64_t>(inner);
                outside_of[sn] = static_cast<std::int64_t>(v);
                held += place_side(inner, sn, vec3_scale(out_dir, -1.0), size[v]);
                held += place_side(v, sn, out_dir, size[v]);
            }
        }
    if (robust)
        log::debug("{}{} skin nodes placed along their normal, {} held by the thickness",
                   kOctreeMeshPrefix, along_normal, held);
    octmh_log_stage(robust ? "robust placement" : "body-fitted placement", pts, rIn.rMesh.mHexes);
    out.mFeature.assign(pts.size(), 0);
    // The feature snaps: the outer skin's, then each interface's on its own
    // shell. The robust fit leaves out creases shorter than a cell.
    struct Job {
        const FeatureGraph* pGraph;
        const std::vector<std::array<std::int64_t, 4>>* pQuads;
    };
    std::vector<Job> jobs;
    if (rIn.pGraph != nullptr)
        jobs.push_back({rIn.pGraph, &rIn.rQuads});
    if (robust && rIn.pInterfaces != nullptr)
        for (const OctmhInterface& in : *rIn.pInterfaces)
            if (in.mGraph)
                jobs.push_back({&*in.mGraph, &in.mQuads});
    // Every snapped chain, across the jobs, and per point the one a node lies
    // on and what a released node is tethered to (-1: none).
    std::vector<FeatureGraph> kept(jobs.size());
    std::vector<const std::vector<std::int64_t>*> chains;
    std::vector<std::int64_t> node_chain, node_corner, tether_chain, tether_corner;
    if (robust && !jobs.empty()) {
        node_chain.assign(pts.size(), -1);
        node_corner.assign(pts.size(), -1);
        tether_chain.assign(pts.size(), -1);
        tether_corner.assign(pts.size(), -1);
    }
    for (std::size_t j = 0; j < jobs.size(); ++j) {
        const FeatureGraph* p_graph = jobs[j].pGraph;
        std::int64_t short_chains = 0;
        if (robust) {
            kept[j] = octmh_drop_short(
                rs, *p_graph, ro.mHexGrading == OctreeMeshHexGrading::Graded ? 0.0 : size.front(),
                short_chains);
            p_graph = &kept[j];
        }
        OctmhFeatures f = octmh_snap(pts, *jobs[j].pQuads, *p_graph, rs.mSkin.mPoints, size,
                                     out.mFeature, robust);
        // Report against the whole graph: a dropped crease is not recovered.
        out.mFeatures.mCorners += static_cast<std::int64_t>(jobs[j].pGraph->mCorners.size());
        out.mFeatures.mChains += static_cast<std::int64_t>(jobs[j].pGraph->mChains.size());
        out.mFeatures.mCornersSnapped += f.mCornersSnapped;
        out.mFeatures.mChainsMissed += f.mChainsMissed + short_chains;
        out.mFeatures.mMoved += f.mMoved;
        out.mShortChains += short_chains;
        out.mShortCorners +=
            static_cast<std::int64_t>(jobs[j].pGraph->mCorners.size() - p_graph->mCorners.size());
        if (!robust)
            continue;
        const auto offset = static_cast<std::int64_t>(chains.size());
        for (const auto& chain : p_graph->mChains)
            chains.push_back(&chain);
        for (std::size_t p = 0; p < pts.size(); ++p) {
            if (f.mNodeChain[p] >= 0)
                node_chain[p] = f.mNodeChain[p] + offset;
            if (f.mNodeCorner[p] >= 0)
                node_corner[p] = f.mNodeCorner[p];
        }
        const std::vector<std::int64_t> turns =
            octmh_release_turns(*jobs[j].pQuads, f, out.mFeature);
        for (std::size_t p = 0; p < pts.size(); ++p)
            if (turns[p] >= 0)
                tether_chain[p] = turns[p] + offset;
    }
    if (!jobs.empty()) {
        if (robust && ro.mPillow)
            // A snapped node can move far (up to 1.5 sizes for a corner): the
            // nodes a pillow's depth from it follow, so its cells keep their height.
            for (std::size_t p = 0; p < pts.size(); ++p) {
                if (!out.mFeature[p])
                    continue;
                if (outside_of[p] >= 0) {
                    // An interface's shell node: both sides follow it.
                    const auto v = static_cast<std::size_t>(outside_of[p]);
                    place_side(static_cast<std::size_t>(inside_of[p]), p,
                               vec3_scale(outward[p], -1.0), size[v]);
                    place_side(v, p, outward[p], size[v]);
                } else if (inside_of[p] >= 0) {
                    const auto v = static_cast<std::size_t>(inside_of[p]);
                    octmh_inset(rs, pts, v, p, outward[p], size[v], p_shells);
                }
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
            const std::int32_t sh = shell_of.empty() ? 0 : shell_of[Node];
            const DistanceQuery& query = octmh_query(rs, p_shells, sh);
            const SurfaceProjection c = query_surface_projections(query, {rP}).front();
            Vec3 f = c.mPoint;
            octmh_robust_foot(rs, rP, slide[Node], c.mDistance, size[Node], f, p_shells, sh);
            if (tether_chain.empty())
                return f;
            Vec3 held_at = f;
            const double within = kOctmhTether * size[Node];
            if (tether_chain[Node] >= 0) {
                const auto& chain = *chains[static_cast<std::size_t>(tether_chain[Node])];
                held_at = octmh_tether(rs.mSkin.mPoints, chain, f, within);
            } else if (tether_corner[Node] >= 0) {
                held_at = octmh_tether(rs.mSkin.mPoints, {tether_corner[Node], tether_corner[Node]},
                                       f, within);
            }
            return held_at != f ? query_closest_points(query, {held_at}).front().mPoint : f;
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
    if (robust && !jobs.empty() && ro.mUntangleIterations >= 5 && out.mWorst < kOctmhRepairBelow) {
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
                tether_chain[n] = node_chain[n];
                tether_corner[n] = node_corner[n];
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

// The material of a winding number @p W: 1 inside the skin, L inside the
// (L-1)-th nested shell, at most @p Levels.
std::int32_t octmh_material(double W, double Threshold, std::int32_t Levels) {
    return std::clamp<std::int32_t>(static_cast<std::int32_t>(std::floor(W - Threshold)) + 1, 1,
                                    Levels);
}

// The shells of the skin under cleaving: each soup triangle takes the material
// just inside it (the winding number a thousandth of its shortest edge inside
// its centroid), each skin point the shell of its triangles, and each shell a
// closest-point query of its own triangles.
OctmhShells octmh_shells(const OctreeMeshSetup& rS, const OctreeMeshOptions& rO,
                         std::int32_t Levels) {
    const TriangleSoup& soup = rS.mSkin.mSoup;
    const std::size_t nt = soup.NumTriangles();
    std::vector<Vec3> probe(nt);
    for (std::size_t t = 0; t < nt; ++t) {
        const Vec3& a = soup.mCorners[3 * t];
        const Vec3& b = soup.mCorners[3 * t + 1];
        const Vec3& c = soup.mCorners[3 * t + 2];
        const Vec3 centroid = vec3_scale(vec3_add(vec3_add(a, b), c), 1.0 / 3.0);
        const Vec3 n = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
        const double edge = std::min(
            {vec3_norm(vec3_sub(b, a)), vec3_norm(vec3_sub(c, b)), vec3_norm(vec3_sub(a, c))});
        probe[t] =
            vec3_norm(n) > 0.0
                ? vec3_sub(centroid,
                           vec3_scale(vec3_normalize(n),
                                      rS.mWindingSign * std::max(1e-3 * edge, 1e-9 * rS.mDiagonal)))
                : centroid;
    }
    const std::vector<double> w = fast_windings(rS.mBvh, probe, rO.mWindingBeta);
    OctmhShells out;
    out.mTriShell.resize(nt);
    for (std::size_t t = 0; t < nt; ++t)
        out.mTriShell[t] = octmh_material(rS.mWindingSign * w[t], rO.mWindingThreshold, Levels);
    out.mPointShell.assign(rS.mSkin.mPoints.size(), 0);
    for (std::size_t t = 0; t < rS.mSkin.mTris.size(); ++t)
        for (std::int64_t v : rS.mSkin.mTris[t])
            out.mPointShell[static_cast<std::size_t>(v)] = out.mTriShell[t];
    SurfaceDistanceOptions dopt;
    dopt.mWatertightCheck = SdfWatertightCheck::Off;
    out.mSoup.resize(static_cast<std::size_t>(Levels));
    for (std::int32_t level = 1; level <= Levels; ++level) {
        TriangleSoup& sub = out.mSoup[static_cast<std::size_t>(level - 1)];
        sub.mPoints = soup.mPoints;
        for (std::size_t t = 0; t < nt; ++t)
            if (out.mTriShell[t] == level) {
                for (std::size_t k = 0; k < 3; ++k)
                    sub.mCorners.push_back(soup.mCorners[3 * t + k]);
                sub.mSourceCell.push_back(soup.mSourceCell[t]);
                sub.mVertices.push_back(soup.mVertices[t]);
            }
    }
    for (std::int32_t level = 1; level <= Levels; ++level) {
        const TriangleSoup& sub = out.mSoup[static_cast<std::size_t>(level - 1)];
        // A shell with no triangle of its own (a level the cells reach only
        // through the winding number's error) projects onto the whole skin.
        out.mQuery.push_back(build_distance_query(sub.NumTriangles() > 0 ? sub : soup, dopt));
    }
    return out;
}

// Nest the materials of the hexahedra of @p rMesh (one per hexahedron in
// @p rMaterial, 1 to @p Levels): for each level L from 2, the set S_L of the
// hexahedra of material L or more keeps off the nodes of everything outside
// S_{L-1} and the mesh's boundary (so each interface is interior and one
// layer from the next), and is made manifold by adding hexahedra of S_{L-1}
// (the largest @p rScore first), else by removing its own (the smallest
// first). @return the levels left with no hexahedron.
std::int32_t octmh_nest(const HexMesh& rMesh, std::vector<std::int32_t>& rMaterial,
                        const std::vector<double>& rScore, std::int32_t Levels) {
    const std::size_t nh = rMesh.mHexes.size(), np = rMesh.mPoints.size();
    std::vector<std::uint8_t> on_boundary(np, 0);
    for (const auto& q : hex_boundary(rMesh).mQuads)
        for (std::int64_t v : q)
            on_boundary[static_cast<std::size_t>(v)] = 1;
    std::int32_t empty = 0;
    for (std::int32_t level = 2; level <= Levels; ++level) {
        std::vector<std::uint8_t> fenced = on_boundary;
        for (std::size_t h = 0; h < nh; ++h)
            if (rMaterial[h] < level - 1)
                for (std::int64_t v : rMesh.mHexes[h])
                    fenced[static_cast<std::size_t>(v)] = 1;
        auto touches = [&](std::size_t H) {
            for (std::int64_t v : rMesh.mHexes[H])
                if (fenced[static_cast<std::size_t>(v)])
                    return true;
            return false;
        };
        std::vector<std::uint8_t> kept(nh, 0), allowed(nh, 0);
        for (std::size_t h = 0; h < nh; ++h) {
            if (rMaterial[h] < level - 1 || touches(h))
                continue;
            allowed[h] = 1;
            kept[h] = rMaterial[h] >= level;
        }
        hex_make_manifold(rMesh, kept, rScore, &allowed);
        hex_shrink_manifold(rMesh, kept, rScore);
        bool any = false;
        for (std::size_t h = 0; h < nh; ++h) {
            if (kept[h])
                rMaterial[h] = std::max(rMaterial[h], level);
            else if (rMaterial[h] >= level)
                rMaterial[h] = level - 1;
            any = any || kept[h];
        }
        empty += any ? 0 : 1;
    }
    return empty;
}

}  // namespace

void octree_mesh_hexahedra(OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR) {
    octree_mesh_report_tree(rS.mOctree.mTree, rR);
    HexMesh hm;
    std::vector<std::int64_t> cell_level;
    std::vector<double> size;
    OctmhCells cells;
    if (rO.mHexGrading == OctreeMeshHexGrading::Graded)
        octmh_graded(rS, rO, rR, hm, cell_level, size, cells);
    else
        octmh_uniform(rS, rO, rR, hm, cell_level, size, cells);
    const std::size_t num_grid = hm.mHexes.size();

    // Materials (cleaving): one per winding level the cells reach, nested and
    // made manifold, each interface pillowed below.
    std::int32_t levels = 1;
    if (rO.mConformity == OctreeMeshConformity::Cleaving) {
        double wmax = rO.mWindingThreshold;
        for (double w : cells.mWinding)
            wmax = std::max(wmax, w);
        levels = octmh_material(wmax, rO.mWindingThreshold, 16);
    }
    std::vector<std::int32_t> material;
    OctmhShells shells;
    if (levels > 1) {
        shells = octmh_shells(rS, rO, levels);
        material.resize(num_grid);
        for (std::size_t h = 0; h < num_grid; ++h)
            material[h] = octmh_material(cells.mWinding[h], rO.mWindingThreshold, levels);
        if (rO.mHexGrading == OctreeMeshHexGrading::Graded) {
            // As at the skin: a centre within a tenth of its cell of its
            // interface (a shell on a lattice plane) is kept by a hair.
            for (std::size_t h = 0; h < num_grid; ++h) {
                if (material[h] < 2)
                    continue;
                const DistanceQuery& query = octmh_query(rS, &shells, material[h]);
                const double d = query_closest_points(query, {cells.mCentre[h]}).front().mDistance;
                const double cell =
                    rS.mFrame.mRootSize / static_cast<double>(std::int64_t{1} << cell_level[h]);
                if (d < 0.1 * cell)
                    --material[h];
            }
        }
        const std::int32_t empty = octmh_nest(hm, material, cells.mWinding, levels);
        if (empty > 0)
            octree_mesh_warn(
                rR,
                std::string(kOctreeMeshPrefix) + std::to_string(empty) +
                    " nested shells have no hexahedron between their "
                    "interfaces; lower " +
                    (rO.mHexGrading == OctreeMeshHexGrading::Graded ? "min_size" : "surface_size"));
    }
    // Per hexahedron, whether it is a pillow cell (grows with the pillows).
    std::vector<std::uint8_t> is_pillow(num_grid, 0);

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
        for (const auto& side : grid_boundary.mSides) {
            const auto h = static_cast<std::size_t>(side[0]);
            const std::int64_t lvl = cell_level[h];
            cell_level.push_back(lvl);
            is_pillow.push_back(1);
            if (!material.empty()) {
                const std::int32_t mat = material[h];
                material.push_back(mat);
            }
        }
    }
    const HexBoundary boundary = rO.mPillow ? hex_boundary(hm) : grid_boundary;
    // The features of shell @p Shell (0: every feature).
    auto features_of = [&](std::int32_t Shell) -> std::optional<FeatureGraph> {
        if (!rO.mPreserveFeatures || rS.mFeatures.empty())
            return std::nullopt;
        std::vector<std::array<std::int64_t, 2>> edges;
        for (const OctreeFeatureEdge& e : rS.mFeatures)
            if (Shell == 0 || (shells.mPointShell[static_cast<std::size_t>(e.mA)] == Shell &&
                               shells.mPointShell[static_cast<std::size_t>(e.mB)] == Shell))
                edges.push_back({e.mA, e.mB});
        if (edges.empty())
            return std::nullopt;
        return feature_graph(rS.mSkin.mPoints, edges, rO.mCornerAngle);
    };
    const std::optional<FeatureGraph> graph = features_of(levels > 1 ? 1 : 0);
    // The interfaces, each pillowed on both sides (cleaving).
    std::vector<OctmhInterface> interfaces;
    for (std::int32_t level = 2; level <= levels; ++level) {
        std::vector<std::uint8_t> inner(hm.mHexes.size(), 0);
        for (std::size_t h = 0; h < inner.size(); ++h)
            inner[h] = material[h] >= level;
        const std::vector<Vec3> before = hm.mPoints;
        const HexInterface pil = hex_pillow_interface(hm, inner);
        if (pil.mQuads.empty())
            continue;
        OctmhInterface itf;
        itf.mShell = level;
        const std::vector<Vec3> nrm =
            hex_boundary_normals(before, pil.mQuads, kOctmhNormalSmoothing);
        size.resize(hm.mPoints.size());
        for (std::size_t v = 0; v < before.size(); ++v)
            if (pil.mShell[v] >= 0) {
                itf.mNodes.push_back(static_cast<std::int64_t>(v));
                itf.mShellNode.push_back(pil.mShell[v]);
                itf.mInnerNode.push_back(pil.mInner[v]);
                itf.mNormals.push_back(nrm[v]);
                size[static_cast<std::size_t>(pil.mShell[v])] = size[v];
                size[static_cast<std::size_t>(pil.mInner[v])] = size[v];
            }
        for (const auto& q : pil.mQuads) {
            std::array<std::int64_t, 4> on{};
            for (std::size_t c = 0; c < 4; ++c)
                on[c] = pil.mShell[static_cast<std::size_t>(q[c])];
            itf.mQuads.push_back(on);
        }
        for (std::size_t f = 0; f < pil.mQuads.size(); ++f) {
            const std::int64_t lvl = cell_level[static_cast<std::size_t>(pil.mInnerHex[f])];
            for (std::int32_t layer : {level, level - 1}) {
                cell_level.push_back(lvl);
                is_pillow.push_back(1);
                material.push_back(layer);
            }
        }
        itf.mGraph = features_of(level);
        interfaces.push_back(std::move(itf));
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
                           graph ? &*graph : nullptr,
                           levels > 1 ? &shells : nullptr,
                           levels > 1 ? &interfaces : nullptr};
    const char* p_env = std::getenv("MESHIOPLUSPLUS_OCTREE_MESH_HEX_FIT");
    const std::string forced = p_env != nullptr ? std::string(p_env) : std::string();
    OctmhFit fit;
    if (levels > 1) {
        // Materials: the 18.0.0 fit had none to compare with.
        fit = octmh_fit(in, OctmhPlacement::Robust);
    } else if (forced == "body_fitted") {
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
    const bool any_features =
        graph || std::any_of(interfaces.begin(), interfaces.end(),
                             [](const OctmhInterface& rI) { return rI.mGraph.has_value(); });
    if (any_features) {
        const OctmhFeatures& f = fit.mFeatures;
        rR.mNumCorners = f.mCorners;
        rR.mNumCornersSnapped = f.mCornersSnapped;
        rR.mNumFeatureChains = f.mChains;
        rR.mNumFeatureChainsMissed = f.mChainsMissed;
        rR.mNumFeatureMoves = f.mMoved;
        if (f.mCornersSnapped + fit.mShortCorners < f.mCorners)
            octree_mesh_warn(
                rR, std::string(kOctreeMeshPrefix) +
                        std::to_string(f.mCorners - f.mCornersSnapped - fit.mShortCorners) +
                        " of " + std::to_string(f.mCorners) + " feature corners not recovered");
        if (f.mChainsMissed > fit.mShortChains)
            octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) +
                                     std::to_string(f.mChainsMissed - fit.mShortChains) + " of " +
                                     std::to_string(f.mChains) + " feature chains not recovered");
        if (fit.mShortChains > 0)
            octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(fit.mShortChains) +
                                     " of " + std::to_string(f.mChains) +
                                     " feature chains are shorter than a cell and were not "
                                     "snapped, nor the corners only they end at");
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
    if (worst <= 0.0) {
        const char* p_size =
            rO.mHexGrading == OctreeMeshHexGrading::Graded ? "min_size" : "surface_size";
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) +
                                 "untangling left a hexahedron with a non-positive corner "
                                 "Jacobian; raise untangle_iterations or lower " +
                                 p_size);
        // Where the skin is thinner than two cells, no fit within the size
        // floor can help: say how many of those hexahedra lie there.
        std::vector<Vec3> at;
        std::vector<double> cell;
        for (std::size_t h = 0; h < hm.mHexes.size(); ++h) {
            if (hex_min_scaled_jacobian(hm.mPoints, hm.mHexes[h]) > 0.0)
                continue;
            Vec3 c{0.0, 0.0, 0.0};
            for (std::int64_t v : hm.mHexes[h])
                c = vec3_add(c, hm.mPoints[static_cast<std::size_t>(v)]);
            at.push_back(vec3_scale(c, 0.125));
            cell.push_back(rS.mFrame.mRootSize /
                           static_cast<double>(std::int64_t{1} << cell_level[h]));
        }
        const std::vector<SurfaceProjection> near = query_surface_projections(rS.mQuery, at);
        std::int64_t thin = 0;
        for (std::size_t k = 0; k < at.size(); ++k) {
            const Vec3 out_dir = vec3_scale(vec3_normalize(near[k].mNormal), rS.mWindingSign);
            thin += octmh_thickness(rS, near[k].mPoint, out_dir, 2.0 * cell[k]) < 2.0 * cell[k];
        }
        if (thin > 0)
            octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(thin) +
                                     " of the " + std::to_string(at.size()) +
                                     " hexahedra with a non-positive corner lie where the skin "
                                     "is thinner than two cells; lower " +
                                     p_size + " to resolve it");
    } else if (worst < 0.2)
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
        lv.As<std::int64_t>()[c] = cell_level[c];
        pillow.As<std::int64_t>()[c] = is_pillow[c];
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
    if (levels > 1) {
        NDArray mt = NDArray::Uninit(DType::Int64, {hm.mHexes.size()});
        std::vector<std::vector<std::int64_t>> by_material(static_cast<std::size_t>(levels) + 1);
        for (std::size_t c = 0; c < hm.mHexes.size(); ++c) {
            mt.As<std::int64_t>()[c] = material[c];
            by_material[static_cast<std::size_t>(material[c])].push_back(
                static_cast<std::int64_t>(c));
        }
        std::vector<NDArray> a;
        a.push_back(std::move(mt));
        m.AddCellData(kOctreeMeshMaterialName, std::move(a));
        for (std::int32_t mat = 1; mat <= levels; ++mat)
            m.AddRegion(Region("material_" + std::to_string(mat), RegionKind::Cell, 3, mat,
                               octmh_ids(by_material[static_cast<std::size_t>(mat)], 1)));
    }
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
    // The surfaces that stand for the skin: the boundary, and under cleaving
    // the interfaces, which stand for the nested shells.
    std::vector<std::array<std::int64_t, 4>> surface_quads = boundary.mQuads;
    for (const OctmhInterface& itf : interfaces)
        surface_quads.insert(surface_quads.end(), itf.mQuads.begin(), itf.mQuads.end());
    for (const auto& q : surface_quads) {
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
