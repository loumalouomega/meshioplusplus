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
// The feature graph and the feature snapping of detail/feature_snap.hpp.
// Anonymous-namespace helpers are prefixed `featsnap_` (the amalgamation
// concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

// Project includes
#include "bcc_octree_tets.hpp"
#include "feature_snap.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kFeatsnapPi = 3.14159265358979323846;
constexpr double kFeatsnapFloorDeg = 2.0;

double featsnap_volume(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) {
    return triple_product(vec3_sub(rB, rA), vec3_sub(rC, rA), vec3_sub(rD, rA)) / 6.0;
}

// The nearest point of a polyline (skin point ids) to `rP`, and its distance.
std::pair<Vec3, double> featsnap_on_polyline(const std::vector<Vec3>& rSkin,
                                             const std::vector<std::int64_t>& rChain,
                                             const Vec3& rP) {
    Vec3 best = rSkin[static_cast<std::size_t>(rChain[0])];
    double best_d = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i + 1 < rChain.size(); ++i) {
        const Vec3& a = rSkin[static_cast<std::size_t>(rChain[i])];
        const Vec3& b = rSkin[static_cast<std::size_t>(rChain[i + 1])];
        const Vec3 ab = vec3_sub(b, a);
        const double len2 = vec3_norm_sq(ab);
        const double t =
            len2 > 0.0 ? std::clamp(vec3_dot(vec3_sub(rP, a), ab) / len2, 0.0, 1.0) : 0.0;
        const Vec3 q = vec3_add(a, vec3_scale(ab, t));
        const double d = vec3_norm(vec3_sub(rP, q));
        if (d < best_d) {
            best_d = d;
            best = q;
        }
    }
    return {best, best_d};
}

}  // namespace

FeatureGraph feature_graph(const std::vector<Vec3>& rPoints,
                           const std::vector<std::array<std::int64_t, 2>>& rEdges,
                           double CornerAngleDeg) {
    FeatureGraph g;
    const std::size_t np = rPoints.size();
    std::vector<std::vector<std::int64_t>> adj(np);
    std::vector<std::array<std::int64_t, 2>> edges;
    for (const auto& e : rEdges)
        if (e[0] != e[1])
            edges.push_back({std::min(e[0], e[1]), std::max(e[0], e[1])});
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    for (const auto& e : edges) {
        adj[static_cast<std::size_t>(e[0])].push_back(e[1]);
        adj[static_cast<std::size_t>(e[1])].push_back(e[0]);
    }
    for (auto& a : adj)
        std::sort(a.begin(), a.end());
    std::vector<std::uint8_t> corner(np, 0);
    for (std::size_t p = 0; p < np; ++p) {
        const auto& a = adj[p];
        if (a.empty())
            continue;
        if (a.size() != 2) {
            corner[p] = 1;
            continue;
        }
        const Vec3 u = vec3_sub(rPoints[static_cast<std::size_t>(a[0])], rPoints[p]);
        const Vec3 v = vec3_sub(rPoints[static_cast<std::size_t>(a[1])], rPoints[p]);
        const double c = vec3_dot(u, v) / std::sqrt(vec3_norm_sq(u) * vec3_norm_sq(v));
        const double turn = 180.0 - std::acos(std::clamp(c, -1.0, 1.0)) * 180.0 / kFeatsnapPi;
        corner[p] = turn > CornerAngleDeg;
    }
    auto edge_index = [&](std::int64_t A, std::int64_t B) {
        const std::array<std::int64_t, 2> k{std::min(A, B), std::max(A, B)};
        return static_cast<std::size_t>(std::lower_bound(edges.begin(), edges.end(), k) -
                                        edges.begin());
    };
    std::vector<std::uint8_t> used(edges.size(), 0);
    auto walk = [&](std::int64_t Start, std::int64_t Next) {
        std::vector<std::int64_t> chain{Start, Next};
        used[edge_index(Start, Next)] = 1;
        std::int64_t prev = Start, cur = Next;
        while (!corner[static_cast<std::size_t>(cur)]) {
            const auto& a = adj[static_cast<std::size_t>(cur)];
            const std::int64_t nxt = a[0] == prev ? a[1] : a[0];
            if (used[edge_index(cur, nxt)])
                break;
            used[edge_index(cur, nxt)] = 1;
            chain.push_back(nxt);
            prev = cur;
            cur = nxt;
        }
        return chain;
    };
    for (std::size_t p = 0; p < np; ++p)
        if (corner[p])
            for (std::int64_t n : adj[p])
                if (!used[edge_index(static_cast<std::int64_t>(p), n)])
                    g.mChains.push_back(walk(static_cast<std::int64_t>(p), n));
    // Closed loops with no corner: cut at the lowest vertex and half way round.
    for (std::size_t e = 0; e < edges.size(); ++e) {
        if (used[e])
            continue;
        const std::int64_t v0 = edges[e][0];
        corner[static_cast<std::size_t>(v0)] = 1;
        std::vector<std::int64_t> loop = walk(v0, edges[e][1]);
        if (loop.size() < 3 || loop.back() != v0) {
            g.mChains.push_back(std::move(loop));
            continue;
        }
        const std::size_t mid = loop.size() / 2;
        corner[static_cast<std::size_t>(loop[mid])] = 1;
        g.mChains.emplace_back(loop.begin(), loop.begin() + static_cast<std::ptrdiff_t>(mid) + 1);
        g.mChains.emplace_back(loop.begin() + static_cast<std::ptrdiff_t>(mid), loop.end());
    }
    for (std::size_t p = 0; p < np; ++p)
        if (corner[p])
            g.mCorners.push_back(static_cast<std::int64_t>(p));
    return g;
}

FeatureSnapReport feature_snap(std::vector<Vec3>& rPoints,
                               const std::vector<std::array<std::int64_t, 4>>& rTets,
                               const std::vector<std::uint8_t>& rInside, const FeatureGraph& rGraph,
                               const std::vector<Vec3>& rSkinPoints,
                               const std::function<double(const Vec3&)>& rLocalSize,
                               double MinVolume, std::vector<std::uint8_t>& rPinned,
                               std::vector<std::uint8_t>& rFeature) {
    FeatureSnapReport report;
    report.mNumCorners = static_cast<std::int64_t>(rGraph.mCorners.size());
    report.mNumChains = static_cast<std::int64_t>(rGraph.mChains.size());
    const std::size_t np = rPoints.size();
    rFeature.resize(np, 0);
    // The boundary faces of the inside tetrahedra, their vertices and edges.
    constexpr int kFace[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};
    std::vector<std::array<std::int64_t, 3>> faces;
    for (std::size_t t = 0; t < rTets.size(); ++t) {
        if (!rInside[t])
            continue;
        for (const auto& f : kFace) {
            std::array<std::int64_t, 3> k{rTets[t][f[0]], rTets[t][f[1]], rTets[t][f[2]]};
            std::sort(k.begin(), k.end());
            faces.push_back(k);
        }
    }
    std::sort(faces.begin(), faces.end());
    std::vector<std::uint8_t> boundary(np, 0);
    std::vector<std::array<std::int64_t, 2>> bedges;
    for (std::size_t i = 0; i < faces.size();) {
        std::size_t j = i;
        while (j < faces.size() && faces[j] == faces[i])
            ++j;
        if (j - i == 1) {
            const auto& f = faces[i];
            for (std::int64_t v : f)
                boundary[static_cast<std::size_t>(v)] = 1;
            bedges.push_back({f[0], f[1]});
            bedges.push_back({f[1], f[2]});
            bedges.push_back({f[0], f[2]});
        }
        i = j;
    }
    std::sort(bedges.begin(), bedges.end());
    bedges.erase(std::unique(bedges.begin(), bedges.end()), bedges.end());
    std::vector<std::vector<std::int64_t>> badj(np);
    for (const auto& e : bedges) {
        badj[static_cast<std::size_t>(e[0])].push_back(e[1]);
        badj[static_cast<std::size_t>(e[1])].push_back(e[0]);
    }
    // Every tetrahedron around each point.
    std::vector<std::vector<std::int64_t>> inc(np);
    for (std::size_t t = 0; t < rTets.size(); ++t)
        for (std::int64_t v : rTets[t])
            inc[static_cast<std::size_t>(v)].push_back(static_cast<std::int64_t>(t));
    auto min_dihedral = [&](std::int64_t V, const Vec3& rAt, bool& rPositive) {
        double lo = 180.0;
        rPositive = true;
        for (std::int64_t t : inc[static_cast<std::size_t>(V)]) {
            Vec3 p[4];
            for (std::size_t c = 0; c < 4; ++c) {
                const std::int64_t q = rTets[static_cast<std::size_t>(t)][c];
                p[c] = q == V ? rAt : rPoints[static_cast<std::size_t>(q)];
            }
            if (featsnap_volume(p[0], p[1], p[2], p[3]) <= MinVolume) {
                rPositive = false;
                return 0.0;
            }
            lo = std::min(lo, tet_dihedral_angles(p[0], p[1], p[2], p[3])[0]);
        }
        return lo;
    };
    auto try_move = [&](std::int64_t V, const Vec3& rTarget) {
        bool positive = true;
        const double before = min_dihedral(V, rPoints[static_cast<std::size_t>(V)], positive);
        const double after = min_dihedral(V, rTarget, positive);
        if (!positive || after < std::min(before, kFeatsnapFloorDeg))
            return false;
        rPoints[static_cast<std::size_t>(V)] = rTarget;
        ++report.mNumMoved;
        return true;
    };
    // Corners: the nearest boundary vertex within reach that can take it.
    std::vector<std::int64_t> bverts;
    for (std::size_t p = 0; p < np; ++p)
        if (boundary[p])
            bverts.push_back(static_cast<std::int64_t>(p));
    std::vector<std::int64_t> snapped_at(rSkinPoints.size(), -1);
    for (std::int64_t c : rGraph.mCorners) {
        const Vec3& cp = rSkinPoints[static_cast<std::size_t>(c)];
        const double reach = 1.5 * rLocalSize(cp);
        std::vector<std::pair<double, std::int64_t>> cand;
        for (std::int64_t v : bverts) {
            if (rPinned[static_cast<std::size_t>(v)])
                continue;
            const double d = vec3_norm(vec3_sub(rPoints[static_cast<std::size_t>(v)], cp));
            if (d <= reach)
                cand.emplace_back(d, v);
        }
        std::sort(cand.begin(), cand.end());
        for (const auto& [d, v] : cand)
            if (try_move(v, cp)) {
                rPinned[static_cast<std::size_t>(v)] = 1;
                rFeature[static_cast<std::size_t>(v)] = 1;
                snapped_at[static_cast<std::size_t>(c)] = v;
                ++report.mNumCornersSnapped;
                break;
            }
    }
    // Chains: a shortest boundary path in a tube around the polyline.
    for (const auto& chain : rGraph.mChains) {
        const std::int64_t from = snapped_at[static_cast<std::size_t>(chain.front())];
        const std::int64_t to = snapped_at[static_cast<std::size_t>(chain.back())];
        if (from < 0 || to < 0 || from == to) {
            ++report.mNumChainsMissed;
            continue;
        }
        std::vector<double> dist(np, std::numeric_limits<double>::infinity());
        std::vector<std::int64_t> prev(np, -1);
        using Item = std::pair<double, std::int64_t>;
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
        dist[static_cast<std::size_t>(from)] = 0.0;
        heap.emplace(0.0, from);
        while (!heap.empty()) {
            const auto [d, v] = heap.top();
            heap.pop();
            if (d > dist[static_cast<std::size_t>(v)] || v == to)
                continue;
            const Vec3& pv = rPoints[static_cast<std::size_t>(v)];
            const double dv = featsnap_on_polyline(rSkinPoints, chain, pv).second;
            for (std::int64_t u : badj[static_cast<std::size_t>(v)]) {
                if (u != to && rPinned[static_cast<std::size_t>(u)])
                    continue;
                const Vec3& pu = rPoints[static_cast<std::size_t>(u)];
                const double h = rLocalSize(pu);
                const double du = featsnap_on_polyline(rSkinPoints, chain, pu).second;
                if (u != to && du > 1.25 * h)
                    continue;
                const double w = vec3_norm(vec3_sub(pu, pv)) * (1.0 + 4.0 * (du + dv) / h);
                const double nd = d + w;
                const auto ui = static_cast<std::size_t>(u);
                if (nd < dist[ui] || (nd == dist[ui] && v < prev[ui])) {
                    dist[ui] = nd;
                    prev[ui] = v;
                    heap.emplace(nd, u);
                }
            }
        }
        if (prev[static_cast<std::size_t>(to)] < 0) {
            ++report.mNumChainsMissed;
            continue;
        }
        std::vector<std::int64_t> path;
        for (std::int64_t v = prev[static_cast<std::size_t>(to)]; v != from;
             v = prev[static_cast<std::size_t>(v)])
            path.push_back(v);
        std::reverse(path.begin(), path.end());
        bool all = true;
        for (std::int64_t v : path) {
            const Vec3 target =
                featsnap_on_polyline(rSkinPoints, chain, rPoints[static_cast<std::size_t>(v)])
                    .first;
            const Vec3 at = rPoints[static_cast<std::size_t>(v)];
            bool moved = false;
            for (double f : {1.0, 0.5, 0.25})
                if (try_move(v, vec3_add(at, vec3_scale(vec3_sub(target, at), f)))) {
                    moved = f == 1.0;
                    break;
                }
            all = all && moved;
            rPinned[static_cast<std::size_t>(v)] = 1;
            rFeature[static_cast<std::size_t>(v)] = 2;
        }
        if (!all)
            ++report.mNumChainsMissed;
    }
    return report;
}

}  // namespace detail
}  // namespace meshioplusplus
