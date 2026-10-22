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
#include "meshioplusplus/log.hpp"

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

// The arc length along a polyline of its nearest point to `rP`.
double featsnap_arc(const std::vector<Vec3>& rSkin, const std::vector<std::int64_t>& rChain,
                    const Vec3& rP) {
    double best_d = std::numeric_limits<double>::infinity(), best_s = 0.0, s = 0.0;
    for (std::size_t i = 0; i + 1 < rChain.size(); ++i) {
        const Vec3& a = rSkin[static_cast<std::size_t>(rChain[i])];
        const Vec3& b = rSkin[static_cast<std::size_t>(rChain[i + 1])];
        const Vec3 ab = vec3_sub(b, a);
        const double len2 = vec3_norm_sq(ab), len = std::sqrt(len2);
        const double t =
            len2 > 0.0 ? std::clamp(vec3_dot(vec3_sub(rP, a), ab) / len2, 0.0, 1.0) : 0.0;
        const double d = vec3_norm(vec3_sub(rP, vec3_add(a, vec3_scale(ab, t))));
        if (d < best_d) {
            best_d = d;
            best_s = s + t * len;
        }
        s += len;
    }
    return best_s;
}

// The point of a polyline at arc length `S`.
Vec3 featsnap_at_arc(const std::vector<Vec3>& rSkin, const std::vector<std::int64_t>& rChain,
                     double S) {
    double s = 0.0;
    for (std::size_t i = 0; i + 1 < rChain.size(); ++i) {
        const Vec3& a = rSkin[static_cast<std::size_t>(rChain[i])];
        const Vec3& b = rSkin[static_cast<std::size_t>(rChain[i + 1])];
        const double len = vec3_norm(vec3_sub(b, a));
        if (S <= s + len || i + 2 == rChain.size()) {
            const double t = len > 0.0 ? std::clamp((S - s) / len, 0.0, 1.0) : 0.0;
            return vec3_add(a, vec3_scale(vec3_sub(b, a), t));
        }
        s += len;
    }
    return rSkin[static_cast<std::size_t>(rChain.back())];
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
    for (auto& a : badj)
        std::sort(a.begin(), a.end());
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
    // Whether every tetrahedron around `U` stays positive with `V` at `rAt`.
    auto star_positive = [&](std::int64_t U, std::int64_t V, const Vec3& rAt) {
        for (std::int64_t t : inc[static_cast<std::size_t>(U)]) {
            Vec3 p[4];
            for (std::size_t c = 0; c < 4; ++c) {
                const std::int64_t q = rTets[static_cast<std::size_t>(t)][c];
                p[c] = q == V ? rAt : rPoints[static_cast<std::size_t>(q)];
            }
            if (featsnap_volume(p[0], p[1], p[2], p[3]) <= MinVolume)
                return false;
        }
        return true;
    };
    auto try_move = [&](std::int64_t V, const Vec3& rTarget) {
        bool positive = true;
        const double before = min_dihedral(V, rPoints[static_cast<std::size_t>(V)], positive);
        double after = min_dihedral(V, rTarget, positive);
        if (!positive) {
            // Make room: move the free interior nodes of the tetrahedra the
            // move would invert to the centroid of their stars.
            std::vector<std::pair<std::int64_t, Vec3>> undo;
            for (int pass = 0; pass < 2 && !positive; ++pass) {
                for (std::int64_t t : inc[static_cast<std::size_t>(V)])
                    for (std::int64_t u : rTets[static_cast<std::size_t>(t)]) {
                        const auto ui = static_cast<std::size_t>(u);
                        if (u == V || boundary[ui] || rPinned[ui])
                            continue;
                        Vec3 c{0.0, 0.0, 0.0};
                        double w = 0.0;
                        for (std::int64_t s : inc[ui])
                            for (std::int64_t q : rTets[static_cast<std::size_t>(s)])
                                if (q != u) {
                                    c = vec3_add(
                                        c, q == V ? rTarget : rPoints[static_cast<std::size_t>(q)]);
                                    w += 1.0;
                                }
                        const Vec3 old = rPoints[ui];
                        rPoints[ui] = vec3_scale(c, 1.0 / w);
                        if (star_positive(u, V, rTarget) &&
                            star_positive(u, V, rPoints[static_cast<std::size_t>(V)]))
                            undo.emplace_back(u, old);
                        else
                            rPoints[ui] = old;
                    }
                positive = true;
                after = min_dihedral(V, rTarget, positive);
            }
            if (!positive || after < std::min(before, kFeatsnapFloorDeg)) {
                for (auto it = undo.rbegin(); it != undo.rend(); ++it)
                    rPoints[static_cast<std::size_t>(it->first)] = it->second;
                return false;
            }
        } else if (after < std::min(before, kFeatsnapFloorDeg)) {
            return false;
        }
        rPoints[static_cast<std::size_t>(V)] = rTarget;
        ++report.mNumMoved;
        return true;
    };
    // A move in four steps, each making room around it as it goes; a refused
    // step keeps the steps before it.
    auto move_by_steps = [&](std::int64_t V, const Vec3& rTarget) {
        const Vec3 at = rPoints[static_cast<std::size_t>(V)];
        for (int k = 1; k <= 4; ++k)
            if (!try_move(V, vec3_add(at, vec3_scale(vec3_sub(rTarget, at), 0.25 * k))))
                return false;
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
        for (const auto& [d, v] : cand) {
            const Vec3 was = rPoints[static_cast<std::size_t>(v)];
            if (!move_by_steps(v, cp)) {
                rPoints[static_cast<std::size_t>(v)] = was;
                continue;
            }
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
            log::debug("feature_snap: a chain's ends are not snapped ({}, {})", from, to);
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
                const double w = vec3_norm(vec3_sub(pu, pv)) * (1.0 + 16.0 * (du + dv) / h);
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
            log::debug("feature_snap: no path for a chain of {} points", chain.size());
            ++report.mNumChainsMissed;
            continue;
        }
        std::vector<std::int64_t> path;
        for (std::int64_t v = prev[static_cast<std::size_t>(to)]; v != from;
             v = prev[static_cast<std::size_t>(v)])
            path.push_back(v);
        std::reverse(path.begin(), path.end());
        // A node whose neighbours on the path share an edge is a detour
        // through a triangle (projecting all three onto the chain would make
        // the triangle flat): drop it from the path.
        {
            std::vector<std::int64_t> full{from};
            full.insert(full.end(), path.begin(), path.end());
            full.push_back(to);
            for (bool changed = true; changed;) {
                changed = false;
                for (std::size_t i = 1; i + 1 < full.size(); ++i) {
                    const auto& a = badj[static_cast<std::size_t>(full[i - 1])];
                    if (std::binary_search(a.begin(), a.end(), full[i + 1])) {
                        full.erase(full.begin() + static_cast<std::ptrdiff_t>(i));
                        changed = true;
                        break;
                    }
                }
            }
            path.assign(full.begin() + 1, full.end() - 1);
        }
        // Targets in order along the chain: each node's nearest arc length,
        // made increasing; when two land closer than a fifth of the local
        // size, or out of order, the nodes are spread evenly instead.
        double total = 0.0;
        for (std::size_t i = 0; i + 1 < chain.size(); ++i)
            total += vec3_norm(vec3_sub(rSkinPoints[static_cast<std::size_t>(chain[i + 1])],
                                        rSkinPoints[static_cast<std::size_t>(chain[i])]));
        std::vector<double> arc(path.size());
        bool even = false;
        double last = 0.0;
        for (std::size_t i = 0; i < path.size(); ++i) {
            const Vec3& pv = rPoints[static_cast<std::size_t>(path[i])];
            arc[i] = featsnap_arc(rSkinPoints, chain, pv);
            const double gap = 0.2 * rLocalSize(pv);
            if (arc[i] < last + gap || arc[i] > total - gap)
                even = true;
            last = arc[i];
        }
        if (even)
            for (std::size_t i = 0; i < path.size(); ++i)
                arc[i] = total * static_cast<double>(i + 1) / static_cast<double>(path.size() + 1);
        bool all = true;
        std::int64_t refused = 0;
        for (std::size_t i = 0; i < path.size(); ++i) {
            const std::int64_t v = path[i];
            const Vec3 target = featsnap_at_arc(rSkinPoints, chain, arc[i]);
            const Vec3 at = rPoints[static_cast<std::size_t>(v)];
            (void)at;
            const bool moved = move_by_steps(v, target);
            all = all && moved;
            refused += !moved;
            rPinned[static_cast<std::size_t>(v)] = 1;
            rFeature[static_cast<std::size_t>(v)] = 2;
        }
        if (!all) {
            log::debug("feature_snap: a chain of {} points ({} nodes) kept {} off it", chain.size(),
                       path.size(), refused);
            ++report.mNumChainsMissed;
        }
    }
    return report;
}

}  // namespace detail
}  // namespace meshioplusplus
