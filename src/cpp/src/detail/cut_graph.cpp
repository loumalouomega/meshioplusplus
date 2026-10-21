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

// System includes
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <utility>

// Project includes
#include "cut_graph.hpp"
#include "param_disc.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

using CgEdge = std::pair<std::size_t, std::size_t>;

CgEdge cg_edge(std::size_t A, std::size_t B) {
    return {std::min(A, B), std::max(A, B)};
}

double cg_dist(const PatchVec& rA, const PatchVec& rB) {
    const double dx = rA[0] - rB[0], dy = rA[1] - rB[1], dz = rA[2] - rB[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Per vertex its neighbours, ascending (indices up to the largest one used).
std::vector<std::vector<std::size_t>> cg_neighbours(const std::vector<PatchTri>& rTris,
                                                    std::size_t NumPts) {
    std::vector<std::vector<std::size_t>> nb(NumPts);
    for (const PatchTri& t : rTris)
        for (int k = 0; k < 3; ++k) {
            nb[t[static_cast<std::size_t>(k)]].push_back(t[static_cast<std::size_t>((k + 1) % 3)]);
            nb[t[static_cast<std::size_t>((k + 1) % 3)]].push_back(t[static_cast<std::size_t>(k)]);
        }
    for (auto& v : nb) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    return nb;
}

std::size_t cg_num_points(const std::vector<PatchTri>& rTris) {
    std::size_t n = 0;
    for (const PatchTri& t : rTris)
        n = std::max({n, t[0] + 1, t[1] + 1, t[2] + 1});
    return n;
}

// Dijkstra from every source at once; ties go to the lower vertex.
void cg_dijkstra(const std::vector<PatchVec>& rPts,
                 const std::vector<std::vector<std::size_t>>& rNb,
                 const std::vector<std::size_t>& rSources,
                 const std::vector<std::uint8_t>& rBlocked, std::vector<double>& rDist,
                 std::vector<std::size_t>& rPrev) {
    const std::size_t n = rNb.size();
    rDist.assign(n, std::numeric_limits<double>::infinity());
    rPrev.assign(n, n);
    using Item = std::pair<double, std::size_t>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
    for (std::size_t s : rSources)
        if (s < n && rDist[s] != 0.0) {
            rDist[s] = 0.0;
            queue.emplace(0.0, s);
        }
    while (!queue.empty()) {
        const auto [d, v] = queue.top();
        queue.pop();
        if (d > rDist[v])
            continue;
        for (std::size_t w : rNb[v]) {
            if (w < rBlocked.size() && rBlocked[w])
                continue;
            const double nd = d + cg_dist(rPts[v], rPts[w]);
            if (nd < rDist[w]) {
                rDist[w] = nd;
                rPrev[w] = v;
                queue.emplace(nd, w);
            }
        }
    }
}

// Per triangle its edge-connected component across the edges not in `rCut`
// (and not on the boundary), numbered by lowest triangle.
std::vector<std::int64_t> cg_components(const std::vector<PatchTri>& rTris,
                                        const std::vector<CgEdge>& rCut, std::size_t& rNum) {
    std::map<CgEdge, std::vector<std::size_t>> edge_tris;
    for (std::size_t t = 0; t < rTris.size(); ++t)
        for (int k = 0; k < 3; ++k)
            edge_tris[cg_edge(rTris[t][static_cast<std::size_t>(k)],
                              rTris[t][static_cast<std::size_t>((k + 1) % 3)])]
                .push_back(t);
    std::vector<std::int64_t> label(rTris.size(), -1);
    rNum = 0;
    std::vector<std::size_t> stack;
    for (std::size_t seed = 0; seed < rTris.size(); ++seed) {
        if (label[seed] >= 0)
            continue;
        const std::int64_t id = static_cast<std::int64_t>(rNum++);
        label[seed] = id;
        stack.assign(1, seed);
        while (!stack.empty()) {
            const std::size_t t = stack.back();
            stack.pop_back();
            for (int k = 0; k < 3; ++k) {
                const CgEdge e = cg_edge(rTris[t][static_cast<std::size_t>(k)],
                                         rTris[t][static_cast<std::size_t>((k + 1) % 3)]);
                if (std::binary_search(rCut.begin(), rCut.end(), e))
                    continue;
                for (std::size_t o : edge_tris[e])
                    if (label[o] < 0) {
                        label[o] = id;
                        stack.push_back(o);
                    }
            }
        }
    }
    return label;
}

// Arc length positions along a closed vertex loop, and its length.
std::vector<double> cg_arc(const std::vector<PatchVec>& rPts, const std::vector<std::size_t>& rLoop,
                           double& rLength) {
    std::vector<double> s(rLoop.size(), 0.0);
    rLength = 0.0;
    for (std::size_t i = 0; i < rLoop.size(); ++i) {
        s[i] = rLength;
        rLength += cg_dist(rPts[rLoop[i]], rPts[rLoop[(i + 1) % rLoop.size()]]);
    }
    return s;
}

// The loop's vertex half way round (by arc length) from `From`.
std::size_t cg_opposite(const std::vector<PatchVec>& rPts, const std::vector<std::size_t>& rLoop,
                        std::size_t From) {
    double len = 0.0;
    const std::vector<double> s = cg_arc(rPts, rLoop, len);
    const std::size_t i0 =
        static_cast<std::size_t>(std::find(rLoop.begin(), rLoop.end(), From) - rLoop.begin());
    if (i0 >= rLoop.size())
        return rLoop.front();
    std::size_t best = i0;
    double best_gap = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < rLoop.size(); ++i) {
        double d = std::fabs(s[i] - s[i0]);
        d = std::min(d, len - d);
        const double gap = std::fabs(d - 0.5 * len);
        if (gap < best_gap) {
            best_gap = gap;
            best = i;
        }
    }
    return rLoop[best];
}

// Every piece a disc of at least `MinTriangles` triangles.
bool cg_pieces_are_discs(const std::vector<PatchTri>& rTris,
                         const std::vector<std::int64_t>& rLabel, std::size_t NumPieces,
                         std::size_t MinTriangles) {
    std::vector<std::vector<PatchTri>> piece(NumPieces);
    for (std::size_t t = 0; t < rTris.size(); ++t)
        piece[static_cast<std::size_t>(rLabel[t])].push_back(rTris[t]);
    for (const auto& p : piece)
        if (p.size() < MinTriangles || !patch_topology(p).IsDisc())
            return false;
    return true;
}

}  // namespace

PatchTopology patch_topology(const std::vector<PatchTri>& rTris) {
    PatchTopology topo;
    topo.mNumTriangles = rTris.size();
    std::map<CgEdge, int> edges;
    std::vector<std::size_t> verts;
    for (const PatchTri& t : rTris)
        for (int k = 0; k < 3; ++k) {
            verts.push_back(t[static_cast<std::size_t>(k)]);
            ++edges[cg_edge(t[static_cast<std::size_t>(k)],
                            t[static_cast<std::size_t>((k + 1) % 3)])];
        }
    std::sort(verts.begin(), verts.end());
    verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
    topo.mNumVertices = verts.size();
    topo.mNumEdges = edges.size();
    for (const auto& kv : edges)
        if (kv.second > 2)
            topo.mManifold = false;
    bool ok = true;
    topo.mNumBoundaryLoops = param_boundary_loops(rTris, ok).size();
    topo.mManifold = topo.mManifold && ok;
    topo.mEuler = static_cast<std::int64_t>(topo.mNumVertices) -
                  static_cast<std::int64_t>(topo.mNumEdges) +
                  static_cast<std::int64_t>(topo.mNumTriangles);
    topo.mGenus = (2 - topo.mEuler - static_cast<std::int64_t>(topo.mNumBoundaryLoops)) / 2;
    return topo;
}

std::vector<std::size_t> patch_shortest_path(const std::vector<PatchVec>& rPts,
                                             const std::vector<PatchTri>& rTris,
                                             const std::vector<std::size_t>& rSources,
                                             const std::vector<std::size_t>& rTargets,
                                             const std::vector<std::uint8_t>& rForbidden) {
    const std::size_t n = std::max(cg_num_points(rTris), rPts.size());
    const auto nb = cg_neighbours(rTris, n);
    std::vector<std::uint8_t> blocked(n, 0);
    for (std::size_t v = 0; v < std::min(n, rForbidden.size()); ++v)
        blocked[v] = rForbidden[v];
    for (std::size_t v : rSources)
        if (v < n)
            blocked[v] = 0;
    for (std::size_t v : rTargets)
        if (v < n)
            blocked[v] = 0;
    std::vector<double> dist;
    std::vector<std::size_t> prev;
    cg_dijkstra(rPts, nb, rSources, blocked, dist, prev);
    std::size_t best = n;
    for (std::size_t t : rTargets)
        if (t < n && std::isfinite(dist[t]) &&
            (best == n || dist[t] < dist[best] || (dist[t] == dist[best] && t < best)))
            best = t;
    std::vector<std::size_t> path;
    for (std::size_t v = best; v != n; v = prev[v])
        path.push_back(v);
    std::reverse(path.begin(), path.end());
    return path;
}

std::vector<std::int64_t> patch_split_to_discs(const std::vector<PatchVec>& rPts,
                                               const std::vector<PatchTri>& rTris,
                                               std::size_t MinTriangles, std::size_t& rNumPieces) {
    rNumPieces = 0;
    const PatchTopology topo = patch_topology(rTris);
    if (!topo.mManifold || topo.mGenus != 0 || topo.IsDisc() || rTris.size() < 2 * MinTriangles)
        return {};
    const std::size_t n = std::max(cg_num_points(rTris), rPts.size());
    std::vector<CgEdge> cut;

    if (topo.mNumBoundaryLoops == 0) {
        // Two caps: the median level of the distance difference to the two
        // most distant vertices (the lowest-indexed vertex, the farthest from
        // it, and the farthest from that).
        const auto nb = cg_neighbours(rTris, n);
        std::size_t a = n;
        for (const PatchTri& t : rTris)
            a = std::min({a, t[0], t[1], t[2]});
        std::vector<double> da, db;
        std::vector<std::size_t> prev;
        const std::vector<std::uint8_t> none;
        const auto farthest = [&](const std::vector<double>& rD) {
            std::size_t f = a;
            for (std::size_t v = 0; v < n; ++v)
                if (std::isfinite(rD[v]) && rD[v] > rD[f])
                    f = v;
            return f;
        };
        cg_dijkstra(rPts, nb, {a}, none, da, prev);
        a = farthest(da);
        cg_dijkstra(rPts, nb, {a}, none, da, prev);
        const std::size_t b = farthest(da);
        cg_dijkstra(rPts, nb, {b}, none, db, prev);
        std::vector<double> key(rTris.size()), area(rTris.size());
        for (std::size_t t = 0; t < rTris.size(); ++t) {
            const PatchTri& tri = rTris[t];
            key[t] =
                (da[tri[0]] + da[tri[1]] + da[tri[2]] - db[tri[0]] - db[tri[1]] - db[tri[2]]) / 3.0;
            const PatchVec& p = rPts[tri[0]];
            const PatchVec& q = rPts[tri[1]];
            const PatchVec& r = rPts[tri[2]];
            const double u[3] = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
            const double v[3] = {r[0] - p[0], r[1] - p[1], r[2] - p[2]};
            const double c[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                 u[0] * v[1] - u[1] * v[0]};
            area[t] = 0.5 * std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        }
        std::vector<std::size_t> order(rTris.size());
        for (std::size_t t = 0; t < order.size(); ++t)
            order[t] = t;
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t A, std::size_t B) { return key[A] < key[B]; });
        double total = 0.0;
        for (double x : area)
            total += x;
        std::vector<std::uint8_t> side(rTris.size(), 1);
        double acc = 0.0;
        for (std::size_t t : order) {
            if (acc >= 0.5 * total)
                break;
            side[t] = 0;
            acc += area[t];
        }
        // The cut: every edge between the two sides.
        std::map<CgEdge, std::vector<std::size_t>> edge_tris;
        for (std::size_t t = 0; t < rTris.size(); ++t)
            for (int k = 0; k < 3; ++k)
                edge_tris[cg_edge(rTris[t][static_cast<std::size_t>(k)],
                                  rTris[t][static_cast<std::size_t>((k + 1) % 3)])]
                    .push_back(t);
        for (const auto& kv : edge_tris)
            if (kv.second.size() == 2 && side[kv.second[0]] != side[kv.second[1]])
                cut.push_back(kv.first);
    } else {
        // A cycle of shortest paths through every loop.
        bool ok = true;
        const std::vector<std::vector<std::size_t>> loops = param_boundary_loops(rTris, ok);
        const std::size_t b = loops.size();
        std::vector<std::uint8_t> forbidden(n, 0);
        std::size_t start = n, arrived = n;
        for (std::size_t i = 0; i < b; ++i) {
            const std::vector<std::size_t>& from = loops[i];
            const std::vector<std::size_t>& to = loops[(i + 1) % b];
            std::vector<std::size_t> sources;
            if (i == 0)
                sources = from;
            else
                sources = {cg_opposite(rPts, from, arrived)};
            std::vector<std::size_t> targets = to;
            if (i + 1 == b) {
                // Back to the first loop, well away from where the cycle left it.
                double len = 0.0;
                const std::vector<double> s = cg_arc(rPts, to, len);
                const std::size_t i0 =
                    static_cast<std::size_t>(std::find(to.begin(), to.end(), start) - to.begin());
                targets.clear();
                for (std::size_t k = 0; k < to.size(); ++k) {
                    double d = std::fabs(s[k] - s[i0]);
                    d = std::min(d, len - d);
                    if (d >= 0.25 * len)
                        targets.push_back(to[k]);
                }
            }
            // Loop vertices are never interior to a path: it would run along
            // the boundary.
            std::vector<std::uint8_t> blocked = forbidden;
            for (const auto& loop : loops)
                for (std::size_t v : loop)
                    blocked[v] = 1;
            const std::vector<std::size_t> path =
                patch_shortest_path(rPts, rTris, sources, targets, blocked);
            if (path.size() < 2)
                return {};
            if (i == 0)
                start = path.front();
            arrived = path.back();
            for (std::size_t k = 0; k + 1 < path.size(); ++k)
                cut.push_back(cg_edge(path[k], path[k + 1]));
            for (std::size_t v : path)
                forbidden[v] = 1;
        }
    }
    std::sort(cut.begin(), cut.end());
    cut.erase(std::unique(cut.begin(), cut.end()), cut.end());
    std::size_t num = 0;
    std::vector<std::int64_t> label = cg_components(rTris, cut, num);
    if (num < 2 || !cg_pieces_are_discs(rTris, label, num, MinTriangles))
        return {};
    rNumPieces = num;
    return label;
}

}  // namespace detail
}  // namespace meshioplusplus
