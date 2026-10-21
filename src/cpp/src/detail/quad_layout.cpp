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
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>

// Project includes
#include "meshioplusplus/log.hpp"

// Project includes (private, not installed)
#include "cross_field.hpp"
#include "cut_graph.hpp"
#include "param_disc.hpp"
#include "quad_layout.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kQlPi = 3.14159265358979323846;
constexpr std::size_t kQlNone = std::numeric_limits<std::size_t>::max();
// About the triangles a cell of the layout needs for its patch to fit.
constexpr std::size_t kQlStripTriangles = 48;

PatchVec ql_sub(const PatchVec& rA, const PatchVec& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}
double ql_dot(const PatchVec& rA, const PatchVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
PatchVec ql_cross(const PatchVec& rA, const PatchVec& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
double ql_norm(const PatchVec& rA) {
    return std::sqrt(ql_dot(rA, rA));
}
PatchVec ql_unit(const PatchVec& rA) {
    const double l = ql_norm(rA);
    return l > 0.0 ? PatchVec{rA[0] / l, rA[1] / l, rA[2] / l} : PatchVec{0.0, 0.0, 0.0};
}
// The angle between two vectors, in degrees.
double ql_angle_deg(const PatchVec& rA, const PatchVec& rB) {
    const double la = ql_norm(rA), lb = ql_norm(rB);
    if (!(la > 0.0) || !(lb > 0.0))
        return 180.0;
    const double c = std::max(-1.0, std::min(1.0, ql_dot(rA, rB) / (la * lb)));
    return std::acos(c) * 180.0 / kQlPi;
}

// The region's half-edges: per triangle and local edge k (from vertex k to
// k + 1), the triangle across it, or -1 on the boundary.
struct QlMesh {
    std::size_t mNumPts = 0;
    std::vector<std::array<std::int64_t, 3>> mNbr;
    std::vector<std::vector<std::size_t>> mVertTris;
    std::vector<std::uint8_t> mOnBoundary;  ///< per vertex
    bool mManifold = true;
};

QlMesh ql_mesh(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris) {
    QlMesh m;
    m.mNumPts = rPts.size();
    for (const PatchTri& t : rTris)
        m.mNumPts = std::max({m.mNumPts, t[0] + 1, t[1] + 1, t[2] + 1});
    const std::size_t nt = rTris.size();
    m.mNbr.assign(nt, {-1, -1, -1});
    m.mVertTris.assign(m.mNumPts, {});
    m.mOnBoundary.assign(m.mNumPts, 0);
    std::unordered_map<std::uint64_t, std::size_t> half;
    half.reserve(3 * nt);
    const auto key = [&](std::size_t A, std::size_t B) {
        return static_cast<std::uint64_t>(A) * static_cast<std::uint64_t>(m.mNumPts) +
               static_cast<std::uint64_t>(B);
    };
    for (std::size_t f = 0; f < nt; ++f)
        for (std::size_t k = 0; k < 3; ++k) {
            m.mVertTris[rTris[f][k]].push_back(f);
            if (!half.emplace(key(rTris[f][k], rTris[f][(k + 1) % 3]), f * 3 + k).second)
                m.mManifold = false;
        }
    for (std::size_t f = 0; f < nt; ++f)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::size_t a = rTris[f][k], b = rTris[f][(k + 1) % 3];
            const auto it = half.find(key(b, a));
            if (it == half.end()) {
                m.mOnBoundary[a] = m.mOnBoundary[b] = 1;
            } else {
                m.mNbr[f][k] = static_cast<std::int64_t>(it->second / 3);
            }
        }
    return m;
}

std::vector<PatchTri> ql_local(const std::vector<PatchTri>& rTris,
                               const std::vector<std::size_t>& rIdx) {
    std::vector<PatchTri> out;
    out.reserve(rIdx.size());
    for (std::size_t f : rIdx)
        out.push_back(rTris[f]);
    return out;
}

// The edge-connected components of `rIdx` (ascending) that `rCut(f, k)` does
// not separate, each ascending, numbered by their lowest triangle.
template <class Cut>
std::vector<std::vector<std::size_t>> ql_components(const QlMesh& rMesh,
                                                    const std::vector<std::size_t>& rIdx,
                                                    Cut rCut) {
    std::unordered_map<std::size_t, std::size_t> slot;
    slot.reserve(rIdx.size());
    for (std::size_t i = 0; i < rIdx.size(); ++i)
        slot[rIdx[i]] = i;
    std::vector<std::size_t> comp(rIdx.size(), kQlNone);
    std::vector<std::vector<std::size_t>> out;
    for (std::size_t s = 0; s < rIdx.size(); ++s) {
        if (comp[s] != kQlNone)
            continue;
        const std::size_t c = out.size();
        out.emplace_back();
        std::vector<std::size_t> stack{s};
        comp[s] = c;
        while (!stack.empty()) {
            const std::size_t i = stack.back();
            stack.pop_back();
            const std::size_t f = rIdx[i];
            out[c].push_back(f);
            for (std::size_t k = 0; k < 3; ++k) {
                const std::int64_t g = rMesh.mNbr[f][k];
                if (g < 0 || rCut(f, k))
                    continue;
                const auto it = slot.find(static_cast<std::size_t>(g));
                if (it == slot.end() || comp[it->second] != kQlNone)
                    continue;
                comp[it->second] = c;
                stack.push_back(it->second);
            }
        }
        std::sort(out[c].begin(), out[c].end());
    }
    return out;
}

// Group `rIdx` by per-triangle labels from a split (`-1` triangles are
// grouped by edge connectivity too), pieces by their lowest triangle.
std::vector<std::vector<std::size_t>> ql_groups(const QlMesh& rMesh,
                                                const std::vector<std::size_t>& rIdx,
                                                const std::vector<std::int64_t>& rLabels) {
    std::unordered_map<std::size_t, std::int64_t> lab;
    lab.reserve(rIdx.size());
    for (std::size_t i = 0; i < rIdx.size(); ++i)
        lab[rIdx[i]] = rLabels[i];
    return ql_components(rMesh, rIdx, [&](std::size_t F, std::size_t K) {
        const auto it = lab.find(static_cast<std::size_t>(rMesh.mNbr[F][K]));
        return it == lab.end() || it->second != lab.at(F);
    });
}

// A disc piece's boundary, a vertex cycle in the triangles' orientation,
// starting at its smallest vertex; empty when the boundary is not one cycle.
std::vector<std::size_t> ql_boundary(const QlMesh& rMesh, const std::vector<PatchTri>& rTris,
                                     const std::vector<std::size_t>& rIdx,
                                     const std::vector<std::uint8_t>& rIn) {
    std::unordered_map<std::size_t, std::size_t> next;
    for (std::size_t f : rIdx)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t g = rMesh.mNbr[f][k];
            if (g >= 0 && rIn[static_cast<std::size_t>(g)])
                continue;
            if (!next.emplace(rTris[f][k], rTris[f][(k + 1) % 3]).second)
                return {};
        }
    if (next.empty())
        return {};
    std::size_t start = kQlNone;
    for (const auto& kv : next)
        start = std::min(start, kv.first);
    std::vector<std::size_t> cycle{start};
    for (std::size_t v = next.at(start); v != start; v = next.at(v)) {
        if (cycle.size() > next.size())
            return {};
        cycle.push_back(v);
        if (next.find(v) == next.end())
            return {};
    }
    return cycle.size() == next.size() ? cycle : std::vector<std::size_t>{};
}

struct QlPiece {
    std::vector<std::size_t> mTris;
    std::vector<std::size_t> mForced;
    int mDepth = 0;
    bool mNoSplit = false;
};

// The layout's working state.
struct QlState {
    const std::vector<PatchVec>& mPts;
    const std::vector<PatchTri>& mTris;
    const QuadLayoutOptions& mOptions;
    QlMesh mMesh;
    CrossField mField;
    std::vector<std::uint8_t> mIn;      ///< scratch: per triangle, in the current piece
    std::vector<std::uint8_t> mCorner;  ///< per vertex: a corner some piece chose
};

void ql_mark(QlState& rS, const std::vector<std::size_t>& rIdx, std::uint8_t Value) {
    for (std::size_t f : rIdx)
        rS.mIn[f] = Value;
}

// Per vertex of the piece, its ring distance from the piece's boundary
// (`kQlNone` for vertices outside it).
std::vector<std::size_t> ql_rings(const QlState& rS, const std::vector<std::size_t>& rIdx,
                                  const std::vector<std::size_t>& rBoundary) {
    std::vector<std::size_t> ring(rS.mMesh.mNumPts, kQlNone);
    std::deque<std::size_t> queue;
    for (std::size_t v : rBoundary) {
        ring[v] = 0;
        queue.push_back(v);
    }
    // Adjacency through the piece's triangles only.
    std::unordered_map<std::size_t, std::vector<std::size_t>> adj;
    for (std::size_t f : rIdx)
        for (std::size_t k = 0; k < 3; ++k) {
            adj[rS.mTris[f][k]].push_back(rS.mTris[f][(k + 1) % 3]);
            adj[rS.mTris[f][(k + 1) % 3]].push_back(rS.mTris[f][k]);
        }
    for (auto& kv : adj)
        std::sort(kv.second.begin(), kv.second.end());
    while (!queue.empty()) {
        const std::size_t v = queue.front();
        queue.pop_front();
        for (std::size_t u : adj[v])
            if (ring[u] == kQlNone) {
                ring[u] = ring[v] + 1;
                queue.push_back(u);
            }
    }
    return ring;
}

// The fan of an interior vertex: its triangles in rotation order and, for
// each, the neighbour at the start of its sector.
bool ql_fan(const QlState& rS, std::size_t V, std::vector<std::size_t>& rFan,
            std::vector<std::size_t>& rNbrs) {
    const auto& ring = rS.mMesh.mVertTris[V];
    if (ring.size() < 3)
        return false;
    std::size_t f = *std::min_element(ring.begin(), ring.end());
    rFan.clear();
    rNbrs.clear();
    for (std::size_t step = 0; step < ring.size(); ++step) {
        const PatchTri& t = rS.mTris[f];
        const std::size_t iv = t[0] == V ? 0 : (t[1] == V ? 1 : 2);
        rFan.push_back(f);
        rNbrs.push_back(t[(iv + 1) % 3]);
        // Across the edge (prev, V): local edge (iv + 2) % 3 runs prev -> V.
        const std::int64_t g = rS.mMesh.mNbr[f][(iv + 2) % 3];
        if (g < 0)
            return false;
        f = static_cast<std::size_t>(g);
        if (f == rFan.front())
            return step + 1 == ring.size();
    }
    return false;
}

double ql_tri_angle(const QlState& rS, std::size_t F, std::size_t V) {
    const PatchTri& t = rS.mTris[F];
    const std::size_t iv = t[0] == V ? 0 : (t[1] == V ? 1 : 2);
    return ql_angle_deg(ql_sub(rS.mPts[t[(iv + 1) % 3]], rS.mPts[V]),
                        ql_sub(rS.mPts[t[(iv + 2) % 3]], rS.mPts[V])) *
           kQlPi / 180.0;
}

// Cut a piece along vertex paths: its components with the path edges as cuts.
std::vector<std::vector<std::size_t>> ql_cut_along(
    const QlState& rS, const std::vector<std::size_t>& rIdx,
    const std::vector<std::vector<std::size_t>>& rPaths) {
    std::unordered_map<std::uint64_t, std::uint8_t> cut;
    const auto key = [&](std::size_t A, std::size_t B) {
        const std::size_t lo = std::min(A, B), hi = std::max(A, B);
        return static_cast<std::uint64_t>(lo) * static_cast<std::uint64_t>(rS.mMesh.mNumPts) +
               static_cast<std::uint64_t>(hi);
    };
    for (const auto& p : rPaths)
        for (std::size_t i = 0; i + 1 < p.size(); ++i)
            cut[key(p[i], p[i + 1])] = 1;
    return ql_components(rS.mMesh, rIdx, [&](std::size_t F, std::size_t K) {
        return cut.count(key(rS.mTris[F][K], rS.mTris[F][(K + 1) % 3])) != 0;
    });
}

bool ql_all_discs(const QlState& rS, const std::vector<std::vector<std::size_t>>& rParts) {
    for (const auto& p : rParts)
        if (!patch_topology(ql_local(rS.mTris, p)).IsDisc())
            return false;
    return true;
}

// Split a piece at its one singularity `Sing` along the separatrices traced
// from it (see the header). Appends the pieces to `rOut`.
bool ql_separatrix_split(QlState& rS, const QlPiece& rPiece, std::size_t Sing, int Valence,
                         std::vector<QlPiece>& rOut) {
    std::vector<std::size_t> fan, nbrs;
    if (!ql_fan(rS, Sing, fan, nbrs))
        return false;
    const std::size_t m = fan.size();
    if (m < static_cast<std::size_t>(Valence))
        return false;
    // Polar angle of each neighbour, and its mismatch with the field.
    std::vector<double> phi(m, 0.0), mis(m, 0.0);
    double total = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        phi[i] = total;
        total += ql_tri_angle(rS, fan[i], Sing);
        const PatchVec d = ql_sub(rS.mPts[nbrs[i]], rS.mPts[Sing]);
        const std::size_t f = fan[i];
        const double a = std::atan2(ql_dot(d, rS.mField.mE2[f]), ql_dot(d, rS.mField.mE1[f])) -
                         rS.mField.mTheta[f];
        const double q = 0.5 * kQlPi;
        double w = std::fmod(a, q);
        if (w > 0.5 * q)
            w -= q;
        else if (w <= -0.5 * q)
            w += q;
        mis[i] = std::abs(w);
    }
    const auto nearest = [&](double Target) {
        std::size_t best = 0;
        double bd = 1e300;
        for (std::size_t i = 0; i < m; ++i) {
            double d = std::abs(std::fmod(phi[i] - Target + 3.0 * total, total));
            d = std::min(d, total - d);
            if (d < bd - 1e-12) {
                bd = d;
                best = i;
            }
        }
        return best;
    };
    std::vector<std::size_t> pick;
    double best_cost = 1e300;
    for (std::size_t i0 = 0; i0 < m; ++i0) {
        std::vector<std::size_t> p;
        double cost = 0.0;
        for (int j = 0; j < Valence; ++j) {
            const std::size_t i = nearest(phi[i0] + total * static_cast<double>(j) / Valence);
            if (std::find(p.begin(), p.end(), i) != p.end()) {
                cost = 1e300;
                break;
            }
            p.push_back(i);
            cost += mis[i];
        }
        if (cost < best_cost - 1e-12) {
            best_cost = cost;
            pick = p;
        }
    }
    if (pick.size() != static_cast<std::size_t>(Valence))
        return false;

    // Trace each separatrix to the piece's boundary.
    ql_mark(rS, rPiece.mTris, 1);
    std::vector<std::uint8_t> used(rS.mMesh.mNumPts, 0);
    used[Sing] = 1;
    for (std::size_t i : pick)
        used[nbrs[i]] = 1;
    std::vector<std::vector<std::size_t>> paths;
    bool ok = true;
    for (std::size_t i : pick) {
        std::vector<std::size_t> path{Sing, nbrs[i]};
        PatchVec dir = ql_unit(ql_sub(rS.mPts[nbrs[i]], rS.mPts[Sing]));
        std::size_t v = nbrs[i];
        while (!rS.mMesh.mOnBoundary[v] && ok) {
            // On the piece's boundary?
            bool on_piece_boundary = false;
            for (std::size_t f : rS.mMesh.mVertTris[v])
                if (!rS.mIn[f])
                    on_piece_boundary = true;
            if (on_piece_boundary)
                break;
            if (path.size() > rPiece.mTris.size()) {
                ok = false;
                break;
            }
            // The field's branch closest to the walk's direction, in the
            // triangle round `v` that direction points into.
            std::size_t ft = kQlNone;
            double fbest = -2.0;
            for (std::size_t f : rS.mMesh.mVertTris[v]) {
                const PatchTri& t = rS.mTris[f];
                const PatchVec c{(rS.mPts[t[0]][0] + rS.mPts[t[1]][0] + rS.mPts[t[2]][0]) / 3.0,
                                 (rS.mPts[t[0]][1] + rS.mPts[t[1]][1] + rS.mPts[t[2]][1]) / 3.0,
                                 (rS.mPts[t[0]][2] + rS.mPts[t[1]][2] + rS.mPts[t[2]][2]) / 3.0};
                const double s = ql_dot(ql_unit(ql_sub(c, rS.mPts[v])), dir);
                if (s > fbest) {
                    fbest = s;
                    ft = f;
                }
            }
            PatchVec branch = dir;
            double bbest = -2.0;
            for (int k = 0; k < 4; ++k) {
                const PatchVec b = cross_branch(rS.mField, ft, k);
                const double s = ql_dot(b, dir);
                if (s > bbest) {
                    bbest = s;
                    branch = b;
                }
            }
            // The neighbour best aligned with it.
            std::size_t next = kQlNone;
            double nbest = 0.3;
            for (std::size_t f : rS.mMesh.mVertTris[v])
                for (std::size_t u : rS.mTris[f]) {
                    if (u == v || used[u])
                        continue;
                    const double s = ql_dot(ql_unit(ql_sub(rS.mPts[u], rS.mPts[v])), branch);
                    if (s > nbest + 1e-12 || (std::abs(s - nbest) <= 1e-12 && u < next)) {
                        nbest = s;
                        next = u;
                    }
                }
            if (next == kQlNone) {
                ok = false;
                break;
            }
            dir = ql_unit(ql_sub(rS.mPts[next], rS.mPts[v]));
            used[next] = 1;
            path.push_back(next);
            v = next;
        }
        if (!ok)
            break;
        paths.push_back(std::move(path));
    }
    ql_mark(rS, rPiece.mTris, 0);
    if (!ok)
        return false;
    const auto parts = ql_cut_along(rS, rPiece.mTris, paths);
    if (parts.size() != static_cast<std::size_t>(Valence) || !ql_all_discs(rS, parts))
        return false;
    for (const auto& p : parts)
        if (p.size() < 2)
            return false;
    // Each part lies between two consecutive separatrices: find it by the
    // fan triangle just after separatrix j.
    std::vector<std::size_t> sorted_pick = pick;
    std::sort(sorted_pick.begin(), sorted_pick.end());
    std::vector<QlPiece> made;
    for (std::size_t j = 0; j < sorted_pick.size(); ++j) {
        const std::size_t f = fan[sorted_pick[j]];
        std::size_t part = kQlNone;
        for (std::size_t q = 0; q < parts.size(); ++q)
            if (std::binary_search(parts[q].begin(), parts[q].end(), f))
                part = q;
        if (part == kQlNone)
            return false;
        // The separatrices bounding it: the one through nbrs[sorted_pick[j]]
        // and the next one round.
        const std::size_t ja = sorted_pick[j], jb = sorted_pick[(j + 1) % sorted_pick.size()];
        std::size_t ea = kQlNone, eb = kQlNone;
        for (std::size_t q = 0; q < pick.size(); ++q) {
            if (pick[q] == ja)
                ea = paths[q].back();
            if (pick[q] == jb)
                eb = paths[q].back();
        }
        QlPiece np;
        np.mTris = parts[part];
        np.mForced = {Sing, ea, eb};
        np.mDepth = rPiece.mDepth + 1;
        np.mNoSplit = true;
        made.push_back(std::move(np));
    }
    std::sort(made.begin(), made.end(),
              [](const QlPiece& rA, const QlPiece& rB) { return rA.mTris[0] < rB.mTris[0]; });
    for (std::size_t q = 1; q < made.size(); ++q)
        if (made[q].mTris[0] == made[q - 1].mTris[0])
            return false;  // two separatrix sectors in one part
    for (auto& p : made)
        rOut.push_back(std::move(p));
    return true;
}

// The chord angle (degrees, 360 - it when concave) at each boundary position,
// over `K` vertices either side.
std::vector<double> ql_turns(const QlState& rS, const std::vector<std::size_t>& rB,
                             const std::vector<std::uint8_t>& rIn, std::size_t K) {
    const std::size_t nb = rB.size();
    std::vector<double> beta(nb, 180.0);
    for (std::size_t i = 0; i < nb; ++i) {
        const std::size_t v = rB[i];
        const PatchVec& p = rS.mPts[v];
        const PatchVec& a = rS.mPts[rB[(i + nb - K) % nb]];
        const PatchVec& b = rS.mPts[rB[(i + K) % nb]];
        PatchVec n{0.0, 0.0, 0.0};
        for (std::size_t f : rS.mMesh.mVertTris[v]) {
            if (!rIn[f])
                continue;
            const PatchTri& t = rS.mTris[f];
            const PatchVec c = ql_cross(ql_sub(rS.mPts[t[1]], rS.mPts[t[0]]),
                                        ql_sub(rS.mPts[t[2]], rS.mPts[t[0]]));
            for (std::size_t j = 0; j < 3; ++j)
                n[j] += c[j];
        }
        // In the tangent plane: the chords without their normal parts.
        n = ql_unit(n);
        const auto flat = [&](const PatchVec& rD) {
            const double h = ql_dot(rD, n);
            return PatchVec{rD[0] - h * n[0], rD[1] - h * n[1], rD[2] - h * n[2]};
        };
        const PatchVec da = flat(ql_sub(a, p)), db = flat(ql_sub(b, p));
        double ang = ql_angle_deg(da, db);
        if (ql_dot(ql_cross(PatchVec{-da[0], -da[1], -da[2]}, db), n) < 0.0)
            ang = 360.0 - ang;
        beta[i] = ang;
    }
    return beta;
}

// Split a disc into an n-gon's quads: paths from its innermost vertex to the
// middles of the sides between the corners `rC` (positions in `rB`, sorted).
bool ql_split_ngon(QlState& rS, const QlPiece& rPiece, const std::vector<std::size_t>& rB,
                   const std::vector<std::size_t>& rC, std::vector<QlPiece>& rMore) {
    const std::size_t nb = rB.size();
    std::vector<double> s(nb + 1, 0.0);
    for (std::size_t i = 0; i < nb; ++i)
        s[i + 1] = s[i] + ql_norm(ql_sub(rS.mPts[rB[(i + 1) % nb]], rS.mPts[rB[i]]));
    const double perim = s[nb];
    const auto gap = [&](std::size_t A, std::size_t B) {
        return B >= A ? s[B] - s[A] : perim - s[A] + s[B];
    };
    // An n-gon: quads from its innermost vertex to its sides' middles.
    const std::vector<std::size_t>& c = rC;
    const std::size_t n = c.size();
    const std::vector<std::size_t> ring = ql_rings(rS, rPiece.mTris, rB);
    std::size_t centre = kQlNone, deepest = 0;
    for (std::size_t f : rPiece.mTris)
        for (std::size_t v : rS.mTris[f])
            if (ring[v] != kQlNone &&
                (ring[v] > deepest || (ring[v] == deepest && v < centre && ring[v] > 0))) {
                deepest = ring[v];
                centre = v;
            }
    std::vector<std::size_t> mids;
    bool ok = deepest >= 2;
    for (std::size_t j = 0; j < n && ok; ++j) {
        const std::size_t a = c[j], b = c[(j + 1) % n];
        const double half = 0.5 * gap(a, b);
        std::size_t best = kQlNone;
        double bd = 1e300;
        for (std::size_t i = (a + 1) % nb; i != b; i = (i + 1) % nb) {
            const double d = std::abs(gap(a, i) - half);
            if (d < bd) {
                bd = d;
                best = i;
            }
        }
        if (best == kQlNone)
            ok = false;
        else
            mids.push_back(best);
    }
    if (ok) {
        const std::vector<PatchTri> local = ql_local(rS.mTris, rPiece.mTris);
        std::vector<std::uint8_t> forbid(rS.mMesh.mNumPts, 0);
        for (std::size_t v : rB)
            forbid[v] = 1;
        std::vector<std::vector<std::size_t>> paths;
        for (std::size_t j = 0; j < n && ok; ++j) {
            auto p = patch_shortest_path(rS.mPts, local, {centre}, {rB[mids[j]]}, forbid);
            if (p.size() < 2) {
                ok = false;
                break;
            }
            for (std::size_t v : p)
                if (v != centre)
                    forbid[v] = 1;
            paths.push_back(std::move(p));
        }
        if (ok) {
            const auto parts = ql_cut_along(rS, rPiece.mTris, paths);
            ok = parts.size() == n && ql_all_discs(rS, parts);
            std::vector<QlPiece> made;
            for (std::size_t j = 0; j < n && ok; ++j) {
                // The part holding corner j: its triangles there.
                const std::size_t cv = rB[c[j]];
                std::size_t part = kQlNone;
                for (std::size_t q = 0; q < parts.size() && part == kQlNone; ++q)
                    for (std::size_t f : rS.mMesh.mVertTris[cv])
                        if (std::binary_search(parts[q].begin(), parts[q].end(), f))
                            part = q;
                if (part == kQlNone) {
                    ok = false;
                    break;
                }
                QlPiece np;
                np.mTris = parts[part];
                np.mForced = {centre, rB[mids[(j + n - 1) % n]], cv, rB[mids[j]]};
                np.mDepth = rPiece.mDepth + 1;
                np.mNoSplit = true;
                made.push_back(std::move(np));
            }
            if (ok) {
                std::vector<std::size_t> firsts;
                for (const auto& p : made)
                    firsts.push_back(p.mTris[0]);
                std::sort(firsts.begin(), firsts.end());
                ok = std::adjacent_find(firsts.begin(), firsts.end()) == firsts.end();
            }
            if (ok) {
                for (auto& p : made) {
                    for (std::size_t v : p.mForced)
                        rS.mCorner[v] = 1;
                    rMore.push_back(std::move(p));
                }
                return true;
            }
        }
    }
    return false;
}

// Split a disc whose four corners `rC` (positions in `rB`, sorted) are not
// all real into an O-grid: a centre quad on the square map's inner square
// and four quads round it, so each side of the disc is one quad's and a
// corner on a smooth boundary is two quads' right angles instead of one
// quad's straight one.
bool ql_split_ogrid(QlState& rS, const QlPiece& rPiece, const std::vector<std::size_t>& rB,
                    const std::vector<std::size_t>& rC, std::vector<QlPiece>& rMore) {
    const std::vector<PatchTri> local = ql_local(rS.mTris, rPiece.mTris);
    const std::array<std::size_t, 4> corners{rC[0], rC[1], rC[2], rC[3]};
    const DiscParam par = param_square(rS.mPts, local, rB, corners);
    if (!par.mOk)
        return false;
    const std::vector<std::size_t> ring = ql_rings(rS, rPiece.mTris, rB);
    // The inner square's corners: the deep vertices nearest (1/4, 1/4) and
    // the like.
    static const double kAt[4][2] = {{0.25, 0.25}, {0.75, 0.25}, {0.75, 0.75}, {0.25, 0.75}};
    std::array<std::size_t, 4> inner{kQlNone, kQlNone, kQlNone, kQlNone};
    for (std::size_t k = 0; k < 4; ++k) {
        double bd = 1e300;
        for (std::size_t f : rPiece.mTris)
            for (std::size_t v : rS.mTris[f]) {
                if (ring[v] == kQlNone || ring[v] < 2)
                    continue;
                const double du = par.mU[v] - kAt[k][0], dv = par.mV[v] - kAt[k][1];
                const double d = du * du + dv * dv;
                if (d < bd - 1e-15 || (std::abs(d - bd) <= 1e-15 && v < inner[k])) {
                    bd = d;
                    inner[k] = v;
                }
            }
        if (inner[k] == kQlNone)
            return false;
        for (std::size_t j = 0; j < k; ++j)
            if (inner[j] == inner[k])
                return false;
    }
    std::vector<std::uint8_t> forbid(rS.mMesh.mNumPts, 0);
    for (std::size_t v : rB)
        forbid[v] = 1;
    std::vector<std::vector<std::size_t>> paths;
    for (std::size_t k = 0; k < 4; ++k) {
        auto p = patch_shortest_path(rS.mPts, local, {inner[k]}, {inner[(k + 1) % 4]}, forbid);
        if (p.size() < 2)
            return false;
        for (std::size_t v : p)
            forbid[v] = 1;
        paths.push_back(std::move(p));
    }
    for (std::size_t k = 0; k < 4; ++k) {
        auto p = patch_shortest_path(rS.mPts, local, {rB[rC[k]]}, {inner[k]}, forbid);
        if (p.size() < 2)
            return false;
        for (std::size_t v : p)
            forbid[v] = 1;
        paths.push_back(std::move(p));
    }
    const auto parts = ql_cut_along(rS, rPiece.mTris, paths);
    if (parts.size() != 5 || !ql_all_discs(rS, parts))
        return false;
    const std::size_t nb = rB.size();
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> owner;
    for (std::size_t f : rPiece.mTris)
        for (std::size_t k = 0; k < 3; ++k)
            owner[{rS.mTris[f][k], rS.mTris[f][(k + 1) % 3]}] = f;
    std::vector<QlPiece> made;
    std::vector<std::uint8_t> used(parts.size(), 0);
    for (std::size_t k = 0; k < 4; ++k) {
        // The ring quad on side k: the triangle on its first boundary edge.
        const auto it = owner.find({rB[rC[k]], rB[(rC[k] + 1) % nb]});
        if (it == owner.end())
            return false;
        std::size_t part = kQlNone;
        for (std::size_t q = 0; q < parts.size(); ++q)
            if (std::binary_search(parts[q].begin(), parts[q].end(), it->second))
                part = q;
        if (part == kQlNone || used[part])
            return false;
        used[part] = 1;
        QlPiece np;
        np.mTris = parts[part];
        np.mForced = {rB[rC[k]], rB[rC[(k + 1) % 4]], inner[(k + 1) % 4], inner[k]};
        np.mDepth = rPiece.mDepth + 1;
        np.mNoSplit = true;
        made.push_back(std::move(np));
    }
    for (std::size_t q = 0; q < parts.size(); ++q)
        if (!used[q]) {
            QlPiece np;
            np.mTris = parts[q];
            np.mForced = {inner[0], inner[1], inner[2], inner[3]};
            np.mDepth = rPiece.mDepth + 1;
            np.mNoSplit = true;
            made.push_back(std::move(np));
        }
    for (auto& p : made) {
        for (std::size_t v : p.mForced)
            rS.mCorner[v] = 1;
        rMore.push_back(std::move(p));
    }
    return true;
}

// Choose four corners for a disc piece, or split it into an n-gon's quads
// (appending them to `rMore`). Returns the corner positions in `rB` (sorted),
// empty on failure; `rSplit` set when the piece was split instead.
std::vector<std::size_t> ql_corners(QlState& rS, const QlPiece& rPiece,
                                    const std::vector<std::size_t>& rB,
                                    const std::vector<std::uint8_t>& rJunction,
                                    std::vector<QlPiece>& rMore, bool& rSplit) {
    rSplit = false;
    const std::size_t nb = rB.size();
    if (nb < 4)
        return {};
    std::unordered_map<std::size_t, std::size_t> pos;
    for (std::size_t i = 0; i < nb; ++i)
        pos[rB[i]] = i;
    std::vector<double> s(nb + 1, 0.0);
    for (std::size_t i = 0; i < nb; ++i)
        s[i + 1] = s[i] + ql_norm(ql_sub(rS.mPts[rB[(i + 1) % nb]], rS.mPts[rB[i]]));
    const double perim = s[nb];
    const auto gap = [&](std::size_t A, std::size_t B) {  // arc length from A forward to B
        return B >= A ? s[B] - s[A] : perim - s[A] + s[B];
    };
    const std::size_t k = std::max<std::size_t>(1, std::min<std::size_t>(4, nb / 8));
    ql_mark(rS, rPiece.mTris, 1);
    const std::vector<double> beta = ql_turns(rS, rB, rS.mIn, k);
    // A corner's turn is concentrated: over twice the window it barely grows,
    // where a smoothly curving side's doubles.
    const std::size_t k2 = std::min(2 * k, std::max<std::size_t>(k, nb / 3));
    const std::vector<double> beta2 = k2 > k ? ql_turns(rS, rB, rS.mIn, k2) : beta;
    ql_mark(rS, rPiece.mTris, 0);
    const auto sharp = [&](std::size_t I) {
        const double t = 180.0 - beta[I], t2 = 180.0 - beta2[I];
        return t > 30.0 && t >= 0.7 * t2;
    };

    std::vector<std::size_t> chosen;
    const auto has = [&](std::size_t I) {
        return std::find(chosen.begin(), chosen.end(), I) != chosen.end();
    };
    for (std::size_t v : rPiece.mForced) {
        const auto it = pos.find(v);
        if (it != pos.end() && !has(it->second) && chosen.size() < 4)
            chosen.push_back(it->second);
    }
    // Sharp turns: locally sharpest under 150 degrees.
    std::vector<std::size_t> strong;
    for (std::size_t i = 0; i < nb; ++i) {
        if (!sharp(i))
            continue;
        bool best = true;
        for (std::size_t d = 1; d <= k && best; ++d) {
            const std::size_t l = (i + nb - d) % nb, r = (i + d) % nb;
            if (beta[l] < beta[i] || beta[r] < beta[i] || (beta[l] == beta[i] && l < i) ||
                (beta[r] == beta[i] && r < i))
                best = false;
        }
        if (best)
            strong.push_back(i);
    }
    std::stable_sort(strong.begin(), strong.end(),
                     [&](std::size_t A, std::size_t B) { return beta[A] < beta[B]; });
    const int dev = rS.mOptions.mMaxValenceDeviation;
    if (rPiece.mForced.empty() && !rPiece.mNoSplit && strong.size() >= 5 &&
        static_cast<int>(strong.size()) <= 4 + dev) {
        std::vector<std::size_t> c = strong;
        std::sort(c.begin(), c.end());
        if (ql_split_ngon(rS, rPiece, rB, c, rMore)) {
            rSplit = true;
            return {};
        }
    }
    const double min_gap = 0.05 * perim;
    const auto far_enough = [&](std::size_t I) {
        for (std::size_t c : chosen)
            if (std::min(gap(c, I), gap(I, c)) < min_gap || c == I)
                return false;
        return true;
    };
    for (std::size_t i : strong)
        if (chosen.size() < 4 && !has(i) && far_enough(i))
            chosen.push_back(i);
    // Junctions and the neighbours' corners, spread out.
    std::vector<std::size_t> inherited;
    for (std::size_t i = 0; i < nb; ++i)
        if ((rJunction[rB[i]] || rS.mCorner[rB[i]]) && !has(i))
            inherited.push_back(i);
    while (chosen.size() < 4 && !inherited.empty()) {
        std::size_t best = kQlNone;
        double bd = -1.0;
        for (std::size_t i : inherited) {
            if (has(i) || !far_enough(i))
                continue;
            double d = perim;
            for (std::size_t c : chosen)
                d = std::min({d, gap(c, i), gap(i, c)});
            // Junctions first, then spread.
            const double score = d + (rJunction[rB[i]] ? perim : 0.0);
            if (score > bd + 1e-12) {
                bd = score;
                best = i;
            }
        }
        if (best == kQlNone)
            break;
        chosen.push_back(best);
    }
    // Fill the longest gaps.
    if (chosen.empty()) {
        std::size_t sharp = 0;
        for (std::size_t i = 1; i < nb; ++i)
            if (beta[i] < beta[sharp] - 1e-9)
                sharp = i;
        chosen.push_back(sharp);
    }
    while (chosen.size() < 4) {
        std::vector<std::size_t> c = chosen;
        std::sort(c.begin(), c.end());
        std::size_t ga = kQlNone, gb = kQlNone;
        double glen = -1.0;
        for (std::size_t j = 0; j < c.size(); ++j) {
            const std::size_t a = c[j], b = c[(j + 1) % c.size()];
            const double len = c.size() == 1 ? perim : gap(a, b);
            if (len > glen + 1e-12) {
                glen = len;
                ga = a;
                gb = b;
            }
        }
        const double half = 0.5 * glen;
        std::size_t best = kQlNone;
        double bd = 1e300;
        for (std::size_t t = 1; t < nb; ++t) {
            const std::size_t i = (ga + t) % nb;
            if (i == gb && c.size() > 1)
                break;
            const double d = std::abs(gap(ga, i) - half);
            if (d < bd) {
                bd = d;
                best = i;
            }
        }
        if (best == kQlNone || has(best))
            return {};
        chosen.push_back(best);
    }
    std::sort(chosen.begin(), chosen.end());
    if (std::adjacent_find(chosen.begin(), chosen.end()) != chosen.end())
        return {};
    // Corners that are no turn at all make degenerate quad corners: three
    // real ones make a triangle's three quads, fewer an O-grid.
    if (rPiece.mForced.empty() && !rPiece.mNoSplit) {
        std::vector<std::size_t> real;
        for (std::size_t i : chosen)
            if (sharp(i))
                real.push_back(i);
        bool done = false;
        if (real.size() == 3 && dev >= 1)
            done = ql_split_ngon(rS, rPiece, rB, real, rMore);
        else if (real.size() <= 2 && rPiece.mTris.size() >= 5 * kQlStripTriangles)
            done = ql_split_ogrid(rS, rPiece, rB, chosen, rMore);
        if (done) {
            rSplit = true;
            return {};
        }
    }
    for (std::size_t i : chosen)
        rS.mCorner[rB[i]] = 1;
    return chosen;
}

struct QlUnionFind {
    std::vector<std::size_t> mParent;
    explicit QlUnionFind(std::size_t N) : mParent(N) {
        std::iota(mParent.begin(), mParent.end(), std::size_t{0});
    }
    std::size_t Find(std::size_t A) {
        while (mParent[A] != A)
            A = mParent[A] = mParent[mParent[A]];
        return A;
    }
    void Unite(std::size_t A, std::size_t B) {
        A = Find(A);
        B = Find(B);
        if (A != B)
            mParent[std::max(A, B)] = std::min(A, B);
    }
};

// A finished piece: its triangles, its boundary from its first corner, its
// corners' positions there and, per direction (sides 0 and 2, sides 1 and
// 3), its strip count.
struct QlQuad {
    std::vector<std::size_t> mTris, mBoundary;
    std::array<std::size_t, 4> mCorners{};
    std::array<std::size_t, 2> mStrips{1, 1};
};

constexpr std::size_t kQlMaxStrips = 16;

// The vertex of a boundary stretch nearest arc length `At` from its start:
// `rStretch` positions into `rB`, ties to the smaller vertex.
std::size_t ql_at_length(const QlState& rS, const std::vector<std::size_t>& rB,
                         const std::vector<std::size_t>& rStretch, const std::vector<double>& rCum,
                         double At) {
    std::size_t best = kQlNone;
    double bd = 1e300;
    for (std::size_t t = 1; t + 1 < rStretch.size(); ++t) {
        const double d = std::abs(rCum[t] - At);
        if (d < bd - 1e-12 * (1.0 + rCum.back()) ||
            (std::abs(d - bd) <= 1e-12 * (1.0 + rCum.back()) &&
             rB[rStretch[t]] < rB[rStretch[best]])) {
            bd = d;
            best = t;
        }
    }
    (void)rS;
    return best;
}

// Split quad `rQ` into `N` strips across direction `D` (its sides D and D + 2
// cut at equal arc-length fractions, joined by shortest paths); false, and
// nothing appended, when it does not split cleanly.
bool ql_split_strips(QlState& rS, const QlQuad& rQ, std::size_t D, std::size_t N,
                     std::vector<QlQuad>& rOut) {
    const std::vector<std::size_t>& b = rQ.mBoundary;
    const std::size_t nb = b.size();
    const auto fail = [&](const char* pWhy) {
        log::debug("quad_layout: {} strips across direction {} of a piece of {} triangles: {}", N,
                   D, rQ.mTris.size(), pWhy);
        return false;
    };
    const auto stretch = [&](std::size_t From, std::size_t To) {  // positions From..To forward
        std::vector<std::size_t> out{From};
        for (std::size_t t = From; t != To;) {
            t = (t + 1) % nb;
            out.push_back(t);
        }
        return out;
    };
    const auto cumulative = [&](const std::vector<std::size_t>& rSt) {
        std::vector<double> c(rSt.size(), 0.0);
        for (std::size_t t = 1; t < rSt.size(); ++t)
            c[t] = c[t - 1] + ql_norm(ql_sub(rS.mPts[b[rSt[t]]], rS.mPts[b[rSt[t - 1]]]));
        return c;
    };
    const std::size_t c0 = rQ.mCorners[D], c1 = rQ.mCorners[(D + 1) % 4];
    const std::size_t c2 = rQ.mCorners[(D + 2) % 4], c3 = rQ.mCorners[(D + 3) % 4];
    std::vector<std::size_t> sa = stretch(c0, c1), sc = stretch(c2, c3);
    std::reverse(sc.begin(), sc.end());  // from corner D + 3 back to D + 2
    if (sa.size() < N + 1 || sc.size() < N + 1)
        return fail("a side has fewer edges than strips");
    const std::vector<double> ca = cumulative(sa), cc = cumulative(sc);
    std::vector<std::size_t> pa{0}, pc{0};
    for (std::size_t k = 1; k < N; ++k) {
        const std::size_t ta = ql_at_length(rS, b, sa, ca, ca.back() * static_cast<double>(k) / N);
        const std::size_t tc = ql_at_length(rS, b, sc, cc, cc.back() * static_cast<double>(k) / N);
        if (ta == kQlNone || tc == kQlNone || ta <= pa.back() || tc <= pc.back())
            return fail("two cut points coincide");
        pa.push_back(ta);
        pc.push_back(tc);
    }
    pa.push_back(sa.size() - 1);
    pc.push_back(sc.size() - 1);
    // Cut paths between the matching points, through the interior.
    const std::vector<PatchTri> local = ql_local(rS.mTris, rQ.mTris);
    std::vector<std::uint8_t> forbid(rS.mMesh.mNumPts, 0);
    for (std::size_t v : b)
        forbid[v] = 1;
    std::vector<std::vector<std::size_t>> paths;
    for (std::size_t k = 1; k < N; ++k) {
        const std::size_t va = b[sa[pa[k]]], vc = b[sc[pc[k]]];
        auto path = patch_shortest_path(rS.mPts, local, {va}, {vc}, forbid);
        if (path.size() < 2)
            return fail("no cut path joins two cut points");
        for (std::size_t v : path)
            forbid[v] = 1;
        paths.push_back(std::move(path));
    }
    const auto parts = ql_cut_along(rS, rQ.mTris, paths);
    if (parts.size() != N || !ql_all_discs(rS, parts))
        return fail("the cuts do not give discs");
    // The triangle on each half-edge of the boundary, to find the strips.
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> owner;
    for (std::size_t f : rQ.mTris)
        for (std::size_t k = 0; k < 3; ++k)
            owner[{rS.mTris[f][k], rS.mTris[f][(k + 1) % 3]}] = f;
    std::vector<QlQuad> made;
    for (std::size_t k = 0; k < N; ++k) {
        const std::size_t va = b[sa[pa[k]]];
        const auto it = owner.find({va, b[(sa[pa[k]] + 1) % nb]});
        if (it == owner.end())
            return fail("a strip has no first edge");
        std::size_t part = kQlNone;
        for (std::size_t q = 0; q < parts.size(); ++q)
            if (std::binary_search(parts[q].begin(), parts[q].end(), it->second))
                part = q;
        if (part == kQlNone)
            return fail("a strip was not found");
        ql_mark(rS, parts[part], 1);
        const std::vector<std::size_t> sb = ql_boundary(rS.mMesh, rS.mTris, parts[part], rS.mIn);
        ql_mark(rS, parts[part], 0);
        // Its corners in the rotated frame, then in the quad's numbering:
        // corner j is rotated corner (j - D) mod 4.
        const std::array<std::size_t, 4> rot{va, b[sa[pa[k + 1]]], b[sc[pc[k + 1]]], b[sc[pc[k]]]};
        std::array<std::size_t, 4> cv{};
        for (std::size_t j = 0; j < 4; ++j)
            cv[j] = rot[(j + 4 - D) % 4];
        std::array<std::size_t, 4> pos{};
        for (std::size_t j = 0; j < 4; ++j) {
            const auto f = std::find(sb.begin(), sb.end(), cv[j]);
            if (f == sb.end())
                return fail("a corner is off its strip");
            pos[j] = static_cast<std::size_t>(f - sb.begin());
        }
        QlQuad s;
        s.mTris = parts[part];
        s.mBoundary.assign(sb.begin() + static_cast<std::ptrdiff_t>(pos[0]), sb.end());
        s.mBoundary.insert(s.mBoundary.end(), sb.begin(),
                           sb.begin() + static_cast<std::ptrdiff_t>(pos[0]));
        for (std::size_t j = 0; j < 4; ++j)
            s.mCorners[j] = (pos[j] + sb.size() - pos[0]) % sb.size();
        for (std::size_t j = 1; j < 4; ++j)
            if (s.mCorners[j] <= s.mCorners[j - 1])
                return fail("a strip's corners are out of order");
        s.mStrips = rQ.mStrips;
        s.mStrips[D % 2] = 1;
        made.push_back(std::move(s));
    }
    for (auto& m : made)
        rOut.push_back(std::move(m));
    return true;
}

// Refine the quads to the target size without T-junctions between them: a
// sheet is the set of directions (piece, sides 0/2 or 1/3) linked across
// shared sides; each sheet takes the strip count its longest side needs, and
// every piece is cut into that many strips along each direction at the same
// arc-length fractions its neighbours cut the shared sides at. A piece that
// does not split cleanly stays whole (its neighbours' cuts end on its side).
std::vector<QlQuad> ql_refine(QlState& rS, std::vector<QlQuad> Quads, double Target) {
    const std::size_t nq = Quads.size();
    std::map<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t>> half;
    std::vector<std::array<double, 4>> len(nq, {0.0, 0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < nq; ++i) {
        const QlQuad& q = Quads[i];
        const std::size_t nb = q.mBoundary.size();
        for (std::size_t j = 0; j < 4; ++j) {
            const std::size_t a = q.mCorners[j], e = j + 1 < 4 ? q.mCorners[j + 1] : nb;
            for (std::size_t t = a; t < e; ++t) {
                const std::size_t u = q.mBoundary[t], v = q.mBoundary[(t + 1) % nb];
                half[{u, v}] = {i, j};
                len[i][j] += ql_norm(ql_sub(rS.mPts[v], rS.mPts[u]));
            }
        }
    }
    QlUnionFind uf(2 * nq);
    for (const auto& kv : half) {
        const auto twin = half.find({kv.first.second, kv.first.first});
        if (twin != half.end())
            uf.Unite(2 * kv.second.first + kv.second.second % 2,
                     2 * twin->second.first + twin->second.second % 2);
    }
    std::vector<std::size_t> count(2 * nq, 1), allowed(2 * nq, kQlMaxStrips);
    for (std::size_t i = 0; i < nq; ++i) {
        // A strip needs data to fit: about kQlStripTriangles triangles a cell.
        const std::size_t most = std::max<std::size_t>(
            1, static_cast<std::size_t>(std::sqrt(static_cast<double>(Quads[i].mTris.size()) /
                                                  static_cast<double>(kQlStripTriangles))));
        for (std::size_t j = 0; j < 4; ++j) {
            const double want = std::ceil(len[i][j] / Target - 1e-9);
            const std::size_t n = static_cast<std::size_t>(
                std::max(1.0, std::min(static_cast<double>(kQlMaxStrips), want)));
            const std::size_t sheet = uf.Find(2 * i + j % 2);
            count[sheet] = std::max(count[sheet], n);
            allowed[sheet] = std::min(allowed[sheet], most);
        }
    }
    for (std::size_t k = 0; k < 2 * nq; ++k)
        count[k] = std::max<std::size_t>(1, std::min(count[k], allowed[k]));
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t d = 0; d < 2; ++d)
            Quads[i].mStrips[d] = count[uf.Find(2 * i + d)];
    for (std::size_t d = 0; d < 2; ++d) {
        std::vector<QlQuad> next;
        for (QlQuad& q : Quads) {
            const std::size_t n = q.mStrips[d];
            if (n > 1 && ql_split_strips(rS, q, d, n, next))
                continue;
            q.mStrips[d] = 1;
            next.push_back(std::move(q));
        }
        Quads = std::move(next);
    }
    std::sort(Quads.begin(), Quads.end(),
              [](const QlQuad& rA, const QlQuad& rB) { return rA.mTris[0] < rB.mTris[0]; });
    return Quads;
}

}  // namespace

QuadLayout quad_layout(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       const QuadLayoutOptions& rOptions) {
    QuadLayout out;
    const std::size_t nt = rTris.size();
    if (nt == 0) {
        out.mFailure = QuadLayoutFailure::CutFailed;
        return out;
    }
    QlState st{rPts, rTris, rOptions, ql_mesh(rPts, rTris), {}, {}, {}};
    if (!st.mMesh.mManifold) {
        out.mFailure = QuadLayoutFailure::NonManifold;
        return out;
    }
    st.mIn.assign(nt, 0);
    st.mCorner.assign(st.mMesh.mNumPts, 0);

    // 1. Features: creases and label changes are cuts.
    std::vector<PatchVec> normal(nt);
    for (std::size_t f = 0; f < nt; ++f)
        normal[f] = ql_unit(ql_cross(ql_sub(rPts[rTris[f][1]], rPts[rTris[f][0]]),
                                     ql_sub(rPts[rTris[f][2]], rPts[rTris[f][0]])));
    const bool labelled = rOptions.mLabels.size() == nt;
    std::vector<std::array<std::uint8_t, 3>> feature(nt, {0, 0, 0});
    CrossFieldOptions co;
    for (std::size_t f = 0; f < nt; ++f)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t g = st.mMesh.mNbr[f][k];
            if (g < 0)
                continue;
            const std::size_t gg = static_cast<std::size_t>(g);
            bool cut = labelled && rOptions.mLabels[f] != rOptions.mLabels[gg];
            if (rOptions.mSnapToFeatures &&
                ql_angle_deg(normal[f], normal[gg]) > rOptions.mFeatureAngleDeg)
                cut = true;
            if (cut)
                feature[f][k] = 1;
        }
    // A crease must run somewhere: erode feature edges that end inside the
    // region (a noisy skin's isolated sharp edges), keeping closed curves and
    // curves that end on the boundary.
    {
        std::vector<int> degree(st.mMesh.mNumPts, 0);
        for (std::size_t f = 0; f < nt; ++f)
            for (std::size_t k = 0; k < 3; ++k)
                if (feature[f][k] && static_cast<std::size_t>(st.mMesh.mNbr[f][k]) > f) {
                    ++degree[rTris[f][k]];
                    ++degree[rTris[f][(k + 1) % 3]];
                }
        const auto dangling = [&](std::size_t V) {
            return degree[V] == 1 && !st.mMesh.mOnBoundary[V];
        };
        for (bool changed = true; changed;) {
            changed = false;
            for (std::size_t f = 0; f < nt; ++f)
                for (std::size_t k = 0; k < 3; ++k) {
                    if (!feature[f][k])
                        continue;
                    const std::size_t a = rTris[f][k], b = rTris[f][(k + 1) % 3];
                    if (!dangling(a) && !dangling(b))
                        continue;
                    const std::size_t g = static_cast<std::size_t>(st.mMesh.mNbr[f][k]);
                    for (std::size_t kk = 0; kk < 3; ++kk)
                        if (rTris[g][kk] == b && rTris[g][(kk + 1) % 3] == a)
                            feature[g][kk] = 0;
                    feature[f][k] = 0;
                    --degree[a];
                    --degree[b];
                    changed = true;
                }
        }
    }
    for (std::size_t f = 0; f < nt; ++f)
        for (std::size_t k = 0; k < 3; ++k)
            if (feature[f][k] && static_cast<std::size_t>(st.mMesh.mNbr[f][k]) > f)
                co.mFeatureEdges.push_back({rTris[f][k], rTris[f][(k + 1) % 3]});
    std::vector<std::size_t> all(nt);
    std::iota(all.begin(), all.end(), std::size_t{0});
    std::deque<QlPiece> work;
    for (auto& c : ql_components(
             st.mMesh, all, [&](std::size_t F, std::size_t K) { return feature[F][K] != 0; })) {
        QlPiece p;
        p.mTris = std::move(c);
        work.push_back(std::move(p));
    }
    if (rOptions.mUseCrossField)
        st.mField = cross_field(rPts, rTris, co);

    // 2-3. Discs and singularities.
    std::vector<QlPiece> final;
    const auto push_groups = [&](const QlPiece& rP, const std::vector<std::int64_t>& rLab) {
        for (auto& g : ql_groups(st.mMesh, rP.mTris, rLab)) {
            QlPiece np;
            np.mTris = std::move(g);
            np.mDepth = rP.mDepth + 1;
            for (std::size_t v : rP.mForced)
                np.mForced.push_back(v);
            work.push_back(std::move(np));
        }
    };
    while (!work.empty()) {
        QlPiece p = std::move(work.front());
        work.pop_front();
        const std::vector<PatchTri> local = ql_local(rTris, p.mTris);
        const PatchTopology topo = patch_topology(local);
        if (!topo.IsDisc()) {
            // Not a disc, or pinched at a vertex: cut it (a pinch by halving).
            if (p.mDepth >= rOptions.mMaxDepth) {
                out.mFailure = QuadLayoutFailure::CutFailed;
                return out;
            }
            std::size_t np = 0;
            std::vector<std::int64_t> lab;
            if (topo.mManifold && topo.mGenus == 0)
                lab = patch_split_to_discs(rPts, local, rOptions.mMinTriangles, np);
            if (np < 2)
                lab = patch_bisect(rPts, local, 1, np);
            if (np < 2 || lab.size() != local.size()) {
                out.mFailure = QuadLayoutFailure::CutFailed;
                return out;
            }
            push_groups(p, lab);
            continue;
        }
        if (st.mField.mOk && !p.mNoSplit) {
            ql_mark(st, p.mTris, 1);
            const std::vector<std::size_t> b = ql_boundary(st.mMesh, rTris, p.mTris, st.mIn);
            ql_mark(st, p.mTris, 0);
            const std::vector<std::size_t> ring = ql_rings(st, p.mTris, b);
            std::vector<std::size_t> sing;
            for (std::size_t v : st.mField.mSingular)
                if (ring[v] != kQlNone && ring[v] >= 2 && !st.mMesh.mOnBoundary[v])
                    sing.push_back(v);
            if (sing.size() == 1) {
                const int valence = cross_valence(st.mField, sing[0]);
                if (valence >= 3 && valence != 4 &&
                    std::abs(valence - 4) <= rOptions.mMaxValenceDeviation) {
                    std::vector<QlPiece> made;
                    if (ql_separatrix_split(st, p, sing[0], valence, made)) {
                        ++out.mNumSeparatrixSplits;
                        for (auto& m : made)
                            work.push_back(std::move(m));
                        continue;
                    }
                }
            }
        }
        final.push_back(std::move(p));
    }

    // A piece under `mMinTriangles` (a sliver between creases or cuts) joins
    // the neighbour it shares most edges with, where the union is a disc.
    for (bool changed = true; changed;) {
        changed = false;
        std::vector<std::int64_t> owner(nt, -1);
        for (std::size_t i = 0; i < final.size(); ++i)
            for (std::size_t f : final[i].mTris)
                owner[f] = static_cast<std::int64_t>(i);
        std::vector<std::size_t> order(final.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::stable_sort(order.begin(), order.end(), [&](std::size_t A, std::size_t B) {
            return final[A].mTris.size() < final[B].mTris.size();
        });
        for (std::size_t i : order) {
            if (final[i].mTris.size() >= rOptions.mMinTriangles || final[i].mTris.empty())
                continue;
            std::map<std::int64_t, std::size_t> shared;
            for (std::size_t f : final[i].mTris)
                for (std::size_t k = 0; k < 3; ++k) {
                    const std::int64_t g = st.mMesh.mNbr[f][k];
                    if (g >= 0 &&
                        owner[static_cast<std::size_t>(g)] != static_cast<std::int64_t>(i))
                        ++shared[owner[static_cast<std::size_t>(g)]];
                }
            std::vector<std::pair<std::size_t, std::int64_t>> cand;
            for (const auto& kv : shared)
                cand.push_back({kv.second, kv.first});
            std::stable_sort(cand.begin(), cand.end(),
                             [](const auto& rA, const auto& rB) { return rA.first > rB.first; });
            for (const auto& c : cand) {
                const std::size_t j = static_cast<std::size_t>(c.second);
                std::vector<std::size_t> merged = final[j].mTris;
                merged.insert(merged.end(), final[i].mTris.begin(), final[i].mTris.end());
                std::sort(merged.begin(), merged.end());
                if (!patch_topology(ql_local(rTris, merged)).IsDisc())
                    continue;
                final[j].mTris = std::move(merged);
                final[i].mTris.clear();
                changed = true;
                break;
            }
            if (changed)
                break;
        }
        final.erase(std::remove_if(final.begin(), final.end(),
                                   [](const QlPiece& rP) { return rP.mTris.empty(); }),
                    final.end());
    }

    // Smooth the cuts: a triangle with two sides on the same neighbour (an
    // ear, whose apex has no vertex inside its piece to start a cut from)
    // moves to that neighbour, where both stay discs.
    {
        std::vector<std::int64_t> owner(nt, -1);
        for (std::size_t i = 0; i < final.size(); ++i)
            for (std::size_t f : final[i].mTris)
                owner[f] = static_cast<std::int64_t>(i);
        const auto members = [&](std::int64_t I) {
            std::vector<std::size_t> m;
            for (std::size_t f = 0; f < nt; ++f)
                if (owner[f] == I)
                    m.push_back(f);
            return m;
        };
        for (int pass = 0; pass < 4; ++pass) {
            bool changed = false;
            for (std::size_t f = 0; f < nt; ++f) {
                const std::int64_t i = owner[f];
                std::int64_t other = -1;
                int across = 0, inside = 0;
                for (std::size_t k = 0; k < 3; ++k) {
                    const std::int64_t g = st.mMesh.mNbr[f][k];
                    if (g < 0)
                        continue;
                    const std::int64_t o = owner[static_cast<std::size_t>(g)];
                    if (o == i) {
                        ++inside;
                    } else if (other < 0 || other == o) {
                        other = o;
                        ++across;
                    } else {
                        across = -9;  // two different neighbours: a junction
                    }
                }
                if (across < 2 || inside > 1 || other < 0)
                    continue;
                owner[f] = other;
                std::vector<std::size_t> from = members(i), to = members(other);
                if ((!from.empty() && !patch_topology(ql_local(rTris, from)).IsDisc()) ||
                    !patch_topology(ql_local(rTris, to)).IsDisc()) {
                    owner[f] = i;
                    continue;
                }
                changed = true;
            }
            if (!changed)
                break;
        }
        for (std::size_t i = 0; i < final.size(); ++i)
            final[i].mTris = members(static_cast<std::int64_t>(i));
        final.erase(std::remove_if(final.begin(), final.end(),
                                   [](const QlPiece& rP) { return rP.mTris.empty(); }),
                    final.end());
    }

    // 4. Corners. Junctions: vertices three or more pieces (or the region's
    // boundary) meet at.
    std::vector<std::int64_t> label(nt, -1);
    for (std::size_t i = 0; i < final.size(); ++i)
        for (std::size_t f : final[i].mTris)
            label[f] = static_cast<std::int64_t>(i);
    std::vector<std::uint8_t> junction(st.mMesh.mNumPts, 0);
    for (std::size_t v = 0; v < st.mMesh.mNumPts; ++v) {
        std::vector<std::int64_t> ids;
        for (std::size_t f : st.mMesh.mVertTris[v])
            ids.push_back(label[f]);
        std::sort(ids.begin(), ids.end());
        const std::size_t distinct =
            static_cast<std::size_t>(std::unique(ids.begin(), ids.end()) - ids.begin());
        if (distinct + (st.mMesh.mOnBoundary[v] ? 1 : 0) >= 3)
            junction[v] = 1;
    }
    for (const QlPiece& p : final)
        for (std::size_t v : p.mForced)
            st.mCorner[v] = 1;
    std::vector<QlQuad> done;
    for (std::size_t i = 0; i < final.size(); ++i) {
        const QlPiece p = final[i];  // `final` may grow below
        ql_mark(st, p.mTris, 1);
        const std::vector<std::size_t> b = ql_boundary(st.mMesh, rTris, p.mTris, st.mIn);
        ql_mark(st, p.mTris, 0);
        if (b.size() < 4) {
            log::debug("quad_layout: piece of {} triangles has a boundary of {} vertices",
                       p.mTris.size(), b.size());
            out.mFailure = QuadLayoutFailure::CornersFailed;
            return out;
        }
        std::vector<QlPiece> more;
        bool split = false;
        const std::vector<std::size_t> c = ql_corners(st, p, b, junction, more, split);
        if (split) {
            for (auto& m : more) {
                for (std::size_t v : m.mForced)
                    junction[v] = 1;
                final.push_back(std::move(m));
            }
            continue;
        }
        if (c.size() != 4) {
            log::debug("quad_layout: no four corners on a piece of {} triangles ({} forced)",
                       p.mTris.size(), p.mForced.size());
            out.mFailure = QuadLayoutFailure::CornersFailed;
            return out;
        }
        QlQuad d;
        d.mTris = p.mTris;
        // Rotate so the first corner leads.
        d.mBoundary.assign(b.begin() + static_cast<std::ptrdiff_t>(c[0]), b.end());
        d.mBoundary.insert(d.mBoundary.end(), b.begin(),
                           b.begin() + static_cast<std::ptrdiff_t>(c[0]));
        for (std::size_t j = 0; j < 4; ++j)
            d.mCorners[j] = c[j] - c[0];
        done.push_back(std::move(d));
    }
    std::sort(done.begin(), done.end(),
              [](const QlQuad& rA, const QlQuad& rB) { return rA.mTris[0] < rB.mTris[0]; });

    // Refine to the target size, conformingly (see `ql_refine`).
    if (rOptions.mTargetSize > 0.0)
        done = ql_refine(st, std::move(done), rOptions.mTargetSize);

    // 5. Assemble: pieces, arcs, chains, sides.
    out.mPatchOfTri.assign(nt, -1);
    out.mPatches.resize(done.size());
    std::vector<std::uint8_t> brk(st.mMesh.mNumPts, 0);
    for (std::size_t i = 0; i < done.size(); ++i) {
        QuadPatch& q = out.mPatches[i];
        q.mTris = done[i].mTris;
        q.mBoundary = done[i].mBoundary;
        for (std::size_t j = 0; j < 4; ++j) {
            q.mCornerPos[j] = done[i].mCorners[j];
            brk[q.mBoundary[q.mCornerPos[j]]] = 1;
        }
        for (std::size_t f : q.mTris)
            out.mPatchOfTri[f] = static_cast<std::int64_t>(i);
    }
    for (std::size_t v = 0; v < st.mMesh.mNumPts; ++v) {
        std::vector<std::int64_t> ids;
        for (std::size_t f : st.mMesh.mVertTris[v])
            ids.push_back(out.mPatchOfTri[f]);
        std::sort(ids.begin(), ids.end());
        const std::size_t distinct =
            static_cast<std::size_t>(std::unique(ids.begin(), ids.end()) - ids.begin());
        if (distinct + (st.mMesh.mOnBoundary[v] ? 1 : 0) >= 3)
            brk[v] = 1;
    }
    // Arcs, keyed by their vertices from the smaller end.
    struct QlUse {
        std::size_t mPatch, mSide;
        bool mForward;
    };
    std::map<std::vector<std::size_t>, std::size_t> arc_id;
    std::vector<std::vector<std::size_t>> arcs;
    std::vector<std::vector<QlUse>> uses;
    // Per patch side, its arcs in order and whether each runs forward.
    std::vector<std::array<std::vector<std::pair<std::size_t, bool>>, 4>> side_arcs(done.size());
    for (std::size_t i = 0; i < out.mPatches.size(); ++i) {
        const QuadPatch& q = out.mPatches[i];
        const std::size_t nb = q.mBoundary.size();
        for (std::size_t j = 0; j < 4; ++j) {
            const std::size_t a = q.mCornerPos[j];
            const std::size_t b = j + 1 < 4 ? q.mCornerPos[j + 1] : nb;
            std::vector<std::size_t> cur{q.mBoundary[a]};
            for (std::size_t t = a + 1; t <= b; ++t) {
                const std::size_t v = q.mBoundary[t % nb];
                cur.push_back(v);
                if (!brk[v] && t != b)
                    continue;
                std::vector<std::size_t> key = cur;
                const bool forward = key.front() < key.back();
                if (!forward)
                    std::reverse(key.begin(), key.end());
                auto it = arc_id.find(key);
                if (it == arc_id.end()) {
                    it = arc_id.emplace(key, arcs.size()).first;
                    arcs.push_back(key);
                    uses.emplace_back();
                }
                uses[it->second].push_back({i, j, forward});
                side_arcs[i][j].push_back({it->second, forward});
                cur.assign(1, v);
            }
        }
    }
    for (std::size_t a = 0; a < arcs.size(); ++a) {
        if (uses[a].size() == 2 && uses[a][0].mForward != uses[a][1].mForward &&
            uses[a][0].mPatch != uses[a][1].mPatch) {
            ++out.mNumInteriorArcs;
            continue;
        }
        bool boundary = uses[a].size() == 1;
        for (std::size_t t = 0; boundary && t + 1 < arcs[a].size(); ++t)
            boundary = st.mMesh.mOnBoundary[arcs[a][t]] && st.mMesh.mOnBoundary[arcs[a][t + 1]];
        if (!boundary) {
            out.mFailure = QuadLayoutFailure::Inconsistent;
            return out;
        }
    }
    // Chains: the arcs a side runs over belong together.
    QlUnionFind uf(arcs.size());
    for (const auto& sides : side_arcs)
        for (const auto& s : sides)
            for (std::size_t t = 1; t < s.size(); ++t)
                uf.Unite(s[0].first, s[t].first);
    std::map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t a = 0; a < arcs.size(); ++a)
        groups[uf.Find(a)].push_back(a);
    std::vector<std::size_t> chain_of(arcs.size(), 0);
    std::vector<std::unordered_map<std::size_t, std::size_t>> chain_pos;
    // The crease edges, for the chains along them: the cuts' features and
    // every edge sharper than the feature angle, eroded or not (a sharp edge
    // the erosion dropped as a cut is still no smooth curve).
    std::set<std::pair<std::size_t, std::size_t>> crease;
    for (std::size_t f = 0; f < nt; ++f)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t g = st.mMesh.mNbr[f][k];
            const bool sharp =
                g >= 0 && ql_angle_deg(normal[f], normal[static_cast<std::size_t>(g)]) >
                              rOptions.mFeatureAngleDeg;
            if (feature[f][k] || sharp)
                crease.insert({std::min(rTris[f][k], rTris[f][(k + 1) % 3]),
                               std::max(rTris[f][k], rTris[f][(k + 1) % 3])});
        }
    // The arcs `rG` as one simple path, or false.
    const auto make_chain = [&](const std::vector<std::size_t>& rG) {
        std::map<std::size_t, std::vector<std::size_t>> ends;  // vertex -> arcs ending there
        for (std::size_t a : rG) {
            ends[arcs[a].front()].push_back(a);
            ends[arcs[a].back()].push_back(a);
        }
        std::size_t start = kQlNone;
        for (const auto& e : ends) {
            if (e.second.size() > 2)
                return false;
            if (e.second.size() == 1 && start == kQlNone)
                start = e.first;
        }
        if (start == kQlNone)
            return false;
        QuadChain ch;
        ch.mVerts.push_back(start);
        ch.mOnBoundary = true;
        std::vector<std::uint8_t> taken(rG.size(), 0);
        std::vector<std::size_t> order;
        std::size_t at = start;
        for (std::size_t n = 0; n < rG.size(); ++n) {
            std::size_t pick = kQlNone;
            for (std::size_t a : ends[at]) {
                const std::size_t gi =
                    static_cast<std::size_t>(std::find(rG.begin(), rG.end(), a) - rG.begin());
                if (!taken[gi]) {
                    pick = gi;
                    break;
                }
            }
            if (pick == kQlNone)
                return false;
            taken[pick] = 1;
            const std::vector<std::size_t>& av = arcs[rG[pick]];
            if (av.front() == at)
                ch.mVerts.insert(ch.mVerts.end(), av.begin() + 1, av.end());
            else
                ch.mVerts.insert(ch.mVerts.end(), av.rbegin() + 1, av.rend());
            at = ch.mVerts.back();
            if (uses[rG[pick]].size() != 1)
                ch.mOnBoundary = false;
            order.push_back(rG[pick]);
        }
        std::unordered_map<std::size_t, std::size_t> where;
        for (std::size_t t = 0; t < ch.mVerts.size(); ++t)
            if (!where.emplace(ch.mVerts[t], t).second)
                return false;
        ch.mSmooth = !ch.mOnBoundary;
        for (std::size_t t = 0; ch.mSmooth && t + 1 < ch.mVerts.size(); ++t)
            if (crease.count({std::min(ch.mVerts[t], ch.mVerts[t + 1]),
                              std::max(ch.mVerts[t], ch.mVerts[t + 1])}) != 0)
                ch.mSmooth = false;
        for (std::size_t a : order)
            chain_of[a] = out.mChains.size();
        chain_pos.push_back(std::move(where));
        out.mChains.push_back(std::move(ch));
        return true;
    };
    for (const auto& kv : groups) {
        if (make_chain(kv.second))
            continue;
        for (std::size_t a : kv.second)
            if (!make_chain({a})) {
                out.mFailure = QuadLayoutFailure::CyclicChain;
                return out;
            }
    }
    out.mSides.resize(out.mPatches.size());
    for (std::size_t i = 0; i < out.mPatches.size(); ++i)
        for (std::size_t j = 0; j < 4; ++j) {
            std::vector<QuadSegment>& segs = out.mSides[i][j].mSegments;
            std::size_t first = kQlNone, last = kQlNone, chain = kQlNone;
            for (const auto& af : side_arcs[i][j]) {
                const std::vector<std::size_t>& av = arcs[af.first];
                const std::size_t a = af.second ? av.front() : av.back();
                const std::size_t b = af.second ? av.back() : av.front();
                if (chain_of[af.first] != chain) {
                    if (chain != kQlNone)
                        segs.push_back(
                            {chain, chain_pos[chain].at(first), chain_pos[chain].at(last)});
                    chain = chain_of[af.first];
                    first = a;
                }
                last = b;
            }
            segs.push_back({chain, chain_pos[chain].at(first), chain_pos[chain].at(last)});
        }
    // Extraordinary vertices and T-junctions.
    std::vector<int> corners(st.mMesh.mNumPts, 0), through(st.mMesh.mNumPts, 0);
    for (std::size_t i = 0; i < out.mPatches.size(); ++i) {
        const QuadPatch& q = out.mPatches[i];
        std::vector<std::uint8_t> is_corner(q.mBoundary.size(), 0);
        for (std::size_t j = 0; j < 4; ++j) {
            is_corner[q.mCornerPos[j]] = 1;
            ++corners[q.mBoundary[q.mCornerPos[j]]];
        }
        for (std::size_t t = 0; t < q.mBoundary.size(); ++t)
            if (!is_corner[t] && brk[q.mBoundary[t]])
                ++through[q.mBoundary[t]];
    }
    for (std::size_t v = 0; v < st.mMesh.mNumPts; ++v) {
        if (corners[v] > 0 && through[v] > 0)
            ++out.mNumTJunctions;
        if (st.mMesh.mOnBoundary[v] || corners[v] == 0 || through[v] > 0 || corners[v] == 4)
            continue;
        out.mExtraordinary.push_back(v);
        out.mValence.push_back(corners[v]);
    }
    out.mFailure = quad_layout_check(out, rTris);
    out.mOk = out.mFailure == QuadLayoutFailure::None;
    return out;
}

QuadLayoutFailure quad_layout_check(const QuadLayout& rLayout, const std::vector<PatchTri>& rTris) {
    const std::size_t nt = rTris.size();
    if (rLayout.mPatchOfTri.size() != nt || rLayout.mSides.size() != rLayout.mPatches.size())
        return QuadLayoutFailure::Inconsistent;
    std::size_t npts = 0;
    for (const PatchTri& t : rTris)
        npts = std::max({npts, t[0] + 1, t[1] + 1, t[2] + 1});
    // Directed boundary half-edges of every piece.
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> half;
    std::size_t covered = 0;
    for (std::size_t i = 0; i < rLayout.mPatches.size(); ++i) {
        const QuadPatch& q = rLayout.mPatches[i];
        covered += q.mTris.size();
        for (std::size_t f : q.mTris)
            if (f >= nt || rLayout.mPatchOfTri[f] != static_cast<std::int64_t>(i))
                return QuadLayoutFailure::Inconsistent;
        const std::vector<PatchTri> local = ql_local(rTris, q.mTris);
        if (!patch_topology(local).IsDisc())
            return QuadLayoutFailure::CutFailed;
        const std::size_t nb = q.mBoundary.size();
        if (nb < 4 || q.mCornerPos[0] != 0)
            return QuadLayoutFailure::CornersFailed;
        for (std::size_t j = 1; j < 4; ++j)
            if (q.mCornerPos[j] <= q.mCornerPos[j - 1] || q.mCornerPos[j] >= nb)
                return QuadLayoutFailure::CornersFailed;
        // The boundary must be the piece's: every triangle half-edge without
        // its twin in the piece, once.
        std::map<std::pair<std::size_t, std::size_t>, int> own;
        for (const PatchTri& t : local)
            for (std::size_t k = 0; k < 3; ++k)
                own[{t[k], t[(k + 1) % 3]}] += 1;
        std::size_t nbound = 0;
        for (const auto& kv : own)
            if (own.find({kv.first.second, kv.first.first}) == own.end())
                ++nbound;
        if (nbound != nb)
            return QuadLayoutFailure::Inconsistent;
        for (std::size_t t = 0; t < nb; ++t) {
            const std::pair<std::size_t, std::size_t> e{q.mBoundary[t], q.mBoundary[(t + 1) % nb]};
            if (own.find(e) == own.end() || own.find({e.second, e.first}) != own.end())
                return QuadLayoutFailure::Inconsistent;
            if (!half.emplace(e, i).second)
                return QuadLayoutFailure::Inconsistent;
        }
    }
    if (covered != nt)
        return QuadLayoutFailure::Inconsistent;
    // Region boundary half-edges.
    std::map<std::pair<std::size_t, std::size_t>, int> region;
    for (const PatchTri& t : rTris)
        for (std::size_t k = 0; k < 3; ++k)
            region[{t[k], t[(k + 1) % 3]}] += 1;
    for (const auto& kv : half) {
        const bool on_region = region.find({kv.first.second, kv.first.first}) == region.end();
        const auto twin = half.find({kv.first.second, kv.first.first});
        if (on_region ? twin != half.end() : (twin == half.end() || twin->second == kv.second))
            return QuadLayoutFailure::Inconsistent;
    }
    // Chains are simple, and every side is the stretch of its chain between
    // its corners.
    for (const QuadChain& c : rLayout.mChains) {
        std::vector<std::size_t> v = c.mVerts;
        std::sort(v.begin(), v.end());
        if (v.size() < 2 || std::adjacent_find(v.begin(), v.end()) != v.end())
            return QuadLayoutFailure::CyclicChain;
    }
    for (std::size_t i = 0; i < rLayout.mPatches.size(); ++i) {
        const QuadPatch& q = rLayout.mPatches[i];
        const std::size_t nb = q.mBoundary.size();
        for (std::size_t j = 0; j < 4; ++j) {
            const std::size_t a = q.mCornerPos[j];
            const std::size_t b = j + 1 < 4 ? q.mCornerPos[j + 1] : nb;
            std::size_t at = a;
            if (rLayout.mSides[i][j].mSegments.empty())
                return QuadLayoutFailure::Inconsistent;
            for (const QuadSegment& s : rLayout.mSides[i][j].mSegments) {
                if (s.mChain >= rLayout.mChains.size())
                    return QuadLayoutFailure::Inconsistent;
                const std::vector<std::size_t>& cv = rLayout.mChains[s.mChain].mVerts;
                if (s.mBegin >= cv.size() || s.mEnd >= cv.size() || s.mBegin == s.mEnd)
                    return QuadLayoutFailure::Inconsistent;
                const std::size_t span = s.mBegin < s.mEnd ? s.mEnd - s.mBegin : s.mBegin - s.mEnd;
                if (at + span > b)
                    return QuadLayoutFailure::Inconsistent;
                for (std::size_t t = 0; t <= span; ++t) {
                    const std::size_t cp = s.mBegin < s.mEnd ? s.mBegin + t : s.mBegin - t;
                    if (cv[cp] != q.mBoundary[(at + t) % nb])
                        return QuadLayoutFailure::Inconsistent;
                }
                at += span;
            }
            if (at != b)
                return QuadLayoutFailure::Inconsistent;
        }
    }
    (void)npts;
    return QuadLayoutFailure::None;
}

}  // namespace detail
}  // namespace meshioplusplus
