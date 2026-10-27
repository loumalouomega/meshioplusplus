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
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/parallel.hpp"
#include "quad_recombine.hpp"
#include "weighted_matching.hpp"

namespace meshioplusplus::detail {

namespace {

using QrV3 = std::array<double, 3>;

QrV3 qr_at(const double* pX, std::int64_t Id) {
    const auto i = static_cast<std::size_t>(Id);
    return {pX[3 * i], pX[3 * i + 1], pX[3 * i + 2]};
}
QrV3 qr_sub(const QrV3& rA, const QrV3& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}
QrV3 qr_cross(const QrV3& rA, const QrV3& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
double qr_dot(const QrV3& rA, const QrV3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
double qr_norm(const QrV3& rA) {
    return std::sqrt(qr_dot(rA, rA));
}

/// The quad's unit normal (its vector area normalised); zero when none.
QrV3 qr_normal(const std::array<QrV3, 4>& rP) {
    QrV3 n{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < 4; ++i) {
        const QrV3 c = qr_cross(rP[i], rP[(i + 1) % 4]);
        for (std::size_t k = 0; k < 3; ++k)
            n[k] += c[k];
    }
    const double l = qr_norm(n);
    if (!(l > 0.0))
        return {0.0, 0.0, 0.0};
    return {n[0] / l, n[1] / l, n[2] / l};
}

std::array<QrV3, 4> qr_points(const double* pX, const std::array<std::int64_t, 4>& rQ) {
    return {qr_at(pX, rQ[0]), qr_at(pX, rQ[1]), qr_at(pX, rQ[2]), qr_at(pX, rQ[3])};
}

}  // namespace

double quad_scaled_jacobian(const double* pX, const std::array<std::int64_t, 4>& rQ) {
    const auto p = qr_points(pX, rQ);
    const QrV3 n = qr_normal(p);
    double sj = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < 4; ++i) {
        const QrV3 e1 = qr_sub(p[(i + 1) % 4], p[i]);
        const QrV3 e2 = qr_sub(p[(i + 3) % 4], p[i]);
        const double l = qr_norm(e1) * qr_norm(e2);
        if (!(l > 0.0))
            return 0.0;
        sj = std::min(sj, qr_dot(qr_cross(e1, e2), n) / l);
    }
    return sj;
}

namespace {

/// `quad_recombine_quality` of a quad given by its corners.
double qr_quality(const std::array<QrV3, 4>& rP) {
    const auto& p = rP;
    const QrV3 n = qr_normal(p);
    double q = 1.0;
    for (std::size_t i = 0; i < 4; ++i) {
        const QrV3 e1 = qr_sub(p[(i + 1) % 4], p[i]);
        const QrV3 e2 = qr_sub(p[(i + 3) % 4], p[i]);
        const double s = qr_dot(qr_cross(e1, e2), n);
        if (!(s > 0.0))
            return 0.0;  // a reflex or flat corner
        const double deg = std::atan2(s, qr_dot(e1, e2)) * (180.0 / 3.14159265358979323846);
        q = std::min(q, 1.0 - std::fabs(deg - 90.0) / 90.0);
    }
    return std::max(0.0, q);
}

}  // namespace

double quad_recombine_quality(const double* pX, const std::array<std::int64_t, 4>& rQ) {
    return qr_quality(qr_points(pX, rQ));
}

namespace {

using QrEdge = MatchingEdge;

/// Whether quad `rQ`, split through its edges' midpoints `rMid` and their
/// mean into four, has all four meet the floor `MinQuality`: a corner the
/// floor admits near a straight angle (150 degrees on the sphere) turns
/// past it once the midpoints lie on a curved surface and the centre is
/// lifted onto it (all-quad, roadmap 6.8). A side without a midpoint leaves
/// the quad as it is.
bool qr_splits_well(const double* pX, const std::array<std::int64_t, 4>& rQ,
                    const std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t>& rMid,
                    double MinQuality) {
    std::array<QrV3, 4> mid{};
    QrV3 centre{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < 4; ++i) {
        const auto it = rMid.find(std::minmax(rQ[i], rQ[(i + 1) % 4]));
        if (it == rMid.end())
            return true;
        mid[i] = qr_at(pX, it->second);
        for (std::size_t k = 0; k < 3; ++k)
            centre[k] += 0.25 * mid[i][k];
    }
    for (std::size_t i = 0; i < 4; ++i)
        if (!(qr_quality({qr_at(pX, rQ[i]), mid[i], centre, mid[(i + 3) % 4]}) >= MinQuality))
            return false;
    return true;
}

// Triangles a face is matched in at once: the matching is exact inside a
// patch, and its cost grows faster than linearly with the patch (measured on
// the sphere at h = 0.025, 45,000 triangles: 45 s whole, 2.6 s in patches of
// 3,000 with a band of 4 rings, leaving 3,074 triangles where the exact
// matching leaves 2,624 and the greedy one 4,796).
constexpr std::int64_t kQrPatch = 3000;
// Rings about the triangles left over that are matched again together.
constexpr int kQrBandRings = 4;

// The matching of the induced subgraph on `rIds` (ascending local ids),
// written into `rMate`; edges leaving it are ignored.
void qr_match_subset(const std::vector<std::int64_t>& rIds,
                     const std::vector<std::vector<std::size_t>>& rIncident,
                     const std::vector<QrEdge>& rEdges, std::vector<std::int64_t>& rMate) {
    std::unordered_map<std::int64_t, std::int64_t> sub;
    for (std::size_t i = 0; i < rIds.size(); ++i)
        sub.emplace(rIds[i], static_cast<std::int64_t>(i));
    std::vector<MatchingEdge> e;
    for (const std::int64_t v : rIds)
        for (const std::size_t k : rIncident[static_cast<std::size_t>(v)]) {
            const QrEdge& r_e = rEdges[k];
            const std::int64_t o = r_e.mA == v ? r_e.mB : r_e.mA;
            if (o > v) {
                const auto it = sub.find(o);
                if (it != sub.end())
                    e.push_back({sub.at(v), it->second, r_e.mWeight});
            }
        }
    const auto m = max_weight_matching(static_cast<std::int64_t>(rIds.size()), e, true);
    for (std::size_t i = 0; i < rIds.size(); ++i)
        rMate[static_cast<std::size_t>(rIds[i])] =
            m[i] < 0 ? -1 : rIds[static_cast<std::size_t>(m[i])];
}

// A maximum-cardinality, maximum-weight matching, exact up to `kQrPatch`
// triangles; above, exact in connected patches of that size grown
// breadth-first in index order, then exact again on the band of
// `kQrBandRings` rings about the triangles the patches left unmatched (with
// the pairs inside it released), which finds the short augmenting paths
// across the patches' seams.
std::vector<std::int64_t> qr_patched_matching(std::int64_t N, const std::vector<QrEdge>& rEdges) {
    if (N <= kQrPatch)
        return max_weight_matching(N, rEdges, true);
    std::vector<std::vector<std::size_t>> incident(static_cast<std::size_t>(N));
    for (std::size_t k = 0; k < rEdges.size(); ++k) {
        incident[static_cast<std::size_t>(rEdges[k].mA)].push_back(k);
        incident[static_cast<std::size_t>(rEdges[k].mB)].push_back(k);
    }
    const auto other = [&](std::size_t K, std::int64_t V) {
        return rEdges[K].mA == V ? rEdges[K].mB : rEdges[K].mA;
    };
    std::vector<std::int64_t> mate(static_cast<std::size_t>(N), -1);
    std::vector<char> seen(static_cast<std::size_t>(N), 0);
    for (std::int64_t start = 0; start < N; ++start) {
        if (seen[static_cast<std::size_t>(start)])
            continue;
        std::vector<std::int64_t> patch{start};
        seen[static_cast<std::size_t>(start)] = 1;
        for (std::size_t q = 0;
             q < patch.size() && static_cast<std::int64_t>(patch.size()) < kQrPatch; ++q)
            for (const std::size_t k : incident[static_cast<std::size_t>(patch[q])]) {
                const std::int64_t o = other(k, patch[q]);
                if (!seen[static_cast<std::size_t>(o)] &&
                    static_cast<std::int64_t>(patch.size()) < kQrPatch) {
                    seen[static_cast<std::size_t>(o)] = 1;
                    patch.push_back(o);
                }
            }
        std::sort(patch.begin(), patch.end());
        qr_match_subset(patch, incident, rEdges, mate);
    }
    // The band about the triangles left over that could still pair.
    std::vector<int> ring(static_cast<std::size_t>(N), -1);
    std::vector<std::int64_t> frontier;
    for (std::int64_t v = 0; v < N; ++v)
        if (mate[static_cast<std::size_t>(v)] < 0 &&
            !incident[static_cast<std::size_t>(v)].empty()) {
            ring[static_cast<std::size_t>(v)] = 0;
            frontier.push_back(v);
        }
    for (int r = 1; r <= kQrBandRings && !frontier.empty(); ++r) {
        std::vector<std::int64_t> next;
        for (const std::int64_t v : frontier)
            for (const std::size_t k : incident[static_cast<std::size_t>(v)]) {
                const std::int64_t o = other(k, v);
                if (ring[static_cast<std::size_t>(o)] < 0) {
                    ring[static_cast<std::size_t>(o)] = r;
                    next.push_back(o);
                }
            }
        frontier = std::move(next);
    }
    // A pair half inside the band stays; both halves inside are released.
    std::vector<std::int64_t> band;
    for (std::int64_t v = 0; v < N; ++v) {
        if (ring[static_cast<std::size_t>(v)] < 0)
            continue;
        const std::int64_t w = mate[static_cast<std::size_t>(v)];
        if (w < 0 || ring[static_cast<std::size_t>(w)] >= 0)
            band.push_back(v);
    }
    if (!band.empty())
        qr_match_subset(band, incident, rEdges, mate);
    return mate;
}

}  // namespace

QuadRecombineResult quad_recombine(
    const std::vector<std::array<std::int64_t, 3>>& rTris, const std::vector<std::int64_t>& rGroup,
    const double* pX, double MinQuality, QuadPairing Pairing,
    const std::vector<std::array<std::int64_t, 3>>* pSideCurve,
    const std::set<std::pair<std::int64_t, std::int64_t>>* pFixed,
    const std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t>* pMidpoint) {
    QuadRecombineResult out;
    // Each directed edge's triangle; an interior edge has both directions.
    std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> of_edge;
    for (std::size_t t = 0; t < rTris.size(); ++t)
        for (std::size_t k = 0; k < 3; ++k)
            of_edge.emplace(std::make_pair(rTris[t][k], rTris[t][(k + 1) % 3]), t);
    const auto normal_of = [&](std::size_t T) {
        const QrV3 a = qr_at(pX, rTris[T][0]);
        return qr_cross(qr_sub(qr_at(pX, rTris[T][1]), a), qr_sub(qr_at(pX, rTris[T][2]), a));
    };
    struct Candidate {
        double mQuality;
        std::int64_t mA, mB;
        std::size_t mT1, mT2;
        std::array<std::int64_t, 4> mQuad;
    };
    std::vector<Candidate> cand;
    for (const auto& [e, t1] : of_edge) {
        const std::int64_t a = e.first, b = e.second;
        if (a > b)
            continue;  // each interior edge once, from its lower end
        const auto it = of_edge.find({b, a});
        if (it == of_edge.end())
            continue;
        const std::size_t t2 = it->second;
        if (t1 == t2 || rGroup[t1] != rGroup[t2])
            continue;
        if (pFixed && pFixed->count({a, b}) > 0)
            continue;
        // The third points: c of (a, b, c), d of (b, a, d).
        std::int64_t c = -1, d = -1;
        for (std::size_t k = 0; k < 3; ++k) {
            if (rTris[t1][k] != a && rTris[t1][k] != b)
                c = rTris[t1][k];
            if (rTris[t2][k] != a && rTris[t2][k] != b)
                d = rTris[t2][k];
        }
        if (c < 0 || d < 0 || c == d)
            continue;
        if (pSideCurve) {
            // The curve of triangle T's side from P to Q.
            const auto curve = [&](std::size_t T, std::int64_t P, std::int64_t Q) {
                for (std::size_t k = 0; k < 3; ++k)
                    if (rTris[T][k] == P && rTris[T][(k + 1) % 3] == Q)
                        return (*pSideCurve)[T][k];
                return std::int64_t{-1};
            };
            // Corner a: sides (a, d) of t2 and (c, a) of t1; corner b: (d, b)
            // of t2 and (b, c) of t1.
            const std::int64_t ad = curve(t2, a, d), ca = curve(t1, c, a);
            const std::int64_t db = curve(t2, d, b), bc = curve(t1, b, c);
            if ((ad >= 0 && ad == ca) || (db >= 0 && db == bc))
                continue;
        }
        const std::array<std::int64_t, 4> quad{a, d, b, c};
        const QrV3 n1 = normal_of(t1), n2 = normal_of(t2);
        const double l = qr_norm(n1) * qr_norm(n2);
        if (!(l > 0.0))
            continue;
        const double q = quad_recombine_quality(pX, quad) * std::max(0.0, qr_dot(n1, n2) / l);
        if (!(q >= MinQuality && q > 0.0))
            continue;
        if (pMidpoint && !qr_splits_well(pX, quad, *pMidpoint, MinQuality))
            continue;
        cand.push_back({q, a, b, std::min(t1, t2), std::max(t1, t2), quad});
    }
    std::sort(cand.begin(), cand.end(), [](const Candidate& rX, const Candidate& rY) {
        if (rX.mQuality != rY.mQuality)
            return rX.mQuality > rY.mQuality;
        return std::tie(rX.mA, rX.mB) < std::tie(rY.mA, rY.mB);
    });
    std::vector<char> used(rTris.size(), 0);
    std::vector<std::pair<std::size_t, std::size_t>> chosen;  // (first triangle, candidate)
    if (Pairing == QuadPairing::Matching) {
        // Each group's triangles and candidates (a pair of triangles once,
        // its best), matched group by group: the groups are independent.
        std::map<std::int64_t, std::size_t> group_index;
        for (std::size_t t = 0; t < rTris.size(); ++t)
            group_index.emplace(rGroup[t], group_index.size());
        const std::size_t ng = group_index.size();
        std::vector<std::vector<std::size_t>> tris_of(ng), cands_of(ng);
        std::vector<std::int64_t> local(rTris.size(), -1);
        for (std::size_t t = 0; t < rTris.size(); ++t) {
            auto& r_t = tris_of[group_index.at(rGroup[t])];
            local[t] = static_cast<std::int64_t>(r_t.size());
            r_t.push_back(t);
        }
        std::set<std::pair<std::size_t, std::size_t>> seen;
        for (std::size_t i = 0; i < cand.size(); ++i)
            if (seen.emplace(cand[i].mT1, cand[i].mT2).second)
                cands_of[group_index.at(rGroup[cand[i].mT1])].push_back(i);
        std::vector<std::vector<std::size_t>> picked(ng);
        parallel_for(
            ng,
            [&](std::size_t G) {
                std::vector<QrEdge> edges;
                for (const std::size_t i : cands_of[G])
                    edges.push_back(
                        {local[cand[i].mT1], local[cand[i].mT2],
                         static_cast<std::int64_t>(std::llround(1e6 * cand[i].mQuality))});
                const auto mate =
                    qr_patched_matching(static_cast<std::int64_t>(tris_of[G].size()), edges);
                for (const std::size_t i : cands_of[G])
                    if (mate[static_cast<std::size_t>(local[cand[i].mT1])] == local[cand[i].mT2])
                        picked[G].push_back(i);
            },
            1);
        for (const auto& r_p : picked)
            for (const std::size_t i : r_p) {
                used[cand[i].mT1] = used[cand[i].mT2] = 1;
                chosen.emplace_back(cand[i].mT1, i);
            }
    } else {
        for (std::size_t i = 0; i < cand.size(); ++i) {
            const Candidate& r_c = cand[i];
            if (used[r_c.mT1] || used[r_c.mT2])
                continue;
            used[r_c.mT1] = used[r_c.mT2] = 1;
            chosen.emplace_back(r_c.mT1, i);
        }
    }
    std::sort(chosen.begin(), chosen.end());
    for (const auto& [t, i] : chosen) {
        out.mQuads.push_back(cand[i].mQuad);
        // The pair as (a, b, c), (b, a, d): the triangle holding a -> b first.
        const std::size_t ta = of_edge.at({cand[i].mA, cand[i].mB});
        out.mPairs.push_back({ta, ta == cand[i].mT1 ? cand[i].mT2 : cand[i].mT1});
    }
    for (std::size_t t = 0; t < rTris.size(); ++t)
        if (!used[t])
            out.mLeft.push_back(t);
    return out;
}

}  // namespace meshioplusplus::detail
