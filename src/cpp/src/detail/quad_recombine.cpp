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
#include <tuple>
#include <utility>
#include <vector>

// Project includes
#include "quad_recombine.hpp"

namespace meshioplusplus::detail {

namespace {

using QrV3 = std::array<double, 3>;

QrV3 qr_at(const double* pX, std::int64_t Id) {
    const auto i = static_cast<std::size_t>(Id);
    return {pX[3 * i], pX[3 * i + 1], pX[3 * i + 2]};
}
QrV3 qr_sub(const QrV3& rA, const QrV3& rB) { return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]}; }
QrV3 qr_cross(const QrV3& rA, const QrV3& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
double qr_dot(const QrV3& rA, const QrV3& rB) { return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2]; }
double qr_norm(const QrV3& rA) { return std::sqrt(qr_dot(rA, rA)); }

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

double quad_recombine_quality(const double* pX, const std::array<std::int64_t, 4>& rQ) {
    const auto p = qr_points(pX, rQ);
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

QuadRecombineResult quad_recombine(const std::vector<std::array<std::int64_t, 3>>& rTris,
                                   const std::vector<std::int64_t>& rGroup, const double* pX,
                                   double MinQuality) {
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
        const std::array<std::int64_t, 4> quad{a, d, b, c};
        const QrV3 n1 = normal_of(t1), n2 = normal_of(t2);
        const double l = qr_norm(n1) * qr_norm(n2);
        if (!(l > 0.0))
            continue;
        const double q = quad_recombine_quality(pX, quad) * std::max(0.0, qr_dot(n1, n2) / l);
        if (q >= MinQuality && q > 0.0)
            cand.push_back({q, a, b, std::min(t1, t2), std::max(t1, t2), quad});
    }
    std::sort(cand.begin(), cand.end(), [](const Candidate& rX, const Candidate& rY) {
        if (rX.mQuality != rY.mQuality)
            return rX.mQuality > rY.mQuality;
        return std::tie(rX.mA, rX.mB) < std::tie(rY.mA, rY.mB);
    });
    std::vector<char> used(rTris.size(), 0);
    std::vector<std::pair<std::size_t, std::size_t>> chosen;  // (first triangle, candidate)
    for (std::size_t i = 0; i < cand.size(); ++i) {
        const Candidate& r_c = cand[i];
        if (used[r_c.mT1] || used[r_c.mT2])
            continue;
        used[r_c.mT1] = used[r_c.mT2] = 1;
        chosen.emplace_back(r_c.mT1, i);
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
