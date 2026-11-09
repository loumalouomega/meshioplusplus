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
// Minimum-zone and envelope fits. See detail/zone_fit.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

// Project includes
#include "zone_fit.hpp"
#include "dense_lp.hpp"
#include "lsq_fit.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

using ZfV3 = std::array<double, 3>;

double zf_dot(const ZfV3& rA, const ZfV3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

ZfV3 zf_unit(ZfV3 V) {
    const double l = std::sqrt(zf_dot(V, V));
    if (l > 0.0)
        for (double& r_x : V)
            r_x /= l;
    return V;
}

// The fit table's sign: the largest-magnitude component positive.
ZfV3 zf_canonical(ZfV3 V) {
    std::size_t big = 0;
    for (std::size_t k = 1; k < 3; ++k)
        if (std::fabs(V[k]) > std::fabs(V[big]))
            big = k;
    if (V[big] < 0.0)
        for (double& r_x : V)
            r_x = -r_x;
    return V;
}

// The surface's current geometry: a point (plane point, sphere centre,
// cylinder axis point) and a direction (plane normal, cylinder axis).
struct ZfState {
    ZfV3 mPoint{};
    ZfV3 mDir{0.0, 0.0, 1.0};
};

struct ZfProblem {
    std::int64_t mType = 0;
    ZoneCriterion mCriterion = ZoneCriterion::MinimumZone;
    const double* mPts = nullptr;
    std::size_t mN = 0;
    double mExtent = 1.0;

    ZfV3 P(std::size_t I) const { return {mPts[3 * I], mPts[3 * I + 1], mPts[3 * I + 2]}; }

    // The signed distance of point I from the plane, or its radius about the
    // centre or the axis.
    double Dist(const ZfState& rS, std::size_t I) const {
        const ZfV3 w{mPts[3 * I] - rS.mPoint[0], mPts[3 * I + 1] - rS.mPoint[1],
                     mPts[3 * I + 2] - rS.mPoint[2]};
        if (mType == 1)
            return zf_dot(w, rS.mDir);
        if (mType == 2)
            return std::sqrt(zf_dot(w, w));
        const double h = zf_dot(w, rS.mDir);
        return std::sqrt(std::max(0.0, zf_dot(w, w) - h * h));
    }

    // The criterion's exact value, to minimise.
    double Objective(const ZfState& rS) const {
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (std::size_t i = 0; i < mN; ++i) {
            const double d = Dist(rS, i);
            lo = std::min(lo, d);
            hi = std::max(hi, d);
        }
        if (mCriterion == ZoneCriterion::MinimumZone)
            return hi - lo;
        return mCriterion == ZoneCriterion::MinCircumscribed ? hi : -lo;
    }

    // The number of geometric step variables.
    std::size_t NumGeo() const { return mType == 1 ? 2 : (mType == 2 ? 3 : 4); }

    // The linearised distance of point I: d0 + g . step (step in lengths).
    double Linear(const ZfState& rS, std::size_t I, const ZfV3& rU, const ZfV3& rV,
                  double* pGrad) const {
        const ZfV3 w{mPts[3 * I] - rS.mPoint[0], mPts[3 * I + 1] - rS.mPoint[1],
                     mPts[3 * I + 2] - rS.mPoint[2]};
        if (mType == 1) {
            // Tilting the normal by (a, b) / extent about the point.
            pGrad[0] = zf_dot(w, rU) / mExtent;
            pGrad[1] = zf_dot(w, rV) / mExtent;
            return zf_dot(w, rS.mDir);
        }
        if (mType == 2) {
            const double r = std::sqrt(zf_dot(w, w));
            for (std::size_t k = 0; k < 3; ++k)
                pGrad[k] = r > 0.0 ? -w[k] / r : 0.0;
            return r;
        }
        const double h = zf_dot(w, rS.mDir);
        ZfV3 q{w[0] - h * rS.mDir[0], w[1] - h * rS.mDir[1], w[2] - h * rS.mDir[2]};
        const double rho = std::sqrt(zf_dot(q, q));
        const ZfV3 e = rho > 0.0 ? ZfV3{q[0] / rho, q[1] / rho, q[2] / rho} : ZfV3{0.0, 0.0, 0.0};
        const double eu = zf_dot(e, rU), ev = zf_dot(e, rV);
        pGrad[0] = -eu;  // moving the axis point
        pGrad[1] = -ev;
        pGrad[2] = -h * eu / mExtent;  // tilting the axis
        pGrad[3] = -h * ev / mExtent;
        return rho;
    }

    ZfState Apply(const ZfState& rS, const ZfV3& rU, const ZfV3& rV, const double* pStep,
                  double Offset) const {
        ZfState n = rS;
        if (mType == 1) {
            for (std::size_t k = 0; k < 3; ++k)
                n.mDir[k] = rS.mDir[k] + (pStep[0] * rU[k] + pStep[1] * rV[k]) / mExtent;
            n.mDir = zf_unit(n.mDir);
            for (std::size_t k = 0; k < 3; ++k)
                n.mPoint[k] = rS.mPoint[k] + Offset * rS.mDir[k];
        } else if (mType == 2) {
            for (std::size_t k = 0; k < 3; ++k)
                n.mPoint[k] = rS.mPoint[k] + pStep[k];
        } else {
            for (std::size_t k = 0; k < 3; ++k) {
                n.mPoint[k] = rS.mPoint[k] + pStep[0] * rU[k] + pStep[1] * rV[k];
                n.mDir[k] = rS.mDir[k] + (pStep[2] * rU[k] + pStep[3] * rV[k]) / mExtent;
            }
            n.mDir = zf_unit(n.mDir);
        }
        return n;
    }
};

}  // namespace

ZoneFitResult zone_fit(std::int64_t Type, const double* pStart, const double* pPts, std::size_t N,
                       ZoneCriterion Criterion) {
    ZoneFitResult out;
    if (Type < 1 || Type > 3) {
        out.mWhy = "a minimum-zone or envelope fit applies to planes, spheres and cylinders";
        return out;
    }
    if (Type == 1 && Criterion != ZoneCriterion::MinimumZone) {
        out.mWhy = "a plane has a minimum zone, not an envelope";
        return out;
    }
    if (N < 4) {
        out.mWhy = "too few points";
        return out;
    }
    ZfProblem pb;
    pb.mType = Type;
    pb.mCriterion = Criterion;
    pb.mPts = pPts;
    pb.mN = N;
    ZfV3 lo{pPts[0], pPts[1], pPts[2]}, hi = lo, mean{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], pPts[3 * i + k]);
            hi[k] = std::max(hi[k], pPts[3 * i + k]);
            mean[k] += pPts[3 * i + k] / static_cast<double>(N);
        }
    const ZfV3 ext{hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
    pb.mExtent = std::sqrt(zf_dot(ext, ext));
    if (!(pb.mExtent > 0.0) || !std::isfinite(pb.mExtent)) {
        out.mWhy = "the points have no extent";
        return out;
    }
    ZfState s;
    s.mPoint = {pStart[0], pStart[1], pStart[2]};
    if (Type != 2)
        s.mDir = zf_unit({pStart[3], pStart[4], pStart[5]});
    const ZfV3 start_dir = s.mDir;
    double f = pb.Objective(s);
    double trust = 0.1 * pb.mExtent;
    const std::size_t g = pb.NumGeo();
    // Variables: the geometric step, then (plane) the offset and the width,
    // or (minimum zone) the inner and outer radius, or the envelope radius.
    const std::size_t nrad = Criterion == ZoneCriterion::MinimumZone ? 2 : 1;
    const std::size_t k = g + nrad;
    std::vector<double> grad(g);
    for (std::int64_t round = 0; round < 200 && trust > 1e-14 * pb.mExtent; ++round) {
        ++out.mRounds;
        ZfV3 u{}, v{};
        lsq_frame(s.mDir.data(), u.data(), v.data());
        std::vector<double> c(k, 0.0), a, b;
        a.reserve((2 * N + 2 * g) * k);
        // Minimum zone: minimise outer - inner (a plane: the width, with the
        // offset as the "inner" variable and lo = offset - w/2 recast as
        // inner <= d <= outer); envelopes: minimise R, or maximise it.
        if (Criterion == ZoneCriterion::MinimumZone) {
            c[g] = -1.0;
            c[g + 1] = 1.0;
        } else {
            c[g] = Criterion == ZoneCriterion::MinCircumscribed ? 1.0 : -1.0;
        }
        for (std::size_t i = 0; i < N; ++i) {
            const double d0 = pb.Linear(s, i, u, v, grad.data());
            const bool upper = Criterion != ZoneCriterion::MaxInscribed;
            const bool lower = Criterion != ZoneCriterion::MinCircumscribed;
            if (upper) {  // d0 + grad.step <= outer
                for (std::size_t j = 0; j < g; ++j)
                    a.push_back(grad[j]);
                if (Criterion == ZoneCriterion::MinimumZone) {
                    a.push_back(0.0);
                    a.push_back(-1.0);
                } else {
                    a.push_back(-1.0);
                }
                b.push_back(-d0);
            }
            if (lower) {  // inner <= d0 + grad.step
                for (std::size_t j = 0; j < g; ++j)
                    a.push_back(-grad[j]);
                a.push_back(1.0);
                if (Criterion == ZoneCriterion::MinimumZone)
                    a.push_back(0.0);
                b.push_back(d0);
            }
        }
        for (std::size_t j = 0; j < g; ++j)
            for (const double sgn : {1.0, -1.0}) {
                for (std::size_t q = 0; q < k; ++q)
                    a.push_back(q == j ? sgn : 0.0);
                b.push_back(trust);
            }
        const DenseLpResult lp = dense_lp_inequality(c, a, b);
        if (lp.mStatus != DenseLpStatus::Optimal) {
            trust *= 0.25;
            continue;
        }
        // A plane's offset: the zone's middle moves by (inner + outer) / 2.
        const double offset = Type == 1 ? 0.5 * (lp.mX[g] + lp.mX[g + 1]) : 0.0;
        const ZfState cand = pb.Apply(s, u, v, lp.mX.data(), offset);
        const double fc = pb.Objective(cand);
        double step = 0.0;
        for (std::size_t j = 0; j < g; ++j)
            step = std::max(step, std::fabs(lp.mX[j]));
        if (fc < f) {
            s = cand;
            f = fc;
            if (step < 1e-14 * pb.mExtent)
                break;
        } else {
            trust = std::min(trust, step) * 0.25;
        }
    }

    // The fit table's parameters.
    double dlo = std::numeric_limits<double>::infinity(), dhi = -dlo;
    if (Type == 1) {
        if (zf_dot(s.mDir, start_dir) < 0.0)
            for (double& r_x : s.mDir)
                r_x = -r_x;
        for (std::size_t i = 0; i < N; ++i) {
            const ZfV3 p = pb.P(i);
            const double d = zf_dot({p[0] - mean[0], p[1] - mean[1], p[2] - mean[2]}, s.mDir);
            dlo = std::min(dlo, d);
            dhi = std::max(dhi, d);
        }
        const double mid = 0.5 * (dlo + dhi);
        for (std::size_t q = 0; q < 3; ++q) {
            out.mParams[q] = mean[q] + mid * s.mDir[q];
            out.mParams[3 + q] = s.mDir[q];
        }
        out.mValue = dhi - dlo;
        out.mOk = true;
        return out;
    }
    if (Type == 3) {
        s.mDir = zf_canonical(s.mDir);
        const double t =
            zf_dot({mean[0] - s.mPoint[0], mean[1] - s.mPoint[1], mean[2] - s.mPoint[2]}, s.mDir);
        for (std::size_t q = 0; q < 3; ++q)
            s.mPoint[q] += t * s.mDir[q];
    }
    for (std::size_t i = 0; i < N; ++i) {
        const double d = pb.Dist(s, i);
        dlo = std::min(dlo, d);
        dhi = std::max(dhi, d);
    }
    for (std::size_t q = 0; q < 3; ++q) {
        out.mParams[q] = s.mPoint[q];
        out.mParams[3 + q] = Type == 3 ? s.mDir[q] : 0.0;
    }
    if (Criterion == ZoneCriterion::MinimumZone) {
        out.mParams[6] = 0.5 * (dlo + dhi);
        out.mValue = dhi - dlo;
    } else if (Criterion == ZoneCriterion::MinCircumscribed) {
        out.mParams[6] = dhi;
        out.mValue = dhi;
    } else {
        out.mParams[6] = dlo;
        out.mValue = dlo;
    }
    out.mOk = true;
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
