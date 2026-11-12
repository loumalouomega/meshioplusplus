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
#include <limits>
#include <numeric>
#include <string>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"
#include "cut_graph.hpp"
#include "freeform_adaptive.hpp"
#include "kinematic_fit.hpp"
#include "lsq_fit.hpp"

namespace meshioplusplus::detail {

namespace {

constexpr double kKfPi = 3.14159265358979323846;

using KfXy = std::array<double, 2>;

// The patch's used vertices, their area weights and centre.
struct KfPatch {
    std::vector<std::size_t> mLive;
    std::vector<double> mArea;
    std::vector<PatchVec> mNormals;
    PatchVec mCentre{0.0, 0.0, 0.0};
    double mScale = 0.0;
};

KfPatch kf_patch(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris) {
    KfPatch p;
    std::vector<std::uint8_t> interior;
    feature_vertex_normals(rPts, rTris, p.mNormals, interior);
    p.mArea.assign(rPts.size(), 0.0);
    std::vector<std::uint8_t> used(rPts.size(), 0);
    for (const PatchTri& t : rTris) {
        const double a = 0.5 * vec3_norm(vec3_cross(vec3_sub(rPts[t[1]], rPts[t[0]]),
                                                    vec3_sub(rPts[t[2]], rPts[t[0]])));
        for (std::size_t v : t) {
            p.mArea[v] += a / 3.0;
            used[v] = 1;
        }
    }
    double wsum = 0.0;
    for (std::size_t v = 0; v < rPts.size(); ++v)
        if (used[v]) {
            p.mLive.push_back(v);
            wsum += p.mArea[v];
            for (std::size_t c = 0; c < 3; ++c)
                p.mCentre[c] += p.mArea[v] * rPts[v][c];
        }
    if (wsum > 0.0)
        for (double& c : p.mCentre)
            c /= wsum;
    for (std::size_t v : p.mLive)
        p.mScale = std::max(p.mScale, vec3_norm(vec3_sub(rPts[v], p.mCentre)));
    return p;
}

// The length of the polyline through the bins' means of the points ordered
// by `rKey`, and those means.
double kf_binned(const std::vector<KfXy>& rXy, const std::vector<double>& rKey,
                 std::vector<KfXy>& rMeans) {
    std::vector<std::size_t> order(rXy.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t A, std::size_t B) { return rKey[A] < rKey[B]; });
    const std::size_t nb = std::clamp<std::size_t>(rXy.size() / 24, 8, 96);
    rMeans.clear();
    for (std::size_t b = 0; b < nb; ++b) {
        const std::size_t i0 = b * order.size() / nb, i1 = (b + 1) * order.size() / nb;
        if (i1 <= i0)
            continue;
        KfXy m{0.0, 0.0};
        for (std::size_t i = i0; i < i1; ++i) {
            m[0] += rXy[order[i]][0];
            m[1] += rXy[order[i]][1];
        }
        m[0] /= static_cast<double>(i1 - i0);
        m[1] /= static_cast<double>(i1 - i0);
        rMeans.push_back(m);
    }
    double len = 0.0;
    for (std::size_t i = 1; i < rMeans.size(); ++i)
        len += std::hypot(rMeans[i][0] - rMeans[i - 1][0], rMeans[i][1] - rMeans[i - 1][1]);
    return len;
}

// A planar profile through scattered points (see the file comment); false,
// with the reason, for a line, a circular arc or a profile no curve
// follows.
bool kf_profile(const std::vector<KfXy>& rXy, double Tol, BsplineCurve& rOut, std::string& rWhy) {
    const std::size_t n = rXy.size();
    if (n < 16) {
        rWhy = "too few points";
        return false;
    }
    KfXy c{0.0, 0.0};
    for (const KfXy& p : rXy) {
        c[0] += p[0];
        c[1] += p[1];
    }
    c[0] /= static_cast<double>(n);
    c[1] /= static_cast<double>(n);
    // The principal direction, and the line test.
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (const KfXy& p : rXy) {
        const double dx = p[0] - c[0], dy = p[1] - c[1];
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    const double ang = 0.5 * std::atan2(2.0 * sxy, sxx - syy);
    const double ux = std::cos(ang), uy = std::sin(ang);
    double line = 0.0;
    std::vector<double> along(n), around(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = rXy[i][0] - c[0], dy = rXy[i][1] - c[1];
        along[i] = dx * ux + dy * uy;
        line = std::max(line, std::abs(-dx * uy + dy * ux));
        around[i] = std::atan2(dy, dx);
    }
    if (line <= 0.4 * Tol) {
        rWhy = "the profile is a line";
        return false;
    }
    {
        std::vector<double> xs(n), ys(n);
        for (std::size_t i = 0; i < n; ++i) {
            xs[i] = rXy[i][0];
            ys[i] = rXy[i][1];
        }
        double cx = 0.0, cy = 0.0, r = 0.0;
        if (lsq_circle_kasa(xs, ys, cx, cy, r) && lsq_circle_geometric(xs, ys, cx, cy, r)) {
            double worst = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                worst = std::max(worst, std::abs(std::hypot(xs[i] - cx, ys[i] - cy) - r));
            if (worst <= 0.4 * Tol) {
                rWhy = "the profile is a circular arc";
                return false;
            }
        }
    }
    // The angle about the centre, unwrapped at its largest gap.
    {
        std::vector<double> sorted = around;
        std::sort(sorted.begin(), sorted.end());
        double gap = sorted.front() + 2.0 * kKfPi - sorted.back(), cut = sorted.front();
        for (std::size_t i = 1; i < sorted.size(); ++i)
            if (sorted[i] - sorted[i - 1] > gap) {
                gap = sorted[i] - sorted[i - 1];
                cut = sorted[i];
            }
        for (double& a : around)
            if (a < cut)
                a += 2.0 * kKfPi;
    }
    std::vector<KfXy> by_along, by_around;
    const double la = kf_binned(rXy, along, by_along);
    const double lr = kf_binned(rXy, around, by_around);
    std::vector<KfXy>& means = lr < la ? by_around : by_along;
    if (means.size() < 4) {
        rWhy = "too few bins";
        return false;
    }
    // Carried a tenth of its length past each end, so the face runs on past
    // its patch and the end points project inside.
    const double len = std::min(la, lr);
    const auto extend = [&](const KfXy& rEnd, const KfXy& rNext) {
        const double dx = rEnd[0] - rNext[0], dy = rEnd[1] - rNext[1];
        const double l = std::hypot(dx, dy);
        return l > 0.0 ? KfXy{rEnd[0] + 0.1 * len * dx / l, rEnd[1] + 0.1 * len * dy / l} : rEnd;
    };
    std::vector<PatchVec> poly;
    const KfXy head = extend(means.front(), means[1]);
    const KfXy tail = extend(means.back(), means[means.size() - 2]);
    poly.push_back({head[0], head[1], 0.0});
    for (const KfXy& m : means)
        poly.push_back({m[0], m[1], 0.0});
    poly.push_back({tail[0], tail[1], 0.0});
    CurveFitOptions co;
    co.mDegree = 3;
    co.mTolerance = 0.2 * Tol;
    co.mMaxPoles = 48;
    CurveFitResult fit = bspline_curve_fit(poly, co);
    if (!fit.mOk) {
        co.mTolerance = 0.35 * Tol;
        fit = bspline_curve_fit(poly, co);
    }
    if (fit.mCurve.mPoles.size() < 4) {
        rWhy = "no profile curve";
        return false;
    }
    // Every point near the curve, not just the bins' means.
    double worst = 0.0;
    for (const KfXy& p : rXy) {
        double t = 0.0;
        worst = std::max(worst, bspline_curve_project(fit.mCurve, {p[0], p[1], 0.0}, t));
    }
    if (!(worst <= 0.5 * Tol)) {
        rWhy = "the profile curve misses its points";
        return false;
    }
    rOut = std::move(fit.mCurve);
    return true;
}

// A patch in a revolution's meridian half-plane: `(rho, z)` per used point,
// its angles unwrapped at their largest gap, and their range.
struct KfMeridian {
    std::vector<KfXy> mRz;
    double mA0 = 0.0, mA1 = 0.0, mGap = 0.0, mRhoMin = 0.0;
};

KfMeridian kf_meridian(const std::vector<PatchVec>& rPts, const std::vector<std::size_t>& rLive,
                       const KinematicFit& rF) {
    KfMeridian m;
    std::vector<double> theta;
    m.mRhoMin = std::numeric_limits<double>::infinity();
    for (std::size_t v : rLive) {
        const PatchVec r = vec3_sub(rPts[v], rF.mOrigin);
        const double z = vec3_dot(r, rF.mAxis);
        const double x = vec3_dot(r, rF.mE1), y = vec3_dot(r, rF.mE2);
        m.mRz.push_back({std::hypot(x, y), z});
        m.mRhoMin = std::min(m.mRhoMin, std::hypot(x, y));
        theta.push_back(std::atan2(y, x));
    }
    std::vector<double> sorted = theta;
    std::sort(sorted.begin(), sorted.end());
    double cut = sorted.front();
    m.mGap = sorted.front() + 2.0 * kKfPi - sorted.back();
    for (std::size_t i = 1; i < sorted.size(); ++i)
        if (sorted[i] - sorted[i - 1] > m.mGap) {
            m.mGap = sorted[i] - sorted[i - 1];
            cut = sorted[i];
        }
    m.mA0 = std::numeric_limits<double>::infinity();
    m.mA1 = -m.mA0;
    for (double a : theta) {
        if (a < cut)
            a += 2.0 * kKfPi;
        m.mA0 = std::min(m.mA0, a);
        m.mA1 = std::max(m.mA1, a);
    }
    return m;
}

// The line-complex axis is biased by the skin's discretised normals (a
// thousandth of a radian on a coarse skin, enough to smear the meridian
// profile past the tolerance), so it is refined against the profile: the
// profile fitted, the axis (two tilts and two offsets) moved by
// Levenberg-Marquardt to bring a fixed subsample of the points onto it,
// twice. Fills `rF.mProfile` and `rMer`; false when no profile fits.
bool kf_refine_axis(const std::vector<PatchVec>& rPts, const KfPatch& rP, double Tol,
                    KinematicFit& rF, KfMeridian& rMer) {
    const std::size_t stride = std::max<std::size_t>(1, rP.mLive.size() / 400);
    std::vector<std::size_t> sub;
    for (std::size_t i = 0; i < rP.mLive.size(); i += stride)
        sub.push_back(rP.mLive[i]);
    for (int round = 0; round < 2; ++round) {
        if (!kf_profile(rMer.mRz, round == 0 ? 4.0 * Tol : Tol, rF.mProfile, rF.mWhy))
            return false;
        const KinematicFit base = rF;
        const auto moved = [&](const std::vector<double>& rTheta, KinematicFit& rOut) {
            rOut.mOrigin = base.mOrigin;
            PatchVec a = base.mAxis;
            for (std::size_t k = 0; k < 3; ++k) {
                a[k] += rTheta[0] * base.mE1[k] + rTheta[1] * base.mE2[k];
                rOut.mOrigin[k] += rTheta[2] * base.mE1[k] + rTheta[3] * base.mE2[k];
            }
            const double l = vec3_norm(a);
            for (double& c : a)
                c /= l;
            rOut.mAxis = a;
            lsq_frame(rOut.mAxis.data(), rOut.mE1.data(), rOut.mE2.data());
        };
        std::vector<double> theta(4, 0.0);
        const std::vector<double> typical{1e-3, 1e-3, 1e-3 * rP.mScale, 1e-3 * rP.mScale};
        lsq_levenberg_marquardt(
            theta, typical, sub.size(),
            [&](const std::vector<double>& rTheta, double* pOut) {
                KinematicFit g = base;
                moved(rTheta, g);
                for (std::size_t i = 0; i < sub.size(); ++i) {
                    const PatchVec r = vec3_sub(rPts[sub[i]], g.mOrigin);
                    const double z = vec3_dot(r, g.mAxis);
                    const double rho = std::hypot(vec3_dot(r, g.mE1), vec3_dot(r, g.mE2));
                    double t = 0.0;
                    const double d = bspline_curve_project(base.mProfile, {rho, z, 0.0}, t, 16);
                    PatchVec d1, d2;
                    const PatchVec c = bspline_curve_eval2(base.mProfile, t, d1, d2);
                    const double side = d1[0] * (z - c[1]) - d1[1] * (rho - c[0]);
                    pOut[i] = side < 0.0 ? -d : d;
                }
            },
            10);
        moved(theta, rF);
        rMer = kf_meridian(rPts, rP.mLive, rF);
    }
    return kf_profile(rMer.mRz, Tol, rF.mProfile, rF.mWhy);
}

// The parameters of every used point on the surface: the nearest node of a
// grid, then Newton.
void kf_params(const BsplineSurface& rS, const std::vector<PatchVec>& rPts, const KfPatch& rP,
               std::vector<double>& rU, std::vector<double>& rV) {
    constexpr int kG = 24;
    std::vector<PatchVec> node;
    std::vector<KfXy> uv;
    for (int j = 0; j <= kG; ++j)
        for (int i = 0; i <= kG; ++i) {
            const double u = rS.U0() + (rS.U1() - rS.U0()) * i / kG;
            const double v = rS.V0() + (rS.V1() - rS.V0()) * j / kG;
            BsplineEval e;
            bspline_eval(rS, u, v, 0, e);
            node.push_back({e.mS[0], e.mS[1], e.mS[2]});
            uv.push_back({u, v});
        }
    rU.assign(rPts.size(), 0.0);
    rV.assign(rPts.size(), 0.0);
    for (std::size_t v : rP.mLive) {
        std::size_t best = 0;
        double bd = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0; k < node.size(); ++k) {
            const double d = vec3_norm_sq(vec3_sub(node[k], rPts[v]));
            if (d < bd) {
                bd = d;
                best = k;
            }
        }
        rU[v] = uv[best][0];
        rV[v] = uv[best][1];
        bspline_project(rS, rPts[v].data(), rU[v], rV[v]);
    }
}

// The parameters of every used point on a helical fit's surface, from the
// screw itself: the angle unwrapped along the patch's edges (a thread winds
// round many times, where the nearest node of a grid can lie a turn away),
// the profile parameter from the point's radius or height; then Newton.
void kf_helical_params(const KinematicFit& rF, const std::vector<PatchVec>& rPts,
                       const std::vector<PatchTri>& rTris, const KfPatch& rP,
                       std::vector<double>& rU, std::vector<double>& rV) {
    const HelicalFit& h = rF.mHelical;
    const std::size_t n = rPts.size();
    std::vector<std::vector<std::size_t>> adj(n);
    for (const PatchTri& t : rTris)
        for (std::size_t k = 0; k < 3; ++k) {
            adj[t[k]].push_back(t[(k + 1) % 3]);
            adj[t[(k + 1) % 3]].push_back(t[k]);
        }
    for (auto& a : adj) {
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
    }
    std::vector<double> raw(n, 0.0), theta(n, 0.0);
    for (std::size_t v : rP.mLive) {
        const PatchVec r = vec3_sub(rPts[v], h.mOrigin);
        raw[v] = std::atan2(vec3_dot(r, h.mE2), vec3_dot(r, h.mE1));
    }
    std::vector<std::uint8_t> seen(n, 0);
    double lo = std::numeric_limits<double>::infinity();
    for (std::size_t seed : rP.mLive) {
        if (seen[seed])
            continue;
        std::vector<std::size_t> queue{seed};
        seen[seed] = 1;
        theta[seed] = raw[seed];
        for (std::size_t q = 0; q < queue.size(); ++q) {
            const std::size_t a = queue[q];
            lo = std::min(lo, theta[a]);
            for (std::size_t b : adj[a])
                if (!seen[b]) {
                    seen[b] = 1;
                    theta[b] = theta[a] + std::remainder(raw[b] - theta[a], 2.0 * kKfPi);
                    queue.push_back(b);
                }
        }
    }
    const double shift = 2.0 * kKfPi * std::round((h.mThetaMin - lo) / (2.0 * kKfPi));
    rU.assign(n, 0.0);
    rV.assign(n, 0.0);
    for (std::size_t v : rP.mLive) {
        const double th = theta[v] + shift;
        const PatchVec r = vec3_sub(rPts[v], h.mOrigin);
        const double g = h.Advance(th);
        const double rho = std::hypot(vec3_dot(r, h.mE1), vec3_dot(r, h.mE2)) - h.mTaper * g;
        const double zeta = vec3_dot(r, h.mAxis) - g;
        rU[v] = std::clamp(h.mRhoOfZeta ? zeta : rho, rF.mSurface.U0(), rF.mSurface.U1());
        rV[v] = std::clamp(th, rF.mSurface.V0(), rF.mSurface.V1());
        bspline_project(rF.mSurface, rPts[v].data(), rU[v], rV[v]);
    }
}

// A cubic B-spline over [U0, U1] x [V0, V1] through the exact surface
// `rAt`, within a twentieth of the tolerance where it can be.
template <class F>
BsplineSurface kf_approximate(F&& rAt, double U0, double U1, double V0, double V1, std::size_t Nu,
                              std::size_t Nv, double Tol) {
    BsplineSurface best;
    for (int attempt = 0; attempt < 3; ++attempt) {
        const std::size_t nu = std::min<std::size_t>(Nu << attempt, 256);
        const std::size_t nv = std::min<std::size_t>(Nv << attempt, 512);
        BsplineSurface s = bspline_shape(3, 3, nu, nv, U0, U1, V0, V1);
        const std::size_t mu = 4 * nu + 1, mv = 4 * nv + 1;
        std::vector<double> us, vs, xyz;
        for (std::size_t j = 0; j < mv; ++j)
            for (std::size_t i = 0; i < mu; ++i) {
                const double u =
                    U0 + (U1 - U0) * static_cast<double>(i) / static_cast<double>(mu - 1);
                const double v =
                    V0 + (V1 - V0) * static_cast<double>(j) / static_cast<double>(mv - 1);
                const PatchVec p = rAt(u, v);
                us.push_back(u);
                vs.push_back(v);
                xyz.insert(xyz.end(), p.begin(), p.end());
            }
        if (!bspline_lsq_fit(s, us, vs, xyz, nullptr, 0.0))
            break;
        best = s;
        // Checked half-way between the samples.
        double worst = 0.0;
        for (std::size_t j = 0; j + 1 < mv; ++j)
            for (std::size_t i = 0; i + 1 < mu; ++i) {
                const double u =
                    U0 + (U1 - U0) * (static_cast<double>(i) + 0.5) / static_cast<double>(mu - 1);
                const double v =
                    V0 + (V1 - V0) * (static_cast<double>(j) + 0.5) / static_cast<double>(mv - 1);
                BsplineEval e;
                bspline_eval(s, u, v, 0, e);
                const PatchVec p = rAt(u, v);
                worst =
                    std::max(worst, vec3_norm({e.mS[0] - p[0], e.mS[1] - p[1], e.mS[2] - p[2]}));
            }
        if (worst <= 0.05 * Tol)
            break;
    }
    return best;
}

// The acceptance of a fit: `freeform_accept` on its surface.
bool kf_accept(KinematicFit& rFit, const std::vector<PatchVec>& rPts,
               const std::vector<PatchTri>& rTris, const KfPatch& rP, double Tol) {
    if (rFit.mSurface.NumPoles() == 0) {
        rFit.mWhy = "no surface";
        return false;
    }
    if (rFit.mKind == KinematicKind::Helical)
        kf_helical_params(rFit, rPts, rTris, rP, rFit.mU, rFit.mV);
    else
        kf_params(rFit.mSurface, rPts, rP, rFit.mU, rFit.mV);
    const FreeformAcceptance acc =
        freeform_accept(rFit.mSurface, rPts, rTris, rFit.mU, rFit.mV, Tol, 0.0, true);
    rFit.mMaxDist = acc.mMaxDist;
    rFit.mRmsDist = acc.mRmsDist;
    rFit.mMaxCentroidDist = acc.mMaxCentroidDist;
    rFit.mReverseMaxDist = acc.mReverseMaxDist;
    if (!acc.mOk)
        rFit.mWhy = "the surface misses the patch (" + std::to_string(acc.mMaxDist) + ", " +
                    std::to_string(acc.mMaxCentroidDist) + ", " +
                    std::to_string(acc.mReverseMaxDist) + ")";
    return acc.mOk;
}

}  // namespace

PatchVec KinematicFit::ProfilePoint(double T) const {
    const PatchVec q = bspline_curve_eval(mProfile, T);
    const PatchVec& second = mKind == KinematicKind::Revolution ? mAxis : mE2;
    return {mOrigin[0] + q[0] * mE1[0] + q[1] * second[0],
            mOrigin[1] + q[0] * mE1[1] + q[1] * second[1],
            mOrigin[2] + q[0] * mE1[2] + q[1] * second[2]};
}

BsplineCurve KinematicFit::SweptCurve() const {
    BsplineCurve c = mProfile;
    const PatchVec& second = mKind == KinematicKind::Revolution ? mAxis : mE2;
    for (PatchVec& q : c.mPoles) {
        const PatchVec p = q;
        for (std::size_t k = 0; k < 3; ++k)
            q[k] = mOrigin[k] + p[0] * mE1[k] + p[1] * second[k];
    }
    return c;
}

PatchVec KinematicFit::At(double T, double S) const {
    switch (mKind) {
        case KinematicKind::Extrusion: {
            const PatchVec p = ProfilePoint(T);
            return {p[0] + S * mAxis[0], p[1] + S * mAxis[1], p[2] + S * mAxis[2]};
        }
        case KinematicKind::Revolution: {
            const PatchVec q = bspline_curve_eval(mProfile, T);
            const double c = std::cos(S), s = std::sin(S);
            PatchVec out{};
            for (std::size_t k = 0; k < 3; ++k)
                out[k] = mOrigin[k] + q[1] * mAxis[k] + q[0] * (c * mE1[k] + s * mE2[k]);
            return out;
        }
        case KinematicKind::Helical:
            return mHelical.At(T, S);
        default:
            return mOrigin;
    }
}

bool kinematic_is_band(const std::vector<PatchTri>& rTris) {
    const PatchTopology t = patch_topology(rTris);
    return t.mNumBoundaryLoops == 2 && t.mEuler == 0;
}

KinematicFit kinematic_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                           const KinematicOptions& rOptions) {
    KinematicFit out;
    const double tol = rOptions.mTolerance;
    if (!(tol > 0.0) || rTris.size() < rOptions.mMinTriangles) {
        out.mWhy = "too small";
        return out;
    }
    if (kinematic_is_band(rTris) || patch_topology(rTris).mNumBoundaryLoops != 1) {
        out.mWhy = "not a disc";
        return out;
    }
    const KfPatch patch = kf_patch(rPts, rTris);
    if (patch.mLive.size() < 16 || !(patch.mScale > 0.0)) {
        out.mWhy = "degenerate";
        return out;
    }
    std::string why;

    // 1. Extrusion: every normal perpendicular to one direction.
    if (rOptions.mExtrusion) {
        double C[9] = {0.0}, wsum = 0.0;
        for (std::size_t v : patch.mLive) {
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    C[i * 3 + j] += patch.mArea[v] *
                                    patch.mNormals[v][static_cast<std::size_t>(i)] *
                                    patch.mNormals[v][static_cast<std::size_t>(j)];
            wsum += patch.mArea[v];
        }
        for (double& c : C)
            c /= wsum;
        double w[3], vec[9];
        lsq_sym3_eigen(C, w, vec);
        if (w[0] <= 3e-3 && w[1] >= 1e-2) {
            KinematicFit f;
            f.mKind = KinematicKind::Extrusion;
            f.mAxis = {vec[0], vec[1], vec[2]};
            lsq_frame(f.mAxis.data(), f.mE1.data(), f.mE2.data());
            f.mOrigin = patch.mCentre;
            std::vector<KfXy> xy;
            double z0 = std::numeric_limits<double>::infinity(), z1 = -z0;
            for (std::size_t v : patch.mLive) {
                const PatchVec r = vec3_sub(rPts[v], f.mOrigin);
                xy.push_back({vec3_dot(r, f.mE1), vec3_dot(r, f.mE2)});
                const double z = vec3_dot(r, f.mAxis);
                z0 = std::min(z0, z);
                z1 = std::max(z1, z);
            }
            if (kf_profile(xy, tol, f.mProfile, f.mWhy) && z1 > z0) {
                const double m = rOptions.mMargin * (z1 - z0);
                f.mV0 = z0 - m;
                f.mV1 = z1 + m;
                // Exactly: the profile's poles at both heights, degree 1 along.
                BsplineSurface& s = f.mSurface;
                s.mDegreeU = f.mProfile.mDegree;
                s.mDegreeV = 1;
                s.mKnotsU = f.mProfile.mKnots;
                s.mKnotsV = {f.mV0, f.mV0, f.mV1, f.mV1};
                s.mNu = f.mProfile.mPoles.size();
                s.mNv = 2;
                s.mPoles.assign(s.mNu * 2 * 3, 0.0);
                for (std::size_t j = 0; j < 2; ++j)
                    for (std::size_t i = 0; i < s.mNu; ++i) {
                        const PatchVec& q = f.mProfile.mPoles[i];
                        const double h = j == 0 ? f.mV0 : f.mV1;
                        for (std::size_t k = 0; k < 3; ++k)
                            s.mPoles[3 * (i + s.mNu * j) + k] =
                                f.mOrigin[k] + q[0] * f.mE1[k] + q[1] * f.mE2[k] + h * f.mAxis[k];
                    }
                if (kf_accept(f, rPts, rTris, patch, tol)) {
                    f.mOk = true;
                    return f;
                }
            }
            why += "extrusion: " + f.mWhy + "; ";
        } else {
            why += "extrusion: the normals turn about no one direction; ";
        }
    }

    // 2. Revolution: a screw with no advance.
    if (rOptions.mRevolution) {
        const ScrewEstimate est = screw_estimate(rPts, patch.mNormals, patch.mArea, patch.mLive,
                                                 patch.mCentre, patch.mScale);
        // The advance the estimate finds over a turn is the discretisation's
        // too (a few thousandths of the size on a coarse skin), so a screw
        // is tried as a revolution when its advance over the patch's angles
        // stays within the tolerance; the acceptance has the last word.
        if (est.mOk && std::abs(est.mK) * 2.0 * kKfPi <= 0.05 * patch.mScale) {
            KinematicFit f;
            f.mKind = KinematicKind::Revolution;
            f.mAxis = est.mAxis;
            // The origin at the foot of the centre on the axis.
            const PatchVec d = vec3_sub(patch.mCentre, est.mOrigin);
            const double along = vec3_dot(d, f.mAxis);
            for (std::size_t k = 0; k < 3; ++k)
                f.mOrigin[k] = est.mOrigin[k] + along * f.mAxis[k];
            lsq_frame(f.mAxis.data(), f.mE1.data(), f.mE2.data());
            KfMeridian mer = kf_meridian(rPts, patch.mLive, f);
            if (std::abs(est.mK) * (mer.mA1 - mer.mA0) > tol) {
                f.mWhy = "the patch advances along the axis";
            } else if (!(mer.mRhoMin > 1e-3 * patch.mScale)) {
                f.mWhy = "the patch reaches the axis";
            } else if (mer.mGap < kKfPi / 12.0) {
                f.mWhy = "the patch winds round the axis";
            } else if (kf_refine_axis(rPts, patch, tol, f, mer)) {
                const double a0 = mer.mA0, a1 = mer.mA1, gap = mer.mGap;
                const double m = std::min(rOptions.mMargin * (a1 - a0), 0.5 * (gap - 1e-3));
                // The frame turned so the angles run about pi, inside
                // (0, 2 pi): a CAD kernel's periodic angle then needs no
                // shift past its seam.
                const double turn = 0.5 * (a0 + a1) - kKfPi;
                const double ct = std::cos(turn), st = std::sin(turn);
                const PatchVec e1 = f.mE1, e2 = f.mE2;
                for (std::size_t k = 0; k < 3; ++k) {
                    f.mE1[k] = ct * e1[k] + st * e2[k];
                    f.mE2[k] = -st * e1[k] + ct * e2[k];
                }
                f.mV0 = a0 - turn - m;
                f.mV1 = a1 - turn + m;
                const std::size_t nu = f.mProfile.mPoles.size() + 2;
                const std::size_t nv = std::max<std::size_t>(
                    4, static_cast<std::size_t>(std::ceil((f.mV1 - f.mV0) / (kKfPi / 6.0))) + 3);
                f.mSurface =
                    kf_approximate([&](double U, double V) { return f.At(U, V); }, f.mProfile.T0(),
                                   f.mProfile.T1(), f.mV0, f.mV1, nu, nv, tol);
                if (kf_accept(f, rPts, rTris, patch, tol)) {
                    f.mOk = true;
                    return f;
                }
            }
            why += "revolution: " + f.mWhy + "; ";
        } else {
            why += "revolution: the normals meet no one axis; ";
        }
    }

    // 3. Helical: an advancing screw.
    if (rOptions.mHelical) {
        KinematicFit f;
        f.mKind = KinematicKind::Helical;
        f.mHelical = helical_fit(rPts, rTris, nullptr, nullptr, tol, true);
        if (f.mHelical.mOk) {
            const HelicalFit& h = f.mHelical;
            const double ds = rOptions.mMargin * (h.mS1 - h.mS0);
            const double dt = 0.02 * (h.mThetaMax - h.mThetaMin);
            f.mOrigin = h.mOrigin;
            f.mAxis = h.mAxis;
            f.mE1 = h.mE1;
            f.mE2 = h.mE2;
            f.mV0 = h.mThetaMin - dt;
            f.mV1 = h.mThetaMax + dt;
            const std::size_t nv = std::min<std::size_t>(
                400, static_cast<std::size_t>(std::ceil((f.mV1 - f.mV0) / (kKfPi / 6.0))) + 3);
            f.mSurface = kf_approximate([&](double U, double V) { return h.At(U, V); }, h.mS0 - ds,
                                        h.mS1 + ds, f.mV0, f.mV1, h.mStraight ? 4 : 8, nv, tol);
            if (kf_accept(f, rPts, rTris, patch, tol)) {
                f.mOk = true;
                return f;
            }
            why += "helical: " + f.mWhy;
        } else {
            why += "helical: no screw fits";
        }
    }
    out.mWhy = why;
    return out;
}

}  // namespace meshioplusplus::detail
