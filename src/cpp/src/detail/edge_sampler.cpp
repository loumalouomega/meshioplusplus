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
// Surface and trim evaluation, arc length, sampling and projection of
// detail/edge_sampler.hpp. Anonymous-namespace helpers are prefixed `edgs_`
// (the amalgamation concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "edge_sampler.hpp"
#include "nurbs_eval.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kEdgsGaussX[5] = {
    -0.906179845938663992797626878299, -0.538469310105683091036314420700, 0.0,
    0.538469310105683091036314420700, 0.906179845938663992797626878299};
constexpr double kEdgsGaussW[5] = {
    0.236926885056189087514264040720, 0.478628670499366468041291514836,
    0.568888888888888888888888888889, 0.478628670499366468041291514836,
    0.236926885056189087514264040720};

[[noreturn]] void edgs_bad(const std::string& rWhat) {
    throw std::invalid_argument("meshio++: mesh_cad_surface: " + rWhat);
}

double edgs_dist(const EdgeVec3& rA, const EdgeVec3& rB) {
    return std::sqrt((rA[0] - rB[0]) * (rA[0] - rB[0]) + (rA[1] - rB[1]) * (rA[1] - rB[1]) +
                     (rA[2] - rB[2]) * (rA[2] - rB[2]));
}

/// The distance from `rP` to the segment `rA rB`.
double edgs_segment_distance(const EdgeVec3& rP, const EdgeVec3& rA, const EdgeVec3& rB) {
    EdgeVec3 d{rB[0] - rA[0], rB[1] - rA[1], rB[2] - rA[2]};
    const double len2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    double t = 0.0;
    if (len2 > 0.0)
        t = std::clamp(
            ((rP[0] - rA[0]) * d[0] + (rP[1] - rA[1]) * d[1] + (rP[2] - rA[2]) * d[2]) / len2, 0.0,
            1.0);
    const EdgeVec3 q{rA[0] + t * d[0], rA[1] + t * d[1], rA[2] + t * d[2]};
    return edgs_dist(rP, q);
}

}  // namespace

BrepSurfaceEval::BrepSurfaceEval(const NurbsPatch& rPatch, bool Rational) : pPatch(&rPatch) {
    if (rPatch.mDegrees.size() != 2 || rPatch.mKnots.size() != 2)
        edgs_bad("a face's surface is not a two-parameter patch");
    std::size_t total = 1;
    for (std::size_t d = 0; d < 2; ++d) {
        const int p = rPatch.mDegrees[d];
        const std::size_t nk = rPatch.mKnots[d].size();
        if (p < 1 || p > 32 || nk < 2 * static_cast<std::size_t>(p) + 2)
            edgs_bad("a face's surface has an invalid degree or knot vector");
        for (std::size_t i = 1; i < nk; ++i)
            if (!(rPatch.mKnots[d][i] >= rPatch.mKnots[d][i - 1]))
                edgs_bad("a face's surface has a decreasing knot vector");
        mCounts[d] = nk - static_cast<std::size_t>(p) - 1;
        total *= mCounts[d];
        mBox[2 * d] = rPatch.mKnots[d][static_cast<std::size_t>(p)];
        mBox[2 * d + 1] = rPatch.mKnots[d][mCounts[d]];
        if (!(mBox[2 * d + 1] > mBox[2 * d]))
            edgs_bad("a face's surface has an empty parameter range");
    }
    const NDArray& r_cp = rPatch.mControlPoints;
    if (r_cp.Dtype() != DType::Float64 || r_cp.Shape().size() != 2 || r_cp.Shape()[0] != total ||
        r_cp.Shape()[1] != 3)
        edgs_bad("a face's surface control points are not Float64 (n, 3)");
    const bool weighted =
        Rational && rPatch.mWeights.Size() == total && rPatch.mWeights.Dtype() == DType::Float64;
    mNet.mDim = 2;
    for (std::size_t d = 0; d < 2; ++d) {
        mNet.mDegrees[d] = rPatch.mDegrees[d];
        mNet.mKnots[d] = &rPatch.mKnots[d];
        mNet.mCounts[d] = mCounts[d];
    }
    mNet.mValues = r_cp.As<double>();
    mNet.mColumns = 3;
    mNet.mWeights = weighted ? rPatch.mWeights.As<double>() : nullptr;
    for (std::size_t i = 0; i < 3 * total; ++i)
        if (std::isfinite(mNet.mValues[i]))
            mScale = std::max(mScale, std::fabs(mNet.mValues[i]));
}

void BrepSurfaceEval::Eval(double U, double V, double* pX, double* pSu, double* pSv) const {
    const double u[2] = {std::clamp(U, mBox[0], mBox[1]), std::clamp(V, mBox[2], mBox[3])};
    std::size_t span[2];
    for (std::size_t d = 0; d < 2; ++d)
        span[d] = nurbs_find_span(pPatch->mDegrees[d], pPatch->mKnots[d], mCounts[d], u[d]);
    double der[6];
    const bool want = pSu != nullptr || pSv != nullptr;
    iga_eval_point(mNet, span, u, want ? 1 : 0, pX, want ? der : nullptr);
    if (pSu)
        for (int k = 0; k < 3; ++k)
            pSu[k] = der[k];
    if (pSv)
        for (int k = 0; k < 3; ++k)
            pSv[k] = der[3 + k];
}

void BrepSurfaceEval::Eval2(double U, double V, double* pX, double* pSu, double* pSv, double* pSuu,
                            double* pSuv, double* pSvv) const {
    const double u[2] = {std::clamp(U, mBox[0], mBox[1]), std::clamp(V, mBox[2], mBox[3])};
    std::size_t span[2];
    for (std::size_t d = 0; d < 2; ++d)
        span[d] = nurbs_find_span(pPatch->mDegrees[d], pPatch->mKnots[d], mCounts[d], u[d]);
    double d1[6], d2[9];
    iga_eval_surface_d2(mNet, span, u, pX, d1, d2);
    for (int k = 0; k < 3; ++k) {
        pSu[k] = d1[k];
        pSv[k] = d1[3 + k];
        pSuu[k] = d2[k];
        pSuv[k] = d2[3 + k];
        pSvv[k] = d2[6 + k];
    }
}

double BrepSurfaceEval::MaxCurvature(double U, double V) const {
    double x[3], su[3], sv[3], suu[3], suv[3], svv[3];
    Eval2(U, V, x, su, sv, suu, suv, svv);
    const auto dot = [](const double* pA, const double* pB) {
        return pA[0] * pB[0] + pA[1] * pB[1] + pA[2] * pB[2];
    };
    double n[3] = {su[1] * sv[2] - su[2] * sv[1], su[2] * sv[0] - su[0] * sv[2],
                   su[0] * sv[1] - su[1] * sv[0]};
    const double len = std::sqrt(dot(n, n));
    const double e = dot(su, su), f = dot(su, sv), g = dot(sv, sv);
    const double det = e * g - f * f;
    if (!(len > 0.0) || !(det > 1e-14 * e * g))
        return 0.0;
    for (double& r_c : n)
        r_c /= len;
    const double l = dot(suu, n), m = dot(suv, n), nn = dot(svv, n);
    // The shape operator I^-1 II: mean and Gaussian curvature, then the
    // principal curvatures H +- sqrt(H^2 - K).
    const double mean = (e * nn - 2.0 * f * m + g * l) / (2.0 * det);
    const double gauss = (l * nn - m * m) / det;
    return std::abs(mean) + std::sqrt(std::max(0.0, mean * mean - gauss));
}

BrepCurveEval::BrepCurveEval(const BrepCurve& rCurve) : mDegree(rCurve.mDegree) {
    mCount = rCurve.mControlPoints.Shape().empty() ? 0 : rCurve.mControlPoints.Shape()[0];
    const auto p = static_cast<std::size_t>(mDegree);
    if (mDegree < 1 || mDegree > 32 || rCurve.mControlPoints.Dtype() != DType::Float64 ||
        rCurve.mControlPoints.Shape().size() != 2 || rCurve.mControlPoints.Shape()[1] != 3 ||
        mCount < p + 1 || rCurve.mKnots.size() != mCount + p + 1)
        edgs_bad("a trimming curve has an invalid degree, control points or knot vector");
    mKnots = rCurve.mKnots;
    for (std::size_t i = 1; i < mKnots.size(); ++i)
        if (!(mKnots[i] >= mKnots[i - 1]))
            edgs_bad("a trimming curve has a decreasing knot vector");
    const double* cp = rCurve.mControlPoints.As<double>();
    mValues.assign(cp, cp + 3 * mCount);
    if (rCurve.mRational && rCurve.mWeights.Size() == mCount &&
        rCurve.mWeights.Dtype() == DType::Float64) {
        const double* w = rCurve.mWeights.As<double>();
        mWeights.assign(w, w + mCount);
    }
    const double k0 = mKnots[p];
    const double k1 = mKnots[mCount];
    mT0 = k0;
    mT1 = k1;
    if (rCurve.mActiveRange.size() == 2) {
        mT0 = std::max(k0, std::min(rCurve.mActiveRange[0], rCurve.mActiveRange[1]));
        mT1 = std::min(k1, std::max(rCurve.mActiveRange[0], rCurve.mActiveRange[1]));
    }
    if (!(mT1 > mT0))
        edgs_bad("a trimming curve has an empty parameter range");
    for (std::size_t i = 0; i < mKnots.size(); ++i) {
        const double k = mKnots[i];
        if (!(k > mT0 && k < mT1) || (!mBreaks.empty() && mBreaks.back() == k))
            continue;
        mBreaks.push_back(k);
        std::size_t mult = 0;
        for (std::size_t j = i; j < mKnots.size() && mKnots[j] == k; ++j)
            ++mult;
        if (mult >= p)
            mCorners.push_back(k);
    }
    mNet.mDim = 1;
    mNet.mDegrees[0] = mDegree;
    mNet.mKnots[0] = &mKnots;
    mNet.mCounts[0] = mCount;
    mNet.mValues = mValues.data();
    mNet.mColumns = 3;
    mNet.mWeights = mWeights.empty() ? nullptr : mWeights.data();
}

void BrepCurveEval::Eval(double T, double* pX, double* pD) const {
    // The net points into this object's own vectors, which a copy or move
    // does not carry: rebuild the pointers on the stack.
    IgaEvalNet net = mNet;
    net.mKnots[0] = &mKnots;
    net.mValues = mValues.data();
    net.mWeights = mWeights.empty() ? nullptr : mWeights.data();
    const double t = std::clamp(T, mT0, mT1);
    const std::size_t span = nurbs_find_span(mDegree, mKnots, mCount, t);
    double der[3];
    iga_eval_point(net, &span, &t, pD ? 1 : 0, pX, pD ? der : nullptr);
    if (pD)
        for (int k = 0; k < 3; ++k)
            pD[k] = der[k];
}

TrimTrace::TrimTrace(const BrepSurfaceEval* pSurface, BrepCurveEval Curve, bool Forward)
    : pSurface(pSurface), mCurve(std::move(Curve)), mForward(Forward) {
    std::vector<double> breaks{mCurve.T0()};
    breaks.insert(breaks.end(), mCurve.Breaks().begin(), mCurve.Breaks().end());
    breaks.push_back(mCurve.T1());
    if (!mForward)
        std::reverse(breaks.begin(), breaks.end());
    constexpr int kStart = 4;
    // The speed carries round-off relative to the coordinates, not to the
    // trace: on a short trace far from the origin, or one whose image is a
    // point (a pole's trim), no relative tolerance is met and the pieces
    // halve to the depth limit. A piece is accurate enough within 1e-14 of
    // the whole trace's length, or shorter than round-off at the surface's
    // size.
    double estimate = 0.0;
    for (std::size_t i = 0; i + 1 < breaks.size(); ++i)
        for (int j = 0; j < kStart; ++j) {
            const double sa = breaks[i] + (breaks[i + 1] - breaks[i]) * (double(j) / kStart);
            const double sb = j + 1 == kStart ? breaks[i + 1]
                                              : breaks[i] + (breaks[i + 1] - breaks[i]) *
                                                                (double(j + 1) / kStart);
            estimate += Gauss(sa, sb);
        }
    mNegligible = 1e-15 * pSurface->Scale();
    mAbsTolerance = 1e-14 * estimate;
    mT.push_back(breaks.front());
    mCum.push_back(0.0);
    for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
        const double a = breaks[i];
        const double b = breaks[i + 1];
        for (int j = 0; j < kStart; ++j) {
            const double sa = a + (b - a) * (double(j) / kStart);
            const double sb = j + 1 == kStart ? b : a + (b - a) * (double(j + 1) / kStart);
            Refine(sa, sb, Gauss(sa, sb), mCum.back(), 0, mT, mCum);
        }
    }
}

EdgeVec2 TrimTrace::Uv(double T) const {
    double x[3];
    mCurve.Eval(T, x, nullptr);
    return {x[0], x[1]};
}

EdgeVec3 TrimTrace::Point(double T) const {
    const EdgeVec2 uv = Uv(T);
    EdgeVec3 p{};
    pSurface->Eval(uv[0], uv[1], p.data(), nullptr, nullptr);
    return p;
}

EdgeVec3 TrimTrace::Velocity(double T) const {
    double c[3];
    double dc[3];
    mCurve.Eval(T, c, dc);
    double x[3];
    double su[3];
    double sv[3];
    pSurface->Eval(c[0], c[1], x, su, sv);
    return {su[0] * dc[0] + sv[0] * dc[1], su[1] * dc[0] + sv[1] * dc[1],
            su[2] * dc[0] + sv[2] * dc[1]};
}

double TrimTrace::Speed(double T) const {
    const EdgeVec3 v = Velocity(T);
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

double TrimTrace::Gauss(double A, double B) const {
    const double h = 0.5 * (B - A);
    const double c = 0.5 * (A + B);
    double s = 0.0;
    for (int i = 0; i < 5; ++i)
        s += kEdgsGaussW[i] * Speed(c + h * kEdgsGaussX[i]);
    return s * std::fabs(h);
}

void TrimTrace::Refine(double A, double B, double Whole, double La, int Depth,
                       std::vector<double>& rT, std::vector<double>& rL) const {
    const double m = 0.5 * (A + B);
    const double left = Gauss(A, m);
    const double right = Gauss(m, B);
    const double halves = left + right;
    // And never past a table of 2^15 entries, where heavier round-off (a
    // short trace on a patch far from the origin) would still halve on:
    // depth first, so the cap falls the same way every time.
    if (Depth < 30 && std::fabs(Whole - halves) > std::max(1e-12 * halves, mAbsTolerance) &&
        halves > mNegligible && rT.size() < (std::size_t(1) << 15) && m != A && m != B) {
        Refine(A, m, left, La, Depth + 1, rT, rL);
        Refine(m, B, right, rL.back(), Depth + 1, rT, rL);
        return;
    }
    rT.push_back(m);
    rL.push_back(La + left);
    rT.push_back(B);
    rL.push_back(La + halves);
}

double TrimTrace::ParamAtLength(double L) const {
    const double total = Length();
    if (!(total > 0.0) || L <= 0.0)
        return mT.front();
    if (L >= total)
        return mT.back();
    const auto it = std::upper_bound(mCum.begin(), mCum.end(), L);
    const auto j = static_cast<std::size_t>(std::max<std::ptrdiff_t>(1, it - mCum.begin())) - 1;
    const double ta = mT[j];
    const double tb = mT[j + 1];
    const double ca = mCum[j];
    const double cb = mCum[j + 1];
    if (!(cb > ca))
        return ta;
    double s = (L - ca) / (cb - ca);
    double lo = 0.0;
    double hi = 1.0;
    for (int it2 = 0; it2 < 40; ++it2) {
        const double t = ta + s * (tb - ta);
        const double g = ca + Gauss(ta, t) - L;
        if (std::fabs(g) <= 1e-14 * total)
            break;
        if (g > 0.0)
            hi = s;
        else
            lo = s;
        const double dg = Speed(t) * std::fabs(tb - ta);
        double next = dg > 0.0 ? s - g / dg : 0.5 * (lo + hi);
        if (!(next > lo && next < hi))
            next = 0.5 * (lo + hi);
        if (next == s)
            break;
        s = next;
    }
    return ta + s * (tb - ta);
}

double TrimTrace::LengthAt(double T) const {
    if (mT.size() < 2)
        return 0.0;
    const bool up = mT.back() > mT.front();
    if (up ? T <= mT.front() : T >= mT.front())
        return 0.0;
    if (up ? T >= mT.back() : T <= mT.back())
        return Length();
    const auto it = up ? std::upper_bound(mT.begin(), mT.end(), T)
                       : std::upper_bound(mT.begin(), mT.end(), T, std::greater<double>());
    const auto j = static_cast<std::size_t>(std::max<std::ptrdiff_t>(1, it - mT.begin())) - 1;
    return mCum[j] + Gauss(mT[j], T);
}

std::vector<double> TrimTrace::CornerLengths() const {
    // Every knot is a node of the table, so a corner's arc length is exact.
    std::vector<double> out;
    for (const double c : mCurve.Corners())
        for (std::size_t i = 0; i < mT.size(); ++i)
            if (mT[i] == c) {
                out.push_back(mCum[i]);
                break;
            }
    std::sort(out.begin(), out.end());
    return out;
}

namespace {

/// One stretch `[L0, L1]` of a trace cut at equal increments of the integral
/// of ds / h(s) (roadmap 6.3.4): the density 1 / h is integrated by adaptive
/// Simpson from pieces no longer than the target, the stretch gets
/// `ceil(N)` pieces (and the minimum the options ask for), and each piece's
/// end is found by interpolating the cumulative integral, then handed to
/// `Split` for the chordal cap.
template <class TSplit>
void edgs_sized_stretch(const TrimTrace& rTrace, const EdgeSampleOptions& rOptions, double L0,
                        double L1, double Length, std::int64_t& rBudget, std::vector<double>& rEnds,
                        TSplit&& rSplit) {
    const double part = L1 - L0;
    const auto density = [&](double L) {
        const double hx = rOptions.mSize(rTrace.Point(rTrace.ParamAtLength(L)));
        return 1.0 / std::min(rOptions.mTargetSize, std::max(hx, 1e-12 * rOptions.mTargetSize));
    };
    std::vector<double> at{L0}, cum{0.0};
    const auto simpson = [](double A, double B, double Fa, double Fm, double Fb) {
        return (B - A) / 6.0 * (Fa + 4.0 * Fm + Fb);
    };
    const auto adapt = [&](auto&& rSelf, double A, double B, double Fa, double Fm, double Fb,
                           double Whole, int Depth) -> void {
        const double m = 0.5 * (A + B);
        const double lm = 0.5 * (A + m), rm = 0.5 * (m + B);
        const double flm = density(lm), frm = density(rm);
        const double left = simpson(A, m, Fa, flm, Fm);
        const double right = simpson(m, B, Fm, frm, Fb);
        if (Depth >= 24 || std::abs(left + right - Whole) <= 1e-3 * (left + right)) {
            at.push_back(m);
            cum.push_back(cum.back() + left);
            at.push_back(B);
            cum.push_back(cum.back() + right);
            return;
        }
        rSelf(rSelf, A, m, Fa, flm, Fm, left, Depth + 1);
        rSelf(rSelf, m, B, Fm, frm, Fb, right, Depth + 1);
    };
    const auto first = static_cast<std::int64_t>(std::ceil(part / rOptions.mTargetSize));
    double fa = density(L0);
    for (std::int64_t k = 0; k < first; ++k) {
        const double a = L0 + part * double(k) / double(first);
        const double b = k + 1 == first ? L1 : L0 + part * double(k + 1) / double(first);
        const double fm = density(0.5 * (a + b)), fb = density(b);
        adapt(adapt, a, b, fa, fm, fb, simpson(a, b, fa, fm, fb), 0);
        fa = fb;
    }
    const double total = cum.back();
    if (!(total <= double(rBudget)))
        throw std::length_error("meshio++: mesh_cad_surface: an edge needs more than " +
                                std::to_string(rOptions.mMaxSegments) + " segments");
    const std::int64_t n = std::max<std::int64_t>(
        {std::int64_t{1}, static_cast<std::int64_t>(std::ceil(total - 1e-9)),
         static_cast<std::int64_t>(std::ceil(double(rOptions.mMinSegments) * part / Length))});
    double la = L0;
    std::size_t j = 0;
    for (std::int64_t k = 1; k <= n; ++k) {
        double lb = L1;
        if (k < n) {
            const double target = total * double(k) / double(n);
            while (j + 1 < cum.size() && cum[j + 1] < target)
                ++j;
            const double span = cum[j + 1] - cum[j];
            const double f = span > 0.0 ? (target - cum[j]) / span : 0.0;
            lb = at[j] + f * (at[j + 1] - at[j]);
        }
        rSplit(rSplit, la, lb, 0);
        la = lb;
    }
}

}  // namespace

EdgeSamples edge_sample(const TrimTrace& rTrace, const EdgeSampleOptions& rOptions) {
    EdgeSamples out;
    const double length = rTrace.Length();
    out.mLength = length;
    const auto push = [&](double T, double F) {
        out.mParams.push_back(T);
        out.mPoints.push_back(rTrace.Point(T));
        out.mFractions.push_back(F);
    };
    if (!(length > rOptions.mTolerance)) {
        out.mDegenerate = true;
        push(rTrace.Start(), 0.0);
        push(rTrace.End(), 1.0);
        return out;
    }
    if (!(rOptions.mTargetSize > 0.0))
        edgs_bad("target_size must be positive");
    const double pieces = std::ceil(length / rOptions.mTargetSize);
    if (!(pieces <= double(rOptions.mMaxSegments)))
        throw std::length_error("meshio++: mesh_cad_surface: an edge needs more than " +
                                std::to_string(rOptions.mMaxSegments) + " segments");
    // The stretches between the corners, each in equal parts.
    std::vector<double> cuts{0.0};
    for (const double c : rTrace.CornerLengths())
        if (c > cuts.back() + rOptions.mTolerance && c < length - rOptions.mTolerance)
            cuts.push_back(c);
    cuts.push_back(length);
    std::vector<double> ends;  // arc lengths of the pieces' ends, after the first point
    const double dev = rOptions.mMaxChordalDeviation;
    std::int64_t budget = rOptions.mMaxSegments;
    // Halve a piece while its sagitta exceeds the cap (depth-first, so the
    // ends come out in order).
    const auto split = [&](auto&& rSelf, double La, double Lb, int Depth) -> void {
        if (dev > 0.0 && Depth < 40) {
            const EdgeVec3 a = rTrace.Point(rTrace.ParamAtLength(La));
            const EdgeVec3 b = rTrace.Point(rTrace.ParamAtLength(Lb));
            double sag = 0.0;
            for (const double f : {0.25, 0.5, 0.75})
                sag = std::max(sag,
                               edgs_segment_distance(
                                   rTrace.Point(rTrace.ParamAtLength(La + f * (Lb - La))), a, b));
            if (sag > dev) {
                if (--budget < 0)
                    throw std::length_error(
                        "meshio++: mesh_cad_surface: the chordal deviation needs more than " +
                        std::to_string(rOptions.mMaxSegments) + " segments on an edge");
                const double m = 0.5 * (La + Lb);
                rSelf(rSelf, La, m, Depth + 1);
                rSelf(rSelf, m, Lb, Depth + 1);
                return;
            }
        }
        ends.push_back(Lb);
    };
    for (std::size_t c = 0; c + 1 < cuts.size(); ++c) {
        const double l0 = cuts[c];
        const double l1 = cuts[c + 1];
        const double part = l1 - l0;
        if (rOptions.mSize) {
            edgs_sized_stretch(rTrace, rOptions, l0, l1, length, budget, ends, split);
            continue;
        }
        const std::int64_t n = std::max<std::int64_t>(
            {std::int64_t{1}, static_cast<std::int64_t>(std::ceil(part / rOptions.mTargetSize)),
             static_cast<std::int64_t>(std::ceil(double(rOptions.mMinSegments) * part / length))});
        for (std::int64_t k = 0; k < n; ++k) {
            const double la = l0 + part * double(k) / double(n);
            const double lb = k + 1 == n ? l1 : l0 + part * double(k + 1) / double(n);
            split(split, la, lb, 0);
        }
    }
    push(rTrace.Start(), 0.0);
    for (std::size_t i = 0; i + 1 < ends.size(); ++i)
        push(rTrace.ParamAtLength(ends[i]), ends[i] / length);
    push(rTrace.End(), 1.0);
    return out;
}

void edge_refine(const TrimTrace& rTrace, EdgeSamples& rSamples,
                 const std::vector<std::size_t>& rPieces) {
    edge_refine(rTrace, rSamples, rPieces, std::vector<double>(rPieces.size(), 0.5));
}

void edge_refine(const TrimTrace& rTrace, EdgeSamples& rSamples,
                 const std::vector<std::size_t>& rPieces, const std::vector<double>& rAt) {
    if (rPieces.empty() || rSamples.mDegenerate)
        return;
    const double length = rSamples.mLength;
    EdgeSamples out;
    out.mLength = length;
    std::size_t next = 0;
    for (std::size_t k = 0; k < rSamples.mParams.size(); ++k) {
        out.mParams.push_back(rSamples.mParams[k]);
        out.mPoints.push_back(rSamples.mPoints[k]);
        out.mFractions.push_back(rSamples.mFractions[k]);
        if (next < rPieces.size() && rPieces[next] == k && k + 1 < rSamples.mParams.size()) {
            const double a = rAt[next];
            const double f = (1.0 - a) * rSamples.mFractions[k] + a * rSamples.mFractions[k + 1];
            const double t = rTrace.ParamAtLength(f * length);
            out.mParams.push_back(t);
            out.mPoints.push_back(rTrace.Point(t));
            out.mFractions.push_back(f);
            ++next;
        }
    }
    rSamples = std::move(out);
}

std::vector<std::array<double, 3>> edge_piece_points(const TrimTrace& rTrace,
                                                     const EdgeSamples& rSamples, int Order) {
    std::vector<std::array<double, 3>> out;
    if (rSamples.mDegenerate || Order < 2)
        return out;
    const double length = rSamples.mLength;
    for (std::size_t k = 0; k + 1 < rSamples.mParams.size(); ++k)
        for (int j = 1; j < Order; ++j) {
            const double a = double(j) / Order;
            const double f = (1.0 - a) * rSamples.mFractions[k] + a * rSamples.mFractions[k + 1];
            out.push_back(rTrace.Point(rTrace.ParamAtLength(f * length)));
        }
    return out;
}

bool edge_same_direction(const TrimTrace& rA, const TrimTrace& rB, double StartLengthB) {
    const double la = rA.Length();
    const double lb = rB.Length();
    if (StartLengthB >= 0.0) {
        // A rotated closed trace: its ends are the master's start, so only
        // the quarter points tell the directions apart.
        const auto at_b = [&](double F) {
            double s = StartLengthB + F * lb;
            if (s >= lb)
                s -= lb;
            return rB.Point(rB.ParamAtLength(s));
        };
        const EdgeVec3 a1 = rA.Point(rA.ParamAtLength(0.25 * la));
        const EdgeVec3 a3 = rA.Point(rA.ParamAtLength(0.75 * la));
        const EdgeVec3 b1 = at_b(0.25);
        const EdgeVec3 b3 = at_b(0.75);
        return edgs_dist(a1, b1) + edgs_dist(a3, b3) <= edgs_dist(a1, b3) + edgs_dist(a3, b1);
    }
    double same = edgs_dist(rA.Point(rA.Start()), rB.Point(rB.Start())) +
                  edgs_dist(rA.Point(rA.End()), rB.Point(rB.End()));
    double opposite = edgs_dist(rA.Point(rA.Start()), rB.Point(rB.End())) +
                      edgs_dist(rA.Point(rA.End()), rB.Point(rB.Start()));
    if (la > 0.0 && lb > 0.0) {
        const EdgeVec3 a1 = rA.Point(rA.ParamAtLength(0.25 * la));
        const EdgeVec3 a3 = rA.Point(rA.ParamAtLength(0.75 * la));
        const EdgeVec3 b1 = rB.Point(rB.ParamAtLength(0.25 * lb));
        const EdgeVec3 b3 = rB.Point(rB.ParamAtLength(0.75 * lb));
        same += edgs_dist(a1, b1) + edgs_dist(a3, b3);
        opposite += edgs_dist(a1, b3) + edgs_dist(a3, b1);
    }
    return same <= opposite;
}

double edge_closest_length(const TrimTrace& rTrace, const EdgeVec3& rP) {
    const double length = rTrace.Length();
    if (!(length > 0.0))
        return 0.0;
    constexpr int kSteps = 64;
    double best_l = 0.0;
    double best_d = edgs_dist(rTrace.Point(rTrace.Start()), rP);
    for (int i = 1; i < kSteps; ++i) {
        const double l = length * double(i) / double(kSteps);
        const double d = edgs_dist(rTrace.Point(rTrace.ParamAtLength(l)), rP);
        if (d < best_d) {
            best_d = d;
            best_l = l;
        }
    }
    // Newton on (X - p) . X' = 0 within the two neighbouring steps.
    const double step = length / double(kSteps);
    const double t_lo = rTrace.ParamAtLength(std::max(0.0, best_l - step));
    const double t_hi = rTrace.ParamAtLength(std::min(length, best_l + step));
    const double lo = std::min(t_lo, t_hi);
    const double hi = std::max(t_lo, t_hi);
    double t = rTrace.ParamAtLength(best_l);
    for (int it = 0; it < 30; ++it) {
        const EdgeVec3 x = rTrace.Point(t);
        const EdgeVec3 v = rTrace.Velocity(t);
        const double vv = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        if (!(vv > 0.0))
            break;
        const double dt =
            -((x[0] - rP[0]) * v[0] + (x[1] - rP[1]) * v[1] + (x[2] - rP[2]) * v[2]) / vv;
        const double next = std::clamp(t + dt, lo, hi);
        if (std::fabs(next - t) <= 1e-15 * (hi - lo)) {
            t = next;
            break;
        }
        t = next;
    }
    return edgs_dist(rTrace.Point(t), rP) <= best_d ? rTrace.LengthAt(t) : best_l;
}

EdgeProjection edge_project(const TrimTrace& rOther, const EdgeSamples& rMaster, bool SameDirection,
                            double StartLength) {
    EdgeProjection out;
    const std::size_t m = rMaster.mParams.size();
    out.mParams.assign(m, 0.0);
    if (m == 0)
        return out;
    const auto point_at = [&](std::size_t K) -> const EdgeVec3& {
        return SameDirection ? rMaster.mPoints[K] : rMaster.mPoints[m - 1 - K];
    };
    const auto fraction_at = [&](std::size_t K) {
        return SameDirection ? rMaster.mFractions[K] : 1.0 - rMaster.mFractions[m - 1 - K];
    };
    const double length = rOther.Length();
    const bool rotated = StartLength >= 0.0;
    const double s0 = rotated ? std::clamp(StartLength, 0.0, length) : 0.0;
    // A rotated trace starts and ends at s0, and runs past its own end and on
    // from its own start: two monotone pieces.
    const double start = rotated ? rOther.ParamAtLength(s0) : rOther.Start();
    const double end = rotated ? start : rOther.End();
    const double span = std::fabs(rOther.End() - rOther.Start());
    out.mParams.front() = start;
    out.mParams.back() = end;
    out.mMaxGap = std::max(edgs_dist(rOther.Point(start), point_at(0)),
                           edgs_dist(rOther.Point(end), point_at(m - 1)));
    double prev = start;
    bool wrapped = false;
    for (std::size_t k = 1; k + 1 < m; ++k) {
        const EdgeVec3& p = point_at(k);
        double s = s0 + fraction_at(k) * length;
        double piece_end = rOther.End();
        if (rotated && s >= length) {
            s -= length;
            piece_end = start;
            if (!wrapped) {
                wrapped = true;
                prev = rOther.Start();
            }
        } else if (!rotated) {
            piece_end = end;
        }
        double t = rOther.ParamAtLength(s);
        const double lo = std::min(prev, piece_end);
        const double hi = std::max(prev, piece_end);
        t = std::clamp(t, lo, hi);
        for (int it = 0; it < 30; ++it) {
            const EdgeVec3 x = rOther.Point(t);
            const EdgeVec3 v = rOther.Velocity(t);
            const double vv = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            if (!(vv > 0.0))
                break;
            const double dt =
                -((x[0] - p[0]) * v[0] + (x[1] - p[1]) * v[1] + (x[2] - p[2]) * v[2]) / vv;
            const double next = std::clamp(t + dt, lo, hi);
            if (std::fabs(next - t) <= 1e-15 * span) {
                t = next;
                break;
            }
            t = next;
        }
        out.mParams[k] = t;
        out.mMaxGap = std::max(out.mMaxGap, edgs_dist(rOther.Point(t), p));
        prev = t;
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
