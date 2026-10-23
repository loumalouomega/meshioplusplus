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
    mT.push_back(breaks.front());
    mCum.push_back(0.0);
    for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
        const double a = breaks[i];
        const double b = breaks[i + 1];
        constexpr int kStart = 4;
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
    if (Depth < 30 && std::fabs(Whole - halves) > 1e-12 * halves && m != A && m != B) {
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

bool edge_same_direction(const TrimTrace& rA, const TrimTrace& rB) {
    const double la = rA.Length();
    const double lb = rB.Length();
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

EdgeProjection edge_project(const TrimTrace& rOther, const EdgeSamples& rMaster,
                            bool SameDirection) {
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
    const double start = rOther.Start();
    const double end = rOther.End();
    out.mParams.front() = start;
    out.mParams.back() = end;
    out.mMaxGap = std::max(edgs_dist(rOther.Point(start), point_at(0)),
                           edgs_dist(rOther.Point(end), point_at(m - 1)));
    const double length = rOther.Length();
    const double span = std::fabs(end - start);
    double prev = start;
    for (std::size_t k = 1; k + 1 < m; ++k) {
        const EdgeVec3& p = point_at(k);
        double t = rOther.ParamAtLength(fraction_at(k) * length);
        const double lo = std::min(prev, end);
        const double hi = std::max(prev, end);
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
