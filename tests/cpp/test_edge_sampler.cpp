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
// detail/edge_sampler: trims composed with their surfaces, sampled by arc
// length and carried to the other uses of an edge. A straight line on a
// scaled plane and a rational quarter circle (whose length and chords are
// known exactly) check the arc length, the equal spacing and the chordal cap;
// the same circle written on a second plane with another parametrisation,
// run backwards, checks that every use lands on the master's points; a
// patch with a collapsed row checks a zero-length edge.

// System includes
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/nurbs.hpp"
#include "../../src/cpp/src/detail/edge_sampler.hpp"
#include "brep_solids.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::detail::BrepCurveEval;
using mio::detail::BrepSurfaceEval;
using mio::detail::EdgeSampleOptions;
using mio::detail::EdgeVec3;
using mio::detail::TrimTrace;

// A bilinear patch over [0, 1]^2: X = (X0 + Sx u, Y0 + Sy v, 0); `Pole`
// collapses the v = 0 row onto (X0, Y0, 0).
mio::NurbsPatch plane(double X0, double Y0, double Sx, double Sy, bool Pole = false) {
    mio::NurbsPatch p;
    p.mDegrees = {1, 1};
    p.mKnots = {{0, 0, 1, 1}, {0, 0, 1, 1}};
    p.mControlPoints = mio::NDArray(mio::DType::Float64, {4, 3});
    double* c = p.mControlPoints.As<double>();
    const double pts[4][3] = {
        {X0, Y0, 0}, {Pole ? X0 : X0 + Sx, Y0, 0}, {X0, Y0 + Sy, 0}, {X0 + Sx, Y0 + Sy, 0}};
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 3; ++k)
            c[i * 3 + k] = pts[i][k];
    p.mWeights = mio::NDArray(mio::DType::Float64, {4});
    for (int i = 0; i < 4; ++i)
        p.mWeights.As<double>()[i] = 1.0;
    return p;
}

mio::BrepCurve curve(int Degree, const std::vector<double>& rKnots,
                     const std::vector<std::array<double, 2>>& rPoles,
                     const std::vector<double>& rWeights = {}) {
    mio::BrepCurve c;
    c.mDegree = Degree;
    c.mKnots = rKnots;
    c.mRational = !rWeights.empty();
    c.mControlPoints = mio::NDArray(mio::DType::Float64, {rPoles.size(), 3});
    c.mWeights = mio::NDArray(mio::DType::Float64, {rPoles.size()});
    for (std::size_t i = 0; i < rPoles.size(); ++i) {
        c.mControlPoints.As<double>()[i * 3] = rPoles[i][0];
        c.mControlPoints.As<double>()[i * 3 + 1] = rPoles[i][1];
        c.mControlPoints.As<double>()[i * 3 + 2] = 0.0;
        c.mWeights.As<double>()[i] = rWeights.empty() ? 1.0 : rWeights[i];
    }
    return c;
}

// A counter-clockwise quarter circle in (u, v) about (Cu, Cv), from angle 0,
// over the parameter range [0, T].
mio::BrepCurve quarter(double Cu, double Cv, double R, double T = 1.0) {
    const double w = std::sqrt(0.5);
    return curve(2, {0, 0, 0, T, T, T}, {{Cu + R, Cv}, {Cu + R, Cv + R}, {Cu, Cv + R}}, {1, w, 1});
}

// A full circle of four rational quarters, from angle 0, over [0, 4 T].
mio::BrepCurve circle(double Cu, double Cv, double R, double T = 1.0) {
    const double w = std::sqrt(0.5);
    return curve(2, {0, 0, 0, T, T, 2 * T, 2 * T, 3 * T, 3 * T, 4 * T, 4 * T, 4 * T},
                 {{Cu + R, Cv},
                  {Cu + R, Cv + R},
                  {Cu, Cv + R},
                  {Cu - R, Cv + R},
                  {Cu - R, Cv},
                  {Cu - R, Cv - R},
                  {Cu, Cv - R},
                  {Cu + R, Cv - R},
                  {Cu + R, Cv}},
                 {1, w, 1, w, 1, w, 1, w, 1});
}

// The same circle started a number of quarter turns later.
mio::BrepCurve circle_from(double Cu, double Cv, double R, int Quarters, double T = 1.0) {
    const double w = std::sqrt(0.5);
    std::vector<std::array<double, 2>> poles = {
        {Cu + R, Cv}, {Cu + R, Cv + R}, {Cu, Cv + R}, {Cu - R, Cv + R},
        {Cu - R, Cv}, {Cu - R, Cv - R}, {Cu, Cv - R}, {Cu + R, Cv - R}};
    std::rotate(poles.begin(), poles.begin() + 2 * Quarters, poles.end());
    poles.push_back(poles.front());
    return curve(2, {0, 0, 0, T, T, 2 * T, 2 * T, 3 * T, 3 * T, 4 * T, 4 * T, 4 * T}, poles,
                 {1, w, 1, w, 1, w, 1, w, 1});
}

double dist(const EdgeVec3& rA, const EdgeVec3& rB) {
    return std::sqrt((rA[0] - rB[0]) * (rA[0] - rB[0]) + (rA[1] - rB[1]) * (rA[1] - rB[1]) +
                     (rA[2] - rB[2]) * (rA[2] - rB[2]));
}

TEST(EdgeSampler, LineOnAScaledPlane) {
    const mio::NurbsPatch p = plane(0, 0, 2, 3);
    const BrepSurfaceEval s(p, false);
    const TrimTrace t(&s, BrepCurveEval(curve(1, {0, 0, 1, 1}, {{0, 0}, {1, 0}})), true);
    EXPECT_NEAR(t.Length(), 2.0, 1e-14);
    EdgeSampleOptions o;
    o.mTargetSize = 0.3;
    const auto e = mio::detail::edge_sample(t, o);
    ASSERT_EQ(e.mParams.size(), 8u);  // ceil(2 / 0.3) = 7 pieces
    for (std::size_t k = 0; k + 1 < e.mPoints.size(); ++k) {
        EXPECT_NEAR(dist(e.mPoints[k], e.mPoints[k + 1]), 2.0 / 7.0, 1e-12);
        EXPECT_NEAR(e.mFractions[k], double(k) / 7.0, 1e-14);
    }
    o.mTargetSize = 10.0;
    o.mMinSegments = 5;
    EXPECT_EQ(mio::detail::edge_sample(t, o).mParams.size(), 6u);
}

TEST(EdgeSampler, RationalQuarterCircle) {
    // The plane doubles (u, v): the arc of radius 0.375 is one of radius 0.75.
    const mio::NurbsPatch p = plane(0, 0, 2, 2);
    const BrepSurfaceEval s(p, false);
    const TrimTrace t(&s, BrepCurveEval(quarter(0.5, 0.5, 0.375)), true);
    const double r = 0.75;
    EXPECT_NEAR(t.Length(), 0.5 * M_PI * r, 1e-12 * r);
    EdgeSampleOptions o;
    o.mTargetSize = 0.1;
    const auto e = mio::detail::edge_sample(t, o);
    const auto n = static_cast<double>(e.mPoints.size() - 1);
    EXPECT_EQ(n, std::ceil(0.5 * M_PI * r / 0.1));
    const double chord = 2.0 * r * std::sin(0.25 * M_PI / n);
    for (std::size_t k = 0; k < e.mPoints.size(); ++k) {
        const EdgeVec3& q = e.mPoints[k];
        EXPECT_NEAR(std::hypot(q[0] - 1.0, q[1] - 1.0), r, 1e-12);
        if (k + 1 < e.mPoints.size()) {
            EXPECT_NEAR(dist(q, e.mPoints[k + 1]), chord, 1e-10);
        }
    }
}

TEST(EdgeSampler, RefineSplitsPiecesByArcLengthAtAChosenFraction) {
    // Halving (a crossing to separate) and a split a quarter of the way along
    // (a concentric shell, a power of two from a sharp corner) both land at
    // that fraction of the piece's arc length, on the rational arc.
    const mio::NurbsPatch p = plane(0, 0, 2, 2);
    const BrepSurfaceEval s(p, false);
    const TrimTrace t(&s, BrepCurveEval(quarter(0.5, 0.5, 0.375)), true);
    const double r = 0.75;
    EdgeSampleOptions o;
    o.mTargetSize = 0.4;
    auto e = mio::detail::edge_sample(t, o);
    const std::size_t n = e.mPoints.size() - 1;
    ASSERT_GE(n, 3u);
    const double piece = 0.5 * M_PI * r / double(n);
    mio::detail::edge_refine(t, e, {0, n - 1}, {0.5, 0.25});
    ASSERT_EQ(e.mPoints.size(), n + 3);
    const auto arc = [&](std::size_t A, std::size_t B) {
        // The angle between two points on the arc, times its radius.
        const EdgeVec3& a = e.mPoints[A];
        const EdgeVec3& b = e.mPoints[B];
        const double ca = std::atan2(a[1] - 1.0, a[0] - 1.0);
        const double cb = std::atan2(b[1] - 1.0, b[0] - 1.0);
        return r * std::fabs(cb - ca);
    };
    EXPECT_NEAR(arc(0, 1), 0.5 * piece, 1e-10);
    EXPECT_NEAR(arc(1, 2), 0.5 * piece, 1e-10);
    EXPECT_NEAR(arc(n, n + 1), 0.25 * piece, 1e-10);
    EXPECT_NEAR(arc(n + 1, n + 2), 0.75 * piece, 1e-10);
    for (std::size_t k = 0; k + 1 < e.mFractions.size(); ++k)
        EXPECT_LT(e.mFractions[k], e.mFractions[k + 1]);
}

TEST(EdgeSampler, ChordalCapIsHonoured) {
    const mio::NurbsPatch p = plane(0, 0, 2, 2);
    const BrepSurfaceEval s(p, false);
    const TrimTrace t(&s, BrepCurveEval(quarter(0.5, 0.5, 0.375)), true);
    const double r = 0.75;
    for (const double dev : {1e-2, 1e-3, 1e-4}) {
        EdgeSampleOptions o;
        o.mTargetSize = 10.0;  // one piece by length
        o.mMaxChordalDeviation = dev;
        const auto e = mio::detail::edge_sample(t, o);
        for (std::size_t k = 0; k + 1 < e.mPoints.size(); ++k) {
            const double half = std::asin(0.5 * dist(e.mPoints[k], e.mPoints[k + 1]) / r);
            EXPECT_LE(r * (1.0 - std::cos(half)), dev * (1.0 + 1e-9));
        }
        const double optimal = std::ceil(0.5 * M_PI / (2.0 * std::acos(1.0 - dev / r)));
        EXPECT_LE(double(e.mPoints.size() - 1), 2.0 * optimal);
    }
}

TEST(EdgeSampler, AnotherUseLandsOnTheMasterPoints) {
    // The same 3-D arc (centre (1, 1), radius 0.75) as a trim of a second
    // plane X = (-1 + 4u, -1 + 4v), over the parameter range [0, 7] and run
    // backwards along its loop.
    const mio::NurbsPatch pa = plane(0, 0, 2, 2);
    const mio::NurbsPatch pb = plane(-1, -1, 4, 4);
    const BrepSurfaceEval sa(pa, false);
    const BrepSurfaceEval sb(pb, false);
    const TrimTrace master(&sa, BrepCurveEval(quarter(0.5, 0.5, 0.375)), true);
    const TrimTrace other(&sb, BrepCurveEval(quarter(0.5, 0.5, 0.1875, 7.0)), false);
    EXPECT_NEAR(master.Length(), other.Length(), 1e-12);
    EXPECT_FALSE(mio::detail::edge_same_direction(master, other));
    EdgeSampleOptions o;
    o.mTargetSize = 0.05;
    const auto e = mio::detail::edge_sample(master, o);
    const auto proj = mio::detail::edge_project(other, e, false);
    ASSERT_EQ(proj.mParams.size(), e.mParams.size());
    EXPECT_LT(proj.mMaxGap, 1e-12);
    const std::size_t m = e.mParams.size();
    for (std::size_t k = 0; k < m; ++k) {
        EXPECT_LT(dist(other.Point(proj.mParams[k]), e.mPoints[m - 1 - k]), 1e-12);
        if (k > 0) {
            EXPECT_LT(proj.mParams[k], proj.mParams[k - 1]);  // runs backwards from 7 to 0
        }
    }
}

TEST(EdgeSampler, ClosedEdgeDirections) {
    const mio::NurbsPatch pa = plane(0, 0, 2, 2);
    const mio::NurbsPatch pb = plane(-1, -1, 4, 4);
    const BrepSurfaceEval sa(pa, true);
    const BrepSurfaceEval sb(pb, true);
    const TrimTrace master(&sa, BrepCurveEval(circle(0.5, 0.5, 0.375)), true);
    const TrimTrace same(&sb, BrepCurveEval(circle(0.5, 0.5, 0.1875, 2.5)), true);
    const TrimTrace back(&sb, BrepCurveEval(circle(0.5, 0.5, 0.1875, 2.5)), false);
    EXPECT_NEAR(master.Length(), 2.0 * M_PI * 0.75, 1e-12);
    EXPECT_TRUE(mio::detail::edge_same_direction(master, same));
    EXPECT_FALSE(mio::detail::edge_same_direction(master, back));
    EdgeSampleOptions o;
    o.mTargetSize = 0.2;
    const auto e = mio::detail::edge_sample(master, o);
    EXPECT_LT(mio::detail::edge_project(same, e, true).mMaxGap, 1e-12);
    EXPECT_LT(mio::detail::edge_project(back, e, false).mMaxGap, 1e-12);
}

TEST(EdgeSampler, LengthAtInvertsParamAtLength) {
    const mio::NurbsPatch pa = plane(0, 0, 2, 2);
    const BrepSurfaceEval sa(pa, true);
    for (const bool forward : {true, false}) {
        const TrimTrace t(&sa, BrepCurveEval(circle(0.5, 0.5, 0.375, 2.5)), forward);
        for (int i = 0; i <= 16; ++i) {
            const double l = t.Length() * double(i) / 16.0;
            EXPECT_NEAR(t.LengthAt(t.ParamAtLength(l)), l, 1e-12) << forward << " " << i;
        }
        EXPECT_EQ(t.LengthAt(t.Start()), 0.0);
        EXPECT_DOUBLE_EQ(t.LengthAt(t.End()), t.Length());
    }
}

TEST(EdgeSampler, AClosedUseStartedElsewhereIsRotatedOntoTheMaster) {
    // The master circle starts at angle 0; the other use is the same 3-D
    // circle, its trim on another plane started a quarter, a half or three
    // quarters of a turn later, run either way. Rotated to start where the
    // master does, it lands on the master's points in two monotone pieces.
    const mio::NurbsPatch pa = plane(0, 0, 2, 2);
    const mio::NurbsPatch pb = plane(-1, -1, 4, 4);
    const BrepSurfaceEval sa(pa, true);
    const BrepSurfaceEval sb(pb, true);
    const TrimTrace master(&sa, BrepCurveEval(circle(0.5, 0.5, 0.375)), true);
    EdgeSampleOptions o;
    o.mTargetSize = 0.15;
    const auto e = mio::detail::edge_sample(master, o);
    const std::size_t m = e.mParams.size();
    const EdgeVec3 start = master.Point(master.Start());
    for (int q = 1; q < 4; ++q)
        for (const bool forward : {true, false}) {
            const TrimTrace other(&sb, BrepCurveEval(circle_from(0.5, 0.5, 0.1875, q, 2.5)),
                                  forward);
            EXPECT_GT(dist(other.Point(other.Start()), start), 0.5);
            const double s0 = mio::detail::edge_closest_length(other, start);
            EXPECT_LT(dist(other.Point(other.ParamAtLength(s0)), start), 1e-12);
            const double turn = forward ? (4 - q) / 4.0 : q / 4.0;
            EXPECT_NEAR(s0, turn * other.Length(), 1e-10) << q << " " << forward;
            const bool same = mio::detail::edge_same_direction(master, other, s0);
            EXPECT_EQ(same, forward);
            const auto proj = mio::detail::edge_project(other, e, same, s0);
            EXPECT_LT(proj.mMaxGap, 1e-12) << q << " " << forward;
            EXPECT_EQ(proj.mParams.front(), proj.mParams.back());
            int wraps = 0;
            for (std::size_t k = 0; k < m; ++k) {
                const EdgeVec3& want = same ? e.mPoints[k] : e.mPoints[m - 1 - k];
                EXPECT_LT(dist(other.Point(proj.mParams[k]), want), 1e-12);
                if (k > 0 && k + 1 < m) {
                    const double d = proj.mParams[k] - proj.mParams[k - 1];
                    if ((d < 0.0) == forward)
                        ++wraps;
                }
            }
            EXPECT_EQ(wraps, 1) << "one wrap past the trace's end";
        }
}

TEST(EdgeSampler, ZeroLengthEdge) {
    const mio::NurbsPatch p = plane(0, 0, 1, 1, /*Pole=*/true);
    const BrepSurfaceEval s(p, false);
    const TrimTrace t(&s, BrepCurveEval(curve(1, {0, 0, 1, 1}, {{0, 0}, {1, 0}})), true);
    EXPECT_EQ(t.Length(), 0.0);
    EdgeSampleOptions o;
    o.mTargetSize = 0.1;
    o.mTolerance = 1e-9;
    const auto e = mio::detail::edge_sample(t, o);
    EXPECT_TRUE(e.mDegenerate);
    ASSERT_EQ(e.mParams.size(), 2u);
    EXPECT_EQ(dist(e.mPoints[0], e.mPoints[1]), 0.0);
}

// A short trace on a patch far from the origin: its speed carries round-off
// of 1e-6 of itself, which no relative tolerance meets. Its table halved to
// the depth limit (2^30 pieces per span, hours) before v17.18.0; now it stops
// at the cap and its length is still right to the round-off.
TEST(EdgeSampler, RoundOffInTheSpeedDoesNotHalveForever) {
    const mio::NurbsPatch p = plane(1e6, 1e6, 1e-4, 1e-4);
    const BrepSurfaceEval s(p, false);
    const TrimTrace t(&s, BrepCurveEval(curve(1, {0, 0, 1, 1}, {{0.1, 0.2}, {0.7, 0.9}})), true);
    const double length = 1e-4 * std::hypot(0.6, 0.7);
    EXPECT_NEAR(t.Length(), length, length * 1e-5);
    EdgeSampleOptions o;
    o.mTargetSize = 0.25 * length * (1.0 + 1e-3);
    EXPECT_EQ(mio::detail::edge_sample(t, o).mParams.size(), 5u);
}

TEST(EdgeSampler, CornersAreKept) {
    // A polyline (0, 0) -> (1, 0) -> (1, 0.3): its corner at the knot 0.5 is
    // a sample whatever the spacing, so no chord cuts it.
    const mio::NurbsPatch p = plane(0, 0, 1, 1);
    const BrepSurfaceEval s(p, false);
    const TrimTrace t(&s, BrepCurveEval(curve(1, {0, 0, 0.5, 1, 1}, {{0, 0}, {1, 0}, {1, 0.3}})),
                      true);
    EXPECT_NEAR(t.Length(), 1.3, 1e-14);
    EdgeSampleOptions o;
    o.mTargetSize = 0.4;
    const auto e = mio::detail::edge_sample(t, o);
    // ceil(1 / 0.4) = 3 parts before the corner, ceil(0.3 / 0.4) = 1 after.
    ASSERT_EQ(e.mPoints.size(), 5u);
    EXPECT_NEAR(e.mPoints[3][0], 1.0, 1e-14);
    EXPECT_NEAR(e.mPoints[3][1], 0.0, 1e-14);
}

TEST(EdgeSampler, Refusals) {
    const mio::NurbsPatch p = plane(0, 0, 1, 1);
    const BrepSurfaceEval s(p, false);
    const TrimTrace t(&s, BrepCurveEval(curve(1, {0, 0, 1, 1}, {{0, 0}, {1, 0}})), true);
    EdgeSampleOptions o;
    o.mTargetSize = 1e-9;
    o.mMaxSegments = 1000;
    EXPECT_THROW(mio::detail::edge_sample(t, o), std::length_error);
    o.mTargetSize = 0.0;
    EXPECT_THROW(mio::detail::edge_sample(t, o), std::invalid_argument);
    EXPECT_THROW(BrepCurveEval(curve(1, {0, 0, 1}, {{0, 0}, {1, 0}})), std::invalid_argument);
    mio::NurbsPatch bad = p;
    bad.mKnots[0] = {0, 0, 0, 1};
    EXPECT_THROW(BrepSurfaceEval(bad, false), std::invalid_argument);
}

}  // namespace

// Second derivatives and curvature (roadmap 6.3.2): the rational surfaces of
// the hand-built solids against central differences of the first
// derivatives, and their principal curvatures against the analytic ones
// (A2.3 and A4.4, Piegl and Tiller 1997).
namespace {

void expect_second_derivatives_match(const mio::BrepInfo& rInfo) {
    for (const mio::BrepFace& r_f : rInfo.mFaces) {
        const mio::NurbsPatch& r_p = rInfo.mPatches[static_cast<std::size_t>(r_f.mPatch)];
        const BrepSurfaceEval s(r_p, r_f.mRational);
        const auto& box = s.Box();
        // Off the knots (the circles' double knots at the quarters, where the
        // second derivative jumps and a central difference straddles both).
        for (const double fu : {0.13, 0.37, 0.61, 0.89})
            for (const double fv : {0.21, 0.43, 0.79}) {
                const double u = box[0] + fu * (box[1] - box[0]);
                const double v = box[2] + fv * (box[3] - box[2]);
                double x[3], su[3], sv[3], suu[3], suv[3], svv[3];
                s.Eval2(u, v, x, su, sv, suu, suv, svv);
                double x1[3], su1[3], sv1[3];
                s.Eval(u, v, x1, su1, sv1);
                const double hu = 1e-5 * (box[1] - box[0]), hv = 1e-5 * (box[3] - box[2]);
                double xa[3], ua[3], va[3], xb[3], ub[3], vb[3];
                s.Eval(u + hu, v, xa, ua, va);
                s.Eval(u - hu, v, xb, ub, vb);
                double xc[3], uc[3], vc[3], xd[3], ud[3], vd[3];
                s.Eval(u, v + hv, xc, uc, vc);
                s.Eval(u, v - hv, xd, ud, vd);
                for (int k = 0; k < 3; ++k) {
                    EXPECT_NEAR(x[k], x1[k], 1e-12);
                    EXPECT_NEAR(su[k], su1[k], 1e-10);
                    EXPECT_NEAR(sv[k], sv1[k], 1e-10);
                    const double scale_u = 1.0 + std::abs(suu[k]);
                    const double scale_v = 1.0 + std::abs(svv[k]);
                    EXPECT_NEAR(suu[k], (ua[k] - ub[k]) / (2 * hu), 1e-5 * scale_u);
                    EXPECT_NEAR(suv[k], (va[k] - vb[k]) / (2 * hu),
                                1e-5 * (1.0 + std::abs(suv[k])));
                    EXPECT_NEAR(suv[k], (uc[k] - ud[k]) / (2 * hv),
                                1e-5 * (1.0 + std::abs(suv[k])));
                    EXPECT_NEAR(svv[k], (vc[k] - vd[k]) / (2 * hv), 1e-5 * scale_v);
                }
            }
    }
}

}  // namespace

TEST(EdgeSampler, SecondDerivativesMatchCentralDifferences) {
    expect_second_derivatives_match(brep_solids::sphere(1.5));
    expect_second_derivatives_match(brep_solids::capped_cylinder(0.7, 2.0));
    expect_second_derivatives_match(brep_solids::box());
}

TEST(EdgeSampler, PrincipalCurvaturesOfTheAnalyticSurfaces) {
    const auto max_curvature = [](const mio::BrepInfo& rInfo, std::size_t Face, double Fu,
                                  double Fv) {
        const mio::BrepFace& r_f = rInfo.mFaces[Face];
        const BrepSurfaceEval s(rInfo.mPatches[static_cast<std::size_t>(r_f.mPatch)],
                                r_f.mRational);
        const auto& box = s.Box();
        return s.MaxCurvature(box[0] + Fu * (box[1] - box[0]), box[2] + Fv * (box[3] - box[2]));
    };
    const mio::BrepInfo sphere = brep_solids::sphere(2.0);
    for (const double fu : {0.1, 0.4, 0.9})
        for (const double fv : {0.2, 0.45, 0.7})
            // A sphere is all umbilics, where sqrt(H^2 - K) of two equal
            // curvatures lifts rounding to about sqrt(eps).
            EXPECT_NEAR(max_curvature(sphere, 0, fu, fv), 0.5, 1e-7) << fu << " " << fv;
    const mio::BrepInfo cylinder = brep_solids::capped_cylinder(0.25, 1.0);
    // The band is face 0, the discs are flat.
    EXPECT_NEAR(max_curvature(cylinder, 0, 0.3, 0.6), 4.0, 1e-9);
    EXPECT_NEAR(max_curvature(cylinder, 1, 0.3, 0.6), 0.0, 1e-12);
    EXPECT_NEAR(max_curvature(brep_solids::box(), 2, 0.5, 0.5), 0.0, 1e-12);
}
