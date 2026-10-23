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
