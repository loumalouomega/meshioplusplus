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

// Kernel-free trimmed NURBS faces (detail/brep_nurbs).

// System includes
#include <cmath>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/cad_json.hpp"

// Project includes (private, not installed)
#include "../../src/cpp/src/detail/brep_nurbs.hpp"
#include "../../src/cpp/src/detail/edge_sampler.hpp"

namespace {

using meshioplusplus::BrepCurve;
using meshioplusplus::BrepInfo;
using meshioplusplus::detail::BrepNurbsFace;
using meshioplusplus::detail::BrepNurbsSurface;
using meshioplusplus::detail::BrepNurbsTrim;
using meshioplusplus::detail::PatchVec;

constexpr double kPi = 3.14159265358979323846;

// A circle of radius R at height Z about the z axis, from angle A0 by Sweep.
BrepNurbsTrim circle_trim(double R, double Z, double A0, double Sweep) {
    BrepNurbsTrim t;
    t.mAt = [=](double T) {
        const double a = A0 + Sweep * T;
        return PatchVec{R * std::cos(a), R * std::sin(a), Z};
    };
    return t;
}

// The distance of every point of every trim's parameter curve, mapped onto
// the patch, from the given analytic surface.
template <class F>
double trim_deviation(const BrepNurbsFace& rF, F&& rDistance) {
    double worst = 0.0;
    for (const auto& loop : rF.mFace.mLoops)
        for (const auto& trim : loop.mTrims) {
            const BrepCurve& c = trim.mCurve;
            const double a = c.mKnots.front(), b = c.mKnots.back();
            for (int i = 0; i <= 64; ++i) {
                const PatchVec uv =
                    meshioplusplus::detail::brep_curve_point(c, a + (b - a) * i / 64.0);
                const PatchVec p =
                    meshioplusplus::detail::nurbs_patch_point(rF.mPatch, uv[0], uv[1]);
                worst = std::max(worst, rDistance(p));
            }
        }
    return worst;
}

// The largest distance from the circle of radius R about the origin of a
// trim's image through the face's patch, evaluated as the mesher does (its
// parameters clamped to the patch).
double circle_deviation(const BrepNurbsFace& rF, const meshioplusplus::BrepTrim& rTrim, double R) {
    const meshioplusplus::detail::BrepSurfaceEval s(rF.mPatch, rF.mFace.mRational);
    const meshioplusplus::detail::BrepCurveEval c(rTrim.mCurve);
    double worst = 0.0;
    for (int i = 0; i <= 2000; ++i) {
        double uv[3], d[3], x[3];
        c.Eval(c.T0() + (c.T1() - c.T0()) * i / 2000.0, uv, d);
        s.Eval(uv[0], uv[1], x, nullptr, nullptr);
        worst = std::max(worst, std::fabs(std::hypot(x[0], x[1]) - R));
    }
    return worst;
}

// The area of a sampled B-rep's face cells (by their corners).
double sampled_area(const BrepInfo& rInfo) {
    meshioplusplus::CadJsonReadOptions o;
    o.mSubdivisions = 48;
    o.mEdges = false;
    const meshioplusplus::Mesh m = meshioplusplus::sample_brep(rInfo, o);
    const double* xyz = m.Points().As<double>();
    double area = 0.0;
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        const auto cells = m.Cells(b);
        const meshioplusplus::NDArray& conn = cells.Conn();
        const std::size_t width = cells.NodesPerCell();
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            // The four corners of a Lagrange quadrilateral come first.
            const double* p[4];
            for (std::size_t k = 0; k < 4; ++k)
                p[k] = xyz + 3 * meshioplusplus::detail::read_int(conn, c * width + k);
            for (int t = 0; t < 2; ++t) {
                const double* a = p[0];
                const double* q = p[t + 1];
                const double* r = p[t + 2];
                const double u[3] = {q[0] - a[0], q[1] - a[1], q[2] - a[2]};
                const double v[3] = {r[0] - a[0], r[1] - a[1], r[2] - a[2]};
                const double x[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                     u[0] * v[1] - u[1] * v[0]};
                area += 0.5 * std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
            }
        }
    }
    return area;
}

// A one-face B-rep (no edges) for sampling.
BrepInfo one_face(const BrepNurbsFace& rF) {
    BrepInfo info;
    info.mPatches.push_back(rF.mPatch);
    meshioplusplus::BrepFace f = rF.mFace;
    f.mId = 1;
    f.mPatch = 0;
    info.mFaces.push_back(f);
    meshioplusplus::Brep b;
    b.mId = 2;
    b.mFaces = {0};
    info.mBreps.push_back(b);
    return info;
}

}  // namespace

TEST(BrepNurbs, UnitArcInvertsItsParameter) {
    const double ranges[][2] = {{0.0, kPi / 3}, {-1.0, 5.0}, {0.0, 2 * kPi}, {0.3, 0.2 + 4 * kPi}};
    for (const auto& r : ranges) {
        const meshioplusplus::detail::UnitArc arc = meshioplusplus::detail::unit_arc(r[0], r[1]);
        const BrepCurve c = meshioplusplus::detail::brep_curve_arc(
            {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0, r[0], r[1]);
        for (int i = 0; i <= 50; ++i) {
            const double a = r[0] + (r[1] - r[0]) * i / 50.0;
            const PatchVec p = meshioplusplus::detail::brep_curve_point(
                c, meshioplusplus::detail::unit_arc_param(arc, a));
            EXPECT_NEAR(p[0], std::cos(a), 1e-12);
            EXPECT_NEAR(p[1], std::sin(a), 1e-12);
        }
    }
    // Backwards, the same points in reverse.
    const BrepCurve f =
        meshioplusplus::detail::brep_curve_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 2.0, 0.0, 1.0);
    const BrepCurve b =
        meshioplusplus::detail::brep_curve_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 2.0, 1.0, 0.0);
    const PatchVec pf = meshioplusplus::detail::brep_curve_point(f, 0.0);
    const PatchVec pb = meshioplusplus::detail::brep_curve_point(b, 1.0);
    EXPECT_NEAR(pf[0], pb[0], 1e-14);
    EXPECT_NEAR(pf[1], pb[1], 1e-14);
}

TEST(BrepNurbs, CylinderBandClosesThroughASeam) {
    BrepNurbsSurface s;
    s.mKind = 3;
    s.mRadius = 2.0;
    // The bottom rim runs anticlockwise, the top clockwise, from other angles.
    std::vector<std::vector<BrepNurbsTrim>> loops(2);
    loops[0].push_back(circle_trim(2.0, 0.0, 0.4, 2 * kPi));
    loops[0].back().mEdge = 0;
    loops[1].push_back(circle_trim(2.0, 3.0, 2.0, -2 * kPi));
    loops[1].back().mEdge = 1;
    BrepNurbsFace face;
    std::string why;
    ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why)) << why;
    ASSERT_EQ(face.mFace.mLoops.size(), 1u);
    EXPECT_EQ(face.mFace.mLoops[0].mType, "outer");
    ASSERT_EQ(face.mFace.mLoops[0].mTrims.size(), 4u);
    ASSERT_EQ(face.mSeams.size(), 1u);
    ASSERT_EQ(face.mRefs.size(), 4u);
    EXPECT_EQ(face.mRefs[0].mEdge, 0);
    EXPECT_EQ(face.mRefs[1].mSeam, 0);
    EXPECT_TRUE(face.mRefs[1].mAlong);
    EXPECT_EQ(face.mRefs[3].mSeam, 0);
    EXPECT_FALSE(face.mRefs[3].mAlong);
    // The rims are parameter lines: two points each.
    EXPECT_EQ(face.mFace.mLoops[0].mTrims[0].mCurve.mControlPoints.Shape()[0], 2u);
    const double dev = trim_deviation(
        face, [](const PatchVec& rP) { return std::fabs(std::hypot(rP[0], rP[1]) - 2.0); });
    EXPECT_LT(dev, 1e-12);
    // The face is the whole band, once.
    EXPECT_NEAR(sampled_area(one_face(face)), 2 * kPi * 2.0 * 3.0, 0.01 * 2 * kPi * 6.0);
}

TEST(BrepNurbs, SphereCapClosesThroughThePole) {
    BrepNurbsSurface s;
    s.mKind = 2;
    s.mRadius = 1.0;
    const double z = 0.5, r = std::sqrt(1.0 - z * z);
    std::vector<std::vector<BrepNurbsTrim>> loops(1);
    loops[0].push_back(circle_trim(r, z, 1.0, 2 * kPi));
    BrepNurbsFace face;
    std::string why;
    ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why)) << why;
    ASSERT_EQ(face.mSeams.size(), 2u);
    EXPECT_TRUE(face.mSeamHasCurve[0]);
    EXPECT_FALSE(face.mSeamHasCurve[1]);
    const double dev = trim_deviation(face, [](const PatchVec& rP) {
        return std::fabs(std::sqrt(rP[0] * rP[0] + rP[1] * rP[1] + rP[2] * rP[2]) - 1.0);
    });
    EXPECT_LT(dev, 1e-12);
    // A cap of height 1 - z has area 2 pi R h.
    EXPECT_NEAR(sampled_area(one_face(face)), 2 * kPi * (1.0 - z), 0.01 * 2 * kPi * (1.0 - z));
}

TEST(BrepNurbs, ConeFrustumAndTorusBands) {
    {
        BrepNurbsSurface s;
        s.mKind = 4;
        s.mTanHalfAngle = 0.5;
        std::vector<std::vector<BrepNurbsTrim>> loops(2);
        loops[0].push_back(circle_trim(0.5, 1.0, 0.0, 2 * kPi));
        loops[1].push_back(circle_trim(1.0, 2.0, 3.0, -2 * kPi));
        BrepNurbsFace face;
        std::string why;
        ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why))
            << why;
        const double dev = trim_deviation(face, [](const PatchVec& rP) {
            return std::fabs(std::hypot(rP[0], rP[1]) - 0.5 * rP[2]);
        });
        EXPECT_LT(dev, 1e-12);
        const double slant = std::sqrt(1.0 + 0.25);
        EXPECT_NEAR(sampled_area(one_face(face)), kPi * (0.5 + 1.0) * slant, 0.02);
    }
    {
        // A torus band between two parallels of its tube, round the axis.
        BrepNurbsSurface s;
        s.mKind = 5;
        s.mRadius = 3.0;
        s.mMinor = 1.0;
        const auto parallel = [](double Phi, double A0, double Sweep) {
            BrepNurbsTrim t;
            t.mAt = [=](double T) {
                const double a = A0 + Sweep * T, rho = 3.0 + std::cos(Phi);
                return PatchVec{rho * std::cos(a), rho * std::sin(a), std::sin(Phi)};
            };
            return t;
        };
        std::vector<std::vector<BrepNurbsTrim>> loops(2);
        loops[0].push_back(parallel(-0.5, 0.0, 2 * kPi));
        loops[1].push_back(parallel(0.5, 1.0, -2 * kPi));
        BrepNurbsFace face;
        std::string why;
        ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why))
            << why;
        const double dev = trim_deviation(face, [](const PatchVec& rP) {
            const double rho = std::hypot(rP[0], rP[1]) - 3.0;
            return std::fabs(std::hypot(rho, rP[2]) - 1.0);
        });
        EXPECT_LT(dev, 1e-12);
        // The outer band: area 2 pi R r (2 sin 0.5) + 2 pi r^2 ... by integration.
        const double area = 2 * kPi * 1.0 * (3.0 * 1.0 + 2.0 * std::sin(0.5));
        EXPECT_NEAR(sampled_area(one_face(face)), area, 0.01 * area);
    }
}

TEST(BrepNurbs, PlaneWithAHoleKeepsItsCircleExact) {
    BrepNurbsSurface s;
    s.mKind = 1;
    std::vector<std::vector<BrepNurbsTrim>> loops(2);
    const PatchVec c[4] = {{-2, -2, 0}, {2, -2, 0}, {2, 2, 0}, {-2, 2, 0}};
    for (int i = 0; i < 4; ++i) {
        BrepNurbsTrim t;
        t.mPoints = {c[i], c[(i + 1) % 4]};
        t.mEdge = i;
        loops[0].push_back(t);
    }
    // The hole runs clockwise.
    BrepNurbsTrim hole = circle_trim(1.0, 0.0, 0.0, -2 * kPi);
    hole.mCircle = true;
    hole.mRadius = 1.0;
    loops[1].push_back(hole);
    BrepNurbsFace face;
    std::string why;
    ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why)) << why;
    ASSERT_EQ(face.mFace.mLoops.size(), 2u);
    EXPECT_EQ(face.mFace.mLoops[0].mType, "outer");
    EXPECT_EQ(face.mFace.mLoops[1].mType, "inner");
    const BrepCurve& arc = face.mFace.mLoops[1].mTrims[0].mCurve;
    EXPECT_TRUE(arc.mRational);
    EXPECT_EQ(arc.mDegree, 2);
    for (int i = 0; i <= 20; ++i) {
        const PatchVec p = meshioplusplus::detail::brep_curve_point(
            arc, arc.mKnots.front() + (arc.mKnots.back() - arc.mKnots.front()) * i / 20.0);
        EXPECT_NEAR(std::hypot(p[0], p[1]), 1.0, 1e-12);
    }
    EXPECT_NEAR(sampled_area(one_face(face)), 16.0 - kPi, 0.01 * (16.0 - kPi));
}

TEST(BrepNurbs, PlaneCircleStartedAnywhereStaysOnItsPatch) {
    // Started between the 16 samples' axis extremes, the exact circle still
    // lies within the plane's patch, which clamps its parameters.
    const double a0 = std::sqrt(2.0);
    BrepNurbsSurface s;
    s.mKind = 1;
    std::string why;
    {
        // A disc.
        std::vector<std::vector<BrepNurbsTrim>> loops(1);
        loops[0].push_back(circle_trim(1.0, 0.0, a0, 2 * kPi));
        loops[0][0].mCircle = true;
        loops[0][0].mRadius = 1.0;
        BrepNurbsFace face;
        ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why))
            << why;
        EXPECT_LT(circle_deviation(face, face.mFace.mLoops[0].mTrims[0], 1.0), 1e-12);
    }
    {
        // A square with a clockwise hole.
        std::vector<std::vector<BrepNurbsTrim>> loops(2);
        const PatchVec c[4] = {{-2, -2, 0}, {2, -2, 0}, {2, 2, 0}, {-2, 2, 0}};
        for (int i = 0; i < 4; ++i) {
            BrepNurbsTrim t;
            t.mPoints = {c[i], c[(i + 1) % 4]};
            loops[0].push_back(t);
        }
        loops[1].push_back(circle_trim(1.0, 0.0, a0, -2 * kPi));
        loops[1][0].mCircle = true;
        loops[1][0].mRadius = 1.0;
        BrepNurbsFace face;
        ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why))
            << why;
        ASSERT_EQ(face.mFace.mLoops.size(), 2u);
        EXPECT_LT(circle_deviation(face, face.mFace.mLoops[1].mTrims[0], 1.0), 1e-12);
    }
    {
        // Three quarters of a disc: an arc across three axis extremes, closed
        // by a chord.
        const double sweep = 1.5 * kPi;
        std::vector<std::vector<BrepNurbsTrim>> loops(1);
        loops[0].push_back(circle_trim(2.0, 0.0, a0, sweep));
        loops[0][0].mCircle = true;
        loops[0][0].mRadius = 2.0;
        BrepNurbsTrim chord;
        chord.mPoints = {{2.0 * std::cos(a0 + sweep), 2.0 * std::sin(a0 + sweep), 0.0},
                         {2.0 * std::cos(a0), 2.0 * std::sin(a0), 0.0}};
        loops[0].push_back(chord);
        BrepNurbsFace face;
        ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why))
            << why;
        EXPECT_LT(circle_deviation(face, face.mFace.mLoops[0].mTrims[0], 2.0), 1e-12);
        // A circular segment: R^2 (theta - sin theta) / 2.
        const double area = 2.0 * (sweep - std::sin(sweep));
        EXPECT_NEAR(sampled_area(one_face(face)), area, 0.01 * area);
    }
}

TEST(BrepNurbs, InwardNormalSwapsTheLoops) {
    // The same disc seen with the plane's normal pointing in.
    BrepNurbsSurface s;
    s.mKind = 1;
    s.mZ = {0.0, 0.0, -1.0};
    s.mY = {0.0, -1.0, 0.0};
    std::vector<std::vector<BrepNurbsTrim>> loops(1);
    BrepNurbsTrim t = circle_trim(1.0, 0.0, 0.0, 2 * kPi);
    t.mCircle = true;
    t.mRadius = 1.0;
    loops[0].push_back(t);
    BrepNurbsFace face;
    std::string why;
    ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, false, 1e-4, face, why)) << why;
    EXPECT_TRUE(face.mFace.mSwappedNormal);
    EXPECT_FALSE(face.mRefs[0].mAlong);
    EXPECT_NEAR(sampled_area(one_face(face)), kPi, 0.02);
}

TEST(BrepNurbs, BsplineFaceTrimsByProjection) {
    meshioplusplus::detail::BsplineSurface bs;
    bs.mDegreeU = bs.mDegreeV = 3;
    bs.mNu = bs.mNv = 4;
    bs.mKnotsU = bs.mKnotsV = {0, 0, 0, 0, 1, 1, 1, 1};
    for (std::size_t j = 0; j < 4; ++j)
        for (std::size_t i = 0; i < 4; ++i) {
            const double x = i / 3.0, y = j / 3.0;
            bs.mPoles.push_back(x);
            bs.mPoles.push_back(y);
            bs.mPoles.push_back(0.1 * x * y);
        }
    BrepNurbsSurface s;
    s.mKind = 6;
    s.pBspline = &bs;
    // The boundary of the parameter square, as points on the surface.
    const auto at = [&](double U, double V) {
        meshioplusplus::detail::BsplineEval e;
        meshioplusplus::detail::bspline_eval(bs, U, V, 0, e);
        return PatchVec{e.mS[0], e.mS[1], e.mS[2]};
    };
    std::vector<std::vector<BrepNurbsTrim>> loops(1);
    const double c[5][2] = {{0.1, 0.1}, {0.9, 0.1}, {0.9, 0.9}, {0.1, 0.9}, {0.1, 0.1}};
    for (int i = 0; i < 4; ++i) {
        BrepNurbsTrim t;
        for (int k = 0; k <= 10; ++k) {
            const double s0 = k / 10.0;
            t.mPoints.push_back(
                at(c[i][0] + (c[i + 1][0] - c[i][0]) * s0, c[i][1] + (c[i + 1][1] - c[i][1]) * s0));
        }
        loops[0].push_back(t);
    }
    BrepNurbsFace face;
    std::string why;
    ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why)) << why;
    // Each side is a parameter line.
    for (const auto& trim : face.mFace.mLoops[0].mTrims)
        EXPECT_EQ(trim.mCurve.mControlPoints.Shape()[0], 2u);
    EXPECT_FALSE(face.mFace.mRational);
}

TEST(BrepNurbs, HalfTorusWindsRoundTheTube) {
    // Half a torus, bounded by the tube's circles at angles 0 and pi.
    BrepNurbsSurface s;
    s.mKind = 5;
    s.mRadius = 3.0;
    s.mMinor = 1.0;
    const auto meridian = [](double Theta, double P0, double Sweep) {
        BrepNurbsTrim t;
        t.mAt = [=](double T) {
            const double p = P0 + Sweep * T, rho = 3.0 + std::cos(p);
            return PatchVec{rho * std::cos(Theta), rho * std::sin(Theta), std::sin(p)};
        };
        return t;
    };
    std::vector<std::vector<BrepNurbsTrim>> loops(2);
    loops[0].push_back(meridian(0.0, 0.3, -2 * kPi));
    loops[1].push_back(meridian(kPi, 1.0, 2 * kPi));
    BrepNurbsFace face;
    std::string why;
    ASSERT_TRUE(meshioplusplus::detail::brep_nurbs_face(s, loops, true, 1e-4, face, why)) << why;
    ASSERT_EQ(face.mSeams.size(), 1u);
    const double area = 2 * kPi * kPi * 3.0 * 1.0;
    EXPECT_NEAR(sampled_area(one_face(face)), area, 0.01 * area);
}
