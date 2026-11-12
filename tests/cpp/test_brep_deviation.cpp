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

// detail/brep_deviation: the two-sided distance and the normal deviation
// between a trimmed NURBS B-rep and a skin (fit_cad's validation report,
// doc/fit_cad_freeform_plan.md step 2), and the report as fit_cad fills it.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "../../src/cpp/src/detail/brep_deviation.hpp"
#include "analytic_surfaces.hpp"
#include "fit_cad_fixtures.hpp"
#include "thread_sweep.hpp"

using namespace meshioplusplus;

namespace {

// A bilinear patch over [0, S] x [0, S] at height Z (control points
// u-fastest), optionally trimmed to the loop through `rLoop` (u, v) corners.
BrepInfo bd_plane(double S, double Z, const std::vector<std::array<double, 2>>& rLoop = {}) {
    BrepInfo info;
    NurbsPatch p;
    p.mDegrees = {1, 1};
    p.mKnots = {{0.0, 0.0, S, S}, {0.0, 0.0, S, S}};
    p.mControlPoints = NDArray(DType::Float64, {4, 3});
    const double xyz[12] = {0, 0, Z, S, 0, Z, 0, S, Z, S, S, Z};
    for (int i = 0; i < 12; ++i)
        p.mControlPoints.As<double>()[i] = xyz[i];
    p.mWeights = NDArray(DType::Float64, {4});
    for (int i = 0; i < 4; ++i)
        p.mWeights.As<double>()[i] = 1.0;
    info.mPatches.push_back(std::move(p));
    BrepFace f;
    f.mPatch = 0;
    if (!rLoop.empty()) {
        f.mTrimmed = true;
        BrepLoop loop;
        for (std::size_t i = 0; i < rLoop.size(); ++i) {
            const auto& a = rLoop[i];
            const auto& b = rLoop[(i + 1) % rLoop.size()];
            BrepTrim t;
            t.mCurve.mDegree = 1;
            t.mCurve.mKnots = {0.0, 0.0, 1.0, 1.0};
            t.mCurve.mControlPoints = NDArray(DType::Float64, {2, 3});
            const double uv[6] = {a[0], a[1], 0.0, b[0], b[1], 0.0};
            for (int k = 0; k < 6; ++k)
                t.mCurve.mControlPoints.As<double>()[k] = uv[k];
            loop.mTrims.push_back(std::move(t));
        }
        f.mLoops.push_back(std::move(loop));
    }
    info.mFaces.push_back(std::move(f));
    return info;
}

// The unit square at z = 0 as two triangles, optionally its far corner raised.
void bd_square(std::vector<double>& rXyz, std::vector<std::int64_t>& rTri, double Raise = 0.0) {
    rXyz = {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, Raise};
    rTri = {0, 1, 3, 0, 3, 2};
}

}  // namespace

TEST(BrepDeviation, PlaneOnItsSkinIsZero) {
    std::vector<double> xyz;
    std::vector<std::int64_t> tri;
    bd_square(xyz, tri);
    const detail::BrepDeviation d = detail::brep_deviation(bd_plane(1.0, 0.0), xyz, tri);
    ASSERT_TRUE(d.mOk) << d.mWhy;
    EXPECT_LE(d.mCadToMeshMax, 1e-14);
    EXPECT_LE(d.mMeshToCadMax, 1e-14);
    EXPECT_LE(d.mNormalMaxDeg, 1e-6);
    ASSERT_EQ(d.mFaces.size(), 1u);
    EXPECT_EQ(d.mFaces[0].mNumSamples, 25 * 25);
}

TEST(BrepDeviation, ShiftedPlaneGivesTheOffsetBothWays) {
    std::vector<double> xyz;
    std::vector<std::int64_t> tri;
    bd_square(xyz, tri);
    const detail::BrepDeviation d = detail::brep_deviation(bd_plane(1.0, 0.01), xyz, tri);
    ASSERT_TRUE(d.mOk);
    EXPECT_NEAR(d.mCadToMeshMax, 0.01, 1e-12);
    EXPECT_NEAR(d.mCadToMeshP95, 0.01, 1e-12);
    EXPECT_NEAR(d.mMeshToCadMax, 0.01, 1e-12);
    EXPECT_NEAR(d.mMeshToCadP99, 0.01, 1e-12);
    EXPECT_LE(d.mNormalMaxDeg, 1e-6);
}

TEST(BrepDeviation, TrimsBoundTheSamples) {
    std::vector<double> xyz;
    std::vector<std::int64_t> tri;
    bd_square(xyz, tri);
    // The plane runs to (2, 2): untrimmed, its far corner is sqrt(2) from the
    // skin; trimmed to the unit square, nothing of it is.
    const detail::BrepDeviation whole = detail::brep_deviation(bd_plane(2.0, 0.0), xyz, tri);
    ASSERT_TRUE(whole.mOk);
    EXPECT_NEAR(whole.mCadToMeshMax, std::sqrt(2.0), 1e-12);
    const detail::BrepDeviation trimmed =
        detail::brep_deviation(bd_plane(2.0, 0.0, {{0, 0}, {1, 0}, {1, 1}, {0, 1}}), xyz, tri);
    ASSERT_TRUE(trimmed.mOk);
    EXPECT_LE(trimmed.mCadToMeshMax, 1e-14);
    EXPECT_LE(trimmed.mMeshToCadMax, 1e-14);
    EXPECT_LT(trimmed.mFaces[0].mNumSamples, whole.mFaces[0].mNumSamples);
}

TEST(BrepDeviation, NormalDeviationSeesATiltedFacet) {
    std::vector<double> xyz;
    std::vector<std::int64_t> tri;
    bd_square(xyz, tri, 0.2);  // the far corner raised: both facets tilt
    const detail::BrepDeviation d = detail::brep_deviation(bd_plane(1.0, 0.0), xyz, tri);
    ASSERT_TRUE(d.mOk);
    // Facet (0, 1, 3) has normal (0, -0.2, 1): an 11.3 degree tilt.
    EXPECT_NEAR(d.mNormalMaxDeg, std::atan(0.2) * 180.0 / 3.14159265358979323846, 1e-9);
    EXPECT_LE(d.mNormalP95Deg, d.mNormalMaxDeg);
    EXPECT_NEAR(d.mCadToMeshMax, 0.2 / std::sqrt(1.04), 1e-9);
}

// A face with a hole and a second face filling it: a skin point over the
// hole, whose nearest grid triangles belong to the holed face (its cells
// overhang the hole), is measured to the face that holds it.
TEST(BrepDeviation, PointOverAHoleFindsTheFaceFillingIt) {
    std::vector<double> xyz;
    std::vector<std::int64_t> tri;
    bd_square(xyz, tri);
    BrepInfo ring = bd_plane(1.0, 0.0, {{0, 0}, {1, 0}, {1, 1}, {0, 1}});
    {
        BrepLoop hole;
        const std::vector<std::array<double, 2>> c = {
            {0.3, 0.3}, {0.3, 0.7}, {0.7, 0.7}, {0.7, 0.3}};
        for (std::size_t i = 0; i < c.size(); ++i) {
            BrepTrim t;
            t.mCurve.mDegree = 1;
            t.mCurve.mKnots = {0.0, 0.0, 1.0, 1.0};
            t.mCurve.mControlPoints = NDArray(DType::Float64, {2, 3});
            const auto& a = c[i];
            const auto& b = c[(i + 1) % c.size()];
            const double uv[6] = {a[0], a[1], 0.0, b[0], b[1], 0.0};
            for (int k = 0; k < 6; ++k)
                t.mCurve.mControlPoints.As<double>()[k] = uv[k];
            hole.mTrims.push_back(std::move(t));
        }
        hole.mType = "inner";
        ring.mFaces[0].mLoops.push_back(std::move(hole));
    }
    // The filling face: the same plane trimmed to the hole.
    BrepInfo fill = bd_plane(1.0, 0.0, {{0.3, 0.3}, {0.7, 0.3}, {0.7, 0.7}, {0.3, 0.7}});
    ring.mPatches.push_back(fill.mPatches[0]);
    fill.mFaces[0].mPatch = 1;
    ring.mFaces.push_back(fill.mFaces[0]);
    // Skin points over the hole: the centre of each facet, and one more.
    xyz.insert(xyz.end(), {0.5, 0.5, 0.0});
    const std::int64_t c = 4;
    tri = {0, 1, c, 1, 3, c, 3, 2, c, 2, 0, c};
    const detail::BrepDeviation d = detail::brep_deviation(ring, xyz, tri);
    ASSERT_TRUE(d.mOk);
    EXPECT_LE(d.mMeshToCadMax, 1e-12);
    EXPECT_LE(d.mCadToMeshMax, 1e-12);
    EXPECT_LE(d.mNormalMaxDeg, 1e-6);
}

TEST(BrepDeviation, NothingToMeasure) {
    std::vector<double> xyz;
    std::vector<std::int64_t> tri;
    bd_square(xyz, tri);
    EXPECT_FALSE(detail::brep_deviation(BrepInfo{}, xyz, tri).mOk);
    EXPECT_FALSE(detail::brep_deviation(bd_plane(1.0, 0.0), xyz, {}).mOk);
}

// fit_cad's report: off, nothing is measured and the solid is the same; on,
// the native capped cylinder's analytic faces lie within the tolerance of
// the skin both ways (their own chordal sag at most), every face is sampled,
// and the measures do not depend on the thread count.
TEST(BrepDeviation, FitCadValidationReport) {
    const Mesh fitted = fcfix::fitted(analytic::capped_cylinder());
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    const FitCadResult plain = fit_cad(fitted, "", o);
    EXPECT_FALSE(plain.mValidationReported);
    EXPECT_TRUE(std::isnan(plain.mCadToMeshMax));
    EXPECT_FALSE(plain.mMesh.HasFieldData(kFitCadFaceDeviationName));
    o.mValidationReport = true;
    // The one-thread values (a result holds a mesh, which a KRATOS build cannot copy).
    double ref_c2m = 0.0, ref_m2c = 0.0, ref_normal = 0.0;
    for (int threads : {1, 4}) {
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
        mio_test::set_threads(threads);
#endif
        const FitCadResult r = fit_cad(fitted, "", o);
        ASSERT_TRUE(r.mValidationReported);
        EXPECT_EQ(r.mCounts.mNumFaces, plain.mCounts.mNumFaces);
        // The wall's facets sag 1 - cos(pi / 32) from the cylinder.
        const double sag = 1.05 * (1.0 - std::cos(3.14159265358979323846 / 32.0));
        EXPECT_LE(r.mCadToMeshMax, std::max(r.mTolerance, sag));
        EXPECT_LE(r.mMeshToCadMax, std::max(r.mTolerance, sag));
        EXPECT_LE(r.mCadToMeshP95, r.mCadToMeshMax);
        EXPECT_LE(r.mMeshToCadP95, r.mMeshToCadMax);
        EXPECT_LT(r.mNormalDeviationMaxDeg, 15.0);
        EXPECT_FALSE(r.mStepRoundTripChecked);  // no file written
        ASSERT_TRUE(r.mMesh.HasFieldData(kFitCadFaceDeviationName));
        const NDArray& table = r.mMesh.FieldData(kFitCadFaceDeviationName);
        EXPECT_EQ(table.Shape()[0], static_cast<std::size_t>(r.mCounts.mNumFaces));
        EXPECT_EQ(table.Shape()[1], 5u);
        if (threads == 1) {
            ref_c2m = r.mCadToMeshMax;
            ref_m2c = r.mMeshToCadMax;
            ref_normal = r.mNormalDeviationMaxDeg;
        } else {
            EXPECT_EQ(r.mCadToMeshMax, ref_c2m);
            EXPECT_EQ(r.mMeshToCadMax, ref_m2c);
            EXPECT_EQ(r.mNormalDeviationMaxDeg, ref_normal);
        }
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(mio_test::num_procs());
#endif
}
