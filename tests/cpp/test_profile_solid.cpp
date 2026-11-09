// profile_solid through a provider (roadmap 7.7.4).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/deviation.hpp"
#include "meshioplusplus/operations/fit_curves.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/profile_solid.hpp"

using namespace meshioplusplus;

namespace {

// A closed polyline through rPts (xyz triples) as a chained line mesh, then
// its fit_curves sketch.
Mesh ps_sketch(const std::vector<double>& rCorners) {
    // Four points an edge, so no circle passes through them all.
    std::vector<double> rPts;
    const std::size_t nc = rCorners.size() / 3;
    for (std::size_t i = 0; i < nc; ++i)
        for (int s = 0; s < 4; ++s)
            for (std::size_t k = 0; k < 3; ++k)
                rPts.push_back(rCorners[3 * i + k] +
                               (rCorners[3 * ((i + 1) % nc) + k] - rCorners[3 * i + k]) * s / 4.0);
    const std::size_t n = rPts.size() / 3;
    Mesh m;
    NDArray p = NDArray::Uninit(DType::Float64, {n, 3});
    std::copy(rPts.begin(), rPts.end(), p.As<double>());
    m.AssignPoints(std::move(p));
    NDArray c = NDArray::Uninit(DType::Int64, {n, 2});
    NDArray chain = NDArray::Uninit(DType::Int64, {n});
    for (std::size_t i = 0; i < n; ++i) {
        c.As<std::int64_t>()[2 * i] = static_cast<std::int64_t>(i);
        c.As<std::int64_t>()[2 * i + 1] = static_cast<std::int64_t>((i + 1) % n);
        chain.As<std::int64_t>()[i] = 0;
    }
    m.AddCellBlock("line", std::move(c));
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(chain));
    m.AddCellData("section:chain", std::move(blocks));
    FitCurvesOptions o;
    o.mTolerance = 1e-6;
    return fit_curves(m, o).mMesh;
}

// A 2 x 1 rectangle in z = 0, wound clockwise (the core rewinds it).
Mesh ps_rectangle() {
    return ps_sketch({0, 0, 0, 0, 1, 0, 2, 1, 0, 2, 0, 0});
}

// The half-section of a stepped shaft in y = 0, against the z axis.
Mesh ps_shaft_section() {
    return ps_sketch({0, 0, 0, 1, 0, 0, 1, 0, 2, 0.6, 0, 2, 0.6, 0, 3, 0, 0, 3});
}

}  // namespace

TEST(ProfileSolid, TheRequestAProviderSees) {
    ProfileSolidRequest seen;
    add_profile_solid_provider({"test", [&](const ProfileSolidRequest& rReq) {
                                    seen = rReq;
                                    ProfileSolidResult r;
                                    r.mValid = true;
                                    return r;
                                }});
    ProfileSolidOptions o;
    o.mProvider = "test";
    o.mFeatures = profile_features_from_text("extrude:0;fillet:0.1");
    // An extrusion of 0 runs to the scan's extent along the normal.
    const Mesh scan = icosa_sphere(1.0, 1, {1.0, 0.5, 2.0});
    const ProfileSolidResult r = profile_solid(ps_rectangle(), "", o, &scan);
    EXPECT_EQ(r.mProvider, "test");
    ASSERT_EQ(seen.mLoops.size(), 1u);
    EXPECT_EQ(seen.mLoops[0].size(), 4u);
    for (std::size_t k = 0; k < 4; ++k) {
        const auto& r_a = seen.mLoops[0][k].mEnd;
        const auto& r_b = seen.mLoops[0][(k + 1) % 4].mStart;
        EXPECT_NEAR(std::hypot(r_a[0] - r_b[0], r_a[1] - r_b[1], r_a[2] - r_b[2]), 0.0, 1e-9);
    }
    // Counter-clockwise about the normal.
    double area = 0.0;
    for (const ProfilePiece& r_p : seen.mLoops[0])
        area +=
            0.5 * ((r_p.mStart[0] * r_p.mEnd[1] - r_p.mEnd[0] * r_p.mStart[1]) * seen.mNormal[2] +
                   (r_p.mStart[1] * r_p.mEnd[2] - r_p.mEnd[1] * r_p.mStart[2]) * seen.mNormal[0] +
                   (r_p.mStart[2] * r_p.mEnd[0] - r_p.mEnd[2] * r_p.mStart[0]) * seen.mNormal[1]);
    EXPECT_NEAR(area, 2.0, 1e-9);
    EXPECT_NEAR(std::fabs(seen.mFeatures[0].mDistance), 3.0, 1e-9);
    EXPECT_EQ(seen.mFeatures[1].mKind, ProfileFeatureKind::Fillet);
    EXPECT_EQ(seen.mFeatures[1].mRadius, 0.1);
    remove_profile_solid_provider("test");
}

TEST(ProfileSolid, TheTextFormAndTheRefusals) {
    const std::vector<ProfileFeature> f =
        profile_features_from_text("revolve:0,0,0:0,0,1:180;chamfer:0.2;shell:0.1");
    ASSERT_EQ(f.size(), 3u);
    EXPECT_EQ(f[0].mAngleDeg, 180.0);
    EXPECT_EQ(f[0].mAxisDirection[2], 1.0);
    EXPECT_EQ(f[2].mKind, ProfileFeatureKind::Shell);
    EXPECT_TRUE(profile_features_from_text("extrude:2:symmetric")[0].mSymmetric);
    EXPECT_THROW(profile_features_from_text("extrude:x"), std::invalid_argument);
    EXPECT_THROW(profile_features_from_text("loft:1"), std::invalid_argument);
    ProfileSolidOptions o;
    o.mFeatures = profile_features_from_text("fillet:0.1");
    EXPECT_THROW(profile_solid(ps_rectangle(), "", o), std::invalid_argument);
    o.mFeatures = profile_features_from_text("extrude:1");
    EXPECT_THROW(profile_solid(icosa_sphere(1.0, 1), "", o), std::invalid_argument);
    o.mProvider = "nowhere";
    EXPECT_THROW(profile_solid(ps_rectangle(), "", o), Unsupported);
}

TEST(ProfileSolid, OpenCascadeBuildsTheSolids) {
    if (!profile_solid_provider("occt"))
        GTEST_SKIP() << "no OpenCASCADE profile-solid provider in this build";
    ProfileSolidOptions o;
    o.mFeatures = profile_features_from_text("extrude:3");
    const ProfileSolidResult box = profile_solid(ps_rectangle(), "", o);
    EXPECT_TRUE(box.mValid);
    EXPECT_EQ(box.mNumFaces, 6);
    EXPECT_NEAR(box.mVolume, 6.0, 1e-9);
    o.mFeatures = profile_features_from_text("extrude:3;fillet:0.1");
    const ProfileSolidResult round = profile_solid(ps_rectangle(), "", o);
    EXPECT_TRUE(round.mValid);
    EXPECT_LT(round.mVolume, 6.0);
    EXPECT_GT(round.mNumFaces, 6);
    o.mFeatures = profile_features_from_text("extrude:3;shell:0.1");
    const ProfileSolidResult cup = profile_solid(ps_rectangle(), "", o);
    EXPECT_TRUE(cup.mValid);
    // The walls and the bottom, the top face removed.
    EXPECT_NEAR(cup.mVolume, 6.0 - 1.8 * 0.8 * 2.9, 1e-6);
    // Two sections of a turned shaft revolved about their axis.
    o.mFeatures = profile_features_from_text("revolve:0,0,0:0,0,1");
    const std::string step = ::testing::TempDir() + "shaft.step";
    const ProfileSolidResult shaft = profile_solid(ps_shaft_section(), step, o);
    EXPECT_TRUE(shaft.mValid);
    const double pi = std::acos(-1.0);
    EXPECT_NEAR(shaft.mVolume, pi * (2.0 + 0.36), 1e-6);
    // A scan of the shaft (points on its two turned faces) lies within 0.05
    // of the solid on 95% of it.
    std::vector<double> xyz;
    for (int k = 0; k < 64; ++k)
        for (int j = 0; j < 9; ++j) {
            const double t = 2.0 * pi * k / 64.0;
            const bool upper = j == 8;
            const double r = upper ? 0.6 : 1.0;
            const double z = upper ? 2.5 : 0.1 + 0.225 * j;
            xyz.insert(xyz.end(), {r * std::cos(t), r * std::sin(t), z});
        }
    Mesh scan;
    NDArray sp = NDArray::Uninit(DType::Float64, {xyz.size() / 3, 3});
    std::copy(xyz.begin(), xyz.end(), sp.As<double>());
    scan.AssignPoints(std::move(sp));
    NDArray sv = NDArray::Uninit(DType::Int64, {xyz.size() / 3, 1});
    for (std::size_t i = 0; i < xyz.size() / 3; ++i)
        sv.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    scan.AddCellBlock("vertex", std::move(sv));
    DeviationOptions d;
    d.mTolerance = 0.05;
    const DeviationResult dev = deviation(scan, shaft.mMesh, d);
    EXPECT_GE(dev.mShareWithin, 0.95);
}
