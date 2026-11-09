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
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/deviation.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// The closed, outward-oriented quad skin of the box [rLo, rHi].
Mesh dev_box_skin(const std::array<double, 3>& rLo, const std::array<double, 3>& rHi) {
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> quads;
    for (int a = 0; a < 3; ++a)
        for (int side = 0; side < 2; ++side) {
            const int b = (a + 1) % 3, c = (a + 2) % 3;
            const double uv[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
            std::vector<std::int64_t> q;
            for (int k = 0; k < 4; ++k) {
                const int kk = side == 1 ? k : 3 - k;
                std::vector<double> p(3);
                p[static_cast<std::size_t>(a)] = side == 1 ? rHi[a] : rLo[a];
                p[static_cast<std::size_t>(b)] = uv[kk][0] > 0 ? rHi[b] : rLo[b];
                p[static_cast<std::size_t>(c)] = uv[kk][1] > 0 ? rHi[c] : rLo[c];
                q.push_back(static_cast<std::int64_t>(pts.size()));
                pts.push_back(p);
            }
            quads.push_back(q);
        }
    return mt::make_mesh(pts, "quad", quads);
}

Mesh dev_scaled(const Mesh& rMesh, double S) {
    Mesh out = rMesh;
    NDArray pts = NDArray::Uninit(DType::Float64, rMesh.Points().Shape());
    for (std::size_t i = 0; i < pts.Size(); ++i)
        pts.As<double>()[i] = S * rMesh.Points().As<double>()[i];
    out.AssignPoints(std::move(pts));
    return out;
}

const double* dev_values(const DeviationResult& rR) {
    return rR.mMesh.PointData(kDeviationDistanceName).As<double>();
}

}  // namespace

TEST(Deviation, AnOffsetSphereReadsItsOffset) {
    const Mesh ref = icosa_sphere(1.0, 3);
    const Mesh scan = dev_scaled(ref, 1.01);
    DeviationOptions o;
    o.mTolerance = 0.02;
    const DeviationResult r = deviation(scan, ref, o);
    EXPECT_EQ(r.mNumMeasured, static_cast<std::int64_t>(scan.NumPoints()));
    EXPECT_NEAR(r.mMean, 0.01, 1e-9);
    EXPECT_NEAR(r.mMin, 0.01, 1e-9);
    EXPECT_NEAR(r.mMax, 0.01, 1e-9);
    EXPECT_NEAR(r.mSigma, 0.0, 1e-9);
    EXPECT_DOUBLE_EQ(r.mShareWithin, 1.0);
    EXPECT_NEAR(r.mSuggestedRange, 0.01, 1e-12);
    // A shrunken copy is inside: negative (its vertices are nearest the
    // faces, not the vertices, so a little less than the offset).
    const DeviationResult in = deviation(dev_scaled(ref, 0.99), ref, o);
    EXPECT_LT(in.mMax, 0.0);
    EXPECT_NEAR(in.mMean, -0.01, 1e-4);
}

TEST(Deviation, TheFarSideRuleKeepsAThinWallsFrontFace) {
    // A 0.04 thick plate, scanned 0.03 too high: the scan's bottom face lies
    // nearer the reference's top face (0.01) than its own bottom face (0.03).
    const Mesh ref = dev_box_skin({0, 0, 0}, {1, 1, 0.04});
    const Mesh scan = dev_box_skin({0, 0, 0.03}, {1, 1, 0.07});
    const DeviationResult with = deviation(scan, ref);
    const DeviationResult without = deviation(scan, ref, [] {
        DeviationOptions o;
        o.mMaxNormalAngleDeg = 180.0;
        return o;
    }());
    const double* dw = dev_values(with);
    const double* dn = dev_values(without);
    double worst_with = 0.0, worst_without = 0.0;
    // dev_box_skin does not weld: points 16 to 19 are the bottom face's.
    for (std::size_t p = 16; p < 20; ++p) {
        worst_with = std::max(worst_with, std::fabs(std::fabs(dw[p]) - 0.03));
        worst_without = std::max(worst_without, std::fabs(std::fabs(dn[p]) - 0.03));
    }
    EXPECT_LT(worst_with, 1e-12);
    EXPECT_GT(worst_without, 0.01);
}

TEST(Deviation, AgainstAFittedCylinder) {
    const Mesh scan = cylinder(2.1, 1.0, 16, 4);
    std::array<double, 10> params{{0, 0, 0, 0, 0, 1, 2.0, 0, 0, 0}};
    const DeviationResult r = deviation_to_primitive(scan, 3, params);
    // Every point off the axis reads 0.1; the cylinder's cap centres sit on it.
    const double* d = dev_values(r);
    std::size_t on_wall = 0;
    for (std::size_t p = 0; p < scan.NumPoints(); ++p) {
        const double x = scan.Points().As<double>()[p * 3],
                     y = scan.Points().As<double>()[p * 3 + 1];
        if (std::fabs(std::hypot(x, y) - 2.1) < 1e-9) {
            EXPECT_NEAR(d[p], 0.1, 1e-12);
            ++on_wall;
        }
    }
    EXPECT_GT(on_wall, 0u);
    EXPECT_THROW(deviation_to_primitive(scan, 7, params), std::invalid_argument);
}

TEST(Deviation, HistogramAndPercentiles) {
    const Mesh ref = dev_box_skin({0, 0, 0}, {1, 1, 1});
    // Points at heights 1 + k / 10 above the top face: deviations 0.0 .. 0.9.
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> verts;
    for (int k = 0; k < 10; ++k) {
        pts.push_back({0.5, 0.5, 1.0 + 0.1 * k});
        verts.push_back({k});
    }
    const Mesh cloud = mt::make_mesh(pts, "vertex", verts);
    DeviationOptions o;
    o.mBins = 4;
    o.mRangeMin = 0.0;
    o.mRangeMax = 0.8;
    const DeviationResult r = deviation(cloud, ref, o);
    EXPECT_EQ(r.mNumMeasured, 10);
    EXPECT_EQ(r.mNumOutsideRange, 1);  // 0.9
    const std::int64_t total =
        std::accumulate(r.mHistogramCounts.begin(), r.mHistogramCounts.end(), std::int64_t{0});
    EXPECT_EQ(total, 9);
    ASSERT_EQ(r.mHistogramEdges.size(), 5u);
    EXPECT_DOUBLE_EQ(r.mHistogramEdges.back(), 0.8);
    EXPECT_NEAR(r.mPercentiles[2], 0.45, 1e-12);
    EXPECT_TRUE(std::isnan(r.mShareWithin));
    EXPECT_GT(r.mNumWithoutNormal, 0);
}

TEST(Deviation, BadOptionsAreRefused) {
    const Mesh ref = dev_box_skin({0, 0, 0}, {1, 1, 1});
    DeviationOptions o;
    o.mBins = 0;
    EXPECT_THROW(deviation(ref, ref, o), std::invalid_argument);
    o = DeviationOptions{};
    o.mMaxNormalAngleDeg = 200.0;
    EXPECT_THROW(deviation(ref, ref, o), std::invalid_argument);
}

namespace {

// A 10 x 10 grid of unit quads in z = 0 with `deviation:distance` +2 on a
// 4 x 4 block of points (a bump) and -1 on another (a dent).
Mesh dev_bumped_plate() {
    std::vector<std::vector<double>> pts;
    std::vector<double> d;
    for (int j = 0; j <= 10; ++j)
        for (int i = 0; i <= 10; ++i) {
            pts.push_back({double(i), double(j), 0.0});
            const bool bump = i >= 1 && i <= 4 && j >= 1 && j <= 4;
            const bool dent = i >= 6 && i <= 9 && j >= 6 && j <= 9;
            d.push_back(bump ? 2.0 : (dent ? -1.0 : 0.0));
        }
    std::vector<std::vector<std::int64_t>> quads;
    for (int j = 0; j < 10; ++j)
        for (int i = 0; i < 10; ++i)
            quads.push_back({j * 11 + i, j * 11 + i + 1, (j + 1) * 11 + i + 1, (j + 1) * 11 + i});
    Mesh m = mt::make_mesh(pts, "quad", quads);
    NDArray a = NDArray::Uninit(DType::Float64, {d.size()});
    std::copy(d.begin(), d.end(), a.As<double>());
    m.AddPointData(kDeviationDistanceName, std::move(a));
    return m;
}

}  // namespace

TEST(DeviationRegions, ABumpAndADentAreTwoRegions) {
    DeviationRegionsOptions o;
    o.mTolerance = 0.5;
    const DeviationRegionsResult r = deviation_regions(dev_bumped_plate(), o);
    ASSERT_EQ(r.mRegions.size(), 2u);
    EXPECT_EQ(r.mNumExcess, 1);
    EXPECT_EQ(r.mNumMissing, 1);
    // The bump holds twice the deviation over the same footprint: first.
    EXPECT_EQ(r.mRegions[0].mName, "deviation:excess_0");
    EXPECT_EQ(r.mRegions[0].mSign, 1);
    EXPECT_DOUBLE_EQ(r.mRegions[0].mExtreme, 2.0);
    EXPECT_EQ(r.mRegions[1].mName, "deviation:missing_0");
    EXPECT_EQ(r.mRegions[1].mSign, -1);
    EXPECT_DOUBLE_EQ(r.mRegions[1].mExtreme, -1.0);
    EXPECT_GT(r.mRegions[0].mVolume, 0.0);
    EXPECT_LT(r.mRegions[1].mVolume, 0.0);
    // Its 3 x 3 inner cells read 2; the ring around reads 1 or 1.5 (above 0.5).
    EXPECT_EQ(r.mRegions[0].mNumCells, 21);
    EXPECT_NEAR(r.mRegions[0].mCentroid[0], 2.5, 1e-12);
    EXPECT_NEAR(r.mRegions[0].mCentroid[1], 2.5, 1e-12);
    EXPECT_TRUE(r.mMesh.HasRegion("deviation:excess_0", RegionKind::Cell));
    EXPECT_TRUE(r.mMesh.HasCellData(kDeviationRegionName));
    // A region smaller than the minimum area is dropped.
    o.mMinArea = 1000.0;
    const DeviationRegionsResult none = deviation_regions(dev_bumped_plate(), o);
    EXPECT_TRUE(none.mRegions.empty());
    EXPECT_EQ(none.mNumDropped, 2);
    o = DeviationRegionsOptions{};
    EXPECT_THROW(deviation_regions(dev_bumped_plate(), o), std::invalid_argument);  // tolerance 0
    o.mTolerance = 0.5;
    o.mField = "nope";
    EXPECT_THROW(deviation_regions(dev_bumped_plate(), o), std::invalid_argument);
}

TEST(Exaggerate, AnOffsetSphereGrowsByTheScaledOffset) {
    const Mesh ref = icosa_sphere(1.0, 3);
    const DeviationResult d = deviation(dev_scaled(ref, 1.01), ref);
    ExaggerateOptions o;
    o.mScale = 20.0;
    const ExaggerateResult r = exaggerate(d.mMesh, o);
    EXPECT_NEAR(r.mClamp, 0.01, 1e-12);
    EXPECT_EQ(r.mNumMoved, static_cast<std::int64_t>(ref.NumPoints()));
    const double* p = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        const double rad = std::sqrt(p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 1] +
                                     p[3 * i + 2] * p[3 * i + 2]);
        EXPECT_NEAR(rad, 1.01 + 20.0 * 0.01, 2e-3);
    }
    // The measured scan is untouched.
    EXPECT_NEAR(d.mMesh.Points().As<double>()[0] * d.mMesh.Points().As<double>()[0] +
                    d.mMesh.Points().As<double>()[1] * d.mMesh.Points().As<double>()[1] +
                    d.mMesh.Points().As<double>()[2] * d.mMesh.Points().As<double>()[2],
                1.01 * 1.01, 1e-12);
    o.mScale = std::nan("");
    EXPECT_THROW(exaggerate(d.mMesh, o), std::invalid_argument);
}
