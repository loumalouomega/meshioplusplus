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
// The OpenCASCADE-free half of fit_cad's fillets and threads (roadmap 6.3.1).

#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/feature_fit.hpp"

using namespace meshioplusplus::detail;

namespace {

constexpr double kPi = 3.14159265358979323846;

// (I + 1) x (J + 1) points through `rAt(u, v)`, u and v in [0, 1], two
// triangles per cell wound so the normal is d/du x d/dv.
template <class F>
void grid(int I, int J, F&& rAt, std::vector<PatchVec>& rPts, std::vector<PatchTri>& rTris) {
    rPts.clear();
    rTris.clear();
    for (int j = 0; j <= J; ++j)
        for (int i = 0; i <= I; ++i)
            rPts.push_back(rAt(static_cast<double>(i) / I, static_cast<double>(j) / J));
    auto id = [&](int i, int j) { return static_cast<std::size_t>(j * (I + 1) + i); };
    for (int j = 0; j < J; ++j)
        for (int i = 0; i < I; ++i) {
            rTris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            rTris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
}

}  // namespace

TEST(FeatureFit, PrimitiveDistances) {
    const double plane[10] = {0, 0, 1, 0, 0, 1, 0, 0, 0, 0};
    EXPECT_NEAR(primitive_distance(1, plane, {3, 4, -1}), 2.0, 1e-12);
    const double sphere[10] = {1, 0, 0, 0, 0, 0, 2, 0, 0, 0};
    EXPECT_NEAR(primitive_distance(2, sphere, {1, 0, 3}), 1.0, 1e-12);
    const double cyl[10] = {0, 0, 0, 0, 0, 1, 1, 0, 0, 0};
    EXPECT_NEAR(primitive_distance(3, cyl, {0, 3, 7}), 2.0, 1e-12);
    const double cone[10] = {0, 0, 0, 0, 0, 1, 0, 0, 45, 0};
    EXPECT_NEAR(primitive_distance(4, cone, {1, 0, 1}), 0.0, 1e-12);
    EXPECT_NEAR(primitive_distance(4, cone, {0, 0, -2}), 2.0, 1e-12);  // behind the apex
    const double torus[10] = {0, 0, 0, 0, 0, 1, 3, 1, 0, 0};
    EXPECT_NEAR(primitive_distance(5, torus, {3, 0, 0.5}), 0.5, 1e-12);
}

TEST(FeatureFit, ConvexRoundedEdgeBetweenTwoPlanes) {
    // The edge of x <= 1, y <= 1 rounded with r = 0.2: centres inside.
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    grid(
        12, 8,
        [](double U, double V) {
            const double a = 0.5 * kPi * U;
            return PatchVec{0.8 + 0.2 * std::cos(a), 0.8 + 0.2 * std::sin(a), 2.0 * V};
        },
        pts, tris);
    const double px[10] = {1, 0, 0, 1, 0, 0, 0, 0, 0, 0};
    const double py[10] = {0, 1, 0, 0, 1, 0, 0, 0, 0, 0};
    const FilletFit f = fillet_fit(pts, tris, 1, px, 1, py, 1e-3, 10.0);
    ASSERT_TRUE(f.mOk);
    EXPECT_EQ(f.mSide, 1);
    EXPECT_NEAR(f.mRadius, 0.2, 2e-3);
    EXPECT_LT(f.mRms, 2e-3);
    for (const PatchVec& c : f.mCentres) {
        EXPECT_NEAR(c[0], 0.8, 5e-3);
        EXPECT_NEAR(c[1], 0.8, 5e-3);
    }
}

TEST(FeatureFit, ConcaveFilletBetweenAPlaneAndACylinder) {
    // A boss of radius 1 on the plane z = 0, its foot rounded with r = 0.3:
    // the ball rolls outside the material, centres on the ring rho = 1.3, z = 0.3.
    const double r = 0.3;
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    grid(
        48, 10,
        [r](double U, double V) {
            const double a = 2.0 * kPi * U * 0.75, s = 0.5 * kPi * V;
            // From the plane contact (straight below the centre) to the cylinder
            // contact (towards the axis), around the corner side.
            const double rho = 1.0 + r - r * std::sin(s), z = r - r * std::cos(s);
            return PatchVec{rho * std::cos(a), rho * std::sin(a), z};
        },
        pts, tris);
    const double plane[10] = {0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
    const double cyl[10] = {0, 0, 0, 0, 0, 1, 1, 0, 0, 0};
    const FilletFit f = fillet_fit(pts, tris, 1, plane, 3, cyl, 1e-3, 10.0);
    ASSERT_TRUE(f.mOk);
    EXPECT_EQ(f.mSide, -1);
    EXPECT_NEAR(f.mRadius, r, 3e-3);
}

TEST(FeatureFit, HelicoidRecoversPitchSlopeAndOffset) {
    const double k = 0.25 / (2.0 * kPi), t = std::tan(30.0 * kPi / 180.0), c = 0.1;
    std::mt19937 rng(7);
    std::normal_distribution<double> noise(0.0, 1e-4);
    for (double sigma : {0.0, 1.0}) {
        std::vector<PatchVec> pts;
        std::vector<PatchTri> tris;
        grid(
            4, 96,
            [&](double U, double V) {
                const double rho = 0.8 + 0.2 * U, th = 4.0 * kPi * V;
                return PatchVec{rho * std::cos(th) + sigma * noise(rng),
                                rho * std::sin(th) + sigma * noise(rng),
                                c + k * th + t * rho + sigma * noise(rng)};
            },
            pts, tris);
        const HelicoidFit h = helicoid_fit(pts, tris, {0, 0, 0}, {0, 0, 1});
        ASSERT_TRUE(h.mOk);
        EXPECT_NEAR(h.Pitch(), 0.25, 1e-3);
        EXPECT_NEAR(h.mT, t, 5e-3);
        EXPECT_NEAR(h.mThetaMax - h.mThetaMin, 4.0 * kPi, sigma > 0 ? 1e-2 : 1e-9);
        EXPECT_LT(h.mRms, sigma > 0 ? 5e-4 : 1e-12);
        // The surface reproduces a vertex.
        const PatchVec q = h.At(0.9, 1.0);
        EXPECT_NEAR(q[2], c + k * 1.0 + t * 0.9, sigma > 0 ? 1e-3 : 1e-12);
    }
}

TEST(FeatureFit, HelicoidNeedsHalfATurnAndAnAxisOffThePatch) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    grid(
        4, 16,
        [](double U, double V) {
            const double rho = 0.8 + 0.2 * U, th = 0.8 * kPi * V;
            return PatchVec{rho * std::cos(th), rho * std::sin(th), 0.05 * th};
        },
        pts, tris);
    EXPECT_FALSE(helicoid_fit(pts, tris, {0, 0, 0}, {0, 0, 1}).mOk);
    // An axis through the patch.
    EXPECT_FALSE(helicoid_fit(pts, tris, pts[5], {0, 0, 1}).mOk);
}
