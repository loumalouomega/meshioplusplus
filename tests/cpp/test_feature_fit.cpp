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

namespace {

// The convex edge of x <= 1, y <= 1 along z in [0, 2], rounded by a ball whose
// radius r(z) = 0.25 + 0.1 cos(pi z / 2) varies along it: the exact envelope
// of the balls, each row a characteristic circle from the foot on x = 1 to
// the foot on y = 1.
double variable_radius(double Z) {
    return 0.25 + 0.1 * std::cos(0.5 * kPi * Z);
}

PatchVec variable_edge(double U, double Z) {
    const double r = variable_radius(Z), dr = -0.05 * kPi * std::sin(0.5 * kPi * Z);
    const PatchVec c{1.0 - r, 1.0 - r, Z};
    const double len = std::sqrt(1.0 + 2.0 * dr * dr);
    const PatchVec t{-dr / len, -dr / len, 1.0 / len};
    const double delta = -r * dr / len, rho = std::sqrt(r * r - delta * delta);
    const PatchVec cc{c[0] + delta * t[0], c[1] + delta * t[1], c[2] + delta * t[2]};
    const PatchVec f1{1.0, 1.0 - r, Z}, f2{1.0 - r, 1.0, Z};
    PatchVec a{}, b{};
    for (int k = 0; k < 3; ++k)
        a[k] = (f1[k] - cc[k]) / rho;
    double ab = 0.0;
    for (int k = 0; k < 3; ++k)
        ab += (f2[k] - cc[k]) * a[k];
    double bl = 0.0;
    for (int k = 0; k < 3; ++k) {
        b[k] = f2[k] - cc[k] - ab * a[k];
        bl += b[k] * b[k];
    }
    bl = std::sqrt(bl);
    const double w = std::atan2(bl, ab) * U;
    PatchVec p{};
    for (int k = 0; k < 3; ++k)
        p[k] = cc[k] + rho * (std::cos(w) * a[k] + std::sin(w) * b[k] / bl);
    return p;
}

bool foot_x1(const PatchVec& rX, PatchVec& rF) {
    rF = {1.0, rX[1], rX[2]};
    return true;
}

bool foot_y1(const PatchVec& rX, PatchVec& rF) {
    rF = {rX[0], 1.0, rX[2]};
    return true;
}

// A helical surface about the axis through rO along the unit rA: the crest of
// a rounded thread, rho = 0.9 + sqrt(0.1^2 - zeta^2), swept by a screw of
// pitch 0.5 over two turns.
PatchVec rounded_crest(double U, double V, const PatchVec& rO, const PatchVec& rA) {
    PatchVec e1{}, e2{};
    {
        PatchVec h{1.0, 0.0, 0.0};
        double d = h[0] * rA[0];
        double l = 0.0;
        for (int k = 0; k < 3; ++k) {
            e1[k] = h[k] - d * rA[k];
            l += e1[k] * e1[k];
        }
        for (int k = 0; k < 3; ++k)
            e1[k] /= std::sqrt(l);
        e2 = {rA[1] * e1[2] - rA[2] * e1[1], rA[2] * e1[0] - rA[0] * e1[2],
              rA[0] * e1[1] - rA[1] * e1[0]};
    }
    const double zeta = -0.07 + 0.14 * U, th = 4.0 * kPi * V;
    const double rho = 0.9 + std::sqrt(0.01 - zeta * zeta), z = zeta + 0.5 / (2.0 * kPi) * th;
    PatchVec p{};
    for (int k = 0; k < 3; ++k)
        p[k] = rO[k] + rho * (std::cos(th) * e1[k] + std::sin(th) * e2[k]) + z * rA[k];
    return p;
}

}  // namespace

TEST(FeatureFit, LawsFitPolynomialsAndPeriodicSeries) {
    std::vector<double> t, y, yp;
    for (int i = 0; i <= 40; ++i) {
        const double x = 3.0 * i / 40.0;
        t.push_back(x);
        y.push_back(1.0 + 0.5 * x - 0.25 * x * x);
        yp.push_back(2.0 + std::cos(2.0 * kPi * x / 3.0) - 0.5 * std::sin(4.0 * kPi * x / 3.0));
    }
    FeatureLaw law;
    ASSERT_TRUE(feature_law_fit(t, y, false, 0.0, 3.0, 2, law));
    EXPECT_NEAR(law.At(1.7), 1.0 + 0.5 * 1.7 - 0.25 * 1.7 * 1.7, 1e-12);
    EXPECT_NEAR(law.Slope(1.7), 0.5 - 0.5 * 1.7, 1e-12);
    ASSERT_TRUE(feature_law_fit(t, yp, true, 0.0, 3.0, 2, law));
    EXPECT_NEAR(law.At(0.4),
                2.0 + std::cos(2.0 * kPi * 0.4 / 3.0) - 0.5 * std::sin(4.0 * kPi * 0.4 / 3.0),
                1e-12);
    EXPECT_NEAR(law.At(0.4 + 3.0), law.At(0.4), 1e-12);
    EXPECT_FALSE(feature_law_fit({0.0, 1.0}, {0.0, 1.0}, false, 0.0, 1.0, 3, law));
}

TEST(FeatureFit, VariableRadiusBlendBetweenTwoPlanes) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    grid(12, 32, variable_edge, pts, tris);
    const RollingBallFit f = rolling_ball_fit(pts, tris, foot_x1, foot_y1, 5e-4, 10.0, 1e-3);
    ASSERT_TRUE(f.mOk);
    EXPECT_EQ(f.mSide, 1);
    EXPECT_TRUE(f.mVariable);
    EXPECT_FALSE(f.mClosed);
    EXPECT_LT(f.mRms, 2.5e-4);
    // Every section is a ball between the planes, of the radius at its height.
    for (int i = 0; i <= 10; ++i) {
        const double t = f.mLength * i / 10.0;
        const RollingBallSection s =
            rolling_ball_section(f, foot_x1, foot_y1, t, f.RailTangent(t), 2.0);
        ASSERT_TRUE(s.mOk) << t;
        EXPECT_NEAR(s.mCentre[0], 1.0 - s.mRadius, 1e-9);
        EXPECT_NEAR(s.mCentre[1], 1.0 - s.mRadius, 1e-9);
        EXPECT_NEAR(s.mRadius, variable_radius(s.mCentre[2]), 1e-3);
    }
    // A constant radius is no variable blend.
    grid(
        12, 16,
        [](double U, double V) {
            const double a = 0.5 * kPi * U;
            return PatchVec{0.8 + 0.2 * std::cos(a), 0.8 + 0.2 * std::sin(a), 2.0 * V};
        },
        pts, tris);
    const RollingBallFit c = rolling_ball_fit(pts, tris, foot_x1, foot_y1, 5e-4, 10.0, 1e-3);
    ASSERT_TRUE(c.mOk);
    EXPECT_FALSE(c.mVariable);
    EXPECT_NEAR(c.mRadius.At(0.5), 0.2, 1e-3);
}

TEST(FeatureFit, NoBlendOnAPatchNoBallFits) {
    // A flat strip between the two planes: no ball touches both from it.
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    grid(
        8, 16,
        [](double U, double V) { return PatchVec{1.0 - 0.3 * U, 0.7 + 0.3 * U * U, 2.0 * V}; }, pts,
        tris);
    EXPECT_FALSE(rolling_ball_fit(pts, tris, foot_x1, foot_y1, 5e-4, 10.0, 1e-3).mOk);
}

TEST(FeatureFit, HelicalFitFindsTheScrewOfARoundedProfile) {
    const double al = std::sqrt(0.01 + 0.04 + 1.0);
    const PatchVec o{0.3, -0.2, 0.1}, a{0.1 / al, 0.2 / al, 1.0 / al};
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    grid(8, 192, [&](double U, double V) { return rounded_crest(U, V, o, a); }, pts, tris);
    const double tol = 1e-3;
    for (bool seeded : {false, true}) {
        // Seeded with a slightly wrong axis, as a fitted cylinder gives one.
        const PatchVec so{o[0] + 1e-3, o[1], o[2]}, sa{a[0] + 1e-3, a[1], a[2]};
        const HelicalFit h = seeded ? helical_fit(pts, tris, &so, &sa, tol)
                                    : helical_fit(pts, tris, nullptr, nullptr, tol);
        ASSERT_TRUE(h.mOk) << seeded;
        EXPECT_FALSE(h.mStraight);
        EXPECT_TRUE(h.mRhoOfZeta);
        EXPECT_NEAR(std::fabs(h.Pitch()), 0.5, 1e-4);
        EXPECT_NEAR(std::fabs(h.mAxis[0] * a[0] + h.mAxis[1] * a[1] + h.mAxis[2] * a[2]), 1.0,
                    1e-9);
        EXPECT_LT(h.mMax, 0.5 * tol);
        EXPECT_NEAR(h.mThetaMax - h.mThetaMin, 4.0 * kPi, 1e-4);
        // The surface passes through the patch.
        const PatchVec p = rounded_crest(0.5, 0.3, o, a);
        double rho = 0.0, z = 0.0;
        {
            PatchVec d{p[0] - h.mOrigin[0], p[1] - h.mOrigin[1], p[2] - h.mOrigin[2]};
            z = d[0] * h.mAxis[0] + d[1] * h.mAxis[1] + d[2] * h.mAxis[2];
            for (int k = 0; k < 3; ++k)
                d[k] -= z * h.mAxis[k];
            rho = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        }
        EXPECT_NEAR(rho, 1.0, 1e-3);
    }
}

TEST(FeatureFit, HelicalFitFindsAStraightFlankWithoutAnAxis) {
    const double k = 0.25 / (2.0 * kPi), t = std::tan(30.0 * kPi / 180.0);
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    grid(
        4, 96,
        [&](double U, double V) {
            const double rho = 0.8 + 0.2 * U, th = 4.0 * kPi * V;
            return PatchVec{rho * std::cos(th), rho * std::sin(th), 0.1 + k * th + t * rho};
        },
        pts, tris);
    const HelicalFit h = helical_fit(pts, tris, nullptr, nullptr, 1e-3);
    ASSERT_TRUE(h.mOk);
    EXPECT_TRUE(h.mStraight);
    EXPECT_NEAR(std::fabs(h.Pitch()), 0.25, 1e-6);
    EXPECT_NEAR(std::fabs(h.mAxis[2]), 1.0, 1e-9);
}

TEST(FeatureFit, HelicalFitRefusesCylindersAndSurfacesOfRevolution) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    // A cylinder: every screw about its axis leaves it unchanged.
    grid(
        48, 8,
        [](double U, double V) {
            const double th = 1.5 * kPi * U;
            return PatchVec{std::cos(th), std::sin(th), 2.0 * V};
        },
        pts, tris);
    EXPECT_FALSE(helical_fit(pts, tris, nullptr, nullptr, 1e-3).mOk);
    // A band of a torus: a rotation, with no advance.
    grid(
        48, 8,
        [](double U, double V) {
            const double th = 1.5 * kPi * U, s = -0.8 + 1.6 * V;
            const double rho = 1.0 + 0.3 * std::cos(s);
            return PatchVec{rho * std::cos(th), rho * std::sin(th), 0.3 * std::sin(s)};
        },
        pts, tris);
    EXPECT_FALSE(helical_fit(pts, tris, nullptr, nullptr, 1e-3).mOk);
}
