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
// Rigid registration by iterative closest points (roadmap 7.5.3).

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/operations/registration.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

using RgV = std::array<double, 3>;

// A lumpy closed cloud: no symmetry for ICP to settle into.
std::vector<RgV> rg_cloud(std::size_t N) {
    std::mt19937_64 rng(5);
    std::normal_distribution<double> g;
    std::vector<RgV> out(N);
    for (auto& q : out) {
        double x = g(rng), y = g(rng), z = g(rng);
        const double r = std::sqrt(x * x + y * y + z * z);
        x /= r;
        y /= r;
        z /= r;
        const double s = 1.0 + 0.3 * std::sin(3 * x) * std::cos(2 * y) + 0.2 * z * z;
        q = {1.6 * s * x, 1.0 * s * y, 0.7 * s * z};
    }
    return out;
}

Mesh rg_mesh(const std::vector<RgV>& rPts) {
    Mesh m;
    NDArray p(DType::Float64, {rPts.size(), std::size_t{3}});
    for (std::size_t i = 0; i < rPts.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            p.As<double>()[i * 3 + k] = rPts[i][k];
    m.AssignPoints(std::move(p));
    return m;
}

// x -> R x + t for a rotation about (1, 2, 3) by Angle.
std::array<double, 16> rg_motion(double Angle, RgV T) {
    RgV a{1, 2, 3};
    const double l = std::sqrt(14.0);
    a = {a[0] / l, a[1] / l, a[2] / l};
    const double c = std::cos(Angle), s = std::sin(Angle), C = 1 - c;
    return {c + a[0] * a[0] * C,
            a[0] * a[1] * C - a[2] * s,
            a[0] * a[2] * C + a[1] * s,
            T[0],
            a[1] * a[0] * C + a[2] * s,
            c + a[1] * a[1] * C,
            a[1] * a[2] * C - a[0] * s,
            T[1],
            a[2] * a[0] * C - a[1] * s,
            a[2] * a[1] * C + a[0] * s,
            c + a[2] * a[2] * C,
            T[2],
            0,
            0,
            0,
            1};
}

std::vector<RgV> rg_moved(const std::vector<RgV>& rP, const std::array<double, 16>& rM) {
    std::vector<RgV> out(rP.size());
    for (std::size_t i = 0; i < rP.size(); ++i)
        for (std::size_t r = 0; r < 3; ++r)
            out[i][r] = rM[4 * r] * rP[i][0] + rM[4 * r + 1] * rP[i][1] + rM[4 * r + 2] * rP[i][2] +
                        rM[4 * r + 3];
    return out;
}

}  // namespace

TEST(Registration, RotatedShiftedCopyRegistersBack) {
    const auto target = rg_cloud(3000);
    const auto motion = rg_motion(0.15, {0.2, -0.1, 0.05});
    const auto source = rg_moved(target, motion);  // source = M target
    for (const RegistrationMethod method :
         {RegistrationMethod::PointToPoint, RegistrationMethod::PointToPlane}) {
        RegistrationOptions o;
        o.mMethod = method;
        Mesh tmesh = rg_mesh(target);
        if (method == RegistrationMethod::PointToPlane) {
            // the exact normals are not needed: any normals make a plane per point
            NDArray nrm(DType::Float64, {target.size(), std::size_t{3}});
            for (std::size_t i = 0; i < target.size(); ++i) {
                const double l =
                    std::sqrt(target[i][0] * target[i][0] / 2.56 + target[i][1] * target[i][1] +
                              target[i][2] * target[i][2] / 0.49);
                nrm.As<double>()[3 * i] = target[i][0] / 2.56 / l;
                nrm.As<double>()[3 * i + 1] = target[i][1] / l;
                nrm.As<double>()[3 * i + 2] = target[i][2] / 0.49 / l;
            }
            tmesh.AddPointData("normals", std::move(nrm));
        }
        const RegistrationResult r = register_rigid(rg_mesh(source), tmesh, o);
        EXPECT_TRUE(r.mConverged) << static_cast<int>(method);
        EXPECT_LT(r.mRms, 1e-9);
        EXPECT_GT(r.mInitialRms, 0.05);
        // r.mMatrix is the inverse motion.
        const auto back = rg_moved(source, r.mMatrix);
        double worst = 0.0;
        for (std::size_t i = 0; i < target.size(); ++i)
            for (std::size_t a = 0; a < 3; ++a)
                worst = std::max(worst, std::abs(back[i][a] - target[i][a]));
        EXPECT_LT(worst, 1e-9) << static_cast<int>(method);
        EXPECT_EQ(r.mNumPairs, 3000);
    }
}

TEST(Registration, PointToPlaneOntoASurface) {
    // The source: the tet fixture's points moved slightly; the target: the
    // cube surface they sample.
    Mesh cube = mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "triangle",
        {{0, 2, 1},
         {0, 3, 2},
         {4, 5, 6},
         {4, 6, 7},
         {0, 1, 5},
         {0, 5, 4},
         {1, 2, 6},
         {1, 6, 5},
         {2, 3, 7},
         {2, 7, 6},
         {3, 0, 4},
         {3, 4, 7}});
    std::vector<RgV> samples;
    std::mt19937_64 rng(2);
    std::uniform_real_distribution<double> u(0.05, 0.95);
    for (int f = 0; f < 6; ++f)
        for (int i = 0; i < 80; ++i) {
            const double a = u(rng), b = u(rng);
            const double side = f % 2;
            RgV p = f / 2 == 0 ? RgV{side, a, b} : f / 2 == 1 ? RgV{a, side, b} : RgV{a, b, side};
            samples.push_back(p);
        }
    const auto source = rg_moved(samples, rg_motion(0.05, {0.03, 0.02, -0.04}));
    RegistrationOptions o;
    o.mMethod = RegistrationMethod::PointToPlane;
    const RegistrationResult r = register_rigid(rg_mesh(source), cube, o);
    EXPECT_TRUE(r.mConverged);
    EXPECT_LT(r.mRms, 1e-9);
    const auto back = rg_moved(source, r.mMatrix);
    double worst = 0.0;
    for (std::size_t i = 0; i < samples.size(); ++i)
        for (std::size_t a = 0; a < 3; ++a)
            worst = std::max(worst, std::abs(back[i][a] - samples[i][a]));
    EXPECT_LT(worst, 1e-9);
}

TEST(Registration, Refusals) {
    const auto pts = rg_cloud(50);
    RegistrationOptions o;
    o.mMethod = RegistrationMethod::PointToPlane;
    EXPECT_THROW(register_rigid(rg_mesh(pts), rg_mesh(pts), o), std::invalid_argument);
    o = RegistrationOptions{};
    o.mMaxIterations = 0;
    EXPECT_THROW(register_rigid(rg_mesh(pts), rg_mesh(pts), o), std::invalid_argument);
    o = RegistrationOptions{};
    o.mMaxDistance = 1e-6;
    auto far = pts;
    for (auto& q : far)
        q[0] += 100.0;
    EXPECT_THROW(register_rigid(rg_mesh(far), rg_mesh(pts), o), std::invalid_argument);
    EXPECT_THROW(register_rigid(rg_mesh(pts), mt::tet_mesh()), std::invalid_argument);
    EXPECT_THROW(registration_method_from_name("affine"), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Alignment beyond local ICP (roadmap 7.3.1, 7.3.2, 7.3.3, 7.3.5, 7.3.6).
// ---------------------------------------------------------------------------

#include "meshioplusplus/operations/auto_align.hpp"
#include "meshioplusplus/operations/deviation.hpp"
#include "meshioplusplus/operations/transform.hpp"

namespace {

// Append an axis-aligned box's surface: K x K quads per face, two outward
// triangles each.
void rg_box(const RgV& rLo, const RgV& rHi, int K, std::vector<RgV>& rPts,
            std::vector<std::array<std::int64_t, 3>>& rTris, bool SkipZHi = false) {
    for (int a = 0; a < 3; ++a)
        for (int side = 0; side < 2; ++side) {
            if (SkipZHi && a == 2 && side == 1)
                continue;
            const int b = (a + 1) % 3, c = (a + 2) % 3;
            const std::int64_t base = static_cast<std::int64_t>(rPts.size());
            for (int i = 0; i <= K; ++i)
                for (int j = 0; j <= K; ++j) {
                    RgV p{};
                    p[static_cast<std::size_t>(a)] = side == 0 ? rLo[static_cast<std::size_t>(a)]
                                                               : rHi[static_cast<std::size_t>(a)];
                    p[static_cast<std::size_t>(b)] =
                        rLo[static_cast<std::size_t>(b)] +
                        (rHi[static_cast<std::size_t>(b)] - rLo[static_cast<std::size_t>(b)]) * i /
                            K;
                    p[static_cast<std::size_t>(c)] =
                        rLo[static_cast<std::size_t>(c)] +
                        (rHi[static_cast<std::size_t>(c)] - rLo[static_cast<std::size_t>(c)]) * j /
                            K;
                    rPts.push_back(p);
                }
            const auto id = [&](int i, int j) { return base + i * (K + 1) + j; };
            for (int i = 0; i < K; ++i)
                for (int j = 0; j < K; ++j) {
                    std::array<std::int64_t, 3> t1{id(i, j), id(i + 1, j), id(i + 1, j + 1)};
                    std::array<std::int64_t, 3> t2{id(i, j), id(i + 1, j + 1), id(i, j + 1)};
                    if (side == 0) {
                        std::swap(t1[1], t1[2]);
                        std::swap(t2[1], t2[2]);
                    }
                    rTris.push_back(t1);
                    rTris.push_back(t2);
                }
        }
}

Mesh rg_tri_mesh(const std::vector<RgV>& rPts,
                 const std::vector<std::array<std::int64_t, 3>>& rTris) {
    Mesh m = rg_mesh(rPts);
    NDArray c(DType::Int64, {rTris.size(), std::size_t{3}});
    for (std::size_t t = 0; t < rTris.size(); ++t)
        for (std::size_t k = 0; k < 3; ++k)
            c.As<std::int64_t>()[3 * t + k] = rTris[t][k];
    m.AddCellBlock("triangle", std::move(c));
    return m;
}

// An L bracket with a rib in one corner: no proper symmetry.
Mesh rg_bracket(std::vector<RgV>* pPts = nullptr) {
    std::vector<RgV> pts;
    std::vector<std::array<std::int64_t, 3>> tris;
    rg_box({0, 0, 0}, {4, 1, 0.5}, 8, pts, tris);
    rg_box({0, 0, 0.5}, {0.5, 1, 3}, 8, pts, tris);
    rg_box({3, 0.6, 0.5}, {4, 1, 1.2}, 4, pts, tris);
    if (pPts != nullptr)
        *pPts = pts;
    return rg_tri_mesh(pts, tris);
}

double rg_max_error(const std::vector<RgV>& rA, const std::vector<RgV>& rB) {
    double e = 0.0;
    for (std::size_t i = 0; i < rA.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            e = std::max(e, std::abs(rA[i][k] - rB[i][k]));
    return e;
}

std::vector<RgV> rg_points_of(const Mesh& rMesh) {
    std::vector<RgV> out(rMesh.NumPoints());
    for (std::size_t i = 0; i < out.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            out[i][k] = rMesh.Points().As<double>()[3 * i + k];
    return out;
}

}  // namespace

TEST(RegistrationAuto, LargeRotationRegistersBackFromThePrincipalFrames) {
    std::vector<RgV> pts;
    const Mesh target = rg_bracket(&pts);
    const auto motion = rg_motion(150.0 * 3.14159265358979323846 / 180.0, {2.0, -1.0, 0.5});
    Mesh source = detail::clone_mesh(target);
    source.AssignPoints(rg_mesh(rg_moved(pts, motion)).Points());
    RegistrationOptions o;
    o.mMethod = RegistrationMethod::PointToPlane;
    o.mMaxIterations = 200;
    o.mTolerance = 1e-14;
    // The identity start settles in a wrong local optimum.
    const RegistrationResult local = register_rigid(source, target, o);
    EXPECT_GT(rg_max_error(rg_moved(rg_moved(pts, motion), local.mMatrix), pts), 1e-3);
    o.mInitialMode = RegistrationInitial::Auto;
    const RegistrationResult r = register_rigid(source, target, o);
    ASSERT_EQ(r.mCandidateRms.size(), 24u);
    EXPECT_GE(r.mChosenCandidate, 0);
    EXPECT_LT(rg_max_error(rg_moved(rg_moved(pts, motion), r.mMatrix), pts), 1e-9);
    EXPECT_TRUE(r.mWarnings.empty());
    EXPECT_EQ(r.mNumRejected, 0);
    // Deterministic: a second run is bit-identical.
    const RegistrationResult again = register_rigid(source, target, o);
    EXPECT_EQ(again.mMatrix, r.mMatrix);
    EXPECT_EQ(again.mCandidateRms, r.mCandidateRms);
}

TEST(RegistrationAuto, SymmetricCheckWarnsOnABox) {
    std::vector<RgV> pts;
    std::vector<std::array<std::int64_t, 3>> tris;
    rg_box({0, 0, 0}, {3, 2, 1}, 6, pts, tris);
    const Mesh target = rg_tri_mesh(pts, tris);
    Mesh source = detail::clone_mesh(target);
    source.AssignPoints(rg_mesh(rg_moved(pts, rg_motion(1.0, {0.5, 0, 0}))).Points());
    RegistrationOptions o;
    o.mInitialMode = RegistrationInitial::Auto;
    o.mSymmetricCheck = true;
    const RegistrationResult r = register_rigid(source, target, o);
    ASSERT_FALSE(r.mWarnings.empty());
    EXPECT_NE(r.mWarnings[0].find("symmetric"), std::string::npos);
    EXPECT_LT(r.mRms, 1e-6);
}

TEST(RegistrationRegions, FixturePlateIsExcludedByTheSourceRegion) {
    std::vector<RgV> pts;
    std::vector<std::array<std::int64_t, 3>> tris;
    rg_box({0, 0, 0}, {4, 1, 0.5}, 8, pts, tris);
    rg_box({0, 0, 0.5}, {0.5, 1, 3}, 8, pts, tris);
    rg_box({3, 0.6, 0.5}, {4, 1, 1.2}, 4, pts, tris);
    const std::size_t part_cells = tris.size();
    const Mesh target = rg_tri_mesh(pts, tris);
    // The scan: the part plus a fixture plate the nominal lacks, slightly off.
    rg_box({-1, -1, -0.6}, {5, 2, -0.1}, 10, pts, tris);
    const auto motion = rg_motion(0.05, {0.05, 0.02, -0.03});
    Mesh scan = rg_tri_mesh(rg_moved(pts, motion), tris);
    NDArray entries(DType::Int64, {part_cells});
    for (std::size_t c = 0; c < part_cells; ++c)
        entries.As<std::int64_t>()[c] = static_cast<std::int64_t>(c);
    scan.AddRegion(meshioplusplus::Region("part", RegionKind::Cell, std::move(entries)));
    RegistrationOptions o;
    o.mMethod = RegistrationMethod::PointToPlane;
    o.mMaxIterations = 200;
    o.mTolerance = 1e-14;
    const auto moved_back = [&](const RegistrationResult& rR) {
        return rg_max_error(rg_moved(rg_moved(pts, motion), rR.mMatrix), pts);
    };
    const RegistrationResult pulled = register_rigid(scan, target, o);
    EXPECT_GT(moved_back(pulled), 1e-3);
    o.mSourceRegion = "part";
    const RegistrationResult r = register_rigid(scan, target, o);
    EXPECT_LT(moved_back(r), 1e-6);
    // Adaptive trimming rejects the plate's far pairs without a region.
    o.mSourceRegion.clear();
    o.mAdaptiveTrim = 3.0;
    const RegistrationResult t = register_rigid(scan, target, o);
    EXPECT_GT(t.mNumRejected, 0);
    EXPECT_LT(moved_back(t), moved_back(pulled));
    o.mSourceRegion = "nope";
    EXPECT_THROW(register_rigid(scan, target, o), std::invalid_argument);
    o.mSourceRegion.clear();
    o.mTargetRegion = "nope";
    EXPECT_THROW(register_rigid(scan, target, o), std::invalid_argument);
    o.mTargetRegion.clear();
    o.mAdaptiveTrim = -1.0;
    EXPECT_THROW(register_rigid(scan, target, o), std::invalid_argument);
}

TEST(RegisterLandmarks, ThreePairsRecoverATransformAndTwoAreRefused) {
    const auto motion = rg_motion(2.0, {1.0, -2.0, 3.0});
    const std::vector<RgV> src{{0, 0, 0}, {1, 0, 0}, {0, 2, 0}};
    const std::vector<RgV> dst = rg_moved(src, motion);
    const LandmarkResult r = register_landmarks(src, dst);
    for (std::size_t k = 0; k < 16; ++k)
        EXPECT_NEAR(r.mMatrix[k], motion[k], 1e-12);
    ASSERT_EQ(r.mResiduals.size(), 3u);
    EXPECT_LT(r.mRms, 1e-12);
    EXPECT_THROW(register_landmarks({src[0], src[1]}, {dst[0], dst[1]}), std::invalid_argument);
    EXPECT_THROW(
        register_landmarks({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}, {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}),
        std::invalid_argument);
    EXPECT_THROW(register_landmarks(src, dst, {1.0, -1.0, 1.0}), std::invalid_argument);
    // A zero weight leaves an outlier out.
    std::vector<RgV> src4 = src, dst4 = dst;
    src4.push_back({5, 5, 5});
    dst4.push_back({100, 0, 0});
    const LandmarkResult w = register_landmarks(src4, dst4, {1, 1, 1, 0});
    for (std::size_t k = 0; k < 16; ++k)
        EXPECT_NEAR(w.mMatrix[k], motion[k], 1e-12);
    EXPECT_GT(w.mResiduals[3], 10.0);
}

TEST(SymmetryPlane, TiltedPlaneOfANoisyBracketIsFound) {
    std::vector<RgV> pts;
    std::vector<std::array<std::int64_t, 3>> tris;
    // Symmetric about x = 0.
    rg_box({-2, 0, 0}, {2, 1, 0.5}, 8, pts, tris);
    rg_box({-0.25, 0, 0.5}, {0.25, 1, 3}, 8, pts, tris);
    rg_box({1.5, 0.6, 0.5}, {2, 1, 1.2}, 4, pts, tris);
    rg_box({-2, 0.6, 0.5}, {-1.5, 1, 1.2}, 4, pts, tris);
    const auto motion = rg_motion(0.6, {0.3, -0.2, 0.1});
    std::vector<RgV> moved = rg_moved(pts, motion);
    std::mt19937_64 rng(11);
    std::normal_distribution<double> g(0.0, 2e-4);
    for (RgV& p : moved)
        for (double& v : p)
            v += g(rng);
    const Mesh mesh = rg_tri_mesh(moved, tris);
    const SymmetryPlaneResult r = find_symmetry_plane(mesh);
    const RgV n{motion[0], motion[4], motion[8]};  // R e_x
    const double c = std::abs(r.mNormal[0] * n[0] + r.mNormal[1] * n[1] + r.mNormal[2] * n[2]);
    EXPECT_LT(std::acos(std::min(1.0, c)), 1e-4);
    // The plane passes through the moved origin's image of x = 0: t.
    const double off = (r.mOrigin[0] - motion[3]) * n[0] + (r.mOrigin[1] - motion[7]) * n[1] +
                       (r.mOrigin[2] - motion[11]) * n[2];
    EXPECT_LT(std::abs(off), 1e-4 * std::sqrt(16.0 + 1.0 + 9.0));
    EXPECT_GT(r.mShare, 0.9);
    EXPECT_EQ(r.mCandidateRms.size(), 3u);
    // The bracket has no such plane.
    SymmetryPlaneOptions strict;
    strict.mMinShare = 0.999;
    EXPECT_THROW(find_symmetry_plane(rg_bracket(), strict), std::invalid_argument);
}

TEST(AutoAlign, RotatedBlockComesBackAxisAligned) {
    std::vector<RgV> pts;
    std::vector<std::array<std::int64_t, 3>> tris;
    rg_box({0, 0, 0}, {3, 2, 1}, 6, pts, tris);
    const Mesh block = rg_tri_mesh(rg_moved(pts, rg_motion(0.9, {4, 5, 6})), tris);
    const AutoAlignResult r = auto_align(block);
    EXPECT_EQ(r.mMode, "face");
    EXPECT_GT(r.mShare, 0.9);
    const std::vector<RgV> q = rg_points_of(r.mMesh);
    RgV lo = q[0], hi = q[0];
    for (const RgV& p : q)
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    EXPECT_NEAR(lo[2], 0.0, 1e-12);
    EXPECT_NEAR(hi[2] - lo[2], 1.0, 1e-9);  // lies on a 3 x 2 face
    EXPECT_NEAR(hi[0] - lo[0], 2.0, 1e-9);  // the larger side faces +-X
    EXPECT_NEAR(hi[1] - lo[1], 3.0, 1e-9);
    EXPECT_NEAR(lo[0] + hi[0], 0.0, 1e-9);
    EXPECT_DOUBLE_EQ(r.mOpeningArea, 0.0);
}

TEST(AutoAlign, TubeStandsOnItsAxisAndAnOpenTrayOnItsOpening) {
    // An open tube along (1, 1, 0) / sqrt 2.
    std::vector<RgV> pts;
    std::vector<std::array<std::int64_t, 3>> tris;
    const int ns = 64, nz = 10;
    for (int k = 0; k <= nz; ++k)
        for (int s = 0; s < ns; ++s) {
            const double a = 2 * 3.14159265358979323846 * s / ns;
            pts.push_back({std::cos(a), std::sin(a), 4.0 * k / nz});
        }
    for (int k = 0; k < nz; ++k)
        for (int s = 0; s < ns; ++s) {
            const std::int64_t a = k * ns + s, b = k * ns + (s + 1) % ns;
            tris.push_back({a, b, b + ns});
            tris.push_back({a, b + ns, a + ns});
        }
    const auto motion = rg_motion(1.1, {0, 0, 0});
    const AutoAlignResult r = auto_align(rg_tri_mesh(rg_moved(pts, motion), tris));
    EXPECT_EQ(r.mMode, "axis");
    const RgV axis{motion[2], motion[6], motion[10]};
    EXPECT_NEAR(
        std::abs(r.mPrimary[0] * axis[0] + r.mPrimary[1] * axis[1] + r.mPrimary[2] * axis[2]), 1.0,
        1e-9);
    EXPECT_GT(r.mShare, 0.9);
    // A tray open at its top: it lies on the opening.
    pts.clear();
    tris.clear();
    rg_box({0, 0, 0}, {3, 2, 1}, 4, pts, tris, true);
    const auto m2 = rg_motion(0.4, {1, 2, 3});
    const AutoAlignResult t = auto_align(rg_tri_mesh(rg_moved(pts, m2), tris));
    EXPECT_NEAR(t.mOpeningArea, 6.0, 1e-9);
    const RgV up{m2[2], m2[6], m2[10]};  // the opening's outward normal
    // The bottom and the opening have the same area: either way the tray
    // lies on a 3 x 2 side.
    EXPECT_NEAR(std::abs(t.mPrimary[0] * up[0] + t.mPrimary[1] * up[1] + t.mPrimary[2] * up[2]),
                1.0, 1e-9);
    EXPECT_THROW(auto_align(rg_mesh(pts)), std::invalid_argument);
    AutoAlignOptions bad;
    bad.mAngleToleranceDeg = 50.0;
    EXPECT_THROW(auto_align(rg_tri_mesh(pts, tris), bad), std::invalid_argument);
}

TEST(DeviationRegister, ScanIsRegisteredBeforeItIsMeasured) {
    std::vector<RgV> pts;
    const Mesh reference = rg_bracket(&pts);
    Mesh scan = detail::clone_mesh(reference);
    scan.AssignPoints(rg_mesh(rg_moved(pts, rg_motion(2.5, {1, 2, 3}))).Points());
    DeviationOptions o;
    o.mTolerance = 1e-6;
    const DeviationResult before = deviation(scan, reference, o);
    EXPECT_LT(before.mShareWithin, 0.5);
    o.mRegister = true;
    const DeviationResult r = deviation(scan, reference, o);
    EXPECT_TRUE(r.mRegistered);
    EXPECT_GT(r.mShareWithin, 0.999);
    EXPECT_LT(r.mRegistrationRms, 1e-9);
    EXPECT_THROW(deviation_to_primitive(scan, 1, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, o),
                 std::invalid_argument);
}
