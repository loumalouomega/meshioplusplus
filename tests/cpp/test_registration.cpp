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
