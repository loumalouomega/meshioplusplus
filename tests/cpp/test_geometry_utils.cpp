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
// Unit tests for detail/geometry_utils.hpp -- GeometryUtils, the triangle
// measures that replaced hand-written copies across the operations.

// System includes
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/geometry_utils.hpp"

using meshioplusplus::detail::GeometryUtils;
using meshioplusplus::detail::Vec3;

namespace {

std::uint64_t bits(double X) {
    return std::bit_cast<std::uint64_t>(X);
}

std::vector<std::array<Vec3, 3>> sample_triangles() {
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> u(-10.0, 10.0);
    std::vector<std::array<Vec3, 3>> out;
    for (int i = 0; i < 100; ++i)
        out.push_back({Vec3{u(rng), u(rng), u(rng)}, Vec3{u(rng), u(rng), u(rng)},
                       Vec3{u(rng), u(rng), u(rng)}});
    // Degenerate: collinear and coincident corners.
    out.push_back({Vec3{0, 0, 0}, Vec3{1, 1, 1}, Vec3{2, 2, 2}});
    out.push_back({Vec3{1, 2, 3}, Vec3{1, 2, 3}, Vec3{1, 2, 3}});
    return out;
}

}  // namespace

TEST(GeometryUtils, TriangleKnownValues) {
    const Vec3 a{0, 0, 0}, b{2, 0, 0}, c{0, 3, 0};
    EXPECT_EQ(GeometryUtils::TriangleNormal(a, b, c), (Vec3{0.0, 0.0, 6.0}));
    EXPECT_EQ(GeometryUtils::TriangleNormal(a, c, b), (Vec3{0.0, 0.0, -6.0}));
    EXPECT_DOUBLE_EQ(GeometryUtils::TriangleTwiceArea(a, b, c), 6.0);
    EXPECT_DOUBLE_EQ(GeometryUtils::TriangleArea(a, b, c), 3.0);
    EXPECT_EQ(GeometryUtils::TriangleArea(a, a, c), 0.0);
}

TEST(GeometryUtils, TriangleAnglesDegSumTo180) {
    const auto ang = GeometryUtils::TriangleAnglesDeg(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0});
    EXPECT_NEAR(ang[0], 90.0, 1e-12);
    EXPECT_NEAR(ang[1], 45.0, 1e-12);
    EXPECT_NEAR(ang[2], 45.0, 1e-12);
    for (const auto& t : sample_triangles()) {
        const auto a = GeometryUtils::TriangleAnglesDeg(t[0], t[1], t[2]);
        if (GeometryUtils::TriangleTwiceArea(t[0], t[1], t[2]) > 1e-6) {
            EXPECT_NEAR(a[0] + a[1] + a[2], 180.0, 1e-9);
        }
    }
}

TEST(GeometryUtils, MatchesTheHandWrittenCopiesBitForBit) {
    for (const auto& t : sample_triangles()) {
        const Vec3 &a = t[0], &b = t[1], &c = t[2];
        // The copy every replaced site wrote out.
        const Vec3 u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const Vec3 v{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const Vec3 n{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                     u[0] * v[1] - u[1] * v[0]};
        const double twice = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        const Vec3 got = GeometryUtils::TriangleNormal(a, b, c);
        for (std::size_t k = 0; k < 3; ++k)
            EXPECT_EQ(bits(got[k]), bits(n[k]));
        EXPECT_EQ(bits(GeometryUtils::TriangleTwiceArea(a, b, c)), bits(twice));
        EXPECT_EQ(bits(GeometryUtils::TriangleArea(a, b, c)), bits(0.5 * twice));
        // The flat-triple overloads agree with the Vec3 ones.
        const Vec3 gp = GeometryUtils::TriangleNormal(a.data(), b.data(), c.data());
        for (std::size_t k = 0; k < 3; ++k)
            EXPECT_EQ(bits(gp[k]), bits(n[k]));
        EXPECT_EQ(bits(GeometryUtils::TriangleArea(a.data(), b.data(), c.data())),
                  bits(0.5 * twice));
        EXPECT_EQ(bits(GeometryUtils::TriangleTwiceArea(a.data(), b.data(), c.data())),
                  bits(twice));
    }
}

TEST(GeometryUtils, SignedTwiceArea2D) {
    EXPECT_DOUBLE_EQ(GeometryUtils::SignedTwiceArea2D(0, 0, 1, 0, 0, 1), 1.0);
    EXPECT_DOUBLE_EQ(GeometryUtils::SignedTwiceArea2D(0, 0, 0, 1, 1, 0), -1.0);
    EXPECT_EQ(GeometryUtils::SignedTwiceArea2D(0, 0, 1, 1, 2, 2), 0.0);
    std::mt19937_64 rng(3);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (int i = 0; i < 200; ++i) {
        const double u0 = u(rng), v0 = u(rng), u1 = u(rng), v1 = u(rng), u2 = u(rng), v2 = u(rng);
        const double copy = (u1 - u0) * (v2 - v0) - (u2 - u0) * (v1 - v0);
        EXPECT_EQ(bits(GeometryUtils::SignedTwiceArea2D(u0, v0, u1, v1, u2, v2)), bits(copy));
    }
}
