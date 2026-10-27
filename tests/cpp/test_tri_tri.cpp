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
// Tests for the exact triangle-triangle kernel (roadmap 6.9.1).

#include <array>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/tri_tri.hpp"

namespace {

namespace mdt = meshioplusplus::detail;
using V = std::array<double, 3>;

bool meet(const V& rA, const V& rB, const V& rC, const V& rD, const V& rE, const V& rF) {
    const bool x =
        mdt::tri_tri_intersect(rA.data(), rB.data(), rC.data(), rD.data(), rE.data(), rF.data());
    // Symmetric in the two triangles and in their corners' order.
    EXPECT_EQ(x, mdt::tri_tri_intersect(rD.data(), rE.data(), rF.data(), rA.data(), rB.data(),
                                        rC.data()));
    EXPECT_EQ(x, mdt::tri_tri_intersect(rB.data(), rA.data(), rC.data(), rF.data(), rD.data(),
                                        rE.data()));
    return x;
}

}  // namespace

TEST(TriTri, CrossingApartAndTouching) {
    const V a{0, 0, 0}, b{1, 0, 0}, c{0, 1, 0};
    // A triangle through the first's inside, standing on it.
    EXPECT_TRUE(meet(a, b, c, V{0.2, 0.2, -1}, V{0.2, 0.2, 1}, V{0.3, 0.25, 1}));
    // Above it.
    EXPECT_FALSE(meet(a, b, c, V{0, 0, 1}, V{1, 0, 1}, V{0, 1, 1}));
    // Through its plane but beside it.
    EXPECT_FALSE(meet(a, b, c, V{2, 2, -1}, V{2, 2, 1}, V{3, 2, 1}));
    // Touching at a corner, and along an edge.
    EXPECT_TRUE(meet(a, b, c, V{1, 0, 0}, V{2, 0, 1}, V{2, 1, 1}));
    EXPECT_TRUE(meet(a, b, c, V{0, 0, 0}, V{1, 0, 0}, V{0.5, 0, 1}));
    // An edge passing exactly through the first's edge.
    EXPECT_TRUE(meet(a, b, c, V{0.5, 0, -1}, V{0.5, 0, 1}, V{0.5, -1, 0}));
}

TEST(TriTri, InOnePlane) {
    const V a{0, 0, 0}, b{1, 0, 0}, c{0, 1, 0};
    EXPECT_TRUE(meet(a, b, c, V{0.1, 0.1, 0}, V{0.2, 0.1, 0}, V{0.1, 0.2, 0}));  // inside
    EXPECT_TRUE(meet(a, b, c, V{0.5, -0.5, 0}, V{0.5, 0.5, 0}, V{2, 0, 0}));     // crossing
    EXPECT_FALSE(meet(a, b, c, V{1, 1, 0}, V{2, 1, 0}, V{1, 2, 0}));             // apart
    EXPECT_TRUE(meet(a, b, c, V{1, 0, 0}, V{2, 0, 0}, V{1, 1, 0}));              // a corner
    // Tilted out of the axes' planes: still exact.
    const V p{0, 0, 0}, q{1, 1, 1}, r{1, 0, 1};
    EXPECT_TRUE(meet(p, q, r, V{0.8, 0.4, 0.8}, V{2, 1, 2}, V{2, 0, 2}));
    EXPECT_FALSE(meet(p, q, r, V{2, 2, 2}, V{3, 3, 3}, V{3, 2, 3}));
}

TEST(TriTri, NearlyTouchingIsDecidedExactly) {
    // The second triangle lies 1e-300 under the first's plane: apart.
    const V a{0, 0, 0}, b{1, 0, 0}, c{0, 1, 0};
    EXPECT_FALSE(meet(a, b, c, V{0.2, 0.2, -1e-300}, V{0.3, 0.2, -1}, V{0.2, 0.3, -1}));
    // A segment through the triangle's corner exactly.
    const V s0{0, 0, -1}, s1{0, 0, 1};
    EXPECT_TRUE(mdt::segment_tri_intersect(s0.data(), s1.data(), a.data(), b.data(), c.data()));
    const V t0{-1e-17, 0, -1}, t1{-1e-17, 0, 1};
    EXPECT_FALSE(mdt::segment_tri_intersect(t0.data(), t1.data(), a.data(), b.data(), c.data()));
}

// A regression: two triangles side by side in one plane, separated by a
// third, nearly coplanar to rounding (the cube's offset caps): apart. A
// segment of one lies in the other's plane, collinear with a corner of it
// but beyond its extent.
TEST(TriTri, SideBySideInOnePlaneAreApart) {
    const V a{0.75000000000000011, 0.13397459621556138, -0.061870000000000008};
    const V b{0.75, -0.061870000000000008, -0.061870000000000008};
    const V c{0.5, 0.13397459621556138, -0.061870000000000008};
    const V d{0.25000000000000006, 0.13397459621556138, -0.061870000000000008};
    const V e{0.49999999999999994, -0.061870000000000008, -0.061870000000000015};
    const V f{0.25, -0.061870000000000008, -0.061870000000000015};
    EXPECT_FALSE(meet(a, b, c, d, e, f));
}
