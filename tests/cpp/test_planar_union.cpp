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
// DXF R2000 drawings (roadmap 7.6.4).
// detail/planar_union: the region a soup of planar triangles covers.

#include <array>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/planar_union.hpp"

using meshioplusplus::detail::planar_loop_twice_area;
using meshioplusplus::detail::planar_region;
using meshioplusplus::detail::PlanarRegion;
using meshioplusplus::detail::PlanarSoup;

namespace {

PlanarSoup pu_square(double X, double Y, double S) {
    return PlanarSoup{{{X, Y}, {X + S, Y}, {X + S, Y + S}, {X, Y + S}}, {{0, 1, 2}, {0, 2, 3}}};
}

void pu_append(PlanarSoup& rInto, const PlanarSoup& rFrom) {
    const auto base = static_cast<std::int64_t>(rInto.mPoints.size());
    rInto.mPoints.insert(rInto.mPoints.end(), rFrom.mPoints.begin(), rFrom.mPoints.end());
    for (const auto& r_t : rFrom.mTriangles)
        rInto.mTriangles.push_back({r_t[0] + base, r_t[1] + base, r_t[2] + base});
}

double pu_area(const PlanarRegion& rR) {
    double a = 0.0;
    for (const auto& r_l : rR.mLoops)
        a += 0.5 * planar_loop_twice_area(rR.mPoints, r_l);
    return a;
}

}  // namespace

TEST(PlanarUnion, OverlappingSquaresAreOneLoop) {
    PlanarSoup s = pu_square(0, 0, 2);
    pu_append(s, pu_square(1, 1, 2));
    const PlanarRegion r = planar_region(s, {}, 0.0);
    ASSERT_TRUE(r.mOk);
    ASSERT_EQ(r.mLoops.size(), 1u);
    EXPECT_EQ(r.mLoops[0].size(), 8u);  // two crossing points added
    EXPECT_DOUBLE_EQ(pu_area(r), 7.0);
}

TEST(PlanarUnion, SquaresTouchingAtACornerStayTwoLoops) {
    PlanarSoup s = pu_square(0, 0, 1);
    pu_append(s, pu_square(1, 1, 1));
    const PlanarRegion r = planar_region(s, {}, 0.0);
    ASSERT_TRUE(r.mOk);
    EXPECT_EQ(r.mLoops.size(), 2u);
    EXPECT_DOUBLE_EQ(pu_area(r), 2.0);
}

TEST(PlanarUnion, SharedAndCollinearEdgesMerge) {
    PlanarSoup s = pu_square(0, 0, 1);
    pu_append(s, pu_square(1, 0, 1));       // shares an edge, through distinct points
    pu_append(s, pu_square(0.5, 0, 0.25));  // inside, on the bottom edge
    pu_append(s, pu_square(0, 0, 1));       // a duplicate, reversed below
    std::swap(s.mTriangles.back()[1], s.mTriangles.back()[2]);
    const PlanarRegion r = planar_region(s, {}, 0.0);
    ASSERT_TRUE(r.mOk);
    ASSERT_EQ(r.mLoops.size(), 1u);
    EXPECT_DOUBLE_EQ(pu_area(r), 2.0);
}

TEST(PlanarUnion, SubtractingMakesAHole) {
    const PlanarRegion r = planar_region(pu_square(0, 0, 3), pu_square(1, 1, 1), 0.0);
    ASSERT_TRUE(r.mOk);
    ASSERT_EQ(r.mLoops.size(), 2u);
    int holes = 0;
    for (const auto& r_l : r.mLoops)
        holes += planar_loop_twice_area(r.mPoints, r_l) < 0.0 ? 1 : 0;
    EXPECT_EQ(holes, 1);
    EXPECT_DOUBLE_EQ(pu_area(r), 8.0);
    EXPECT_TRUE(planar_region(pu_square(0, 0, 1), pu_square(-1, -1, 3), 0.0).mLoops.empty());
}

TEST(PlanarUnion, MergesWithinTheTolerance) {
    PlanarSoup s = pu_square(0, 0, 1);
    pu_append(s, pu_square(1.0 + 1e-13, 0, 1));
    EXPECT_EQ(planar_region(s, {}, 0.0).mLoops.size(), 2u);
    EXPECT_EQ(planar_region(s, {}, 1e-9).mLoops.size(), 1u);
    EXPECT_FALSE(planar_region(PlanarSoup{{{0, 0}, {1e300, 0}, {0, 1}}, {{0, 1, 2}}}, {}, 0.0).mOk);
}
