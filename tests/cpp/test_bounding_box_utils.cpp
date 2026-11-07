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
// Unit tests for detail/bounding_box_utils.hpp -- BoundingBoxUtils, the
// first-point-seeded axis-aligned boxes that replaced hand-written loops.

// System includes
#include <array>
#include <cmath>
#include <limits>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/bounding_box_utils.hpp"

using meshioplusplus::detail::BoundingBox3;
using meshioplusplus::detail::BoundingBoxUtils;
using meshioplusplus::detail::Vec3;

TEST(BoundingBoxUtils, FromPointsAndMeasures) {
    const std::vector<Vec3> p{{1, 2, 3}, {-1, 5, 0}, {4, -2, 1}};
    const BoundingBox3 box = BoundingBoxUtils::FromPoints(p);
    EXPECT_EQ(box.mMin, (Vec3{-1, -2, 0}));
    EXPECT_EQ(box.mMax, (Vec3{4, 5, 3}));
    EXPECT_EQ(BoundingBoxUtils::Extent(box), (Vec3{5, 7, 3}));
    EXPECT_EQ(BoundingBoxUtils::MaxExtent(box), 7.0);
    EXPECT_EQ(BoundingBoxUtils::HalfSurfaceArea(box), 5.0 * 7.0 + 7.0 * 3.0 + 3.0 * 5.0);
    EXPECT_DOUBLE_EQ(BoundingBoxUtils::Diagonal(box), std::sqrt(25.0 + 49.0 + 9.0));
}

TEST(BoundingBoxUtils, SinglePointAndEmptySet) {
    const BoundingBox3 one = BoundingBoxUtils::FromPoints({Vec3{1, 2, 3}});
    EXPECT_EQ(one.mMin, (Vec3{1, 2, 3}));
    EXPECT_EQ(one.mMax, (Vec3{1, 2, 3}));
    EXPECT_EQ(BoundingBoxUtils::MaxExtent(one), 0.0);
    const BoundingBox3 none = BoundingBoxUtils::FromPoints({});
    EXPECT_EQ(none.mMin, (Vec3{0, 0, 0}));
    EXPECT_EQ(none.mMax, (Vec3{0, 0, 0}));
}

TEST(BoundingBoxUtils, NaNCoordinatesKeepTheSeed) {
    // std::min(lo, NaN) keeps lo: the copies' behaviour, which a loop seeded
    // with +-infinity does not share.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const BoundingBox3 box = BoundingBoxUtils::FromPoints({Vec3{0, 0, 0}, Vec3{nan, 1, -1}});
    EXPECT_EQ(box.mMin, (Vec3{0, 0, -1}));
    EXPECT_EQ(box.mMax, (Vec3{0, 1, 0}));
}

TEST(BoundingBoxUtils, FromCoordsMatchesFromPoints) {
    const std::vector<Vec3> p{{3, -1, 2}, {0, 0, 0}, {-5, 9, 1}, {2, 2, -7}};
    double lo[3], hi[3];
    BoundingBoxUtils::FromCoords(p.size(), [&](std::size_t i) { return p[i]; }, lo, hi);
    const BoundingBox3 box = BoundingBoxUtils::FromPoints(p);
    for (int d = 0; d < 3; ++d) {
        EXPECT_EQ(lo[d], box.mMin[static_cast<std::size_t>(d)]);
        EXPECT_EQ(hi[d], box.mMax[static_cast<std::size_t>(d)]);
    }
    // One coordinate: the box is the point.
    BoundingBoxUtils::FromCoords(1, [&](std::size_t) { return Vec3{4, 5, 6}; }, lo, hi);
    EXPECT_EQ(lo[0], 4.0);
    EXPECT_EQ(hi[2], 6.0);
}

TEST(BoundingBoxUtils, ExpandWidensOnlyWhereNeeded) {
    BoundingBox3 box{{0, 0, 0}, {1, 1, 1}};
    BoundingBoxUtils::Expand(box, Vec3{0.5, -2, 3});
    EXPECT_EQ(box.mMin, (Vec3{0, -2, 0}));
    EXPECT_EQ(box.mMax, (Vec3{1, 1, 3}));
}
