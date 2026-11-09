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
#include <cmath>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/colormap.hpp"
#include "meshioplusplus/operations/colorize.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// Ten points on a line carrying `v` = 0, 1, ..., 9 (and one NaN at the end).
Mesh col_line() {
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> lines;
    for (int i = 0; i < 11; ++i)
        pts.push_back({double(i), 0.0, 0.0});
    for (int i = 0; i < 10; ++i)
        lines.push_back({i, i + 1});
    Mesh m = mt::make_mesh(pts, "line", lines);
    NDArray v = NDArray::Uninit(DType::Float64, {11});
    for (int i = 0; i < 10; ++i)
        v.As<double>()[i] = i;
    v.As<double>()[10] = std::nan("");
    m.AddPointData("v", std::move(v));
    return m;
}

}  // namespace

TEST(Colorize, ContinuousColoursFollowTheColormap) {
    ColorizeOptions o;
    o.mField = "v";
    const ColorizeResult r = colorize(col_line(), o);
    EXPECT_EQ(r.mLocation, "point");
    EXPECT_DOUBLE_EQ(r.mVMin, 0.0);
    EXPECT_DOUBLE_EQ(r.mVMax, 9.0);
    const std::uint8_t* red = r.mMesh.PointData(kColorRedName).As<std::uint8_t>();
    const std::uint8_t* green = r.mMesh.PointData(kColorGreenName).As<std::uint8_t>();
    const detail::Rgb lo = detail::colormap_lookup(detail::colormap_table("viridis"), 0.0);
    const detail::Rgb hi = detail::colormap_lookup(detail::colormap_table("viridis"), 1.0);
    EXPECT_EQ(red[0], lo.mR);
    EXPECT_EQ(green[9], hi.mG);
    EXPECT_EQ(red[10], 128);  // NaN: the nan colour
    EXPECT_EQ(r.mMesh.PointData(kColorRedName).Dtype(), DType::UInt8);
}

TEST(Colorize, BandsGiveFlatSteps) {
    ColorizeOptions o;
    o.mField = "v";
    o.mSteps = 3;
    const ColorizeResult r = colorize(col_line(), o);
    std::set<std::tuple<int, int, int>> seen;
    for (int i = 0; i < 10; ++i)
        seen.insert({r.mMesh.PointData(kColorRedName).As<std::uint8_t>()[i],
                     r.mMesh.PointData(kColorGreenName).As<std::uint8_t>()[i],
                     r.mMesh.PointData(kColorBlueName).As<std::uint8_t>()[i]});
    EXPECT_EQ(seen.size(), 3u);
}

TEST(Colorize, SymmetricRangeAndBadOptions) {
    ColorizeOptions o;
    o.mField = "v";
    o.mSymmetric = true;
    const ColorizeResult r = colorize(col_line(), o);
    EXPECT_DOUBLE_EQ(r.mVMin, -9.0);
    EXPECT_DOUBLE_EQ(r.mVMax, 9.0);
    o.mField = "nope";
    EXPECT_THROW(colorize(col_line(), o), std::invalid_argument);
    o.mField = "v";
    o.mCmap = "nope";
    EXPECT_THROW(colorize(col_line(), o), std::invalid_argument);
    o.mCmap = "viridis";
    o.mNanColor = "red";
    EXPECT_THROW(colorize(col_line(), o), std::invalid_argument);
}
