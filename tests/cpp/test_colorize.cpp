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
#include "meshioplusplus/detail/value_io.hpp"
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
    // read_int: the NATIVE backend stores every integer array as Int64.
    const NDArray& red = r.mMesh.PointData(kColorRedName);
    const NDArray& green = r.mMesh.PointData(kColorGreenName);
    const detail::Rgb lo = detail::colormap_lookup(detail::colormap_table("viridis"), 0.0);
    const detail::Rgb hi = detail::colormap_lookup(detail::colormap_table("viridis"), 1.0);
    EXPECT_EQ(detail::read_int(red, 0), lo.mR);
    EXPECT_EQ(detail::read_int(green, 9), hi.mG);
    EXPECT_EQ(detail::read_int(red, 10), 128);  // NaN: the nan colour
#if defined(MESHIOPLUSPLUS_MESH_BACKEND_MESHIO)
    EXPECT_EQ(red.Dtype(), DType::UInt8);
#endif
}

TEST(Colorize, BandsGiveFlatSteps) {
    ColorizeOptions o;
    o.mField = "v";
    o.mSteps = 3;
    const ColorizeResult r = colorize(col_line(), o);
    std::set<std::tuple<int, int, int>> seen;
    for (int i = 0; i < 10; ++i)
        seen.insert({static_cast<int>(detail::read_int(r.mMesh.PointData(kColorRedName), i)),
                     static_cast<int>(detail::read_int(r.mMesh.PointData(kColorGreenName), i)),
                     static_cast<int>(detail::read_int(r.mMesh.PointData(kColorBlueName), i))});
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
