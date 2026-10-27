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
// Tests for error_size (roadmap 6.10.2): the Zienkiewicz-Zhu (1987) size
// field from a per-cell error indicator.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/error_size.hpp"
#include "meshioplusplus/region.hpp"

namespace {

namespace mio = meshioplusplus;

/// Two unit right triangles of the unit square, with indicator values.
mio::Mesh square(double E0, double E1) {
    mio::Mesh m;
    mio::NDArray p(mio::DType::Float64, {4, 3});
    const double x[12] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    std::copy(x, x + 12, p.As<double>());
    m.AssignPoints(std::move(p));
    mio::NDArray c(mio::DType::Int64, {2, 3});
    const std::int64_t t[6] = {0, 1, 2, 0, 2, 3};
    std::copy(t, t + 6, c.As<std::int64_t>());
    m.AddCellBlock("triangle", std::move(c));
    mio::NDArray e(mio::DType::Float64, {2});
    e.As<double>()[0] = E0;
    e.As<double>()[1] = E1;
    std::vector<mio::NDArray> blocks;
    blocks.push_back(std::move(e));
    m.AddCellData("error:zz", std::move(blocks));
    return m;
}

}  // namespace

// Each cell asks h_K (eta / sqrt(N) / eta_K)^(1/p); a point takes the least
// of its cells'. Diameter sqrt 2, N = 2, eta = 1: a cell with eta_K = 1 asks
// sqrt 2 / sqrt 2 = 1, one with eta_K = 4 a quarter (p = 1).
TEST(ErrorSize, TheZienkiewiczZhuSizes) {
    mio::ErrorSizeOptions o;
    o.mTargetError = 1.0;
    const auto r = mio::error_size(square(1.0, 4.0), o);
    EXPECT_EQ(r.mNumCells, 2);
    EXPECT_EQ(r.mNumSkipped, 0);
    const double* s = r.mMesh.PointData("size").As<double>();
    EXPECT_NEAR(s[1], 1.0, 1e-14);   // the first cell's alone
    EXPECT_NEAR(s[3], 0.25, 1e-14);  // the second's alone
    EXPECT_NEAR(s[0], 0.25, 1e-14);  // both: the least
    EXPECT_NEAR(r.mMinSize, 0.25, 1e-14);
    EXPECT_NEAR(r.mMaxSize, 1.0, 1e-14);
    // p = 2: sqrt 2 (1 / sqrt 2 / 4)^(1/2) = 2^(-3/4).
    o.mOrder = 2.0;
    EXPECT_NEAR(mio::error_size(square(1.0, 4.0), o).mMesh.PointData("size").As<double>()[3],
                std::pow(2.0, -0.75), 1e-14);
    // Bounds clamp and are counted.
    o.mOrder = 1.0;
    o.mMinSize = 0.5;
    o.mMaxSize = 0.8;
    const auto c = mio::error_size(square(1.0, 4.0), o);
    EXPECT_EQ(c.mNumClamped, 2);
    EXPECT_NEAR(c.mMesh.PointData("size").As<double>()[1], 0.8, 1e-14);
    EXPECT_NEAR(c.mMesh.PointData("size").As<double>()[3], 0.5, 1e-14);
}

// A cell with no finite indicator is skipped; a point only it uses reads NaN.
// A zero indicator asks for no bound, so the cell's size is the maximum or
// it is skipped. Regions and arrays pass through.
TEST(ErrorSize, SkippedCellsAndPassThrough) {
    mio::Mesh m = square(std::numeric_limits<double>::quiet_NaN(), 1.0);
    mio::NDArray idx(mio::DType::Int64, {1});
    idx.As<std::int64_t>()[0] = 0;
    m.AddRegion(mio::Region("first", mio::RegionKind::Cell, 2, -1, std::move(idx)));
    mio::ErrorSizeOptions o;
    o.mTargetError = 1.0;
    const auto r = mio::error_size(m, o);
    EXPECT_EQ(r.mNumCells, 1);
    EXPECT_EQ(r.mNumSkipped, 1);
    const double* s = r.mMesh.PointData("size").As<double>();
    EXPECT_TRUE(std::isnan(s[1]));
    EXPECT_FALSE(std::isnan(s[3]));
    EXPECT_EQ(r.mMesh.NumRegions(), 1u);
    EXPECT_TRUE(r.mMesh.HasCellData("error:zz"));
    const auto z = mio::error_size(square(0.0, 1.0), o);
    EXPECT_EQ(z.mNumSkipped, 1);
    o.mMaxSize = 3.0;
    const auto zb = mio::error_size(square(0.0, 1.0), o);
    EXPECT_EQ(zb.mNumSkipped, 0);
    EXPECT_NEAR(zb.mMesh.PointData("size").As<double>()[1], 3.0, 1e-14);
}

TEST(ErrorSize, RefusalsAreByName) {
    const mio::Mesh m = square(1.0, 1.0);
    mio::ErrorSizeOptions o;
    EXPECT_THROW((void)mio::error_size(m, o), std::invalid_argument);  // no target
    o.mTargetError = 1.0;
    o.mOrder = 0.0;
    EXPECT_THROW((void)mio::error_size(m, o), std::invalid_argument);
    o.mOrder = 1.0;
    o.mMinSize = 2.0;
    o.mMaxSize = 1.0;
    EXPECT_THROW((void)mio::error_size(m, o), std::invalid_argument);
    o.mMinSize = 0.0;
    o.mMaxSize = 0.0;
    o.mErrorArray = "nope";
    EXPECT_THROW((void)mio::error_size(m, o), std::invalid_argument);
    o.mErrorArray = "error:zz";
    const mio::Mesh once = mio::error_size(m, o).mMesh;
    EXPECT_THROW((void)mio::error_size(once, o), std::invalid_argument);  // "size" taken
    o.mOverwrite = true;
    EXPECT_NO_THROW((void)mio::error_size(once, o));
}
