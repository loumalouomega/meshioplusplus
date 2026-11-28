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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/operations/subsample_points.hpp"

using namespace meshioplusplus;

namespace {

Mesh make_grid_cloud(std::size_t nx, std::size_t ny) {
    const std::size_t n = nx * ny;
    NDArray pts = NDArray::Uninit(DType::Float64, {n, std::size_t{3}});
    auto* p = pts.As<double>();
    NDArray scalar = NDArray::Uninit(DType::Float64, {n});
    auto* s = scalar.As<double>();
    std::size_t k = 0;
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            p[3 * k + 0] = static_cast<double>(i);
            p[3 * k + 1] = static_cast<double>(j);
            p[3 * k + 2] = 0.0;
            s[k] = static_cast<double>(k) * 2.0;
            ++k;
        }
    }
    Mesh m;
    m.AssignPoints(std::move(pts));
    m.AddPointData("u", std::move(scalar));
    return m;
}

TEST(SubsamplePoints, FarthestGridAndRandomExactCount) {
    Mesh cloud = make_grid_cloud(6, 5);  // 30 points
    for (SubsampleMethod method :
         {SubsampleMethod::Farthest, SubsampleMethod::Grid, SubsampleMethod::Random}) {
        SubsamplePointsOptions opts;
        opts.mCount = 10;
        opts.mMethod = method;
        opts.mSeed = 7;
        opts.mRecordIds = true;
        auto res = subsample_points(cloud, opts);
        EXPECT_EQ(res.mMesh.NumPoints(), 10u);
        EXPECT_EQ(res.mPointMap.size(), 10u);
        std::set<std::int64_t> uniq(res.mPointMap.begin(), res.mPointMap.end());
        EXPECT_EQ(uniq.size(), 10u);
        ASSERT_TRUE(res.mMesh.HasPointData("u"));
        ASSERT_TRUE(res.mMesh.HasPointData("budget:original_point_id"));
        const auto* u = res.mMesh.PointData("u").As<double>();
        const auto* orig = res.mMesh.PointData("budget:original_point_id").As<std::int64_t>();
        for (std::size_t i = 0; i < 10; ++i) {
            EXPECT_EQ(orig[i], res.mPointMap[i]);
            EXPECT_DOUBLE_EQ(u[i], static_cast<double>(res.mPointMap[i]) * 2.0);
        }
    }
}

TEST(SubsamplePoints, BoundsFilteringAndValidation) {
    Mesh cloud = make_grid_cloud(5, 5);  // x,y in [0, 4]
    SubsamplePointsOptions opts;
    opts.mCount = 4;
    opts.mMethod = SubsampleMethod::Farthest;
    opts.mHasBounds = true;
    opts.mBounds = {0.0, 0.0, -1.0, 2.0, 2.0, 1.0};  // 3x3 = 9 candidates
    opts.mStart = 0;
    auto res = subsample_points(cloud, opts);
    EXPECT_EQ(res.mMesh.NumPoints(), 4u);
    const auto* pts = res.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_LE(pts[3 * i + 0], 2.0);
        EXPECT_LE(pts[3 * i + 1], 2.0);
    }

    opts.mCount = 0;
    EXPECT_THROW(subsample_points(cloud, opts), std::invalid_argument);
}

}  // namespace
