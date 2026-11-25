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
// Ray casting (roadmap 9.3): analytic sphere and box, seams, limits, regions
// and the refusals.

// System includes
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/cast_rays.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::cast_rays;
using meshioplusplus::CastRaysOptions;
using meshioplusplus::Mesh;

// The unit cube [0, 1]^3 as 12 outward-wound triangles.
Mesh cr_cube() {
    const std::vector<std::vector<std::int64_t>> tris = {
        {0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
        {1, 2, 6}, {1, 6, 5}, {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    return mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1},
                          {1, 1, 1}, {0, 1, 1}},
                         "triangle", tris);
}

}  // namespace

TEST(CastRays, HitsTheFacesOfACubeFromOutside) {
    const Mesh cube = cr_cube();
    // From the left, the right, below and above, towards the centre of a face.
    const std::vector<double> o = {-1, 0.5, 0.5, 2, 0.5, 0.5, 0.25, 0.25, -3, 0.5, 0.5, 5};
    const std::vector<double> d = {1, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, -1};
    const auto r = cast_rays(cube, o, d);
    ASSERT_EQ(r.mHit.size(), 4u);
    EXPECT_EQ(r.mNumHits, 4);
    const double dist[4] = {1.0, 1.0, 3.0, 4.0};
    const double nx[4][3] = {{-1, 0, 0}, {1, 0, 0}, {0, 0, -1}, {0, 0, 1}};
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(r.mHit[i], 1) << i;
        EXPECT_NEAR(r.mDistance[i], dist[i], 1e-12) << i;
        for (int k = 0; k < 3; ++k)
            EXPECT_NEAR(r.mNormal[3 * i + k], nx[i][k], 1e-12) << i << "," << k;
        EXPECT_GE(r.mCell[i], 0);
        EXPECT_LT(r.mCell[i], 12);
    }
    EXPECT_NEAR(r.mPoint[0], 0.0, 1e-12);
    EXPECT_NEAR(r.mPoint[3], 1.0, 1e-12);
}

TEST(CastRays, SingleDirectionBroadcastsAndNormalises) {
    const Mesh cube = cr_cube();
    const std::vector<double> o = {-2, 0.25, 0.25, -2, 0.75, 0.5};
    const auto r = cast_rays(cube, o, {5, 0, 0});  // not a unit vector
    EXPECT_EQ(r.mNumHits, 2);
    EXPECT_NEAR(r.mDistance[0], 2.0, 1e-12);
    EXPECT_NEAR(r.mDistance[1], 2.0, 1e-12);
}

TEST(CastRays, MissesAndTheDistanceLimit) {
    const Mesh cube = cr_cube();
    const auto miss = cast_rays(cube, {-1, 5, 5}, {1, 0, 0});
    EXPECT_EQ(miss.mHit[0], 0);
    EXPECT_TRUE(std::isnan(miss.mDistance[0]));
    EXPECT_TRUE(std::isnan(miss.mPoint[0]));
    EXPECT_EQ(miss.mCell[0], -1);
    CastRaysOptions near;
    near.mMaxDistance = 0.5;
    EXPECT_EQ(cast_rays(cube, {-1, 0.5, 0.5}, {1, 0, 0}, near).mNumHits, 0);
    near.mMaxDistance = 1.0;  // the limit is inclusive
    EXPECT_EQ(cast_rays(cube, {-1, 0.5, 0.5}, {1, 0, 0}, near).mNumHits, 1);
    // From inside, the nearest hit is the wall in front.
    const auto inside = cast_rays(cube, {0.5, 0.5, 0.5}, {0, 1, 0});
    EXPECT_NEAR(inside.mDistance[0], 0.5, 1e-12);
    // A ray starting on the surface does not hit it at t = 0.
    const auto on = cast_rays(cube, {0, 0.5, 0.5}, {1, 0, 0});
    EXPECT_NEAR(on.mDistance[0], 1.0, 1e-12);
}

TEST(CastRays, SeamsAreHitOnce) {
    const Mesh cube = cr_cube();
    // Through the diagonal shared by the two triangles of the face x = 0, then
    // through the vertex (0, 0, 0) itself and along a cube edge.
    const auto r = cast_rays(cube, {-1, 0.5, 0.5, -1, 0, 0, -1, 0.5, 0}, {1, 0, 0});
    EXPECT_EQ(r.mHit[0], 1);
    EXPECT_NEAR(r.mDistance[0], 1.0, 1e-9);
    EXPECT_EQ(r.mHit[1], 1);
    EXPECT_NEAR(r.mDistance[1], 1.0, 1e-9);
    EXPECT_EQ(r.mHit[2], 1);
    EXPECT_NEAR(r.mDistance[2], 1.0, 1e-9);
}

TEST(CastRays, SphereAgainstTheClosedForm) {
    const Mesh sphere = meshioplusplus::uv_sphere(1.0, 64, 32);
    // Rays from outside towards the centre hit at about the radius.
    std::vector<double> o, d;
    const double angles[4][2] = {{0.3, 0.4}, {1.1, 2.0}, {2.5, 0.9}, {4.0, 1.5}};
    for (const auto& a : angles) {
        const double x = std::sin(a[1]) * std::cos(a[0]);
        const double y = std::sin(a[1]) * std::sin(a[0]);
        const double z = std::cos(a[1]);
        o.insert(o.end(), {3 * x, 3 * y, 3 * z});
        d.insert(d.end(), {-x, -y, -z});
    }
    const auto r = cast_rays(sphere, o, d);
    EXPECT_EQ(r.mNumHits, 4);
    for (int i = 0; i < 4; ++i) {
        // A facet is inscribed in the sphere, so the hit is a little inside it.
        EXPECT_GT(r.mDistance[i], 1.99);
        EXPECT_LT(r.mDistance[i], 2.01);
        const double nx = r.mNormal[3 * i], ny = r.mNormal[3 * i + 1], nz = r.mNormal[3 * i + 2];
        EXPECT_NEAR(nx * nx + ny * ny + nz * nz, 1.0, 1e-12);
        // The outward normal faces back towards the ray's origin.
        EXPECT_LT(nx * d[3 * i] + ny * d[3 * i + 1] + nz * d[3 * i + 2], 0.0);
    }
}

TEST(CastRays, RegionRestrictsTheSurface) {
    Mesh cube = cr_cube();
    meshioplusplus::NDArray ids =
        meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Int64, {std::size_t{2}});
    ids.As<std::int64_t>()[0] = 10;  // the two triangles of the face x = 0
    ids.As<std::int64_t>()[1] = 11;
    cube.AddRegion(meshioplusplus::Region("left", meshioplusplus::RegionKind::Cell, ids));
    CastRaysOptions o;
    o.mRegion = "left";
    // From the right, the face x = 1 is not in the region: the ray reaches x = 0.
    const auto r = cast_rays(cube, {2, 0.5, 0.5}, {-1, 0, 0}, o);
    EXPECT_NEAR(r.mDistance[0], 2.0, 1e-12);
    EXPECT_GE(r.mCell[0], 10);
    o.mRegion = "nope";
    EXPECT_THROW(cast_rays(cube, {0, 0, 0}, {1, 0, 0}, o), std::invalid_argument);
}

TEST(CastRays, Refusals) {
    const Mesh cube = cr_cube();
    EXPECT_THROW(cast_rays(cube, {0, 0}, {1, 0, 0}), std::invalid_argument);
    EXPECT_THROW(cast_rays(cube, {0, 0, 0}, {}), std::invalid_argument);
    EXPECT_THROW(cast_rays(cube, {0, 0, 0}, {0, 0, 0}), std::invalid_argument);
    EXPECT_THROW(cast_rays(cube, {0, 0, 0, 1, 1, 1}, {1, 0, 0, 0, 1, 0, 0, 0, 1}),
                 std::invalid_argument);
    EXPECT_THROW(cast_rays(cube, {std::nan(""), 0, 0}, {1, 0, 0}), std::invalid_argument);
    const Mesh tets = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra",
                                    {{0, 1, 2, 3}});
    EXPECT_THROW(cast_rays(tets, {2, 2, 2}, {-1, -1, -1}), std::invalid_argument);
    EXPECT_TRUE(cast_rays(cube, {}, {1, 0, 0}).mHit.empty());
}
