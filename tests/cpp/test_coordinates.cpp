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
// Cylindrical and spherical components about an axis (roadmap 7.6.1).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/coordinates.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// Points around the z axis carrying a velocity of a rigid rotation about it
// (omega x r): purely tangential, so (0, |r|, 0) in cylindrical components.
Mesh co_ring() {
    std::vector<std::vector<double>> pts;
    for (int k = 0; k < 8; ++k) {
        const double t = k * M_PI / 4;
        pts.push_back({2.0 * std::cos(t), 2.0 * std::sin(t), 0.5 * k});
    }
    Mesh m = mt::make_mesh(pts, "line", {{0, 1}, {2, 3}, {4, 5}, {6, 7}});
    NDArray v(DType::Float64, {8, 3});
    for (std::size_t i = 0; i < 8; ++i) {
        const double x = pts[i][0], y = pts[i][1];
        v.As<double>()[i * 3 + 0] = -y;
        v.As<double>()[i * 3 + 1] = x;
        v.As<double>()[i * 3 + 2] = 0.0;
    }
    m.AddPointData("v", std::move(v));
    return m;
}

}  // namespace

TEST(Coordinates, RigidRotationIsPurelyTangential) {
    CoordinateOptions o;
    o.mArrays = {"v"};
    const CoordinateResult r = convert_coordinates(co_ring(), o);
    EXPECT_EQ(r.mNumArrays, 1);
    EXPECT_EQ(r.mNumSingular, 0);
    const double* v = r.mMesh.PointData("v").As<double>();
    const double* c = r.mMesh.PointData(kCylindricalCoordsName).As<double>();
    for (std::size_t i = 0; i < 8; ++i) {
        EXPECT_NEAR(v[i * 3 + 0], 0.0, 1e-12);
        EXPECT_NEAR(v[i * 3 + 1], 2.0, 1e-12);
        EXPECT_NEAR(v[i * 3 + 2], 0.0, 1e-12);
        EXPECT_NEAR(c[i * 3 + 0], 2.0, 1e-12);
        EXPECT_NEAR(c[i * 3 + 2], 0.5 * static_cast<double>(i), 1e-12);
    }
    EXPECT_NEAR(c[1 * 3 + 1], M_PI / 4, 1e-12);
}

TEST(Coordinates, ForwardThenInverseRoundTrips) {
    for (const CoordinateSystem sys :
         {CoordinateSystem::Cylindrical, CoordinateSystem::Spherical}) {
        CoordinateOptions o;
        o.mSystem = sys;
        o.mAxis = {1.0, 2.0, 0.5};
        o.mOrigin = {0.3, -0.1, 0.2};
        o.mArrays = {"v"};
        const Mesh m = co_ring();
        const CoordinateResult fwd = convert_coordinates(m, o);
        o.mInverse = true;
        o.mRecordCoordinates = false;
        const CoordinateResult back = convert_coordinates(fwd.mMesh, o);
        for (std::size_t k = 0; k < 24; ++k)
            EXPECT_NEAR(back.mMesh.PointData("v").As<double>()[k], m.PointData("v").As<double>()[k],
                        1e-12);
    }
}

TEST(Coordinates, Refusals) {
    CoordinateOptions o;
    o.mAxis = {0.0, 0.0, 0.0};
    EXPECT_THROW(convert_coordinates(co_ring(), o), std::invalid_argument);
    o.mAxis = {0.0, 0.0, 1.0};
    o.mArrays = {"nope"};
    EXPECT_THROW(convert_coordinates(co_ring(), o), std::invalid_argument);
    EXPECT_THROW(coordinate_system_from_name("polar"), std::invalid_argument);
}
