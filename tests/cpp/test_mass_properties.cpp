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
// Mass, centre of mass and inertia (roadmap 9.4): closed forms on a box, a
// sphere's cell-count independent limit, the parallel-axis theorem, shells,
// regions and the refusals.

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
#include "meshioplusplus/operations/mass_properties.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::MassAbout;
using meshioplusplus::mass_properties;
using meshioplusplus::MassPropertiesOptions;
using meshioplusplus::MassPropertiesReport;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;

constexpr double kTol = 1e-9;

// A box [x, x+2] x [y, y+3] x [z, z+4].
std::vector<std::vector<double>> mp_box_points(double X, double Y, double Z) {
    std::vector<std::vector<double>> p;
    const double dx[8] = {0, 2, 2, 0, 0, 2, 2, 0};
    const double dy[8] = {0, 0, 3, 3, 0, 0, 3, 3};
    const double dz[8] = {0, 0, 0, 0, 4, 4, 4, 4};
    for (int v = 0; v < 8; ++v)
        p.push_back({X + dx[v], Y + dy[v], Z + dz[v]});
    return p;
}

Mesh mp_hex_box(double X = 0.0, double Y = 0.0, double Z = 0.0) {
    return mt::make_mesh(mp_box_points(X, Y, Z), "hexahedron",
                                              {{0, 1, 2, 3, 4, 5, 6, 7}});
}

// The same box as 12 outward triangles; Flip winds them inward.
Mesh mp_tri_box(bool Flip, double X = 0.0, double Y = 0.0, double Z = 0.0) {
    const std::vector<std::vector<std::int64_t>> quads = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                                          {3, 7, 6, 2}, {0, 4, 7, 3}, {1, 2, 6, 5}};
    std::vector<std::vector<std::int64_t>> tris;
    for (const auto& q : quads) {
        tris.push_back({q[0], q[1], q[2]});
        tris.push_back({q[0], q[2], q[3]});
    }
    if (Flip)
        for (auto& t : tris)
            std::swap(t[1], t[2]);
    return mt::make_mesh(mp_box_points(X, Y, Z), "triangle", tris);
}

void expect_box(const MassPropertiesReport& rR, double Rho, double X, double Y, double Z) {
    const auto& t = rR.mTotal;
    EXPECT_NEAR(t.mMeasure, 24.0, kTol * 24);
    EXPECT_NEAR(t.mMass, 24.0 * Rho, kTol * 24 * Rho);
    EXPECT_NEAR(t.mCentre[0], X + 1.0, kTol);
    EXPECT_NEAR(t.mCentre[1], Y + 1.5, kTol);
    EXPECT_NEAR(t.mCentre[2], Z + 2.0, kTol);
    const double m = t.mMass;
    EXPECT_NEAR(t.mInertia[0], m * (9.0 + 16.0) / 12.0, kTol * m * 3);
    EXPECT_NEAR(t.mInertia[1], m * (4.0 + 16.0) / 12.0, kTol * m * 3);
    EXPECT_NEAR(t.mInertia[2], m * (4.0 + 9.0) / 12.0, kTol * m * 3);
    EXPECT_NEAR(t.mInertia[3], 0.0, kTol * m * 3);
    EXPECT_NEAR(t.mInertia[4], 0.0, kTol * m * 3);
    EXPECT_NEAR(t.mInertia[5], 0.0, kTol * m * 3);
    // Ascending principal moments of a box with sides 2, 3, 4.
    EXPECT_NEAR(t.mPrincipalMoments[0], m * 13.0 / 12.0, kTol * m * 3);
    EXPECT_NEAR(t.mPrincipalMoments[2], m * 25.0 / 12.0, kTol * m * 3);
    EXPECT_NEAR(t.mGyrationRadii[0], std::sqrt(13.0 / 12.0), 1e-8);
}

}  // namespace

TEST(MassProperties, HexahedronMatchesBoxClosedForm) {
    const auto r = mass_properties(mp_hex_box(), {.mDensity = 2.0});
    EXPECT_EQ(r.mKind, "solid");
    EXPECT_EQ(r.mTotal.mNumCells, 1);
    expect_box(r, 2.0, 0.0, 0.0, 0.0);
}

TEST(MassProperties, FarFromTheOriginKeepsAccuracy) {
    const auto r = mass_properties(mp_hex_box(1e6, -2e6, 3e6), {.mDensity = 7850.0});
    expect_box(r, 7850.0, 1e6, -2e6, 3e6);
}

TEST(MassProperties, ClosedTriangleSurfaceEqualsTheSolid) {
    expect_box(mass_properties(mp_tri_box(false), {.mDensity = 2.0}), 2.0, 0.0, 0.0, 0.0);
}

TEST(MassProperties, InwardWoundSurfaceIsTurnedOutward) {
    const auto r = mass_properties(mp_tri_box(true), {.mDensity = 2.0});
    EXPECT_GT(r.mTotal.mMass, 0.0);
    expect_box(r, 2.0, 0.0, 0.0, 0.0);
}

TEST(MassProperties, ParallelAxisTheoremAboutTheOrigin) {
    MassPropertiesOptions o;
    o.mDensity = 1.5;
    const auto c = mass_properties(mp_hex_box(1.0, 2.0, 3.0), o);
    o.mAbout = MassAbout::Origin;
    const auto g = mass_properties(mp_hex_box(1.0, 2.0, 3.0), o);
    const double m = c.mTotal.mMass;
    const std::array<double, 3> d = c.mTotal.mCentre;
    const double d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    EXPECT_NEAR(g.mTotal.mInertia[0], c.mTotal.mInertia[0] + m * (d2 - d[0] * d[0]), 1e-8);
    EXPECT_NEAR(g.mTotal.mInertia[1], c.mTotal.mInertia[1] + m * (d2 - d[1] * d[1]), 1e-8);
    EXPECT_NEAR(g.mTotal.mInertia[2], c.mTotal.mInertia[2] + m * (d2 - d[2] * d[2]), 1e-8);
    EXPECT_NEAR(g.mTotal.mInertia[3], c.mTotal.mInertia[3] - m * d[0] * d[1], 1e-8);
    EXPECT_NEAR(g.mTotal.mInertia[4], c.mTotal.mInertia[4] - m * d[1] * d[2], 1e-8);
    EXPECT_NEAR(g.mTotal.mInertia[5], c.mTotal.mInertia[5] - m * d[2] * d[0], 1e-8);
    EXPECT_NEAR(g.mTotal.mCentre[0], c.mTotal.mCentre[0], 1e-12);
}

TEST(MassProperties, TetrahedronVolumeAndCentre) {
    const Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                                                      "tetra", {{0, 1, 2, 3}});
    const auto r = mass_properties(m, {.mDensity = 6.0});
    EXPECT_NEAR(r.mTotal.mMeasure, 1.0 / 6.0, kTol);
    EXPECT_NEAR(r.mTotal.mMass, 1.0, kTol);
    EXPECT_NEAR(r.mTotal.mCentre[0], 0.25, kTol);
    // About the origin: Sxx = m/10 for the unit corner tetrahedron.
    MassPropertiesOptions o;
    o.mDensity = 6.0;
    o.mAbout = MassAbout::Origin;
    const auto g = mass_properties(m, o);
    EXPECT_NEAR(g.mTotal.mInertia[0], 2.0 / 10.0, kTol);
    EXPECT_NEAR(g.mTotal.mInertia[3], -1.0 / 20.0, kTol);
}

TEST(MassProperties, ShellOfAPlate) {
    const Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}},
                                                      "quad", {{0, 1, 2, 3}});
    MassPropertiesOptions o;
    o.mDensity = 8.0;
    o.mThickness = 0.5;
    const auto r = mass_properties(m, o);
    EXPECT_EQ(r.mKind, "shell");
    EXPECT_NEAR(r.mTotal.mMeasure, 1.0, kTol);
    EXPECT_NEAR(r.mTotal.mMass, 4.0, kTol);
    EXPECT_NEAR(r.mTotal.mInertia[2], 4.0 / 6.0, kTol);  // m (a^2 + b^2) / 12
    EXPECT_NEAR(r.mTotal.mInertia[0], 4.0 / 12.0, kTol);
}

TEST(MassProperties, OpenSurfaceWithoutThicknessIsRefused) {
    const Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}},
                                                      "quad", {{0, 1, 2, 3}});
    try {
        mass_properties(m);
        FAIL() << "expected an open surface to be refused";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("closed"), std::string::npos);
    }
}

TEST(MassProperties, DensityArrayAndRegions) {
    Mesh m = mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
         {0, 0, 2}, {1, 0, 2}, {1, 1, 2}, {0, 1, 2}},
        "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}, {4, 5, 6, 7, 8, 9, 10, 11}});
    m.AddCellData("rho", {mt::data_array({1.0, 3.0})});
    NDArray ids = NDArray::Uninit(DType::Int64, {std::size_t{1}});
    ids.As<std::int64_t>()[0] = 1;
    m.AddRegion(Region("top", RegionKind::Cell, ids));
    MassPropertiesOptions o;
    o.mDensityArray = "rho";
    const auto r = mass_properties(m, o);
    EXPECT_NEAR(r.mTotal.mMass, 4.0, kTol);
    EXPECT_NEAR(r.mTotal.mCentre[2], (1.0 * 0.5 + 3.0 * 1.5) / 4.0, kTol);
    ASSERT_EQ(r.mRegions.size(), 1u);
    EXPECT_EQ(r.mRegions[0].mName, "top");
    EXPECT_NEAR(r.mRegions[0].mMass, 3.0, kTol);
    EXPECT_NEAR(r.mRegions[0].mCentre[2], 1.5, kTol);
    // A region option restricts the total and lists no sub-regions.
    MassPropertiesOptions only;
    only.mRegion = "top";
    const auto s = mass_properties(m, only);
    EXPECT_NEAR(s.mTotal.mMass, 1.0, kTol);
    EXPECT_TRUE(s.mRegions.empty());
}

TEST(MassProperties, Refusals) {
    const Mesh box = mp_hex_box();
    EXPECT_THROW(mass_properties(box, {.mDensityArray = "nope"}), std::invalid_argument);
    EXPECT_THROW(mass_properties(box, {.mRegion = "nope"}), std::invalid_argument);
    EXPECT_THROW(mass_properties(box, {.mThickness = -1.0}), std::invalid_argument);
    EXPECT_THROW(mass_properties(box, {.mDensity = std::nan("")}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::mass_about_from_name("moon"), std::invalid_argument);
    Mesh lines = mt::make_mesh({{0, 0, 0}, {1, 0, 0}}, "line", {{0, 1}});
    EXPECT_THROW(mass_properties(lines), std::invalid_argument);
}
