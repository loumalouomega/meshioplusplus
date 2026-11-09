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
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/thickness.hpp"

using namespace meshioplusplus;

namespace {

// The cells of the (one-block) surface wound the other way round.
Mesh thk_flipped(const Mesh& rMesh) {
    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rMesh.Points()));
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cells = rMesh.Cells(b);
        const std::size_t n = cells.NumCells(), k = cells.NodesPerCell();
        NDArray conn = detail::data_owned_copy(cells.Conn());
        std::int64_t* p = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < n; ++i)
            std::reverse(p + i * k + 1, p + (i + 1) * k);
        out.AddCellBlock(cells.Type(), std::move(conn));
    }
    return out;
}

const NDArray& thk_dist(const ThicknessResult& rResult, std::size_t Block = 0) {
    return rResult.mMesh.CellData(kThicknessDistanceName, Block);
}

}  // namespace

TEST(Thickness, HollowBoxWallIsOneEverywhere) {
    const Mesh outer = extract_surface(box({{1, 1, 1}}, {{3.0, 3.0, 3.0}}));
    const Mesh inner =
        thk_flipped(extract_surface(box({{1, 1, 1}}, {{1.0, 1.0, 1.0}}, {{1.0, 1.0, 1.0}})));
    const Mesh shell = merge({&outer, &inner}).mMesh;
    const ThicknessResult r = thickness(shell);
    EXPECT_EQ(r.mDirection, ThicknessDirection::Inward);
    EXPECT_EQ(r.mNumCells, 12);
    EXPECT_EQ(r.mNumHit, 12);
    EXPECT_EQ(r.mNumMissed, 0);
    EXPECT_NEAR(r.mMin, 1.0, 1e-12);
    EXPECT_NEAR(r.mMax, 1.0, 1e-12);
    EXPECT_NEAR(r.mMean, 1.0, 1e-12);
    for (std::size_t b = 0; b < r.mMesh.NumCellBlocks(); ++b) {
        const NDArray& d = thk_dist(r, b);
        const NDArray& h = r.mMesh.CellData(kThicknessHitCellName, b);
        for (std::size_t i = 0; i < d.Size(); ++i) {
            EXPECT_NEAR(d.As<double>()[i], 1.0, 1e-12);
            EXPECT_GE(h.As<std::int64_t>()[i], 0);
        }
    }
    // A threshold counts the thin cells.
    ThicknessOptions o;
    o.mThreshold = 1.5;
    EXPECT_EQ(thickness(shell, o).mNumBelow, 12);
    // Outward from the outer box and into the cavity: the cavity's far wall,
    // or nothing.
    o = {};
    o.mDirection = ThicknessDirection::Outward;
    const ThicknessResult gap = thickness(shell, o);
    EXPECT_EQ(gap.mNumHit, 6);
    EXPECT_NEAR(gap.mMin, 1.0, 1e-12);
}

TEST(Thickness, SphereIsAsThickAsItsDiameter) {
    const ThicknessResult r = thickness(icosa_sphere(1.0, 3));
    EXPECT_EQ(r.mNumMissed, 0);
    EXPECT_NEAR(r.mMin, 2.0, 0.02);
    EXPECT_LE(r.mMax, 2.0);
    EXPECT_TRUE(r.mWarnings.empty());
}

TEST(Thickness, ClearanceToASecondPart) {
    const Mesh a = icosa_sphere(1.0, 3), b = icosa_sphere(1.0, 3, {{3.0, 0.0, 0.0}});
    const ThicknessResult r = thickness(a, b);
    EXPECT_EQ(r.mDirection, ThicknessDirection::Outward);
    EXPECT_GT(r.mNumHit, 0);
    EXPECT_GT(r.mNumMissed, 0);
    // Facet centroids lie inside the spheres: a little more than 1.
    EXPECT_NEAR(r.mMin, 1.0, 0.05);
    EXPECT_GE(r.mMin, 1.0);
    EXPECT_GT(r.mMinPoint[0], 0.9);
    // Bounded rays miss.
    ThicknessOptions o;
    o.mMaxDistance = 0.5;
    EXPECT_EQ(thickness(a, b, o).mNumHit, 0);
    // Deterministic.
    const ThicknessResult again = thickness(a, b);
    EXPECT_EQ(again.mMinCell, r.mMinCell);
    EXPECT_EQ(again.mMean, r.mMean);
}

TEST(Thickness, OpenSheetMissesAndBadArguments) {
    const ThicknessResult r = thickness(rectangle({{2, 2}}));
    EXPECT_EQ(r.mNumHit, 0);
    EXPECT_EQ(r.mNumMissed, 4);
    EXPECT_TRUE(std::isnan(thk_dist(r).As<double>()[0]));
    ThicknessOptions o;
    o.mMaxDistance = -1.0;
    EXPECT_THROW(thickness(rectangle(), o), std::invalid_argument);
    EXPECT_THROW(thickness(box()), std::invalid_argument);
    // Lines are not measured and read as misses.
    const ThicknessResult lines = thickness(mt::line_mesh());
    EXPECT_EQ(lines.mNumCells, 0);
    EXPECT_FALSE(lines.mWarnings.empty());
}

// ---------------------------------------------------------------------------
// Beyond one ray (roadmap 7.5): the sphere and cone methods, scan-safe hits.
// ---------------------------------------------------------------------------

namespace {

constexpr double kThkTestPi = 3.14159265358979323846;

// A closed wedge whose bottom (z = 0, x in [0, L], y in [0, W]) and top
// (z = x tan(Angle)) meet along the y axis; the bottom is K x K quads split
// into triangles, the rest coarse. Outward normals.
Mesh thk_wedge(double L, double W, double Angle, int K) {
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> tris;
    const auto add = [&](std::vector<double> p) {
        pts.push_back(std::move(p));
        return static_cast<std::int64_t>(pts.size() - 1);
    };
    std::vector<std::int64_t> grid;
    for (int i = 0; i <= K; ++i)
        for (int j = 0; j <= K; ++j)
            grid.push_back(add({L * i / K, W * j / K, 0.0}));
    const auto g = [&](int i, int j) { return grid[static_cast<std::size_t>(i * (K + 1) + j)]; };
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < K; ++j) {  // normal -z
            tris.push_back({g(i, j), g(i + 1, j + 1), g(i + 1, j)});
            tris.push_back({g(i, j), g(i, j + 1), g(i + 1, j + 1)});
        }
    const double h = L * std::tan(Angle);
    const std::int64_t a0 = add({0, 0, 0}), b0 = add({L, 0, 0}), c0 = add({L, 0, h});
    const std::int64_t a1 = add({0, W, 0}), b1 = add({L, W, 0}), c1 = add({L, W, h});
    tris.push_back({a0, c1, a1});  // top
    tris.push_back({a0, c0, c1});
    tris.push_back({b0, b1, c1});  // end, +x
    tris.push_back({b0, c1, c0});
    tris.push_back({a0, b0, c0});  // y = 0, -y
    tris.push_back({a1, c1, b1});  // y = W, +y
    return mt::make_mesh(pts, "triangle", tris);
}

// The measured distance of the bottom facet nearest (X, Y), and its centroid x.
std::pair<double, double> thk_at(const ThicknessResult& rR, const Mesh& rMesh, double X, double Y) {
    const auto cells = rMesh.Cells(0);
    const double* p = rMesh.Points().As<double>();
    double best = 1e300, value = 0.0, cx = 0.0;
    for (std::size_t c = 0; c + 6 < cells.NumCells(); ++c) {
        double x = 0.0, y = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            const auto v = static_cast<std::size_t>(cells.Conn().As<std::int64_t>()[3 * c + k]);
            x += p[3 * v] / 3.0;
            y += p[3 * v + 1] / 3.0;
        }
        const double d = (x - X) * (x - X) + (y - Y) * (y - Y);
        if (d < best) {
            best = d;
            value = thk_dist(rR).As<double>()[c];
            cx = x;
        }
    }
    return {value, cx};
}

}  // namespace

TEST(ThicknessMethods, AWedgeReadsTheSphereValueAndTheLongerRay) {
    const double angle = 30.0 * kThkTestPi / 180.0;
    const Mesh wedge = thk_wedge(4.0, 8.0, angle, 40);
    const auto [ray, x] = thk_at(thickness(wedge), wedge, 2.0, 4.0);
    const double perpendicular = x * std::sin(angle);  // to the top face
    EXPECT_NEAR(ray, x * std::tan(angle), 1e-9);
    EXPECT_NEAR(ray, perpendicular / std::cos(angle), 1e-9);
    ThicknessOptions o;
    o.mMethod = ThicknessMethod::Sphere;
    const ThicknessResult sr = thickness(wedge, o);
    const double sphere = thk_at(sr, wedge, 2.0, 4.0).first;
    EXPECT_NEAR(sphere, perpendicular, 0.02 * perpendicular);
    EXPECT_LE(sphere, ray);
    EXPECT_EQ(sr.mMethod, ThicknessMethod::Sphere);
    // A cone fan finds the narrow spot the straight ray misses.
    o.mMethod = ThicknessMethod::Cone;
    o.mConeAngleDeg = 35.0;
    o.mConeRays = 128;
    const double cone = thk_at(thickness(wedge, o), wedge, 2.0, 4.0).first;
    EXPECT_GE(cone, perpendicular - 1e-9);
    EXPECT_LT(cone, 1.02 * perpendicular);
    EXPECT_LT(cone, ray);
    // A face more edge-on than the limit is stepped over: the ray then misses.
    o = {};
    o.mMaxNormalDeviationDeg = 20.0;
    EXPECT_TRUE(std::isnan(thk_at(thickness(wedge, o), wedge, 2.0, 4.0).first));
    o.mMaxNormalDeviationDeg = 40.0;
    EXPECT_NEAR(thk_at(thickness(wedge, o), wedge, 2.0, 4.0).first, ray, 1e-12);
}

TEST(ThicknessMethods, TwoSidedReadsAWallWhoseFacetsAreFlipped) {
    const Mesh outer = extract_surface(box({{4, 4, 4}}, {{4.0, 4.0, 4.0}}));
    const Mesh inner =
        thk_flipped(extract_surface(box({{2, 2, 2}}, {{2.0, 2.0, 2.0}}, {{1.0, 1.0, 1.0}})));
    const Mesh shell = merge({&outer, &inner}).mMesh;
    // Flip every third cell, as an unreliable scan winding would.
    Mesh scan;
    scan.AssignPoints(detail::data_owned_copy(shell.Points()));
    for (std::size_t b = 0; b < shell.NumCellBlocks(); ++b) {
        const auto cells = shell.Cells(b);
        const std::size_t k = cells.NodesPerCell();
        NDArray conn = detail::data_owned_copy(cells.Conn());
        std::int64_t* p = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < cells.NumCells(); i += 3)
            std::reverse(p + i * k + 1, p + (i + 1) * k);
        scan.AddCellBlock(cells.Type(), std::move(conn));
    }
    const auto count_ones = [](const ThicknessResult& rR) {
        std::int64_t n = 0;
        for (std::size_t b = 0; b < rR.mMesh.NumCellBlocks(); ++b) {
            const NDArray& d = thk_dist(rR, b);
            for (std::size_t i = 0; i < d.Size(); ++i)
                n += std::abs(d.As<double>()[i] - 1.0) < 1e-9 ? 1 : 0;
        }
        return n;
    };
    // The facets facing the inner box read the wall, 1; the rest see no wall
    // within a fifth of the diagonal.
    const std::int64_t walls = count_ones(thickness(shell));
    EXPECT_EQ(walls, 48);
    const ThicknessResult plain = thickness(scan);
    EXPECT_LT(count_ones(plain), walls);
    ThicknessOptions o;
    o.mTwoSided = true;
    const ThicknessResult two = thickness(scan, o);
    EXPECT_EQ(count_ones(two), walls);
    EXPECT_EQ(two.mNumHit, walls);
    EXPECT_NEAR(two.mMaxDistance, 0.2 * std::sqrt(48.0), 1e-12);
    EXPECT_NEAR(two.mP2, 1.0, 1e-9);
    EXPECT_NEAR(two.mP95, 1.0, 1e-9);
}

TEST(ThicknessMethods, DefaultsAreUnchangedAndBadOptionsAreRefused) {
    const Mesh sphere = icosa_sphere(1.0, 3);
    const ThicknessResult r = thickness(sphere);
    EXPECT_EQ(r.mMethod, ThicknessMethod::Ray);
    EXPECT_GT(r.mP95, r.mP2);
    EXPECT_LE(r.mP95, r.mMax);
    EXPECT_GE(r.mP2, r.mMin);
    ThicknessOptions o;
    o.mMaxNormalDeviationDeg = 0.0;
    EXPECT_THROW(thickness(sphere, o), std::invalid_argument);
    o = {};
    o.mMethod = ThicknessMethod::Cone;
    o.mConeAngleDeg = 95.0;
    EXPECT_THROW(thickness(sphere, o), std::invalid_argument);
    EXPECT_EQ(thickness_method_from_name("sphere"), ThicknessMethod::Sphere);
    EXPECT_THROW(thickness_method_from_name("laser"), std::invalid_argument);
}
