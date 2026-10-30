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
