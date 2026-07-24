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
// Tests for the split operation (by type / component / tag).

// System includes
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;
using meshioplusplus::split;
using meshioplusplus::SplitBy;
using meshioplusplus::SplitResult;

Mesh mixed_tri_quad() {
    Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 3}}));
    m.AddCellBlock("quad", mt::conn_from({{1, 2, 5, 4}}));
    return m;
}

TEST(Split, ByType) {
    SplitResult r = split(mixed_tri_quad(), SplitBy::Type);
    ASSERT_EQ(r.mPieces.size(), 2u);
    std::int64_t total = 0;
    for (const auto& p : r.mPieces) {
        ASSERT_EQ(p.mMesh.NumCellBlocks(), 1u);  // each piece is one type
        total += static_cast<std::int64_t>(p.mMesh.Cells(0).NumCells());
    }
    EXPECT_EQ(total, 2);
    EXPECT_EQ(r.mPieces[0].mKey, "triangle");
    EXPECT_EQ(r.mPieces[1].mKey, "quad");
}

TEST(Split, ByComponentTwoBlobs) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {5, 0, 0}, {6, 0, 0}, {5, 1, 0}},
                           "triangle", {{0, 1, 2}, {3, 4, 5}});
    SplitResult r = split(m, SplitBy::Component);
    ASSERT_EQ(r.mPieces.size(), 2u);
    for (const auto& p : r.mPieces) {
        EXPECT_EQ(p.mMesh.NumPoints(), 3u);
        EXPECT_EQ(p.mMesh.Cells(0).NumCells(), 1u);
    }
}

TEST(Split, ByTag) {
    Mesh m = mixed_tri_quad();
    m.AddCellData("region", {mt::conn_from({{7}}), mt::conn_from({{9}})});
    SplitResult r = split(m, SplitBy::Tag, "region");
    ASSERT_EQ(r.mPieces.size(), 2u);
    EXPECT_EQ(r.mPieces[0].mKey, "7");
    EXPECT_EQ(r.mPieces[1].mKey, "9");
}

TEST(Split, ByRegionOnePiecePerNamedCellRegion) {
    // triangle 3-1-4 is the last triangle in the block below; add a second
    // triangle so the block has two cells, one per named region.
    Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 4}, {1, 2, 5}}));

    NDArray first = NDArray::Uninit(DType::Int64, {1});
    first.As<std::int64_t>()[0] = 0;
    m.AddRegion(Region("left", RegionKind::Cell, 2, 1, std::move(first)));
    NDArray second = NDArray::Uninit(DType::Int64, {1});
    second.As<std::int64_t>()[0] = 1;
    m.AddRegion(Region("right", RegionKind::Cell, 2, 2, std::move(second)));

    SplitResult r = split(m, SplitBy::Region);
    ASSERT_EQ(r.mPieces.size(), 2u);
    // RegionNames() is sorted, so "left" < "right".
    EXPECT_EQ(r.mPieces[0].mKey, "left");
    EXPECT_EQ(r.mPieces[0].mMesh.Cells(0).NumCells(), 1u);
    EXPECT_EQ(r.mPieces[1].mKey, "right");
    EXPECT_EQ(r.mPieces[1].mMesh.Cells(0).NumCells(), 1u);
}

TEST(Split, ByRegionOverlappingRegionsAreNotAPartition) {
    // Unlike every other criterion, a cell in two regions must appear in both
    // output pieces -- split-by-region is not required to partition the mesh.
    Mesh m = mixed_tri_quad();  // triangle block (1 cell), quad block (1 cell)

    NDArray tri = NDArray::Uninit(DType::Int64, {1});
    tri.As<std::int64_t>()[0] = 0;  // global cell 0 = the triangle
    m.AddRegion(Region("a", RegionKind::Cell, std::move(tri)));

    NDArray both = NDArray::Uninit(DType::Int64, {2});
    both.As<std::int64_t>()[0] = 0;  // the same triangle...
    both.As<std::int64_t>()[1] = 1;  // ...plus the quad
    m.AddRegion(Region("b", RegionKind::Cell, std::move(both)));

    SplitResult r = split(m, SplitBy::Region);
    ASSERT_EQ(r.mPieces.size(), 2u);
    // "a" -> just the triangle; "b" -> both cells. Total cells across pieces
    // (2) exceeds the mesh's own cell count (2 blocks x 1 cell), confirming
    // the triangle was NOT consumed by the first piece it matched.
    std::size_t total = 0;
    for (const auto& p : r.mPieces)
        for (std::size_t i = 0; i < p.mMesh.NumCellBlocks(); ++i)
            total += p.mMesh.Cells(i).NumCells();
    EXPECT_EQ(total, 3u);  // "a": 1 cell, "b": 2 cells
}

TEST(Split, ByRegionIgnoresNonCellRegions) {
    // Point/Side regions produce no piece -- split's contract is whole
    // submeshes, and there is no sound default for "these facets alone".
    Mesh m = mixed_tri_quad();
    NDArray pts = NDArray::Uninit(DType::Int64, {1});
    pts.As<std::int64_t>()[0] = 0;
    m.AddRegion(Region("fixed", RegionKind::Point, std::move(pts)));

    SplitResult r = split(m, SplitBy::Region);
    EXPECT_TRUE(r.mPieces.empty());
}

TEST(Split, ByFromNameAcceptsRegionsPlural) {
    EXPECT_EQ(meshioplusplus::split_by_from_name("regions"), SplitBy::Region);
    // The pre-existing singular alias for Tag is untouched.
    EXPECT_EQ(meshioplusplus::split_by_from_name("region"), SplitBy::Tag);
}

}  // namespace
