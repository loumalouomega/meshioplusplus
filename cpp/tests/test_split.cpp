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

namespace {

using meshioplusplus::Mesh;
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

}  // namespace
