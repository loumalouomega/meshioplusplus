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
// Tests for the mesh cleanup operation (weld / prune / de-dup).

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/clean.hpp"

namespace {

using meshioplusplus::clean;
using meshioplusplus::CleanOptions;
using meshioplusplus::CleanResult;
using meshioplusplus::Mesh;

TEST(Clean, WeldReducesPointCount) {
    // node 3 coincides with node 0
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 0}}, "triangle",
                           {{0, 1, 2}, {3, 1, 2}});
    CleanOptions opts;
    opts.weld = true;
    opts.atol = 1e-9;
    opts.remove_orphans = false;
    opts.drop_degenerate = false;
    opts.drop_duplicate_cells = false;
    CleanResult r = clean(m, opts);
    EXPECT_EQ(r.mMesh.NumPoints(), 3u);
    EXPECT_EQ(r.mPointsWelded, 1);
}

TEST(Clean, RemovesOrphans) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {9, 9, 9}}, "triangle", {{0, 1, 2}});
    CleanOptions opts;  // defaults: remove orphans on
    CleanResult r = clean(m, opts);
    EXPECT_EQ(r.mMesh.NumPoints(), 3u);
    EXPECT_EQ(r.mPointsRemovedOrphan, 1);
}

TEST(Clean, DropsDegenerateAndDuplicate) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, "triangle",
                           {{0, 1, 2}, {0, 1, 2}, {0, 0, 1}});
    CleanOptions opts;
    CleanResult r = clean(m, opts);
    EXPECT_EQ(r.mCellsDroppedDuplicate, 1);
    EXPECT_EQ(r.mCellsDroppedDegenerate, 1);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 1u);
}

TEST(Clean, DefaultDoesNotWeld) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 0}}, "triangle",
                           {{0, 1, 2}, {3, 1, 2}});
    CleanOptions opts;  // weld off by default
    CleanResult r = clean(m, opts);
    EXPECT_EQ(r.mPointsWelded, 0);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 2u);
}

}  // namespace
