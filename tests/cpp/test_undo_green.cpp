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
// Tests for green-element undo (`operations/undo_green.hpp`). The core oracle
// (`ExactlyRestoresTheCoarseParentPerGreenGroup`) is a byte-for-byte round
// trip: every green sibling group's substituted row must equal, verbatim, the
// coarse mesh's own row at that group's parent -- a lookup, not an
// approximation, which is what the whole design turns on.

// System includes
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/undo_green.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::refine;
using meshioplusplus::RefineOptions;
using meshioplusplus::RefineResult;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;
using meshioplusplus::undo_green;
using meshioplusplus::UndoGreenResult;

namespace {

// Two triangles per cell of an n x n grid over the unit square -- the same
// fixture test_refine.cpp's own selective suite uses, duplicated here rather
// than shared (that helper is file-private there too).
Mesh tri_grid(std::size_t n) {
    std::vector<std::vector<double>> pts;
    for (std::size_t j = 0; j <= n; ++j)
        for (std::size_t i = 0; i <= n; ++i)
            pts.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});
    std::vector<std::vector<std::int64_t>> cells;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const std::int64_t a = static_cast<std::int64_t>(j * (n + 1) + i);
            const std::int64_t b = a + 1;
            const std::int64_t cc = a + static_cast<std::int64_t>(n) + 2;
            const std::int64_t dd = a + static_cast<std::int64_t>(n) + 1;
            cells.push_back({a, b, cc});
            cells.push_back({a, cc, dd});
        }
    }
    return mt::make_mesh(std::move(pts), "triangle", std::move(cells));
}

RefineOptions cells_opt(std::vector<std::int64_t> cells) {
    RefineOptions o;
    o.mCells = std::move(cells);
    o.mRecordHierarchy = true;
    o.mRecordLevels = true;
    return o;
}

NDArray i64(const std::vector<std::int64_t>& rVals) {
    NDArray a = NDArray::Uninit(meshioplusplus::DType::Int64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

NDArray i64_pairs(const std::vector<std::int64_t>& rFlat) {
    NDArray a = NDArray::Uninit(meshioplusplus::DType::Int64, {rFlat.size() / 2, 2});
    for (std::size_t i = 0; i < rFlat.size(); ++i)
        a.As<std::int64_t>()[i] = rFlat[i];
    return a;
}

bool row_equal(const NDArray& rA, std::size_t RowA, const NDArray& rB, std::size_t RowB,
               std::size_t NodesPerCell) {
    for (std::size_t k = 0; k < NodesPerCell; ++k)
        if (meshioplusplus::detail::read_int(rA, RowA * NodesPerCell + k) !=
            meshioplusplus::detail::read_int(rB, RowB * NodesPerCell + k))
            return false;
    return true;
}

}  // namespace

// --------------------------------------------------------------------------
// The core oracle: byte-for-byte lookup, not reconstruction
// --------------------------------------------------------------------------

TEST(UndoGreen, ExactlyRestoresTheCoarseParentPerGreenGroupAndKeepsRedUnchanged) {
    const Mesh coarse = tri_grid(3);
    const RefineResult fine = refine(coarse, cells_opt({8}));

    ASSERT_TRUE(fine.mMesh.HasCellData("refine:cell_id"));
    ASSERT_TRUE(fine.mMesh.HasCellData("refine:parent_id"));
    ASSERT_TRUE(fine.mMesh.HasCellData("refine:level"));

    const UndoGreenResult undone = undo_green(coarse, fine.mMesh);
    EXPECT_GT(undone.mNumGroupsUndone, 0) << "the fixture must produce at least one green group";
    EXPECT_GT(undone.mNumCellsRemoved, 0);

    // Fewer cells than the fine mesh (green groups collapsed), more than the
    // coarse mesh (the red region is still genuinely refined).
    EXPECT_LT(undone.mMesh.Cells(0).NumCells(), fine.mMesh.Cells(0).NumCells());
    EXPECT_GT(undone.mMesh.Cells(0).NumCells(), coarse.Cells(0).NumCells());

    // Points are never pruned or renumbered.
    EXPECT_EQ(undone.mMesh.NumPoints(), fine.mMesh.NumPoints());

    // The reserved refine:* arrays describe a now-stale hierarchy and are
    // dropped entirely -- four are cell_data, two (entity/hanging) are
    // point_data, so both accessors must be checked.
    for (const char* name :
         {"refine:cell_id", "refine:parent_id", "refine:level", "refine:parent_cell"})
        EXPECT_FALSE(undone.mMesh.HasCellData(name)) << name;
    for (const char* name : {"refine:entity", "refine:hanging"})
        EXPECT_FALSE(undone.mMesh.HasPointData(name)) << name;

    // Classify every fine cell exactly as RedGreenAndUntouchedAreDistinguish...
    // does, then check the undone mesh against the coarse mesh directly for
    // every green group and against the fine mesh for every red child.
    const NDArray& ids = fine.mMesh.CellData("refine:cell_id", 0);
    const NDArray& parents = fine.mMesh.CellData("refine:parent_id", 0);
    const NDArray& levels = fine.mMesh.CellData("refine:level", 0);
    const NDArray& fine_conn = fine.mMesh.Cells(0).Conn();
    const NDArray& coarse_conn = coarse.Cells(0).Conn();
    const NDArray& undone_conn = undone.mMesh.Cells(0).Conn();
    const std::int64_t* cmap = undone.mCellMaps[0].As<std::int64_t>();
    constexpr std::size_t kNpc = 3;  // triangle

    std::size_t red_checked = 0, green_groups_checked = 0;
    std::map<std::int64_t, std::int64_t> green_group_output;  // parent_id -> output row
    for (std::size_t c = 0; c < ids.Size(); ++c) {
        const std::int64_t id = meshioplusplus::detail::read_int(ids, c);
        const std::int64_t parent = meshioplusplus::detail::read_int(parents, c);
        const std::int64_t level = meshioplusplus::detail::read_int(levels, c);
        const std::int64_t out_row = cmap[c];
        if (id == parent) {
            // Untouched: byte-identical to its own fine row.
            EXPECT_TRUE(
                row_equal(fine_conn, c, undone_conn, static_cast<std::size_t>(out_row), kNpc));
            continue;
        }
        if (level == 1) {
            // Red: kept exactly as refine() produced it.
            EXPECT_TRUE(
                row_equal(fine_conn, c, undone_conn, static_cast<std::size_t>(out_row), kNpc));
            ++red_checked;
            continue;
        }
        // Green: the whole group must collapse to ONE output row, and that
        // row must equal the COARSE mesh's own row at `parent` verbatim --
        // parent IS the coarse row index here (tri_grid has no prior
        // refine:cell_id, so ids are implicit == row index).
        auto it = green_group_output.find(parent);
        if (it == green_group_output.end()) {
            green_group_output[parent] = out_row;
            ++green_groups_checked;
        } else {
            EXPECT_EQ(it->second, out_row) << "every member of one green group must map to the "
                                              "same output row";
        }
        EXPECT_TRUE(row_equal(coarse_conn, static_cast<std::size_t>(parent), undone_conn,
                              static_cast<std::size_t>(out_row), kNpc));
    }
    EXPECT_GT(red_checked, 0u);
    EXPECT_EQ(green_groups_checked, static_cast<std::size_t>(undone.mNumGroupsUndone));
}

// --------------------------------------------------------------------------
// Regions: the first genuinely non-injective CellMapKind::Direct use
// --------------------------------------------------------------------------

TEST(UndoGreen, RegionsSurviveTheNonInjectiveCollapseAndSideRegionsDoNotSurvive) {
    const Mesh coarse = tri_grid(3);
    RefineResult fine = refine(coarse, cells_opt({8}));

    // Discover a genuine green pair to build a Cell region spanning two
    // members of the same group (the exact non-injective shape).
    const NDArray& ids = fine.mMesh.CellData("refine:cell_id", 0);
    const NDArray& parents = fine.mMesh.CellData("refine:parent_id", 0);
    const NDArray& levels = fine.mMesh.CellData("refine:level", 0);
    std::vector<std::int64_t> one_green_pair;
    std::map<std::int64_t, std::vector<std::int64_t>> by_parent_green;
    for (std::size_t c = 0; c < ids.Size(); ++c) {
        if (meshioplusplus::detail::read_int(ids, c) ==
            meshioplusplus::detail::read_int(parents, c))
            continue;
        if (meshioplusplus::detail::read_int(levels, c) != 0)
            continue;  // red
        by_parent_green[meshioplusplus::detail::read_int(parents, c)].push_back(
            static_cast<std::int64_t>(c));
    }
    for (auto& kv : by_parent_green)
        if (kv.second.size() >= 2) {
            one_green_pair = {kv.second[0], kv.second[1]};
            break;
        }
    ASSERT_EQ(one_green_pair.size(), 2u) << "the fixture must produce a green group of size >= 2";

    fine.mMesh.AddRegion(Region("greens", RegionKind::Cell, i64(one_green_pair)));
    fine.mMesh.AddRegion(Region("apex", RegionKind::Point, i64({0})));
    fine.mMesh.AddRegion(Region("edge", RegionKind::Side, i64_pairs({0, 0})));

    const UndoGreenResult undone = undo_green(coarse, fine.mMesh);

    const std::size_t greens_idx = undone.mMesh.FindRegion("greens", RegionKind::Cell);
    ASSERT_NE(greens_idx, Mesh::npos);
    const Region& greens = undone.mMesh.Region(greens_idx);
    EXPECT_EQ(greens.NumEntries(), 1u) << "both members of the pair collapse via Canonicalize";

    ASSERT_NE(undone.mMesh.FindRegion("apex", RegionKind::Point), Mesh::npos);
    const Region& apex = undone.mMesh.Region(undone.mMesh.FindRegion("apex", RegionKind::Point));
    EXPECT_EQ(meshioplusplus::detail::read_int(apex.mEntries, 0), 0)
        << "points are never renumbered";

    EXPECT_EQ(undone.mMesh.FindRegion("edge", RegionKind::Side), Mesh::npos)
        << "named Side regions do not survive undo_green at all (mDropSideRegions, since facet "
           "identity is not preserved across a substitution)";
}

// --------------------------------------------------------------------------
// Preconditions and errors
// --------------------------------------------------------------------------

TEST(UndoGreen, MissingHierarchyThrowsNamingTheFix) {
    const Mesh coarse = tri_grid(3);
    RefineOptions o;
    o.mCells = {8};  // no record_hierarchy, no record_levels
    const RefineResult fine = refine(coarse, o);
    EXPECT_THROW(undo_green(coarse, fine.mMesh), std::invalid_argument);
}

TEST(UndoGreen, HierarchyWithoutLevelsStillThrows) {
    const Mesh coarse = tri_grid(3);
    RefineOptions o;
    o.mCells = {8};
    o.mRecordHierarchy = true;  // levels deliberately not requested
    const RefineResult fine = refine(coarse, o);
    EXPECT_THROW(undo_green(coarse, fine.mMesh), std::invalid_argument);
}

TEST(UndoGreen, AnUnrelatedCoarseMeshThrowsRatherThanGuessing) {
    const Mesh coarse = tri_grid(3);
    const RefineResult fine = refine(coarse, cells_opt({8}));
    const Mesh unrelated = tri_grid(1);  // far too few cells to resolve parent ids against
    EXPECT_THROW(undo_green(unrelated, fine.mMesh), std::invalid_argument);
}

TEST(UndoGreen, ACoarseMeshWithMorePointsThanFineThrows) {
    const Mesh coarse = tri_grid(5);
    const RefineResult fine = refine(tri_grid(3), cells_opt({8}));
    EXPECT_THROW(undo_green(coarse, fine.mMesh), std::invalid_argument);
}

TEST(UndoGreen, MultiLevelHierarchyIsRefusedByName) {
    const Mesh coarse = tri_grid(3);
    RefineOptions o = cells_opt({8});
    o.mLevels = 2;
    const RefineResult fine = refine(coarse, o);
    EXPECT_THROW(undo_green(coarse, fine.mMesh), std::invalid_argument);
}

// --------------------------------------------------------------------------
// A no-op case: nothing selected, nothing green
// --------------------------------------------------------------------------

TEST(UndoGreen, UniformRefinementHasNoGreenGroupsAndUndoIsANoOp) {
    // No selector at all -> every cell is refined (uniform mode), so every
    // sibling group is red -- there is nothing to undo.
    const Mesh coarse = tri_grid(2);
    RefineOptions o;
    o.mRecordHierarchy = true;
    o.mRecordLevels = true;
    const RefineResult fine = refine(coarse, o);
    const UndoGreenResult undone = undo_green(coarse, fine.mMesh);
    EXPECT_EQ(undone.mNumGroupsUndone, 0);
    EXPECT_EQ(undone.mNumCellsRemoved, 0);
    EXPECT_EQ(undone.mMesh.Cells(0).NumCells(), fine.mMesh.Cells(0).NumCells());
}
