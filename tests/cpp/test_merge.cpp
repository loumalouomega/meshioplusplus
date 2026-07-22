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
// Unit tests for the merge/weld operation (operations/merge.hpp). Exercised
// under every mesh backend via the cpp-tests CI matrix.

// System includes
#include <cmath>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "mesh_fixtures.hpp"

using meshioplusplus::DType;
using meshioplusplus::MergeDataPolicy;
using meshioplusplus::MergeOptions;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

namespace {

// Two adjacent unit triangles: mesh B is shifted so its left edge coincides
// with mesh A's right edge (points (1,0,0) and (1,1,0) are shared).
Mesh block_a() {
    return mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle",
                         {{0, 1, 2}, {0, 2, 3}});
}
Mesh block_b() {
    return mt::make_mesh({{1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}}, "triangle",
                         {{0, 1, 2}, {0, 2, 3}});
}

std::size_t total_cells(const Mesh& rM) {
    std::size_t n = 0;
    for (const auto cb : rM.CellRange())
        n += cb.NumCells();
    return n;
}

// Maximum node index referenced by any cell (to check connectivity is in range).
std::int64_t max_conn(const Mesh& rM) {
    std::int64_t m = -1;
    for (const auto cb : rM.CellRange()) {
        const NDArray& conn = cb.Conn();
        const std::size_t n = cb.NumCells() * cb.NodesPerCell();
        for (std::size_t i = 0; i < n; ++i)
            m = std::max(m, meshioplusplus::detail::read_int(conn, i));
    }
    return m;
}

}  // namespace

TEST(Merge, ConcatenationCounts) {
    Mesh a = block_a();
    Mesh b = block_b();
    std::vector<const Mesh*> meshes{&a, &b};
    meshioplusplus::MergeResult r = meshioplusplus::merge(meshes);

    EXPECT_EQ(r.mMesh.NumPoints(), a.NumPoints() + b.NumPoints());
    EXPECT_EQ(total_cells(r.mMesh), total_cells(a) + total_cells(b));
    // Same-type blocks merged into one.
    EXPECT_EQ(r.mMesh.NumCellBlocks(), 1u);
    // All connectivity indices stay in [0, total_points).
    EXPECT_LT(max_conn(r.mMesh), static_cast<std::int64_t>(r.mMesh.NumPoints()));
    EXPECT_GE(max_conn(r.mMesh), 0);
}

TEST(Merge, OffsetCorrectness) {
    Mesh a = block_a();
    Mesh b = block_b();
    std::vector<const Mesh*> meshes{&a, &b};
    meshioplusplus::MergeResult r = meshioplusplus::merge(meshes);

    // B's first cell is local {0,1,2} -> offset by 4 -> {4,5,6}. Its first node
    // must land on B's original point 0 = (1,0,0).
    const auto cb = r.mMesh.Cells(0);
    const NDArray& conn = cb.Conn();
    // Block order: A's 2 cells then B's 2 cells.
    const std::int64_t n0 = meshioplusplus::detail::read_int(conn, 2 * cb.NodesPerCell() + 0);
    EXPECT_EQ(n0, 4);
    const NDArray& pts = r.mMesh.Points();
    const std::size_t dim = r.mMesh.PointDim();
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(pts, 4 * dim + 0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(pts, 4 * dim + 1), 0.0);
}

TEST(Merge, SourceTag) {
    Mesh a = block_a();
    Mesh b = block_b();
    std::vector<const Mesh*> meshes{&a, &b};
    meshioplusplus::MergeResult r = meshioplusplus::merge(meshes);

    ASSERT_TRUE(r.mMesh.HasCellData("source_mesh_id"));
    const NDArray& tag = r.mMesh.CellData("source_mesh_id", 0);
    ASSERT_EQ(meshioplusplus::detail::rows(tag), total_cells(r.mMesh));
    // A's two cells tagged 0, B's two cells tagged 1.
    EXPECT_EQ(meshioplusplus::detail::read_int(tag, 0), 0);
    EXPECT_EQ(meshioplusplus::detail::read_int(tag, 1), 0);
    EXPECT_EQ(meshioplusplus::detail::read_int(tag, 2), 1);
    EXPECT_EQ(meshioplusplus::detail::read_int(tag, 3), 1);
}

TEST(Merge, NoSourceTag) {
    Mesh a = block_a();
    Mesh b = block_b();
    std::vector<const Mesh*> meshes{&a, &b};
    MergeOptions opts;
    opts.source_tag = false;
    meshioplusplus::MergeResult r = meshioplusplus::merge(meshes, opts);
    EXPECT_FALSE(r.mMesh.HasCellData("source_mesh_id"));
}

TEST(Merge, WeldSharedInterface) {
    Mesh a = block_a();
    Mesh b = block_b();
    std::vector<const Mesh*> meshes{&a, &b};
    MergeOptions opts;
    opts.weld = true;
    opts.atol = 1e-9;
    meshioplusplus::MergeResult r = meshioplusplus::merge(meshes, opts);

    // 8 input points, 2 coincident pairs -> 6 output points.
    EXPECT_EQ(r.mMesh.NumPoints(), 6u);
    EXPECT_EQ(total_cells(r.mMesh), 4u);
    EXPECT_LT(max_conn(r.mMesh), 6);
}

TEST(Merge, WeldTolerance) {
    Mesh a = block_a();
    // B shifted by 1e-3 in x, so its left edge is *near* but not *at* A's right.
    Mesh b = mt::make_mesh({{1.001, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1.001, 1, 0}}, "triangle",
                           {{0, 1, 2}, {0, 2, 3}});
    std::vector<const Mesh*> meshes{&a, &b};

    // Below tolerance: nothing welds.
    MergeOptions tight;
    tight.weld = true;
    tight.atol = 1e-6;
    EXPECT_EQ(meshioplusplus::merge(meshes, tight).mMesh.NumPoints(), 8u);

    // Above tolerance: the 2 near-coincident pairs weld.
    MergeOptions loose;
    loose.weld = true;
    loose.atol = 1e-2;
    EXPECT_EQ(meshioplusplus::merge(meshes, loose).mMesh.NumPoints(), 6u);
}

TEST(Merge, DropDuplicateCells) {
    // Merge a mesh with itself and weld: with dedup the cells collapse.
    Mesh a = block_a();
    Mesh a2 = block_a();
    std::vector<const Mesh*> meshes{&a, &a2};
    MergeOptions opts;
    opts.weld = true;
    opts.atol = 1e-9;
    opts.drop_duplicate_cells = true;
    meshioplusplus::MergeResult r = meshioplusplus::merge(meshes, opts);
    EXPECT_EQ(r.mMesh.NumPoints(), 4u);
    EXPECT_EQ(total_cells(r.mMesh), 2u);
}

TEST(Merge, DataPolicyIntersection) {
    Mesh a = block_a();
    a.AddPointData("T", mt::points_from({{1}, {2}, {3}, {4}}));
    a.AddPointData("onlyA", mt::points_from({{9}, {8}, {7}, {6}}));
    Mesh b = block_b();
    b.AddPointData("T", mt::points_from({{5}, {6}, {7}, {8}}));
    std::vector<const Mesh*> meshes{&a, &b};

    // Intersection: 'T' survives (present in both), 'onlyA' dropped.
    meshioplusplus::MergeResult ri = meshioplusplus::merge(meshes);
    EXPECT_TRUE(ri.mMesh.HasPointData("T"));
    EXPECT_FALSE(ri.mMesh.HasPointData("onlyA"));
    EXPECT_EQ(meshioplusplus::detail::rows(ri.mMesh.PointData("T")), 8u);

    // Fill: 'onlyA' kept, NaN in B's rows.
    MergeOptions fill;
    fill.data_policy = MergeDataPolicy::Fill;
    meshioplusplus::MergeResult rf = meshioplusplus::merge(meshes, fill);
    ASSERT_TRUE(rf.mMesh.HasPointData("onlyA"));
    const NDArray& oa = rf.mMesh.PointData("onlyA");
    EXPECT_TRUE(std::isnan(meshioplusplus::detail::read_double(oa, 4)));
}

TEST(Merge, PointCellMaps) {
    Mesh a = block_a();
    Mesh b = block_b();
    std::vector<const Mesh*> meshes{&a, &b};
    MergeOptions opts;
    opts.weld = true;
    opts.atol = 1e-9;
    meshioplusplus::MergeResult r = meshioplusplus::merge(meshes, opts);

    ASSERT_EQ(r.mPointMaps.size(), 2u);
    ASSERT_EQ(r.mCellMaps.size(), 2u);
    // Every mapped point index is valid.
    for (const NDArray& pm : r.mPointMaps) {
        for (std::size_t i = 0; i < meshioplusplus::detail::rows(pm); ++i) {
            const std::int64_t v = meshioplusplus::detail::read_int(pm, i);
            EXPECT_GE(v, 0);
            EXPECT_LT(v, static_cast<std::int64_t>(r.mMesh.NumPoints()));
        }
    }
    // Every mapped cell index is valid (no dedup -> none dropped).
    for (const NDArray& cm : r.mCellMaps) {
        for (std::size_t i = 0; i < meshioplusplus::detail::rows(cm); ++i) {
            const std::int64_t v = meshioplusplus::detail::read_int(cm, i);
            EXPECT_GE(v, 0);
            EXPECT_LT(v, static_cast<std::int64_t>(total_cells(r.mMesh)));
        }
    }
}
