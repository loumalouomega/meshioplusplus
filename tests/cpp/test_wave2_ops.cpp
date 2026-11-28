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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/subdivide_surface.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/undo_green.hpp"
#include "meshioplusplus/providers.hpp"

using namespace meshioplusplus;

namespace {

TEST(SubdivideSurface, LoopAndCatmullClark) {
    // 1. Loop on a closed triangulated box
    const Mesh vol = box({1, 1, 1}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");
    Mesh skin = extract_surface(vol);
    const std::size_t n_tris = skin.Cells(0).NumCells();
    NDArray f = NDArray::Uninit(DType::Float64, {skin.NumPoints()});
    for (std::size_t i = 0; i < skin.NumPoints(); ++i)
        f.As<double>()[i] = skin.Points().As<double>()[3 * i];
    skin.AddPointData("x", std::move(f));

    SubdivideSurfaceOptions lopts;
    lopts.mScheme = SubdivisionScheme::Loop;
    lopts.mLevels = 2;
    lopts.mCreaseAngle = 0.0;  // smooth everywhere
    lopts.mRecordParentIds = true;
    auto lres = subdivide_surface(skin, lopts);
    EXPECT_EQ(lres.mMesh.Cells(0).Type(), "triangle");
    EXPECT_EQ(lres.mMesh.Cells(0).NumCells(), n_tris * 16u);
    EXPECT_TRUE(lres.mMesh.HasCellData(kSubdivideParentCellName));
    EXPECT_TRUE(lres.mMesh.HasPointData("x"));

    // 2. Catmull-Clark on a quad cube surface (6 quads -> 24 quads after 1 level)
    const Mesh hbox = box({1, 1, 1}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "hexahedron");
    const Mesh qskin = extract_surface(hbox);
    ASSERT_EQ(qskin.Cells(0).Type(), "quad");
    ASSERT_EQ(qskin.Cells(0).NumCells(), 6u);

    SubdivideSurfaceOptions copts;
    copts.mScheme = SubdivisionScheme::CatmullClark;
    copts.mLevels = 1;
    copts.mCreaseAngle = 60.0;  // cube edges are 90 deg -> all creases
    auto cres = subdivide_surface(qskin, copts);
    EXPECT_EQ(cres.mMesh.Cells(0).Type(), "quad");
    EXPECT_EQ(cres.mMesh.Cells(0).NumCells(), 24u);
    EXPECT_EQ(cres.mNumCreaseEdges, 12);

    // Loop refuses quads
    lopts.mScheme = SubdivisionScheme::Loop;
    EXPECT_THROW(subdivide_surface(qskin, lopts), std::invalid_argument);
}

TEST(ConvertCellsOrder, ElevateBeyondOrder2AndReduceRoundTrip) {
    // 1. Tetra -> tetra20 (order 3) -> reduce to order 1
    Mesh tet = box({1, 1, 1}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");
    const std::size_t n_pts_orig = tet.NumPoints();
    NDArray lin = NDArray::Uninit(DType::Float64, {n_pts_orig});
    for (std::size_t i = 0; i < n_pts_orig; ++i) {
        const double* p = tet.Points().As<double>() + 3 * i;
        lin.As<double>()[i] = 2.0 * p[0] - 3.0 * p[1] + 5.0 * p[2] + 1.0;
    }
    tet.AddPointData("lin", std::move(lin));

    ConvertCellsOptions e3;
    e3.mMode = ConvertCellsMode::Elevate;
    e3.mOrder = 3;
    auto r3 = convert_cells(tet, e3);
    ASSERT_EQ(r3.mMesh.Cells(0).Type(), "tetra20");
    ASSERT_EQ(r3.mMesh.Cells(0).NodesPerCell(), 20u);
    // Linear field must be reproduced exactly at all order-3 nodes.
    const auto* pts3 = r3.mMesh.Points().As<double>();
    const auto* val3 = r3.mMesh.PointData("lin").As<double>();
    for (std::size_t i = 0; i < r3.mMesh.NumPoints(); ++i) {
        const double expected =
            2.0 * pts3[3 * i] - 3.0 * pts3[3 * i + 1] + 5.0 * pts3[3 * i + 2] + 1.0;
        EXPECT_NEAR(val3[i], expected, 1e-12);
    }

    // Reduce order 3 -> order 2 (tetra10) -> order 1 (tetra)
    ConvertCellsOptions red2;
    red2.mMode = ConvertCellsMode::Reduce;
    red2.mOrder = 2;
    auto r2 = convert_cells(r3.mMesh, red2);
    EXPECT_EQ(r2.mMesh.Cells(0).Type(), "tetra10");

    ConvertCellsOptions red1;
    red1.mMode = ConvertCellsMode::Reduce;
    red1.mOrder = 1;
    auto r1 = convert_cells(r2.mMesh, red1);
    EXPECT_EQ(r1.mMesh.Cells(0).Type(), "tetra");
    EXPECT_EQ(r1.mMesh.NumPoints(), n_pts_orig);

    // 2. Hexahedron order=2, full=true -> hexahedron27
    const Mesh hex = box({1, 1, 1}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "hexahedron");
    ConvertCellsOptions ef;
    ef.mMode = ConvertCellsMode::Elevate;
    ef.mOrder = 2;
    ef.mFull = true;
    auto rf = convert_cells(hex, ef);
    EXPECT_EQ(rf.mMesh.Cells(0).Type(), "hexahedron27");
    EXPECT_EQ(rf.mMesh.NumPoints(), 27u);
}

TEST(RefineUndoGreen, SingleMeshUndoAndRepeatedSelectivePasses) {
    const Mesh base = box({2, 2, 2}, {2.0, 2.0, 2.0}, {0.0, 0.0, 0.0}, "tetra");
    RefineOptions ropts;
    ropts.mLevels = 1;
    ropts.mCells = {0};
    ropts.mRecordLevels = true;
    ropts.mRecordHierarchy = true;
    auto pass1 = refine(base, ropts);

    // Single-mesh undo_green without passing `base`
    auto ug = undo_green(pass1.mMesh);
    EXPECT_GT(ug.mNumGroupsUndone, 0);
    EXPECT_LT(ug.mMesh.Cells(0).NumCells(), pass1.mMesh.Cells(0).NumCells());

    // 5 repeated selective passes with mUndoGreenFirst=true keep quality bounded
    Mesh cur = pass1.mMesh;
    for (int iter = 0; iter < 4; ++iter) {
        RefineOptions ro;
        ro.mLevels = 1;
        ro.mCells = {0};
        ro.mRecordLevels = true;
        ro.mRecordHierarchy = true;
        ro.mUndoGreenFirst = true;
        cur = refine(cur, ro).mMesh;
    }
    const Mesh qm = attach_quality(cur);
    const auto& sj = qm.CellData("quality:scaled_jacobian", 0);
    double min_sj = 1.0;
    for (std::size_t i = 0; i < sj.Size(); ++i)
        min_sj = std::min(min_sj, sj.As<double>()[i]);
    EXPECT_GT(min_sj, 0.02);
}

TEST(PartitionProviders, CustomProviderAndConnectedParts) {
    const Mesh grid = box({4, 4, 1}, {4.0, 4.0, 1.0}, {0.0, 0.0, 0.0}, "hexahedron");

    // 1. Register a mock GraphPartitionerProvider
    bool called = false;
    GraphPartitionerProvider mock;
    mock.mName = "mock_checker";
    mock.mVersion = "1.0";
    mock.mPartition = [&](const GraphPartitionerRequest& rReq) {
        called = true;
        std::vector<std::int64_t> out(rReq.mNumVertices);
        for (std::size_t i = 0; i < rReq.mNumVertices; ++i)
            out[i] = static_cast<std::int64_t>(i % static_cast<std::size_t>(rReq.mNParts));
        return out;
    };
    add_graph_partitioner(mock);

    PartitionOptions popts;
    popts.mNParts = 2;
    popts.mMethodName = "mock_checker";
    auto lbls = partition_labels(grid, popts);
    EXPECT_TRUE(called);
    ASSERT_EQ(lbls.size(), 1u);

    // 2. With mConnectedParts = true, the checkerboard (every cell disconnected) merges cleanly
    popts.mConnectedParts = true;
    auto clbls = partition_labels(grid, popts);
    ASSERT_EQ(clbls.size(), 1u);
    EXPECT_EQ(clbls[0].Size(), grid.Cells(0).NumCells());

    remove_graph_partitioner("mock_checker");
    EXPECT_EQ(graph_partitioner("mock_checker"), nullptr);
}

}  // namespace
