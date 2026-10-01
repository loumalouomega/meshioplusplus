// Interfaces and contact operations (roadmap §5.2).

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/interfaces.hpp"
#include "meshioplusplus/region.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

Region cell_region(const std::string& rName, std::vector<std::int64_t> rEntries, int Dim = 3,
                   std::int64_t Tag = -1) {
    Region region;
    region.mName = rName;
    region.mKind = RegionKind::Cell;
    region.mDim = Dim;
    region.mTag = Tag;
    region.mEntries = NDArray::Uninit(DType::Int64, {rEntries.size()});
    std::copy(rEntries.begin(), rEntries.end(), region.mEntries.As<std::int64_t>());
    return region;
}

Region side_region(const std::string& rName,
                   std::vector<std::pair<std::int64_t, std::int64_t>> rEntries) {
    Region region;
    region.mName = rName;
    region.mKind = RegionKind::Side;
    region.mEntries = NDArray::Uninit(DType::Int64, {rEntries.size(), std::size_t{2}});
    std::int64_t* data = region.mEntries.As<std::int64_t>();
    for (std::size_t i = 0; i < rEntries.size(); ++i) {
        data[2 * i] = rEntries[i].first;
        data[2 * i + 1] = rEntries[i].second;
    }
    return region;
}

Mesh two_tets() {
    Mesh mesh = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}}, "tetra",
                              {{0, 1, 2, 3}, {0, 2, 1, 4}});
    mesh.AddRegion(cell_region("upper", {0}));
    mesh.AddRegion(cell_region("lower", {1}));
    return mesh;
}

}  // namespace

TEST(RegionAdjacency, FindsOneConformingSharedTriangleAndItsArea) {
    const Mesh input = two_tets();
    const Mesh out = region_adjacency(input);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "triangle");
    ASSERT_EQ(out.Cells(0).NumCells(), 1u);
    EXPECT_EQ(out.CellData("interface:region_a", 0).As<std::int64_t>()[0], 0);
    EXPECT_EQ(out.CellData("interface:region_b", 0).As<std::int64_t>()[0], 1);
    EXPECT_EQ(out.CellData("interface:parent_cell_a", 0).As<std::int64_t>()[0], 1);
    EXPECT_EQ(out.CellData("interface:parent_cell_b", 0).As<std::int64_t>()[0], 0);
    EXPECT_EQ(out.CellData("interface:shared_count", 0).As<std::int64_t>()[0], 2);
    EXPECT_NEAR(out.CellData("interface:measure", 0).As<double>()[0], 0.5, 1e-12);
    ASSERT_EQ(out.NumRegions(), 1u);
    EXPECT_EQ(out.Region(0).mName, "adjacency:0:1");
    EXPECT_EQ(input.NumCellBlocks(), 1u);
}

TEST(RegionAdjacency, SupportsTwoDimensionalEdgesAndExplicitSelectorOrder) {
    Mesh mesh = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}},
                              "quad", {{0, 1, 4, 3}, {1, 2, 5, 4}});
    mesh.AddRegion(cell_region("left", {0}, 2));
    mesh.AddRegion(cell_region("right", {1}, 2));
    RegionSelector left;
    left.mName = "left";
    RegionSelector right;
    right.mName = "right";
    const Mesh out = region_adjacency(mesh, {right, left});
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "line");
    EXPECT_NEAR(out.CellData("interface:measure", 0).As<double>()[0], 1.0, 1e-12);
    EXPECT_EQ(out.CellData("interface:region_a", 0).As<std::int64_t>()[0], 0);
    EXPECT_EQ(out.CellData("interface:region_b", 0).As<std::int64_t>()[0], 1);
}

TEST(RegionAdjacency, UsesCellBlocksAsGroupsWhenNoCellRegionsExist) {
    Mesh mesh = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}}, "tetra",
                              {{0, 1, 2, 3}});
    NDArray second_block(DType::Int64, {1, std::size_t{4}});
    auto second = second_block.As<std::int64_t>();
    second[0] = 0;
    second[1] = 2;
    second[2] = 1;
    second[3] = 4;
    mesh.AddCellBlock("tetra", std::move(second_block));
    const Mesh out = region_adjacency(mesh);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "triangle");
    ASSERT_EQ(out.Cells(0).NumCells(), 1u);
    EXPECT_EQ(out.Region(0).mName, "adjacency:0:1");
}

TEST(RegionAdjacency, RejectsTooFewOrRepeatedRegionsAndResolvesAmbiguity) {
    Mesh mesh = two_tets();
    EXPECT_THROW(region_adjacency(mesh, {{"upper"}}), std::invalid_argument);
    EXPECT_THROW(region_adjacency(mesh, {{"upper"}, {"upper"}}), std::invalid_argument);
    mesh.AddRegion(cell_region("upper", {0}, 3, 8));
    EXPECT_THROW(region_adjacency(mesh, {{"upper"}, {"lower"}}), std::invalid_argument);
    RegionSelector upper;
    upper.mName = "upper";
    upper.mKind = static_cast<std::int32_t>(RegionKind::Cell);
    upper.mDim = 3;
    upper.mTag = 8;
    RegionSelector lower;
    lower.mName = "lower";
    lower.mKind = static_cast<std::int32_t>(RegionKind::Cell);
    lower.mDim = 3;
    EXPECT_NO_THROW(region_adjacency(mesh, {upper, lower}));
}

TEST(FindInterface, ConformingResultHasBothSideRegionsAndReport) {
    Mesh mesh = two_tets();
    RegionSelector upper;
    upper.mName = "upper";
    RegionSelector lower;
    lower.mName = "lower";
    const FindInterfaceResult result = find_interface(mesh, upper, lower);
    ASSERT_EQ(result.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(result.mMesh.Cells(0).Type(), "triangle");
    ASSERT_EQ(result.mReport.mNumPairs, 1);
    EXPECT_NEAR(result.mReport.mArea, 0.5, 1e-12);
    EXPECT_DOUBLE_EQ(result.mReport.mMaxGap, 0.0);
    EXPECT_EQ(result.mSideA.mKind, RegionKind::Side);
    EXPECT_EQ(result.mSideB.mKind, RegionKind::Side);
    EXPECT_EQ(result.mSideA.mEntries.Shape(), (std::vector<std::size_t>{1, 2}));
    EXPECT_EQ(result.mMesh.CellData("interface:parent_facet", 0).As<std::int64_t>()[0], 3);
}

TEST(FindInterface, ProximityMatchesSeparatelyNumberedBoundaryFaces) {
    Mesh a = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra", {{0, 1, 2, 3}});
    Mesh b = mt::make_mesh({{0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, -1}}, "tetra", {{0, 1, 2, 3}});
    a.AddRegion(cell_region("part", {0}));
    b.AddRegion(cell_region("part", {0}));
    RegionSelector part;
    part.mName = "part";
    FindInterfaceOptions options;
    options.mMode = InterfaceMode::Proximity;
    options.mGapTolerance = 1e-10;
    options.mAngleTolerance = 5.0;
    const auto result = find_interface(a, part, b, part, options);
    EXPECT_EQ(result.mReport.mNumPairs, 1);
    EXPECT_NEAR(result.mReport.mArea, 0.5, 1e-12);
    EXPECT_NEAR(result.mReport.mMaxGap, 0.0, 1e-12);
}

TEST(FindInterface, SignedGapFollowsTheChosenMasterSide) {
    Mesh a = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra", {{0, 1, 2, 3}});
    Mesh b = mt::make_mesh({{0, 0, -0.02}, {0, 1, -0.02}, {1, 0, -0.02}, {0, 0, -1}}, "tetra",
                           {{0, 1, 2, 3}});
    a.AddRegion(cell_region("part", {0}));
    b.AddRegion(cell_region("part", {0}));
    RegionSelector part;
    part.mName = "part";
    FindInterfaceOptions options;
    options.mMode = InterfaceMode::Proximity;
    options.mGapTolerance = 0.03;
    options.mAngleTolerance = 5.0;
    options.mMaster = InterfaceMaster::A;
    const auto result_a = find_interface(a, part, b, part, options);
    options.mMaster = InterfaceMaster::B;
    const auto result_b = find_interface(a, part, b, part, options);
    const double gap_a = result_a.mMesh.CellData("interface:gap", 0).As<double>()[0];
    const double gap_b = result_b.mMesh.CellData("interface:gap", 0).As<double>()[0];
    EXPECT_NEAR(std::fabs(gap_a), 0.02, 1e-12);
    EXPECT_NEAR(std::fabs(gap_b), 0.02, 1e-12);
    // Along each master's own outward normal, a separation is negative.
    EXPECT_LT(gap_a, 0.0);
    EXPECT_LT(gap_b, 0.0);
}

TEST(FindInterface, ProximityMasterBReportsEveryFinerFacet) {
    // A has one coarse triangle on z = 0; B tiles the same triangle with four smaller ones.
    Mesh a = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra", {{0, 1, 2, 3}});
    Mesh b = mt::make_mesh({{0, 0, 0},
                            {1, 0, 0},
                            {0, 1, 0},
                            {0.5, 0, 0},
                            {0.5, 0.5, 0},
                            {0, 0.5, 0},
                            {0.3, 0.3, -1}},
                           "tetra", {{0, 5, 3, 6}, {3, 4, 1, 6}, {5, 2, 4, 6}, {3, 5, 4, 6}});
    a.AddRegion(cell_region("part", {0}));
    b.AddRegion(cell_region("part", {0, 1, 2, 3}));
    RegionSelector part;
    part.mName = "part";
    FindInterfaceOptions options;
    options.mMode = InterfaceMode::Proximity;
    options.mGapTolerance = 1e-10;
    options.mAngleTolerance = 5.0;
    options.mMaster = InterfaceMaster::A;
    const auto result_a = find_interface(a, part, b, part, options);
    EXPECT_EQ(result_a.mReport.mNumPairs, 1);
    EXPECT_NEAR(result_a.mReport.mArea, 0.5, 1e-12);
    options.mMaster = InterfaceMaster::B;
    const auto result_b = find_interface(a, part, b, part, options);
    EXPECT_EQ(result_b.mReport.mNumPairs, 4);
    EXPECT_NEAR(result_b.mReport.mArea, 0.5, 1e-12);
    EXPECT_EQ(result_b.mReport.mUnmatchedB, 0);
}

TEST(ContactPairs, ProjectsPointRegionToClosestMasterFacet) {
    Mesh mesh = two_tets();
    NDArray point_ids = NDArray::Uninit(DType::Int64, {3});
    point_ids.As<std::int64_t>()[0] = 0;
    point_ids.As<std::int64_t>()[1] = 1;
    point_ids.As<std::int64_t>()[2] = 2;
    Region points;
    points.mName = "slave";
    points.mKind = RegionKind::Point;
    points.mEntries = std::move(point_ids);
    mesh.AddRegion(points);
    RegionSelector slave;
    slave.mName = "slave";
    RegionSelector master;
    master.mName = "lower";
    const ContactPairsResult result = contact_pairs(mesh, slave, mesh, master);
    ASSERT_EQ(result.mSlavePoint.Size(), 3u);
    EXPECT_EQ(result.mUnmatched.Size(), 0u);
    for (std::size_t i = 0; i < result.mMasterCell.Size(); ++i)
        EXPECT_EQ(result.mMasterCell.As<std::int64_t>()[i], 1);
    EXPECT_NEAR(result.mGap.As<double>()[0], 0.0, 1e-12);
    EXPECT_EQ(result.mLocalCoordinates.Shape(), (std::vector<std::size_t>{3, 3}));
}

TEST(ContactPairs, ReportsPolygonFanTriangleForLocalCoordinates) {
    Mesh mesh = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.25, 0.75, 0.1}},
                              "quad", {{0, 1, 2, 3}});
    NDArray point_ids = NDArray::Uninit(DType::Int64, {1});
    point_ids.As<std::int64_t>()[0] = 4;
    Region points;
    points.mName = "slave";
    points.mKind = RegionKind::Point;
    points.mEntries = std::move(point_ids);
    mesh.AddRegion(points);
    mesh.AddRegion(cell_region("master", {0}));
    RegionSelector slave;
    slave.mName = "slave";
    RegionSelector master;
    master.mName = "master";
    ContactPairsOptions options;
    options.mTolerance = 0.2;
    const auto result = contact_pairs(mesh, slave, mesh, master, options);
    EXPECT_EQ(result.mMasterSubfacet.As<std::int64_t>()[0], 1);
    EXPECT_NEAR(result.mLocalCoordinates.As<double>()[0], 0.25, 1e-12);
    EXPECT_NEAR(result.mLocalCoordinates.As<double>()[1], 0.25, 1e-12);
    EXPECT_NEAR(result.mLocalCoordinates.As<double>()[2], 0.5, 1e-12);
}

TEST(SplitInterface, DuplicatesOnlyTheOppositePointFanAndAddsWedge) {
    Mesh mesh = two_tets();
    NDArray point_data(DType::Float64, {5});
    for (std::size_t i = 0; i < 5; ++i)
        point_data.As<double>()[i] = static_cast<double>(i);
    mesh.AddPointData("temperature", std::move(point_data));
    mesh.AddRegion(side_region("cut", {{0, 3}}));
    RegionSelector cut;
    cut.mName = "cut";
    SplitInterfaceOptions options;
    options.mAddCohesive = true;
    const SplitInterfaceResult result = split_interface(mesh, cut, options);
    EXPECT_EQ(result.mNumDuplicatedPoints, 3);
    EXPECT_EQ(result.mNumCohesiveCells, 1);
    EXPECT_EQ(result.mMesh.NumPoints(), 8u);
    EXPECT_EQ(result.mMesh.NumCellBlocks(), 2u);
    EXPECT_EQ(result.mMesh.Cells(1).Type(), "wedge");
    EXPECT_EQ(result.mMesh.PointData("temperature").Size(), 8u);
    EXPECT_EQ(result.mMesh.PointData("temperature").As<double>()[7], 2.0);
}

TEST(SplitInterface, EncodesTheSecondTwoDimensionalCohesiveTraceAsCellData) {
    Mesh mesh = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}},
                              "quad", {{0, 1, 4, 3}, {1, 2, 5, 4}});
    mesh.AddRegion(side_region("cut", {{0, 1}}));
    RegionSelector cut;
    cut.mName = "cut";
    SplitInterfaceOptions options;
    options.mAddCohesive = true;
    const auto result = split_interface(mesh, cut, options);
    ASSERT_EQ(result.mNumCohesiveCells, 1);
    ASSERT_EQ(result.mMesh.NumCellBlocks(), 2u);
    EXPECT_EQ(result.mMesh.Cells(1).Type(), "line");
    const NDArray& trace_b = result.mMesh.CellData("cohesive:trace_b", 1);
    EXPECT_EQ(trace_b.Shape(), (std::vector<std::size_t>{1, 2}));
    EXPECT_NE(trace_b.As<std::int64_t>()[0], trace_b.As<std::int64_t>()[1]);
}

TEST(SplitInterface, BothSideEntriesStillCreateOneCohesiveCell) {
    Mesh mesh = two_tets();
    mesh.AddRegion(side_region("cut", {{0, 3}, {1, 3}}));
    RegionSelector cut;
    cut.mName = "cut";
    SplitInterfaceOptions options;
    options.mAddCohesive = true;
    const auto result = split_interface(mesh, cut, options);
    EXPECT_EQ(result.mNumCohesiveCells, 1);
}
