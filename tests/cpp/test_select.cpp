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
// Tests for select_region: every criterion and region kind, their
// intersection, the data-only contract, and the hand-off to find_interface.

// System includes
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "analytic_surfaces.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/interfaces.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/select.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/region.hpp"

using namespace meshioplusplus;

namespace {

// A 2 x 2 x 2 hexahedral cube of side 2, with primitives' side regions.
Mesh cube() {
    return box({{2, 2, 2}}, {{2.0, 2.0, 2.0}});
}

std::vector<std::int64_t> entries(const Region& rRegion) {
    const std::int64_t* p = rRegion.Entries();
    return std::vector<std::int64_t>(p, p + rRegion.mEntries.Size());
}

NDArray f64(const std::vector<double>& rV) {
    NDArray a = NDArray::Uninit(DType::Float64, {rV.size()});
    std::copy(rV.begin(), rV.end(), a.As<double>());
    return a;
}

bool same_bytes(const NDArray& rA, const NDArray& rB) {
    return rA.Dtype() == rB.Dtype() && rA.Shape() == rB.Shape() &&
           std::memcmp(rA.Data(), rB.Data(), rA.Nbytes()) == 0;
}

std::int64_t total_cells(const Mesh& rMesh) {
    std::int64_t n = 0;
    for (const auto cb : rMesh.CellRange())
        n += static_cast<std::int64_t>(cb.NumCells());
    return n;
}

SelectOptions named(const char* pName, RegionKind Kind) {
    SelectOptions o;
    o.mName = pName;
    o.mKind = Kind;
    return o;
}

}  // namespace

TEST(Select, BoxSidesMatchThePrimitiveRegion) {
    const Mesh m = cube();
    SelectOptions o = named("left", RegionKind::Side);
    o.mBox = true;
    o.mBoxLo = {-1.0, -1.0, -1.0};
    o.mBoxHi = {1e-9, 3.0, 3.0};
    const SelectResult r = select_region(m, o);
    EXPECT_EQ(r.mNumCandidates, 24);  // 6 faces x 4 boundary quads
    EXPECT_EQ(r.mNumSelected, 4);
    EXPECT_EQ(r.mRegion.mDim, 2);
    const std::size_t xmin = m.FindRegion("xmin", RegionKind::Side);
    ASSERT_NE(xmin, Mesh::npos);
    EXPECT_EQ(entries(r.mRegion), entries(m.Region(xmin)));
    const std::size_t left = r.mMesh.FindRegion("left", RegionKind::Side);
    ASSERT_NE(left, Mesh::npos);
    EXPECT_TRUE(regions_equal(r.mMesh.Region(left), r.mRegion));
}

TEST(Select, InteriorFacesWithoutBoundaryOnly) {
    SelectOptions o = named("mid", RegionKind::Side);
    o.mBoundaryOnly = false;
    o.mPlane = true;
    o.mPlanePoint = {1.0, 0.0, 0.0};
    o.mPlaneNormal = {1.0, 0.0, 0.0};
    o.mBox = true;
    o.mBoxLo = {-1.0, -1.0, -1.0};
    o.mBoxHi = {1.0, 3.0, 3.0};
    const SelectResult r = select_region(cube(), o);
    EXPECT_EQ(r.mNumCandidates, 8 * 6);
    // The plane x = 1 is four faces, each named from both of its cells.
    EXPECT_EQ(r.mNumSelected, 8);
}

TEST(Select, CellModes) {
    const Mesh m = cube();
    SelectOptions o = named("half", RegionKind::Cell);
    o.mPlane = true;
    o.mPlanePoint = {0.6, 0.0, 0.0};
    o.mPlaneNormal = {1.0, 0.0, 0.0};
    o.mMode = SelectMode::All;
    EXPECT_EQ(select_region(m, o).mNumSelected, 4);
    o.mMode = SelectMode::Any;
    EXPECT_EQ(select_region(m, o).mNumSelected, 8);
    o.mMode = SelectMode::Centroid;
    EXPECT_EQ(select_region(m, o).mNumSelected, 4);
    // A centroid on the plane is inside: the half-space is closed.
    o.mPlanePoint = {1.5, 0.0, 0.0};
    EXPECT_EQ(select_region(m, o).mNumSelected, 4);
    o.mPlanePoint = {1.6, 0.0, 0.0};
    EXPECT_EQ(select_region(m, o).mNumSelected, 0);
    o.mMode = SelectMode::Any;
    EXPECT_EQ(select_region(m, o).mNumSelected, 4);
}

TEST(Select, PointsInASphere) {
    SelectOptions o = named("corner", RegionKind::Point);
    o.mSphere = true;
    o.mSphereCenter = {0.0, 0.0, 0.0};
    o.mSphereRadius = 1.0;
    const SelectResult r = select_region(cube(), o);
    EXPECT_EQ(r.mNumCandidates, 27);
    EXPECT_EQ(r.mNumSelected, 4);  // the corner and its three axis neighbours
    EXPECT_EQ(r.mRegion.mDim, -1);
}

TEST(Select, WhereOnCellAndPointData) {
    Mesh m = cube();
    std::vector<double> mat = {1, 2, 1, 2, 1, 2, 1, 2};
    m.AddCellData("mat", {f64(mat)});
    std::vector<double> x(m.NumPoints());
    for (std::size_t i = 0; i < x.size(); ++i)
        x[i] = m.Points().As<double>()[3 * i];
    m.AddPointData("x", f64(x));

    SelectOptions o = named("two", RegionKind::Cell);
    o.mWhereArray = "mat";
    o.mWhereOp = RefineCompare::Equal;
    o.mWhereValue = 2.0;
    EXPECT_EQ(entries(select_region(m, o).mRegion), (std::vector<std::int64_t>{1, 3, 5, 7}));

    // A side takes its owner's value.
    o.mKind = RegionKind::Side;
    EXPECT_EQ(select_region(m, o).mNumSelected, 12);

    o.mKind = RegionKind::Point;
    o.mWhereArray = "x";
    o.mWhereOp = RefineCompare::GreaterEqual;
    o.mWhereValue = 2.0;
    EXPECT_EQ(select_region(m, o).mNumSelected, 9);

    // The wrong location is refused by name.
    o.mWhereArray = "mat";
    EXPECT_THROW(select_region(m, o), std::invalid_argument);
    o.mKind = RegionKind::Cell;
    o.mWhereArray = "x";
    EXPECT_THROW(select_region(m, o), std::invalid_argument);
    o.mWhereArray = "nope";
    EXPECT_THROW(select_region(m, o), std::invalid_argument);
}

TEST(Select, QualityThreshold) {
    SelectOptions o = named("bad", RegionKind::Cell);
    o.mQualityMetric = "scaled_jacobian";
    o.mQualityOp = RefineCompare::Less;
    o.mQualityValue = 0.5;
    EXPECT_EQ(select_region(cube(), o).mNumSelected, 0);
    o.mQualityOp = RefineCompare::Greater;
    EXPECT_EQ(select_region(cube(), o).mNumSelected, 8);
    o.mQualityMetric = "quality:scaled_jacobian";  // attach_quality's spelling
    EXPECT_EQ(select_region(cube(), o).mNumSelected, 8);
    o.mQualityMetric = "no_such_metric";
    EXPECT_THROW(select_region(cube(), o), std::invalid_argument);
    o.mQualityMetric = "scaled_jacobian";
    o.mKind = RegionKind::Point;
    EXPECT_THROW(select_region(cube(), o), std::invalid_argument);
}

TEST(Select, ComponentMatchesSplit) {
    const Mesh a = cube();
    const Mesh b = box({{1, 1, 1}}, {{1.0, 1.0, 1.0}}, {{5.0, 0.0, 0.0}});
    const Mesh m = merge({&a, &b}).mMesh;
    SelectOptions o = named("second", RegionKind::Cell);
    o.mComponent = 1;
    const SelectResult r = select_region(m, o);
    EXPECT_EQ(r.mNumComponents, 2);
    EXPECT_EQ(r.mNumSelected, 1);
    const std::vector<SplitPiece> pieces = split(m, SplitBy::Component).mPieces;
    ASSERT_EQ(pieces.size(), 2u);
    EXPECT_EQ(total_cells(pieces[1].mMesh), r.mNumSelected);

    o.mKind = RegionKind::Point;
    EXPECT_EQ(select_region(m, o).mNumSelected, 8);
    o.mKind = RegionKind::Side;
    EXPECT_EQ(select_region(m, o).mNumSelected, 6);
    o.mComponent = 2;
    EXPECT_EQ(select_region(m, o).mNumSelected, 0);
}

TEST(Select, GrowStopsAtFeatureEdges) {
    const Mesh m = analytic::box_skin(3);
    SelectOptions o = named("top", RegionKind::Cell);
    o.mSeed = true;
    o.mSeedPoint = {0.5, 0.5, 1.2};
    SelectResult r = select_region(m, o);
    EXPECT_EQ(r.mNumGrown, 9);
    EXPECT_EQ(r.mNumSelected, 9);
    // Every selected quad lies in z = 1.
    for (const std::int64_t c : entries(r.mRegion))
        for (int k = 0; k < 4; ++k) {
            const std::int64_t v = m.Cells(0).Conn().As<std::int64_t>()[4 * c + k];
            EXPECT_DOUBLE_EQ(m.Points().As<double>()[3 * v + 2], 1.0);
        }
    // A feature angle above 90 degrees crosses the cube's edges.
    o.mFeatureAngleDeg = 95.0;
    EXPECT_EQ(select_region(m, o).mNumSelected, 54);
    // A normal tolerance holds the patch to the seed's face again.
    o.mNormalToleranceDeg = 10.0;
    EXPECT_EQ(select_region(m, o).mNumSelected, 9);
    // The points of the patch.
    o.mKind = RegionKind::Point;
    EXPECT_EQ(select_region(m, o).mNumSelected, 16);
}

TEST(Select, GrowSidesOfAVolume) {
    const Mesh m = cube();
    SelectOptions o = named("top", RegionKind::Side);
    o.mSeed = true;
    o.mSeedPoint = {1.0, 1.0, 5.0};
    const SelectResult r = select_region(m, o);
    const std::size_t zmax = m.FindRegion("zmax", RegionKind::Side);
    ASSERT_NE(zmax, Mesh::npos);
    EXPECT_EQ(entries(r.mRegion), entries(m.Region(zmax)));
    // On a volume without surface cells, points grow over the boundary faces.
    o.mKind = RegionKind::Point;
    EXPECT_EQ(select_region(m, o).mNumSelected, 9);
    // A surface mesh has edges for sides: growing them is refused.
    o.mKind = RegionKind::Side;
    EXPECT_THROW(select_region(analytic::box_skin(2), o), std::invalid_argument);
    // So is growing cells of a volume mesh.
    o.mKind = RegionKind::Cell;
    EXPECT_THROW(select_region(m, o), std::invalid_argument);
}

TEST(Select, CriteriaIntersect) {
    Mesh m = cube();
    std::vector<double> mat = {1, 2, 1, 2, 1, 2, 1, 2};
    m.AddCellData("mat", {f64(mat)});
    SelectOptions o = named("both", RegionKind::Side);
    o.mSeed = true;
    o.mSeedPoint = {1.0, 1.0, 5.0};
    o.mWhereArray = "mat";
    o.mWhereOp = RefineCompare::Equal;
    o.mWhereValue = 2.0;
    const SelectResult r = select_region(m, o);
    EXPECT_EQ(r.mNumGrown, 4);
    EXPECT_EQ(r.mNumSelected, 2);
}

TEST(Select, IsDataOnlyAndReplacesItsRegion) {
    const Mesh m = cube();
    SelectOptions o = named("xmin", RegionKind::Cell);
    o.mBox = true;
    o.mBoxLo = {-1.0, -1.0, -1.0};
    o.mBoxHi = {1.0, 3.0, 3.0};
    const SelectResult r = select_region(m, o);
    EXPECT_TRUE(same_bytes(r.mMesh.Points(), m.Points()));
    ASSERT_EQ(r.mMesh.NumCellBlocks(), m.NumCellBlocks());
    EXPECT_TRUE(same_bytes(r.mMesh.Cells(0).Conn(), m.Cells(0).Conn()));
    // The side region of the same name survives; a cell region was added.
    EXPECT_EQ(r.mMesh.NumRegions(), m.NumRegions() + 1);
    EXPECT_NE(r.mMesh.FindRegion("xmin", RegionKind::Side), Mesh::npos);
    // Selecting again under the same name and kind replaces, even with a tag.
    o.mTag = 7;
    o.mBoxHi = {2.5, 3.0, 3.0};
    const SelectResult again = select_region(r.mMesh, o);
    EXPECT_EQ(again.mMesh.NumRegions(), m.NumRegions() + 1);
    const std::size_t i = again.mMesh.FindRegion("xmin", RegionKind::Cell);
    ASSERT_NE(i, Mesh::npos);
    EXPECT_EQ(again.mMesh.Region(i).NumEntries(), 8u);
    EXPECT_EQ(again.mMesh.Region(i).mTag, 7);
}

TEST(Select, FeedsFindInterface) {
    const Mesh m = cube();
    SelectOptions o = named("left", RegionKind::Cell);
    o.mMode = SelectMode::Centroid;
    o.mPlane = true;
    o.mPlanePoint = {1.0, 0.0, 0.0};
    o.mPlaneNormal = {-1.0, 0.0, 0.0};
    Mesh step = select_region(m, o).mMesh;
    o.mName = "right";
    o.mPlaneNormal = {1.0, 0.0, 0.0};
    step = select_region(step, o).mMesh;
    RegionSelector a, b;
    a.mName = "left";
    a.mKind = static_cast<std::int32_t>(RegionKind::Cell);
    b.mName = "right";
    b.mKind = static_cast<std::int32_t>(RegionKind::Cell);
    const FindInterfaceResult fi = find_interface(step, a, b, {});
    EXPECT_EQ(fi.mReport.mNumPairs, 4);
}

TEST(Select, RefusesBadOptions) {
    const Mesh m = cube();
    SelectOptions o;
    o.mBox = true;
    EXPECT_THROW(select_region(m, o), std::invalid_argument);  // no name
    o = named("x", RegionKind::Cell);
    EXPECT_THROW(select_region(m, o), std::invalid_argument);  // no criterion
    o.mSphere = true;
    o.mSphereRadius = -1.0;
    EXPECT_THROW(select_region(m, o), std::invalid_argument);
    o = named("x", RegionKind::Cell);
    o.mPlane = true;
    o.mPlaneNormal = {0.0, 0.0, 0.0};
    EXPECT_THROW(select_region(m, o), std::invalid_argument);
    o = named("x", RegionKind::Side);
    o.mSeed = true;
    o.mFeatureAngleDeg = 200.0;
    EXPECT_THROW(select_region(m, o), std::invalid_argument);
    o.mKind = RegionKind::Side;
    o.mFeatureAngleDeg = 30.0;
    o.mSeed = false;
    o.mBox = true;
    EXPECT_THROW(select_region(line(4), o), std::invalid_argument);  // no faces or edges
    EXPECT_EQ(select_mode_from_name("centroid"), SelectMode::Centroid);
    EXPECT_STREQ(select_mode_name(SelectMode::Any), "any");
    EXPECT_THROW(select_mode_from_name("most"), std::invalid_argument);
}

TEST(Select, PolygonsAndPolyhedra) {
    // Polygon edges are sides of a polygon surface.
    const Mesh hex = honeycomb({{2, 2}});
    SelectOptions o = named("rim", RegionKind::Side);
    o.mBox = true;
    o.mBoxLo = {-100.0, -100.0, -100.0};
    o.mBoxHi = {100.0, 100.0, 100.0};
    const SelectResult r = select_region(hex, o);
    EXPECT_GT(r.mNumSelected, 0);
    EXPECT_LT(r.mNumSelected, r.mNumCandidates + 1);
    EXPECT_EQ(r.mRegion.mDim, 1);
    // Growing a polygon surface takes it whole.
    SelectOptions g = named("all", RegionKind::Cell);
    g.mSeed = true;
    EXPECT_EQ(select_region(hex, g).mNumSelected, total_cells(hex));
}
