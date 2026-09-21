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

// System includes
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/region.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

using mt::icosphere;
using mt::make_mesh;

// A unit cube of six outward-wound quads over 8 points.
Mesh cube_quads() {
    return make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "quad",
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 4, 7, 3}, {1, 2, 6, 5}});
}

NDArray i64(const std::vector<std::int64_t>& rVals) {
    NDArray a = NDArray::Uninit(DType::Int64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

const double* normals_of(const Mesh& rMesh) {
    return rMesh.PointData(kNormalsName).As<double>();
}

double length_of(const double* pN) {
    return std::sqrt(pN[0] * pN[0] + pN[1] * pN[1] + pN[2] * pN[2]);
}

// The normal is one of the six axis directions, exactly.
bool is_axis_aligned(const double* pN) {
    int nonzero = 0;
    for (int k = 0; k < 3; ++k)
        if (pN[k] != 0.0) {
            ++nonzero;
            if (std::abs(pN[k]) != 1.0)
                return false;
        }
    return nonzero == 1;
}

}  // namespace

TEST(Normals, FlatSquareHasUpNormals) {
    NormalsOptions o;
    o.mCellNormals = true;
    const NormalsResult r = compute_normals(mt::tri_mesh(), o);
    const double* n = normals_of(r.mMesh);
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        EXPECT_EQ(n[i * 3 + 0], 0.0);
        EXPECT_EQ(n[i * 3 + 1], 0.0);
        EXPECT_EQ(n[i * 3 + 2], 1.0);
    }
    const NDArray& cn = r.mMesh.CellData(kNormalsName, 0);
    ASSERT_EQ(cn.Size(), 6u);
    EXPECT_EQ(cn.As<double>()[2], 1.0);
    EXPECT_EQ(cn.As<double>()[5], 1.0);
    EXPECT_EQ(r.mNumIsolated, 0);
    EXPECT_EQ(r.mNumAddedPoints, 0);
    EXPECT_TRUE(r.mQuality.mBoundaryEdges > 0);
}

TEST(Normals, TwoDimensionalPointsGetAnUpNormal) {
    const NormalsResult r = compute_normals(mt::tri_mesh_2d());
    const NDArray& a = r.mMesh.PointData(kNormalsName);
    ASSERT_EQ(a.Shape().size(), 2u);
    EXPECT_EQ(a.Shape()[1], 3u);
    EXPECT_EQ(a.As<double>()[2], 1.0);
}

TEST(Normals, SmoothCubeNormalsPointOutwardAndAreUnit) {
    const Mesh cube = cube_quads();
    const NormalsResult r = compute_normals(cube);
    EXPECT_EQ(r.mMesh.NumPoints(), 8u);
    const double* n = normals_of(r.mMesh);
    const NDArray& p = cube.Points();
    for (std::size_t i = 0; i < 8; ++i) {
        EXPECT_NEAR(length_of(n + i * 3), 1.0, 1e-14);
        for (int k = 0; k < 3; ++k)
            EXPECT_GT(n[i * 3 + k] * (2.0 * p.As<double>()[i * 3 + k] - 1.0), 0.0)
                << "point " << i << " component " << k;
    }
}

TEST(Normals, SplittingACubeGivesAxisAlignedNormals) {
    NormalsOptions o;
    o.mSplit = true;
    o.mSplitAngle = 30.0;
    o.mRecordParentIds = true;
    o.mCellNormals = true;
    const NormalsResult r = compute_normals(cube_quads(), o);
    ASSERT_EQ(r.mMesh.NumPoints(), 24u);
    EXPECT_EQ(r.mNumAddedPoints, 16);
    EXPECT_EQ(r.mNumSplitPoints, 8);

    const double* n = normals_of(r.mMesh);
    for (std::size_t i = 0; i < 24; ++i)
        EXPECT_TRUE(is_axis_aligned(n + i * 3)) << "point " << i;

    // Every cell's corners now sit on points whose normal is the cell's own.
    const NDArray& conn = r.mMesh.Cells(0).Conn();
    const double* cn = r.mMesh.CellData(kNormalsName, 0).As<double>();
    for (std::size_t c = 0; c < 6; ++c)
        for (std::size_t k = 0; k < 4; ++k) {
            const std::int64_t id = conn.As<std::int64_t>()[c * 4 + k];
            for (int d = 0; d < 3; ++d)
                EXPECT_EQ(n[static_cast<std::size_t>(id) * 3 + d], cn[c * 3 + d])
                    << "cell " << c << " corner " << k;
        }

    // The originals keep their indices; every copy names its parent.
    const std::int64_t* parent = r.mMesh.PointData(kNormalsParentPointName).As<std::int64_t>();
    for (std::size_t i = 0; i < 8; ++i)
        EXPECT_EQ(parent[i], static_cast<std::int64_t>(i));
    const double* q = r.mMesh.Points().As<double>();
    for (std::size_t i = 8; i < 24; ++i)
        for (int d = 0; d < 3; ++d)
            EXPECT_EQ(q[i * 3 + d], q[static_cast<std::size_t>(parent[i]) * 3 + d]);
}

TEST(Normals, AFullSplitAngleOnAClosedManifoldIsTheUnsplitResult) {
    const Mesh cube = cube_quads();
    const NormalsResult smooth = compute_normals(cube);
    NormalsOptions o;
    o.mSplit = true;
    o.mSplitAngle = 180.0;
    const NormalsResult wide = compute_normals(cube, o);
    ASSERT_EQ(wide.mMesh.NumPoints(), smooth.mMesh.NumPoints());
    EXPECT_EQ(wide.mNumAddedPoints, 0);
    EXPECT_EQ(std::memcmp(normals_of(wide.mMesh), normals_of(smooth.mMesh), 8 * 3 * sizeof(double)),
              0);
}

TEST(Normals, ASphereHasRadialNormalsUnderBothWeights) {
    for (SdfPseudonormalWeight w : {SdfPseudonormalWeight::Angle, SdfPseudonormalWeight::Area}) {
        NormalsOptions o;
        o.mWeight = w;
        const Mesh sphere = icosphere(3, 2.0);
        const NormalsResult r = compute_normals(sphere, o);
        const double* n = normals_of(r.mMesh);
        const double* p = sphere.Points().As<double>();
        for (std::size_t i = 0; i < sphere.NumPoints(); ++i) {
            const double len = length_of(p + i * 3);
            for (int d = 0; d < 3; ++d)
                EXPECT_NEAR(n[i * 3 + d], p[i * 3 + d] / len, 1e-2);
        }
        EXPECT_TRUE(r.mQuality.mWatertight);
    }
}

TEST(Normals, ASphereBelowTheSplitAngleIsNotSplit) {
    NormalsOptions o;
    o.mSplit = true;
    o.mSplitAngle = 30.0;
    const Mesh sphere = icosphere(3, 1.0);
    EXPECT_EQ(compute_normals(sphere, o).mNumAddedPoints, 0);
}

TEST(Normals, TwoTrianglesTouchingAtAVertexSplitThere) {
    const Mesh m = make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}}, "triangle",
                             {{0, 1, 2}, {0, 3, 4}});
    EXPECT_EQ(compute_normals(m).mMesh.NumPoints(), 5u);
    NormalsOptions o;
    o.mSplit = true;
    o.mSplitAngle = 180.0;
    const NormalsResult r = compute_normals(m, o);
    EXPECT_EQ(r.mNumAddedPoints, 1);
    EXPECT_EQ(r.mNumSplitPoints, 1);
}

TEST(Normals, AnInconsistentPairIsReportedAndAlwaysSplit) {
    // Two triangles sharing edge (1, 2), both wound the same way round it.
    const Mesh m =
        make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}, "triangle", {{0, 1, 2}, {3, 1, 2}});
    NormalsOptions o;
    o.mSplit = true;
    o.mSplitAngle = 180.0;
    const NormalsResult r = compute_normals(m, o);
    EXPECT_EQ(r.mQuality.mInconsistentPairs, 1);
    EXPECT_EQ(r.mNumSplitPoints, 2);
    EXPECT_EQ(r.mNumAddedPoints, 2);
}

TEST(Normals, ANonPlanarQuadHasItsNewellNormal) {
    NormalsOptions o;
    o.mCellNormals = true;
    const Mesh m = make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 1}, {0, 1, 0}}, "quad", {{0, 1, 2, 3}});
    const NormalsResult r = compute_normals(m, o);
    const double* c = r.mMesh.CellData(kNormalsName, 0).As<double>();
    const double s = 1.0 / std::sqrt(6.0);
    EXPECT_NEAR(c[0], -s, 1e-15);
    EXPECT_NEAR(c[1], -s, 1e-15);
    EXPECT_NEAR(c[2], 2.0 * s, 1e-15);
}

TEST(Normals, ARaggedPolygonBlockSplitsLikeAnyOther) {
    // A planar pentagon and a triangle folded 90 degrees along the edge (0, 4).
    Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {1.5, 1, 0}, {0.5, 2, 0}, {-0.5, 1, 0}, {0, 0, 1}}));
    m.AddPolygonBlock("polygon", {{0, 1, 2, 3, 4}});
    m.AddCellBlock("triangle", mt::conn_from({{0, 5, 4}}));
    NormalsOptions o;
    o.mSplit = true;
    o.mSplitAngle = 30.0;
    o.mCellNormals = true;
    const NormalsResult r = compute_normals(m, o);
    // Edge (0, 4) is shared and folds by 90 degrees, so points 0 and 4 split;
    // the pentagon's own fan diagonals never do.
    EXPECT_EQ(r.mNumAddedPoints, 2);
    const double* pn = normals_of(r.mMesh);
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
        EXPECT_NEAR(length_of(pn + i * 3), 1.0, 1e-14);
}

TEST(Normals, LinesAndVerticesAreIgnored) {
    NormalsOptions o;
    o.mCellNormals = true;
    const NormalsResult r = compute_normals(mt::line_mesh(), o);
    const double* n = normals_of(r.mMesh);
    EXPECT_TRUE(std::isnan(n[0]));
    EXPECT_EQ(r.mNumIsolated, static_cast<std::int64_t>(r.mMesh.NumPoints()));
    EXPECT_TRUE(std::isnan(r.mMesh.CellData(kNormalsName, 0).As<double>()[0]));
}

TEST(Normals, RefusesVolumesAndHigherOrderSurfaces) {
    try {
        compute_normals(mt::tet_mesh());
        FAIL() << "a tetra block must be refused";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("extract_surface"), std::string::npos);
    }
    try {
        compute_normals(mt::triangle6_mesh());
        FAIL() << "a triangle6 block must be refused";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("linearize"), std::string::npos);
    }
}

TEST(Normals, RejectsABadSplitAngleAndAnUnknownRegion) {
    NormalsOptions o;
    o.mSplit = true;
    o.mSplitAngle = 181.0;
    EXPECT_THROW(compute_normals(mt::tri_mesh(), o), std::invalid_argument);
    o.mSplitAngle = -1.0;
    EXPECT_THROW(compute_normals(mt::tri_mesh(), o), std::invalid_argument);
    NormalsOptions r;
    r.mRegion = "nowhere";
    EXPECT_THROW(compute_normals(mt::tri_mesh(), r), std::invalid_argument);
}

TEST(Normals, ARegionRestrictsTheSurface) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("first", RegionKind::Cell, i64({0})));
    NormalsOptions o;
    o.mRegion = "first";
    const NormalsResult r = compute_normals(m, o);
    const double* n = normals_of(r.mMesh);
    // Cell 0 is (0, 1, 2): point 3 is not touched, so its normal is NaN.
    EXPECT_EQ(n[2], 1.0);
    EXPECT_TRUE(std::isnan(n[3 * 3 + 2]));
    EXPECT_EQ(r.mNumIsolated, 1);
}

TEST(Normals, SplittingCarriesDataAndRegionsAlong) {
    Mesh m = cube_quads();
    NDArray temp(DType::Float64, {std::size_t{8}, std::size_t{1}});
    for (std::size_t i = 0; i < 8; ++i)
        temp.As<double>()[i] = 10.0 * static_cast<double>(i);
    m.AddPointData("T", std::move(temp));
    std::vector<NDArray> tags;
    tags.push_back(NDArray(DType::Float64, {std::size_t{6}}));
    for (std::size_t c = 0; c < 6; ++c)
        tags[0].As<double>()[c] = static_cast<double>(c);
    m.AddCellData("id", std::move(tags));
    m.AddRegion(Region("corner", RegionKind::Point, i64({0})));
    m.AddRegion(Region("top", RegionKind::Cell, i64({1})));

    NormalsOptions o;
    o.mSplit = true;
    o.mRecordParentIds = true;
    const NormalsResult r = compute_normals(m, o);
    ASSERT_EQ(r.mMesh.NumPoints(), 24u);

    const double* t = r.mMesh.PointData("T").As<double>();
    const std::int64_t* parent = r.mMesh.PointData(kNormalsParentPointName).As<std::int64_t>();
    for (std::size_t i = 0; i < 24; ++i)
        EXPECT_EQ(t[i], 10.0 * static_cast<double>(parent[i]));

    // Point 0 and both of its copies belong to the "corner" region.
    const std::size_t idx = r.mMesh.FindRegion("corner", RegionKind::Point);
    ASSERT_NE(idx, Mesh::npos);
    EXPECT_EQ(r.mMesh.Region(idx).NumEntries(), 3u);
    EXPECT_NE(r.mMesh.FindRegion("top", RegionKind::Cell), Mesh::npos);
    EXPECT_EQ(r.mMesh.CellData("id", 0).As<double>()[3], 3.0);
}

TEST(Normals, AnExistingNormalsArrayIsReplaced) {
    Mesh m = mt::tri_mesh();
    NDArray old(DType::Float64, {std::size_t{4}, std::size_t{3}});
    for (std::size_t i = 0; i < 12; ++i)
        old.As<double>()[i] = 7.0;
    m.AddPointData(kNormalsName, std::move(old));
    const NormalsResult r = compute_normals(m);
    EXPECT_EQ(normals_of(r.mMesh)[0], 0.0);
    EXPECT_EQ(normals_of(r.mMesh)[2], 1.0);
}

TEST(Normals, PointNormalsCanBeSwitchedOff) {
    NormalsOptions o;
    o.mPointNormals = false;
    o.mCellNormals = true;
    const NormalsResult r = compute_normals(mt::tri_mesh(), o);
    EXPECT_FALSE(r.mMesh.HasPointData(kNormalsName));
    EXPECT_TRUE(r.mMesh.HasCellData(kNormalsName));
}
