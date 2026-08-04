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
// Tests for the convert_cells operation (linearize / simplexify / elevate),
// including the orientation and volume-conservation invariants that pin the
// decomposition templates.

// System includes
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/stats.hpp"

namespace {

using meshioplusplus::compute_stats;
using meshioplusplus::convert_cells;
using meshioplusplus::convert_cells_mode_from_name;
using meshioplusplus::ConvertCellsMode;
using meshioplusplus::ConvertCellsOptions;
using meshioplusplus::ConvertCellsResult;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::StatsReport;

ConvertCellsOptions opts(ConvertCellsMode mode, bool record = false) {
    ConvertCellsOptions o;
    o.mMode = mode;
    o.mRecordParentIds = record;
    return o;
}

// The signed volume of the tetrahedron (a, b, c, d), positive when d lies on
// the side the (a, b, c) right-hand normal points to.
double signed_tet_volume(const NDArray& rPoints, std::size_t Dim, std::int64_t a, std::int64_t b,
                         std::int64_t c, std::int64_t d) {
    auto at = [&](std::int64_t i, std::size_t k) {
        return rPoints.As<double>()[static_cast<std::size_t>(i) * Dim + k];
    };
    const double ux = at(b, 0) - at(a, 0), uy = at(b, 1) - at(a, 1), uz = at(b, 2) - at(a, 2);
    const double vx = at(c, 0) - at(a, 0), vy = at(c, 1) - at(a, 1), vz = at(c, 2) - at(a, 2);
    const double wx = at(d, 0) - at(a, 0), wy = at(d, 1) - at(a, 1), wz = at(d, 2) - at(a, 2);
    const double cx = uy * vz - uz * vy;
    const double cy = uz * vx - ux * vz;
    const double cz = ux * vy - uy * vx;
    return (cx * wx + cy * wy + cz * wz) / 6.0;
}

// The connectivity rows of a mesh's single cell block.
std::vector<std::vector<std::int64_t>> rows_of(const Mesh& rMesh, std::size_t Block = 0) {
    const auto cb = rMesh.Cells(Block);
    const std::size_t npc = cb.NodesPerCell();
    const std::int64_t* conn = cb.Conn().As<std::int64_t>();
    std::vector<std::vector<std::int64_t>> out;
    for (std::size_t c = 0; c < cb.NumCells(); ++c)
        out.emplace_back(conn + c * npc, conn + (c + 1) * npc);
    return out;
}

// --- mode parsing ------------------------------------------------------------

TEST(ConvertCells, ModeFromName) {
    EXPECT_EQ(convert_cells_mode_from_name("linearize"), ConvertCellsMode::Linearize);
    EXPECT_EQ(convert_cells_mode_from_name("simplexify"), ConvertCellsMode::Simplexify);
    EXPECT_EQ(convert_cells_mode_from_name("elevate"), ConvertCellsMode::Elevate);
    EXPECT_THROW(convert_cells_mode_from_name("nope"), std::invalid_argument);
}

// --- linearize ---------------------------------------------------------------

TEST(ConvertCells, LinearizeTetra10KeepsCornersAndPrunesMidNodes) {
    Mesh m = mt::tet10_mesh();
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Linearize));

    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "tetra");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 1u);
    // The 6 mid-edge nodes are unreferenced and must be pruned.
    EXPECT_EQ(r.mMesh.NumPoints(), 4u);

    // Corner connectivity survives verbatim (through the point remap).
    const std::vector<std::vector<std::int64_t>> conn = rows_of(r.mMesh);
    ASSERT_EQ(conn.size(), 1u);
    EXPECT_EQ(conn[0], (std::vector<std::int64_t>{0, 1, 2, 3}));

    // The point map sends the 4 corners to 0..3 and the mid nodes to -1.
    const std::int64_t* pm = r.mPointMap.As<std::int64_t>();
    for (std::int64_t i = 0; i < 4; ++i)
        EXPECT_EQ(pm[i], i);
    for (std::size_t i = 4; i < 10; ++i)
        EXPECT_EQ(pm[i], -1);
}

TEST(ConvertCells, LinearizeLeavesLinearCellsAlone) {
    Mesh m = mt::hex_mesh();
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Linearize));
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "hexahedron");
    EXPECT_EQ(r.mMesh.NumPoints(), m.NumPoints());
    EXPECT_EQ(rows_of(r.mMesh), rows_of(m));
}

TEST(ConvertCells, LinearizeCarriesPointDataThroughThePrune) {
    Mesh m = mt::tet10_mesh();
    NDArray t = NDArray::Uninit(meshioplusplus::DType::Float64, {10});
    for (std::size_t i = 0; i < 10; ++i)
        t.As<double>()[i] = static_cast<double>(i);
    m.AddPointData("T", std::move(t));

    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Linearize));
    ASSERT_TRUE(r.mMesh.HasPointData("T"));
    const NDArray& got = r.mMesh.PointData("T");
    ASSERT_EQ(got.Size(), 4u);
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_DOUBLE_EQ(got.As<double>()[i], static_cast<double>(i));
}

// --- simplexify --------------------------------------------------------------

TEST(ConvertCells, SimplexifyQuadGivesTwoTriangles) {
    Mesh m = mt::quad_mesh();
    const std::size_t nquads = m.Cells(0).NumCells();
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Simplexify));

    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "triangle");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 2 * nquads);
    // Points are untouched by the decomposition.
    EXPECT_EQ(r.mMesh.NumPoints(), m.NumPoints());
}

TEST(ConvertCells, SimplexifyHexGivesSixTetraConservingVolume) {
    Mesh m = mt::hex_mesh();
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Simplexify));

    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "tetra");
    ASSERT_EQ(r.mMesh.Cells(0).NumCells(), 6u);
    EXPECT_EQ(r.mMesh.NumPoints(), 8u);

    // Orientation invariant: every child is positively oriented, and the six
    // volumes sum to the unit cube's volume.
    const NDArray& pts = r.mMesh.Points();
    const std::size_t dim = pts.Shape()[1];
    double total = 0.0;
    for (const std::vector<std::int64_t>& row : rows_of(r.mMesh)) {
        const double v = signed_tet_volume(pts, dim, row[0], row[1], row[2], row[3]);
        EXPECT_GT(v, 0.0) << "simplexify emitted an inverted tetrahedron";
        total += v;
    }
    EXPECT_NEAR(total, 1.0, 1e-12);
}

TEST(ConvertCells, SimplexifyWedgeAndPyramidAreOrientedAndConserveVolume) {
    struct Case {
        Mesh mMesh;
        std::size_t mExpectedChildren;
        double mExpectedVolume;
    };
    std::vector<Case> cases;
    cases.push_back({mt::wedge_mesh(), 3, 0.5});
    // A unit-base pyramid of height 1: volume = 1/3.
    cases.push_back({mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}},
                                   "pyramid", {{0, 1, 2, 3, 4}}),
                     2, 1.0 / 3.0});

    for (Case& c : cases) {
        ConvertCellsResult r = convert_cells(c.mMesh, opts(ConvertCellsMode::Simplexify));
        ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
        EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "tetra");
        ASSERT_EQ(r.mMesh.Cells(0).NumCells(), c.mExpectedChildren);

        const NDArray& pts = r.mMesh.Points();
        const std::size_t dim = pts.Shape()[1];
        double total = 0.0;
        for (const std::vector<std::int64_t>& row : rows_of(r.mMesh)) {
            const double v = signed_tet_volume(pts, dim, row[0], row[1], row[2], row[3]);
            EXPECT_GT(v, 0.0) << "simplexify emitted an inverted tetrahedron";
            total += v;
        }
        EXPECT_NEAR(total, c.mExpectedVolume, 1e-12);
    }
}

TEST(ConvertCells, SimplexifyReportsNoInvertedCells) {
    Mesh m = mt::hex_mesh();
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Simplexify));
    const StatsReport s = compute_stats(r.mMesh);
    EXPECT_EQ(s.mNumInverted, 0);
    EXPECT_NEAR(s.mUnsignedVolume, 1.0, 1e-12);
}

TEST(ConvertCells, SimplexifyReplicatesCellDataAndRecordsParents) {
    Mesh m = mt::hex_mesh();
    m.AddCellData("mat", {mt::int_data_array({7})});

    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Simplexify, /*record=*/true));

    // Exactly one cell_data array per output block, sized to the output cells.
    ASSERT_TRUE(r.mMesh.HasCellData("mat"));
    ASSERT_EQ(r.mMesh.CellDataNumBlocks("mat"), r.mMesh.NumCellBlocks());
    const NDArray& mat = r.mMesh.CellData("mat", 0);
    ASSERT_EQ(mat.Shape()[0], 6u);
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(mat, i), 7);

    ASSERT_TRUE(r.mMesh.HasCellData("convert:parent_cell"));
    const NDArray& parents = r.mMesh.CellData("convert:parent_cell", 0);
    ASSERT_EQ(parents.Shape()[0], 6u);
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_EQ(parents.As<std::int64_t>()[i], 0);

    // The cell map points at each parent's first child.
    ASSERT_EQ(r.mCellMaps.size(), 1u);
    EXPECT_EQ(r.mCellMaps[0].As<std::int64_t>()[0], 0);
}

TEST(ConvertCells, SimplexifyPolyhedronFansIntoTetrahedra) {
    // Since v9.17.0 a polyhedron decomposes rather than raising: one tet per
    // (face, edge-of-that-face), from the face's corner average to the cell's.
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));
    m.AddPolyhedronBlock("polyhedron4", {{{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}});
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Simplexify));
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "tetra");
    // 4 faces x 3 edges each = 12 children; 1 cell centroid + 4 face centroids.
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 12u);
    EXPECT_EQ(r.mMesh.NumPoints(), 4u + 5u);
    // Block structure stays 1:1 and the map is a contiguous FirstChild run.
    ASSERT_EQ(r.mCellMaps.size(), 1u);
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMaps[0], 0), 0);
}

TEST(ConvertCells, SimplexifyingAPolyhedronConservesVolumeExactly) {
    // The oracle the shared fan buys: the children's total volume is the
    // parent's, both being literally the same sum of the same tetrahedra. A
    // single-centroid decomposition would not conserve it on a warped cell.
    const std::vector<std::vector<double>> pts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0},       {0, 1, 0},
                                                  {0, 0, 1}, {1, 0, 1}, {1.6, 1.4, 1.3}, {0, 1, 1}};
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock(
        "polyhedron8",
        {{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}}});
    const double before = compute_stats(m).mSignedVolume;
    const Mesh simp = convert_cells(m, opts(ConvertCellsMode::Simplexify)).mMesh;
    const double after = compute_stats(simp).mSignedVolume;
    EXPECT_GT(before, 0.0);
    EXPECT_NEAR(after, before, 1e-14) << "the decomposition did not conserve volume";

    // Every child must be positively oriented, or the "same sum" claim is a
    // coincidence of cancelling signs.
    EXPECT_EQ(compute_stats(simp).mNumInverted, 0);
}

TEST(ConvertCells, SimplexifyPolygonFansIntoTriangles) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0.5, 1.5, 0}, {0, 1, 0}}));
    m.AddPolygonBlock("polygon", {{0, 1, 2, 3, 4}});
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Simplexify));
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "triangle");
    // An n-gon fans into n - 2 triangles.
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 3u);
}

TEST(ConvertCells, SimplexifyRectangularPolygonBlockAlsoFans) {
    // A polygon block whose cells all have the same node count is stored
    // rectangularly, not jagged -- it must still fan.
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1.5, 1, 0}, {0.5, 1.5, 0}, {-0.5, 1, 0}},
                           "polygon", {{0, 1, 2, 3, 4}});
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Simplexify));
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "triangle");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 3u);
}

// --- elevate -----------------------------------------------------------------

TEST(ConvertCells, ElevateTriangleAddsEdgeMidpoints) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {2, 0, 0}, {0, 2, 0}}, "triangle", {{0, 1, 2}});
    NDArray t = NDArray::Uninit(meshioplusplus::DType::Float64, {3});
    t.As<double>()[0] = 0.0;
    t.As<double>()[1] = 10.0;
    t.As<double>()[2] = 20.0;
    m.AddPointData("T", std::move(t));

    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Elevate));

    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "triangle6");
    ASSERT_EQ(r.mMesh.NumPoints(), 6u);

    const std::vector<std::vector<std::int64_t>> conn = rows_of(r.mMesh);
    ASSERT_EQ(conn.size(), 1u);
    EXPECT_EQ(conn[0], (std::vector<std::int64_t>{0, 1, 2, 3, 4, 5}));

    // Node 3 = mid(0,1), node 4 = mid(1,2), node 5 = mid(2,0).
    const double* p = r.mMesh.Points().As<double>();
    const std::size_t dim = r.mMesh.Points().Shape()[1];
    EXPECT_DOUBLE_EQ(p[3 * dim + 0], 1.0);
    EXPECT_DOUBLE_EQ(p[3 * dim + 1], 0.0);
    EXPECT_DOUBLE_EQ(p[4 * dim + 0], 1.0);
    EXPECT_DOUBLE_EQ(p[4 * dim + 1], 1.0);
    EXPECT_DOUBLE_EQ(p[5 * dim + 0], 0.0);
    EXPECT_DOUBLE_EQ(p[5 * dim + 1], 1.0);

    // point_data on a new node is the mean of its edge endpoints.
    const NDArray& got = r.mMesh.PointData("T");
    ASSERT_EQ(got.Size(), 6u);
    EXPECT_DOUBLE_EQ(got.As<double>()[3], 5.0);   // mean(0, 10)
    EXPECT_DOUBLE_EQ(got.As<double>()[4], 15.0);  // mean(10, 20)
    EXPECT_DOUBLE_EQ(got.As<double>()[5], 10.0);  // mean(20, 0)
}

TEST(ConvertCells, ElevateSharesMidNodesBetweenAdjacentCells) {
    // Two triangles sharing the edge (1, 2): 5 unique edges -> 5 new nodes.
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}, "triangle",
                           {{0, 1, 2}, {1, 3, 2}});
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Elevate));
    EXPECT_EQ(r.mMesh.NumPoints(), 4u + 5u);

    const std::vector<std::vector<std::int64_t>> conn = rows_of(r.mMesh);
    ASSERT_EQ(conn.size(), 2u);
    // The shared edge (1,2) is cell 0's edge 1 and cell 1's edge 2.
    EXPECT_EQ(conn[0][4], conn[1][5]);
}

TEST(ConvertCells, ElevateHexGivesHexahedron20) {
    Mesh m = mt::hex_mesh();
    ConvertCellsResult r = convert_cells(m, opts(ConvertCellsMode::Elevate));
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "hexahedron20");
    EXPECT_EQ(r.mMesh.Cells(0).NodesPerCell(), 20u);
    EXPECT_EQ(r.mMesh.NumPoints(), 8u + 12u);
}

TEST(ConvertCells, ElevateFullLagrangeTargetThrows) {
    Mesh m = mt::make_mesh({{0, 0, 0},
                            {1, 0, 0},
                            {1, 1, 0},
                            {0, 1, 0},
                            {0.5, 0, 0},
                            {1, 0.5, 0},
                            {0.5, 1, 0},
                            {0, 0.5, 0},
                            {0.5, 0.5, 0}},
                           "quad9", {{0, 1, 2, 3, 4, 5, 6, 7, 8}});
    EXPECT_THROW(convert_cells(m, opts(ConvertCellsMode::Elevate)), std::invalid_argument);
}

// --- round-trip identities ---------------------------------------------------

TEST(ConvertCells, ElevateThenLinearizeReturnsTheOriginal) {
    Mesh m = mt::tet_mesh();
    ConvertCellsResult up = convert_cells(m, opts(ConvertCellsMode::Elevate));
    ConvertCellsResult down = convert_cells(up.mMesh, opts(ConvertCellsMode::Linearize));

    EXPECT_EQ(down.mMesh.NumPoints(), m.NumPoints());
    ASSERT_EQ(down.mMesh.NumCellBlocks(), m.NumCellBlocks());
    EXPECT_EQ(std::string(down.mMesh.Cells(0).Type()), std::string(m.Cells(0).Type()));
    EXPECT_EQ(rows_of(down.mMesh), rows_of(m));

    const double* a = m.Points().As<double>();
    const double* b = down.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < m.Points().Size(); ++i)
        EXPECT_DOUBLE_EQ(a[i], b[i]);
}

TEST(ConvertCells, LinearizeThenElevateRestoresTheTopology) {
    Mesh m = mt::tet10_mesh();
    ConvertCellsResult down = convert_cells(m, opts(ConvertCellsMode::Linearize));
    ConvertCellsResult up = convert_cells(down.mMesh, opts(ConvertCellsMode::Elevate));

    EXPECT_EQ(std::string(up.mMesh.Cells(0).Type()), std::string(m.Cells(0).Type()));
    EXPECT_EQ(up.mMesh.Cells(0).NumCells(), m.Cells(0).NumCells());
    EXPECT_EQ(up.mMesh.NumPoints(), m.NumPoints());
}

// --- determinism -------------------------------------------------------------

TEST(ConvertCells, ElevateNumberingIsStableAcrossRepeatedRuns) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 0, 0}, {2, 1, 0}},
                           "triangle", {{0, 1, 2}, {1, 3, 2}, {1, 4, 3}, {4, 5, 3}});
    ConvertCellsResult a = convert_cells(m, opts(ConvertCellsMode::Elevate));
    for (int i = 0; i < 4; ++i) {
        ConvertCellsResult b = convert_cells(m, opts(ConvertCellsMode::Elevate));
        EXPECT_EQ(rows_of(a.mMesh), rows_of(b.mMesh));
        ASSERT_EQ(a.mMesh.Points().Size(), b.mMesh.Points().Size());
        for (std::size_t k = 0; k < a.mMesh.Points().Size(); ++k)
            EXPECT_DOUBLE_EQ(a.mMesh.Points().As<double>()[k], b.mMesh.Points().As<double>()[k]);
    }
}

}  // namespace
