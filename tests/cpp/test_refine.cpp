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
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/surface.hpp"

namespace {

using meshioplusplus::compute_stats;
using meshioplusplus::extract_surface;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::refine;
using meshioplusplus::RefineOptions;
using meshioplusplus::RefineResult;
using meshioplusplus::StatsReport;

RefineOptions opts(int levels, bool record = false) {
    RefineOptions o;
    o.mLevels = levels;
    o.mRecordParentIds = record;
    return o;
}

// The signed volume of the tetrahedron (a, b, c, d), positive when d lies on
// the side the (a, b, c) right-hand normal points to. Mirrors the helper in
// test_convert_cells.cpp.
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

std::vector<std::vector<std::int64_t>> rows_of(const Mesh& rMesh, std::size_t Block = 0) {
    const auto cb = rMesh.Cells(Block);
    const std::size_t npc = cb.NodesPerCell();
    const std::int64_t* conn = cb.Conn().As<std::int64_t>();
    std::vector<std::vector<std::int64_t>> out;
    for (std::size_t c = 0; c < cb.NumCells(); ++c)
        out.emplace_back(conn + c * npc, conn + (c + 1) * npc);
    return out;
}

std::vector<double> point_at(const Mesh& rMesh, std::int64_t i) {
    const std::size_t dim = rMesh.PointDim();
    const double* p = rMesh.Points().As<double>();
    return std::vector<double>(p + static_cast<std::size_t>(i) * dim,
                               p + (static_cast<std::size_t>(i) + 1) * dim);
}

// --- shared tables -----------------------------------------------------------

// The quad-face rows of cell_subdivision.hpp must agree with cell_faces.hpp's
// own face-centre columns, or refine's new nodes land in the wrong slot of the
// full-Lagrange layout. This is the guard that keeps the two tables from
// drifting (see cell_subdivision.hpp's file comment).
TEST(CellSubdivision, QuadFacesAgreeWithCellFaces) {
    using meshioplusplus::CellType;
    namespace d = meshioplusplus::detail;
    struct Case {
        CellType mLinear;
        CellType mLagrange;
        std::uint8_t mFirstCentre;
    };
    for (const Case& c : std::vector<Case>{{CellType::Hexahedron, CellType::Hexahedron27, 20},
                                           {CellType::Wedge, CellType::Wedge18, 15}}) {
        const auto& quads = d::cell_refine_quad_faces(c.mLinear);
        // Every quad face of the Lagrange type, keyed by its sorted corners.
        std::map<std::set<std::uint8_t>, std::uint8_t> centre_of;
        for (const d::CellFaceDef& f : d::cell_faces(c.mLagrange)) {
            if (f.mNumCorners != 4)
                continue;
            centre_of[{f.mNodes[0], f.mNodes[1], f.mNodes[2], f.mNodes[3]}] = f.mNodes[8];
        }
        ASSERT_EQ(centre_of.size(), quads.size()) << cell_type_name(c.mLinear);
        for (std::size_t k = 0; k < quads.size(); ++k) {
            const std::set<std::uint8_t> key{quads[k][0], quads[k][1], quads[k][2], quads[k][3]};
            ASSERT_TRUE(centre_of.count(key)) << cell_type_name(c.mLinear) << " face " << k;
            EXPECT_EQ(centre_of[key], c.mFirstCentre + k)
                << cell_type_name(c.mLinear) << " face " << k << " maps to the wrong centre node";
        }
    }
}

// --- 2D templates ------------------------------------------------------------

TEST(Refine, TriangleSplitsIntoFourWithSharedMidEdgeNodes) {
    // Two triangles sharing edge (0, 2): 4 points, 5 unique edges.
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle",
                           {{0, 1, 2}, {0, 2, 3}});
    const RefineResult r = refine(m, opts(1));

    EXPECT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "triangle");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 8u);
    // 4 original + 5 unique edges -- NOT 4 + 2*3, which is what a per-cell
    // (undeduped) midpoint would give.
    EXPECT_EQ(r.mMesh.NumPoints(), 9u);

    // The shared edge's midpoint must be one node, referenced by both parents.
    const auto rows = rows_of(r.mMesh);
    std::set<std::int64_t> from_a(rows[0].begin(), rows[0].end());
    for (std::size_t i = 1; i < 4; ++i)
        from_a.insert(rows[i].begin(), rows[i].end());
    std::set<std::int64_t> from_b(rows[4].begin(), rows[4].end());
    for (std::size_t i = 5; i < 8; ++i)
        from_b.insert(rows[i].begin(), rows[i].end());
    std::vector<std::int64_t> shared;
    std::set_intersection(from_a.begin(), from_a.end(), from_b.begin(), from_b.end(),
                          std::back_inserter(shared));
    // corners 0 and 2, plus the midpoint of edge (0, 2).
    EXPECT_EQ(shared.size(), 3u);
}

TEST(Refine, TriangleMidNodesSitAtEdgeMidpoints) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {2, 0, 0}, {0, 2, 0}}, "triangle", {{0, 1, 2}});
    const RefineResult r = refine(m, opts(1));
    ASSERT_EQ(r.mMesh.NumPoints(), 6u);
    // New nodes are appended in cell_refine_edges order: (0,1), (1,2), (2,0).
    EXPECT_EQ(point_at(r.mMesh, 3), (std::vector<double>{1, 0, 0}));
    EXPECT_EQ(point_at(r.mMesh, 4), (std::vector<double>{1, 1, 0}));
    EXPECT_EQ(point_at(r.mMesh, 5), (std::vector<double>{0, 1, 0}));
}

TEST(Refine, TriangleInterpolatesLinearPointDataExactly) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {2, 0, 0}, {0, 2, 0}}, "triangle", {{0, 1, 2}});
    // f(x, y) = 3x + 5y + 1 -- a linear field is reproduced exactly by
    // midpoint averaging.
    m.AddPointData("f", mt::data_array({1.0, 7.0, 11.0}));
    const RefineResult r = refine(m, opts(1));
    const double* f = r.mMesh.PointData("f").As<double>();
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        const std::vector<double> p = point_at(r.mMesh, i);
        EXPECT_DOUBLE_EQ(f[i], 3.0 * p[0] + 5.0 * p[1] + 1.0) << "point " << i;
    }
}

TEST(Refine, TriangleConservesAreaAndCellCount) {
    Mesh m = mt::tri_mesh();
    const StatsReport before = compute_stats(m);
    const RefineResult r = refine(m, opts(2));
    const StatsReport after = compute_stats(r.mMesh);
    EXPECT_EQ(after.mNumCells, before.mNumCells * 16);  // 4^2
    EXPECT_NEAR(after.mTotalArea, before.mTotalArea, 1e-12);
}

TEST(Refine, QuadSplitsIntoFourWithCentreAtCornerMean) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}}, "quad", {{0, 1, 2, 3}});
    const RefineResult r = refine(m, opts(1));
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "quad");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 4u);
    // 4 corners + 4 edge mids + 1 face centre.
    ASSERT_EQ(r.mMesh.NumPoints(), 9u);
    EXPECT_EQ(point_at(r.mMesh, 8), (std::vector<double>{1, 1, 0}));
    const StatsReport after = compute_stats(r.mMesh);
    EXPECT_NEAR(after.mTotalArea, 4.0, 1e-12);
}

TEST(Refine, LineSplitsIntoTwo) {
    Mesh m = mt::line_mesh();
    const auto before = m.Cells(0).NumCells();
    const RefineResult r = refine(m, opts(1));
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "line");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), before * 2);
}

// --- 3D templates ------------------------------------------------------------

TEST(Refine, TetraSplitsIntoEightConservingVolume) {
    Mesh m = mt::tet_mesh();
    const StatsReport before = compute_stats(m);
    const RefineResult r = refine(m, opts(1));
    const StatsReport after = compute_stats(r.mMesh);

    EXPECT_EQ(r.mMesh.Cells(0).Type(), "tetra");
    EXPECT_EQ(after.mNumCells, before.mNumCells * 8);
    // A tetrahedron's subdivision is affine, so volume is conserved exactly.
    EXPECT_NEAR(after.mSignedVolume, before.mSignedVolume, 1e-12);
    EXPECT_EQ(after.mNumInverted, 0);
}

TEST(Refine, TetraChildrenAreAllPositivelyOriented) {
    Mesh m = mt::tet_mesh();
    const RefineResult r = refine(m, opts(2));
    const auto rows = rows_of(r.mMesh);
    const NDArray& pts = r.mMesh.Points();
    const std::size_t dim = r.mMesh.PointDim();
    for (std::size_t c = 0; c < rows.size(); ++c)
        EXPECT_GT(signed_tet_volume(pts, dim, rows[c][0], rows[c][1], rows[c][2], rows[c][3]), 0.0)
            << "child " << c << " is inverted";
}

TEST(Refine, HexSplitsIntoEightConservingVolume) {
    Mesh m = mt::hex_mesh();  // the unit cube: affine, so volume is exact
    const StatsReport before = compute_stats(m);
    const RefineResult r = refine(m, opts(1));
    const StatsReport after = compute_stats(r.mMesh);

    EXPECT_EQ(r.mMesh.Cells(0).Type(), "hexahedron");
    EXPECT_EQ(after.mNumCells, 8);
    // 8 corners + 12 edge mids + 6 face centres + 1 body = 27.
    EXPECT_EQ(r.mMesh.NumPoints(), 27u);
    EXPECT_NEAR(after.mSignedVolume, before.mSignedVolume, 1e-12);
    EXPECT_EQ(after.mNumInverted, 0);
}

TEST(Refine, HexCentreNodeIsTheCornerMean) {
    Mesh m = mt::hex_mesh();
    const RefineResult r = refine(m, opts(1));
    ASSERT_EQ(r.mMesh.NumPoints(), 27u);
    // The body centre is allocated last, after the 12 edge + 6 face nodes.
    EXPECT_EQ(point_at(r.mMesh, 26), (std::vector<double>{0.5, 0.5, 0.5}));
}

TEST(Refine, WedgeSplitsIntoEightConservingVolume) {
    Mesh m = mt::wedge_mesh();  // a right prism: affine, so volume is exact
    const StatsReport before = compute_stats(m);
    const RefineResult r = refine(m, opts(1));
    const StatsReport after = compute_stats(r.mMesh);

    EXPECT_EQ(r.mMesh.Cells(0).Type(), "wedge");
    EXPECT_EQ(after.mNumCells, 8);
    // 6 corners + 9 edge mids + 3 quad-face centres, no body node.
    EXPECT_EQ(r.mMesh.NumPoints(), 18u);
    EXPECT_NEAR(after.mSignedVolume, before.mSignedVolume, 1e-12);
    EXPECT_EQ(after.mNumInverted, 0);
}

// --- conformity --------------------------------------------------------------

// The direct test that quad-face centres are deduped: with a per-cell copy,
// the interior face between two hexahedra would appear on the boundary.
TEST(Refine, AdjacentHexesRefineConformingly) {
    Mesh m = mt::make_mesh({{0, 0, 0},
                            {1, 0, 0},
                            {1, 1, 0},
                            {0, 1, 0},
                            {0, 0, 1},
                            {1, 0, 1},
                            {1, 1, 1},
                            {0, 1, 1},
                            {0, 0, 2},
                            {1, 0, 2},
                            {1, 1, 2},
                            {0, 1, 2}},
                           "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}, {4, 5, 6, 7, 8, 9, 10, 11}});
    const RefineResult r = refine(m, opts(1));
    // 2x2x4 lattice of hexahedra = 3x3x5 nodes.
    EXPECT_EQ(r.mMesh.NumPoints(), 45u);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 16u);

    // A conforming 1x1x2 box refined once has 4*2 side faces per long side and
    // 4 per cap: 6 faces of a 2x2 grid on the caps, 2*(2x4) on the sides.
    const Mesh skin = extract_surface(m);
    const Mesh refined_skin = extract_surface(r.mMesh);
    std::size_t before = 0, after = 0;
    for (const auto cb : skin.CellRange())
        before += cb.NumCells();
    for (const auto cb : refined_skin.CellRange())
        after += cb.NumCells();
    // Every boundary facet splits into 4; no interior facet may appear.
    EXPECT_EQ(after, before * 4);
}

// --- levels, maps, data ------------------------------------------------------

TEST(Refine, TwoLevelsEqualsRefiningTwice) {
    Mesh m = mt::tri_quad_mesh();
    const RefineResult once = refine(m, opts(1));
    const RefineResult twice = refine(once.mMesh, opts(1));
    const RefineResult direct = refine(m, opts(2));

    ASSERT_EQ(direct.mMesh.NumCellBlocks(), twice.mMesh.NumCellBlocks());
    EXPECT_EQ(direct.mMesh.NumPoints(), twice.mMesh.NumPoints());
    for (std::size_t b = 0; b < direct.mMesh.NumCellBlocks(); ++b)
        EXPECT_EQ(rows_of(direct.mMesh, b), rows_of(twice.mMesh, b)) << "block " << b;
    mt::expect_points_close(direct.mMesh, twice.mMesh, 1e-15);
}

TEST(Refine, ZeroLevelsIsAnUnchangedClone) {
    Mesh m = mt::tri_quad_mesh();
    const RefineResult r = refine(m, opts(0));
    mt::expect_mesh_eq(m, r.mMesh);
}

TEST(Refine, PreservesBlockStructureAndReplicatesCellData) {
    Mesh m = mt::tri_quad_mesh();
    m.AddCellData(
        "tag", {mt::int_data_array({10, 20}), mt::int_data_array({30}), mt::int_data_array({40})});
    const RefineResult r = refine(m, opts(1));

    ASSERT_EQ(r.mMesh.NumCellBlocks(), 3u);
    ASSERT_EQ(r.mMesh.CellDataNumBlocks("tag"), 3u);
    // Block 0: two triangles -> 8 cells, each parent's row replicated 4x.
    // Read dtype-agnostically: NATIVE/KRATOS canonicalize integer cell_data to
    // Int64 at ingest, so the array's dtype is backend-dependent even though
    // its values are not.
    const auto& b0 = r.mMesh.CellData("tag", 0);
    ASSERT_EQ(b0.Size(), 8u);
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(b0, i), 10);
    for (std::size_t i = 4; i < 8; ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(b0, i), 20);
}

TEST(Refine, CellMapsPointAtContiguousChildren) {
    Mesh m = mt::tri_mesh();
    const RefineResult r = refine(m, opts(1));
    ASSERT_EQ(r.mCellMaps.size(), 1u);
    const std::int64_t* map = r.mCellMaps[0].As<std::int64_t>();
    for (std::size_t c = 0; c < r.mCellMaps[0].Size(); ++c)
        EXPECT_EQ(map[c], static_cast<std::int64_t>(c * 4));
}

TEST(Refine, PointMapIsTheIdentity) {
    Mesh m = mt::tri_mesh();
    const RefineResult r = refine(m, opts(2));
    ASSERT_EQ(r.mPointMap.Size(), m.NumPoints());
    for (std::size_t i = 0; i < r.mPointMap.Size(); ++i)
        EXPECT_EQ(r.mPointMap.As<std::int64_t>()[i], static_cast<std::int64_t>(i));
}

TEST(Refine, RecordParentIdsNamesTheOriginalAncestor) {
    Mesh m = mt::tri_mesh();  // 2 triangles
    const RefineResult r = refine(m, opts(2, /*record=*/true));
    ASSERT_TRUE(r.mMesh.HasCellData("refine:parent_cell"));
    const NDArray& a = r.mMesh.CellData("refine:parent_cell", 0);
    ASSERT_EQ(a.Size(), 32u);  // 2 parents * 16 grandchildren
    for (std::size_t i = 0; i < 16; ++i)
        EXPECT_EQ(a.As<std::int64_t>()[i], 0);
    for (std::size_t i = 16; i < 32; ++i)
        EXPECT_EQ(a.As<std::int64_t>()[i], 1);
}

// --- rejections --------------------------------------------------------------

TEST(Refine, HigherOrderCellsRaise) {
    EXPECT_THROW(refine(mt::tet10_mesh()), std::invalid_argument);
    EXPECT_THROW(refine(mt::triangle6_mesh()), std::invalid_argument);
    EXPECT_THROW(refine(mt::hex20_mesh()), std::invalid_argument);
}

TEST(Refine, PyramidRaises) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}}, "pyramid",
                           {{0, 1, 2, 3, 4}});
    EXPECT_THROW(refine(m), std::invalid_argument);
}

TEST(Refine, RaggedBlocksRaise) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 1.5, 0}}));
    m.AddPolygonBlock("polygon", {{0, 1, 2, 3}, {3, 2, 4}});
    EXPECT_THROW(refine(m), std::invalid_argument);
}

// --- determinism -------------------------------------------------------------

TEST(Refine, NumberingIsStableAcrossRepeatedRuns) {
    Mesh m = mt::tri_quad_mesh();
    const RefineResult first = refine(m, opts(2));
    for (int run = 0; run < 5; ++run) {
        const RefineResult again = refine(m, opts(2));
        ASSERT_EQ(again.mMesh.NumCellBlocks(), first.mMesh.NumCellBlocks());
        for (std::size_t b = 0; b < first.mMesh.NumCellBlocks(); ++b)
            ASSERT_EQ(rows_of(again.mMesh, b), rows_of(first.mMesh, b))
                << "run " << run << " block " << b;
        mt::expect_points_close(first.mMesh, again.mMesh, 0.0);
    }
}

}  // namespace
