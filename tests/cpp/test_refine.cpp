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
#include <array>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/refine_templates.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/split.hpp"
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

// --- selective refinement ----------------------------------------------------

namespace selective {

using meshioplusplus::CellType;
using meshioplusplus::RefineClosure;
namespace d = meshioplusplus::detail;

RefineOptions cells_opt(std::vector<std::int64_t> cells,
                        RefineClosure closure = RefineClosure::RedGreen) {
    RefineOptions o;
    o.mCells = std::move(cells);
    o.mClosure = closure;
    return o;
}

// The three conformity oracles. They catch different things and none of them is
// redundant. Each is written as a *counting* function so that
// `TheConformityOraclesActuallyFire` below can assert it reports a mesh that is
// genuinely non-conforming -- an oracle nobody tests for firing is exactly how a
// guard ships inert.

// (1) No facet may be shared by more than two cells.
std::size_t count_overshared_facets(const Mesh& rMesh) {
    std::map<std::set<std::int64_t>, int> counts;
    for (const auto cb : rMesh.CellRange()) {
        const CellType type = meshioplusplus::cell_type_from_name(std::string(cb.Type()));
        const std::size_t npc = cb.NodesPerCell();
        const std::int64_t* conn = cb.Conn().As<std::int64_t>();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            if (meshioplusplus::cell_type_dimension(type) == 3) {
                for (const d::CellFaceDef& f : d::cell_faces(type)) {
                    std::set<std::int64_t> key;
                    for (std::size_t i = 0; i < f.mNumCorners; ++i)
                        key.insert(conn[c * npc + f.mNodes[i]]);
                    ++counts[key];
                }
            } else if (meshioplusplus::cell_type_dimension(type) == 2) {
                for (std::size_t i = 0; i < npc; ++i)
                    ++counts[{conn[c * npc + i], conn[c * npc + (i + 1) % npc]}];
            }
        }
    }
    std::size_t bad = 0;
    for (const auto& [facet, count] : counts)
        bad += count > 2 ? 1 : 0;
    return bad;
}

// (2) The one that actually finds a hanging node. A hanging node leaves one big
// facet on one side and two small ones on the other, so (1) still passes with
// every count at 1. Instead: no output node may sit exactly at the midpoint of
// an output cell's edge. New nodes are order-independent corner means, so the
// comparison is exact rather than approximate.
std::size_t count_hanging_nodes(const Mesh& rMesh) {
    std::size_t hanging = 0;
    const std::size_t dim = rMesh.PointDim();
    const double* p = rMesh.Points().As<double>();
    std::map<std::array<double, 3>, std::int64_t> at;
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        std::array<double, 3> key{0, 0, 0};
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            key[k] = p[i * dim + k];
        at.emplace(key, static_cast<std::int64_t>(i));
    }
    for (const auto cb : rMesh.CellRange()) {
        const CellType type = meshioplusplus::cell_type_from_name(std::string(cb.Type()));
        const std::size_t npc = cb.NodesPerCell();
        const std::int64_t* conn = cb.Conn().As<std::int64_t>();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            for (const d::CellEdgePair& e : d::cell_refine_edges(type)) {
                const std::int64_t u = conn[c * npc + e[0]];
                const std::int64_t v = conn[c * npc + e[1]];
                std::array<double, 3> mid{0, 0, 0};
                for (std::size_t k = 0; k < dim && k < 3; ++k)
                    mid[k] = (p[static_cast<std::size_t>(u) * dim + k] +
                              p[static_cast<std::size_t>(v) * dim + k]) *
                             0.5;
                if (at.find(mid) != at.end())
                    ++hanging;
            }
        }
    }
    return hanging;
}

// (3) A conforming refinement does not move the boundary. compute_stats reports
// "area of 2D cells + boundary area of 3D cells", so for a VOLUME mesh an
// unmatched interior facet counts as boundary and the number jumps. On a surface
// mesh the same call only re-checks area conservation, which is worth having but
// is not a conformity test -- (2) is the one that carries that weight in 2D.
void expect_boundary_preserved(const Mesh& rIn, const Mesh& rOut) {
    EXPECT_NEAR(compute_stats(rOut).mTotalArea, compute_stats(rIn).mTotalArea, 1e-12);
}

void expect_conforming(const Mesh& rIn, const Mesh& rOut) {
    EXPECT_EQ(count_overshared_facets(rOut), 0u) << "a facet is shared by more than two cells";
    EXPECT_EQ(count_hanging_nodes(rOut), 0u) << "the mesh has hanging nodes";
    expect_boundary_preserved(rIn, rOut);
}

// An n x n grid of quadrilaterals over the unit square.
Mesh quad_grid(std::size_t n) {
    std::vector<std::vector<double>> pts;
    for (std::size_t j = 0; j <= n; ++j)
        for (std::size_t i = 0; i <= n; ++i)
            pts.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});
    std::vector<std::vector<std::int64_t>> cells;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const std::int64_t a = static_cast<std::int64_t>(j * (n + 1) + i);
            cells.push_back({a, a + 1, a + static_cast<std::int64_t>(n) + 2,
                             a + static_cast<std::int64_t>(n) + 1});
        }
    }
    return mt::make_mesh(std::move(pts), "quad", std::move(cells));
}

// Two triangles per cell of an n x n grid.
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

// An n x n x n grid of hexahedra.
Mesh hex_grid(std::size_t n) {
    const std::size_t s = n + 1;
    std::vector<std::vector<double>> pts;
    for (std::size_t k = 0; k < s; ++k)
        for (std::size_t j = 0; j < s; ++j)
            for (std::size_t i = 0; i < s; ++i)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    const auto id = [s](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<std::int64_t>((k * s + j) * s + i);
    };
    std::vector<std::vector<std::int64_t>> cells;
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                cells.push_back({id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k),
                                 id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1),
                                 id(i, j + 1, k + 1)});
    return mt::make_mesh(std::move(pts), "hexahedron", std::move(cells));
}

// A conformity oracle that cannot fail proves nothing, so build the two failure
// modes by hand and check that each is reported.
TEST(RefineSelective, TheConformityOraclesActuallyFire) {
    // A hanging node: the left quadrilateral is split in two across its mid
    // height, its right-hand neighbour is not, so node 7 at (1, 0.5) sits
    // exactly on the neighbour's edge (1, 4).
    Mesh hanging = mt::make_mesh({{0, 0, 0},
                                  {1, 0, 0},
                                  {2, 0, 0},
                                  {0, 1, 0},
                                  {1, 1, 0},
                                  {2, 1, 0},
                                  {0, 0.5, 0},
                                  {1, 0.5, 0}},
                                 "quad", {{0, 1, 7, 6}, {6, 7, 4, 3}, {1, 2, 5, 4}});
    EXPECT_EQ(count_hanging_nodes(hanging), 1u) << "the hanging-node oracle is inert";
    EXPECT_EQ(count_overshared_facets(hanging), 0u)
        << "and it is the only oracle that sees this failure -- facet counting cannot";

    // Three triangles meeting at one edge: a non-manifold configuration the
    // facet counter is the one that sees.
    Mesh overshared = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}},
                                    "triangle", {{0, 1, 2}, {0, 1, 3}, {0, 1, 4}});
    EXPECT_EQ(count_overshared_facets(overshared), 1u) << "the facet oracle is inert";
}

TEST(RefineSelective, NoSelectorIsTheUniformRefinement) {
    Mesh m = mt::tri_quad_mesh();
    const RefineResult uniform = refine(m, opts(1));
    RefineOptions o;  // every field default: no selector at all
    const RefineResult same = refine(m, o);
    ASSERT_EQ(same.mMesh.NumCellBlocks(), uniform.mMesh.NumCellBlocks());
    for (std::size_t b = 0; b < uniform.mMesh.NumCellBlocks(); ++b)
        EXPECT_EQ(rows_of(same.mMesh, b), rows_of(uniform.mMesh, b));
    mt::expect_points_close(uniform.mMesh, same.mMesh, 0.0);
}

TEST(RefineSelective, ASelectorThatMatchesNothingLeavesTheMeshAlone) {
    Mesh m = mt::tri_quad_mesh();
    RefineOptions o;
    o.mPredicateArray = "";
    o.mCells = {};
    // An empty cell list is "no selector", so use a predicate no cell satisfies.
    Mesh tagged = mt::tri_quad_mesh();
    std::vector<meshioplusplus::NDArray> blocks;
    std::size_t b = 0;
    for (const auto cb : tagged.CellRange()) {
        (void)b;
        NDArray a = NDArray::Uninit(meshioplusplus::DType::Float64, {cb.NumCells(), 1});
        std::fill(a.As<double>(), a.As<double>() + cb.NumCells(), 5.0);
        blocks.push_back(std::move(a));
    }
    tagged.AddCellData("score", std::move(blocks));
    RefineOptions p;
    p.mPredicateArray = "score";
    p.mPredicateOp = meshioplusplus::RefineCompare::Less;
    p.mPredicateValue = 1.0;
    const RefineResult res = refine(tagged, p);
    EXPECT_EQ(res.mMesh.NumPoints(), tagged.NumPoints());
    for (std::size_t i = 0; i < tagged.NumCellBlocks(); ++i)
        EXPECT_EQ(res.mMesh.Cells(i).NumCells(), tagged.Cells(i).NumCells());
}

TEST(RefineSelective, OneTriangleGreensItsNeighboursConformingly) {
    Mesh m = tri_grid(3);
    const std::size_t before = m.Cells(0).NumCells();
    const RefineResult res = refine(m, cells_opt({8}));
    expect_conforming(m, res.mMesh);
    // The selected triangle becomes 4; its three edge-neighbours are split
    // transitionally rather than fully, and nothing else is touched.
    EXPECT_GT(res.mMesh.Cells(0).NumCells(), before);
    EXPECT_LT(res.mMesh.Cells(0).NumCells(), before + 12);
    EXPECT_EQ(res.mCellMaps.size(), 1u);
    const std::int64_t* map = res.mCellMaps[0].As<std::int64_t>();
    EXPECT_EQ(map[8 + 1] - map[8], 4) << "the selected cell must get the full 1-to-4 split";
}

// The test that pins the anisotropic quadrilateral rule. Promoting a single
// split edge to its opposite pair confines the extra refinement to one row, so
// the cost is O(n) rather than the O(n^2) of full propagation.
TEST(RefineSelective, OneQuadRefinesOneRowNotTheWholeGrid) {
    const std::size_t n = 8;
    Mesh m = quad_grid(n);
    const std::size_t before = m.Cells(0).NumCells();  // 64
    const RefineResult res = refine(m, cells_opt({27}));
    expect_conforming(m, res.mMesh);
    EXPECT_EQ(std::string(res.mMesh.Cells(0).Type()), "quad") << "green quads stay quads";
    // The selected cell contributes 4; its row and column each contribute one
    // extra cell per member. Full propagation would give 4 * 64 = 256.
    EXPECT_LT(res.mMesh.Cells(0).NumCells(), before + 6 * n);
    EXPECT_GT(res.mMesh.Cells(0).NumCells(), before);
}

TEST(RefineSelective, OneHexRefinesSheetsNotTheWholeBlock) {
    const std::size_t n = 4;
    Mesh m = hex_grid(n);
    const std::size_t before = m.Cells(0).NumCells();  // 64
    const RefineResult res = refine(m, cells_opt({21}));
    expect_conforming(m, res.mMesh);
    EXPECT_EQ(std::string(res.mMesh.Cells(0).Type()), "hexahedron");
    // Three sheets of n*n cells split in two, plus the selected cell's own
    // 8-way split. Full propagation would give 8 * 64 = 512.
    EXPECT_LT(res.mMesh.Cells(0).NumCells(), before + 4 * n * n);
}

TEST(RefineSelective, OneTetraGreensItsNeighboursConformingly) {
    Mesh m = mt::tet_mesh();
    const RefineResult res = refine(m, cells_opt({0}));
    expect_conforming(m, res.mMesh);
    EXPECT_EQ(std::string(res.mMesh.Cells(0).Type()), "tetra");
}

// Propagation is the oracle: on a connected mesh it must reproduce the uniform
// refinement exactly, which is both the correctness check and the honest
// statement of what that mode costs.
TEST(RefineSelective, PropagateReproducesUniformRefinement) {
    // Built one at a time rather than in a range-for over a braced list: the
    // KRATOS Mesh is deliberately not copy-constructible.
    const auto check = [](Mesh m) {
        const RefineResult uniform = refine(m, opts(1));
        const RefineResult propagated = refine(m, cells_opt({0}, RefineClosure::Propagate));
        ASSERT_EQ(propagated.mMesh.NumCellBlocks(), uniform.mMesh.NumCellBlocks());
        for (std::size_t b = 0; b < uniform.mMesh.NumCellBlocks(); ++b)
            EXPECT_EQ(rows_of(propagated.mMesh, b), rows_of(uniform.mMesh, b));
        mt::expect_points_close(uniform.mMesh, propagated.mMesh, 0.0);
    };
    check(tri_grid(3));
    check(quad_grid(3));
    check(hex_grid(2));
}

// --- the balanced closure: 2:1 balance, hanging nodes kept -------------------

// The largest level difference between two cells sharing a NODE. Node adjacency,
// not edge adjacency: across a hanging interface the coarse cell spans a whole
// edge while the fine cell has only half of it, so an edge-keyed check is blind
// to exactly the coarse/fine adjacency 2:1 balance exists to police.
std::int64_t max_level_gap(const Mesh& rMesh) {
    const NDArray& levels = rMesh.CellData("refine:level", 0);
    std::map<std::int64_t, std::pair<std::int64_t, std::int64_t>> range;
    std::size_t c0 = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t npc = cb.NodesPerCell();
        const std::int64_t* conn = cb.Conn().As<std::int64_t>();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const std::int64_t lv = meshioplusplus::detail::read_int(levels, c0 + c);
            for (std::size_t n = 0; n < npc; ++n) {
                auto it = range.find(conn[c * npc + n]);
                if (it == range.end())
                    range.emplace(conn[c * npc + n], std::make_pair(lv, lv));
                else {
                    it->second.first = std::min(it->second.first, lv);
                    it->second.second = std::max(it->second.second, lv);
                }
            }
        }
        c0 += cb.NumCells();
    }
    std::int64_t gap = 0;
    for (const auto& [node, lohi] : range)
        gap = std::max(gap, lohi.second - lohi.first);
    return gap;
}

TEST(RefineSelective, BalancedKeepsHangingNodesInsteadOfClosing) {
    Mesh m = hex_grid(4);
    RefineOptions o = cells_opt({21}, RefineClosure::Balanced);
    o.mRecordLevels = true;
    const RefineResult res = refine(m, o);
    // One cell split into 8, every other cell untouched: 64 - 1 + 8.
    EXPECT_EQ(res.mMesh.Cells(0).NumCells(), 71u);
    // Deliberately NOT conforming -- and it says so.
    EXPECT_GT(count_hanging_nodes(res.mMesh), 0u);
    EXPECT_TRUE(res.mMesh.HasPointData("refine:hanging"));
    EXPECT_LE(max_level_gap(res.mMesh), 1);
}

TEST(RefineSelective, TheConformingClosuresLeaveNoHangingArray) {
    Mesh m = hex_grid(3);
    for (RefineClosure closure : {RefineClosure::RedGreen, RefineClosure::Propagate}) {
        RefineOptions o = cells_opt({13}, closure);
        o.mRecordLevels = true;
        const RefineResult res = refine(m, o);
        EXPECT_EQ(count_hanging_nodes(res.mMesh), 0u);
        EXPECT_FALSE(res.mMesh.HasPointData("refine:hanging"));
    }
}

// --- multi-pass balanced: tearing, and reporting every constrained node -----
//
// Two more counting oracles, for the two things a *second* balanced pass can get
// wrong. Both are computed from the output's geometry and connectivity alone, so
// neither shares a code path with refine's own bookkeeping, and
// `TheTearAndReportOraclesActuallyFire` asserts each of them fires.

// (4) Distinct node ids at the same position, both referenced by cells: the mesh
// is torn, not merely 1-irregular. New nodes are order-independent corner means,
// so equality here is exact rather than approximate.
std::size_t count_torn_positions(const Mesh& rMesh) {
    const std::size_t dim = rMesh.PointDim();
    const double* p = rMesh.Points().As<double>();
    std::set<std::int64_t> used;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t* conn = cb.Conn().As<std::int64_t>();
        for (std::size_t i = 0; i < cb.NumCells() * cb.NodesPerCell(); ++i)
            used.insert(conn[i]);
    }
    std::map<std::array<double, 3>, std::size_t> counts;
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        if (used.find(static_cast<std::int64_t>(i)) == used.end())
            continue;
        std::array<double, 3> key{0, 0, 0};
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            key[k] = p[i * dim + k];
        ++counts[key];
    }
    std::size_t torn = 0;
    for (const auto& [pos, n] : counts)
        torn += n > 1 ? 1 : 0;
    return torn;
}

// (5) Constrained nodes the `refine:hanging` array fails to name. A node is
// constrained when it sits at an edge midpoint or quad-face centre of a cell
// that does not reference it. Both entity kinds matter: an edge-only oracle
// silently passes a mesh whose face centres are unreported.
std::size_t count_unreported_hanging(const Mesh& rMesh) {
    const std::size_t dim = rMesh.PointDim();
    const double* p = rMesh.Points().As<double>();
    const auto coord = [&](std::int64_t node, std::size_t k) {
        return p[static_cast<std::size_t>(node) * dim + k];
    };
    std::map<std::array<double, 3>, std::int64_t> at;
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        std::array<double, 3> key{0, 0, 0};
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            key[k] = p[i * dim + k];
        at.emplace(key, static_cast<std::int64_t>(i));
    }

    std::set<std::int64_t> constrained;
    for (const auto cb : rMesh.CellRange()) {
        const CellType type = meshioplusplus::cell_type_from_name(std::string(cb.Type()));
        const std::size_t npc = cb.NodesPerCell();
        const std::int64_t* conn = cb.Conn().As<std::int64_t>();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const std::int64_t* row = conn + c * npc;
            std::set<std::int64_t> own(row, row + npc);
            std::vector<std::array<double, 3>> spots;
            for (const d::CellEdgePair& e : d::cell_refine_edges(type)) {
                std::array<double, 3> mid{0, 0, 0};
                for (std::size_t k = 0; k < dim && k < 3; ++k)
                    mid[k] = (coord(row[e[0]], k) + coord(row[e[1]], k)) * 0.5;
                spots.push_back(mid);
            }
            for (const d::CellQuadFace& f : d::cell_refine_quad_faces(type)) {
                std::array<double, 3> ctr{0, 0, 0};
                for (std::size_t k = 0; k < dim && k < 3; ++k)
                    ctr[k] = (coord(row[f[0]], k) + coord(row[f[1]], k) + coord(row[f[2]], k) +
                              coord(row[f[3]], k)) *
                             0.25;
                spots.push_back(ctr);
            }
            for (const std::array<double, 3>& s : spots) {
                auto it = at.find(s);
                if (it != at.end() && own.find(it->second) == own.end())
                    constrained.insert(it->second);
            }
        }
    }

    std::set<std::int64_t> reported;
    if (rMesh.HasPointData("refine:hanging")) {
        const NDArray& a = rMesh.PointData("refine:hanging");
        for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
            if (meshioplusplus::detail::read_int(a, i) != 0)
                reported.insert(static_cast<std::int64_t>(i));
    }
    std::size_t missed = 0;
    for (std::int64_t node : constrained)
        missed += reported.find(node) == reported.end() ? 1 : 0;
    // A node reported but not constrained is just as wrong, and would otherwise
    // let a "mark everything" implementation pass.
    for (std::int64_t node : reported)
        missed += constrained.find(node) == constrained.end() ? 1 : 0;
    return missed;
}

TEST(RefineSelective, TheTearAndReportOraclesActuallyFire) {
    // (4) fires on a mesh carrying two referenced nodes at one position.
    {
        Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 0}}, "triangle",
                               {{0, 1, 2}, {0, 2, 3}, {4, 2, 3}});
        EXPECT_EQ(count_torn_positions(m), 1u) << "the tear oracle does not fire";
        Mesh clean = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle",
                                   {{0, 1, 2}, {0, 2, 3}});
        EXPECT_EQ(count_torn_positions(clean), 0u);
    }
    // (5) fires on a mesh with a genuine hanging node and no array to name it:
    // the big triangle spans the edge whose midpoint node 4 occupies.
    {
        Mesh m = mt::make_mesh({{0, 0, 0}, {2, 0, 0}, {1, 2, 0}, {2, 2, 0}, {1, 0, 0}}, "triangle",
                               {{0, 1, 2}, {1, 3, 2}, {0, 4, 2}});
        EXPECT_GT(count_unreported_hanging(m), 0u) << "the report oracle does not fire";
    }
}

TEST(RefineSelective, MultiPassBalancedDoesNotTearTheMesh) {
    // Refining one cell twice draws its coarse neighbours in by the balance rule.
    // Each of those already carries hanging nodes from the first pass, and their
    // own refinement must reuse them rather than allocate a coincident second
    // node -- which is what refine:entity is for.
    for (int levels : {2, 3}) {
        Mesh m = hex_grid(4);
        RefineOptions o = cells_opt({0}, RefineClosure::Balanced);
        o.mLevels = levels;
        o.mRecordLevels = true;
        const RefineResult res = refine(m, o);
        EXPECT_EQ(count_torn_positions(res.mMesh), 0u)
            << "levels=" << levels << ": the mesh is torn";
        EXPECT_LE(max_level_gap(res.mMesh), 1) << "levels=" << levels;
    }
}

TEST(RefineSelective, MultiPassBalancedReportsEveryConstrainedNode) {
    // The array must name exactly the constrained nodes -- neither a superset nor
    // a subset, which is what doc/refine.md promises. Stating the rule over the
    // *input* cells' entities under-reports here: a cell the balance rule draws
    // in has children whose sub-edges the input cell never had.
    for (int levels : {1, 2, 3}) {
        Mesh m = hex_grid(4);
        RefineOptions o = cells_opt({0}, RefineClosure::Balanced);
        o.mLevels = levels;
        o.mRecordLevels = true;
        const RefineResult res = refine(m, o);
        EXPECT_EQ(count_unreported_hanging(res.mMesh), 0u) << "levels=" << levels;
    }
}

TEST(RefineSelective, EntityKeysDescribeTheNodesTheyName) {
    Mesh m = hex_grid(3);
    RefineOptions o = cells_opt({13}, RefineClosure::Balanced);
    const RefineResult res = refine(m, o);
    ASSERT_TRUE(res.mMesh.HasPointData("refine:entity"));
    const NDArray& keys = res.mMesh.PointData("refine:entity");
    const std::size_t dim = res.mMesh.PointDim();
    const double* p = res.mMesh.Points().As<double>();
    std::size_t keyed = 0;
    for (std::size_t i = 0; i < res.mMesh.NumPoints(); ++i) {
        std::array<std::int64_t, 4> key{};
        for (std::size_t k = 0; k < 4; ++k)
            key[k] = meshioplusplus::detail::read_int(keys, i * 4 + k);
        if (key[3] < 0)
            continue;  // sentinel: an original point or a body centre
        ++keyed;
        const std::size_t first = key[0] < 0 ? 2 : 0;
        for (std::size_t k = 0; k < dim; ++k) {
            double sum = 0.0;
            for (std::size_t c = first; c < 4; ++c)
                sum += p[static_cast<std::size_t>(key[c]) * dim + k];
            EXPECT_DOUBLE_EQ(p[i * dim + k], sum / static_cast<double>(4 - first));
        }
    }
    EXPECT_GT(keyed, 0u) << "no point carries a real entity key";
}

TEST(RefineSelective, TheConformingClosuresAttachNoEntityArray) {
    // The array is what keeps the conforming closures' output byte-identical to
    // what it was before this bookkeeping existed: they cannot tear, so they must
    // not pay for it.
    Mesh m = hex_grid(3);
    for (RefineClosure closure : {RefineClosure::RedGreen, RefineClosure::Propagate}) {
        RefineOptions o = cells_opt({13}, closure);
        o.mLevels = 2;
        EXPECT_FALSE(refine(m, o).mMesh.HasPointData("refine:entity"));
    }
    RefineOptions uniform;  // no selector at all
    uniform.mLevels = 2;
    EXPECT_FALSE(refine(m, uniform).mMesh.HasPointData("refine:entity"));
}

TEST(RefineSelective, AStaleEntityArrayIsIgnoredRatherThanTrusted) {
    // The keys name point indices, which refine never renumbers but other
    // operations do. Trusting a stale one would graft an entity's node onto
    // another's, so a key that no longer reproduces its own point's coordinates
    // invalidates the whole array -- refine falls back to allocating fresh nodes,
    // which is merely the old behaviour rather than a wrong answer.
    Mesh stale = hex_grid(3);
    NDArray bad =
        NDArray::Uninit(meshioplusplus::DType::Int64, {stale.NumPoints(), std::size_t{4}});
    std::int64_t* dst = bad.As<std::int64_t>();
    std::fill(dst, dst + stale.NumPoints() * 4, static_cast<std::int64_t>(-1));
    // Point 0 is an original corner; claiming it is the midpoint of edge (1, 2)
    // is exactly the shape a renumbering leaves behind.
    dst[2] = 1;
    dst[3] = 2;
    stale.AddPointData("refine:entity", std::move(bad));

    RefineOptions o = cells_opt({13}, RefineClosure::Balanced);
    const RefineResult res = refine(stale, o);
    const RefineResult reference = refine(hex_grid(3), o);

    // Dropped, not obeyed: the output is exactly what a mesh carrying no array at
    // all produces.
    ASSERT_EQ(res.mMesh.NumPoints(), reference.mMesh.NumPoints());
    ASSERT_EQ(res.mMesh.Cells(0).NumCells(), reference.mMesh.Cells(0).NumCells());
    const std::int64_t* got = res.mMesh.Cells(0).Conn().As<std::int64_t>();
    const std::int64_t* want = reference.mMesh.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t i = 0; i < reference.mMesh.Cells(0).NumCells() * 8; ++i)
        ASSERT_EQ(got[i], want[i]) << "at " << i;
    EXPECT_EQ(count_unreported_hanging(res.mMesh), 0u);
}

TEST(RefineSelective, BalancedDoesNotPropagateOnAMeshOfUniformLevel) {
    // The whole point: with every cell at the same level, refining one puts
    // nothing else out of balance, so nothing else is touched.
    const auto check = [](Mesh m, std::int64_t selected, std::size_t children) {
        const std::size_t before = m.Cells(0).NumCells();
        RefineOptions o = cells_opt({selected}, RefineClosure::Balanced);
        o.mRecordLevels = true;
        const RefineResult res = refine(m, o);
        EXPECT_EQ(res.mMesh.Cells(0).NumCells(), before - 1 + children);
        EXPECT_LE(max_level_gap(res.mMesh), 1);
    };
    check(tri_grid(4), 8, 4);
    check(quad_grid(5), 7, 4);
    check(hex_grid(4), 21, 8);
}

TEST(RefineSelective, BalancedDrawsInNeighboursThatWouldFallTwoLevelsBehind) {
    Mesh m = hex_grid(4);
    RefineOptions first_opts = cells_opt({21}, RefineClosure::Balanced);
    first_opts.mRecordLevels = true;
    const RefineResult first = refine(m, first_opts);

    // Refining a level-1 child would leave its level-0 neighbours two levels
    // coarser, so balancing must promote them -- and only them.
    const NDArray& levels = first.mMesh.CellData("refine:level", 0);
    std::int64_t child = -1;
    for (std::size_t i = 0; i < levels.Size() && child < 0; ++i)
        if (meshioplusplus::detail::read_int(levels, i) == 1)
            child = static_cast<std::int64_t>(i);
    ASSERT_GE(child, 0);

    RefineOptions second_opts = cells_opt({child}, RefineClosure::Balanced);
    second_opts.mRecordLevels = true;
    const RefineResult second = refine(first.mMesh, second_opts);
    const std::size_t grew = second.mMesh.Cells(0).NumCells() - first.mMesh.Cells(0).NumCells();
    EXPECT_GT(grew, 7u) << "balancing must have drawn in cells beyond the selected one";
    EXPECT_LT(grew, 7u * 30u) << "and must stay local, not reach the whole mesh";
    EXPECT_LE(max_level_gap(second.mMesh), 1);
}

TEST(RefineSelective, BalancedIsFarCheaperThanPropagate) {
    Mesh m = hex_grid(3);
    const RefineResult balanced = refine(m, cells_opt({13}, RefineClosure::Balanced));
    const RefineResult propagated = refine(m, cells_opt({13}, RefineClosure::Propagate));
    EXPECT_LT(balanced.mMesh.Cells(0).NumCells() * 4, propagated.mMesh.Cells(0).NumCells());
}

TEST(RefineSelective, TwoSuccessivePassesStayConforming) {
    Mesh m = tri_grid(4);
    const RefineResult first = refine(m, cells_opt({10}));
    expect_conforming(m, first.mMesh);
    const RefineResult second = refine(first.mMesh, cells_opt({0}));
    expect_conforming(first.mMesh, second.mMesh);
}

TEST(RefineSelective, MultipleLevelsRefineTheChildrenOfTheRedCells) {
    Mesh m = tri_grid(3);
    RefineOptions o = cells_opt({4});
    o.mLevels = 2;
    o.mRecordLevels = true;
    const RefineResult res = refine(m, o);
    expect_conforming(m, res.mMesh);
    ASSERT_TRUE(res.mMesh.HasCellData("refine:level"));
    const NDArray& lv = res.mMesh.CellData("refine:level", 0);
    std::int64_t max_level = 0;
    for (std::size_t i = 0; i < lv.Size(); ++i)
        max_level = std::max(max_level, meshioplusplus::detail::read_int(lv, i));
    EXPECT_EQ(max_level, 2) << "the selected cell's descendants must reach level 2";
}

// --- refine:level ------------------------------------------------------------

TEST(RefineSelective, LevelsAreNotRecordedUnlessAsked) {
    Mesh m = tri_grid(2);
    EXPECT_FALSE(refine(m, cells_opt({0})).mMesh.HasCellData("refine:level"));
}

TEST(RefineSelective, GreenChildrenInheritTheirParentsLevel) {
    Mesh m = tri_grid(3);
    RefineOptions o = cells_opt({8});
    o.mRecordLevels = true;
    const RefineResult res = refine(m, o);
    ASSERT_TRUE(res.mMesh.HasCellData("refine:level"));
    const NDArray& lv = res.mMesh.CellData("refine:level", 0);
    const std::int64_t* map = res.mCellMaps[0].As<std::int64_t>();
    // The selected cell's four children are level 1 ...
    for (std::int64_t c = map[8]; c < map[9]; ++c)
        EXPECT_EQ(meshioplusplus::detail::read_int(lv, static_cast<std::size_t>(c)), 1);
    // ... while a green neighbour's children stay at level 0: a transitional
    // split is a closure, not a refinement.
    std::size_t green_children = 0;
    for (std::size_t p = 0; p < res.mCellMaps[0].Size(); ++p) {
        if (p == 8)
            continue;
        const std::int64_t lo = map[p];
        const std::int64_t hi = p + 1 < res.mCellMaps[0].Size()
                                    ? map[p + 1]
                                    : static_cast<std::int64_t>(res.mMesh.Cells(0).NumCells());
        if (hi - lo <= 1)
            continue;
        ++green_children;
        for (std::int64_t c = lo; c < hi; ++c)
            EXPECT_EQ(meshioplusplus::detail::read_int(lv, static_cast<std::size_t>(c)), 0);
    }
    EXPECT_GT(green_children, 0u) << "the selection must have produced green cells to check";
}

TEST(RefineSelective, AnExistingLevelArrayAccumulatesRatherThanReplicating) {
    Mesh m = tri_grid(3);
    RefineOptions o = cells_opt({4});
    o.mRecordLevels = true;
    const RefineResult first = refine(m, o);
    // The second pass sees the array and must UPDATE it, not replicate it -- the
    // flag only controls creating it.
    RefineOptions second = cells_opt({0});
    const RefineResult res = refine(first.mMesh, second);
    ASSERT_TRUE(res.mMesh.HasCellData("refine:level"));
    const NDArray& lv = res.mMesh.CellData("refine:level", 0);
    std::int64_t max_level = 0;
    for (std::size_t i = 0; i < lv.Size(); ++i)
        max_level = std::max(max_level, meshioplusplus::detail::read_int(lv, i));
    EXPECT_GE(max_level, 1);
}

// --- selectors ---------------------------------------------------------------

TEST(RefineSelective, CellRegionSelectsTheSameCellsAsAnIndexList) {
    Mesh a = tri_grid(3);
    Mesh b = tri_grid(3);
    NDArray entries = NDArray::Uninit(meshioplusplus::DType::Int64, {2});
    entries.As<std::int64_t>()[0] = 4;
    entries.As<std::int64_t>()[1] = 9;
    b.AddRegion(
        meshioplusplus::Region("hot", meshioplusplus::RegionKind::Cell, std::move(entries)));

    const RefineResult by_index = refine(a, cells_opt({4, 9}));
    RefineOptions o;
    o.mRegion = "hot";
    const RefineResult by_region = refine(b, o);
    EXPECT_EQ(rows_of(by_index.mMesh, 0), rows_of(by_region.mMesh, 0));
}

TEST(RefineSelective, PointRegionSelectsEveryCellTouchingIt) {
    Mesh m = tri_grid(3);
    NDArray entries = NDArray::Uninit(meshioplusplus::DType::Int64, {1});
    entries.As<std::int64_t>()[0] = 0;  // the corner node
    m.AddRegion(
        meshioplusplus::Region("corner", meshioplusplus::RegionKind::Point, std::move(entries)));
    RefineOptions o;
    o.mRegion = "corner";
    const RefineResult res = refine(m, o);
    expect_conforming(m, res.mMesh);
    EXPECT_GT(res.mMesh.Cells(0).NumCells(), m.Cells(0).NumCells());
}

TEST(RefineSelective, PredicateSelectsByCellData) {
    Mesh m = tri_grid(3);
    const std::size_t n = m.Cells(0).NumCells();
    NDArray score = NDArray::Uninit(meshioplusplus::DType::Float64, {n, 1});
    for (std::size_t c = 0; c < n; ++c)
        score.As<double>()[c] = c == 5 ? 0.1 : 0.9;
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(score));
    m.AddCellData("quality:scaled_jacobian", std::move(blocks));

    RefineOptions o;
    o.mPredicateArray = "quality:scaled_jacobian";
    o.mPredicateOp = meshioplusplus::RefineCompare::Less;
    o.mPredicateValue = 0.3;
    const RefineResult by_pred = refine(m, o);
    const RefineResult by_index = refine(m, cells_opt({5}));
    EXPECT_EQ(rows_of(by_pred.mMesh, 0), rows_of(by_index.mMesh, 0));
}

// compute_quality reports NaN where a metric does not apply, so a predicate over
// `quality:*` on a mixed mesh is the headline use case. Rejecting the array
// would break it; a non-finite value simply never matches.
TEST(RefineSelective, NonFiniteCellValuesNeverMatchAndAreNotAnError) {
    Mesh m = tri_grid(2);
    const std::size_t n = m.Cells(0).NumCells();
    NDArray score = NDArray::Uninit(meshioplusplus::DType::Float64, {n, 1});
    for (std::size_t c = 0; c < n; ++c)
        score.As<double>()[c] = std::numeric_limits<double>::quiet_NaN();
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(score));
    m.AddCellData("q", std::move(blocks));
    RefineOptions o;
    o.mPredicateArray = "q";
    o.mPredicateOp = meshioplusplus::RefineCompare::Less;
    o.mPredicateValue = 1.0;
    const RefineResult res = refine(m, o);
    EXPECT_EQ(res.mMesh.Cells(0).NumCells(), n) << "no cell should have been selected";
}

// --- validation --------------------------------------------------------------

TEST(RefineSelective, TwoSelectorsRaise) {
    Mesh m = tri_grid(2);
    RefineOptions o;
    o.mCells = {0};
    o.mRegion = "anything";
    EXPECT_THROW(refine(m, o), std::invalid_argument);
}

TEST(RefineSelective, AnOutOfRangeCellIndexRaises) {
    Mesh m = tri_grid(2);
    EXPECT_THROW(refine(m, cells_opt({1000})), std::invalid_argument);
    EXPECT_THROW(refine(m, cells_opt({-1})), std::invalid_argument);
}

TEST(RefineSelective, AnUnknownRegionRaises) {
    Mesh m = tri_grid(2);
    RefineOptions o;
    o.mRegion = "nope";
    EXPECT_THROW(refine(m, o), std::invalid_argument);
}

TEST(RefineSelective, ASideRegionIsNotASelector) {
    Mesh m = mt::tet_mesh();
    NDArray entries = NDArray::Uninit(meshioplusplus::DType::Int64, {1, 2});
    entries.As<std::int64_t>()[0] = 0;
    entries.As<std::int64_t>()[1] = 0;
    m.AddRegion(
        meshioplusplus::Region("wall", meshioplusplus::RegionKind::Side, std::move(entries)));
    RefineOptions o;
    o.mRegion = "wall";
    EXPECT_THROW(refine(m, o), std::invalid_argument);
}

TEST(RefineSelective, AMissingOrVectorPredicateArrayRaises) {
    Mesh m = tri_grid(2);
    RefineOptions o;
    o.mPredicateArray = "absent";
    EXPECT_THROW(refine(m, o), std::invalid_argument);

    const std::size_t n = m.Cells(0).NumCells();
    NDArray vec = NDArray::Uninit(meshioplusplus::DType::Float64, {n, 3});
    std::fill(vec.As<double>(), vec.As<double>() + vec.Size(), 1.0);
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(vec));
    m.AddCellData("v", std::move(blocks));
    RefineOptions p;
    p.mPredicateArray = "v";
    EXPECT_THROW(refine(m, p), std::invalid_argument);
}

TEST(RefineSelective, ClosureAndComparisonNamesParse) {
    using meshioplusplus::refine_closure_from_name;
    using meshioplusplus::refine_compare_from_name;
    EXPECT_EQ(refine_closure_from_name(""), RefineClosure::RedGreen);
    EXPECT_EQ(refine_closure_from_name("green"), RefineClosure::RedGreen);
    EXPECT_EQ(refine_closure_from_name("redgreen"), RefineClosure::RedGreen);
    EXPECT_EQ(refine_closure_from_name("propagate"), RefineClosure::Propagate);
    EXPECT_EQ(refine_closure_from_name("balanced"), RefineClosure::Balanced);
    EXPECT_EQ(refine_closure_from_name("2:1"), RefineClosure::Balanced);
    EXPECT_THROW(refine_closure_from_name("blue"), std::invalid_argument);
    EXPECT_EQ(refine_compare_from_name("<"), meshioplusplus::RefineCompare::Less);
    EXPECT_EQ(refine_compare_from_name(">="), meshioplusplus::RefineCompare::GreaterEqual);
    EXPECT_EQ(refine_compare_from_name("!="), meshioplusplus::RefineCompare::NotEqual);
    EXPECT_THROW(refine_compare_from_name("~="), std::invalid_argument);
}

// --- data and regions through a selective pass -------------------------------

TEST(RefineSelective, CellDataIsGatheredNotBlindlyRepeated) {
    Mesh m = tri_grid(2);
    const std::size_t n = m.Cells(0).NumCells();
    NDArray tag = NDArray::Uninit(meshioplusplus::DType::Int64, {n, 1});
    for (std::size_t c = 0; c < n; ++c)
        tag.As<std::int64_t>()[c] = static_cast<std::int64_t>(c);
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(tag));
    m.AddCellData("tag", std::move(blocks));

    const RefineResult res = refine(m, cells_opt({3}));
    const NDArray& out = res.mMesh.CellData("tag", 0);
    const std::int64_t* map = res.mCellMaps[0].As<std::int64_t>();
    for (std::size_t p = 0; p < n; ++p) {
        const std::int64_t lo = map[p];
        const std::int64_t hi =
            p + 1 < n ? map[p + 1] : static_cast<std::int64_t>(res.mMesh.Cells(0).NumCells());
        for (std::int64_t c = lo; c < hi; ++c)
            EXPECT_EQ(meshioplusplus::detail::read_int(out, static_cast<std::size_t>(c)),
                      static_cast<std::int64_t>(p));
    }
}

TEST(RefineSelective, PointDataStaysExactForALinearField) {
    Mesh m = tri_grid(3);
    NDArray f = NDArray::Uninit(meshioplusplus::DType::Float64, {m.NumPoints(), 1});
    const double* p = m.Points().As<double>();
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        f.As<double>()[i] = 3.0 * p[i * 3] + 5.0 * p[i * 3 + 1] + 1.0;
    m.AddPointData("f", std::move(f));

    const RefineResult res = refine(m, cells_opt({5}));
    const NDArray& out = res.mMesh.PointData("f");
    const double* q = res.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < res.mMesh.NumPoints(); ++i)
        EXPECT_NEAR(out.As<double>()[i], 3.0 * q[i * 3] + 5.0 * q[i * 3 + 1] + 1.0, 1e-12);
}

TEST(RefineSelective, RegionsExpandToTheChildren) {
    Mesh m = tri_grid(2);
    NDArray entries = NDArray::Uninit(meshioplusplus::DType::Int64, {1});
    entries.As<std::int64_t>()[0] = 3;
    m.AddRegion(
        meshioplusplus::Region("patch", meshioplusplus::RegionKind::Cell, std::move(entries)));
    const RefineResult res = refine(m, cells_opt({3}));
    const std::size_t idx = res.mMesh.FindRegion("patch", meshioplusplus::RegionKind::Cell);
    ASSERT_NE(idx, Mesh::npos);
    EXPECT_EQ(res.mMesh.Region(idx).NumEntries(), 4u) << "the refined cell's four children";
}

TEST(RefineSelective, RecordParentIdsStillNamesTheOriginalAncestor) {
    Mesh m = tri_grid(2);
    RefineOptions o = cells_opt({3});
    o.mRecordParentIds = true;
    const RefineResult res = refine(m, o);
    ASSERT_TRUE(res.mMesh.HasCellData("refine:parent_cell"));
    const NDArray& parents = res.mMesh.CellData("refine:parent_cell", 0);
    const std::int64_t* map = res.mCellMaps[0].As<std::int64_t>();
    for (std::int64_t c = map[3]; c < map[4]; ++c)
        EXPECT_EQ(meshioplusplus::detail::read_int(parents, static_cast<std::size_t>(c)), 3);
}

TEST(RefineSelective, SelectiveNumberingIsStableAcrossRepeatedRuns) {
    Mesh m = quad_grid(5);
    const RefineResult first = refine(m, cells_opt({7, 18}));
    for (int run = 0; run < 5; ++run) {
        const RefineResult again = refine(m, cells_opt({7, 18}));
        ASSERT_EQ(rows_of(again.mMesh, 0), rows_of(first.mMesh, 0)) << "run " << run;
        mt::expect_points_close(first.mMesh, again.mMesh, 0.0);
    }
}

// --- refine:cell_id / refine:parent_id ----------------------------------------
//
// The persistent hierarchy: "a link between two meshes, not a tree inside
// one". A multigrid or green-undo caller keeps every pass's OUTPUT mesh; the
// fine mesh's refine:parent_id resolves against the coarse mesh's
// refine:cell_id (implicit, or the coarse mesh's own recorded ids).

TEST(RefineSelective, HierarchyIsNotRecordedUnlessAsked) {
    Mesh m = tri_grid(2);
    const RefineResult res = refine(m, cells_opt({0}));
    EXPECT_FALSE(res.mMesh.HasCellData("refine:cell_id"));
    EXPECT_FALSE(res.mMesh.HasCellData("refine:parent_id"));
}

TEST(RefineSelective, AnUnsplitCellKeepsItsIdAndIsItsOwnParent) {
    Mesh m = tri_grid(3);
    RefineOptions o = cells_opt({8});
    o.mRecordHierarchy = true;
    const RefineResult res = refine(m, o);
    ASSERT_TRUE(res.mMesh.HasCellData("refine:cell_id"));
    ASSERT_TRUE(res.mMesh.HasCellData("refine:parent_id"));
    const NDArray& ids = res.mMesh.CellData("refine:cell_id", 0);
    const NDArray& parents = res.mMesh.CellData("refine:parent_id", 0);
    const std::int64_t* map = res.mCellMaps[0].As<std::int64_t>();
    const std::size_t nparents = res.mCellMaps[0].Size();
    const std::size_t ncells = res.mMesh.Cells(0).NumCells();
    std::size_t untouched_checked = 0;
    for (std::size_t p = 0; p < nparents; ++p) {
        const std::int64_t lo = map[p];
        const std::int64_t hi = p + 1 < nparents ? map[p + 1] : static_cast<std::int64_t>(ncells);
        if (hi - lo != 1)
            continue;  // this parent was split; not what this test checks
        const std::size_t c = static_cast<std::size_t>(lo);
        EXPECT_EQ(meshioplusplus::detail::read_int(ids, c), static_cast<std::int64_t>(p))
            << "an untouched cell must keep its own implicit id";
        EXPECT_EQ(meshioplusplus::detail::read_int(parents, c), static_cast<std::int64_t>(p))
            << "an untouched cell is its own parent";
        ++untouched_checked;
    }
    EXPECT_GT(untouched_checked, 0u) << "the selection must leave some cell untouched to check";
}

TEST(RefineSelective, EveryCellIdIsUniqueAcrossBlocks) {
    Mesh m = mt::tri_quad_mesh();
    RefineOptions o;
    o.mLevels = 2;
    o.mRecordHierarchy = true;
    const RefineResult res = refine(m, o);
    std::set<std::int64_t> seen;
    for (std::size_t b = 0; b < res.mMesh.NumCellBlocks(); ++b) {
        const NDArray& ids = res.mMesh.CellData("refine:cell_id", b);
        for (std::size_t i = 0; i < ids.Size(); ++i) {
            const std::int64_t id = meshioplusplus::detail::read_int(ids, i);
            EXPECT_TRUE(seen.insert(id).second) << "duplicate id " << id << " in block " << b;
        }
    }
}

TEST(RefineSelective, ImplicitIdsAreTheGlobalBlockMajorIndex) {
    // A mixed-block mesh, so implicit numbering across blocks (not just within
    // one) is what is actually being checked -- the same numbering
    // detail::block_bases/partition_labels use.
    Mesh m = mt::tri_quad_mesh();
    RefineOptions o = cells_opt({0});
    o.mRecordHierarchy = true;
    const RefineResult res = refine(m, o);
    const std::vector<std::int64_t> bases = meshioplusplus::detail::block_bases(m);
    for (std::size_t b = 0; b < res.mMesh.NumCellBlocks(); ++b) {
        const NDArray& ids = res.mMesh.CellData("refine:cell_id", b);
        const std::int64_t* map = res.mCellMaps[b].As<std::int64_t>();
        const std::size_t nparents = res.mCellMaps[b].Size();
        const std::size_t ncells = res.mMesh.Cells(b).NumCells();
        for (std::size_t p = 0; p < nparents; ++p) {
            const std::int64_t lo = map[p];
            const std::int64_t hi =
                p + 1 < nparents ? map[p + 1] : static_cast<std::int64_t>(ncells);
            if (hi - lo != 1)
                continue;  // only untouched cells keep their implicit id
            EXPECT_EQ(meshioplusplus::detail::read_int(ids, static_cast<std::size_t>(lo)),
                      bases[b] + static_cast<std::int64_t>(p));
        }
    }
}

TEST(RefineSelective, FreshIdsExceedEveryInputId) {
    Mesh m = tri_grid(3);
    RefineOptions o1 = cells_opt({8});
    o1.mRecordHierarchy = true;
    const RefineResult first = refine(m, o1);
    std::int64_t max_first_id = -1;
    const NDArray& ids1 = first.mMesh.CellData("refine:cell_id", 0);
    for (std::size_t i = 0; i < ids1.Size(); ++i)
        max_first_id = std::max(max_first_id, meshioplusplus::detail::read_int(ids1, i));

    RefineOptions o2 = cells_opt({0});
    o2.mRecordHierarchy = true;
    const RefineResult second = refine(first.mMesh, o2);
    const NDArray& ids2 = second.mMesh.CellData("refine:cell_id", 0);
    const std::int64_t* map2 = second.mCellMaps[0].As<std::int64_t>();
    const std::size_t nparents2 = second.mCellMaps[0].Size();
    const std::size_t ncells2 = second.mMesh.Cells(0).NumCells();
    bool checked_a_fresh_id = false;
    for (std::size_t p = 0; p < nparents2; ++p) {
        const std::int64_t lo = map2[p];
        const std::int64_t hi =
            p + 1 < nparents2 ? map2[p + 1] : static_cast<std::int64_t>(ncells2);
        if (hi - lo == 1)
            continue;  // untouched, kept its old (possibly small) id
        for (std::int64_t c = lo; c < hi; ++c) {
            EXPECT_GT(meshioplusplus::detail::read_int(ids2, static_cast<std::size_t>(c)),
                      max_first_id);
            checked_a_fresh_id = true;
        }
    }
    EXPECT_TRUE(checked_a_fresh_id) << "the second call must have split something to check";
}

TEST(RefineSelective, EveryParentIdResolvesInTheInputMesh) {
    Mesh m = tri_grid(3);
    RefineOptions o = cells_opt({8});
    o.mRecordHierarchy = true;
    const RefineResult res = refine(m, o);
    // m carries no refine:cell_id of its own, so its implicit ids are its
    // global block-major indices: [0, NumCells()).
    const std::int64_t coarse_n = static_cast<std::int64_t>(m.Cells(0).NumCells());
    const NDArray& parents = res.mMesh.CellData("refine:parent_id", 0);
    for (std::size_t i = 0; i < parents.Size(); ++i) {
        const std::int64_t p = meshioplusplus::detail::read_int(parents, i);
        EXPECT_GE(p, 0);
        EXPECT_LT(p, coarse_n) << "every parent_id must resolve to a cell of the coarse mesh";
    }
}

TEST(RefineSelective, RedGreenAndUntouchedAreDistinguishableFromTheTwoMeshes) {
    Mesh m = tri_grid(3);
    RefineOptions o = cells_opt({8});
    o.mRecordHierarchy = true;
    o.mRecordLevels = true;
    const RefineResult res = refine(m, o);
    const NDArray& ids = res.mMesh.CellData("refine:cell_id", 0);
    const NDArray& parents = res.mMesh.CellData("refine:parent_id", 0);
    const NDArray& levels = res.mMesh.CellData("refine:level", 0);
    // m has no PRIOR refine:level, so every parent's own level is implicitly 0
    // and the classification collapses to the rule stated in refine.hpp:
    // untouched: id == parent; green: id != parent && level == 0 (unchanged);
    // red: id != parent && level == 1 (one more than the parent's 0).
    std::size_t untouched = 0, green = 0, red = 0;
    for (std::size_t c = 0; c < ids.Size(); ++c) {
        const std::int64_t id = meshioplusplus::detail::read_int(ids, c);
        const std::int64_t parent = meshioplusplus::detail::read_int(parents, c);
        const std::int64_t level = meshioplusplus::detail::read_int(levels, c);
        if (id == parent) {
            EXPECT_EQ(level, 0) << "an untouched cell must be at level 0 too";
            ++untouched;
        } else if (level == 1) {
            ++red;
        } else {
            EXPECT_EQ(level, 0) << "a green child inherits its parent's level unchanged";
            ++green;
        }
    }
    EXPECT_GT(untouched, 0u) << "no third array is needed only if all three classes are covered";
    EXPECT_GT(green, 0u);
    EXPECT_GT(red, 0u);
}

TEST(RefineSelective, SiblingsGroupUnderSplitByTag) {
    using meshioplusplus::split;
    using meshioplusplus::SplitBy;
    Mesh m = tri_grid(3);
    RefineOptions o = cells_opt({4, 8});
    o.mRecordHierarchy = true;
    const RefineResult res = refine(m, o);
    const NDArray& parents = res.mMesh.CellData("refine:parent_id", 0);
    std::set<std::int64_t> distinct_parents;
    for (std::size_t i = 0; i < parents.Size(); ++i)
        distinct_parents.insert(meshioplusplus::detail::read_int(parents, i));

    const auto grouped = split(res.mMesh, SplitBy::Tag, "refine:parent_id");
    EXPECT_EQ(grouped.mPieces.size(), distinct_parents.size())
        << "one piece per distinct parent -- siblings grouped, untouched cells singletons";
    for (const auto& piece : grouped.mPieces) {
        const std::int64_t key = std::stoll(piece.mKey);
        EXPECT_TRUE(distinct_parents.count(key));
    }
}

TEST(RefineSelective, TwoLevelsNameTheOriginalInputLikeParentCell) {
    Mesh m = tri_grid(3);
    RefineOptions o;
    o.mLevels = 2;
    o.mRecordHierarchy = true;
    o.mRecordParentIds = true;
    const RefineResult res = refine(m, o);
    ASSERT_TRUE(res.mMesh.HasCellData("refine:parent_cell"));
    ASSERT_TRUE(res.mMesh.HasCellData("refine:parent_id"));
    const NDArray& parent_cell = res.mMesh.CellData("refine:parent_cell", 0);
    const NDArray& parent_id = res.mMesh.CellData("refine:parent_id", 0);
    // refine:parent_cell names the original ancestor's ROW INDEX; for a mesh
    // with no prior hierarchy the implicit refine:cell_id of an input cell IS
    // its row index, so the two must induce the identical partition even
    // across several internal levels.
    ASSERT_EQ(parent_cell.Size(), parent_id.Size());
    for (std::size_t i = 0; i < parent_cell.Size(); ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(parent_cell, i),
                  meshioplusplus::detail::read_int(parent_id, i));
}

TEST(RefineSelective, RecordHierarchyAlsoAttachesTheEntityKeys) {
    // RedGreen leaves no hanging nodes, so refine:entity would normally never
    // be attached at all -- but mRecordHierarchy needs it as the multigrid
    // prolongation stencil (each new fine node's coarse corners and weights).
    Mesh m = tri_grid(3);
    RefineOptions o = cells_opt({8}, RefineClosure::RedGreen);
    o.mRecordHierarchy = true;
    const RefineResult res = refine(m, o);
    ASSERT_TRUE(res.mMesh.HasPointData("refine:entity"));
    const NDArray& keys = res.mMesh.PointData("refine:entity");
    const NDArray& pts = res.mMesh.Points();
    const std::size_t dim = res.mMesh.PointDim();
    std::size_t new_nodes_checked = 0;
    for (std::size_t p = 0; p < res.mMesh.NumPoints(); ++p) {
        std::array<std::int64_t, 4> key{};
        for (std::size_t k = 0; k < 4; ++k)
            key[k] = meshioplusplus::detail::read_int(keys, p * 4 + k);
        if (key[3] < 0)
            continue;  // an original point, not one refine created
        const std::size_t first = key[0] < 0 ? 2 : 0;
        for (std::size_t k = 0; k < dim; ++k) {
            double sum = 0.0;
            for (std::size_t c = first; c < 4; ++c)
                sum += pts.As<double>()[static_cast<std::size_t>(key[c]) * dim + k];
            EXPECT_DOUBLE_EQ(pts.As<double>()[p * dim + k], sum / static_cast<double>(4 - first))
                << "point " << p << " must sit at the corner mean of its recorded entity";
        }
        ++new_nodes_checked;
    }
    EXPECT_GT(new_nodes_checked, 0u) << "the selection must have created some new node to check";
}

// A validation guard nobody tests for firing is exactly how an oracle ships
// inert (see TheConformityOraclesActuallyFire above). This proves the
// hierarchy guard checks TWO independent things -- uniqueness and
// non-negativity -- and does NOT reject a genuinely valid array, so the two
// sabotages below are not meaningless (a guard that rejected everything would
// "catch" both for the wrong reason).
TEST(RefineSelective, TheHierarchyOraclesActuallyFire) {
    const auto attach_and_refine = [](std::vector<std::int64_t> ids,
                                      std::vector<std::int64_t> parents) {
        Mesh m = tri_grid(2);
        const std::size_t n = m.Cells(0).NumCells();
        NDArray id_arr = NDArray::Uninit(meshioplusplus::DType::Int64, {n, 1});
        NDArray parent_arr = NDArray::Uninit(meshioplusplus::DType::Int64, {n, 1});
        for (std::size_t c = 0; c < n; ++c) {
            id_arr.As<std::int64_t>()[c] = ids[c];
            parent_arr.As<std::int64_t>()[c] = parents[c];
        }
        std::vector<NDArray> id_blocks, parent_blocks;
        id_blocks.push_back(std::move(id_arr));
        parent_blocks.push_back(std::move(parent_arr));
        m.AddCellData("refine:cell_id", std::move(id_blocks));
        m.AddCellData("refine:parent_id", std::move(parent_blocks));
        return refine(m, cells_opt({0})).mMesh.HasCellData("refine:cell_id");
    };
    const std::size_t n = tri_grid(2).Cells(0).NumCells();
    std::vector<std::int64_t> sequential(n), self_parent(n);
    for (std::size_t c = 0; c < n; ++c)
        sequential[c] = self_parent[c] = static_cast<std::int64_t>(c);

    // Sabotage 1: a duplicated id, otherwise well-formed. A guard checking
    // only shape/coverage would let this through.
    {
        std::vector<std::int64_t> dup = sequential;
        dup[0] = dup[1];
        EXPECT_FALSE(attach_and_refine(dup, self_parent))
            << "the uniqueness check must have fired on the duplicated id";
    }
    // Sabotage 2: a negative id, unique otherwise. A guard checking only
    // uniqueness (e.g. via a hash set) would let this through.
    {
        std::vector<std::int64_t> negative = sequential;
        negative[0] = -1;
        EXPECT_FALSE(attach_and_refine(negative, self_parent))
            << "the non-negativity check must have fired on the negative id";
    }
    // Control: genuinely valid input must survive, or the two sabotages above
    // prove nothing (a guard rejecting everything "catches" both trivially).
    EXPECT_TRUE(attach_and_refine(sequential, self_parent))
        << "a genuinely valid hierarchy must be maintained, not rejected";
}

}  // namespace selective

}  // namespace
