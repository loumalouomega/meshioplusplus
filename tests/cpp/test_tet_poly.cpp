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
// Tests for tetrahedral <-> polyhedral conversion (`operations/tet_poly.hpp`).
//
// The dual's correctness rests on three checkable properties: every dual cell
// is a closed surface already wound outward (orient_rings reports Consistent,
// so the winding is right by construction, not repaired), every face is shared
// by exactly two cells or lies on the boundary (conformity), and the cells'
// volumes sum to the tetrahedral mesh's (planar boundary faces make it exact).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/tet_poly.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::poly_to_tet;
using meshioplusplus::PolyToTetOptions;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;
using meshioplusplus::tet_to_poly;
using meshioplusplus::TetToPolyOptions;
namespace detail = meshioplusplus::detail;

namespace {

NDArray i64(const std::vector<std::int64_t>& rVals) {
    NDArray a = NDArray::Uninit(DType::Int64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

NDArray f64(const std::vector<double>& rVals) {
    NDArray a = NDArray::Uninit(DType::Float64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<double>()[i] = rVals[i];
    return a;
}

double tet_volume_sum(const Mesh& rMesh) {
    double vol = 0.0;
    const NDArray& pts = rMesh.Points();
    for (const auto cb : rMesh.CellRange()) {
        if (std::string(cb.Type()) != "tetra")
            continue;
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            double p[4][3];
            for (int k = 0; k < 4; ++k)
                for (int d = 0; d < 3; ++d)
                    p[k][d] = detail::read_double(
                        pts, static_cast<std::size_t>(detail::read_int(conn, 4 * c + k)) * 3 + d);
            double u[3][3];
            for (int r = 0; r < 3; ++r)
                for (int d = 0; d < 3; ++d)
                    u[r][d] = p[r + 1][d] - p[0][d];
            vol += (u[0][0] * (u[1][1] * u[2][2] - u[1][2] * u[2][1]) -
                    u[0][1] * (u[1][0] * u[2][2] - u[1][2] * u[2][0]) +
                    u[0][2] * (u[1][0] * u[2][1] - u[1][1] * u[2][0])) /
                   6.0;
        }
    }
    return vol;
}

/// Checks every polyhedron is closed and already outward-wound, and returns
/// the total volume.
double dual_volume(const Mesh& rMesh) {
    double vol = 0.0;
    detail::CellRings rings;
    std::vector<detail::Vec3> coords;
    for (const auto cb : rMesh.CellRange()) {
        if (!cb.IsPolyhedron())
            continue;
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            EXPECT_TRUE(detail::cell_rings(cb, c, rMesh.Points(), rMesh.PointDim(), rings, coords));
            EXPECT_EQ(detail::orient_rings(rings, coords.data()),
                      detail::RingOrientation::Consistent)
                << "cell " << c << " is not closed and outward-wound as built";
            const detail::PolyMeasure m = detail::poly_measure(rings, coords.data());
            EXPECT_TRUE(m.mClosed);
            EXPECT_GT(m.mVolume, 0.0);
            vol += m.mVolume;
        }
    }
    return vol;
}

/// Face (as its sorted node set) -> how many polyhedron faces use it, and the
/// number of times it is used forward minus reversed (0 for a shared face).
std::map<std::vector<std::int64_t>, int> face_uses(const Mesh& rMesh) {
    std::map<std::vector<std::int64_t>, int> uses;
    for (const auto cb : rMesh.CellRange()) {
        if (!cb.IsPolyhedron())
            continue;
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                const auto face = cb.Face(c, f);
                std::vector<std::int64_t> key(face.first, face.first + face.second);
                std::sort(key.begin(), key.end());
                ++uses[key];
            }
    }
    return uses;
}

bool polyhedron_references(const Mesh& rMesh, std::size_t Block, std::size_t Cell,
                           std::int64_t Node) {
    const auto cb = rMesh.Cells(Block);
    for (std::size_t f = 0; f < cb.NumFaces(Cell); ++f) {
        const auto face = cb.Face(Cell, f);
        for (std::size_t k = 0; k < face.second; ++k)
            if (face.first[k] == Node)
                return true;
    }
    return false;
}

Mesh single_tet() {
    return mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra", {{0, 1, 2, 3}});
}

Mesh kuhn_box(std::int64_t N) {
    return meshioplusplus::box({{N, N, N}}, {{1.0, 1.0, 1.0}}, {{0.0, 0.0, 0.0}}, "tetra");
}

}  // namespace

// ---------------------------------------------------------------------------
// tet_to_poly
// ---------------------------------------------------------------------------

TEST(TetToPoly, SingleTetHasFourClosedCellsAndConservesVolume) {
    const Mesh m = single_tet();
    const auto r = tet_to_poly(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    const auto cb = r.mMesh.Cells(0);
    EXPECT_TRUE(cb.IsPolyhedron());
    EXPECT_EQ(std::string(cb.Type()), "polyhedron");
    ASSERT_EQ(cb.NumCells(), 4u);
    for (std::size_t c = 0; c < 4; ++c)
        EXPECT_EQ(cb.NumFaces(c), 6u);  // three edge faces, three boundary faces
    EXPECT_EQ(r.mNumBoundaryFaces, 12);
    EXPECT_EQ(r.mNumSplitVertices, 0);
    // 6 edge midpoints + 4 face centroids + 1 tet centroid.
    EXPECT_EQ(r.mNumNewPoints, 11);
    EXPECT_EQ(r.mMesh.NumPoints(), 15u);
    EXPECT_NEAR(dual_volume(r.mMesh), 1.0 / 6.0, 1e-14);
    ASSERT_EQ(r.mVertexMap.Size(), 4u);
    for (std::int64_t c = 0; c < 4; ++c)
        EXPECT_EQ(r.mVertexMap.As<std::int64_t>()[c], c);
}

TEST(TetToPoly, KuhnBoxIsConformingAndConservesVolume) {
    const Mesh m = kuhn_box(3);
    ASSERT_EQ(m.NumPoints(), 64u);
    const auto r = tet_to_poly(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    const auto cb = r.mMesh.Cells(0);
    ASSERT_EQ(cb.NumCells(), 64u);  // one cell per vertex
    EXPECT_NEAR(dual_volume(r.mMesh), tet_volume_sum(m), 1e-13);
    EXPECT_NEAR(dual_volume(r.mMesh), 1.0, 1e-13);

    // Face-interior vertices get one boundary face (a closed run), edge
    // vertices two, corners three: 24 + 2 * 24 + 3 * 8.
    EXPECT_EQ(r.mNumBoundaryFaces, 96);
    int shared = 0;
    int single = 0;
    for (const auto& [key, n] : face_uses(r.mMesh)) {
        EXPECT_LE(n, 2);
        (n == 2 ? shared : single) += 1;
    }
    EXPECT_EQ(single, 96);
    EXPECT_GT(shared, 0);

    // A cell references its own vertex only where a boundary run stops short
    // of going all the way round -- on a box edge or corner -- and never any
    // other input vertex.
    const NDArray& pts = m.Points();
    const std::int64_t* vmap = r.mVertexMap.As<std::int64_t>();
    for (std::size_t c = 0; c < cb.NumCells(); ++c) {
        const std::int64_t v = vmap[c];
        int on_boundary = 0;
        for (int d = 0; d < 3; ++d) {
            const double x = detail::read_double(pts, static_cast<std::size_t>(v) * 3 + d);
            on_boundary += x < 1e-9 || x > 1.0 - 1e-9;
        }
        EXPECT_EQ(polyhedron_references(r.mMesh, 0, c, v), on_boundary >= 2) << "vertex " << v;
        for (std::int64_t u = 0; u < 64; ++u)
            if (u != v)
                EXPECT_FALSE(polyhedron_references(r.mMesh, 0, c, u));
    }
}

TEST(TetToPoly, ConnectivityOrderDoesNotMatter) {
    // Half the tets listed in negative orientation: the dual must not care,
    // since one consistent orientation is propagated across shared faces.
    const Mesh box = kuhn_box(2);
    const auto cb = box.Cells(0);
    std::vector<std::vector<std::int64_t>> rows;
    for (std::size_t c = 0; c < cb.NumCells(); ++c) {
        std::vector<std::int64_t> row(4);
        for (int k = 0; k < 4; ++k)
            row[k] = detail::read_int(cb.Conn(), 4 * c + k);
        if (c % 2 == 1)
            std::swap(row[1], row[2]);
        rows.push_back(row);
    }
    Mesh m;
    m.AssignPoints(detail::data_owned_copy(box.Points()));
    m.AddCellBlock("tetra", mt::conn_from(rows));
    const auto a = tet_to_poly(box);
    const auto b = tet_to_poly(m);
    ASSERT_EQ(b.mMesh.Cells(0).NumCells(), 27u);
    EXPECT_EQ(b.mNumBoundaryFaces, a.mNumBoundaryFaces);
    EXPECT_NEAR(dual_volume(b.mMesh), 1.0, 1e-13);
    for (const auto& [key, n] : face_uses(b.mMesh))
        EXPECT_LE(n, 2);
}

TEST(TetToPoly, LargeFeatureAngleFusesEdgeRunsButKeepsCorners) {
    const Mesh m = kuhn_box(2);
    TetToPolyOptions opts;
    opts.mFeatureAngleDeg = 95.0;  // the box's 90-degree edges no longer break runs
    const auto r = tet_to_poly(m, opts);
    // Every boundary vertex's run then goes all the way round: one face each.
    EXPECT_EQ(r.mNumBoundaryFaces, 26);
    for (const auto cb : r.mMesh.CellRange())
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            detail::CellRings rings;
            std::vector<detail::Vec3> coords;
            ASSERT_TRUE(detail::cell_rings(cb, c, r.mMesh.Points(), 3, rings, coords));
            EXPECT_NE(detail::orient_rings(rings, coords.data()),
                      detail::RingOrientation::Unorientable);
        }
}

TEST(TetToPoly, TwoBlocksSplitInterfaceVerticesAndStayConforming) {
    const Mesh box = kuhn_box(2);
    // Split the tets by centroid x into two blocks.
    const auto cb = box.Cells(0);
    const NDArray& conn = cb.Conn();
    std::vector<std::vector<std::int64_t>> left;
    std::vector<std::vector<std::int64_t>> right;
    for (std::size_t c = 0; c < cb.NumCells(); ++c) {
        std::vector<std::int64_t> row(4);
        double cx = 0.0;
        for (int k = 0; k < 4; ++k) {
            row[k] = detail::read_int(conn, 4 * c + k);
            cx += detail::read_double(box.Points(), static_cast<std::size_t>(row[k]) * 3);
        }
        (cx / 4.0 < 0.5 ? left : right).push_back(row);
    }
    Mesh m;
    m.AssignPoints(detail::data_owned_copy(box.Points()));
    m.AddCellBlock("tetra", mt::conn_from(left));
    m.AddCellBlock("tetra", mt::conn_from(right));

    const auto r = tet_to_poly(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 2u);
    EXPECT_EQ(r.mNumSplitVertices, 9);  // the 3 x 3 vertices on x = 0.5
    EXPECT_EQ(r.mMesh.Cells(0).NumCells() + r.mMesh.Cells(1).NumCells(), 27u + 9u);
    EXPECT_NEAR(dual_volume(r.mMesh), 1.0, 1e-13);
    for (const auto& [key, n] : face_uses(r.mMesh))
        EXPECT_LE(n, 2);
    // Faces used once are the mesh boundary; the remaining boundary faces are
    // interface faces, each matched (reversed) by the other side: the 9
    // interface vertices get one coplanar run per side.
    int single = 0;
    for (const auto& [key, n] : face_uses(r.mMesh))
        single += n == 1;
    EXPECT_EQ(r.mNumBoundaryFaces - single, 18);
}

TEST(TetToPoly, DataRegionsAndTriangleBlocks) {
    Mesh m = mt::tet_mesh();  // two tets of equal volume sharing face (0, 2, 4)
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}, {0, 2, 3}}));  // the base
    m.AddCellData("w", {f64({1.0, 3.0}), f64({10.0, 20.0})});
    m.AddCellData("id", {i64({7, 9}), i64({1, 2})});
    m.AddPointData("x", f64({0.0, 1.0, 1.0, 0.0, 0.5}));
    m.AddRegion(Region("apex", RegionKind::Point, i64({4})));
    m.AddRegion(Region("both", RegionKind::Cell, i64({0, 1})));

    TetToPolyOptions opts;
    opts.mPointDataToCellData = true;
    opts.mRecordVertexIds = true;
    const auto r = tet_to_poly(m, opts);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 2u);
    const auto poly = r.mMesh.Cells(0);
    const auto bnd = r.mMesh.Cells(1);
    EXPECT_TRUE(poly.IsPolyhedron());
    EXPECT_EQ(std::string(bnd.Type()), "polygon");
    ASSERT_EQ(poly.NumCells(), 5u);
    // Vertices 0..3 each touch the base square, split by the diagonal 0-2 into
    // two coplanar triangles: 0 and 2 see both (one fused run), 1 and 3 one.
    EXPECT_EQ(bnd.NumCells(), 4u);

    const NDArray& vid = r.mMesh.CellData("dual:vertex", 0);
    const NDArray& w = r.mMesh.CellData("w", 0);
    const NDArray& id = r.mMesh.CellData("id", 0);
    const NDArray& x = r.mMesh.CellData("x", 0);
    for (std::size_t c = 0; c < poly.NumCells(); ++c) {
        const std::int64_t v = vid.As<std::int64_t>()[c];
        const double wv = w.As<double>()[c];
        const std::int64_t iv = id.As<std::int64_t>()[c];
        if (v == 1) {
            EXPECT_DOUBLE_EQ(wv, 1.0);
            EXPECT_EQ(iv, 7);
        } else if (v == 3) {
            EXPECT_DOUBLE_EQ(wv, 3.0);
            EXPECT_EQ(iv, 9);
        } else {
            EXPECT_DOUBLE_EQ(wv, 2.0);  // equal volumes
            EXPECT_EQ(iv, 7);           // a tie goes to the lower index
        }
        EXPECT_DOUBLE_EQ(x.As<double>()[c], m.PointData("x").As<double>()[v]);
    }
    // Boundary polygons carry the triangle block's data by area.
    const NDArray& wb = r.mMesh.CellData("w", 1);
    const NDArray& vb = r.mMesh.CellData("dual:vertex", 1);
    for (std::size_t c = 0; c < bnd.NumCells(); ++c) {
        const std::int64_t v = vb.As<std::int64_t>()[c];
        const double expect = v == 1 ? 10.0 : (v == 3 ? 20.0 : 15.0);
        EXPECT_DOUBLE_EQ(wb.As<double>()[c], expect);
    }
    // point_data is the source mean on new points: the last point is the
    // second tet's centroid, (0 + 0.5 + 1 + 0) / 4 in x.
    const NDArray& px = r.mMesh.PointData("x");
    ASSERT_EQ(detail::rows(px), r.mMesh.NumPoints());
    EXPECT_DOUBLE_EQ(px.As<double>()[r.mMesh.NumPoints() - 1], (0.0 + 1.0 + 0.0 + 0.5) / 4.0);
    // Regions: the point region survives; the cell region comes back empty.
    ASSERT_NE(r.mMesh.FindRegion("apex", RegionKind::Point), Mesh::npos);
    EXPECT_EQ(r.mMesh.Region(r.mMesh.FindRegion("apex", RegionKind::Point)).NumEntries(), 1u);
    const std::size_t both = r.mMesh.FindRegion("both", RegionKind::Cell);
    if (both != Mesh::npos)
        EXPECT_EQ(r.mMesh.Region(both).NumEntries(), 0u);
    EXPECT_NEAR(dual_volume(r.mMesh), tet_volume_sum(m), 1e-14);
}

TEST(TetToPoly, Refusals) {
    EXPECT_THROW(tet_to_poly(mt::hex_mesh()), std::invalid_argument);
    EXPECT_THROW(tet_to_poly(mt::tet10_mesh()), std::invalid_argument);
    EXPECT_THROW(tet_to_poly(mt::tri_mesh()), std::invalid_argument);
    EXPECT_THROW(tet_to_poly(tet_to_poly(single_tet()).mMesh), std::invalid_argument);
    // Three tets on one face: non-manifold.
    const Mesh nm =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}, {1, 1, 1}}, "tetra",
                      {{0, 1, 2, 3}, {0, 2, 1, 4}, {0, 1, 2, 5}});
    EXPECT_THROW(tet_to_poly(nm), std::invalid_argument);
    TetToPolyOptions bad;
    bad.mFeatureAngleDeg = 180.0;
    EXPECT_THROW(tet_to_poly(single_tet(), bad), std::invalid_argument);
}

TEST(TetToPoly, IsDeterministic) {
    const Mesh m = kuhn_box(3);
    const auto a = tet_to_poly(m);
    const auto b = tet_to_poly(m);
    ASSERT_EQ(a.mMesh.NumPoints(), b.mMesh.NumPoints());
    EXPECT_EQ(
        std::memcmp(a.mMesh.Points().Data(), b.mMesh.Points().Data(), a.mMesh.Points().Nbytes()),
        0);
    const auto ca = a.mMesh.Cells(0);
    const auto cbb = b.mMesh.Cells(0);
    ASSERT_EQ(ca.NumCells(), cbb.NumCells());
    for (std::size_t c = 0; c < ca.NumCells(); ++c) {
        ASSERT_EQ(ca.NumFaces(c), cbb.NumFaces(c));
        for (std::size_t f = 0; f < ca.NumFaces(c); ++f) {
            const auto fa = ca.Face(c, f);
            const auto fb = cbb.Face(c, f);
            ASSERT_EQ(fa.second, fb.second);
            EXPECT_TRUE(std::equal(fa.first, fa.first + fa.second, fb.first));
        }
    }
}

// ---------------------------------------------------------------------------
// poly_to_tet
// ---------------------------------------------------------------------------

TEST(PolyToTet, TetrahedraStoredAsPolyhedraComeBackExactly) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}));
    // Two tets, the second with deliberately inward-wound faces.
    m.AddPolyhedronBlock("polyhedron4", std::vector<std::vector<std::vector<std::int64_t>>>{
                                            {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}},
                                            {{1, 2, 3}, {1, 4, 2}, {2, 4, 3}, {3, 4, 1}}});
    const auto r = poly_to_tet(m);
    EXPECT_EQ(r.mNumExactTets, 2);
    EXPECT_EQ(r.mNumFannedCells, 0);
    ASSERT_EQ(r.mMesh.NumPoints(), 5u);  // no new points
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    const auto cb = r.mMesh.Cells(0);
    EXPECT_EQ(std::string(cb.Type()), "tetra");
    ASSERT_EQ(cb.NumCells(), 2u);
    const double v = tet_volume_sum(r.mMesh);
    // Each tet positive: |volumes| sum equals the signed sum.
    EXPECT_NEAR(v, 1.0 / 6.0 + 1.0 / 3.0, 1e-14);
    const NDArray& conn = cb.Conn();
    std::vector<std::int64_t> first(4);
    for (int k = 0; k < 4; ++k)
        first[k] = detail::read_int(conn, k);
    std::sort(first.begin(), first.end());
    EXPECT_EQ(first, (std::vector<std::int64_t>{0, 1, 2, 3}));
}

TEST(PolyToTet, MixedPolyhedraFanOnlyTheNonTetrahedra) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1},
                                    {2, 0, 0}}));
    m.AddPolyhedronBlock(
        "polyhedron",
        std::vector<std::vector<std::vector<std::int64_t>>>{
            {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}},
            {{1, 8, 2}, {1, 5, 8}, {2, 8, 5}, {1, 2, 5}}});
    PolyToTetOptions opts;
    opts.mRecordParentIds = true;
    const auto r = poly_to_tet(m, opts);
    EXPECT_EQ(r.mNumExactTets, 1);
    EXPECT_EQ(r.mNumFannedCells, 1);
    const auto cb = r.mMesh.Cells(0);
    EXPECT_EQ(cb.NumCells(), 24u + 1u);
    EXPECT_EQ(r.mMesh.NumPoints(), 9u + 7u);  // cube centroid + 6 face centroids
    EXPECT_NEAR(tet_volume_sum(r.mMesh), 1.0 + 1.0 / 6.0, 1e-13);
    ASSERT_EQ(r.mCellMaps.size(), 1u);
    EXPECT_EQ(r.mCellMaps[0].As<std::int64_t>()[1], 24);
}

TEST(PolyToTet, RoundTripThroughTheDualConservesVolume) {
    const Mesh m = kuhn_box(2);
    const auto dual = tet_to_poly(m);
    const auto back = poly_to_tet(dual.mMesh);
    EXPECT_EQ(back.mNumExactTets, 0);
    EXPECT_EQ(back.mNumFannedCells, 27);
    EXPECT_NEAR(tet_volume_sum(back.mMesh), 1.0, 1e-13);
    // Every tet positively oriented: the signed sum equals the absolute sum.
    double abs_sum = 0.0;
    const auto cb = back.mMesh.Cells(0);
    for (std::size_t c = 0; c < cb.NumCells(); ++c) {
        Mesh one;
        one.AssignPoints(detail::data_owned_copy(back.mMesh.Points()));
        std::vector<std::int64_t> row(4);
        for (int k = 0; k < 4; ++k)
            row[k] = detail::read_int(cb.Conn(), 4 * c + k);
        one.AddCellBlock("tetra", mt::conn_from({row}));
        const double v = tet_volume_sum(one);
        EXPECT_GE(v, -1e-15);
        abs_sum += std::abs(v);
    }
    EXPECT_NEAR(abs_sum, 1.0, 1e-13);
}
