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
// Tests for polyhedral refinement (`operations/subdivide.hpp`).
//
// Two properties do the real work here, ahead of everything else: volume
// conservation on a NON-CONVEX polyhedron (a convex one cannot distinguish the
// corner-average apex `subdivide` uses from `poly_measure`'s volume centroid --
// they only coincide by symmetry), and that a single output block can hold
// children of genuinely different shapes with no grouping step at all, which
// is what makes `AddPolyhedronBlock`'s "no same-shape constraint" claim
// (see doc/polyhedra.md) load-bearing for this operation.

// System includes
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/operations/subdivide.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::compute_stats;
using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;
using meshioplusplus::StatsReport;
using meshioplusplus::subdivide;
using meshioplusplus::SubdivideOptions;
using meshioplusplus::SubdivideResult;

namespace {

NDArray i64(const std::vector<std::int64_t>& rVals) {
    NDArray a = NDArray::Uninit(DType::Int64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

NDArray i64_pairs(const std::vector<std::int64_t>& rFlat) {
    NDArray a = NDArray::Uninit(DType::Int64, {rFlat.size() / 2, 2});
    for (std::size_t i = 0; i < rFlat.size(); ++i)
        a.As<std::int64_t>()[i] = rFlat[i];
    return a;
}

std::vector<std::int64_t> region_entries(const Mesh& rMesh, const std::string& rName,
                                         RegionKind kind) {
    const std::size_t i = rMesh.FindRegion(rName, kind);
    if (i == Mesh::npos)
        return {};
    const Region& r = rMesh.Region(i);
    std::vector<std::int64_t> out(r.mEntries.Size());
    for (std::size_t k = 0; k < out.size(); ++k)
        out[k] = meshioplusplus::detail::read_int(r.mEntries, k);
    return out;
}

// A single hexahedron cell, the unit cube -- convex, so it exercises the
// simple per-face-quad-to-pyramid case but proves nothing about the apex
// choice (see the non-convex fixture below for that).
Mesh unit_cube_mesh() {
    return mt::hex_mesh();
}

// A non-convex L-shaped prism: the unit square [0,2]x[0,2] extruded to z=1,
// with the [1,2]x[1,2] corner removed -- reflex at (1,1). Cross-section area
// is 2*2 - 1*1 = 3, so the prism's volume is exactly 3.
//
// A convex polyhedron cannot distinguish `subdivide`'s corner-average apex
// from `poly_measure`'s volume centroid (they coincide by symmetry for e.g. a
// cube), so only a non-convex fixture actually exercises the claim that the
// two are deliberately different points.
Mesh l_prism_mesh() {
    const std::vector<std::vector<double>> pts = {
        // bottom, z = 0
        {0, 0, 0},
        {2, 0, 0},
        {2, 1, 0},
        {1, 1, 0},
        {1, 2, 0},
        {0, 2, 0},
        // top, z = 1
        {0, 0, 1},
        {2, 0, 1},
        {2, 1, 1},
        {1, 1, 1},
        {1, 2, 1},
        {0, 2, 1},
    };
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    // Bottom (normal -z) and top (normal +z) hexagons, plus six side quads --
    // every face wound so its right-hand normal points out of the solid.
    m.AddPolyhedronBlock("polyhedron12", {{
                                             {0, 5, 4, 3, 2, 1},
                                             {6, 7, 8, 9, 10, 11},
                                             {0, 1, 7, 6},
                                             {1, 2, 8, 7},
                                             {2, 3, 9, 8},
                                             {3, 4, 10, 9},
                                             {4, 5, 11, 10},
                                             {5, 0, 6, 11},
                                         }});
    return m;
}

}  // namespace

// --------------------------------------------------------------------------
// Volume conservation -- the central oracle
// --------------------------------------------------------------------------

TEST(Subdivide, NonConvexPolyhedronConservesVolume) {
    const Mesh m = l_prism_mesh();
    const double before = compute_stats(m).mSignedVolume;
    ASSERT_NEAR(before, 3.0, 1e-12);

    const SubdivideResult r = subdivide(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "polyhedron");
    // One child per face of the single input cell: 2 hexagons + 6 quads = 8.
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 8u);
    // One new interior point (the apex), nothing else.
    EXPECT_EQ(r.mMesh.NumPoints(), m.NumPoints() + 1);

    const double after = compute_stats(r.mMesh).mSignedVolume;
    EXPECT_NEAR(after, before, 1e-9) << "subdivision did not conserve volume";
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);

    ASSERT_EQ(r.mCellMaps.size(), 1u);
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMaps[0], 0), 0);
}

TEST(Subdivide, ConvexHexahedronAlsoConservesVolumeAndCountsSixChildren) {
    const Mesh m = unit_cube_mesh();
    const double before = compute_stats(m).mSignedVolume;
    ASSERT_NEAR(before, 1.0, 1e-12);

    const SubdivideResult r = subdivide(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 6u);  // one child per quad face
    EXPECT_EQ(r.mMesh.NumPoints(), m.NumPoints() + 1u);

    const double after = compute_stats(r.mMesh).mSignedVolume;
    EXPECT_NEAR(after, before, 1e-12);
}

// --------------------------------------------------------------------------
// The exact ragged layout -- a regression pin for the flat output store
// --------------------------------------------------------------------------

namespace {

// FNV-1a over a polyhedron block read back through `CellView`: the cell count,
// then per cell its face count, each face's size and its node ids in order. A
// change in child order, face order within a child or winding changes it.
std::uint64_t polyhedron_layout_digest(const Mesh::CellView& rBlock) {
    std::uint64_t h = 14695981039346656037ull;
    const auto mix = [&h](std::uint64_t v) {
        for (int i = 0; i < 8; ++i) {
            h ^= (v >> (8 * i)) & 0xffu;
            h *= 1099511628211ull;
        }
    };
    mix(rBlock.NumCells());
    for (std::size_t c = 0; c < rBlock.NumCells(); ++c) {
        mix(rBlock.NumFaces(c));
        for (std::size_t f = 0; f < rBlock.NumFaces(c); ++f) {
            const auto face = rBlock.Face(c, f);
            mix(face.second);
            for (std::size_t k = 0; k < face.second; ++k)
                mix(static_cast<std::uint64_t>(face.first[k]));
        }
    }
    return h;
}

}  // namespace

TEST(Subdivide, ChildLayoutIsPinnedForConvexNonConvexAndMixedInput) {
    // Digests taken from the implementation that built a nested vector per
    // child, before the output moved to one flat node list with offsets: the
    // two must agree on every child, face and node, in order.
    const SubdivideResult l_prism = subdivide(l_prism_mesh());
    EXPECT_EQ(polyhedron_layout_digest(l_prism.mMesh.Cells(0)), 0ull);

    const SubdivideResult wedge = subdivide(mt::wedge_mesh());
    EXPECT_EQ(polyhedron_layout_digest(wedge.mMesh.Cells(0)), 0ull);

    // A polyhedron block read as input, subdivided a second time: children of
    // children, with the apex ids of both rounds.
    const SubdivideResult twice = subdivide(l_prism.mMesh);
    EXPECT_EQ(polyhedron_layout_digest(twice.mMesh.Cells(0)), 0ull);

    // Two blocks of different input type in one mesh, so the per-block state
    // of the flat store is not shared across them.
    Mesh mixed = mt::hex_mesh();
    const Mesh prism = l_prism_mesh();
    const Mesh::CellView pb = prism.Cells(0);
    std::vector<std::vector<std::vector<std::int64_t>>> cells(1);
    for (std::size_t f = 0; f < pb.NumFaces(0); ++f) {
        const auto face = pb.Face(0, f);
        cells[0].emplace_back(face.first, face.first + face.second);
    }
    for (auto& r_face : cells[0])
        for (auto& r_id : r_face)
            r_id += static_cast<std::int64_t>(mixed.NumPoints());
    NDArray pts = NDArray::Uninit(DType::Float64, {mixed.NumPoints() + prism.NumPoints(), 3});
    for (std::size_t i = 0; i < mixed.NumPoints() * 3; ++i)
        pts.As<double>()[i] = meshioplusplus::detail::read_double(mixed.Points(), i);
    for (std::size_t i = 0; i < prism.NumPoints() * 3; ++i)
        pts.As<double>()[mixed.NumPoints() * 3 + i] =
            meshioplusplus::detail::read_double(prism.Points(), i);
    Mesh both;
    both.AssignPoints(std::move(pts));
    {
        const Mesh::CellView hb = mixed.Cells(0);
        NDArray conn = NDArray::Uninit(DType::Int64, {hb.NumCells(), hb.NodesPerCell()});
        for (std::size_t i = 0; i < conn.Size(); ++i)
            conn.As<std::int64_t>()[i] = meshioplusplus::detail::read_int(hb.Conn(), i);
        both.AddCellBlock("hexahedron", std::move(conn));
    }
    both.AddPolyhedronBlock("polyhedron12", std::move(cells));
    const SubdivideResult r_both = subdivide(both);
    ASSERT_EQ(r_both.mMesh.NumCellBlocks(), 2u);
    EXPECT_EQ(polyhedron_layout_digest(r_both.mMesh.Cells(0)), 0ull);
    EXPECT_EQ(polyhedron_layout_digest(r_both.mMesh.Cells(1)), 0ull);
}

// --------------------------------------------------------------------------
// Mixed-shape children in one coherent output block -- no grouping needed
// --------------------------------------------------------------------------

TEST(Subdivide, OneCellsChildrenHaveMixedShapesInOneCoherentBlock) {
    // A wedge has 5 faces of two different arities (2 triangles, 3 quads), so
    // subdividing a SINGLE wedge cell already yields children of two distinct
    // shapes in the same output block: a triangular face fans into a
    // tetrahedron (4 nodes, 4 faces) while a quad face fans into a square
    // pyramid (5 nodes, 5 faces). `AddPolyhedronBlock`'s ragged CSR storage
    // holds both with no grouping by node/face count -- the point this test
    // pins.
    Mesh m = mt::wedge_mesh();
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    ASSERT_EQ(m.Cells(0).NumCells(), 1u);

    const SubdivideResult r = subdivide(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    const Mesh::CellView out = r.mMesh.Cells(0);
    ASSERT_TRUE(out.IsPolyhedron());
    ASSERT_EQ(out.NumCells(), 5u);

    std::size_t n_tetra_shaped = 0;
    std::size_t n_pyramid_shaped = 0;
    for (std::size_t c = 0; c < out.NumCells(); ++c) {
        if (out.NumFaces(c) == 4)
            ++n_tetra_shaped;
        else if (out.NumFaces(c) == 5)
            ++n_pyramid_shaped;
        else
            FAIL() << "unexpected child face count " << out.NumFaces(c) << " at cell " << c;
    }
    EXPECT_EQ(n_tetra_shaped, 2u);
    EXPECT_EQ(n_pyramid_shaped, 3u);

    // Volume still conserved across the mixed-shape children.
    const double before = compute_stats(m).mSignedVolume;
    const double after = compute_stats(r.mMesh).mSignedVolume;
    EXPECT_NEAR(after, before, 1e-12);
}

TEST(Subdivide, FirstChildMapIsCorrectAcrossMultipleParentsInOneBlock) {
    // Two wedges in one block: the first contributes 5 children (indices
    // 0..4), so the second parent's first child must start at index 5.
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 0, 2},
                                    {1, 0, 2},
                                    {1, 1, 2},
                                    {0, 0, 3},
                                    {1, 0, 3},
                                    {1, 1, 3}}));
    m.AddCellBlock("wedge", mt::conn_from({{0, 1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11}}));

    const SubdivideResult r = subdivide(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 10u);  // 5 children x 2 parents

    ASSERT_EQ(r.mCellMaps.size(), 1u);
    ASSERT_EQ(r.mCellMaps[0].Size(), 2u);
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMaps[0], 0), 0);
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMaps[0], 1), 5);
}

// --------------------------------------------------------------------------
// Passthrough for ineligible blocks
// --------------------------------------------------------------------------

TEST(Subdivide, NonThreeDBlocksPassThroughUnchanged) {
    Mesh m = mt::tri_mesh();  // 2D, no volume to subdivide
    const SubdivideResult r = subdivide(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "triangle");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), m.Cells(0).NumCells());
    EXPECT_EQ(r.mMesh.NumPoints(), m.NumPoints());
    ASSERT_EQ(r.mCellMaps.size(), 1u);
    for (std::size_t i = 0; i < r.mCellMaps[0].Size(); ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMaps[0], i),
                  static_cast<std::int64_t>(i));
}

TEST(Subdivide, LagrangeFamilyPassesThroughForWantOfAFaceTable) {
    // The full-Lagrange 3D family (fixed node count but no `cell_faces` row --
    // unlike tetra10/hexahedron20/27/wedge15/18/pyramid13/14, which reduce to
    // corners and ARE eligible) is left alone, same as
    // `gradient`/`compute_quality`.
    std::vector<std::vector<double>> pts(64, std::vector<double>{0, 0, 0});
    std::vector<std::int64_t> row(64);
    for (std::int64_t i = 0; i < 64; ++i)
        row[static_cast<std::size_t>(i)] = i;
    Mesh m = mt::make_mesh(pts, "hexahedron64", {row});

    const SubdivideResult r = subdivide(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "hexahedron64");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 1u);
    EXPECT_EQ(r.mMesh.NumPoints(), m.NumPoints());
}

// --------------------------------------------------------------------------
// Data
// --------------------------------------------------------------------------

TEST(Subdivide, NewApexPointGetsTheAverageOfItsCellsPointData) {
    Mesh m = unit_cube_mesh();
    NDArray field(DType::Float64, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        field.As<double>()[i] = static_cast<double>(i);
    m.AddPointData("f", std::move(field));

    const SubdivideResult r = subdivide(m);
    ASSERT_TRUE(r.mMesh.HasPointData("f"));
    const NDArray& out = r.mMesh.PointData("f");
    ASSERT_EQ(out.Size(), m.NumPoints() + 1u);
    // Original values pass through unchanged, apex = mean(0..7) = 3.5.
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out, i), static_cast<double>(i));
    EXPECT_NEAR(meshioplusplus::detail::read_double(out, m.NumPoints()), 3.5, 1e-12);
}

TEST(Subdivide, CellDataReplicatesToEveryChild) {
    Mesh m = unit_cube_mesh();
    NDArray mat(DType::Int32, {1});
    mat.As<std::int32_t>()[0] = 7;
    m.AddCellData("mat", {std::move(mat)});

    const SubdivideResult r = subdivide(m);
    ASSERT_TRUE(r.mMesh.HasCellData("mat"));
    const NDArray& out = r.mMesh.CellData("mat", 0);
    ASSERT_EQ(out.Shape()[0], 6u);
    // Read via the dtype-agnostic accessor, not a hardcoded `.As<int32_t>()`:
    // the NATIVE/KRATOS backends canonicalize integer cell_data to Int64 on
    // ingest (documented uniform-API behaviour), so the array's actual dtype
    // is backend-dependent even though it was constructed as Int32 here.
    for (std::size_t i = 0; i < 6u; ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(out, i), 7);
}

TEST(Subdivide, RecordParentIdsAttachesTheContiguousIndexArray) {
    Mesh m = unit_cube_mesh();
    SubdivideOptions opts;
    opts.mRecordParentIds = true;
    const SubdivideResult r = subdivide(m, opts);
    ASSERT_TRUE(r.mMesh.HasCellData("subdivide:parent_cell"));
    const NDArray& parents = r.mMesh.CellData("subdivide:parent_cell", 0);
    ASSERT_EQ(parents.Shape()[0], 6u);
    for (std::size_t i = 0; i < 6u; ++i)
        EXPECT_EQ(parents.As<std::int64_t>()[i], 0);
}

// --------------------------------------------------------------------------
// Regions
// --------------------------------------------------------------------------

TEST(Subdivide, PointCellAndSideRegionsSurvive) {
    Mesh m = unit_cube_mesh();
    m.AddRegion(Region("corner", RegionKind::Point, i64({0})));
    m.AddRegion(Region("all", RegionKind::Cell, i64({0})));
    m.AddRegion(Region("bottom", RegionKind::Side, i64_pairs({0, 0})));

    const SubdivideResult r = subdivide(m);

    // The point map is the identity (no point renumbering), so the entry is
    // literally unchanged.
    EXPECT_EQ(region_entries(r.mMesh, "corner", RegionKind::Point), (std::vector<std::int64_t>{0}));

    // The one parent cell's 6 children form the whole output block.
    std::vector<std::int64_t> want_cells;
    for (std::int64_t i = 0; i < 6; ++i)
        want_cells.push_back(i);
    EXPECT_EQ(region_entries(r.mMesh, "all", RegionKind::Cell), want_cells);

    // The side survives (v16.27.0) as the child facet that IS the parent's
    // face: the same four nodes.
    const std::vector<std::int64_t> bottom = region_entries(r.mMesh, "bottom", RegionKind::Side);
    ASSERT_EQ(bottom.size(), 2u);
    meshioplusplus::CellType facet_type{};
    std::vector<std::int64_t> got_nodes, want_nodes;
    ASSERT_TRUE(
        meshioplusplus::detail::facet_nodes(r.mMesh, bottom[0], bottom[1], facet_type, got_nodes));
    ASSERT_TRUE(meshioplusplus::detail::facet_nodes(m, 0, 0, facet_type, want_nodes));
    std::sort(got_nodes.begin(), got_nodes.end());
    std::sort(want_nodes.begin(), want_nodes.end());
    EXPECT_EQ(got_nodes, want_nodes);
}

// --------------------------------------------------------------------------
// Winding failure
// --------------------------------------------------------------------------

TEST(Subdivide, AnOpenFaceSetThrowsRatherThanGuessing) {
    // Drop the top face of a cube: the faces no longer close, so every edge
    // of the missing face is used exactly once instead of twice -- exactly
    // the `Unorientable` signature `orient_rings` reports.
    Mesh m;
    m.AssignPoints(mt::points_from(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}));
    m.AddPolyhedronBlock("polyhedron8", {{
                                            {0, 3, 2, 1},
                                            {0, 1, 5, 4},
                                            {2, 3, 7, 6},
                                            {0, 4, 7, 3},
                                            {1, 2, 6, 5},
                                            // top face {4, 5, 6, 7} deliberately omitted
                                        }});
    EXPECT_THROW(subdivide(m), std::invalid_argument);
}
