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
// Tests for polyhedral coarsening (`operations/agglomerate.hpp`).
//
// Three properties do the real work here, in order: the compact<->global
// bridge (`GlobalFaces::mCellToGlobal`), previously exercised by NO existing
// caller of `build_global_faces` in any direction, round-trips cell_data and
// regions correctly on a mesh where the two numberings genuinely diverge
// (a 2D block before the volume block); a real merge conserves volume
// EXACTLY (an identity of surviving boundary faces, not a divergence-theorem
// coincidence); and a 3-cell chain proves the internal/external face filter
// distinguishes "shared with a group-mate" from "shared with a different
// group" -- a distinction a 2-cell fixture cannot exercise at all.

// System includes
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::agglomerate;
using meshioplusplus::AgglomerateOptions;
using meshioplusplus::AgglomerateResult;
using meshioplusplus::compute_stats;
using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;

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

// Two unit hexahedra sharing exactly one face (x=1 plane), spanning x in
// [0,1] and [1,2] -- the same fixture `test_face_mesh.cpp`'s `two_hexes()`
// uses, duplicated here rather than shared (that helper is file-private).
Mesh two_hexes() {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {0, 1, 0},
                                    {0, 1, 1},
                                    {0, 0, 1},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {1, 1, 1},
                                    {1, 0, 1},
                                    {2, 0, 0},
                                    {2, 1, 0},
                                    {2, 1, 1},
                                    {2, 0, 1}}));
    m.AddCellBlock("hexahedron",
                   mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}, {4, 5, 6, 7, 8, 9, 10, 11}}));
    return m;
}

// Three unit hexahedra in a row, x in [0,1], [1,2], [2,3].
Mesh three_hexes() {
    Mesh m;
    std::vector<std::vector<double>> pts;
    for (int x = 0; x <= 3; ++x) {
        pts.push_back({static_cast<double>(x), 0, 0});
        pts.push_back({static_cast<double>(x), 1, 0});
        pts.push_back({static_cast<double>(x), 1, 1});
        pts.push_back({static_cast<double>(x), 0, 1});
    }
    m.AssignPoints(mt::points_from(pts));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7},
                                                {4, 5, 6, 7, 8, 9, 10, 11},
                                                {8, 9, 10, 11, 12, 13, 14, 15}}));
    return m;
}

// Whether merged cell `Cell` of block `Block` has a face containing exactly
// the given (unordered) node id set -- a winding-independent fingerprint.
bool has_face_with_nodes(const Mesh& rMesh, std::size_t Block, std::size_t Cell,
                         std::vector<std::int64_t> rWant) {
    std::sort(rWant.begin(), rWant.end());
    const auto cb = rMesh.Cells(Block);
    for (std::size_t f = 0; f < cb.NumFaces(Cell); ++f) {
        const auto face = cb.Face(Cell, f);
        std::vector<std::int64_t> got(face.first, face.first + face.second);
        std::sort(got.begin(), got.end());
        if (got == rWant)
            return true;
    }
    return false;
}

}  // namespace

// --------------------------------------------------------------------------
// The compact<->global bridge -- the central de-risking oracle
// --------------------------------------------------------------------------

TEST(Agglomerate, IdentityGroupingRoundTripsThroughTheCompactToGlobalBridge) {
    // A quad (non-volume) block BEFORE the hex block: compact and global cell
    // numbering genuinely diverge here (quad is global cell 0 but is not in
    // the compact space at all; the hexes are global cells 1 and 2 but
    // compact cells 0 and 1) -- the shape that actually exercises
    // mCellToGlobal, which no existing build_global_faces caller does.
    Mesh ordered;
    ordered.AssignPoints(mt::points_from({{0, 0, 0},
                                          {0, 1, 0},
                                          {0, 1, 1},
                                          {0, 0, 1},
                                          {1, 0, 0},
                                          {1, 1, 0},
                                          {1, 1, 1},
                                          {1, 0, 1},
                                          {2, 0, 0},
                                          {2, 1, 0},
                                          {2, 1, 1},
                                          {2, 0, 1}}));
    ordered.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}}));  // global cell 0
    ordered.AddCellBlock("hexahedron",
                         mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}, {4, 5, 6, 7, 8, 9, 10, 11}}));
    ordered.AddCellData("material", {i64({100}), i64({7, 8})});

    AgglomerateOptions opts;
    opts.mTargetGroupSize = 1;  // identity: every cell its own group
    const AgglomerateResult r = agglomerate(ordered, opts);

    ASSERT_EQ(r.mMesh.NumCellBlocks(), 2u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "quad");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(1).Type()), "polyhedron");
    EXPECT_EQ(r.mMesh.Cells(1).NumCells(), 2u);
    // Each hex kept its full 6-face boundary, including the shared face
    // (present in BOTH singleton "groups", exactly as two separate original
    // cells shared it).
    EXPECT_EQ(r.mMesh.Cells(1).NumFaces(0), 6u);
    EXPECT_EQ(r.mMesh.Cells(1).NumFaces(1), 6u);

    EXPECT_NEAR(compute_stats(r.mMesh).mSignedVolume, compute_stats(ordered).mSignedVolume, 1e-12);
    EXPECT_EQ(compute_stats(ordered).mSignedVolume, 2.0);

    // The flat cell map: quad(0)->0, hex0(1)->1, hex1(2)->2 -- a literal
    // identity, since the quad is the only pass-through block (base 0) and
    // the merged block immediately follows it (base 1).
    ASSERT_EQ(r.mCellMap.Size(), 3u);
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMap, 0), 0);
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMap, 1), 1);
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMap, 2), 2);

    // cell_data round-tripped through mCellToGlobal + block/row resolution.
    ASSERT_TRUE(r.mMesh.HasCellData("material"));
    ASSERT_EQ(r.mMesh.CellDataNumBlocks("material"), 2u);
    EXPECT_EQ(r.mMesh.CellData("material", 0).As<std::int64_t>()[0], 100);
    const NDArray& merged_mat = r.mMesh.CellData("material", 1);
    EXPECT_EQ(merged_mat.As<std::int64_t>()[0], 7);
    EXPECT_EQ(merged_mat.As<std::int64_t>()[1], 8);
}

// --------------------------------------------------------------------------
// Exact volume conservation on a real (non-identity) merge
// --------------------------------------------------------------------------

TEST(Agglomerate, TwoAdjacentHexesMergeIntoOnePolyhedronConservingVolumeExactly) {
    const Mesh m = two_hexes();
    AgglomerateOptions opts;
    opts.mTargetGroupSize = 2;
    const AgglomerateResult r = agglomerate(m, opts);

    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "polyhedron");
    ASSERT_EQ(r.mMesh.Cells(0).NumCells(), 1u);
    // 12 total face-references (6 per hex) minus the 2 references to the 1
    // shared, now-internal face.
    EXPECT_EQ(r.mMesh.Cells(0).NumFaces(0), 10u);

    // An identity of surviving boundary faces, not a divergence-theorem
    // coincidence -- exact, not a tolerance.
    EXPECT_EQ(compute_stats(r.mMesh).mSignedVolume, 2.0);
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);

    ASSERT_EQ(r.mCellMap.Size(), 2u);
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMap, 0),
              meshioplusplus::detail::read_int(r.mCellMap, 1))
        << "both hexes must land in the same merged cell";
}

// --------------------------------------------------------------------------
// The internal-vs-external filter, the case two cells cannot exercise
// --------------------------------------------------------------------------

TEST(Agglomerate, AThreeHexChainKeepsTheCrossGroupFaceAndDropsTheInternalOne) {
    const Mesh m = three_hexes();
    AgglomerateOptions opts;
    opts.mTargetGroupSize = 2;  // greedy growth: {hex0,hex1} then {hex2} alone
    const AgglomerateResult r = agglomerate(m, opts);

    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    ASSERT_EQ(r.mMesh.Cells(0).NumCells(), 2u);

    // group0 = {hex0, hex1}: the hex0/hex1 shared face (nodes 4,5,6,7) is
    // internal and dropped from BOTH cells; the hex1/hex2 shared face (nodes
    // 8,9,10,11) crosses into the OTHER group and survives on both sides.
    EXPECT_FALSE(has_face_with_nodes(r.mMesh, 0, 0, {4, 5, 6, 7}));
    EXPECT_FALSE(has_face_with_nodes(r.mMesh, 0, 1, {4, 5, 6, 7}));
    EXPECT_TRUE(has_face_with_nodes(r.mMesh, 0, 0, {8, 9, 10, 11}));
    EXPECT_TRUE(has_face_with_nodes(r.mMesh, 0, 1, {8, 9, 10, 11}));

    EXPECT_EQ(r.mMesh.Cells(0).NumFaces(0), 10u);  // {hex0,hex1}: as the 2-hex test
    EXPECT_EQ(r.mMesh.Cells(0).NumFaces(1), 6u);   // {hex2} alone: its full boundary

    EXPECT_EQ(compute_stats(r.mMesh).mSignedVolume, 3.0);

    ASSERT_EQ(r.mCellMap.Size(), 3u);
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMap, 0),
              meshioplusplus::detail::read_int(r.mCellMap, 1));
    EXPECT_NE(meshioplusplus::detail::read_int(r.mCellMap, 1),
              meshioplusplus::detail::read_int(r.mCellMap, 2));
}

// --------------------------------------------------------------------------
// Options and errors
// --------------------------------------------------------------------------

TEST(Agglomerate, ZeroTargetGroupSizeThrows) {
    AgglomerateOptions opts;
    opts.mTargetGroupSize = 0;
    EXPECT_THROW(agglomerate(two_hexes(), opts), std::invalid_argument);
}

TEST(Agglomerate, ANonManifoldFaceIsRefusedByName) {
    // Three tetrahedra sharing one triangular face -- that face is used by
    // three cells, which orient_rings/build_global_faces reports as
    // non-manifold rather than guessing an owner/neighbour pairing for it.
    Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}, {1, 1, 1}}));
    m.AddCellBlock("tetra", mt::conn_from({{0, 1, 2, 3}, {0, 2, 1, 4}, {0, 1, 2, 5}}));
    EXPECT_THROW(agglomerate(m), std::invalid_argument);
}

TEST(Agglomerate, NonVolumeBlocksPassThroughUnchangedWhenNoVolumeCellsExist) {
    Mesh m = mt::tri_mesh();  // 2D only, no volume cells at all
    const AgglomerateResult r = agglomerate(m);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "triangle");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), m.Cells(0).NumCells());
    ASSERT_EQ(r.mCellMap.Size(), m.Cells(0).NumCells());
    for (std::size_t i = 0; i < r.mCellMap.Size(); ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(r.mCellMap, i), static_cast<std::int64_t>(i));
}

// --------------------------------------------------------------------------
// Regions
// --------------------------------------------------------------------------

TEST(Agglomerate, PointAndCellRegionsSurviveThroughTheGlobalMap) {
    Mesh m = two_hexes();
    m.AddRegion(Region("corner", RegionKind::Point, i64({0})));
    m.AddRegion(Region("both", RegionKind::Cell, i64({0, 1})));
    m.AddRegion(Region("bottom", RegionKind::Side, i64_pairs({0, 0})));

    AgglomerateOptions opts;
    opts.mTargetGroupSize = 2;
    const AgglomerateResult r = agglomerate(m, opts);

    // Points are never renumbered: the region entry is untouched.
    ASSERT_NE(r.mMesh.FindRegion("corner", RegionKind::Point), Mesh::npos);
    const Region& corner = r.mMesh.Region(r.mMesh.FindRegion("corner", RegionKind::Point));
    EXPECT_EQ(meshioplusplus::detail::read_int(corner.mEntries, 0), 0);

    // Both original cells collapsed into the one merged cell (global 0).
    ASSERT_NE(r.mMesh.FindRegion("both", RegionKind::Cell), Mesh::npos);
    const Region& both = r.mMesh.Region(r.mMesh.FindRegion("both", RegionKind::Cell));
    ASSERT_EQ(both.NumEntries(), 1u) << "duplicate global-0 entries collapse via Canonicalize";
    EXPECT_EQ(meshioplusplus::detail::read_int(both.mEntries, 0), 0);

    // A many-to-one collapse can never preserve facet identity -- the merged
    // cell's type ("polyhedron") never matches the original ("hexahedron"),
    // so `remap_region`'s Side branch drops every entry. The region is still
    // *carried*, as a named empty group (the same "the name is information"
    // convention `test_region_remap.cpp`'s
    // `ARegionLeftEmptyIsStillCarriedAsAnEmptyGroup` pins) -- not removed.
    const std::size_t bottom_idx = r.mMesh.FindRegion("bottom", RegionKind::Side);
    ASSERT_NE(bottom_idx, Mesh::npos);
    EXPECT_EQ(r.mMesh.Region(bottom_idx).NumEntries(), 0u);
}
