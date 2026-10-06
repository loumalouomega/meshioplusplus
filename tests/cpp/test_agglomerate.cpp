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
#include "meshioplusplus/detail/facet_index.hpp"
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
// The exact ragged layout -- a regression pin for the flat output store
// --------------------------------------------------------------------------

namespace {

// An nx x ny x nz grid of unit hexahedra, x fastest.
Mesh hex_grid(std::size_t nx, std::size_t ny, std::size_t nz) {
    std::vector<std::vector<double>> pts;
    const auto id = [&](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<std::int64_t>((k * (ny + 1) + j) * (nx + 1) + i);
    };
    for (std::size_t k = 0; k <= nz; ++k)
        for (std::size_t j = 0; j <= ny; ++j)
            for (std::size_t i = 0; i <= nx; ++i)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    std::vector<std::vector<std::int64_t>> cells;
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i)
                cells.push_back({id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k),
                                 id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1),
                                 id(i, j + 1, k + 1)});
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddCellBlock("hexahedron", mt::conn_from(cells));
    return m;
}

// FNV-1a over a polyhedron block read back through `CellView`: the cell count,
// then per cell its face count, each face's size and its node ids in order.
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

TEST(Agglomerate, MergedLayoutIsPinnedWithAndWithoutCoplanarFusion) {
    // Digests taken from the implementation that built a nested vector per
    // external face, before the merged block moved to one flat node list with
    // offsets: the two must agree on every cell, face and node, in order,
    // including the reversed winding of faces the group sees from its far side
    // and the rings of fused coplanar patches.
    struct Case {
        bool mMerge;
        std::size_t mTarget;
        std::uint64_t mDigest;
    };
    const Case cases[] = {
        {false, 2, 15994731423313902787ull}, {false, 4, 10708260906522639862ull},
        {false, 12, 1369574336007378820ull}, {true, 2, 1323374833727562619ull},
        {true, 4, 18196367985544503846ull},  {true, 12, 18108531758129993186ull},
    };
    const Mesh grid = hex_grid(3, 2, 2);
    for (const Case& c : cases) {
        AgglomerateOptions o;
        o.mTargetGroupSize = c.mTarget;
        o.mMergeCoplanarFaces = c.mMerge;
        const AgglomerateResult r = agglomerate(grid, o);
        ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
        EXPECT_EQ(polyhedron_layout_digest(r.mMesh.Cells(0)), c.mDigest)
            << "merge=" << c.mMerge << " target=" << c.mTarget;
    }
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

TEST(Agglomerate, RegionsSurviveThroughTheGlobalMap) {
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

    // The side survives (v16.27.0): the merged polyhedron's face containing
    // hexahedron 0's bottom face, found by its nodes.
    const std::size_t bottom_idx = r.mMesh.FindRegion("bottom", RegionKind::Side);
    ASSERT_NE(bottom_idx, Mesh::npos);
    const Region& bottom = r.mMesh.Region(bottom_idx);
    ASSERT_EQ(bottom.NumEntries(), 1u);
    meshioplusplus::CellType facet_type{};
    std::vector<std::int64_t> want;
    ASSERT_TRUE(meshioplusplus::detail::facet_nodes(m, 0, 0, facet_type, want));
    const auto face = r.mMesh.Cells(0).Face(
        0, static_cast<std::size_t>(meshioplusplus::detail::read_int(bottom.mEntries, 1)));
    std::vector<std::int64_t> got_nodes(face.first, face.first + face.second);
    for (std::int64_t v : want)
        EXPECT_NE(std::find(got_nodes.begin(), got_nodes.end(), v), got_nodes.end()) << v;
}
