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
/**
 * @file test_region_remap.cpp
 * @brief The shared region remapper (`detail/region_remap.hpp`) against
 * hand-written index maps.
 *
 * Driving it with synthetic maps rather than through a real operation is
 * deliberate: it pins the *contract* the operations rely on — what each
 * `CellMapKind` means, when a side entry survives — independently of any one
 * operation's behaviour. The operations' own suites then only have to check
 * that they pass the right maps.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;
using meshioplusplus::detail::CellMapKind;
using meshioplusplus::detail::RegionRemap;

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

/// The flattened entries of the output region named @p rName, or an empty list.
std::vector<std::int64_t> got(const Mesh& rMesh, const std::string& rName, RegionKind kind) {
    const std::size_t i = rMesh.FindRegion(rName, kind);
    if (i == Mesh::npos)
        return {};
    const Region& r = rMesh.Region(i);
    std::vector<std::int64_t> out(r.mEntries.Size());
    for (std::size_t k = 0; k < out.size(); ++k)
        out[k] = meshioplusplus::detail::read_int(r.mEntries, k);
    return out;
}

}  // namespace

// --------------------------------------------------------------------------
// Point regions
// --------------------------------------------------------------------------

TEST(RegionRemap, PointEntriesAreRemappedAndDroppedEntriesVanish) {
    Mesh in = mt::tri_mesh();  // 4 points
    Mesh out = mt::make_mesh({{0, 0, 0}, {1, 0, 0}}, "triangle", {{0, 1, 0}});
    in.AddRegion(Region("fixed", RegionKind::Point, i64({0, 2, 3})));

    // points 0 -> 1, 1 -> 0, 2 and 3 dropped.
    NDArray pm = i64({1, 0, -1, -1});
    RegionRemap maps;
    maps.pPointMap = &pm;
    maps.mOpName = "test";
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(got(out, "fixed", RegionKind::Point), (std::vector<std::int64_t>{1}));
}

TEST(RegionRemap, ARegionLeftEmptyIsStillCarriedAsAnEmptyGroup) {
    // The group's *name* is information: a crop that keeps none of "gone"'s
    // points has still been told that a "gone" exists. This is also what the
    // per-operation Python shims did before the shared helper replaced them.
    Mesh in = mt::tri_mesh();
    Mesh out = mt::tri_mesh();
    in.AddRegion(Region("gone", RegionKind::Point, i64({0, 1})));

    NDArray pm = i64({-1, -1, -1, -1});
    RegionRemap maps;
    maps.pPointMap = &pm;
    meshioplusplus::detail::remap_regions(in, out, maps);

    ASSERT_EQ(out.NumRegions(), 1u);
    EXPECT_EQ(out.Region(0).mName, "gone");
    EXPECT_EQ(out.Region(0).NumEntries(), 0u);
}

TEST(RegionRemap, ANullPointMapIsTheIdentity) {
    Mesh in = mt::tri_mesh();
    Mesh out = mt::tri_mesh();
    in.AddRegion(Region("fixed", RegionKind::Point, i64({1, 3})));

    RegionRemap maps;  // no maps at all
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(got(out, "fixed", RegionKind::Point), (std::vector<std::int64_t>{1, 3}));
}

TEST(RegionRemap, DimAndTagSurviveTheRemap) {
    Mesh in = mt::tri_mesh();
    Mesh out = mt::tri_mesh();
    in.AddRegion(Region("Surface", RegionKind::Cell, /*dim=*/2, /*tag=*/9, i64({0})));

    RegionRemap maps;
    meshioplusplus::detail::remap_regions(in, out, maps);

    const std::size_t i = out.FindRegion("Surface", RegionKind::Cell);
    ASSERT_NE(i, Mesh::npos);
    EXPECT_EQ(out.Region(i).mDim, 2);
    EXPECT_EQ(out.Region(i).mTag, 9);
}

// --------------------------------------------------------------------------
// Cell regions: the three map shapes
// --------------------------------------------------------------------------

TEST(RegionRemap, DirectMapKeepsSurvivorsAndDropsTheRest) {
    Mesh in = mt::tri_quad_mesh();  // triangle block then quad block
    Mesh out = mt::tri_quad_mesh();
    const std::size_t ntri = in.Cells(0).NumCells();
    const std::size_t nquad = in.Cells(1).NumCells();
    ASSERT_GE(ntri, 2u);
    ASSERT_GE(nquad, 1u);

    // Global numbering: triangles [0, ntri), quads [ntri, ntri + nquad).
    in.AddRegion(Region("mixed", RegionKind::Cell, i64({0, static_cast<std::int64_t>(ntri)})));

    std::vector<NDArray> cell_maps;
    std::vector<std::int64_t> tri_map(ntri, -1), quad_map(nquad, -1);
    tri_map[0] = 1;   // triangle 0 -> output triangle 1
    quad_map[0] = 0;  // quad 0 -> output quad 0
    cell_maps.push_back(i64(tri_map));
    cell_maps.push_back(i64(quad_map));

    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::Direct;
    maps.pCellMaps = &cell_maps;
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(got(out, "mixed", RegionKind::Cell),
              (std::vector<std::int64_t>{1, static_cast<std::int64_t>(ntri)}));
}

// Every current Direct-map caller (crop, split, clean, partition, reorder) is
// injective by construction -- confirmed by reading each one -- so a
// genuinely non-injective map (several input cells collapsing to the SAME
// output index, the shape undo_green needs) is an untested code path until
// this test. The mechanism needs no change: rremap_children's Direct branch
// is a plain per-entry lookup with no uniqueness assumption anywhere, and
// duplicate entries are collapsed by Region::Canonicalize's sort+unique,
// which RegionList::Add runs unconditionally on every stored region.
TEST(RegionRemap, ANonInjectiveDirectMapCollapsesDuplicateRegionEntries) {
    Mesh in = mt::tri_mesh();
    Mesh out = mt::tri_mesh();
    const std::size_t ntri = in.Cells(0).NumCells();
    ASSERT_GE(ntri, 2u);

    in.AddRegion(Region("pair", RegionKind::Cell, i64({0, 1})));

    std::vector<NDArray> cell_maps;
    std::vector<std::int64_t> tri_map(ntri, -1);
    tri_map[0] = 0;  // both cells 0 and 1 collapse onto the SAME output cell 0
    tri_map[1] = 0;
    cell_maps.push_back(i64(tri_map));

    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::Direct;
    maps.pCellMaps = &cell_maps;
    meshioplusplus::detail::remap_regions(in, out, maps);

    // Two entries collapsing onto the same output index dedup to one.
    EXPECT_EQ(got(out, "pair", RegionKind::Cell), (std::vector<std::int64_t>{0}));
}

TEST(RegionRemap, FirstChildMapExpandsAParentToItsChildren) {
    // One input triangle block of 2 cells; each cell becomes 3 output cells.
    Mesh in = mt::tri_mesh();
    Mesh out = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle",
                             {{0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 2, 3}, {0, 2, 3}, {0, 2, 3}});
    in.AddRegion(Region("left", RegionKind::Cell, i64({0})));

    std::vector<NDArray> cell_maps{i64({0, 3})};
    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::FirstChild;
    maps.pCellMaps = &cell_maps;
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(got(out, "left", RegionKind::Cell), (std::vector<std::int64_t>{0, 1, 2}));
}

TEST(RegionRemap, FirstChildRunEndsAtTheOutputBlockWhenTheParentIsLast) {
    Mesh in = mt::tri_mesh();
    Mesh out = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle",
                             {{0, 1, 2}, {0, 1, 2}, {0, 2, 3}, {0, 2, 3}});
    in.AddRegion(Region("right", RegionKind::Cell, i64({1})));

    std::vector<NDArray> cell_maps{i64({0, 2})};
    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::FirstChild;
    maps.pCellMaps = &cell_maps;
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(got(out, "right", RegionKind::Cell), (std::vector<std::int64_t>{2, 3}));
}

TEST(RegionRemap, FirstChildSkipsOverAFullyDroppedParent) {
    // Parent 1 collapsed entirely (map -1); parent 0's run must still stop at
    // parent 2's first child rather than swallowing it.
    Mesh in = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle",
                            {{0, 1, 2}, {0, 1, 2}, {0, 2, 3}});
    Mesh out = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle",
                             {{0, 1, 2}, {0, 1, 2}, {0, 2, 3}});
    in.AddRegion(Region("a", RegionKind::Cell, i64({0})));
    in.AddRegion(Region("b", RegionKind::Cell, i64({1})));

    std::vector<NDArray> cell_maps{i64({0, -1, 2})};
    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::FirstChild;
    maps.pCellMaps = &cell_maps;
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(got(out, "a", RegionKind::Cell), (std::vector<std::int64_t>{0, 1}));
    EXPECT_TRUE(got(out, "b", RegionKind::Cell).empty());  // carried, but with no entries
}

TEST(RegionRemap, GlobalMapIsIndexedByInputGlobalCell) {
    Mesh in = mt::tri_quad_mesh();
    Mesh out = mt::tri_quad_mesh();
    const auto ntri = static_cast<std::int64_t>(in.Cells(0).NumCells());
    in.AddRegion(Region("g", RegionKind::Cell, i64({0, ntri})));

    const std::int64_t ncells = ntri + static_cast<std::int64_t>(in.Cells(1).NumCells());
    std::vector<std::int64_t> flat(static_cast<std::size_t>(ncells), -1);
    flat[0] = ncells - 1;                       // first cell -> last output cell
    flat[static_cast<std::size_t>(ntri)] = -1;  // the quad was dropped
    NDArray gmap = i64(flat);

    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::Global;
    maps.pGlobalCellMap = &gmap;
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(got(out, "g", RegionKind::Cell), (std::vector<std::int64_t>{ncells - 1}));
}

TEST(RegionRemap, BlockMapRedirectsToADifferentOutputBlock) {
    // split(drop_empty_blocks) drops the triangle block, so input block 1
    // becomes output block 0.
    Mesh in = mt::tri_quad_mesh();
    Mesh out = mt::quad_mesh();
    const auto ntri = static_cast<std::int64_t>(in.Cells(0).NumCells());
    in.AddRegion(Region("q", RegionKind::Cell, i64({ntri})));

    std::vector<NDArray> cell_maps{
        i64(std::vector<std::int64_t>(static_cast<std::size_t>(ntri), -1)), i64({0, -1})};
    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::Direct;
    maps.pCellMaps = &cell_maps;
    maps.mBlockMap = {meshioplusplus::detail::kBlockDropped, 0};
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(got(out, "q", RegionKind::Cell), (std::vector<std::int64_t>{0}));
}

// --------------------------------------------------------------------------
// Side regions
// --------------------------------------------------------------------------

TEST(RegionRemap, SideEntrySurvivesWhenItsCellAndFacetDo) {
    Mesh in = mt::tet_mesh();  // 2 tetra, 4 faces each
    Mesh out = mt::tet_mesh();
    in.AddRegion(Region("wall", RegionKind::Side, i64_pairs({0, 2, 1, 0})));

    std::vector<NDArray> cell_maps{i64({1, 0})};  // swap the two tets
    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::Direct;
    maps.pCellMaps = &cell_maps;
    meshioplusplus::detail::remap_regions(in, out, maps);

    // (0,2) -> (1,2) and (1,0) -> (0,0), then canonically sorted.
    EXPECT_EQ(got(out, "wall", RegionKind::Side), (std::vector<std::int64_t>{0, 0, 1, 2}));
}

TEST(RegionRemap, SideEntryVanishesWithItsCell) {
    Mesh in = mt::tet_mesh();
    Mesh out = mt::tet_mesh();
    in.AddRegion(Region("wall", RegionKind::Side, i64_pairs({0, 1, 1, 3})));

    std::vector<NDArray> cell_maps{i64({-1, 0})};  // tet 0 dropped
    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::Direct;
    maps.pCellMaps = &cell_maps;
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(got(out, "wall", RegionKind::Side), (std::vector<std::int64_t>{0, 3}));
}

TEST(RegionRemap, AnOutOfRangeFacetIsDropped) {
    Mesh in = mt::tet_mesh();
    Mesh out = mt::tet_mesh();
    // A tetra has 4 faces, so facet 7 cannot exist.
    in.AddRegion(Region("bogus", RegionKind::Side, i64_pairs({0, 7})));

    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::Direct;
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_TRUE(got(out, "bogus", RegionKind::Side).empty());
}

TEST(RegionRemap, SideEntryVanishesWhenTheCellTypeChanges) {
    // A hexahedron's facets have no counterpart on a tetra, even 1:1.
    Mesh in = mt::hex_mesh();
    Mesh out = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 1}}, "tetra", {{0, 1, 2, 3}});
    in.AddRegion(Region("wall", RegionKind::Side, i64_pairs({0, 1})));

    std::vector<NDArray> cell_maps{i64({0})};
    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::Direct;
    maps.pCellMaps = &cell_maps;
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_TRUE(got(out, "wall", RegionKind::Side).empty());
}

TEST(RegionRemap, SideRegionsAreDroppedUnderAnExpandingMap) {
    Mesh in = mt::tet_mesh();
    Mesh out = mt::tet_mesh();
    in.AddRegion(Region("wall", RegionKind::Side, i64_pairs({0, 1})));
    in.AddRegion(Region("keep", RegionKind::Cell, i64({0})));

    std::vector<NDArray> cell_maps{i64({0, 1})};
    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::FirstChild;
    maps.pCellMaps = &cell_maps;
    maps.mOpName = "refine";
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_TRUE(got(out, "wall", RegionKind::Side).empty());
    // ...but the cell region beside it is carried normally.
    EXPECT_EQ(got(out, "keep", RegionKind::Cell), (std::vector<std::int64_t>{0}));
}

TEST(RegionRemap, DropSideRegionsForcesTheDrop) {
    Mesh in = mt::tet_mesh();
    Mesh out = mt::tet_mesh();
    in.AddRegion(Region("wall", RegionKind::Side, i64_pairs({0, 1})));

    RegionRemap maps;
    maps.mCellMapKind = CellMapKind::Direct;
    maps.mDropSideRegions = true;
    maps.mOpName = "test";
    meshioplusplus::detail::remap_regions(in, out, maps);

    EXPECT_EQ(out.NumRegions(), 0u);
}

// --------------------------------------------------------------------------
// Facet counts
// --------------------------------------------------------------------------

TEST(RegionRemap, FacetCountsComeFromTheFaceAndEdgeTables) {
    using meshioplusplus::detail::region_num_facets;
    EXPECT_EQ(region_num_facets("tetra"), 4u);
    EXPECT_EQ(region_num_facets("hexahedron"), 6u);
    EXPECT_EQ(region_num_facets("wedge"), 5u);
    EXPECT_EQ(region_num_facets("triangle"), 3u);
    EXPECT_EQ(region_num_facets("quad"), 4u);
    // A type with no facet table at all: every side entry on it is invalid.
    EXPECT_EQ(region_num_facets("line"), 0u);
}
