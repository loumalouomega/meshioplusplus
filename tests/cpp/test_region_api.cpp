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
 * @file test_region_api.cpp
 * @brief Backend-agnostic tests of the named-region part of the uniform mesh
 * API (`region.hpp` + `mesh_api.hpp`).
 *
 * Every assertion here must hold for ALL mesh backends (MESHIO, NATIVE,
 * KRATOS) — this file is to regions what `test_mesh_api.cpp` is to points and
 * cells, and it runs under each backend leg of CI. Backend-specific behaviour
 * (the KRATOS SubModelPart materialization) lives at the bottom behind the
 * backend macro, mirroring `test_kratos_backend.cpp`'s split.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;

namespace {

/// An Int64 `(n,)` entry array from a plain list.
NDArray entries(const std::vector<std::int64_t>& rVals) {
    NDArray a = NDArray::Uninit(DType::Int64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

/// An Int64 `(n,2)` side-entry array from flat (cell, facet) pairs.
NDArray side_entries(const std::vector<std::int64_t>& rFlat) {
    NDArray a = NDArray::Uninit(DType::Int64, {rFlat.size() / 2, 2});
    for (std::size_t i = 0; i < rFlat.size(); ++i)
        a.As<std::int64_t>()[i] = rFlat[i];
    return a;
}

/// The entries of region @p i, flattened, for comparison against a literal.
std::vector<std::int64_t> flat(const Mesh& rMesh, std::size_t i) {
    const Region& r = rMesh.Region(i);
    std::vector<std::int64_t> out(r.mEntries.Size());
    for (std::size_t k = 0; k < out.size(); ++k)
        out[k] = meshioplusplus::detail::read_int(r.mEntries, k);
    return out;
}

}  // namespace

// --------------------------------------------------------------------------
// The kind vocabulary
// --------------------------------------------------------------------------

TEST(RegionApi, KindNamesRoundTrip) {
    for (RegionKind k : {RegionKind::Point, RegionKind::Cell, RegionKind::Side}) {
        const std::string name = meshioplusplus::region_kind_name(k);
        EXPECT_EQ(meshioplusplus::region_kind_from_name(name), k);
    }
    EXPECT_EQ(std::string(meshioplusplus::region_kind_name(RegionKind::Point)), "point");
    EXPECT_EQ(std::string(meshioplusplus::region_kind_name(RegionKind::Cell)), "cell");
    EXPECT_EQ(std::string(meshioplusplus::region_kind_name(RegionKind::Side)), "side");
    EXPECT_THROW(meshioplusplus::region_kind_from_name("facet"), std::invalid_argument);
}

// --------------------------------------------------------------------------
// Canonicalization: sorted + de-duplicated, on every backend
// --------------------------------------------------------------------------

TEST(RegionApi, EntriesAreSortedAndDeduplicated) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("corners", RegionKind::Point, entries({2, 0, 2, 1, 0})));
    ASSERT_EQ(m.NumRegions(), 1u);
    EXPECT_EQ(flat(m, 0), (std::vector<std::int64_t>{0, 1, 2}));
    EXPECT_EQ(m.Region(0).NumEntries(), 3u);
    EXPECT_EQ(m.Region(0).Stride(), 1u);
}

TEST(RegionApi, SideEntriesSortLexicographicallyOnThePair) {
    Mesh m = mt::tet_mesh();
    // Deliberately unsorted, with one exact duplicate pair.
    m.AddRegion(Region("wall", RegionKind::Side, side_entries({1, 2, 0, 3, 1, 0, 0, 3})));
    ASSERT_EQ(m.NumRegions(), 1u);
    EXPECT_EQ(flat(m, 0), (std::vector<std::int64_t>{0, 3, 1, 0, 1, 2}));
    EXPECT_EQ(m.Region(0).NumEntries(), 3u);
    EXPECT_EQ(m.Region(0).Stride(), 2u);
    ASSERT_EQ(m.Region(0).mEntries.Ndim(), 2u);
    EXPECT_EQ(m.Region(0).mEntries.Shape()[1], 2u);
}

TEST(RegionApi, CanonicalizeIsIdempotentAndWidensNarrowDtypes) {
    NDArray narrow = NDArray::Uninit(DType::Int32, {4});
    narrow.As<std::int32_t>()[0] = 7;
    narrow.As<std::int32_t>()[1] = 3;
    narrow.As<std::int32_t>()[2] = 7;
    narrow.As<std::int32_t>()[3] = 1;

    Region r("s", RegionKind::Cell, std::move(narrow));
    r.Canonicalize();
    EXPECT_EQ(r.mEntries.Dtype(), DType::Int64);
    ASSERT_EQ(r.NumEntries(), 3u);
    EXPECT_EQ(r.Entries()[0], 1);
    EXPECT_EQ(r.Entries()[2], 7);

    const std::vector<std::int64_t> before(r.Entries(), r.Entries() + r.NumEntries());
    r.Canonicalize();
    const std::vector<std::int64_t> after(r.Entries(), r.Entries() + r.NumEntries());
    EXPECT_EQ(before, after);
}

TEST(RegionApi, EmptyRegionIsWellFormed) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("empty", RegionKind::Cell, entries({})));
    ASSERT_EQ(m.NumRegions(), 1u);
    EXPECT_EQ(m.Region(0).NumEntries(), 0u);
    EXPECT_EQ(m.Region(0).Entries(), nullptr);
}

// --------------------------------------------------------------------------
// Storage: insert-or-replace on the key, deterministic order, sorted names
// --------------------------------------------------------------------------

TEST(RegionApi, SameKeyReplacesRatherThanDuplicates) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("wall", RegionKind::Point, entries({0, 1})));
    m.AddRegion(Region("wall", RegionKind::Point, entries({2})));
    ASSERT_EQ(m.NumRegions(), 1u);
    EXPECT_EQ(flat(m, 0), (std::vector<std::int64_t>{2}));
}

TEST(RegionApi, NameAndKindTogetherFormTheKey) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("wall", RegionKind::Point, entries({0})));
    m.AddRegion(Region("wall", RegionKind::Cell, entries({0})));
    EXPECT_EQ(m.NumRegions(), 2u);
    // ...and so do dim and tag: gmsh physical groups of different dimensions
    // legitimately share a name.
    m.AddRegion(Region("wall", RegionKind::Cell, /*dim=*/2, /*tag=*/7, entries({0})));
    EXPECT_EQ(m.NumRegions(), 3u);
}

TEST(RegionApi, RegionOrderIsDeterministicRegardlessOfInsertionOrder) {
    Mesh a = mt::tri_mesh();
    a.AddRegion(Region("zeta", RegionKind::Cell, entries({0})));
    a.AddRegion(Region("alpha", RegionKind::Point, entries({1})));
    a.AddRegion(Region("beta", RegionKind::Side, side_entries({0, 0})));

    Mesh b = mt::tri_mesh();
    b.AddRegion(Region("beta", RegionKind::Side, side_entries({0, 0})));
    b.AddRegion(Region("zeta", RegionKind::Cell, entries({0})));
    b.AddRegion(Region("alpha", RegionKind::Point, entries({1})));

    ASSERT_EQ(a.NumRegions(), b.NumRegions());
    for (std::size_t i = 0; i < a.NumRegions(); ++i)
        EXPECT_TRUE(meshioplusplus::regions_equal(a.Region(i), b.Region(i)))
            << "region " << i << " differs";
    // Point < Cell < Side, then by name.
    EXPECT_EQ(a.Region(0).mName, "alpha");
    EXPECT_EQ(a.Region(1).mName, "zeta");
    EXPECT_EQ(a.Region(2).mName, "beta");
}

TEST(RegionApi, RegionNamesAreSortedAndDeduplicated) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("wall", RegionKind::Cell, entries({0})));
    m.AddRegion(Region("inlet", RegionKind::Point, entries({0})));
    m.AddRegion(Region("wall", RegionKind::Point, entries({1})));
    EXPECT_EQ(m.RegionNames(), (std::vector<std::string>{"inlet", "wall"}));
}

TEST(RegionApi, HasAndFind) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("wall", RegionKind::Cell, entries({0})));
    EXPECT_TRUE(m.HasRegion("wall"));
    EXPECT_TRUE(m.HasRegion("wall", RegionKind::Cell));
    EXPECT_FALSE(m.HasRegion("wall", RegionKind::Point));
    EXPECT_FALSE(m.HasRegion("floor"));
    EXPECT_EQ(m.FindRegion("wall", RegionKind::Cell), 0u);
    EXPECT_EQ(m.FindRegion("wall", RegionKind::Point), Mesh::npos);
}

TEST(RegionApi, DimAndTagAreCarried) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("Surface", RegionKind::Cell, /*dim=*/2, /*tag=*/12, entries({0, 1})));
    ASSERT_EQ(m.NumRegions(), 1u);
    EXPECT_EQ(m.Region(0).mDim, 2);
    EXPECT_EQ(m.Region(0).mTag, 12);
    // Defaults say "unspecified" rather than 0, so a real tag of 0 is expressible.
    Region plain("p", RegionKind::Point, entries({0}));
    EXPECT_EQ(plain.mDim, -1);
    EXPECT_EQ(plain.mTag, -1);
}

TEST(RegionApi, RegionsEqualComparesKeyAndEntries) {
    Region a("w", RegionKind::Cell, entries({1, 0}));
    Region b("w", RegionKind::Cell, entries({0, 1}));
    a.Canonicalize();
    b.Canonicalize();
    EXPECT_TRUE(meshioplusplus::regions_equal(a, b));

    Region c("w", RegionKind::Cell, entries({0, 2}));
    c.Canonicalize();
    EXPECT_FALSE(meshioplusplus::regions_equal(a, c));

    Region d("w", RegionKind::Point, entries({0, 1}));
    d.Canonicalize();
    EXPECT_FALSE(meshioplusplus::regions_equal(a, d));
}

TEST(RegionApi, MeshWithoutRegionsReportsNone) {
    const Mesh m = mt::tri_mesh();
    EXPECT_EQ(m.NumRegions(), 0u);
    EXPECT_TRUE(m.RegionNames().empty());
    EXPECT_FALSE(m.HasRegion("anything"));
}

// --------------------------------------------------------------------------
// The global (block-major) cell index the Cell/Side kinds are defined against
// --------------------------------------------------------------------------

TEST(RegionApi, BlockBasesCoverEveryBlockIncludingEmptyOnes) {
    Mesh m = mt::tri_quad_mesh();
    const auto bases = meshioplusplus::detail::block_bases(m);
    ASSERT_EQ(bases.size(), m.NumCellBlocks() + 1u);
    EXPECT_EQ(bases.front(), 0);

    std::int64_t total = 0;
    for (const auto cb : m.CellRange())
        total += static_cast<std::int64_t>(cb.NumCells());
    EXPECT_EQ(meshioplusplus::detail::total_cells(bases), total);
}

TEST(RegionApi, GlobalCellIndexRoundTrips) {
    Mesh m = mt::tri_quad_mesh();
    const auto bases = meshioplusplus::detail::block_bases(m);
    for (std::size_t b = 0; b + 1 < bases.size(); ++b) {
        for (std::int64_t row = 0; row < bases[b + 1] - bases[b]; ++row) {
            const std::int64_t g = meshioplusplus::detail::block_row_to_global(bases, b, row);
            const auto [gb, grow] = meshioplusplus::detail::global_to_block_row(bases, g);
            EXPECT_EQ(gb, b);
            EXPECT_EQ(grow, row);
        }
    }
    constexpr std::size_t npos = static_cast<std::size_t>(-1);
    EXPECT_EQ(meshioplusplus::detail::global_to_block_row(bases, -1).first, npos);
    EXPECT_EQ(meshioplusplus::detail::global_to_block_row(bases, bases.back()).first, npos);
}

// --------------------------------------------------------------------------
// KRATOS: regions materialize as SubModelParts and win their names
// --------------------------------------------------------------------------

#if defined(MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)

// NOTE: SubModelParts share the root's entity containers, so `Elements()` /
// `Nodes()` on one returns the *root's* container. The per-part membership is
// `NumberOfElements()` / `NumberOfNodes()`, which is what these assert.

TEST(RegionApiKratos, CellRegionBecomesASubModelPart) {
    Mesh m = mt::tri_mesh();  // 2 triangles over 4 points
    m.AddRegion(Region("Surface", RegionKind::Cell, /*dim=*/2, /*tag=*/1, entries({0})));

    meshioplusplus::ModelPart& r_mp = m.GetModelPart();
    ASSERT_TRUE(r_mp.HasSubModelPart("Surface"));
    const meshioplusplus::ModelPart& r_smp = r_mp.GetSubModelPart("Surface");
    EXPECT_EQ(r_smp.NumberOfElements(), 1u);
    EXPECT_EQ(r_smp.NumberOfConditions(), 0u);
    // Kratos convention: a sub model part carries its entities' nodes.
    EXPECT_EQ(r_smp.NumberOfNodes(), 3u);
}

TEST(RegionApiKratos, PointRegionBecomesANodeOnlySubModelPart) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("Clamped", RegionKind::Point, entries({0, 2})));

    meshioplusplus::ModelPart& r_mp = m.GetModelPart();
    ASSERT_TRUE(r_mp.HasSubModelPart("Clamped"));
    const meshioplusplus::ModelPart& r_smp = r_mp.GetSubModelPart("Clamped");
    EXPECT_EQ(r_smp.NumberOfNodes(), 2u);
    EXPECT_EQ(r_smp.NumberOfElements(), 0u);
    EXPECT_EQ(r_smp.NumberOfConditions(), 0u);
}

TEST(RegionApiKratos, SideRegionIsKeptButNotMaterialized) {
    Mesh m = mt::tet_mesh();
    m.AddRegion(Region("wall", RegionKind::Side, side_entries({0, 1})));

    meshioplusplus::ModelPart& r_mp = m.GetModelPart();
    EXPECT_FALSE(r_mp.HasSubModelPart("wall"));
    // The region itself survives — only the ModelPart projection is skipped.
    EXPECT_TRUE(m.HasRegion("wall", RegionKind::Side));
}

TEST(RegionApiKratos, ExplicitRegionWinsOverTheInferredTagSubModelPart) {
    Mesh m = mt::tri_mesh();  // 2 triangles
    // The tag pass would name this SubModelPart "gmsh_physical_5" and put both
    // elements in it; the explicit region claims the name first.
    NDArray tags = NDArray::Uninit(DType::Int64, {2});
    tags.As<std::int64_t>()[0] = 5;
    tags.As<std::int64_t>()[1] = 5;
    m.AddCellData("gmsh:physical", {std::move(tags)});
    m.AddRegion(Region("gmsh_physical_5", RegionKind::Point, entries({1})));

    meshioplusplus::ModelPart& r_mp = m.GetModelPart();
    ASSERT_TRUE(r_mp.HasSubModelPart("gmsh_physical_5"));
    const meshioplusplus::ModelPart& r_smp = r_mp.GetSubModelPart("gmsh_physical_5");
    // The region's node-only projection, not the tag pass's two elements.
    EXPECT_EQ(r_smp.NumberOfElements(), 0u);
    EXPECT_EQ(r_smp.NumberOfNodes(), 1u);
}

TEST(RegionApiKratos, RegionsSurviveAModelPartInvalidation) {
    Mesh m = mt::tri_mesh();
    m.AddRegion(Region("Surface", RegionKind::Cell, entries({0})));
    (void)m.GetModelPart();
    m.InvalidateBlocks();  // stage is rebuilt from the ModelPart on next use
    EXPECT_TRUE(m.HasRegion("Surface", RegionKind::Cell));
    EXPECT_EQ(m.NumRegions(), 1u);
}

#endif  // MESHIOPLUSPLUS_MESH_BACKEND_KRATOS
