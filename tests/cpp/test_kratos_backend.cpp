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
 * @file test_kratos_backend.cpp
 * @brief KRATOS-backend-specific tests: lazy ModelPart materialization, the
 * Elements/Conditions split, automatic tags -> SubModelParts, ragged
 * pass-through, and the mutate + `InvalidateBlocks()` rebuild path.
 *
 * The whole file compiles away under other backends.
 */

#ifdef MESHIOPLUSPLUS_MESH_BACKEND_KRATOS

// System includes
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/mesh.hpp"

using meshioplusplus::CellType;
using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::ModelPart;
using meshioplusplus::NDArray;

namespace {

NDArray int64_array(const std::vector<std::int64_t>& rVals) {
    NDArray a = NDArray::Uninit(DType::Int64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

// A tri (2D) + tet (3D) mesh: tets become Elements, triangles Conditions.
Mesh mixed_dim_mesh() {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 0.5}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}, {0, 2, 3}}));
    m.AddCellBlock("tetra", mt::conn_from({{0, 1, 2, 4}, {0, 2, 3, 4}}));
    return m;
}

}  // namespace

TEST(KratosBackend, BackendName) {
    EXPECT_STREQ(meshioplusplus::mesh_backend_name(), "kratos");
}

TEST(KratosBackend, LazyMaterialization) {
    Mesh m = mt::tet_mesh();
    EXPECT_FALSE(m.IsMaterialized());
    EXPECT_EQ(m.NumPoints(), 5u);  // accessors do not materialize
    EXPECT_FALSE(m.IsMaterialized());
    ModelPart& r_mp = m.GetModelPart();
    EXPECT_TRUE(m.IsMaterialized());
    EXPECT_EQ(r_mp.NumberOfNodes(), 5u);
    EXPECT_EQ(r_mp.NumberOfElements(), 2u);
    EXPECT_EQ(r_mp.NumberOfConditions(), 0u);
    // Node ids are index + 1 with the staged coordinates.
    EXPECT_DOUBLE_EQ(r_mp.GetNode(2).X(), 1.0);
    EXPECT_EQ(r_mp.GetElement(1).Type(), CellType::Tetra);
    EXPECT_EQ(r_mp.GetElement(2).NodeIds(), (std::vector<meshioplusplus::IndexType>{1, 3, 4, 5}));
}

TEST(KratosBackend, ElementsConditionsSplitByDimension) {
    Mesh m = mixed_dim_mesh();
    ModelPart& r_mp = m.GetModelPart();
    EXPECT_EQ(r_mp.NumberOfElements(), 2u);    // the tets (dim 3 == mesh dim)
    EXPECT_EQ(r_mp.NumberOfConditions(), 2u);  // the triangles (dim 2)
    EXPECT_EQ(r_mp.GetCondition(1).Type(), CellType::Triangle);
    EXPECT_EQ(r_mp.GetElement(1).Type(), CellType::Tetra);
    // Writer path is unaffected: blocks and order are preserved.
    EXPECT_EQ(m.NumCellBlocks(), 2u);
    EXPECT_EQ(m.Cells(0).Type(), "triangle");
    EXPECT_EQ(m.Cells(1).Type(), "tetra");
}

TEST(KratosBackend, PointAndCellDataBecomeVariables) {
    Mesh m = mixed_dim_mesh();
    m.AddPointData("temp", mt::points_from({{1}, {2}, {3}, {4}, {5}}));
    m.AddCellData("gmsh:physical", {int64_array({10, 20}), int64_array({7, 7})});
    ModelPart& r_mp = m.GetModelPart();
    ASSERT_TRUE(r_mp.HasNodalData("temp"));
    EXPECT_DOUBLE_EQ(r_mp.GetNodalValue("temp", 3), 3.0);
    // Split per kind: triangle block rows -> conditional, tet block -> elemental.
    ASSERT_TRUE(r_mp.HasElementalData("gmsh:physical"));
    ASSERT_TRUE(r_mp.HasConditionalData("gmsh:physical"));
    EXPECT_EQ(r_mp.GetElementalData("gmsh:physical").As<std::int64_t>()[0], 7);
    EXPECT_EQ(r_mp.GetConditionalData("gmsh:physical").As<std::int64_t>()[1], 20);
    // The tag array is still cell_data (round-trip fidelity).
    EXPECT_TRUE(m.HasCellData("gmsh:physical"));
}

TEST(KratosBackend, TagsBecomeSubModelParts) {
    Mesh m = mixed_dim_mesh();
    m.AddCellData("gmsh:physical", {int64_array({10, 20}), int64_array({7, 7})});
    ModelPart& r_mp = m.GetModelPart();
    ASSERT_TRUE(r_mp.HasSubModelPart("gmsh_physical_10"));
    ASSERT_TRUE(r_mp.HasSubModelPart("gmsh_physical_20"));
    ASSERT_TRUE(r_mp.HasSubModelPart("gmsh_physical_7"));
    const ModelPart& r_smp7 = r_mp.GetSubModelPart("gmsh_physical_7");
    EXPECT_EQ(r_smp7.NumberOfElements(), 2u);
    EXPECT_EQ(r_smp7.NumberOfConditions(), 0u);
    EXPECT_EQ(r_smp7.NumberOfNodes(), 5u);  // union of the two tets' nodes
    const ModelPart& r_smp10 = r_mp.GetSubModelPart("gmsh_physical_10");
    EXPECT_EQ(r_smp10.NumberOfConditions(), 1u);
    EXPECT_EQ(r_smp10.NumberOfNodes(), 3u);
    EXPECT_TRUE(r_smp10.HasCondition(1));
}

TEST(KratosBackend, TagsToSubModelPartsCanBeDisabled) {
    Mesh m = mixed_dim_mesh();
    m.AddCellData("gmsh:physical", {int64_array({10, 20}), int64_array({7, 7})});
    m.SetBuildSubModelPartsFromTags(false);
    ModelPart& r_mp = m.GetModelPart();
    EXPECT_EQ(r_mp.NumberOfSubModelParts(), 0u);
    // The data itself still lands as elemental/conditional variables.
    EXPECT_TRUE(r_mp.HasElementalData("gmsh:physical"));
}

TEST(KratosBackend, RaggedBlocksPassThrough) {
    Mesh m;
    m.AssignPoints(NDArray(DType::Float64, {6, 3}));
    m.AddPolygonBlock("polygon", {{0, 1, 2}, {1, 3, 4, 2, 5}});
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    ModelPart& r_mp = m.GetModelPart();
    // The polygon block creates no entities but stays readable for writers.
    EXPECT_EQ(r_mp.NumberOfElements(), 1u);
    ASSERT_EQ(m.NumCellBlocks(), 2u);
    EXPECT_TRUE(m.Cells(0).IsRagged());
    EXPECT_EQ(m.Cells(0).RowSize(1), 5u);
}

TEST(KratosBackend, InvalidateBlocksRebuildsFromModelPart) {
    Mesh m = mt::tet_mesh();
    ModelPart& r_mp = m.GetModelPart();
    // Mutate the ModelPart directly: add a node and a triangle condition.
    r_mp.CreateNewNode(6, 2.0, 0.0, 0.0);
    r_mp.CreateNewElement("Element3D4N", 3, {2, 3, 4, 6});
    m.InvalidateBlocks();

    EXPECT_EQ(m.NumPoints(), 6u);
    ASSERT_EQ(m.NumCellBlocks(), 1u);  // consecutive tetras regroup into one block
    EXPECT_EQ(m.Cells(0).Type(), "tetra");
    EXPECT_EQ(m.Cells(0).NumCells(), 3u);
    // New element's connectivity is re-expressed as 0-based point indices.
    const auto conn = m.Cells(0).Conn();
    EXPECT_EQ(conn.As<std::int64_t>()[2 * 4 + 3], 5);
}

TEST(KratosBackend, GetModelPartIsReachableThroughAConstMesh) {
    // The whole point: a wrapper whose API takes `const Mesh&` must be able to
    // reach the ModelPart. Before v9.0.0 there was no const overload at all.
    const Mesh m = mt::tet_mesh();
    const ModelPart& r_mp = m.GetModelPart();
    EXPECT_EQ(r_mp.NumberOfNodes(), 5u);
    EXPECT_EQ(r_mp.NumberOfElements(), 2u);
    // Materialization is lazy, so the const call had to do the work itself.
    EXPECT_TRUE(m.IsMaterialized());
}

TEST(KratosBackend, InvalidateBlocksKeepsRaggedBlocksInPlace) {
    // A ragged block never becomes an entity, so a ModelPart edit cannot have
    // touched it -- and the rebuild must not drop it (it used to).
    Mesh m;
    m.AssignPoints(mt::points_from(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}, {0, 0, 1}}));
    m.AddPolygonBlock("polygon", {{0, 1, 2, 3}, {1, 4, 5, 2, 0}});
    m.AddCellBlock("tetra", mt::conn_from({{0, 1, 3, 6}}));
    m.AppendCellData("mark", mt::data_array({7.0, 8.0}));  // ragged block's slice
    m.AppendCellData("mark", mt::data_array({9.0}));       // tetra block's slice

    ModelPart& r_mp = m.GetModelPart();
    r_mp.CreateNewNode(8, 3.0, 0.0, 0.0);
    m.InvalidateBlocks();

    ASSERT_EQ(m.NumCellBlocks(), 2u);
    // Position preserved: the ragged block was block 0 and still is.
    EXPECT_TRUE(m.Cells(0).IsRagged());
    EXPECT_EQ(m.Cells(0).NumCells(), 2u);
    EXPECT_EQ(m.Cells(0).RowSize(0), 4u);
    EXPECT_EQ(m.Cells(0).RowSize(1), 5u);
    EXPECT_EQ(m.Cells(0).Row(1)[2], 5);
    EXPECT_EQ(m.Cells(1).Type(), "tetra");
    // ...and its cell_data slice came back with it, still block-aligned.
    ASSERT_EQ(m.CellDataNumBlocks("mark"), 2u);
    EXPECT_EQ(m.CellData("mark", 0).As<double>()[1], 8.0);
    EXPECT_EQ(m.CellData("mark", 1).As<double>()[0], 9.0);
}

TEST(KratosBackend, InvalidateBlocksReadsSubModelPartsBackAsRegions) {
    Mesh m = mt::tet_mesh();
    ModelPart& r_mp = m.GetModelPart();
    ModelPart& r_smp = r_mp.CreateSubModelPart("loaded_face");
    r_smp.AddElements({1});
    m.InvalidateBlocks();

    ASSERT_TRUE(m.HasRegion("loaded_face", meshioplusplus::RegionKind::Cell));
    const meshioplusplus::Region& r_region =
        m.Region(m.FindRegion("loaded_face", meshioplusplus::RegionKind::Cell));
    ASSERT_EQ(r_region.mEntries.Size(), 1u);
    EXPECT_EQ(r_region.mEntries.As<std::int64_t>()[0], 0);  // element 1 -> global cell 0
    // A SubModelPart carries no dim/tag, and no staged region supplied one.
    EXPECT_EQ(r_region.mDim, -1);
    EXPECT_EQ(r_region.mTag, -1);
}

TEST(KratosBackend, InvalidateBlocksKeepsRegionsWithNoSubModelPart) {
    // Side regions are deliberately never materialized as SubModelParts, so the
    // rebuild has to carry them through from staging rather than lose them.
    Mesh m = mt::tet_mesh();
    meshioplusplus::Region side;
    side.mName = "wall";
    side.mKind = meshioplusplus::RegionKind::Side;
    side.mEntries = NDArray(DType::Int64, {1, 2});
    side.mEntries.As<std::int64_t>()[0] = 0;
    side.mEntries.As<std::int64_t>()[1] = 2;
    m.AddRegion(std::move(side));

    ModelPart& r_mp = m.GetModelPart();
    r_mp.CreateNewNode(9, 5.0, 0.0, 0.0);
    m.InvalidateBlocks();

    EXPECT_TRUE(m.HasRegion("wall", meshioplusplus::RegionKind::Side));
}

TEST(KratosBackend, RoundTripThroughFormatStillWorks) {
    // A full write -> read through a real format under the KRATOS backend.
    mt::roundtrip([](const std::string& rPath,
                     const Mesh& rMesh) { meshioplusplus::write_su2(rPath, rMesh); },
                  [](const std::string& rPath) { return meshioplusplus::read_su2(rPath); },
                  mixed_dim_mesh(), ".su2");
}

// ---------------------------------------------------------------------------
// v9.1.0: entity names, properties ids, nested SubModelParts
// ---------------------------------------------------------------------------

TEST(KratosBackend, BlockEntityNamesReachTheModelPart) {
    Mesh m = mt::tet_mesh();
    EXPECT_EQ(m.BlockEntityName(0), "");  // unset means "derive from the type"
    m.SetBlockEntityName(0, "SmallDisplacementElement3D4N");
    const ModelPart& r_mp = m.GetModelPart();
    ASSERT_EQ(r_mp.NumberOfElements(), 2u);
    for (const meshioplusplus::Element& r_e : r_mp.Elements()) {
        EXPECT_TRUE(r_e.HasName());
        EXPECT_EQ(r_e.Name(), "SmallDisplacementElement3D4N");
    }
}

TEST(KratosBackend, EntityNamesSurviveTheStagingRoundTrip) {
    Mesh m = mt::tet_mesh();
    m.SetBlockEntityName(0, "TotalLagrangianElement3D4N");
    m.GetModelPart();  // materialize
    m.InvalidateBlocks();
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.BlockEntityName(0), "TotalLagrangianElement3D4N");
}

TEST(KratosBackend, MaterializeHonoursPropertyIds) {
    // Before v9.1.0 every entity landed on properties 0, silently dropping the
    // material assignment the file carried.
    Mesh m = mt::tet_mesh();
    NDArray tags(DType::Int64, {2});
    tags.As<std::int64_t>()[0] = 3;
    tags.As<std::int64_t>()[1] = 7;
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(tags));
    m.AddCellData("gmsh:physical", std::move(blocks));

    const ModelPart& r_mp = m.GetModelPart();
    ASSERT_EQ(r_mp.NumberOfElements(), 2u);
    EXPECT_EQ(r_mp.GetElement(1).PropertiesId(), 3u);
    EXPECT_EQ(r_mp.GetElement(2).PropertiesId(), 7u);
    // ...and the referenced blocks exist.
    EXPECT_TRUE(r_mp.HasProperties(3));
    EXPECT_TRUE(r_mp.HasProperties(7));
    EXPECT_FALSE(r_mp.HasProperties(5));
}

TEST(KratosBackend, NestedSubModelPartsRoundTripThroughSlashRegions) {
    Mesh m = mt::tet_mesh();
    NDArray entries(DType::Int64, {1});
    entries.As<std::int64_t>()[0] = 0;
    m.AddRegion(meshioplusplus::Region("Structure/Loads", meshioplusplus::RegionKind::Cell,
                                       std::move(entries)));

    ModelPart& r_mp = m.GetModelPart();
    // The '/' path became real nesting rather than one flatly-named part.
    ASSERT_TRUE(r_mp.HasSubModelPart("Structure"));
    ModelPart& r_parent = r_mp.GetSubModelPart("Structure");
    ASSERT_TRUE(r_parent.HasSubModelPart("Loads"));
    EXPECT_EQ(r_parent.GetSubModelPart("Loads").NumberOfElements(), 1u);

    // ...and reading it back reproduces the same single flattened name, which
    // is what makes the whole trip lossless.
    m.InvalidateBlocks();
    bool found = false;
    for (std::size_t i = 0; i < m.NumRegions(); ++i)
        if (m.Region(i).mName == "Structure/Loads" &&
            m.Region(i).mKind == meshioplusplus::RegionKind::Cell)
            found = true;
    EXPECT_TRUE(found) << "nested SubModelPart did not come back as 'Structure/Loads'";
}

TEST(KratosBackend, RegionNameThatCannotBeASubModelPartIsWarnedNotThrown) {
    // A '.' is reserved by ModelPart::FullName. Materializing must not let an
    // exception escape a lazy GetModelPart(); the region stays on the mesh.
    Mesh m = mt::tet_mesh();
    NDArray entries(DType::Int64, {1});
    entries.As<std::int64_t>()[0] = 0;
    m.AddRegion(
        meshioplusplus::Region("bad.name", meshioplusplus::RegionKind::Cell, std::move(entries)));
    ModelPart& r_mp = m.GetModelPart();
    EXPECT_FALSE(r_mp.HasSubModelPart("bad.name"));
    EXPECT_TRUE(m.HasRegion("bad.name"));
}

// --- v9.2.0: per-key control over the automatic tag pass -------------------
//
// A `.mdpa` properties id arrives as a `gmsh:physical` cell tag, so the tag
// pass synthesized a `gmsh_physical_<id>` SubModelPart beside the file's real
// ones -- material assignment surfacing as a group, which a consumer then had
// to filter by name. For a genuine gmsh file that inference is wanted, so the
// key stays in KnownTagKeys and the caller says which meaning applies.

TEST(KratosBackend, TagSubModelPartKeysDefaultsToEveryKnownKey) {
    Mesh m = mixed_dim_mesh();
    EXPECT_EQ(m.TagSubModelPartKeys(), Mesh::KnownTagKeys());
}

TEST(KratosBackend, ExcludingOneTagKeyLeavesTheOthersWorking) {
    Mesh m = mixed_dim_mesh();
    m.AddCellData("gmsh:physical", {int64_array({10, 20}), int64_array({7, 7})});
    m.AddCellData("cell_tags", {int64_array({3, 3}), int64_array({4, 4})});
    m.ExcludeTagSubModelPartKey("gmsh:physical");

    const std::vector<std::string> keys = m.TagSubModelPartKeys();
    EXPECT_EQ(std::find(keys.begin(), keys.end(), "gmsh:physical"), keys.end());
    EXPECT_NE(std::find(keys.begin(), keys.end(), "cell_tags"), keys.end());

    ModelPart& r_mp = m.GetModelPart();
    EXPECT_FALSE(r_mp.HasSubModelPart("gmsh_physical_7"));
    EXPECT_FALSE(r_mp.HasSubModelPart("gmsh_physical_10"));
    EXPECT_TRUE(r_mp.HasSubModelPart("cell_tags_3"));
    EXPECT_TRUE(r_mp.HasSubModelPart("cell_tags_4"));
    // The tag array itself is untouched, so writer round-trips are unaffected.
    EXPECT_TRUE(m.HasCellData("gmsh:physical"));
}

TEST(KratosBackend, AnInclusionListRestrictsToJustThoseKeys) {
    Mesh m = mixed_dim_mesh();
    m.AddCellData("gmsh:physical", {int64_array({10, 20}), int64_array({7, 7})});
    m.AddCellData("cell_tags", {int64_array({3, 3}), int64_array({4, 4})});
    m.SetTagSubModelPartKeys({"cell_tags"});

    ModelPart& r_mp = m.GetModelPart();
    EXPECT_TRUE(r_mp.HasSubModelPart("cell_tags_3"));
    EXPECT_FALSE(r_mp.HasSubModelPart("gmsh_physical_7"));
}

TEST(KratosBackend, AnEmptyInclusionListRestoresTheDefault) {
    Mesh m = mixed_dim_mesh();
    m.SetTagSubModelPartKeys({"cell_tags"});
    m.SetTagSubModelPartKeys({});
    EXPECT_EQ(m.TagSubModelPartKeys(), Mesh::KnownTagKeys());
}

TEST(KratosBackend, ExcludingEveryKeyMatchesDisablingTheTagPass) {
    Mesh m = mixed_dim_mesh();
    m.AddCellData("gmsh:physical", {int64_array({10, 20}), int64_array({7, 7})});
    for (const std::string& r_key : Mesh::KnownTagKeys())
        m.ExcludeTagSubModelPartKey(r_key);
    EXPECT_TRUE(m.TagSubModelPartKeys().empty());
    EXPECT_EQ(m.GetModelPart().NumberOfSubModelParts(), 0u);
}

TEST(KratosBackend, ChangingTheTagKeysAfterMaterializationReMaterializes) {
    Mesh m = mixed_dim_mesh();
    m.AddCellData("gmsh:physical", {int64_array({10, 20}), int64_array({7, 7})});
    ASSERT_TRUE(m.GetModelPart().HasSubModelPart("gmsh_physical_7"));
    m.ExcludeTagSubModelPartKey("gmsh:physical");
    EXPECT_FALSE(m.GetModelPart().HasSubModelPart("gmsh_physical_7"));
}
#endif  // MESHIOPLUSPLUS_MESH_BACKEND_KRATOS
