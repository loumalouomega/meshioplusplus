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

TEST(KratosBackend, RoundTripThroughFormatStillWorks) {
    // A full write -> read through a real format under the KRATOS backend.
    mt::roundtrip([](const std::string& rPath,
                     const Mesh& rMesh) { meshioplusplus::write_su2(rPath, rMesh); },
                  [](const std::string& rPath) { return meshioplusplus::read_su2(rPath); },
                  mixed_dim_mesh(), ".su2");
}

#endif  // MESHIOPLUSPLUS_MESH_BACKEND_KRATOS
