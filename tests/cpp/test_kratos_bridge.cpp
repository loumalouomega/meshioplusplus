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
 * @file test_kratos_bridge.cpp
 * @brief Backend-independent tests of `meshioplusplus::ModelPart` semantics
 * and the templated Kratos bridge (`kratos_bridge.hpp`), exercised against a
 * mock Kratos-like model part class.
 *
 * `model_part.hpp` / `kratos_bridge.hpp` never include `mesh.hpp`, so these
 * tests compile and run identically under every mesh backend. CI
 * additionally compile-checks the bridge against the real CoSimIO
 * `ModelPart` headers (see .github/workflows/ci.yml).
 */

// System includes
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/backends/kratos_names.hpp"
#include "meshioplusplus/backends/model_part.hpp"
#include "meshioplusplus/kratos_bridge.hpp"

using meshioplusplus::CellType;
using meshioplusplus::IndexType;
using meshioplusplus::ModelPart;

// ---------------------------------------------------------------------------
// ModelPart semantics
// ---------------------------------------------------------------------------

TEST(ModelPart, NodeCreationAndLookup) {
    ModelPart mp("Main");
    mp.CreateNewNode(1, 0.0, 0.0, 0.0);
    mp.CreateNewNode(2, 1.0, 0.5, 0.25);
    EXPECT_EQ(mp.NumberOfNodes(), 2u);
    EXPECT_TRUE(mp.HasNode(2));
    EXPECT_FALSE(mp.HasNode(3));
    EXPECT_DOUBLE_EQ(mp.GetNode(2).Y(), 0.5);
    EXPECT_THROW(mp.CreateNewNode(2, 9, 9, 9), std::invalid_argument);  // duplicate id
    EXPECT_THROW(mp.CreateNewNode(0, 0, 0, 0), std::invalid_argument);  // ids are 1-based
    EXPECT_THROW(mp.GetNode(99), std::out_of_range);
}

TEST(ModelPart, ElementCreationValidatesNodesAndNames) {
    ModelPart mp;
    for (IndexType i = 1; i <= 4; ++i)
        mp.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
    // Kratos element name, Kratos geometry name, and meshio name all resolve.
    mp.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4});
    mp.CreateNewElement("Tetrahedra3D4", 2, {1, 2, 3, 4});
    mp.CreateNewElement("tetra", 3, {1, 2, 3, 4});
    EXPECT_EQ(mp.NumberOfElements(), 3u);
    EXPECT_EQ(mp.GetElement(2).Type(), CellType::Tetra);
    // Node-count mismatch and unknown node ids throw.
    EXPECT_THROW(mp.CreateNewElement("Element3D4N", 4, {1, 2, 3}), std::invalid_argument);
    EXPECT_THROW(mp.CreateNewElement("Element3D4N", 4, {1, 2, 3, 99}), std::invalid_argument);
    EXPECT_THROW(mp.CreateNewElement("NoSuchThing42", 4, {1, 2, 3, 4}), std::invalid_argument);
}

TEST(ModelPart, SubModelPartsShareEntitiesWithRoot) {
    ModelPart root("Main");
    ModelPart& r_smp = root.CreateSubModelPart("inlet");
    ModelPart& r_nested = r_smp.CreateSubModelPart("lip");

    // Creation on a nested part inserts into the root and records
    // membership in the whole ancestor chain.
    r_nested.CreateNewNode(1, 0, 0, 0);
    EXPECT_EQ(root.NumberOfNodes(), 1u);
    EXPECT_EQ(r_smp.NumberOfNodes(), 1u);
    EXPECT_EQ(r_nested.NumberOfNodes(), 1u);
    EXPECT_TRUE(r_nested.HasNode(1));

    root.CreateNewNode(2, 1, 0, 0);
    root.CreateNewNode(3, 0, 1, 0);
    EXPECT_EQ(root.NumberOfNodes(), 3u);
    EXPECT_EQ(r_smp.NumberOfNodes(), 1u);  // membership stays local

    // AddNodes of existing root entities; unknown ids throw.
    r_smp.AddNodes({2, 3});
    EXPECT_EQ(r_smp.NumberOfNodes(), 3u);
    EXPECT_EQ(r_nested.NumberOfNodes(), 1u);
    EXPECT_THROW(r_smp.AddNodes({99}), std::invalid_argument);

    EXPECT_EQ(r_nested.FullName(), "Main.inlet.lip");
    EXPECT_EQ(&r_nested.GetRootModelPart(), &root);
    EXPECT_TRUE(r_nested.IsSubModelPart());
    EXPECT_FALSE(root.IsSubModelPart());
    EXPECT_THROW(root.CreateSubModelPart("inlet"), std::invalid_argument);  // duplicate
    EXPECT_THROW(root.CreateSubModelPart("a.b"), std::invalid_argument);    // '.' reserved
    EXPECT_THROW(root.GetSubModelPart("outlet"), std::out_of_range);
}

TEST(ModelPart, MoveKeepsParentPointersValid) {
    ModelPart root("Main");
    root.CreateSubModelPart("part");
    ModelPart moved = std::move(root);
    ModelPart& r_smp = moved.GetSubModelPart("part");
    EXPECT_EQ(&r_smp.GetRootModelPart(), &moved);
    r_smp.CreateNewNode(1, 0, 0, 0);  // must land in `moved`, not the husk
    EXPECT_EQ(moved.NumberOfNodes(), 1u);
}

// ---------------------------------------------------------------------------
// Bridge: meshioplusplus::ModelPart -> mock Kratos-like class and back
// ---------------------------------------------------------------------------

namespace {

// Minimal Kratos-shaped destination: name-string creation API, integer
// properties, nested sub model parts. Intentionally NOT meshioplusplus types.
class MockKratosModelPart {
public:
    struct Entity {
        std::string mName;
        IndexType mId;
        std::vector<IndexType> mNodes;
        IndexType mProps;
    };

    explicit MockKratosModelPart(std::string name = "Root") : mName(std::move(name)) {}

    void CreateNewNode(IndexType id, double x, double y, double z) { mNodes[id] = {x, y, z}; }
    void CreateNewElement(const std::string& rName, IndexType id, std::vector<IndexType> nodes,
                          IndexType props) {
        mElements.push_back({rName, id, std::move(nodes), props});
    }
    void CreateNewCondition(const std::string& rName, IndexType id, std::vector<IndexType> nodes,
                            IndexType props) {
        mConditions.push_back({rName, id, std::move(nodes), props});
    }
    MockKratosModelPart& CreateSubModelPart(const std::string& rName) {
        mSubs.emplace_back(new MockKratosModelPart(rName));
        return *mSubs.back();
    }
    void AddNodes(const std::vector<IndexType>& rIds) {
        mMemberNodes.insert(mMemberNodes.end(), rIds.begin(), rIds.end());
    }
    void AddElements(const std::vector<IndexType>& rIds) {
        mMemberElems.insert(mMemberElems.end(), rIds.begin(), rIds.end());
    }
    void AddConditions(const std::vector<IndexType>& rIds) {
        mMemberConds.insert(mMemberConds.end(), rIds.begin(), rIds.end());
    }

    std::string mName;
    std::map<IndexType, std::array<double, 3>> mNodes;
    std::vector<Entity> mElements, mConditions;
    std::vector<IndexType> mMemberNodes, mMemberElems, mMemberConds;
    std::vector<std::unique_ptr<MockKratosModelPart>> mSubs;
};

ModelPart sample_model_part() {
    ModelPart mp("Main");
    mp.CreateNewNode(1, 0, 0, 0);
    mp.CreateNewNode(2, 1, 0, 0);
    mp.CreateNewNode(3, 1, 1, 0);
    mp.CreateNewNode(4, 0, 1, 0);
    mp.CreateNewNode(5, 0.5, 0.5, 0.5);
    mp.CreateNewElement("Element3D4N", 1, {1, 2, 3, 5}, 7);
    mp.CreateNewElement("Element3D4N", 2, {1, 3, 4, 5}, 7);
    mp.CreateNewCondition("SurfaceCondition3D3N", 1, {1, 2, 3});
    ModelPart& r_smp = mp.CreateSubModelPart("skin");
    r_smp.AddConditions({1});
    r_smp.AddNodes({1, 2, 3});
    return mp;
}

}  // namespace

TEST(KratosBridge, ToModelPartPopulatesMock) {
    ModelPart src = sample_model_part();
    MockKratosModelPart dest("Kratos");
    meshioplusplus::to_model_part(src, dest);

    EXPECT_EQ(dest.mNodes.size(), 5u);
    EXPECT_DOUBLE_EQ(dest.mNodes.at(5)[2], 0.5);
    ASSERT_EQ(dest.mElements.size(), 2u);
    EXPECT_EQ(dest.mElements[0].mName, "Element3D4N");
    EXPECT_EQ(dest.mElements[0].mNodes, (std::vector<IndexType>{1, 2, 3, 5}));
    EXPECT_EQ(dest.mElements[0].mProps, 7u);
    ASSERT_EQ(dest.mConditions.size(), 1u);
    EXPECT_EQ(dest.mConditions[0].mName, "SurfaceCondition3D3N");

    ASSERT_EQ(dest.mSubs.size(), 1u);
    EXPECT_EQ(dest.mSubs[0]->mName, "skin");
    EXPECT_EQ(dest.mSubs[0]->mMemberConds, (std::vector<IndexType>{1}));
    EXPECT_EQ(dest.mSubs[0]->mMemberNodes, (std::vector<IndexType>{1, 2, 3}));
}

TEST(KratosBridge, PropertiesGetterIsApplied) {
    ModelPart src = sample_model_part();
    MockKratosModelPart dest;
    // Getter maps the id (e.g. to a Properties::Pointer in real Kratos);
    // here it shifts the id so we can observe it was routed through.
    meshioplusplus::to_model_part(src, dest, [](IndexType pid) { return pid + 100; });
    EXPECT_EQ(dest.mElements[0].mProps, 107u);
    EXPECT_EQ(dest.mConditions[0].mProps, 100u);
}

TEST(KratosBridge, FromModelPartRoundTrip) {
    ModelPart src = sample_model_part();
    // Our own ModelPart satisfies the default bridge_traits shape, so a
    // ModelPart -> ModelPart round-trip exercises from_model_part fully.
    ModelPart copy = meshioplusplus::from_model_part(src);
    EXPECT_EQ(copy.NumberOfNodes(), 5u);
    EXPECT_EQ(copy.NumberOfElements(), 2u);
    EXPECT_EQ(copy.NumberOfConditions(), 1u);
    EXPECT_DOUBLE_EQ(copy.GetNode(5).Z(), 0.5);
    EXPECT_EQ(copy.GetElement(2).NodeIds(), (std::vector<IndexType>{1, 3, 4, 5}));
    EXPECT_EQ(copy.GetElement(1).PropertiesId(), 7u);
    ASSERT_TRUE(copy.HasSubModelPart("skin"));
    EXPECT_EQ(copy.GetSubModelPart("skin").NumberOfConditions(), 1u);
    EXPECT_EQ(copy.GetSubModelPart("skin").NumberOfNodes(), 3u);
}

TEST(KratosBridge, NameTablesRoundTrip) {
    // Every default element/condition name must resolve back to its type.
    using meshioplusplus::cell_type_from_kratos_name;
    using meshioplusplus::kratos_condition_name;
    using meshioplusplus::kratos_element_name;
    for (CellType t : {CellType::Vertex, CellType::Line, CellType::Line3, CellType::Triangle,
                       CellType::Triangle6, CellType::Quad, CellType::Quad8, CellType::Quad9,
                       CellType::Tetra, CellType::Tetra10, CellType::Pyramid, CellType::Pyramid13,
                       CellType::Wedge, CellType::Wedge15, CellType::Hexahedron,
                       CellType::Hexahedron20, CellType::Hexahedron27}) {
        EXPECT_EQ(cell_type_from_kratos_name(kratos_element_name(t)), t) << kratos_element_name(t);
        EXPECT_EQ(cell_type_from_kratos_name(kratos_condition_name(t)), t)
            << kratos_condition_name(t);
    }
}

// v9.2.0: the applier overload iterates rSource.Properties(), which held
// id-only sets before property sets could ride on the Mesh -- so for a mesh
// read from a file it transferred nothing at all, silently.
TEST(KratosBridge, ApplierFiresForPropertiesThatCameFromTheMesh) {
    meshioplusplus::ModelPart source("Main");
    meshioplusplus::PropertySet ps;
    ps.mId = 3;
    meshioplusplus::PropertyValue density;
    density.mKey = "DENSITY";
    density.mValues = meshioplusplus::NDArray(meshioplusplus::DType::Float64, {1});
    density.mValues.As<double>()[0] = 7850.0;
    ps.mValues.push_back(std::move(density));
    source.CreateNewProperties(3).mValues = ps.mValues;

    meshioplusplus::ModelPart dest("Main");
    std::vector<std::string> applied_keys;
    meshioplusplus::to_model_part(
        source, dest, [](meshioplusplus::IndexType id) { return id; },
        [&](meshioplusplus::IndexType, const meshioplusplus::PropertyValue& rValue) {
            applied_keys.push_back(rValue.mKey);
        });
    ASSERT_EQ(applied_keys.size(), 1u);
    EXPECT_EQ(applied_keys[0], "DENSITY");
}
