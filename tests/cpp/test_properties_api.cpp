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

/**
 * @file test_properties_api.cpp
 * @brief The uniform mesh API's property-set surface, which every backend must
 * implement identically -- the `test_region_api.cpp` counterpart for
 * `properties.hpp`.
 *
 * Runs under MESHIO, NATIVE and KRATOS alike (the `cpp-tests` matrix), which is
 * what makes "the three backends cannot drift on ordering or on the
 * replace-by-id rule" a check rather than a claim.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/properties.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::PropertySet;
using meshioplusplus::PropertyValue;

/** @brief A set with one numeric entry, so tests read at a glance. */
PropertySet prop_set(std::int64_t Id, const std::string& rKey, double Value) {
    PropertySet ps;
    ps.mId = Id;
    PropertyValue pv;
    pv.mKey = rKey;
    pv.mValues = meshioplusplus::NDArray(meshioplusplus::DType::Float64, {1});
    pv.mValues.As<double>()[0] = Value;
    ps.mValues.push_back(std::move(pv));
    return ps;
}

}  // namespace

TEST(PropertiesApi, AFreshMeshHasNone) {
    Mesh m = mt::tri_mesh();
    EXPECT_EQ(m.NumPropertySets(), 0u);
    EXPECT_FALSE(m.HasPropertySet(0));
    EXPECT_EQ(m.FindPropertySet(0), Mesh::npos);
}

TEST(PropertiesApi, StoredSetsComeBackAscendingByIdWhateverTheInsertionOrder) {
    Mesh m = mt::tri_mesh();
    m.AddPropertySet(prop_set(7, "DENSITY", 7850.0));
    m.AddPropertySet(prop_set(2, "DENSITY", 2700.0));
    m.AddPropertySet(prop_set(5, "DENSITY", 5000.0));

    ASSERT_EQ(m.NumPropertySets(), 3u);
    EXPECT_EQ(m.GetPropertySet(0).mId, 2);
    EXPECT_EQ(m.GetPropertySet(1).mId, 5);
    EXPECT_EQ(m.GetPropertySet(2).mId, 7);
}

TEST(PropertiesApi, AddingTheSameIdReplacesRatherThanDuplicates) {
    Mesh m = mt::tri_mesh();
    m.AddPropertySet(prop_set(3, "DENSITY", 1.0));
    m.AddPropertySet(prop_set(3, "DENSITY", 2.0));
    ASSERT_EQ(m.NumPropertySets(), 1u);
    ASSERT_EQ(m.GetPropertySet(0).mValues.size(), 1u);
    EXPECT_DOUBLE_EQ(m.GetPropertySet(0).mValues[0].mValues.As<double>()[0], 2.0);
}

TEST(PropertiesApi, FindAndHasAgreeWithTheStoredOrder) {
    Mesh m = mt::tri_mesh();
    m.AddPropertySet(prop_set(10, "A", 1.0));
    m.AddPropertySet(prop_set(4, "B", 2.0));

    EXPECT_TRUE(m.HasPropertySet(4));
    EXPECT_TRUE(m.HasPropertySet(10));
    EXPECT_FALSE(m.HasPropertySet(5));
    EXPECT_EQ(m.FindPropertySet(4), 0u);
    EXPECT_EQ(m.FindPropertySet(10), 1u);
    EXPECT_EQ(m.FindPropertySet(5), Mesh::npos);
}

TEST(PropertiesApi, TextValuesSurviveVerbatim) {
    // The whole reason properties do not go through field_data: NDArray has no
    // string dtype, so a constitutive-law name has nowhere else to live.
    Mesh m = mt::tri_mesh();
    PropertySet ps;
    ps.mId = 1;
    PropertyValue pv;
    pv.mKey = "CONSTITUTIVE_LAW";
    pv.mText = "LinearElastic3DLaw";
    ps.mValues.push_back(std::move(pv));
    m.AddPropertySet(std::move(ps));

    ASSERT_EQ(m.NumPropertySets(), 1u);
    const PropertyValue& r_v = m.GetPropertySet(0).mValues[0];
    EXPECT_TRUE(r_v.IsText());
    EXPECT_EQ(r_v.mText, "LinearElastic3DLaw");
}

TEST(PropertiesApi, SetsSurviveLaterGeometryIngestion) {
    // Property sets are keyed by id, never by entity index, so adding cells
    // after them must neither drop nor renumber anything.
    Mesh m = mt::tri_mesh();
    m.AddPropertySet(prop_set(1, "DENSITY", 7850.0));
    const Mesh other = mt::tri_mesh();
    for (const auto cb : other.CellRange()) {
        meshioplusplus::NDArray conn(cb.Conn().Dtype(), cb.Conn().Shape());
        std::memcpy(conn.Data(), cb.Conn().Data(), cb.Conn().Nbytes());
        m.AddCellBlock(cb.Type(), std::move(conn));
    }
    ASSERT_EQ(m.NumPropertySets(), 1u);
    EXPECT_DOUBLE_EQ(m.GetPropertySet(0).mValues[0].mValues.As<double>()[0], 7850.0);
}
