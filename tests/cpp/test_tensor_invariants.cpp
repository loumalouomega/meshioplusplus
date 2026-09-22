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
// Tests for the tensor_invariants data operation.

// System includes
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/tensor_invariants.hpp"

namespace {

using meshioplusplus::data_num_components;
using meshioplusplus::DataLocation;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::tensor_invariant_from_name;
using meshioplusplus::TensorInvariant;
using meshioplusplus::tensor_invariants;
using meshioplusplus::TensorInvariantsOptions;
using meshioplusplus::detail::read_double;

Mesh tensor_mesh() {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    // Row 0: a known tensor; row 1: hydrostatic (isotropic); row 2: has a NaN.
    m.AddPointData("s", mt::data_array({1.0, 2.0, 3.0, 0.5, 0.6, 0.7,
                                        2.0, 2.0, 2.0, 0.0, 0.0, 0.0,
                                        std::numeric_limits<double>::quiet_NaN(), 1.0, 1.0, 0.0,
                                        0.0, 0.0},
                                       6));
    m.AddPointData("g9", mt::data_array({1.0, 0.0, 0.0, 0.0, 2.0, 0.5, 0.0, 0.5, 3.0,
                                         2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0,
                                         1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                                        9));
    return m;
}

TEST(TensorInvariants, MisesMatchesClosedForm) {
    Mesh in = tensor_mesh();
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;
    o.names = {"s"};
    o.outputs = TensorInvariant::Mises;
    Mesh out = tensor_invariants(in, o);
    const double xx = 1, yy = 2, zz = 3, xy = 0.5, yz = 0.6, zx = 0.7;
    const double expect = std::sqrt(0.5 * ((xx - yy) * (xx - yy) + (yy - zz) * (yy - zz) +
                                           (zz - xx) * (zz - xx) +
                                           6.0 * (xy * xy + yz * yz + zx * zx)));
    EXPECT_NEAR(read_double(out.PointData("s_mises"), 0), expect, 1e-12);
    // Isotropic tensor: zero deviatoric stress, zero von Mises.
    EXPECT_NEAR(read_double(out.PointData("s_mises"), 1), 0.0, 1e-12);
    EXPECT_TRUE(std::isnan(read_double(out.PointData("s_mises"), 2)));
    mt::expect_same_geometry(in, out);
}

TEST(TensorInvariants, PrincipalAscendingAndSumsToTrace) {
    Mesh in = tensor_mesh();
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;
    o.names = {"s"};
    o.outputs = TensorInvariant::Principal;
    Mesh out = tensor_invariants(in, o);
    const NDArray& p = out.PointData("s_principal");
    ASSERT_EQ(p.Shape().size(), 2u);
    ASSERT_EQ(p.Shape()[1], 3u);
    const double p0 = read_double(p, 0), p1 = read_double(p, 1), p2 = read_double(p, 2);
    EXPECT_LE(p0, p1);
    EXPECT_LE(p1, p2);
    EXPECT_NEAR(p0 + p1 + p2, 1.0 + 2.0 + 3.0, 1e-9);  // trace is invariant
    // Isotropic tensor: all three eigenvalues equal 2.
    EXPECT_NEAR(read_double(p, 3), 2.0, 1e-12);
    EXPECT_NEAR(read_double(p, 4), 2.0, 1e-12);
    EXPECT_NEAR(read_double(p, 5), 2.0, 1e-12);
}

TEST(TensorInvariants, HydrostaticIsMeanDiagonal) {
    Mesh in = tensor_mesh();
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;
    o.names = {"s"};
    o.outputs = TensorInvariant::Hydrostatic;
    Mesh out = tensor_invariants(in, o);
    EXPECT_NEAR(read_double(out.PointData("s_hydrostatic"), 0), (1.0 + 2.0 + 3.0) / 3.0, 1e-12);
    EXPECT_NEAR(read_double(out.PointData("s_hydrostatic"), 1), 2.0, 1e-12);
}

TEST(TensorInvariants, DeviatoricSubtractsHydrostaticFromDiagonalOnly) {
    Mesh in = tensor_mesh();
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;
    o.names = {"s"};
    o.outputs = TensorInvariant::Deviatoric;
    Mesh out = tensor_invariants(in, o);
    const NDArray& d = out.PointData("s_deviatoric");
    const double h = (1.0 + 2.0 + 3.0) / 3.0;
    EXPECT_NEAR(read_double(d, 0), 1.0 - h, 1e-12);  // xx
    EXPECT_NEAR(read_double(d, 1), 2.0 - h, 1e-12);  // yy
    EXPECT_NEAR(read_double(d, 2), 3.0 - h, 1e-12);  // zz
    EXPECT_NEAR(read_double(d, 3), 0.5, 1e-12);       // xy unchanged
    EXPECT_NEAR(read_double(d, 4), 0.6, 1e-12);       // yz unchanged
    EXPECT_NEAR(read_double(d, 5), 0.7, 1e-12);       // zx unchanged
    // Isotropic tensor: deviatoric part is zero.
    EXPECT_NEAR(read_double(d, 6), 0.0, 1e-12);
    EXPECT_NEAR(read_double(d, 7), 0.0, 1e-12);
    EXPECT_NEAR(read_double(d, 8), 0.0, 1e-12);
}

TEST(TensorInvariants, NineComponentUsesSymmetricPartForMisesAndPrincipal) {
    Mesh in = tensor_mesh();
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;
    o.names = {"g9"};
    o.outputs = TensorInvariant::All;
    Mesh out = tensor_invariants(in, o);
    // g9 row 0 is symmetric already: xx=1 yy=2 zz=3 xy=0 yz=0.5 zx=0.
    const double expect = std::sqrt(0.5 * ((1 - 2) * (1 - 2) + (2 - 3) * (2 - 3) +
                                           (3 - 1) * (3 - 1) + 6.0 * (0.5 * 0.5)));
    EXPECT_NEAR(read_double(out.PointData("g9_mises"), 0), expect, 1e-12);
    const NDArray& dev = out.PointData("g9_deviatoric");
    ASSERT_EQ(data_num_components(dev), 9u);
    // Off-diagonal entries pass through unchanged, including the asymmetric ones.
    EXPECT_NEAR(read_double(dev, 1), 0.0, 1e-12);  // xy unchanged
    EXPECT_NEAR(read_double(dev, 5), 0.5, 1e-12);  // yz unchanged
    EXPECT_NEAR(read_double(dev, 7), 0.5, 1e-12);  // zy unchanged (not symmetrized)
}

TEST(TensorInvariants, RejectsFieldLocation) {
    Mesh in = tensor_mesh();
    TensorInvariantsOptions o;
    o.location = DataLocation::Field;
    EXPECT_THROW(tensor_invariants(in, o), std::invalid_argument);
}

TEST(TensorInvariants, RejectsUnknownName) {
    Mesh in = tensor_mesh();
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;
    o.names = {"nope"};
    EXPECT_THROW(tensor_invariants(in, o), std::invalid_argument);
}

TEST(TensorInvariants, RejectsWrongComponentCountWhenNamedExplicitly) {
    Mesh in = tensor_mesh();
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;
    o.names = {"T"};  // scalar-shaped, from a manually added array below
    in.AddPointData("T", mt::data_array({0.0, 1.0, 2.0}));
    EXPECT_THROW(tensor_invariants(in, o), std::invalid_argument);
}

TEST(TensorInvariants, AutoSelectSkipsNonTensorShapedArrays) {
    Mesh in = tensor_mesh();
    in.AddPointData("T", mt::data_array({0.0, 1.0, 2.0}));
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;  // names empty: auto-select
    o.outputs = TensorInvariant::Mises;
    Mesh out = tensor_invariants(in, o);
    EXPECT_TRUE(out.HasPointData("s_mises"));
    EXPECT_TRUE(out.HasPointData("g9_mises"));
    EXPECT_FALSE(out.HasPointData("T_mises"));
}

TEST(TensorInvariants, OverwriteFalseRefusesCollision) {
    Mesh in = tensor_mesh();
    in.AddPointData("s_mises", mt::data_array({0.0, 0.0, 0.0}));
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;
    o.names = {"s"};
    o.outputs = TensorInvariant::Mises;
    o.overwrite = false;
    EXPECT_THROW(tensor_invariants(in, o), std::invalid_argument);
}

TEST(TensorInvariants, CellDataOneOutputPerBlock) {
    Mesh in;
    in.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {2, 0, 0}}));
    in.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    in.AddCellBlock("triangle", mt::conn_from({{0, 2, 3}}));
    std::vector<NDArray> blocks;
    blocks.push_back(mt::data_array({1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0}, 9));
    blocks.push_back(mt::data_array({2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0}, 9));
    in.AddCellData("g", std::move(blocks));
    TensorInvariantsOptions o;
    o.location = DataLocation::Cell;
    o.names = {"g"};
    o.outputs = TensorInvariant::Hydrostatic;
    Mesh out = tensor_invariants(in, o);
    ASSERT_EQ(out.CellDataNumBlocks("g_hydrostatic"), 2u);
    EXPECT_NEAR(read_double(out.CellData("g_hydrostatic", 0), 0), 2.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("g_hydrostatic", 1), 0), 2.0, 1e-12);
    mt::expect_same_geometry(in, out);
}

TEST(TensorInvariants, PrefixAndSuffix) {
    Mesh in = tensor_mesh();
    TensorInvariantsOptions o;
    o.location = DataLocation::Point;
    o.names = {"s"};
    o.outputs = TensorInvariant::Mises;
    o.prefix = "p_";
    o.suffix = "_q";
    Mesh out = tensor_invariants(in, o);
    EXPECT_TRUE(out.HasPointData("p_s_mises_q"));
}

TEST(TensorInvariantsParse, NameParsing) {
    EXPECT_EQ(static_cast<unsigned>(tensor_invariant_from_name("mises")),
             static_cast<unsigned>(TensorInvariant::Mises));
    EXPECT_EQ(static_cast<unsigned>(tensor_invariant_from_name("mises,principal")),
             static_cast<unsigned>(TensorInvariant::Mises | TensorInvariant::Principal));
    EXPECT_EQ(static_cast<unsigned>(tensor_invariant_from_name("all")),
             static_cast<unsigned>(TensorInvariant::All));
    EXPECT_THROW(tensor_invariant_from_name("bogus"), std::invalid_argument);
    EXPECT_THROW(tensor_invariant_from_name(""), std::invalid_argument);
}

}  // namespace
