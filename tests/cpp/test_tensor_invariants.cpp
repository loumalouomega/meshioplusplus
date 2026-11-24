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
//                   SPDX-License-Identifier: MIT
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
#include "meshioplusplus/detail/sym3_eigen.hpp"
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

// A mesh with the tensors the algebraic flags are checked against:
//   row 0: diag(2, 3, 5)                    (9 components)
//   row 1: [[1,2,3],[0,1,4],[5,6,0]]        (det 1, an all-integer inverse)
//   row 2: a singular matrix                 (det 0)
//   row 3: one NaN component
Mesh algebra_mesh() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    m.AddPointData("a", mt::data_array({2, 0, 0, 0, 3, 0, 0, 0, 5, 1,   2, 3, 0, 1, 4, 5, 6, 0,
                                        1, 2, 3, 2, 4, 6, 0, 0, 1, nan, 0, 0, 0, 1, 0, 0, 0, 1},
                                       9));
    // The same diag(2, 3, 5) as a symmetric 6-component tensor.
    m.AddPointData(
        "s", mt::data_array(
                 {2, 3, 5, 0, 0, 0, 1, 1, 1, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0}, 6));
    return m;
}

Mesh run_algebra(const std::string& rName, TensorInvariant flag) {
    TensorInvariantsOptions o;
    o.names = {rName};
    o.outputs = flag;
    return tensor_invariants(algebra_mesh(), o);
}

TEST(TensorInvariants, TraceAndDeterminantMatchClosedForms) {
    Mesh t = run_algebra("a", TensorInvariant::Trace | TensorInvariant::Determinant);
    const NDArray& tr = t.PointData("a_trace");
    const NDArray& det = t.PointData("a_determinant");
    EXPECT_NEAR(read_double(tr, 0), 10.0, 1e-12);
    EXPECT_NEAR(read_double(det, 0), 30.0, 1e-12);
    EXPECT_NEAR(read_double(tr, 1), 2.0, 1e-12);
    EXPECT_NEAR(read_double(det, 1), 1.0, 1e-12);
    EXPECT_NEAR(read_double(det, 2), 0.0, 1e-12);
    EXPECT_TRUE(std::isnan(read_double(tr, 3)));
    EXPECT_TRUE(std::isnan(read_double(det, 3)));
}

TEST(TensorInvariants, InvariantsAreTraceSecondAndDeterminant) {
    Mesh t = run_algebra("a", TensorInvariant::Invariants);
    const NDArray& inv = t.PointData("a_invariants");
    EXPECT_EQ(data_num_components(inv), 3u);
    EXPECT_NEAR(read_double(inv, 0 * 3 + 0), 10.0, 1e-12);
    EXPECT_NEAR(read_double(inv, 0 * 3 + 1), 31.0, 1e-12);  // 6 + 15 + 10
    EXPECT_NEAR(read_double(inv, 0 * 3 + 2), 30.0, 1e-12);
    EXPECT_NEAR(read_double(inv, 1 * 3 + 0), 2.0, 1e-12);
    EXPECT_NEAR(read_double(inv, 1 * 3 + 1), -38.0, 1e-12);
    EXPECT_NEAR(read_double(inv, 1 * 3 + 2), 1.0, 1e-12);
    EXPECT_TRUE(std::isnan(read_double(inv, 3 * 3 + 1)));
}

TEST(TensorInvariants, InverseOfAGeneralMatrixIsExact) {
    Mesh t = run_algebra("a", TensorInvariant::Inverse);
    const NDArray& inv = t.PointData("a_inverse");
    ASSERT_EQ(data_num_components(inv), 9u);
    const double diag[9] = {0.5, 0, 0, 0, 1.0 / 3.0, 0, 0, 0, 0.2};
    const double gen[9] = {-24, 18, 5, 20, -15, -4, -5, 4, 1};
    for (int k = 0; k < 9; ++k) {
        EXPECT_NEAR(read_double(inv, k), diag[k], 1e-14) << k;
        EXPECT_NEAR(read_double(inv, 9 + k), gen[k], 1e-12) << k;
    }
    // A singular row and a row with a NaN are NaN throughout.
    for (int k = 0; k < 9; ++k) {
        EXPECT_TRUE(std::isnan(read_double(inv, 2 * 9 + k))) << k;
        EXPECT_TRUE(std::isnan(read_double(inv, 3 * 9 + k))) << k;
    }
}

TEST(TensorInvariants, InverseOfASymmetricInputKeepsSixComponents) {
    Mesh t = run_algebra("s", TensorInvariant::Inverse);
    const NDArray& inv = t.PointData("s_inverse");
    ASSERT_EQ(data_num_components(inv), 6u);
    const double expect[6] = {0.5, 1.0 / 3.0, 0.2, 0, 0, 0};
    for (int k = 0; k < 6; ++k)
        EXPECT_NEAR(read_double(inv, k), expect[k], 1e-14) << k;
}

TEST(TensorInvariants, TheAlgebraicFlagsAreOptInAndReadTheTensorAsGiven) {
    // The default is still the four original arrays.
    TensorInvariantsOptions o;
    o.names = {"a"};
    Mesh out = tensor_invariants(algebra_mesh(), o);
    EXPECT_TRUE(out.HasPointData("a_mises"));
    EXPECT_FALSE(out.HasPointData("a_trace"));
    EXPECT_FALSE(out.HasPointData("a_inverse"));
    // A nonsymmetric 9-component input: the determinant is that of the matrix
    // itself, not of its symmetric part (which would be 0.25 * ... != 1).
    Mesh d = run_algebra("a", TensorInvariant::Determinant);
    EXPECT_NEAR(read_double(d.PointData("a_determinant"), 1), 1.0, 1e-12);
}

TEST(TensorInvariants, InverseTimesTheMatrixIsTheIdentity) {
    Mesh in;
    in.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}));
    in.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    const double a[9] = {4, -2, 1, 3, 6, -5, 2, 7, 8};
    in.AddPointData("t", mt::data_array({4,  -2, 1, 3, 6, -5, 2, 7, 8, 4,  -2, 1, 3, 6,
                                         -5, 2,  7, 8, 4, -2, 1, 3, 6, -5, 2,  7, 8},
                                        9));
    TensorInvariantsOptions o;
    o.names = {"t"};
    o.outputs = TensorInvariant::Inverse;
    const Mesh out = tensor_invariants(in, o);
    const NDArray& inv = out.PointData("t_inverse");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double acc = 0.0;
            for (int k = 0; k < 3; ++k)
                acc += a[i * 3 + k] * read_double(inv, k * 3 + j);
            EXPECT_NEAR(acc, i == j ? 1.0 : 0.0, 1e-12);
        }
    }
}

TEST(TensorInvariants, AlgebraicFlagsOnCellData) {
    Mesh in;
    in.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {2, 0, 0}}));
    in.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    in.AddCellBlock("triangle", mt::conn_from({{0, 2, 3}}));
    std::vector<NDArray> blocks;
    blocks.push_back(mt::data_array({2, 0, 0, 0, 3, 0, 0, 0, 5}, 9));
    blocks.push_back(mt::data_array({1, 2, 3, 0, 1, 4, 5, 6, 0}, 9));
    in.AddCellData("g", std::move(blocks));
    TensorInvariantsOptions o;
    o.location = DataLocation::Cell;
    o.names = {"g"};
    o.outputs = TensorInvariant::Determinant | TensorInvariant::Inverse;
    Mesh out = tensor_invariants(in, o);
    ASSERT_EQ(out.CellDataNumBlocks("g_determinant"), 2u);
    ASSERT_EQ(out.CellDataNumBlocks("g_inverse"), 2u);
    EXPECT_NEAR(read_double(out.CellData("g_determinant", 0), 0), 30.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("g_determinant", 1), 0), 1.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("g_inverse", 1), 0), -24.0, 1e-12);
    EXPECT_EQ(data_num_components(out.CellData("g_inverse", 1)), 9u);
}

// A symmetric tensor with three distinct eigenvalues and no zero component:
// the eigenvector check holds for any basis of an eigenspace, so the rows use
// distinct values (a repeated one has no unique eigenvectors).
Mesh eigen_mesh() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    m.AddPointData("s", mt::data_array({1.0, 2.0, 3.0, 0.5, 0.6, 0.7,
                                        4.0, 1.0, 2.0, -0.3, 0.2, 0.1,
                                        2.0, 3.0, 5.0, 0.0, 0.0, 0.0,
                                        nan, 1.0, 1.0, 0.0, 0.0, 0.0},
                                       6));
    return m;
}

TEST(TensorInvariants, EigenvectorsSolveTheEigenproblemAndAreOrthonormal) {
    Mesh in = eigen_mesh();
    TensorInvariantsOptions o;
    o.names = {"s"};
    o.outputs = TensorInvariant::Principal | TensorInvariant::Eigenvectors;
    const Mesh out = tensor_invariants(in, o);
    const NDArray& vec = out.PointData("s_eigenvectors");
    const NDArray& val = out.PointData("s_principal");
    ASSERT_EQ(data_num_components(vec), 9u);
    const NDArray& src = in.PointData("s");
    for (int r = 0; r < 3; ++r) {
        const double xx = read_double(src, r * 6 + 0), yy = read_double(src, r * 6 + 1),
                     zz = read_double(src, r * 6 + 2), xy = read_double(src, r * 6 + 3),
                     yz = read_double(src, r * 6 + 4), zx = read_double(src, r * 6 + 5);
        const double A[3][3] = {{xx, xy, zx}, {xy, yy, yz}, {zx, yz, zz}};
        for (int k = 0; k < 3; ++k) {
            double v[3];
            for (int i = 0; i < 3; ++i)
                v[i] = read_double(vec, r * 9 + i * 3 + k);
            const double lambda = read_double(val, r * 3 + k);
            for (int i = 0; i < 3; ++i) {
                double av = 0.0;
                for (int j = 0; j < 3; ++j)
                    av += A[i][j] * v[j];
                EXPECT_NEAR(av, lambda * v[i], 1e-12) << r << " " << k << " " << i;
            }
            EXPECT_NEAR(v[0] * v[0] + v[1] * v[1] + v[2] * v[2], 1.0, 1e-12);
            // The sign rule: the first component above 1e-8 is positive.
            for (int i = 0; i < 3; ++i) {
                if (std::fabs(v[i]) > 1e-8) {
                    EXPECT_GT(v[i], 0.0);
                    break;
                }
            }
            for (int l = k + 1; l < 3; ++l) {
                double dot = 0.0;
                for (int i = 0; i < 3; ++i)
                    dot += v[i] * read_double(vec, r * 9 + i * 3 + l);
                EXPECT_NEAR(dot, 0.0, 1e-12);
            }
        }
    }
    for (int k = 0; k < 9; ++k)
        EXPECT_TRUE(std::isnan(read_double(vec, 3 * 9 + k))) << k;
}

TEST(TensorInvariants, EigenvectorsOfADiagonalTensorAreTheAxes) {
    // diag(2, 3, 5): the eigenvectors are the axes, in ascending value order.
    Mesh t = eigen_mesh();
    TensorInvariantsOptions o;
    o.names = {"s"};
    o.outputs = TensorInvariant::Eigenvectors;
    const Mesh out = tensor_invariants(t, o);
    const NDArray& vec = out.PointData("s_eigenvectors");
    const double identity[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for (int k = 0; k < 9; ++k)
        EXPECT_NEAR(read_double(vec, 2 * 9 + k), identity[k], 1e-14) << k;
}

TEST(TensorInvariants, EigenvaluesOfTheVectorSolverMatchPrincipalBitForBit) {
    const double t[6] = {1.0, 2.0, 3.0, 0.5, 0.6, 0.7};
    double a[3];
    double b[3];
    double v[9];
    meshioplusplus::detail::sym3_principal(t, a);
    meshioplusplus::detail::sym3_eigen(t, b, v);
    for (int k = 0; k < 3; ++k)
        EXPECT_EQ(a[k], b[k]) << k;
}

TEST(TensorInvariantsParse, NameParsing) {
    EXPECT_EQ(static_cast<unsigned>(tensor_invariant_from_name("mises")),
             static_cast<unsigned>(TensorInvariant::Mises));
    EXPECT_EQ(static_cast<unsigned>(tensor_invariant_from_name("mises,principal")),
             static_cast<unsigned>(TensorInvariant::Mises | TensorInvariant::Principal));
    EXPECT_EQ(static_cast<unsigned>(tensor_invariant_from_name("all")),
             static_cast<unsigned>(TensorInvariant::All));
    EXPECT_EQ(static_cast<unsigned>(tensor_invariant_from_name("trace,determinant,invariants,inverse")),
              static_cast<unsigned>(TensorInvariant::Trace | TensorInvariant::Determinant |
                                    TensorInvariant::Invariants | TensorInvariant::Inverse));
    EXPECT_EQ(static_cast<unsigned>(tensor_invariant_from_name("eigenvectors")),
              static_cast<unsigned>(TensorInvariant::Eigenvectors));
    EXPECT_EQ(static_cast<unsigned>(TensorInvariant::Eigenvectors), 256u);
    // `all` is still the four original flags.
    EXPECT_EQ(static_cast<unsigned>(TensorInvariant::All), 15u);
    EXPECT_THROW(tensor_invariant_from_name("bogus"), std::invalid_argument);
    EXPECT_THROW(tensor_invariant_from_name(""), std::invalid_argument);
}

}  // namespace
