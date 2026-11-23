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
// Tests for the elementwise data-expression evaluator.

// System includes
#include <cmath>
#include <stdexcept>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/data_calc.hpp"
#include "meshioplusplus/operations/gradient.hpp"

namespace {

using meshioplusplus::data_calc;
using meshioplusplus::DataCalcOptions;
using meshioplusplus::DataLocation;
using meshioplusplus::Mesh;
using meshioplusplus::detail::read_double;

DataCalcOptions point_opts(const std::string& rOutput) {
    DataCalcOptions o;
    o.location = DataLocation::Point;
    o.output = rOutput;
    return o;
}

TEST(DataCalc, NormEqualsEuclideanMagnitude) {
    Mesh in = mt::data_mesh();
    Mesh out = data_calc(in, "norm(v)", point_opts("speed"));
    ASSERT_TRUE(out.HasPointData("speed"));
    // v rows: (1,0,0) (0,1,0) (0,0,1) (1,1,0) (2,0,0) (0,2,0)
    EXPECT_NEAR(read_double(out.PointData("speed"), 0), 1.0, 1e-12);
    EXPECT_NEAR(read_double(out.PointData("speed"), 3), std::sqrt(2.0), 1e-12);
    EXPECT_NEAR(read_double(out.PointData("speed"), 4), 2.0, 1e-12);
    // norm collapses to a scalar, so the result is 1-D.
    EXPECT_EQ(out.PointData("speed").Shape().size(), 1u);
    mt::expect_same_geometry(in, out);
}

TEST(DataCalc, LinearCombination) {
    Mesh in = mt::data_mesh();
    Mesh out = data_calc(in, "2 * T + 1", point_opts("y"));
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_NEAR(read_double(out.PointData("y"), i),
                    2.0 * read_double(in.PointData("T"), i) + 1.0, 1e-12);
}

TEST(DataCalc, OperatorPrecedenceAndParens) {
    Mesh in = mt::data_mesh();
    Mesh a = data_calc(in, "1 + 2 * 3", point_opts("a"));
    EXPECT_NEAR(read_double(a.PointData("a"), 0), 7.0, 1e-12);
    Mesh b = data_calc(in, "(1 + 2) * 3", point_opts("b"));
    EXPECT_NEAR(read_double(b.PointData("b"), 0), 9.0, 1e-12);
}

TEST(DataCalc, UnaryMinus) {
    Mesh in = mt::data_mesh();
    Mesh a = data_calc(in, "-T", point_opts("a"));
    EXPECT_NEAR(read_double(a.PointData("a"), 1), -1.0, 1e-12);
    Mesh b = data_calc(in, "2 * -3", point_opts("b"));
    EXPECT_NEAR(read_double(b.PointData("b"), 0), -6.0, 1e-12);
    Mesh c = data_calc(in, "--T", point_opts("c"));
    EXPECT_NEAR(read_double(c.PointData("c"), 1), 1.0, 1e-12);
}

TEST(DataCalc, Functions) {
    Mesh in = mt::data_mesh();
    Mesh a = data_calc(in, "abs(0 - T)", point_opts("a"));
    EXPECT_NEAR(read_double(a.PointData("a"), 2), 11.0, 1e-12);
    Mesh b = data_calc(in, "sqrt(T)", point_opts("b"));
    EXPECT_NEAR(read_double(b.PointData("b"), 3), std::sqrt(10.0), 1e-12);
    Mesh c = data_calc(in, "min(T, 5)", point_opts("c"));
    EXPECT_NEAR(read_double(c.PointData("c"), 2), 5.0, 1e-12);
    Mesh d = data_calc(in, "max(T, 5)", point_opts("d"));
    EXPECT_NEAR(read_double(d.PointData("d"), 0), 5.0, 1e-12);
}

TEST(DataCalc, MathFunctionsMatchTheLibm) {
    Mesh in = mt::data_mesh();  // T = {0, 1, 11, 10, 2, 12}
    const double t[6] = {0.0, 1.0, 11.0, 10.0, 2.0, 12.0};
    struct Case {
        const char* mExpr;
        double (*mFn)(double);
    };
    const Case cases[] = {
        {"exp(T / 8)", [](double x) { return std::exp(x / 8); }},
        {"log(T + 1)", [](double x) { return std::log(x + 1); }},
        {"log10(T + 1)", [](double x) { return std::log10(x + 1); }},
        {"sin(T)", [](double x) { return std::sin(x); }},
        {"cos(T)", [](double x) { return std::cos(x); }},
        {"tan(T / 8)", [](double x) { return std::tan(x / 8); }},
        {"atan(T)", [](double x) { return std::atan(x); }},
        {"sinh(T / 8)", [](double x) { return std::sinh(x / 8); }},
        {"cosh(T / 8)", [](double x) { return std::cosh(x / 8); }},
        {"tanh(T / 8)", [](double x) { return std::tanh(x / 8); }},
        {"floor(T / 3)", [](double x) { return std::floor(x / 3); }},
        {"ceil(T / 3)", [](double x) { return std::ceil(x / 3); }},
        {"sign(T - 5)", [](double x) { return x - 5 > 0 ? 1.0 : (x - 5 < 0 ? -1.0 : 0.0); }},
    };
    for (const Case& c : cases) {
        const Mesh out = data_calc(in, c.mExpr, point_opts("o"));
        for (int i = 0; i < 6; ++i)
            EXPECT_NEAR(read_double(out.PointData("o"), i), c.mFn(t[i]), 1e-12) << c.mExpr << i;
    }
    // asin and acos on [-1, 1]; outside it the answer is NaN, not an error.
    const Mesh as = data_calc(in, "asin(T / 12)", point_opts("o"));
    EXPECT_NEAR(read_double(as.PointData("o"), 3), std::asin(10.0 / 12), 1e-12);
    const Mesh ac = data_calc(in, "acos(T / 6)", point_opts("o"));
    EXPECT_TRUE(std::isnan(read_double(ac.PointData("o"), 2)));
    EXPECT_NEAR(read_double(ac.PointData("o"), 4), std::acos(2.0 / 6), 1e-12);
}

TEST(DataCalc, LogOfZeroIsAnInfinityAndOfANegativeANaN) {
    Mesh in = mt::data_mesh();
    EXPECT_TRUE(
        std::isinf(read_double(data_calc(in, "log(T)", point_opts("o")).PointData("o"), 0)));
    EXPECT_TRUE(
        std::isnan(read_double(data_calc(in, "log(-T)", point_opts("o")).PointData("o"), 1)));
}

TEST(DataCalc, BinaryMathFunctionsBroadcastAScalar) {
    Mesh in = mt::data_mesh();
    const Mesh p = data_calc(in, "pow(T, 2)", point_opts("o"));
    EXPECT_NEAR(read_double(p.PointData("o"), 2), 121.0, 1e-12);
    const Mesh q = data_calc(in, "pow(2, T / 4)", point_opts("o"));
    EXPECT_NEAR(read_double(q.PointData("o"), 4), std::pow(2.0, 0.5), 1e-12);
    const Mesh a = data_calc(in, "atan2(v, 2)", point_opts("o"));
    EXPECT_EQ(meshioplusplus::data_num_components(a.PointData("o")), 3u);
    EXPECT_NEAR(read_double(a.PointData("o"), 1 * 3 + 1), std::atan2(1.0, 2.0), 1e-12);
    // A fractional power of a negative base is a domain error: NaN.
    EXPECT_TRUE(std::isnan(
        read_double(data_calc(in, "pow(0 - T, 0.5)", point_opts("o")).PointData("o"), 1)));
}

Mesh tuple_mesh() {
    Mesh in = mt::data_mesh();  // v: (1,0,0) (0,1,0) (0,0,1) (1,1,0) (2,0,0) (0,2,0)
    in.AddPointData("t9", mt::data_array({1, 2, 3, 0,  1, 4, 5, 6, 0,  //
                                          2, 0, 0, 0,  3, 0, 0, 0, 5,  //
                                          1, 0, 0, 0,  1, 0, 0, 0, 1,  //
                                          0, 1, 0, -1, 0, 0, 0, 0, 2,  //
                                          4, 1, 2, 1,  5, 3, 2, 3, 6,  //
                                          1, 1, 1, 1,  1, 1, 1, 1, 1},
                                         9));
    in.AddPointData("t6", mt::data_array({2, 3, 5, 0.5, 0.25, 0.75,  //
                                          2, 3, 5, 0.0, 0.0,  0.0,   //
                                          1, 1, 1, 0.0, 0.0,  0.0,   //
                                          1, 2, 3, 0.1, 0.2,  0.3,   //
                                          4, 5, 6, 1.0, 1.0,  1.0,   //
                                          1, 1, 1, 1.0, 1.0,  1.0},
                                         6));
    return in;
}

TEST(DataCalc, DotAndCross) {
    Mesh in = tuple_mesh();
    const Mesh d = data_calc(in, "dot(v, v * 2)", point_opts("o"));
    EXPECT_EQ(meshioplusplus::data_num_components(d.PointData("o")), 1u);
    EXPECT_NEAR(read_double(d.PointData("o"), 4), 8.0, 1e-12);  // v = (2,0,0)
    const Mesh c = data_calc(in, "cross(v, 1 + v)", point_opts("o"));
    EXPECT_EQ(meshioplusplus::data_num_components(c.PointData("o")), 3u);
    // v = (1,0,0): (1,0,0) x (2,1,1) = (0*1-0*1, 0*2-1*1, 1*1-0*2) = (0,-1,1)
    EXPECT_NEAR(read_double(c.PointData("o"), 0), 0.0, 1e-12);
    EXPECT_NEAR(read_double(c.PointData("o"), 1), -1.0, 1e-12);
    EXPECT_NEAR(read_double(c.PointData("o"), 2), 1.0, 1e-12);
}

TEST(DataCalc, DetTraceAndTranspose) {
    Mesh in = tuple_mesh();
    const Mesh det = data_calc(in, "det(t9)", point_opts("o"));
    EXPECT_NEAR(read_double(det.PointData("o"), 0), 1.0, 1e-12);   // [[1,2,3],[0,1,4],[5,6,0]]
    EXPECT_NEAR(read_double(det.PointData("o"), 1), 30.0, 1e-12);  // diag(2,3,5)
    EXPECT_NEAR(read_double(det.PointData("o"), 5), 0.0, 1e-12);   // all ones
    const Mesh tr = data_calc(in, "trace(t9)", point_opts("o"));
    EXPECT_NEAR(read_double(tr.PointData("o"), 0), 2.0, 1e-12);
    EXPECT_NEAR(read_double(tr.PointData("o"), 1), 10.0, 1e-12);
    const Mesh tp = data_calc(in, "transpose(t9)", point_opts("o"));
    EXPECT_EQ(meshioplusplus::data_num_components(tp.PointData("o")), 9u);
    // Row 0 is [[1,2,3],[0,1,4],[5,6,0]]; its transpose is [[1,0,5],[2,1,6],[3,4,0]].
    const double expect[9] = {1, 0, 5, 2, 1, 6, 3, 4, 0};
    for (int k = 0; k < 9; ++k)
        EXPECT_NEAR(read_double(tp.PointData("o"), k), expect[k], 1e-12) << k;
    // A symmetric 6-component tensor: diag(2,3,5) has det 30 and trace 10.
    const Mesh d6 = data_calc(in, "det(t6)", point_opts("o"));
    EXPECT_NEAR(read_double(d6.PointData("o"), 1), 30.0, 1e-12);
    const Mesh t6 = data_calc(in, "trace(t6)", point_opts("o"));
    EXPECT_NEAR(read_double(t6.PointData("o"), 1), 10.0, 1e-12);
    // mag is norm.
    const Mesh mg = data_calc(in, "mag(v)", point_opts("o"));
    EXPECT_NEAR(read_double(mg.PointData("o"), 4), 2.0, 1e-12);
}

TEST(DataCalc, MatmulOuterAndInv) {
    Mesh in = tuple_mesh();  // t9 row 0 = [[1,2,3],[0,1,4],[5,6,0]], det 1
    const Mesh inv = data_calc(in, "inv(t9)", point_opts("o"));
    ASSERT_EQ(meshioplusplus::data_num_components(inv.PointData("o")), 9u);
    const double expect[9] = {-24, 18, 5, 20, -15, -4, -5, 4, 1};
    for (int k = 0; k < 9; ++k)
        EXPECT_NEAR(read_double(inv.PointData("o"), k), expect[k], 1e-12) << k;
    // The singular row (all ones) is NaN throughout.
    for (int k = 0; k < 9; ++k)
        EXPECT_TRUE(std::isnan(read_double(inv.PointData("o"), 5 * 9 + k))) << k;

    // matmul of a matrix and its inverse is the identity on a regular row.
    const Mesh id = data_calc(in, "matmul(t9, inv(t9))", point_opts("o"));
    const double eye[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for (int k = 0; k < 9; ++k)
        EXPECT_NEAR(read_double(id.PointData("o"), k), eye[k], 1e-12) << k;

    // matmul by a vector: row 4 of v is (2,0,0), so t9 row 4 times it is its first column * 2.
    const Mesh mv = data_calc(in, "matmul(t9, v)", point_opts("o"));
    ASSERT_EQ(meshioplusplus::data_num_components(mv.PointData("o")), 3u);
    EXPECT_NEAR(read_double(mv.PointData("o"), 4 * 3 + 0), 8.0, 1e-12);  // 4 * 2
    EXPECT_NEAR(read_double(mv.PointData("o"), 4 * 3 + 1), 2.0, 1e-12);  // 1 * 2
    EXPECT_NEAR(read_double(mv.PointData("o"), 4 * 3 + 2), 4.0, 1e-12);  // 2 * 2

    // outer: v = (1,1,0) at row 3; outer(v, v) has ones in the top-left 2x2.
    const Mesh ou = data_calc(in, "outer(v, v)", point_opts("o"));
    ASSERT_EQ(meshioplusplus::data_num_components(ou.PointData("o")), 9u);
    const double ones[9] = {1, 1, 0, 1, 1, 0, 0, 0, 0};
    for (int k = 0; k < 9; ++k)
        EXPECT_NEAR(read_double(ou.PointData("o"), 3 * 9 + k), ones[k], 1e-12) << k;

    // A symmetric 6-component tensor keeps six components: diag(2,3,5).
    const Mesh i6 = data_calc(in, "inv(t6)", point_opts("o"));
    ASSERT_EQ(meshioplusplus::data_num_components(i6.PointData("o")), 6u);
    EXPECT_NEAR(read_double(i6.PointData("o"), 1 * 6 + 0), 0.5, 1e-14);
    EXPECT_NEAR(read_double(i6.PointData("o"), 1 * 6 + 1), 1.0 / 3.0, 1e-14);
    EXPECT_NEAR(read_double(i6.PointData("o"), 1 * 6 + 2), 0.2, 1e-14);
}

TEST(DataCalc, CoordinatesAtThePoints) {
    Mesh in = mt::data_mesh();  // T = x + 10 y
    const Mesh d = data_calc(in, "x + 10 * y - T", point_opts("o"));
    for (int i = 0; i < 6; ++i)
        EXPECT_NEAR(read_double(d.PointData("o"), i), 0.0, 1e-12) << i;
    const Mesh c = data_calc(in, "coords", point_opts("o"));
    ASSERT_EQ(meshioplusplus::data_num_components(c.PointData("o")), 3u);
    EXPECT_NEAR(read_double(c.PointData("o"), 5 * 3 + 0), 2.0, 1e-12);
    EXPECT_NEAR(read_double(c.PointData("o"), 5 * 3 + 1), 1.0, 1e-12);
    const Mesh n = data_calc(in, "norm(coords)", point_opts("o"));
    EXPECT_NEAR(read_double(n.PointData("o"), 5), std::sqrt(5.0), 1e-12);
}

TEST(DataCalc, CoordinatesAtTheCellsAreTheCentroids) {
    Mesh in = mt::data_mesh();
    DataCalcOptions o;
    o.location = DataLocation::Cell;
    o.output = "c";
    const Mesh out = data_calc(in, "coords", o);
    ASSERT_EQ(out.CellDataNumBlocks("c"), 2u);
    // The first triangle (0,1,2): centroid (2/3, 1/3, 0); the quad (1,4,5,2): (1.5, 0.5, 0).
    EXPECT_NEAR(read_double(out.CellData("c", 0), 0), 2.0 / 3.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("c", 0), 1), 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("c", 1), 0), 1.5, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("c", 1), 1), 0.5, 1e-12);
    const Mesh s = data_calc(in, "x + y", o);
    EXPECT_NEAR(read_double(s.CellData("c", 0), 1), 1.0, 1e-12);
    EXPECT_NEAR(read_double(s.CellData("c", 1), 0), 2.0, 1e-12);
}

TEST(DataCalc, ADataArrayNamedLikeACoordinateShadowsIt) {
    Mesh in = mt::data_mesh();
    in.AddPointData("x", mt::data_array({100, 200, 300, 400, 500, 600}));
    const Mesh out = data_calc(in, "x + y", point_opts("o"));
    EXPECT_NEAR(read_double(out.PointData("o"), 1), 200.0, 1e-12);  // 200 + y(0)
    EXPECT_NEAR(read_double(out.PointData("o"), 2), 301.0, 1e-12);  // 300 + y(1)
}

TEST(DataCalc, GradDivCurlMatchTheGradientOperation) {
    Mesh in = mt::data_mesh();  // T = x + 10 y: grad T = (1, 10, 0) in every cell
    DataCalcOptions cell;
    cell.location = DataLocation::Cell;
    cell.output = "o";
    const Mesh g = data_calc(in, "grad(T)", cell);
    ASSERT_EQ(g.CellDataNumBlocks("o"), 2u);
    ASSERT_EQ(meshioplusplus::data_num_components(g.CellData("o", 0)), 3u);
    for (std::size_t b = 0; b < 2; ++b) {
        EXPECT_NEAR(read_double(g.CellData("o", b), 0), 1.0, 1e-10);
        EXPECT_NEAR(read_double(g.CellData("o", b), 1), 10.0, 1e-10);
        EXPECT_NEAR(read_double(g.CellData("o", b), 2), 0.0, 1e-10);
    }
    // div and curl of the vector field v, against the operation itself.
    meshioplusplus::GradientOptions go;
    go.mArrayName = "v";
    go.mLocation = DataLocation::Cell;
    go.mOperator = meshioplusplus::GradientOperator::Divergence;
    const Mesh ref = meshioplusplus::gradient(in, go).mMesh;
    const Mesh d = data_calc(in, "div(v)", cell);
    for (std::size_t b = 0; b < 2; ++b)
        for (std::size_t i = 0; i < meshioplusplus::detail::rows(ref.CellData("v:divergence", b));
             ++i)
            EXPECT_NEAR(read_double(d.CellData("o", b), i),
                        read_double(ref.CellData("v:divergence", b), i), 1e-12);
    go.mOperator = meshioplusplus::GradientOperator::Curl;
    const Mesh cref = meshioplusplus::gradient(in, go).mMesh;
    const Mesh c = data_calc(in, "curl(v)", cell);
    ASSERT_EQ(meshioplusplus::data_num_components(c.CellData("o", 0)), 3u);
    EXPECT_NEAR(read_double(c.CellData("o", 0), 2), read_double(cref.CellData("v:curl", 0), 2),
                1e-12);
    // At the point location the call is the gradient averaged to the points.
    const Mesh p = data_calc(in, "dot(grad(T), grad(T))", point_opts("o"));
    for (int i = 0; i < 6; ++i)
        EXPECT_NEAR(read_double(p.PointData("o"), i), 101.0, 1e-9) << i;
    // A call combines with an array of the expression's own location.
    const Mesh m = data_calc(in, "mat * norm(grad(T))", cell);
    EXPECT_NEAR(read_double(m.CellData("o", 0), 1), 2.0 * std::sqrt(101.0), 1e-9);
    EXPECT_NEAR(read_double(m.CellData("o", 1), 0), 3.0 * std::sqrt(101.0), 1e-9);
}

TEST(DataCalc, AreaVolumeAndSizeOfTheCells) {
    Mesh in = mt::data_mesh();  // two triangles of area 1/2 and a unit quad
    DataCalcOptions cell;
    cell.location = DataLocation::Cell;
    cell.output = "o";
    const Mesh a = data_calc(in, "area()", cell);
    EXPECT_NEAR(read_double(a.CellData("o", 0), 1), 0.5, 1e-12);
    EXPECT_NEAR(read_double(a.CellData("o", 1), 0), 1.0, 1e-12);
    const Mesh s = data_calc(in, "size() * 2", cell);
    EXPECT_NEAR(read_double(s.CellData("o", 1), 0), 2.0, 1e-12);
    const Mesh v = data_calc(in, "volume()", cell);
    EXPECT_TRUE(std::isnan(read_double(v.CellData("o", 0), 0)));  // no 3-D cell here
    Mesh hex = mt::make_mesh(
        {{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 1}, {2, 1, 1}, {0, 1, 1}},
        "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    EXPECT_NEAR(read_double(data_calc(hex, "volume()", cell).CellData("o", 0), 0), 2.0, 1e-12);
    EXPECT_TRUE(std::isnan(read_double(data_calc(hex, "area()", cell).CellData("o", 0), 0)));
}

TEST(DataCalc, ScalarBroadcastAgainstAVector) {
    Mesh in = mt::data_mesh();
    Mesh out = data_calc(in, "v * 2", point_opts("v2"));
    ASSERT_EQ(out.PointData("v2").Shape().size(), 2u);
    EXPECT_EQ(out.PointData("v2").Shape()[1], 3u);
    EXPECT_NEAR(read_double(out.PointData("v2"), 0), 2.0, 1e-12);
    // The other operand order broadcasts too.
    Mesh out2 = data_calc(in, "2 * v", point_opts("v3"));
    EXPECT_NEAR(read_double(out2.PointData("v3"), 0), 2.0, 1e-12);
}

TEST(DataCalc, VectorMinusVectorStaysAVector) {
    Mesh in = mt::data_mesh();
    Mesh out = data_calc(in, "v - v", point_opts("z"));
    ASSERT_EQ(out.PointData("z").Shape().size(), 2u);
    EXPECT_EQ(out.PointData("z").Shape()[1], 3u);
    for (std::size_t i = 0; i < out.PointData("z").Size(); ++i)
        EXPECT_NEAR(read_double(out.PointData("z"), i), 0.0, 1e-12);
}

TEST(DataCalc, ColonAndBacktickIdentifiers) {
    Mesh in = mt::data_mesh();
    in.AddPointData("gmsh:physical", mt::data_array({1, 2, 3, 4, 5, 6}));
    in.AddPointData("with space", mt::data_array({1, 1, 1, 1, 1, 1}));

    Mesh a = data_calc(in, "gmsh:physical + 1", point_opts("a"));
    EXPECT_NEAR(read_double(a.PointData("a"), 0), 2.0, 1e-12);

    Mesh b = data_calc(in, "`with space` * 3", point_opts("b"));
    EXPECT_NEAR(read_double(b.PointData("b"), 0), 3.0, 1e-12);
}

TEST(DataCalc, CellLocationEvaluatesPerBlock) {
    Mesh in = mt::data_mesh();
    DataCalcOptions o;
    o.location = DataLocation::Cell;
    o.output = "mat2";
    Mesh out = data_calc(in, "mat * 10", o);
    ASSERT_TRUE(out.HasCellData("mat2"));
    // Exactly one array per cell block.
    ASSERT_EQ(out.CellDataNumBlocks("mat2"), out.NumCellBlocks());
    EXPECT_NEAR(read_double(out.CellData("mat2", 0), 0), 10.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("mat2", 0), 1), 20.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("mat2", 1), 0), 30.0, 1e-12);
}

TEST(DataCalc, FieldLocation) {
    Mesh in = mt::data_mesh();
    DataCalcOptions o;
    o.location = DataLocation::Field;
    o.output = "meta2";
    Mesh out = data_calc(in, "meta + 1", o);
    ASSERT_TRUE(out.HasFieldData("meta2"));
    EXPECT_NEAR(read_double(out.FieldData("meta2"), 0), 2.0, 1e-12);
    EXPECT_NEAR(read_double(out.FieldData("meta2"), 2), 4.0, 1e-12);
}

TEST(DataCalc, DivisionByZeroIsNotAnError) {
    Mesh in = mt::data_mesh();
    Mesh out = data_calc(in, "1 / 0", point_opts("inf"));
    EXPECT_TRUE(std::isinf(read_double(out.PointData("inf"), 0)));
}

// --- error paths ----------------------------------------------------------

void expect_message(const std::string& rExpr, const std::string& rNeedle) {
    Mesh in = mt::data_mesh();
    try {
        data_calc(in, rExpr, point_opts("out"));
        FAIL() << "expected an exception for '" << rExpr << "'";
    } catch (const std::invalid_argument& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find(rNeedle), std::string::npos)
            << "expression '" << rExpr << "' produced: " << msg;
    }
}

TEST(DataCalc, UnknownArrayNameListsAvailableKeys) {
    expect_message("temp + 1", "unknown point_data array 'temp'");
    expect_message("temp + 1", "available:");
    // field_data has no geometry: a coordinate name is just an unknown array.
    DataCalcOptions fo;
    fo.location = DataLocation::Field;
    fo.output = "o";
    EXPECT_THROW(data_calc(mt::data_mesh(), "x + 1", fo), std::invalid_argument);
}

TEST(DataCalc, UnknownFunctionListsKnownOnes) {
    expect_message("foo(T)", "unknown function 'foo'");
    expect_message("foo(T)", "abs, sqrt, min, max, norm, exp, log");
}

TEST(DataCalc, ArityErrors) {
    expect_message("norm(v, T)", "takes exactly 1 argument");
    expect_message("min(T)", "takes exactly 2 arguments");
    expect_message("pow(T)", "takes exactly 2 arguments");
    expect_message("sin(T, T)", "takes exactly 1 argument");
    expect_message("dot(v, T)", "'dot' needs two operands of the same width");
    expect_message("cross(v, T)", "'cross' needs two 3-component operands");
    expect_message("det(v)", "'det' needs a 6-component");
    expect_message("transpose(v)", "'transpose' needs a 9-component");
    expect_message("matmul(v, v)", "'matmul' needs a 9-component");
    expect_message("outer(v, T)", "'outer' needs two 3-component operands");
    expect_message("inv(v)", "'inv' needs a 6-component");
    expect_message("grad(T + 1)", "needs a point_data array name");
    expect_message("grad(grad(T))", "needs a point_data array name");
    expect_message("grad(mat)", "is not one");
    expect_message("area()", "per-cell quantity");
    expect_message("area(T)", "takes exactly 0 arguments");
    DataCalcOptions mo;
    mo.location = DataLocation::Field;
    mo.output = "o";
    EXPECT_THROW(data_calc(mt::data_mesh(), "grad(T)", mo), std::invalid_argument);
}

TEST(DataCalc, SyntaxErrors) {
    expect_message("T +", "unexpected end of expression");
    expect_message("(T", "expected ')'");
    expect_message("T T", "trailing input");
    expect_message("T # 1", "unexpected character '#'");
    expect_message("", "the expression is empty");
}

TEST(DataCalc, WidthMismatchIsDiagnosed) {
    Mesh in = mt::data_mesh();
    in.AddPointData("t9", mt::data_array(std::vector<double>(6 * 9, 1.0), 9));
    try {
        data_calc(in, "v * t9", point_opts("bad"));
        FAIL() << "expected an exception";
    } catch (const std::invalid_argument& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("cannot combine"), std::string::npos) << msg;
    }
}

TEST(DataCalc, DepthGuard) {
    Mesh in = mt::data_mesh();
    std::string expr(200, '(');
    expr += "T";
    expr.append(200, ')');
    try {
        data_calc(in, expr, point_opts("deep"));
        FAIL() << "expected an exception";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("nests deeper"), std::string::npos);
    }
}

TEST(DataCalc, ExistingOutputNameNeedsOverwrite) {
    Mesh in = mt::data_mesh();
    EXPECT_THROW(data_calc(in, "T * 2", point_opts("T")), std::invalid_argument);
    DataCalcOptions o = point_opts("T");
    o.overwrite = true;
    Mesh out = data_calc(in, "T * 2", o);
    EXPECT_NEAR(read_double(out.PointData("T"), 1), 2.0, 1e-12);
}

TEST(DataCalc, ErrorMessagesCarryAPosition) {
    expect_message("T + nope", "at position 4");
}

}  // namespace
