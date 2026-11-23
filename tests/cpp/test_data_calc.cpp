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
