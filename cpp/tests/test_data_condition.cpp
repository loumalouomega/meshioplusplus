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
// Tests for data value conditioning (clamp / normalize / standardize).

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
#include "meshioplusplus/operations/data_condition.hpp"

namespace {

using meshioplusplus::ConditionMode;
using meshioplusplus::ConditionScope;
using meshioplusplus::data_condition;
using meshioplusplus::DataConditionOptions;
using meshioplusplus::DataLocation;
using meshioplusplus::Mesh;
using meshioplusplus::NanPolicy;
using meshioplusplus::detail::read_double;

DataConditionOptions point_opts(ConditionMode mode) {
    DataConditionOptions o;
    o.location = DataLocation::Point;
    o.mode = mode;
    return o;
}

TEST(DataCondition, ClampBoundsValues) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o = point_opts(ConditionMode::Clamp);
    o.names = {"T"};
    o.lo = 1.0;
    o.hi = 10.0;
    Mesh out = data_condition(in, o);
    // T = {0, 1, 11, 10, 2, 12}
    EXPECT_DOUBLE_EQ(read_double(out.PointData("T"), 0), 1.0);   // clamped up
    EXPECT_DOUBLE_EQ(read_double(out.PointData("T"), 1), 1.0);   // unchanged
    EXPECT_DOUBLE_EQ(read_double(out.PointData("T"), 2), 10.0);  // clamped down
    EXPECT_DOUBLE_EQ(read_double(out.PointData("T"), 3), 10.0);  // unchanged
    EXPECT_DOUBLE_EQ(read_double(out.PointData("T"), 5), 10.0);  // clamped down
    mt::expect_same_geometry(in, out);
}

TEST(DataCondition, NormalizeMapsMinToZeroAndMaxToOne) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o = point_opts(ConditionMode::Normalize);
    o.names = {"T"};
    Mesh out = data_condition(in, o);
    // T min = 0 (index 0), max = 12 (index 5)
    EXPECT_NEAR(read_double(out.PointData("T"), 0), 0.0, 1e-12);
    EXPECT_NEAR(read_double(out.PointData("T"), 5), 1.0, 1e-12);
    EXPECT_NEAR(read_double(out.PointData("T"), 1), 1.0 / 12.0, 1e-12);
}

TEST(DataCondition, NormalizeToACustomRange) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o = point_opts(ConditionMode::Normalize);
    o.names = {"T"};
    o.lo = -1.0;
    o.hi = 1.0;
    Mesh out = data_condition(in, o);
    EXPECT_NEAR(read_double(out.PointData("T"), 0), -1.0, 1e-12);
    EXPECT_NEAR(read_double(out.PointData("T"), 5), 1.0, 1e-12);
}

TEST(DataCondition, StandardizeGivesZeroMeanUnitStd) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o = point_opts(ConditionMode::Standardize);
    o.names = {"T"};
    Mesh out = data_condition(in, o);

    double sum = 0.0;
    double sumsq = 0.0;
    const std::size_t n = out.PointData("T").Size();
    for (std::size_t i = 0; i < n; ++i) {
        const double v = read_double(out.PointData("T"), i);
        sum += v;
        sumsq += v * v;
    }
    const double mean = sum / static_cast<double>(n);
    const double var = sumsq / static_cast<double>(n) - mean * mean;
    EXPECT_NEAR(mean, 0.0, 1e-12);
    EXPECT_NEAR(std::sqrt(var), 1.0, 1e-12);
}

TEST(DataCondition, ConstantArrayNormalizesToTheTargetLowerBound) {
    Mesh in = mt::data_mesh();
    in.AddPointData("k", mt::data_array({5, 5, 5, 5, 5, 5}));
    DataConditionOptions o = point_opts(ConditionMode::Normalize);
    o.names = {"k"};
    Mesh out = data_condition(in, o);
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_DOUBLE_EQ(read_double(out.PointData("k"), i), 0.0);
}

TEST(DataCondition, ComponentScopeConditionsEachComponentIndependently) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o = point_opts(ConditionMode::Normalize);
    o.names = {"v"};
    Mesh out = data_condition(in, o);
    // Component 0 of v is {1,0,0,1,2,0}: min 0, max 2 -> normalized {0.5,0,0,0.5,1,0}
    EXPECT_NEAR(read_double(out.PointData("v"), 0 * 3 + 0), 0.5, 1e-12);
    EXPECT_NEAR(read_double(out.PointData("v"), 4 * 3 + 0), 1.0, 1e-12);
    EXPECT_NEAR(read_double(out.PointData("v"), 1 * 3 + 0), 0.0, 1e-12);
}

TEST(DataCondition, MagnitudeScopePreservesDirection) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o = point_opts(ConditionMode::Normalize);
    o.names = {"v"};
    o.scope = ConditionScope::Magnitude;
    Mesh out = data_condition(in, o);
    // Row 3 of v is (1,1,0): after rescaling it must still point along (1,1,0).
    const double x = read_double(out.PointData("v"), 3 * 3 + 0);
    const double y = read_double(out.PointData("v"), 3 * 3 + 1);
    const double z = read_double(out.PointData("v"), 3 * 3 + 2);
    EXPECT_NEAR(x, y, 1e-12);
    EXPECT_NEAR(z, 0.0, 1e-12);
}

TEST(DataCondition, ClampPreservesDtypeKind) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o;
    o.location = DataLocation::Cell;
    o.names = {"tag"};  // Int32 on the MESHIO backend
    o.mode = ConditionMode::Clamp;
    o.lo = 0.0;
    o.hi = 25.0;
    Mesh out = data_condition(in, o);
    // Assert the KIND, not the width: NATIVE/KRATOS canonicalize Int32 -> Int64.
    EXPECT_FALSE(meshioplusplus::detail::is_float_dtype(out.CellData("tag", 0).Dtype()));
    EXPECT_DOUBLE_EQ(read_double(out.CellData("tag", 1), 0), 25.0);
}

TEST(DataCondition, NormalizeAlwaysProducesFloat) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o;
    o.location = DataLocation::Cell;
    o.names = {"tag"};
    o.mode = ConditionMode::Normalize;
    Mesh out = data_condition(in, o);
    EXPECT_TRUE(meshioplusplus::detail::is_float_dtype(out.CellData("tag", 0).Dtype()));
}

TEST(DataCondition, CellDataStatisticsAreJointAcrossBlocks) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o;
    o.location = DataLocation::Cell;
    o.names = {"mat"};
    o.mode = ConditionMode::Normalize;
    Mesh out = data_condition(in, o);
    // mat spans both blocks: {1,2} and {3}. Joint min 1, max 3.
    EXPECT_NEAR(read_double(out.CellData("mat", 0), 0), 0.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("mat", 0), 1), 0.5, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("mat", 1), 0), 1.0, 1e-12);
    ASSERT_EQ(out.CellDataNumBlocks("mat"), out.NumCellBlocks());
}

TEST(DataCondition, NaNIsExcludedFromReductionsAndPassedThrough) {
    Mesh in = mt::data_mesh();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    in.AddPointData("n", mt::data_array({0.0, nan, 4.0, 2.0, 1.0, 3.0}));
    DataConditionOptions o = point_opts(ConditionMode::Normalize);
    o.names = {"n"};
    Mesh out = data_condition(in, o);
    // min 0, max 4 (NaN ignored) -> 4 maps to 1.
    EXPECT_NEAR(read_double(out.PointData("n"), 0), 0.0, 1e-12);
    EXPECT_NEAR(read_double(out.PointData("n"), 2), 1.0, 1e-12);
    // Ignore policy passes the NaN through.
    EXPECT_TRUE(std::isnan(read_double(out.PointData("n"), 1)));
}

TEST(DataCondition, NaNReplacePolicy) {
    Mesh in = mt::data_mesh();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    in.AddPointData("n", mt::data_array({0.0, nan, 4.0, 2.0, 1.0, 3.0}));
    DataConditionOptions o = point_opts(ConditionMode::Normalize);
    o.names = {"n"};
    o.nan_policy = NanPolicy::Replace;
    o.nan_replacement = -1.0;
    Mesh out = data_condition(in, o);
    EXPECT_DOUBLE_EQ(read_double(out.PointData("n"), 1), -1.0);
}

TEST(DataCondition, NaNFailPolicyThrows) {
    Mesh in = mt::data_mesh();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    in.AddPointData("n", mt::data_array({0.0, nan, 4.0, 2.0, 1.0, 3.0}));
    DataConditionOptions o = point_opts(ConditionMode::Normalize);
    o.names = {"n"};
    o.nan_policy = NanPolicy::Fail;
    EXPECT_THROW(data_condition(in, o), std::invalid_argument);
}

TEST(DataCondition, SuffixLeavesTheOriginalAlone) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o = point_opts(ConditionMode::Normalize);
    o.names = {"T"};
    o.suffix = "_n";
    Mesh out = data_condition(in, o);
    ASSERT_TRUE(out.HasPointData("T_n"));
    EXPECT_NEAR(read_double(out.PointData("T_n"), 5), 1.0, 1e-12);
    // Original values untouched.
    EXPECT_DOUBLE_EQ(read_double(out.PointData("T"), 5), 12.0);
}

TEST(DataCondition, InvertedBoundsThrow) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o = point_opts(ConditionMode::Clamp);
    o.names = {"T"};
    o.lo = 10.0;
    o.hi = 1.0;
    EXPECT_THROW(data_condition(in, o), std::invalid_argument);
}

TEST(DataCondition, UnknownNameThrows) {
    Mesh in = mt::data_mesh();
    DataConditionOptions o = point_opts(ConditionMode::Clamp);
    o.names = {"nope"};
    EXPECT_THROW(data_condition(in, o), std::invalid_argument);
}

}  // namespace
