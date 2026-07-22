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
// Tests for point <-> cell data averaging.

// System includes
#include <cmath>
#include <stdexcept>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/data_average.hpp"

namespace {

using meshioplusplus::cell_data_to_point_data;
using meshioplusplus::CellPointWeight;
using meshioplusplus::DataAverageOptions;
using meshioplusplus::Mesh;
using meshioplusplus::point_data_to_cell_data;
using meshioplusplus::detail::read_double;

TEST(DataAverage, PointToCellIsTheMeanOfTheCellsNodes) {
    Mesh in = mt::data_mesh();
    DataAverageOptions opts;
    opts.names = {"T"};
    Mesh out = point_data_to_cell_data(in, opts);

    ASSERT_TRUE(out.HasCellData("T"));
    ASSERT_EQ(out.CellDataNumBlocks("T"), out.NumCellBlocks());

    // T = {0, 1, 11, 10, 2, 12}; triangle cells {0,1,2} and {0,2,3}.
    EXPECT_NEAR(read_double(out.CellData("T", 0), 0), (0.0 + 1.0 + 11.0) / 3.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("T", 0), 1), (0.0 + 11.0 + 10.0) / 3.0, 1e-12);
    // quad cell {1,4,5,2}.
    EXPECT_NEAR(read_double(out.CellData("T", 1), 0), (1.0 + 2.0 + 12.0 + 11.0) / 4.0, 1e-12);

    mt::expect_same_geometry(in, out);
}

TEST(DataAverage, PointToCellHandlesVectorArrays) {
    Mesh in = mt::data_mesh();
    DataAverageOptions opts;
    opts.names = {"v"};
    Mesh out = point_data_to_cell_data(in, opts);
    ASSERT_TRUE(out.HasCellData("v"));
    // 3 components preserved.
    ASSERT_EQ(out.CellData("v", 0).Shape().size(), 2u);
    EXPECT_EQ(out.CellData("v", 0).Shape()[1], 3u);
    // triangle {0,1,2}: v = (1,0,0), (0,1,0), (0,0,1) -> mean (1/3, 1/3, 1/3)
    EXPECT_NEAR(read_double(out.CellData("v", 0), 0), 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("v", 0), 1), 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(read_double(out.CellData("v", 0), 2), 1.0 / 3.0, 1e-12);
}

TEST(DataAverage, CellToPointIsTheIncidentCellMean) {
    Mesh in = mt::data_mesh();
    DataAverageOptions opts;
    opts.names = {"mat"};
    Mesh out = cell_data_to_point_data(in, opts);
    ASSERT_TRUE(out.HasPointData("mat"));

    // mat: triangle cells = {1, 2}, quad cell = {3}.
    // Point 0 belongs to triangles 0 and 1 -> (1+2)/2 = 1.5
    EXPECT_NEAR(read_double(out.PointData("mat"), 0), 1.5, 1e-12);
    // Point 4 belongs only to the quad -> 3
    EXPECT_NEAR(read_double(out.PointData("mat"), 4), 3.0, 1e-12);
    // Point 2 belongs to both triangles and the quad -> (1+2+3)/3 = 2
    EXPECT_NEAR(read_double(out.PointData("mat"), 2), 2.0, 1e-12);
    // Point 3 belongs only to triangle 1 -> 2
    EXPECT_NEAR(read_double(out.PointData("mat"), 3), 2.0, 1e-12);

    mt::expect_same_geometry(in, out);
}

TEST(DataAverage, WeightedDiffersFromUniformOnANonUniformMesh) {
    // The two triangles each have area 1/2; the quad has area 1. Point 2 is
    // shared by all three, so the measure-weighted mean differs from the plain
    // one.
    Mesh in = mt::data_mesh();
    DataAverageOptions uniform;
    uniform.names = {"mat"};
    DataAverageOptions weighted;
    weighted.names = {"mat"};
    weighted.weight = CellPointWeight::Measure;

    Mesh u = cell_data_to_point_data(in, uniform);
    Mesh w = cell_data_to_point_data(in, weighted);

    const double uv = read_double(u.PointData("mat"), 2);
    const double wv = read_double(w.PointData("mat"), 2);
    EXPECT_NEAR(uv, (1.0 + 2.0 + 3.0) / 3.0, 1e-12);
    // (0.5*1 + 0.5*2 + 1*3) / (0.5 + 0.5 + 1) = 4.5 / 2 = 2.25
    EXPECT_NEAR(wv, 2.25, 1e-12);
    EXPECT_NE(uv, wv);
}

TEST(DataAverage, ConstantFieldSurvivesTheRoundTripExactly) {
    Mesh in = mt::data_mesh();
    // Overwrite T with a constant.
    in.AddPointData("T", mt::data_array({7.0, 7.0, 7.0, 7.0, 7.0, 7.0}));

    DataAverageOptions to_cell;
    to_cell.names = {"T"};
    Mesh mid = point_data_to_cell_data(in, to_cell);

    DataAverageOptions to_point;
    to_point.names = {"T"};
    Mesh back = cell_data_to_point_data(mid, to_point);

    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_NEAR(read_double(back.PointData("T"), i), 7.0, 1e-12);
}

TEST(DataAverage, OutputIsAlwaysFloat64EvenForIntegerInput) {
    Mesh in = mt::data_mesh();
    DataAverageOptions opts;
    opts.names = {"tag"};  // Int32 cell_data
    Mesh out = cell_data_to_point_data(in, opts);
    ASSERT_TRUE(out.HasPointData("tag"));
    // A mean is not an integer: assert the dtype KIND, since the NATIVE and
    // KRATOS backends canonicalize widths on ingest.
    EXPECT_TRUE(meshioplusplus::detail::is_float_dtype(out.PointData("tag").Dtype()));
}

TEST(DataAverage, SuffixNamesTheOutputWithoutTouchingTheInput) {
    Mesh in = mt::data_mesh();
    DataAverageOptions opts;
    opts.names = {"T"};
    opts.suffix = "_c";
    Mesh out = point_data_to_cell_data(in, opts);
    EXPECT_TRUE(out.HasCellData("T_c"));
    EXPECT_TRUE(out.HasPointData("T"));  // source untouched
    EXPECT_FALSE(out.HasCellData("T"));
}

TEST(DataAverage, IsolatedPointYieldsNaN) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {5, 5, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    std::vector<meshioplusplus::NDArray> mat;
    mat.push_back(mt::data_array({4.0}));
    m.AddCellData("mat", std::move(mat));

    Mesh out = cell_data_to_point_data(m);
    EXPECT_NEAR(read_double(out.PointData("mat"), 0), 4.0, 1e-12);
    // Point 3 touches no cell.
    EXPECT_TRUE(std::isnan(read_double(out.PointData("mat"), 3)));
}

TEST(DataAverage, UnknownNameThrows) {
    Mesh in = mt::data_mesh();
    DataAverageOptions opts;
    opts.names = {"nope"};
    EXPECT_THROW(point_data_to_cell_data(in, opts), std::invalid_argument);
}

TEST(DataAverage, DeterministicAcrossRuns) {
    // The cell->point accumulation is deliberately serial because FP addition
    // is not associative; two runs must agree bit for bit.
    Mesh in = mt::data_mesh();
    Mesh a = cell_data_to_point_data(in);
    Mesh b = cell_data_to_point_data(in);
    ASSERT_EQ(a.PointData("mat").Size(), b.PointData("mat").Size());
    for (std::size_t i = 0; i < a.PointData("mat").Size(); ++i)
        EXPECT_EQ(read_double(a.PointData("mat"), i), read_double(b.PointData("mat"), i));
}

}  // namespace
