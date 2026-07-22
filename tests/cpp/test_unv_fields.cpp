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
// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/dex.hpp"
#include "meshioplusplus/formats/ip.hpp"
#include "meshioplusplus/formats/mff.hpp"
#include "meshioplusplus/formats/unv.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

namespace {

NDArray col(std::vector<double> v) {
    NDArray a(DType::Float64, {v.size()});
    for (std::size_t i = 0; i < v.size(); ++i)
        a.As<double>()[i] = v[i];
    return a;
}

NDArray mat(std::size_t rows, std::size_t cols, std::vector<double> v) {
    NDArray a(DType::Float64, {rows, cols});
    for (std::size_t i = 0; i < v.size(); ++i)
        a.As<double>()[i] = v[i];
    return a;
}

// Build tri_mesh (4 points, 2 triangles) with node + element fields.
Mesh field_mesh() {
    Mesh m = mt::tri_mesh();
    m.AddPointData("temp", col({1.0, 2.0, 3.0, 4.0}));
    m.AddPointData("disp", mat(4, 3, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    std::vector<NDArray> stress;
    stress.push_back(mat(2, 6, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    m.AddCellData("stress", std::move(stress));
    return m;
}

void expect_fields(const Mesh& m) {
    ASSERT_TRUE(m.HasPointData("temp"));
    ASSERT_TRUE(m.HasPointData("disp"));
    ASSERT_TRUE(m.HasCellData("stress"));
    const NDArray& temp = m.PointData("temp");
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(temp, i), double(i + 1));
    const NDArray& disp = m.PointData("disp");
    for (std::size_t i = 0; i < 12; ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(disp, i), double(i));
    const NDArray& stress = m.CellData("stress", 0);
    for (std::size_t i = 0; i < 12; ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(stress, i), double(i));
}

std::string tmp(const std::string& suffix) {
    return mt::temp_path(suffix);
}

}  // namespace

TEST(UnvField, Dataset2414) {
    std::string p = tmp(".unv");
    meshioplusplus::write_unv(p, field_mesh());
    expect_fields(meshioplusplus::read_unv(p));
    std::filesystem::remove(p);
}

TEST(UnvField, CodeAster5557) {
    std::string p = tmp(".unv");
    meshioplusplus::write_unv(p, field_mesh(), /*code_aster=*/true);
    expect_fields(meshioplusplus::read_unv(p));
    std::filesystem::remove(p);
}

TEST(UnvField, Wedge15) {
    Mesh m = mt::make_mesh({{0, 0, 0},
                            {1, 0, 0},
                            {0, 1, 0},
                            {0, 0, 1},
                            {1, 0, 1},
                            {0, 1, 1},
                            {0.5, 0, 0},
                            {0.5, 0.5, 0},
                            {0, 0.5, 0},
                            {0, 0, 0.5},
                            {1, 0, 0.5},
                            {0, 1, 0.5},
                            {0.5, 0, 1},
                            {0.5, 0.5, 1},
                            {0, 0.5, 1}},
                           "wedge15", {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}});
    std::string p = tmp(".unv");
    meshioplusplus::write_unv(p, m);
    Mesh out = meshioplusplus::read_unv(p);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "wedge15");
    std::filesystem::remove(p);
}

TEST(UnvGroup, PointAndCellSets) {
    Mesh m = mt::tri_mesh();
    meshioplusplus::UnvInfo in;
    in.mPointSets["corners"] = {0, 2};
    in.mCellSets["all"] = {{0, 1}};
    std::string p = tmp(".unv");
    meshioplusplus::write_unv(p, m, in);
    meshioplusplus::UnvInfo out;
    meshioplusplus::read_unv(p, out);
    ASSERT_EQ(out.mPointSets.count("corners"), 1u);
    EXPECT_EQ(out.mPointSets["corners"], (std::vector<std::int64_t>{0, 2}));
    ASSERT_EQ(out.mCellSets.count("all"), 1u);
    ASSERT_EQ(out.mCellSets["all"].size(), 1u);
    EXPECT_EQ(out.mCellSets["all"][0], (std::vector<std::int64_t>{0, 1}));
    std::filesystem::remove(p);
}

TEST(Mff, ValuesRoundtrip) {
    Mesh m;
    m.AssignPoints(NDArray(DType::Float64, {5, 0}));
    m.AddPointData("mff:field", col({1.5, -2.25, 3.0, 4.0, 5.0}));
    std::string p = tmp(".mff");
    meshioplusplus::write_mff(p, m);
    Mesh out = meshioplusplus::read_mff(p);
    ASSERT_TRUE(out.HasPointData("mff:field"));
    const NDArray& v = out.PointData("mff:field");
    std::vector<double> want = {1.5, -2.25, 3.0, 4.0, 5.0};
    for (std::size_t i = 0; i < want.size(); ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(v, i), want[i]);
    std::filesystem::remove(p);
}

TEST(Dex, CoordsAndValues) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}));
    m.AddPointData("mGradT", mat(3, 3, {0, 1, 2, 3, 4, 5, 6, 7, 8}));
    std::string p = tmp(".dex");
    meshioplusplus::write_dex(p, m);
    Mesh out = meshioplusplus::read_dex(p);
    ASSERT_EQ(out.NumPoints(), 3u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out.Points(), 3), 1.0);  // point 1, x
    const NDArray& v = out.PointData("mGradT");
    for (std::size_t i = 0; i < 9; ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(v, i), double(i));
    std::filesystem::remove(p);
}

TEST(Ip, CoordsAndFields) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0}, {1, 0}, {0, 1}, {1, 1}}));
    m.AddPointData("pressure", col({10, 20, 30, 40}));
    m.AddPointData("x-velocity", col({1, 2, 3, 4}));
    std::string p = tmp(".ip");
    meshioplusplus::write_ip(p, m);
    Mesh out = meshioplusplus::read_ip(p);
    ASSERT_EQ(out.NumPoints(), 4u);
    const NDArray& pr = out.PointData("pressure");
    std::vector<double> want = {10, 20, 30, 40};
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(pr, i), want[i]);
    std::filesystem::remove(p);
}
