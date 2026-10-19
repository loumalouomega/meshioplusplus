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

// Tecplot SZL (.szplt) through TecIO: synthetic files written with TecIO's own
// writer must read as the same data in ASCII does.
#ifdef MESHIOPLUSPLUS_HAS_TECIO

// System includes
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

#include "TECIO.h"

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/szplt.hpp"
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;

// Two transient quad zones (strand 1, t = 0 and 1): the second shares X, Y and
// the connectivity of the first; P nodal, Q cell-centred.
void szplt_write(const std::string& rPath) {
    void* h = nullptr;
    ASSERT_EQ(tecFileWriterOpen(rPath.c_str(), "t", "X,Y,P,Q", 1, 0, 2, nullptr, &h), 0);
    std::vector<int32_t> types(4, 2), loc = {1, 1, 1, 0}, passive(4, 0);
    for (int s = 0; s < 2; ++s) {
        std::vector<int32_t> share(4, 0);
        if (s == 1)
            share[0] = share[1] = 1;
        int32_t zone = 0;
        ASSERT_EQ(
            tecZoneCreateFE(h, s == 0 ? "first" : "second", 3, 6, 2, types.data(), share.data(),
                            loc.data(), passive.data(), s == 1 ? 1 : 0, 0, 0, &zone),
            0);
        ASSERT_EQ(tecZoneSetUnsteadyOptions(h, zone, static_cast<double>(s), 1), 0);
        if (s == 0) {
            const std::vector<double> x = {0, 1, 2, 0, 1, 2}, y = {0, 0, 0, 1, 1, 1};
            ASSERT_EQ(tecZoneVarWriteDoubleValues(h, zone, 1, 0, 6, x.data()), 0);
            ASSERT_EQ(tecZoneVarWriteDoubleValues(h, zone, 2, 0, 6, y.data()), 0);
        }
        std::vector<double> p;
        for (int n = 0; n < 6; ++n)
            p.push_back(n + 10.0 * s);
        const std::vector<double> q = {1.0 + s, 2.0 + s};
        ASSERT_EQ(tecZoneVarWriteDoubleValues(h, zone, 3, 0, 6, p.data()), 0);
        ASSERT_EQ(tecZoneVarWriteDoubleValues(h, zone, 4, 0, 2, q.data()), 0);
        if (s == 0) {
            const std::vector<int32_t> nodes = {1, 2, 5, 4, 2, 3, 6, 5};
            ASSERT_EQ(tecZoneNodeMapWrite32(h, zone, 0, 1, 8, nodes.data()), 0);
        }
    }
    ASSERT_EQ(tecFileWriterClose(&h), 0);
}

const char* const kAscii = R"(VARIABLES = "X" "Y" "P" "Q"
ZONE T="first", ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK, NODES=6, ELEMENTS=2, VARLOCATION=([4]=CELLCENTERED), SOLUTIONTIME=0, STRANDID=1
0 1 2 0 1 2
0 0 0 1 1 1
0 1 2 3 4 5
1 2
1 2 5 4
2 3 6 5
ZONE T="second", ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK, NODES=6, ELEMENTS=2, VARLOCATION=([4]=CELLCENTERED), VARSHARELIST=([1,2]=1), CONNECTIVITYSHAREZONE=1, SOLUTIONTIME=1, STRANDID=1
10 11 12 13 14 15
2 3
)";

TEST(Szplt, ReadsAsTheSameDataInAscii) {
    const std::string szl = mt::temp_path(".szplt");
    const std::string dat = mt::temp_path(".dat");
    szplt_write(szl);
    {
        std::ofstream out(dat);
        out << kAscii;
    }
    EXPECT_EQ(meshioplusplus::szplt_time_values(szl), (std::vector<double>{0.0, 1.0}));
    for (int step = 0; step < 2; ++step) {
        ReadOptions opts;
        opts.mTimeStep = step;
        const Mesh a = meshioplusplus::read_szplt(szl, opts);
        const Mesh b = meshioplusplus::read_tecplot(dat, opts);
        mt::expect_mesh_eq(a, b, 0.0);
        const auto* pa = a.PointData("P").As<double>();
        const auto* pb = b.PointData("P").As<double>();
        for (std::size_t i = 0; i < 6; ++i)
            EXPECT_EQ(pa[i], pb[i]);
        EXPECT_EQ(a.CellData("Q", 0).As<double>()[1], 2.0 + step);
    }
    ReadOptions opts;
    opts.mTimeStep = 2;
    EXPECT_THROW(meshioplusplus::read_szplt(szl, opts), ReadError);
    EXPECT_EQ(meshioplusplus::sniff_format(szl), "szplt");
    EXPECT_EQ(meshioplusplus::resolve_format(szl, ""), "szplt");
    std::filesystem::remove(szl);
    std::filesystem::remove(dat);
}

TEST(SzpltWriter, StaticFileAndSeriesReadBack) {
    // Two blocks (a tetrahedron and a triangle) sharing nodes: each zone gets the points
    // its cells use.
    Mesh mesh = mt::tet_mesh();
    {
        meshioplusplus::NDArray tri(meshioplusplus::DType::Int64, {std::size_t{1}, std::size_t{3}});
        tri.As<std::int64_t>()[0] = 0;
        tri.As<std::int64_t>()[1] = 1;
        tri.As<std::int64_t>()[2] = 2;
        mesh.AddCellBlock("triangle", std::move(tri));
        meshioplusplus::NDArray p(meshioplusplus::DType::Float64, {mesh.NumPoints()});
        for (std::size_t i = 0; i < mesh.NumPoints(); ++i)
            p.As<double>()[i] = 3.0 * static_cast<double>(i);
        mesh.AddPointData("P", std::move(p));
    }
    const std::string path = mt::temp_path(".szplt");
    meshioplusplus::write_szplt(path, mesh);
    const Mesh back = meshioplusplus::read_szplt(path);
    ASSERT_EQ(back.NumCellBlocks(), 2u);
    EXPECT_EQ(back.Cells(0).Type(), "tetra");
    EXPECT_EQ(back.Cells(1).Type(), "triangle");
    EXPECT_EQ(back.NumPoints(), 8u);  // 5 for the tetrahedra, 3 repeated for the triangle
    EXPECT_EQ(back.PointData("P").As<double>()[6], 3.0);  // the triangle's second node

    // A lone block keeps every point; a series shares the first step's mesh.
    const Mesh tet = mt::tet_mesh();
    const std::string series = mt::temp_path(".szplt");
    {
        meshioplusplus::SzpltSeriesWriter writer(series);
        writer.Write(0.0, tet);
        writer.Write(0.5, tet);
        writer.Finalize();
        EXPECT_EQ(writer.NumSteps(), 2u);
    }
    EXPECT_EQ(meshioplusplus::szplt_time_values(series), (std::vector<double>{0.0, 0.5}));
    ReadOptions opts;
    opts.mTimeStep = 1;
    const Mesh second = meshioplusplus::read_szplt(series, opts);
    EXPECT_EQ(second.NumPoints(), tet.NumPoints());
    EXPECT_EQ(second.Cells(0).Type(), "tetra");

    meshioplusplus::SzpltSeriesWriter bad(mt::temp_path(".szplt"));
    bad.Write(0.0, tet);
    Mesh moved = mt::tet_mesh();
    {
        meshioplusplus::NDArray pts(meshioplusplus::DType::Float64,
                                    {moved.NumPoints(), moved.PointDim()});
        const auto* src = tet.Points().As<double>();
        for (std::size_t i = 0; i < pts.Size(); ++i)
            pts.As<double>()[i] = src[i];
        pts.As<double>()[0] = 9.0;
        moved.AssignPoints(std::move(pts));
    }
    EXPECT_THROW(bad.Write(1.0, moved), meshioplusplus::WriteError);
}

// The type TecIO stores a variable in, read through its own reader.
std::vector<int32_t> szplt_variable_types(const std::string& rPath) {
    void* h = nullptr;
    EXPECT_EQ(tecFileReaderOpen(rPath.c_str(), &h), 0);
    int32_t count = 0;
    EXPECT_EQ(tecDataSetGetNumVars(h, &count), 0);
    std::vector<int32_t> out;
    for (int32_t v = 1; v <= count; ++v) {
        int32_t type = 0;
        EXPECT_EQ(tecZoneVarGetType(h, 1, v, &type), 0);
        out.push_back(type);
    }
    tecFileReaderClose(&h);
    return out;
}

TEST(SzpltWriter, VariablesKeepTheTypeThatHoldsTheirDType) {
    Mesh mesh = mt::tet_mesh();
    const std::size_t n = mesh.NumPoints();
    meshioplusplus::NDArray f32(meshioplusplus::DType::Float32, {n});
    meshioplusplus::NDArray i16(meshioplusplus::DType::Int16, {n});
    meshioplusplus::NDArray u8(meshioplusplus::DType::UInt8, {n});
    meshioplusplus::NDArray i64(meshioplusplus::DType::Int64, {n});
    for (std::size_t i = 0; i < n; ++i) {
        f32.As<float>()[i] = 0.25f * static_cast<float>(i);
        i16.As<std::int16_t>()[i] = static_cast<std::int16_t>(-3 * static_cast<int>(i));
        u8.As<std::uint8_t>()[i] = static_cast<std::uint8_t>(200 + i);
        i64.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    }
    mesh.AddPointData("a_f32", std::move(f32));
    mesh.AddPointData("b_i16", std::move(i16));
    mesh.AddPointData("c_u8", std::move(u8));
    mesh.AddPointData("d_i64", std::move(i64));
    const std::string path = mt::temp_path(".szplt");
    meshioplusplus::write_szplt(path, mesh);
    // X, Y, Z in the points' double precision, then the four arrays in their own types
    EXPECT_EQ(szplt_variable_types(path), (std::vector<int32_t>{2, 2, 2, 1, 4, 5, 2}));
    const Mesh back = meshioplusplus::read_szplt(path);
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQ(back.PointData("a_f32").As<double>()[i], 0.25 * static_cast<double>(i));
        EXPECT_EQ(back.PointData("b_i16").As<double>()[i], -3.0 * static_cast<double>(i));
        EXPECT_EQ(back.PointData("c_u8").As<double>()[i], 200.0 + static_cast<double>(i));
    }
    std::filesystem::remove(path);
}

TEST(SzpltWriter, CellsWithPartsAreWrittenAsPartitions) {
    // A row of four cubes along x, in parts 0 1 0 1: every cell keeps its geometry.
    Mesh mesh;
    meshioplusplus::NDArray pts(meshioplusplus::DType::Float64, {std::size_t{20}, std::size_t{3}});
    const double yz[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    for (std::size_t i = 0; i < 5; ++i)
        for (std::size_t c = 0; c < 4; ++c) {
            double* row = pts.As<double>() + (i * 4 + c) * 3;
            row[0] = static_cast<double>(i);
            row[1] = yz[c][0];
            row[2] = yz[c][1];
        }
    mesh.AssignPoints(std::move(pts));
    meshioplusplus::NDArray conn(meshioplusplus::DType::Int64, {std::size_t{4}, std::size_t{8}});
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t c = 0; c < 8; ++c)
            conn.As<std::int64_t>()[i * 8 + c] = static_cast<std::int64_t>(i * 4 + c);
    mesh.AddCellBlock("hexahedron", std::move(conn));
    meshioplusplus::NDArray part(meshioplusplus::DType::Int32, {std::size_t{4}});
    meshioplusplus::NDArray rho(meshioplusplus::DType::Float64, {std::size_t{4}});
    for (std::size_t i = 0; i < 4; ++i) {
        part.As<std::int32_t>()[i] = static_cast<std::int32_t>(i % 2);
        rho.As<double>()[i] = 1.0 + static_cast<double>(i);
    }
    mesh.AddCellData("partition:part", {std::move(part)});
    mesh.AddCellData("rho", {std::move(rho)});
    const std::string path = mt::temp_path(".szplt");
    meshioplusplus::write_szplt(path, mesh);
    const Mesh back = meshioplusplus::read_szplt(path);
    ASSERT_EQ(back.NumPoints(), 20u);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    ASSERT_EQ(back.Cells(0).NumCells(), 4u);
    // The cells come back grouped by part (0 2 | 1 3), each one a cube at x = i .. i + 1.
    const std::int64_t* cells = back.Cells(0).Conn().As<std::int64_t>();
    const double* rho_back = back.CellData("rho", 0).As<double>();
    const double expected[4] = {1.0, 3.0, 2.0, 4.0};
    for (std::size_t c = 0; c < 4; ++c) {
        EXPECT_EQ(rho_back[c], expected[c]);
        const double x0 = expected[c] - 1.0;
        double lo = 1e9, hi = -1e9;
        for (std::size_t j = 0; j < 8; ++j) {
            const double x =
                back.Points().As<double>()[static_cast<std::size_t>(cells[c * 8 + j]) * 3];
            lo = std::min(lo, x);
            hi = std::max(hi, x);
        }
        EXPECT_EQ(lo, x0);
        EXPECT_EQ(hi, x0 + 1.0);
    }
    std::filesystem::remove(path);
}

TEST(Szplt, MissingFileIsAReadError) {
    EXPECT_THROW(meshioplusplus::read_szplt("/nonexistent/x.szplt"), ReadError);
}

}  // namespace

#endif  // MESHIOPLUSPLUS_HAS_TECIO
