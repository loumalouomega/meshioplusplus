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
 * @file test_patran.cpp
 * @brief Patran 2 neutral reader/writer: hand-written hex20 and wedge15 cards in
 *        Patran's node order, components and the property fallback, the error
 *        paths, and a round trip through the writer.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <array>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/patran.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::ReadError;
using meshioplusplus::RegionKind;
namespace detail = meshioplusplus::detail;

std::string write_file(const std::string& rBody) {
    const std::string path = mt::temp_path(".pat");
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

std::string read_text(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::string header(int It, long long Id, int Iv, int Kc) {
    char buf[96];
    std::snprintf(buf, sizeof(buf), "%2d%8lld%8d%8d%8d%8d%8d%8d%8d\n", It, Id, Iv, Kc, 0, 0, 0, 0,
                  0);
    return buf;
}

std::string node(long long Id, double X, double Y, double Z) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%16.9E%16.9E%16.9E\n", X, Y, Z);
    return header(1, Id, 0, 2) + buf + "1G       6       0       0  000000\n";
}

std::string element(long long Id, int Shape, int Pid, const std::vector<long long>& rNodes) {
    std::string out = header(2, Id, Shape, 1 + static_cast<int>((rNodes.size() + 9) / 10));
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%8zu%8d%8d%8d%16.9E%16.9E%16.9E\n", rNodes.size(), 0, Pid, 0,
                  0.0, 0.0, 0.0);
    out += buf;
    for (std::size_t k = 0; k < rNodes.size(); ++k) {
        std::snprintf(buf, sizeof(buf), "%8lld", rNodes[k]);
        out += buf;
        if (k % 10 == 9 || k + 1 == rNodes.size())
            out += '\n';
    }
    return out;
}

std::array<double, 3> point(const Mesh& rMesh, std::int64_t i) {
    std::array<double, 3> p{0, 0, 0};
    for (std::size_t d = 0; d < 3; ++d)
        p[d] = detail::read_double(rMesh.Points(), static_cast<std::size_t>(i) * 3 + d);
    return p;
}

std::vector<std::int64_t> entries(const Mesh& rMesh, std::size_t r) {
    const auto& reg = rMesh.Region(r);
    return {reg.Entries(), reg.Entries() + reg.NumEntries()};
}

// A unit hex20 in Patran order: corners 1-8, bottom ring (1-2, 2-3, 3-4, 4-1),
// vertical edges (1-5 ... 4-8), top ring (5-6 ... 8-5); node ids step by 10.
// Then a wedge15: corners, bottom ring, verticals, top ring. Components: the
// hex as "BLOCK" (component 7) plus its four bottom nodes; the wedge is in none.
std::string hex20_and_wedge15() {
    const double hex[20][3] = {{0, 0, 0},  {1, 0, 0},  {1, 1, 0},  {0, 1, 0},  {0, 0, 1},
                               {1, 0, 1},  {1, 1, 1},  {0, 1, 1},  {.5, 0, 0}, {1, .5, 0},
                               {.5, 1, 0}, {0, .5, 0}, {0, 0, .5}, {1, 0, .5}, {1, 1, .5},
                               {0, 1, .5}, {.5, 0, 1}, {1, .5, 1}, {.5, 1, 1}, {0, .5, 1}};
    const double wedge[15][3] = {{2, 0, 0},  {3, 0, 0},   {2, 1, 0},    {2, 0, 1},    {3, 0, 1},
                                 {2, 1, 1},  {2.5, 0, 0}, {2.5, .5, 0}, {2, .5, 0},   {2, 0, .5},
                                 {3, 0, .5}, {2, 1, .5},  {2.5, 0, 1},  {2.5, .5, 1}, {2, .5, 1}};
    std::string out = header(25, 0, 0, 1) + "hand-written\n";
    std::vector<long long> hex_ids, wedge_ids;
    for (int k = 0; k < 20; ++k) {
        hex_ids.push_back(10 * (k + 1));
        out += node(hex_ids.back(), hex[k][0], hex[k][1], hex[k][2]);
    }
    for (int k = 0; k < 15; ++k) {
        wedge_ids.push_back(1000 + k);
        out += node(wedge_ids.back(), wedge[k][0], wedge[k][1], wedge[k][2]);
    }
    out += element(55, 8, 3, hex_ids);
    out += element(66, 7, 4, wedge_ids);
    out += header(21, 7, 10, 2) + "BLOCK\n";
    out += "      12      55       5      10       5      20       5      30       5      40\n";
    out += header(99, 0, 0, 1);
    return out;
}

// meshio++ edge order: hexahedron20 bottom, top, verticals; wedge15 likewise.
const int kHexEdges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                              {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
const int kWedgeEdges[9][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5},
                               {5, 3}, {0, 3}, {1, 4}, {2, 5}};

template <std::size_t N>
void expect_midpoints(const Mesh& rMesh, std::size_t Block, std::size_t Corners,
                      const int (&rEdges)[N][2]) {
    const auto& conn = rMesh.Cells(Block).Conn();
    for (std::size_t e = 0; e < N; ++e) {
        const auto a = point(rMesh, detail::read_int(conn, rEdges[e][0]));
        const auto b = point(rMesh, detail::read_int(conn, rEdges[e][1]));
        const auto m = point(rMesh, detail::read_int(conn, Corners + e));
        for (int d = 0; d < 3; ++d)
            EXPECT_DOUBLE_EQ(m[d], 0.5 * (a[d] + b[d])) << "edge " << e;
    }
}

}  // namespace

TEST(Patran, ReadsHex20AndWedge15IntoMeshioOrder) {
    const Mesh mesh = meshioplusplus::read_patran(write_file(hex20_and_wedge15()));
    ASSERT_EQ(mesh.NumPoints(), 35u);
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    ASSERT_EQ(mesh.Cells(0).Type(), "hexahedron20");
    ASSERT_EQ(mesh.Cells(1).Type(), "wedge15");
    expect_midpoints(mesh, 0, 8, kHexEdges);
    expect_midpoints(mesh, 1, 6, kWedgeEdges);
    ASSERT_TRUE(mesh.HasCellData("patran:property"));
    EXPECT_EQ(detail::read_int(mesh.CellData("patran:property", 0), 0), 3);
    EXPECT_EQ(detail::read_int(mesh.CellData("patran:property", 1), 0), 4);
}

TEST(Patran, ComponentsBecomeRegionsAndPropertiesTheFallback) {
    const Mesh mesh = meshioplusplus::read_patran(write_file(hex20_and_wedge15()));
    const std::size_t block = mesh.FindRegion("BLOCK", RegionKind::Cell);
    const std::size_t nodes = mesh.FindRegion("BLOCK", RegionKind::Point);
    const std::size_t prop = mesh.FindRegion("property_4", RegionKind::Cell);
    ASSERT_NE(block, Mesh::npos);
    ASSERT_NE(nodes, Mesh::npos);
    ASSERT_NE(prop, Mesh::npos);
    EXPECT_EQ(mesh.FindRegion("property_3", RegionKind::Cell), Mesh::npos);
    EXPECT_EQ(entries(mesh, block), (std::vector<std::int64_t>{0}));
    EXPECT_EQ(mesh.Region(block).mTag, 7);
    EXPECT_EQ(mesh.Region(block).mDim, 3);
    EXPECT_EQ(entries(mesh, nodes), (std::vector<std::int64_t>{0, 1, 2, 3}));
    EXPECT_EQ(entries(mesh, prop), (std::vector<std::int64_t>{1}));
    EXPECT_EQ(mesh.Region(prop).mTag, 4);
}

TEST(Patran, RoundTripsThroughItsOwnWriter) {
    const Mesh mesh = meshioplusplus::read_patran(write_file(hex20_and_wedge15()));
    const std::string out = mt::temp_path(".pat");
    meshioplusplus::write_patran(out, mesh);
    const std::string text = read_text(out);
    EXPECT_EQ(text.rfind("25       0       0       1", 0), 0u);
    EXPECT_NE(text.find("\n26       0       0       1      35       2       0       2"),
              std::string::npos);
    EXPECT_NE(text.find("\n99       0       0       1"), std::string::npos);
    const Mesh back = meshioplusplus::read_patran(out);
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    ASSERT_EQ(back.NumCellBlocks(), 2u);
    for (std::size_t b = 0; b < 2; ++b) {
        const auto& x = mesh.Cells(b).Conn();
        const auto& y = back.Cells(b).Conn();
        for (std::size_t k = 0; k < x.Size(); ++k)
            EXPECT_EQ(detail::read_int(x, k), detail::read_int(y, k));
    }
    for (std::size_t k = 0; k < mesh.NumPoints() * 3; ++k)
        EXPECT_DOUBLE_EQ(detail::read_double(mesh.Points(), k),
                         detail::read_double(back.Points(), k));
    ASSERT_EQ(back.NumRegions(), mesh.NumRegions());
    for (std::size_t r = 0; r < mesh.NumRegions(); ++r) {
        EXPECT_EQ(back.Region(r).mName, mesh.Region(r).mName);
        EXPECT_EQ(back.Region(r).mTag, mesh.Region(r).mTag);
        EXPECT_EQ(entries(back, r), entries(mesh, r));
    }
}

TEST(Patran, ErrorsNameTheCulprit) {
    const std::string n1 = node(1, 0, 0, 0);
    EXPECT_THROW(meshioplusplus::read_patran(write_file(n1 + element(7, 2, 1, {1, 9}))), ReadError);
    EXPECT_THROW(meshioplusplus::read_patran(write_file(n1 + n1)), ReadError);
    EXPECT_THROW(meshioplusplus::read_patran(write_file(header(1, 1, 0, 5))), ReadError);
    try {
        meshioplusplus::read_patran(write_file(n1 + element(7, 2, 1, {1, 9})));
        FAIL() << "no error";
    } catch (const ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("undefined node 9"), std::string::npos) << e.what();
    }
}

TEST(Patran, UnknownShapesArePassedOver) {
    const std::string body = node(1, 0, 0, 0) + node(2, 1, 0, 0) + element(5, 2, 1, {1, 2, 1, 2}) +
                             element(6, 2, 1, {1, 2}) + header(99, 0, 0, 1);
    const Mesh mesh = meshioplusplus::read_patran(write_file(body));
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "line");
}

TEST(Patran, WriterDropsCellTypesWithoutAShape) {
    Mesh mesh;
    NDArray pts(meshioplusplus::DType::Float64, {3, 3});
    double* p = pts.As<double>();
    const double xyz[9] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    std::copy(xyz, xyz + 9, p);
    mesh.AssignPoints(std::move(pts));
    NDArray vtx(meshioplusplus::DType::Int64, {1, 1});
    vtx.As<std::int64_t>()[0] = 0;
    mesh.AddCellBlock("vertex", std::move(vtx));
    NDArray tri(meshioplusplus::DType::Int64, {1, 3});
    for (int k = 0; k < 3; ++k)
        tri.As<std::int64_t>()[k] = k;
    mesh.AddCellBlock("triangle", std::move(tri));
    const std::string out = mt::temp_path(".pat");
    meshioplusplus::write_patran(out, mesh);
    const Mesh back = meshioplusplus::read_patran(out);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "triangle");
}

TEST(Patran, Quad9AndTriangle7RoundTrip) {
    // Patran's QUAD9 and TRI7 are corners, mid-edges, then the centre: the
    // meshio++ order, so the slots pass through unchanged.
    Mesh mesh;
    NDArray pts(meshioplusplus::DType::Float64, {9, 3});
    const double xy[9][2] = {{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 0},
                             {2, 1}, {1, 2}, {0, 1}, {1, 1}};
    for (int k = 0; k < 9; ++k) {
        pts.As<double>()[3 * k] = xy[k][0];
        pts.As<double>()[3 * k + 1] = xy[k][1];
        pts.As<double>()[3 * k + 2] = 0.0;
    }
    mesh.AssignPoints(std::move(pts));
    NDArray q9(meshioplusplus::DType::Int64, {1, 9});
    for (int k = 0; k < 9; ++k)
        q9.As<std::int64_t>()[k] = k;
    mesh.AddCellBlock("quad9", std::move(q9));
    NDArray t7(meshioplusplus::DType::Int64, {1, 7});
    const std::int64_t tri[7] = {0, 1, 2, 4, 5, 8, 7};  // slots only; geometry unchecked
    std::copy(tri, tri + 7, t7.As<std::int64_t>());
    mesh.AddCellBlock("triangle7", std::move(t7));
    const std::string out = mt::temp_path(".pat");
    meshioplusplus::write_patran(out, mesh);
    const Mesh back = meshioplusplus::read_patran(out);
    ASSERT_EQ(back.NumCellBlocks(), 2u);
    EXPECT_EQ(back.Cells(0).Type(), "quad9");
    EXPECT_EQ(back.Cells(1).Type(), "triangle7");
    for (int k = 0; k < 9; ++k)
        EXPECT_EQ(
            meshioplusplus::detail::read_int(back.Cells(0).Conn(), static_cast<std::size_t>(k)), k);
    for (int k = 0; k < 7; ++k)
        EXPECT_EQ(
            meshioplusplus::detail::read_int(back.Cells(1).Conn(), static_cast<std::size_t>(k)),
            tri[k]);
}

namespace {

// One quad (element 7, nodes 1-4) with a force (07), a constraint (08), a node
// and an element temperature (10, 11) and a pressure on its edge 2 (06).
std::string quad_with_loads() {
    std::string s = header(25, 0, 0, 1) + "loads\n";
    s += node(1, 0, 0, 0) + node(2, 1, 0, 0) + node(3, 1, 1, 0) + node(4, 0, 1, 0);
    s += element(7, 4, 1, {1, 2, 3, 4});
    s += header(7, 3, 2, 2) + "       9001100\n" + " 5.000000000E+00-1.000000000E+00\n";
    s += header(8, 1, 1, 2) + "       0111000\n" +
         " 0.000000000E+00 0.000000000E+00 0.000000000E+00\n";
    char buf[96];
    std::snprintf(buf, sizeof(buf), "%2d%8d%8d%8d%8d%8d%8d%8d%8d\n", 10, 2, 4, 1, 1, 0, 0, 0, 0);
    s += std::string(buf) + " 3.000000000E+02\n";
    std::snprintf(buf, sizeof(buf), "%2d%8d%8d%8d%8d%8d%8d%8d%8d\n", 11, 7, 4, 1, 1, 0, 0, 0, 0);
    s += std::string(buf) + " 2.500000000E+01\n";
    // a line load (LTYPE 0) at the edge's two nodes (GFLAG), component 3
    s += header(6, 7, 2, 2) + "00100100011000000 2\n" + " 1.000000000E+00 2.000000000E+00\n";
    s += header(99, 0, 0, 1);
    return s;
}

std::string result_file(const std::string& rBody, const std::string& rSuffix) {
    const std::string path = mt::temp_path(rSuffix);
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

// A Fortran unformatted record: its length, the bytes, the length again.
void fortran_record(std::string& rOut, const std::string& rBytes, bool Big) {
    std::uint32_t n = static_cast<std::uint32_t>(rBytes.size());
    char len[4];
    for (int k = 0; k < 4; ++k)
        len[Big ? 3 - k : k] = static_cast<char>((n >> (8 * k)) & 0xFF);
    rOut.append(len, 4);
    rOut += rBytes;
    rOut.append(len, 4);
}

std::string word(std::uint32_t V, bool Big) {
    char b[4];
    for (int k = 0; k < 4; ++k)
        b[Big ? 3 - k : k] = static_cast<char>((V >> (8 * k)) & 0xFF);
    return std::string(b, 4);
}

std::string word(float V, bool Big) {
    std::uint32_t u;
    std::memcpy(&u, &V, 4);
    return word(u, Big);
}

}  // namespace

TEST(Patran, LoadsAndBoundaryConditions) {
    const Mesh mesh = meshioplusplus::read_patran(write_file(quad_with_loads()));
    const NDArray& f = mesh.PointData("patran:force:2");
    EXPECT_DOUBLE_EQ(detail::read_double(f, 2 * 6 + 2), 5.0);  // node 3, component 3
    EXPECT_DOUBLE_EQ(detail::read_double(f, 2 * 6 + 3), -1.0);
    EXPECT_TRUE(std::isnan(detail::read_double(f, 0)));
    EXPECT_EQ(detail::read_int(mesh.PointData("patran:force_frame:2"), 2), 9);
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.PointData("patran:displacement:1"), 1), 0.0);
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.PointData("patran:temperature:4"), 1), 300.0);
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.CellData("patran:element_temperature:4", 0), 0),
                     25.0);
    const NDArray& d = mesh.FieldData("patran:distributed_load");
    EXPECT_EQ(detail::read_int(d, 1), 2);   // set
    EXPECT_EQ(detail::read_int(d, 19), 2);  // edge 2
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.FieldData("patran:distributed_load_values"), 1), 2.0);
    // written back as the same packets
    const std::string out = mt::temp_path(".pat");
    meshioplusplus::write_patran(out, mesh);
    const Mesh back = meshioplusplus::read_patran(out);
    EXPECT_EQ(detail::read_int(back.PointData("patran:force_frame:2"), 2), 9);
    EXPECT_EQ(detail::read_int(back.FieldData("patran:distributed_load"), 19), 2);
    EXPECT_NE(read_text(out).find("00100100011000000 2"), std::string::npos);
}

TEST(Patran, TextAndBinaryResultFiles) {
    const std::string pat = write_file(quad_with_loads());
    // text: nodal (NWIDTH 2) and element (NWIDTH 7: two lines)
    const std::string nod = result_file(
        "title\n        4        4   0.000000E+00        0        2\nsub\nsub\n"
        "       1 0.100000E+01-0.200000E+01\n       3 0.300000E+01 0.400000E+01\n",
        ".nod");
    const std::string els = result_file(
        "title\n    7\nsub\nsub\n       7       4\n"
        " 0.100000E+01 0.200000E+01 0.300000E+01 0.400000E+01 0.500000E+01 0.600000E+01\n"
        " 0.700000E+01\n",
        ".els");
    for (const bool big : {false, true}) {
        std::string bin, rec(80 * 4, ' ');
        rec += word(std::uint32_t{4}, big) + word(std::uint32_t{4}, big) + word(0.0F, big) +
               word(std::uint32_t{0}, big) + word(std::uint32_t{1}, big);
        fortran_record(bin, rec, big);
        fortran_record(bin, std::string(80 * 4, ' '), big);
        fortran_record(bin, std::string(80 * 4, ' '), big);
        fortran_record(bin, word(std::uint32_t{2}, big) + word(7.5F, big), big);
        const std::string dis = result_file(bin, ".dis");
        const Mesh mesh = meshioplusplus::read_patran(pat, {{"n", nod}, {"e", els}, {"d", dis}});
        const NDArray& n = mesh.PointData("n");
        EXPECT_DOUBLE_EQ(detail::read_double(n, 1), -2.0);
        EXPECT_DOUBLE_EQ(detail::read_double(n, 2 * 2 + 1), 4.0);
        EXPECT_TRUE(std::isnan(detail::read_double(n, 2)));  // node 2: no record
        EXPECT_DOUBLE_EQ(detail::read_double(mesh.CellData("e", 0), 6), 7.0);
        const NDArray& dd = mesh.PointData("d");  // one column: 1-D
        ASSERT_EQ(dd.Shape().size(), 1u);
        EXPECT_DOUBLE_EQ(detail::read_double(dd, 1), 7.5);
    }
    EXPECT_THROW(meshioplusplus::read_patran(pat, {{"x", result_file("t\n1\n", ".nod")}}),
                 ReadError);
}
