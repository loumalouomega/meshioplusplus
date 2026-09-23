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
 * @file test_code_aster.cpp
 * @brief Code_Aster `.mail` reader/writer: the syntax (comments, commas, wrapped
 *        records, multi-line headers, the 80-column limit), every element keyword,
 *        groups, the error paths, and a HEXA20 through MED and back.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/code_aster.hpp"
#include "meshioplusplus/region.hpp"
#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/formats/med.hpp"
#endif

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::RegionKind;
using meshioplusplus::WriteError;
namespace detail = meshioplusplus::detail;

std::string write_file(const std::string& rBody) {
    const std::string path = mt::temp_path(".mail");
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

std::string read_text(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::string without_comments(const std::string& rText) {
    std::istringstream lines(rText);
    std::string line, out;
    while (std::getline(lines, line))
        if (line.rfind('%', 0) != 0)
            out += line + "\n";
    return out;
}

std::array<double, 3> point(const Mesh& rMesh, std::int64_t i) {
    std::array<double, 3> p{0, 0, 0};
    const std::size_t dim = rMesh.PointDim();
    for (std::size_t d = 0; d < dim; ++d)
        p[d] = detail::read_double(rMesh.Points(), static_cast<std::size_t>(i) * dim + d);
    return p;
}

std::size_t find_region(const Mesh& rMesh, const std::string& rName, RegionKind Kind) {
    const std::size_t r = rMesh.FindRegion(rName, Kind);
    EXPECT_NE(r, Mesh::npos) << rName;
    return r;
}

std::vector<std::int64_t> entries(const Mesh& rMesh, std::size_t r) {
    const auto& reg = rMesh.Region(r);
    return {reg.Entries(), reg.Entries() + reg.NumEntries()};
}

// A unit hexahedron's 20 nodes named in Code_Aster's HEXA20 order: corners
// 1-8, then the bottom-ring mid-edges (1-2, 2-3, 3-4, 4-1), the vertical ones
// (1-5, 2-6, 3-7, 4-8) and the top ring (5-6, 6-7, 7-8, 8-5).
const char* kHexa20 =
    "% a unit HEXA20 in Code_Aster order\n"
    "TITRE\n"
    " a title line, with a comma and COOR_3D in it\n"
    "FINSF\n"
    "COOR_3D NOM=INDEFINI\n"
    "        NUMIN=1 NUMAX=20\n"
    " A1 0 0 0\n A2 1 0 0\n A3 1 1 0\n A4 0 1 0\n"
    " A5 0 0 1\n A6 1 0 1\n A7 1 1 1\n A8 0 1 1\n"
    " B1 0.5 0 0\n B2 1 0.5 0\n B3 0.5 1 0\n B4 0 0.5 0\n"
    " V1 0 0 0.5\n V2 1 0 0.5\n V3 1 1 0.5\n V4 0 1 0.5\n"
    " T1 0.5 0 1\n T2 1 0.5 1\n T3 0.5 1 1\n T4 0 0.5 1\n"
    "FINSF\n"
    "HEXA20\n"
    " CUBE A1 A2 A3 A4 A5 A6 A7 A8\n"
    "      B1,B2,B3,B4, V1 V2 V3 V4 % a comment in a record\n"
    "      T1 T2 T3 T4\n"
    "FINSF\n"
    "QUAD8 NOM = BOTTOM\n"
    " FLOOR A1 A4 A3 A2 B4 B3 B2 B1\n"
    "FINSF\n"
    "GROUP_MA NOM=SOLID\n CUBE\nFINSF\n"
    "GROUP_NO\n BASE A1 A2 A3 A4\nFINSF\n"
    "GROUP_MA\n BASE FLOOR\nFINSF\n"
    "FIN\n";

// meshio++ hexahedron20 edge order: bottom ring, top ring, verticals.
const int kHexEdges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                              {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

}  // namespace

TEST(CodeAster, ReadsAHandWrittenHexa20IntoMeshioOrder) {
    const Mesh mesh = meshioplusplus::read_code_aster(write_file(kHexa20));
    ASSERT_EQ(mesh.NumPoints(), 20u);
    ASSERT_EQ(mesh.PointDim(), 3u);
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    ASSERT_EQ(mesh.Cells(0).Type(), "hexahedron20");
    ASSERT_EQ(mesh.Cells(1).Type(), "quad8");
    const auto& conn = mesh.Cells(0).Conn();
    for (int e = 0; e < 12; ++e) {
        const auto a = point(mesh, detail::read_int(conn, kHexEdges[e][0]));
        const auto b = point(mesh, detail::read_int(conn, kHexEdges[e][1]));
        const auto m = point(mesh, detail::read_int(conn, 8 + e));
        for (int d = 0; d < 3; ++d)
            EXPECT_DOUBLE_EQ(m[d], 0.5 * (a[d] + b[d])) << "edge " << e;
    }
    EXPECT_EQ(entries(mesh, find_region(mesh, "SOLID", RegionKind::Cell)),
              std::vector<std::int64_t>{0});
    EXPECT_EQ(mesh.Region(find_region(mesh, "SOLID", RegionKind::Cell)).mDim, 3);
    // One name, both kinds: two regions.
    EXPECT_EQ(entries(mesh, find_region(mesh, "BASE", RegionKind::Point)),
              (std::vector<std::int64_t>{0, 1, 2, 3}));
    EXPECT_EQ(entries(mesh, find_region(mesh, "BASE", RegionKind::Cell)),
              std::vector<std::int64_t>{1});
    EXPECT_EQ(mesh.Region(find_region(mesh, "BASE", RegionKind::Cell)).mDim, 2);
}

TEST(CodeAster, Hexa20RoundTripsThroughItsOwnWriter) {
    const Mesh mesh = meshioplusplus::read_code_aster(write_file(kHexa20));
    const std::string out = mt::temp_path(".mail");
    meshioplusplus::write_code_aster(out, mesh);
    const Mesh back = meshioplusplus::read_code_aster(out);
    ASSERT_EQ(back.NumCellBlocks(), 2u);
    for (std::size_t b = 0; b < 2; ++b) {
        const auto& c0 = mesh.Cells(b).Conn();
        const auto& c1 = back.Cells(b).Conn();
        ASSERT_EQ(c0.Size(), c1.Size());
        for (std::size_t i = 0; i < c0.Size(); ++i)
            EXPECT_EQ(detail::read_int(c0, i), detail::read_int(c1, i));
    }
    EXPECT_EQ(back.NumRegions(), mesh.NumRegions());
    // A second write is byte-identical: names, order and layout are canonical.
    const std::string again = mt::temp_path(".mail");
    meshioplusplus::write_code_aster(again, back);
    EXPECT_EQ(read_text(out), read_text(again));
}

TEST(CodeAster, EveryWrittenLineFitsIn80Columns) {
    Mesh mesh = meshioplusplus::read_code_aster(write_file(kHexa20));
    const std::string out = mt::temp_path(".mail");
    meshioplusplus::write_code_aster(out, mesh);
    std::istringstream lines(read_text(out));
    std::string line;
    std::size_t count = 0;
    while (std::getline(lines, line)) {
        EXPECT_LE(line.size(), 80u) << line;
        ++count;
    }
    EXPECT_GT(count, 20u);
    EXPECT_NE(read_text(out).find("\nFIN\n"), std::string::npos);
}

TEST(CodeAster, ReadsEveryElementKeyword) {
    struct Case {
        const char* mKeyword;
        const char* mType;
        int mNodes;
    };
    const Case cases[] = {
        {"POI1", "vertex", 1},          {"SEG2", "line", 2},
        {"SEG3", "line3", 3},           {"SEG4", "line4", 4},
        {"TRIA3", "triangle", 3},       {"TRIA6", "triangle6", 6},
        {"TRIA7", "triangle7", 7},      {"QUAD4", "quad", 4},
        {"QUAD8", "quad8", 8},          {"QUAD9", "quad9", 9},
        {"TETRA4", "tetra", 4},         {"TETRA10", "tetra10", 10},
        {"PENTA6", "wedge", 6},         {"PENTA15", "wedge15", 15},
        {"PENTA18", "wedge18", 18},     {"PYRAM5", "pyramid", 5},
        {"PYRAM13", "pyramid13", 13},   {"HEXA8", "hexahedron", 8},
        {"HEXA20", "hexahedron20", 20}, {"HEXA27", "hexahedron27", 27},
    };
    std::string body = "COOR_3D\n";
    for (int p = 1; p <= 27; ++p)
        body += " N" + std::to_string(p) + " " + std::to_string(p) + " 0 0\n";
    body += "FINSF\n";
    for (const Case& c : cases) {
        body += std::string(c.mKeyword) + "\n E" + c.mKeyword;
        for (int k = 1; k <= c.mNodes; ++k)
            body += " N" + std::to_string(k) + (k % 6 == 0 ? "\n" : "");
        body += "\nFINSF\n";
    }
    body += "FIN\n";
    const Mesh mesh = meshioplusplus::read_code_aster(write_file(body));
    ASSERT_EQ(mesh.NumCellBlocks(), std::size(cases));
    for (std::size_t b = 0; b < std::size(cases); ++b)
        EXPECT_EQ(mesh.Cells(b).Type(), cases[b].mType);
    // And back out: every keyword survives the writer.
    const std::string out = mt::temp_path(".mail");
    meshioplusplus::write_code_aster(out, mesh);
    const Mesh back = meshioplusplus::read_code_aster(out);
    ASSERT_EQ(back.NumCellBlocks(), std::size(cases));
    for (std::size_t b = 0; b < std::size(cases); ++b) {
        const auto& c0 = mesh.Cells(b).Conn();
        const auto& c1 = back.Cells(b).Conn();
        for (std::size_t i = 0; i < c0.Size(); ++i)
            EXPECT_EQ(detail::read_int(c0, i), detail::read_int(c1, i)) << cases[b].mKeyword;
    }
}

TEST(CodeAster, Coor2dGivesTwoDimensionalPoints) {
    const Mesh mesh = meshioplusplus::read_code_aster(
        write_file("coor_2d\n N1 0. 0.\n N2 1.D0 0\n N3 0 1.0E+0\nFINSF\n"
                   "tria3\n M1 N1 N2 N3\nFINSF\nFIN\n"));
    EXPECT_EQ(mesh.PointDim(), 2u);
    EXPECT_DOUBLE_EQ(point(mesh, 1)[0], 1.0);
    const std::string out = mt::temp_path(".mail");
    meshioplusplus::write_code_aster(out, mesh);
    EXPECT_NE(read_text(out).find("COOR_2D\n"), std::string::npos);
}

TEST(CodeAster, CoordinatesKeepFullPrecision) {
    const Mesh mesh = meshioplusplus::read_code_aster(
        write_file("COOR_3D\n N1 0.1 -2.2250738585072014E-300 3.141592653589793\n"
                   "FINSF\nPOI1\n M1 N1\nFINSF\nFIN\n"));
    const std::string out = mt::temp_path(".mail");
    meshioplusplus::write_code_aster(out, mesh);
    const Mesh back = meshioplusplus::read_code_aster(out);
    for (int d = 0; d < 3; ++d)
        EXPECT_EQ(point(back, 0)[d], point(mesh, 0)[d]);
    // A three-digit exponent pushes the record past 80 columns: it wraps.
    std::istringstream lines(read_text(out));
    std::string line;
    while (std::getline(lines, line))
        EXPECT_LE(line.size(), 80u);
}

TEST(CodeAster, ColumnsPast80AreIgnoredAsCodeAsterDoes) {
    const std::string pad(80, ' ');
    const Mesh mesh = meshioplusplus::read_code_aster(
        write_file("COOR_3D\n N1 0 0 0\n N2 1 0 0\n N3 0 1 0\nFINSF\n"
                   "SEG2\n M1 N1 N2\n" +
                   pad + "  M2 N2 N3\nFINSF\nFIN\n"));
    EXPECT_EQ(mesh.Cells(0).NumCells(), 1u);
}

TEST(CodeAster, WriterSanitisesGroupNames) {
    Mesh mesh = mt::tri_mesh();
    meshioplusplus::NDArray e(meshioplusplus::DType::Int64, {1});
    e.As<std::int64_t>()[0] = 0;
    mesh.AddRegion(
        meshioplusplus::Region("a very long group name with spaces", RegionKind::Cell, 2, -1, e));
    mesh.AddRegion(meshioplusplus::Region("a very long group name with spaces too",
                                          RegionKind::Cell, 2, -1, e));
    const std::string out = mt::temp_path(".mail");
    meshioplusplus::write_code_aster(out, mesh);
    const Mesh back = meshioplusplus::read_code_aster(out);
    std::vector<std::string> names;
    for (std::size_t r = 0; r < back.NumRegions(); ++r) {
        names.push_back(back.Region(r).mName);
        EXPECT_LE(names.back().size(), 24u);
        EXPECT_EQ(names.back().find(' '), std::string::npos);
    }
    ASSERT_EQ(names.size(), 2u);
    EXPECT_NE(names[0], names[1]);
}

TEST(CodeAster, ErrorsNameTheCulprit) {
    auto read = [](const std::string& rBody) {
        return meshioplusplus::read_code_aster(write_file(rBody));
    };
    EXPECT_THROW(read("COOR_3D\n N1 0 0 0\nFINSF\nSEG2\n M1 N1 N9\nFINSF\nFIN\n"), ReadError);
    EXPECT_THROW(read("COOR_3D\n N1 0 0 0\n N1 1 0 0\nFINSF\nFIN\n"), ReadError);
    EXPECT_THROW(read("COOR_3D\n N1 0 0 0\nFINSF\nGROUP_NO\n G N2\nFINSF\nFIN\n"), ReadError);
    EXPECT_THROW(read("COOR_3D\n N1 0 0 zero\nFINSF\nFIN\n"), ReadError);
    EXPECT_THROW(read("COOR_3D\n N1 0 0 0\n"), ReadError);
    EXPECT_THROW(read("COOR_3D\n N1 0 0 0\nFINSF\nSEG2\n M1 N1\nFINSF\nFIN\n"), ReadError);
    try {
        read("COOR_3D\n N1 0 0 0\nFINSF\nSEG2\n M1 N1 N9\nFINSF\nFIN\n");
    } catch (const ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("N9"), std::string::npos) << e.what();
    }
}

TEST(CodeAster, UnsupportedKeywordsAreSkipped) {
    const Mesh mesh = meshioplusplus::read_code_aster(
        write_file("COOR_3D\n N1 0 0 0\n N2 1 0 0\nFINSF\n"
                   "SEG22\n M9 N1 N2 N1 N2\nFINSF\nSYS_COOR\n anything\nFINSF\n"
                   "SEG2\n M1 N1 N2\nFINSF\nFIN\n"));
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "line");
}

TEST(CodeAster, WriterRefusesTypesWithNoKeyword) {
    Mesh mesh = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 0, 0}}, "polygon",
                              {{0, 1, 3, 2, 4}});
    EXPECT_THROW(meshioplusplus::write_code_aster(mt::temp_path(".mail"), mesh), WriteError);
}

#ifdef MESHIOPLUSPLUS_HAS_HDF5
TEST(CodeAster, Hexa20RoundTripsThroughMed) {
    const Mesh mesh = meshioplusplus::read_code_aster(write_file(kHexa20));
    const std::string first = mt::temp_path(".mail");
    meshioplusplus::write_code_aster(first, mesh);
    const std::string med = mt::temp_path(".med");
    meshioplusplus::MedInfo info;
    meshioplusplus::write_med(med, mesh, info);
    meshioplusplus::MedInfo back_info;
    const Mesh via_med = meshioplusplus::read_med(med, back_info);
    const std::string second = mt::temp_path(".mail");
    meshioplusplus::write_code_aster(second, via_med);
    // The MED reader adds its family arrays, which the .mail writer drops with a
    // provenance note: compare everything but the `%` comment lines.
    EXPECT_EQ(without_comments(read_text(first)), without_comments(read_text(second)));
}
#endif
