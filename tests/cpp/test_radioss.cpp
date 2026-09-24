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
 * @file test_radioss.cpp
 * @brief OpenRadioss starter decks: degenerate bricks, the BRIC20 node order, an
 *        inverted TETRA4, parts, groups, surfaces, includes, old 8/16-column
 *        fields and the refusals.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/radioss.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::RegionKind;
namespace detail = meshioplusplus::detail;
namespace fs = std::filesystem;

std::string ints(const std::vector<long long>& rValues, int Width = 10) {
    std::string out;
    char buf[32];
    for (long long v : rValues) {
        std::snprintf(buf, sizeof(buf), "%*lld", Width, v);
        out += buf;
    }
    return out + "\n";
}

std::string node(long long Id, double X, double Y, double Z, int Iw = 10, int Rw = 20) {
    char buf[96];
    std::snprintf(buf, sizeof(buf), "%*lld%*.6f%*.6f%*.6f\n", Iw, Id, Rw, X, Rw, Y, Rw, Z);
    return buf;
}

// Nodes 1..8 a unit cube, 9..28 a BRIC20 at x in [2, 3] (Radioss order: bottom
// ring, verticals, top ring), 29..32 an inverted tetra at x = 5.
std::string deck(const fs::path& rDir) {
    std::string s = "#RADIOSS STARTER\n/BEGIN\nrun\n" + ints({2019, 0}) + "\n\n";
    s += "#include part2.inc\n/NODE\n";
    const double c[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    for (int k = 0; k < 8; ++k)
        s += node(k + 1, c[k][0], c[k][1], c[k][2]);
    for (int k = 0; k < 8; ++k)
        s += node(9 + k, c[k][0] + 2, c[k][1], c[k][2]);
    const int ring[3][4][2] = {{{0, 1}, {1, 2}, {2, 3}, {3, 0}},
                               {{0, 4}, {1, 5}, {2, 6}, {3, 7}},
                               {{4, 5}, {5, 6}, {6, 7}, {7, 4}}};
    long long id = 17;
    for (const auto& r : ring)
        for (const auto& e : r)
            s += node(id++, (c[e[0]][0] + c[e[1]][0]) / 2 + 2, (c[e[0]][1] + c[e[1]][1]) / 2,
                      (c[e[0]][2] + c[e[1]][2]) / 2);
    s += node(29, 5, 0, 0) + node(30, 5, 1, 0) + node(31, 6, 0, 0) + node(32, 5, 0, 1);
    s += "/PART/1\nbricks\n" + ints({1, 1, 0});
    s += "/BRICK/1\n" + ints({1, 1, 2, 3, 4, 5, 6, 7, 8}) + ints({2, 1, 2, 3, 1, 5, 6, 7, 5});
    s += "/BRIC20/1\n" + ints({3, 9, 10, 11, 12, 13, 14, 15, 16}) +
         ints({17, 18, 19, 20, 21, 22, 23, 24}) + ints({25, 26, 27, 28});
    s += "/TETRA4/1\n" + ints({4, 29, 30, 31, 32});
    s += "/SHELL/3\n" + ints({1, 5, 6, 7, 8});
    s += "/GRBRIC/BRIC/10\nsome\n" + ints({1, 2, 3, -2});
    s += "/GRNOD/PART/11\npart nodes\n" + ints({1});
    s += "/SURF/SEG/12\nbottom\n" + ints({1, 1, 4, 3, 2});
    s += "/END\n/NODE\n" + node(99, 9, 9, 9);
    std::ofstream(rDir / "part2.inc")
        << "/PART/3\nskin\n" + ints({3, 1, 0}) + "#enddata\n/PART/4\nx\n";
    return s;
}

std::string write_deck(const std::string& rBody, fs::path* pDir = nullptr) {
    // A fresh directory per deck (its #include sits beside it); unique across
    // ctest processes, whose `mt::temp_path` counters all start at 0.
    static std::atomic<unsigned> counter{0};
    const fs::path dir =
        fs::temp_directory_path() /
        ("meshio_rad_" + std::to_string(std::random_device{}()) + "_" + std::to_string(counter++));
    fs::remove_all(dir);
    fs::create_directories(dir);
    std::string body = rBody.empty() ? deck(dir) : rBody;
    const fs::path path = dir / "deck_0000.rad";
    std::ofstream(path, std::ios::binary) << body;
    if (pDir)
        *pDir = dir;
    return path.string();
}

double coord(const Mesh& rMesh, std::size_t Block, std::size_t Slot, std::size_t Dim) {
    const auto n = static_cast<std::size_t>(detail::read_int(rMesh.Cells(Block).Conn(), Slot));
    return detail::read_double(rMesh.Points(), n * 3 + Dim);
}

std::size_t block_of(const Mesh& rMesh, const std::string& rType) {
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        if (rMesh.Cells(b).Type() == rType)
            return b;
    return rMesh.NumCellBlocks();
}

}  // namespace

TEST(Radioss, ReadsEveryCardAndRegion) {
    const Mesh mesh = meshioplusplus::read_radioss(write_deck(""));
    EXPECT_EQ(mesh.NumPoints(), 32u);  // node 99 is after /END
    ASSERT_LT(block_of(mesh, "wedge"), mesh.NumCellBlocks());
    const std::size_t hex20 = block_of(mesh, "hexahedron20");
    ASSERT_LT(hex20, mesh.NumCellBlocks());
    // meshio++ slot 12 is the top edge 4-5 (z = 1), slot 16 the vertical 0-4.
    EXPECT_EQ(coord(mesh, hex20, 12, 2), 1.0);
    EXPECT_EQ(coord(mesh, hex20, 16, 2), 0.5);
    // The inverted TETRA4 comes out positive: node 1 and 2 swapped back.
    const std::size_t tet = block_of(mesh, "tetra");
    ASSERT_LT(tet, mesh.NumCellBlocks());
    EXPECT_EQ(coord(mesh, tet, 1, 0), 6.0);
    const std::size_t shells = mesh.FindRegion("skin", RegionKind::Cell);  // from the include
    ASSERT_NE(shells, Mesh::npos);
    EXPECT_EQ(mesh.Region(shells).NumEntries(), 1u);
    EXPECT_FALSE(mesh.HasRegion("x"));  // after #enddata
    const std::size_t some = mesh.FindRegion("some", RegionKind::Cell);
    ASSERT_NE(some, Mesh::npos);
    EXPECT_EQ(mesh.Region(some).NumEntries(), 2u);  // bricks 1 and 3; 2 removed
    const std::size_t nodes = mesh.FindRegion("part nodes", RegionKind::Point);
    ASSERT_NE(nodes, Mesh::npos);
    EXPECT_EQ(mesh.Region(nodes).NumEntries(), 32u);  // cube, BRIC20 and tetra nodes
    const std::size_t bottom = mesh.FindRegion("bottom", RegionKind::Side);
    ASSERT_NE(bottom, Mesh::npos);
    EXPECT_EQ(mesh.Region(bottom).NumEntries(), 1u);
    EXPECT_EQ(detail::read_int(mesh.FieldData("radioss:version"), 0), 2019);
}

TEST(Radioss, OldFormatColumnsAndRefusals) {
    std::string old = "#RADIOSS STARTER\n/BEGIN\nold\n" + ints({44, 0}, 8) + "\n\n/NODE\n";
    const double c[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    for (int k = 0; k < 8; ++k)
        old += node(k + 1, c[k][0], c[k][1], c[k][2], 8, 16);
    old +=
        "/PART/1\np\n" + ints({1, 1, 0}, 8) + "/BRICK/1\n" + ints({1, 1, 2, 3, 4, 5, 6, 7, 8}, 8);
    const Mesh mesh = meshioplusplus::read_radioss(write_deck(old + "/END\n"));
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron");
    EXPECT_EQ(detail::read_int(mesh.FieldData("radioss:version"), 0), 44);

    EXPECT_THROW(meshioplusplus::read_radioss(write_deck("#RADIOSS ENGINE\n/RUN/r/1\n")),
                 ReadError);
    const std::string missing = "#RADIOSS STARTER\n/BEGIN\nr\n" + ints({2019, 0}) + "\n\n/NODE\n" +
                                node(1, 0, 0, 0) + "/TRUSS/1\n" + ints({1, 1, 2}) + "/END\n";
    EXPECT_THROW(meshioplusplus::read_radioss(write_deck(missing)), ReadError);
}

TEST(Radioss, UnitsBoxesGeneratorsAndSurfaces) {
    // Two bricks (part 1) along x in millimetres, metres of work units.
    std::string body = "#RADIOSS STARTER\n/BEGIN\nfeatures\n" + ints({2022, 0});
    body += "                  kg                  mm                   s\n";
    body += "                  kg                   m                   s\n/NODE\n";
    auto id = [](int i, int j, int k) { return 1 + i + 3 * j + 6 * k; };
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 3; ++i)
                body += node(id(i, j, k), 1000.0 * i, 1000.0 * j, 1000.0 * k);
    body += "/BRICK/1\n";
    for (int i = 0; i < 2; ++i)
        body += ints({101 + i, id(i, 0, 0), id(i + 1, 0, 0), id(i + 1, 1, 0), id(i, 1, 0),
                      id(i, 0, 1), id(i + 1, 0, 1), id(i + 1, 1, 1), id(i, 1, 1)});
    body += "/PART/1\nsolid\n" + ints({1, 5});
    body += "/BOX/RECTA/1\nnear half\n" + ints({0, 0, 0});
    body += node(0, -1.0, -1.0, -1.0, 0).substr(1) + node(0, 1001.0, 1001.0, 1001.0, 0).substr(1);
    body += "/GRNOD/BOX/10\nin box\n" + ints({1});
    body += "/GRBRIC/BOX/11\nbricks inside\n" + ints({1});
    body += "/GRBRIC/BOX2/12\nbricks touching\n" + ints({1});
    body += "/GRNOD/GENE/13\nranges\n" + ints({1, 3, 7, 8});
    body += "/GRNOD/NODE/14/1\nunit suffix\n" + ints({5});
    body += "/SURF/PART/EXT/20\noutside\n" + ints({1});
    body += "/SURF/GRBRIC/FREE/21\nfirst free\n" + ints({11});
    body += "/SURF/SURF/22\nsecond free\n" + ints({20, -21});
    body += "/END\n";
    const Mesh mesh = meshioplusplus::read_radioss(write_deck(body));
    EXPECT_DOUBLE_EQ(mesh.FieldData("radioss:length_scale").As<double>()[0], 1e-3);
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.Points(), 3 * 2), 2.0);  // node 3: x = 2 m
    auto entries = [&](const char* pName, RegionKind Kind) {
        const std::size_t r = mesh.FindRegion(pName, Kind);
        return r == Mesh::npos ? std::size_t(999) : mesh.Region(r).NumEntries();
    };
    EXPECT_EQ(entries("in box", RegionKind::Point), 8u);
    EXPECT_EQ(entries("bricks inside", RegionKind::Cell), 1u);
    EXPECT_EQ(entries("bricks touching", RegionKind::Cell), 2u);
    EXPECT_EQ(entries("ranges", RegionKind::Point), 5u);
    const std::size_t suffix = mesh.FindRegion("unit suffix", RegionKind::Point);
    ASSERT_NE(suffix, Mesh::npos);
    EXPECT_EQ(mesh.Region(suffix).mTag, 14);
    EXPECT_EQ(entries("outside", RegionKind::Side), 10u);
    EXPECT_EQ(entries("first free", RegionKind::Side), 5u);
    EXPECT_EQ(entries("second free", RegionKind::Side), 5u);
}
