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
 * @file test_nastran.cpp
 * @brief Nastran/OptiStruct bulk-data reader and writer: the three field
 *        layouts and their continuations, quadratic node orders, HyperMesh
 *        components ($HMMOVE / $HMNAME COMP, the property rule), OptiStruct
 *        SET cards, skipped cards, the error paths, and regions written back
 *        as components. The Python twin is tested in tests/python/test_nastran.py,
 *        which also checks both engines agree on the real decks.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <array>
#include <cstdint>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;
using meshioplusplus::WriteError;
namespace detail = meshioplusplus::detail;

std::string nas_write_file(const std::string& rBody) {
    const std::string path = mt::temp_path(".fem");
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

std::string nas_read_text(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// (kind, name, dim, tag) -> entries
using RegionMap =
    std::map<std::tuple<int, std::string, int, std::int64_t>, std::vector<std::int64_t>>;

RegionMap nas_regions(const Mesh& rMesh) {
    RegionMap out;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        std::vector<std::int64_t> e(reg.Entries(), reg.Entries() + reg.NumEntries());
        out[{static_cast<int>(reg.mKind), reg.mName, reg.mDim, reg.mTag}] = e;
    }
    return out;
}

std::array<double, 3> nas_point(const Mesh& rMesh, std::int64_t i) {
    std::array<double, 3> p{0, 0, 0};
    for (std::size_t d = 0; d < 3; ++d)
        p[d] = detail::read_double(rMesh.Points(), static_cast<std::size_t>(i) * 3 + d);
    return p;
}

std::int64_t nas_node(const Mesh& rMesh, std::size_t Block, std::size_t Slot) {
    return detail::read_int(rMesh.Cells(Block).Conn(), Slot);
}

}  // namespace

TEST(Nastran, ReadsSmallLargeAndFreeFieldGrids) {
    const std::string path = nas_write_file(
        "SOL 101\nCEND\nBEGIN BULK\n"
        "GRID    10              1.0     2.0     3.0\n"
        "GRID*   20                              5.-1            -2.5E+00\n"
        "*       7.0\n"
        "GRID,30,4,0.,1.+1,-3.\n"
        "CTRIA3  1       2       10      20      30\n"
        "ENDDATA\n");
    const Mesh mesh = meshioplusplus::read_nastran(path);
    ASSERT_EQ(mesh.NumPoints(), 3u);
    EXPECT_EQ(nas_point(mesh, 0), (std::array<double, 3>{1, 2, 3}));
    EXPECT_EQ(nas_point(mesh, 1), (std::array<double, 3>{0.5, -2.5, 7}));
    EXPECT_EQ(nas_point(mesh, 2), (std::array<double, 3>{0, 10, -3}));
    // nastran:ref: the CP field of GRID 30, blanks as 0; the PID of the element
    EXPECT_EQ(detail::read_int(mesh.PointData("nastran:ref"), 2), 4);
    EXPECT_EQ(detail::read_int(mesh.PointData("nastran:ref"), 0), 0);
    EXPECT_EQ(detail::read_int(mesh.CellData("nastran:ref", 0), 0), 2);
    EXPECT_EQ(mesh.Cells(0).Type(), "triangle");
}

TEST(Nastran, QuadraticSolidsAndContinuations) {
    // CHEXA 20 in Nastran order: corners, bottom mid-edges, VERTICAL mid-edges,
    // top mid-edges; meshio puts the top ring before the vertical edges.
    const std::vector<std::array<double, 3>> c = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                                  {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5},
                              {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}};
    std::string deck = "BEGIN BULK\n";
    char buf[128];
    std::vector<std::array<double, 3>> pts(c);
    for (const auto& e : edges)
        pts.push_back({(c[e[0]][0] + c[e[1]][0]) / 2, (c[e[0]][1] + c[e[1]][1]) / 2,
                       (c[e[0]][2] + c[e[1]][2]) / 2});
    for (std::size_t k = 0; k < pts.size(); ++k) {
        std::snprintf(buf, sizeof(buf), "GRID,%zu,,%g,%g,%g\n", k + 1, pts[k][0], pts[k][1],
                      pts[k][2]);
        deck += buf;
    }
    deck +=
        "CHEXA   7       1       1       2       3       4       5       6       +A\n"
        "+A      7       8       9       10      11      12      13      14      +B\n"
        "+B      15      16      17      18      19      20\n"
        // an implicit continuation: blank tenth field, blank first field
        "CTETRA  8       1       1       2       4       5       9       14      "
        "        \n"
        "        12      13      15      16\n"
        "ENDDATA\n";
    const Mesh mesh = meshioplusplus::read_nastran(nas_write_file(deck));
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron20");
    EXPECT_EQ(mesh.Cells(1).Type(), "tetra10");
    const int meshio_edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                     {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (int k = 0; k < 12; ++k) {
        const auto m = nas_point(mesh, nas_node(mesh, 0, 8 + k));
        const auto a = nas_point(mesh, nas_node(mesh, 0, meshio_edges[k][0]));
        const auto b = nas_point(mesh, nas_node(mesh, 0, meshio_edges[k][1]));
        for (int d = 0; d < 3; ++d)
            EXPECT_DOUBLE_EQ(m[d], (a[d] + b[d]) / 2) << "edge " << k;
    }
    // CTETRA node ids come straight through (the tet10 order is VTK's)
    EXPECT_EQ(nas_node(mesh, 1, 9), 15);
}

TEST(Nastran, HyperMeshComponentsSetsAndProperties) {
    const std::string path = nas_write_file(
        "BEGIN BULK\n"
        "GRID    1               0.0     0.0     0.0\n"
        "GRID    2               1.0     0.0     0.0\n"
        "GRID    3               0.0     1.0     0.0\n"
        "GRID    4               0.0     0.0     1.0\n"
        "CTRIA3  11      5       1       2       3\n"
        "CTRIA3  12      5       1       2       4\n"
        "CTRIA3  13      5       1       3       4\n"
        "CTETRA  14      6       1       2       3       4\n"
        "CBAR    15      9       1       4       0.0     1.0     0.0\n"
        "$\n$HMMOVE       40\n"
        "$             11THRU          12\n"
        "$             14\n"
        "$$\n"
        "$HMSET        3        1 \"tip\" 18\n"
        "SET     3       GRID    LIST    4\n"
        "SET,8,ELEM,LIST,11,13,\n"
        "+,15,99\n"
        "SET     9       ELEM    RANGE   1       2\n"
        "$HMNAME COMP                  40\"moved\"\n"
        "$HMNAME COMP                  50\"by_pid\"        5 \"PSHELL\" 4\n"
        "$HMNAME COMP                  60\"empty\"\n"
        "ENDDATA\n");
    const Mesh mesh = meshioplusplus::read_nastran(path);
    const int cell = static_cast<int>(RegionKind::Cell);
    const int point = static_cast<int>(RegionKind::Point);
    const RegionMap expected = {
        // $HMMOVE with a THRU range and a single id, named after the elements
        {{cell, "moved", -1, 40}, {0, 1, 3}},
        // no $HMMOVE: the elements whose PID is the component's property,
        // except those moved elsewhere
        {{cell, "by_pid", 2, 50}, {2}},
        {{cell, "empty", -1, 60}, {}},
        {{point, "tip", -1, 3}, {3}},
        // 99 names no element: dropped (with a warning)
        {{cell, "set_8", -1, 8}, {0, 2, 4}},
    };
    EXPECT_EQ(nas_regions(mesh), expected);
}

TEST(Nastran, SkipsUnknownCardsWithoutFailing) {
    const Mesh mesh =
        meshioplusplus::read_nastran(nas_write_file("BEGIN BULK\n"
                                                    "GRID    1               0.0     0.0     0.0\n"
                                                    "DESVAR  1       T1      0.1     0.01    1.0\n"
                                                    "CONTACT 2       FREEZE  1       2\n"
                                                    "CONM2   3       1               2.5\n"
                                                    "RBE2    4       1       123456  1\n"
                                                    "PSHELL  1       1       0.1\n"
                                                    "ENDDATA\n"));
    EXPECT_EQ(mesh.NumPoints(), 1u);
    EXPECT_EQ(mesh.NumCellBlocks(), 0u);
}

TEST(Nastran, ErrorPaths) {
    EXPECT_THROW(meshioplusplus::read_nastran(nas_write_file("GRID    1\n")), ReadError);
    // a CTETRA with neither 4 nor 10 nodes
    EXPECT_THROW(meshioplusplus::read_nastran(nas_write_file(
                     "BEGIN BULK\nCTETRA  1       1       1       2       3\nENDDATA\n")),
                 ReadError);
    // a fixed-count element missing a node
    EXPECT_THROW(meshioplusplus::read_nastran(
                     nas_write_file("BEGIN BULK\nGRID    1               0.0     0.0     0.0\n"
                                    "CQUAD4  1       1       1       1       1\nENDDATA\n")),
                 ReadError);
    // an undefined grid
    EXPECT_THROW(meshioplusplus::read_nastran(nas_write_file(
                     "BEGIN BULK\nCTRIA3  1       1       1       2       3\nENDDATA\n")),
                 ReadError);
    // a garbage field names the card
    try {
        meshioplusplus::read_nastran(
            nas_write_file("BEGIN BULK\nGRID    1               x       0.0     0.0\nENDDATA\n"));
        FAIL();
    } catch (const ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("Nastran: invalid real field 'x' in a GRID card"),
                  std::string::npos);
    }
}

TEST(Nastran, WritesRegionsAsHyperMeshComponents) {
    Mesh mesh = mt::tet_mesh();
    auto entries = [](std::vector<std::int64_t> v) {
        meshioplusplus::NDArray a(meshioplusplus::DType::Int64, {v.size()});
        std::copy(v.begin(), v.end(), a.As<std::int64_t>());
        return a;
    };
    mesh.AddRegion(Region("left", RegionKind::Cell, 3, 7, entries({0})));
    mesh.AddRegion(Region("right", RegionKind::Cell, 3, 9, entries({1})));
    mesh.AddRegion(Region("apex", RegionKind::Point, entries({4})));
    const std::string path = mt::temp_path(".fem");
    meshioplusplus::write_nastran(path, mesh);
    const std::string text = nas_read_text(path);
    EXPECT_EQ(text.find("CTETRA_"), std::string::npos);
    EXPECT_NE(
        text.find("$HMMOVE        7\n$              1\n$HMNAME COMP                   7\"left\"\n"),
        std::string::npos);
    const Mesh back = meshioplusplus::read_nastran(path);
    const int cell = static_cast<int>(RegionKind::Cell);
    const RegionMap expected = {{{cell, "left", 3, 7}, {0}}, {{cell, "right", 3, 9}, {1}}};
    EXPECT_EQ(nas_regions(back), expected);

    Mesh bad = mt::tet_mesh();
    bad.AddRegion(Region("bad", RegionKind::Cell, entries({5})));
    EXPECT_THROW(meshioplusplus::write_nastran(mt::temp_path(".fem"), bad), WriteError);
}

TEST(Nastran, QuadraticRoundTripAndRefs) {
    const Mesh mesh = mt::hex20_mesh();
    const std::string path = mt::temp_path(".bdf");
    meshioplusplus::write_nastran(path, mesh);
    const Mesh back = meshioplusplus::read_nastran(path);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "hexahedron20");
    for (std::size_t j = 0; j < 20; ++j)
        EXPECT_EQ(nas_node(back, 0, j), nas_node(mesh, 0, j));
    // no nastran:ref on the input: every field is blank, so none comes back
    EXPECT_FALSE(back.HasCellData("nastran:ref"));
    EXPECT_FALSE(back.HasPointData("nastran:ref"));
}
