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

// System includes
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/lsdyna.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;
using meshioplusplus::WriteError;

std::string pad(const std::string& rText, std::size_t width) {
    return rText.size() >= width ? rText : std::string(width - rText.size(), ' ') + rText;
}

std::string ipad(std::int64_t value, std::size_t width) {
    return pad(std::to_string(value), width);
}

// One card of right-justified integer fields of `width` columns.
std::string ints(std::initializer_list<std::int64_t> values, std::size_t width = 8) {
    std::string out;
    for (std::int64_t v : values)
        out += ipad(v, width);
    return out + "\n";
}

std::string std_node(std::int64_t id, const char* x, const char* y, const char* z) {
    return ipad(id, 8) + pad(x, 16) + pad(y, 16) + pad(z, 16) + "\n";
}

const char* kCubeNodes[9][3] = {
    {"0.0", "0.0", "0.0"}, {"1.0", "0.0", "0.0"}, {"1.0", "1.0", "0.0"},
    {"0.0", "1.0", "0.0"}, {"0.0", "0.0", "1.0"}, {"1.0", "0.0", "1.0"},
    {"1.0", "1.0", "1.0"}, {"0.0", "1.0", "1.0"}, {"0.5", "0.5", "2.0"}};

std::string cube_nodes() {
    std::string out = "*NODE\n";
    for (int i = 0; i < 9; ++i)
        out += std_node(i + 1, kCubeNodes[i][0], kCubeNodes[i][1], kCubeNodes[i][2]);
    return out;
}

std::string write_deck(const std::string& rSuffix, const std::string& rBody) {
    const std::string path = mt::temp_path(rSuffix);
    std::ofstream out(path, std::ios::binary);
    out << rBody;
    return path;
}

std::vector<std::int64_t> row_of(const Mesh& rMesh, std::size_t block, std::size_t row) {
    const auto cb = rMesh.Cells(block);
    const auto& conn = cb.Conn();
    const std::size_t k = cb.NodesPerCell();
    std::vector<std::int64_t> out;
    for (std::size_t c = 0; c < k; ++c)
        out.push_back(meshioplusplus::detail::read_int(conn, row * k + c));
    return out;
}

std::vector<std::int64_t> entries_of(const Mesh& rMesh, const std::string& rName, RegionKind kind) {
    const std::size_t i = rMesh.FindRegion(rName, kind);
    EXPECT_NE(i, Mesh::npos) << rName;
    std::vector<std::int64_t> out;
    if (i == Mesh::npos)
        return out;
    const Region& r = rMesh.Region(i);
    for (std::size_t k = 0; k < r.NumEntries() * r.Stride(); ++k)
        out.push_back(r.Entries()[k]);
    return out;
}

}  // namespace

using V = std::vector<std::int64_t>;

TEST(Lsdyna, RoundTripsTheStandardMeshes) {
    using namespace meshioplusplus;
    const auto write = [](const std::string& p, const Mesh& m) { write_lsdyna(p, m); };
    const auto read = [](const std::string& p) { return read_lsdyna(p); };
    mt::roundtrip(write, read, mt::line_mesh(), ".k", 1e-8);
    mt::roundtrip(write, read, mt::tri_mesh(), ".k", 1e-8);
    mt::roundtrip(write, read, mt::quad_mesh(), ".k", 1e-8);
    mt::roundtrip(write, read, mt::tri_quad_mesh(), ".k", 1e-8);
    mt::roundtrip(write, read, mt::tet_mesh(), ".k", 1e-8);
    mt::roundtrip(write, read, mt::tet10_mesh(), ".k", 1e-8);
    mt::roundtrip(write, read, mt::hex_mesh(), ".k", 1e-8);
    mt::roundtrip(write, read, mt::wedge_mesh(), ".k", 1e-8);
}

TEST(Lsdyna, CollapsesDegenerateHexahedra) {
    const std::string path = write_deck(
        ".k", cube_nodes() + "*ELEMENT_SOLID\n" + ints({1, 1, 1, 2, 3, 4, 5, 6, 7, 8}) +
                  ints({2, 1, 1, 2, 3, 4, 9, 9, 9, 9}) + ints({3, 1, 1, 2, 3, 3, 5, 6, 7, 7}) +
                  ints({4, 1, 1, 2, 3, 9, 9, 9, 9, 9}) + ints({5, 1, 1, 2, 3, 3, 9, 9, 9, 9}) +
                  ints({6, 1, 1, 5, 2, 4, 9, 6, 9, 3}) + "*ELEMENT_SHELL\n" +
                  ints({1, 1, 1, 2, 3, 3}) + ints({2, 1, 1, 2, 3}) + ints({3, 1, 1, 2, 3, 4}));
    const Mesh mesh = meshioplusplus::read_lsdyna(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 6u);
    EXPECT_EQ(std::string(mesh.Cells(0).Type()), "hexahedron");
    EXPECT_EQ(row_of(mesh, 0, 0), (V{0, 1, 2, 3, 4, 5, 6, 7}));
    EXPECT_EQ(std::string(mesh.Cells(1).Type()), "pyramid");
    EXPECT_EQ(row_of(mesh, 1, 0), (V{0, 1, 2, 3, 8}));
    EXPECT_EQ(std::string(mesh.Cells(2).Type()), "wedge");
    EXPECT_EQ(row_of(mesh, 2, 0), (V{0, 1, 2, 4, 5, 6}));
    EXPECT_EQ(std::string(mesh.Cells(3).Type()), "tetra");
    EXPECT_EQ(row_of(mesh, 3, 0), (V{0, 1, 2, 8}));
    EXPECT_EQ(std::string(mesh.Cells(3).Type()), "tetra");
    EXPECT_EQ(row_of(mesh, 3, 1), (V{0, 1, 2, 8}));  // n3 == n4 and n5..n8 equal
    EXPECT_EQ(std::string(mesh.Cells(4).Type()), "triangle");
    EXPECT_EQ(mesh.Cells(4).NumCells(), 2u);  // n4 == n3, and n4 left blank
    EXPECT_EQ(std::string(mesh.Cells(5).Type()), "quad");
}

TEST(Lsdyna, TheOtherWedgeFormIsReadWithPositiveOrientation) {
    // n1 n2 n3 n4 n5 n5 n6 n6: triangles (n1 n2 n5) and (n4 n3 n6).
    const std::string body = "*NODE\n" + std_node(1, "0.0", "0.0", "0.0") +
                             std_node(2, "1.0", "0.0", "0.0") + std_node(3, "1.0", "1.0", "0.0") +
                             std_node(4, "0.0", "1.0", "0.0") + std_node(5, "0.5", "0.0", "1.0") +
                             std_node(6, "0.5", "1.0", "1.0") + "*ELEMENT_SOLID\n" +
                             ints({1, 1, 1, 2, 3, 4, 5, 5, 6, 6});
    const Mesh mesh = meshioplusplus::read_lsdyna(write_deck(".k", body));
    ASSERT_EQ(std::string(mesh.Cells(0).Type()), "wedge");
    const V w = row_of(mesh, 0, 0);
    const auto& pts = mesh.Points();
    auto at = [&](std::int64_t n, int c) {
        return meshioplusplus::detail::read_double(pts, static_cast<std::size_t>(n) * 3 + c);
    };
    double a[3][3], top[3] = {0, 0, 0}, bot[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int c = 0; c < 3; ++c) {
            a[i][c] = at(w[i], c);
            bot[c] += at(w[i], c) / 3.0;
            top[c] += at(w[3 + i], c) / 3.0;
        }
    const double u[3] = {a[1][0] - a[0][0], a[1][1] - a[0][1], a[1][2] - a[0][2]};
    const double v[3] = {a[2][0] - a[0][0], a[2][1] - a[0][1], a[2][2] - a[0][2]};
    const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                         u[0] * v[1] - u[1] * v[0]};
    const double volume =
        (n[0] * (top[0] - bot[0]) + n[1] * (top[1] - bot[1]) + n[2] * (top[2] - bot[2])) / 2.0;
    EXPECT_NEAR(volume, 0.5, 1e-12);
}

TEST(Lsdyna, ReadsAllFourCardFormatsInOneDeck) {
    std::ostringstream body;
    body << "*KEYWORD\n*NODE\n";
    for (int i = 0; i < 4; ++i)
        body << std_node(i + 1, kCubeNodes[i][0], kCubeNodes[i][1], kCubeNodes[i][2]);
    body << "*NODE+\n";  // long: every field 20 columns
    for (int i = 4; i < 8; ++i)
        body << ipad(i + 1, 20) << pad(kCubeNodes[i][0], 20) << pad(kCubeNodes[i][1], 20)
             << pad(kCubeNodes[i][2], 20) << "\n";
    body << "*NODE\n9, 0.5, 0.5, 2.0\n";  // free format, per line
    body << "*ELEMENT_SOLID\n" << ints({1, 1, 1, 2, 3, 4, 5, 6, 7, 8});
    body << "*ELEMENT_SOLID%\n" << ints({2, 1, 1, 2, 3, 4, 9, 9, 9, 9}, 10);  // i10
    body << "*ELEMENT_SOLID\n3,2,1,2,3,3,5,6,7,7\n";
    const Mesh mesh = meshioplusplus::read_lsdyna(write_deck(".k", body.str()));
    EXPECT_EQ(mesh.NumPoints(), 9u);
    EXPECT_EQ(meshioplusplus::detail::read_double(mesh.Points(), 8 * 3 + 2), 2.0);
    ASSERT_EQ(mesh.NumCellBlocks(), 3u);
}

TEST(Lsdyna, DeckWideLongAndI10) {
    {
        std::ostringstream body;
        body << "*KEYWORD LONG=Y\n*NODE\n";
        for (int i = 0; i < 9; ++i)
            body << ipad(i + 1, 20) << pad(kCubeNodes[i][0], 20) << pad(kCubeNodes[i][1], 20)
                 << pad(kCubeNodes[i][2], 20) << "\n";
        body << "*ELEMENT_SOLID\n"
             << ints({1, 1, 1, 2, 3, 4, 5, 6, 7, 8}, 20) << "*PART\nwide\n"
             << ints({7, 1, 1}, 20);
        const Mesh mesh = meshioplusplus::read_lsdyna(write_deck(".k", body.str()));
        EXPECT_EQ(mesh.NumPoints(), 9u);
        const std::size_t i = mesh.FindRegion("wide", RegionKind::Cell);
        ASSERT_NE(i, Mesh::npos);
        EXPECT_EQ(mesh.Region(i).mTag, 7);
    }
    {
        std::ostringstream body;
        body << "*KEYWORD I10=Y\n*NODE\n";
        for (int i = 0; i < 9; ++i)
            body << ipad(i + 1, 10) << pad(kCubeNodes[i][0], 16) << pad(kCubeNodes[i][1], 16)
                 << pad(kCubeNodes[i][2], 16) << "\n";
        body << "*ELEMENT_SOLID\n" << ints({1, 1, 1, 2, 3, 4, 5, 6, 7, 8}, 10);
        const Mesh mesh = meshioplusplus::read_lsdyna(write_deck(".k", body.str()));
        EXPECT_EQ(row_of(mesh, 0, 0), (V{0, 1, 2, 3, 4, 5, 6, 7}));
    }
}

TEST(Lsdyna, PartsAreRegionsAndSetsCarryTheirIds) {
    const std::string path = write_deck(
        ".k", cube_nodes() + "*ELEMENT_SOLID\n" + ints({1, 1, 1, 2, 3, 4, 5, 6, 7, 8}) +
                  ints({2, 2, 1, 2, 3, 4, 9, 9, 9, 9}) + "*ELEMENT_SHELL\n" +
                  ints({1, 3, 1, 2, 3, 4}) + "*PART\ntop\n" + ints({1, 1, 1}, 10) + "*PART\n\n" +
                  ints({2, 1, 1}, 10) + "*PART\nskin\n" + ints({3, 1, 1}, 10) +
                  "*SET_NODE_LIST_TITLE\nfixed\n" + ints({10}, 10) + ints({1, 2, 3}, 10) +
                  "*SET_NODE_LIST_GENERATE\n" + ints({20}, 10) + ints({1, 3, 7, 9}, 10) +
                  "*SET_SOLID_LIST\n" + ints({30}, 10) + ints({1, 2}, 10) + "*SET_PART_LIST\n" +
                  ints({40}, 10) + ints({1, 3}, 10) + "*SET_SEGMENT_TITLE\nfaces\n" +
                  ints({50}, 10) + ints({1, 2, 3, 4}, 10) + ints({1, 2, 9, 9}, 10) +
                  "*END\n*NODE\n999,junk\n");
    const Mesh mesh = meshioplusplus::read_lsdyna(path);
    EXPECT_EQ(mesh.NumPoints(), 9u);

    const std::size_t top = mesh.FindRegion("top", RegionKind::Cell);
    ASSERT_NE(top, Mesh::npos);
    EXPECT_EQ(mesh.Region(top).mTag, 1);
    EXPECT_EQ(mesh.Region(top).mDim, 3);
    EXPECT_EQ(entries_of(mesh, "top", RegionKind::Cell), (V{0}));
    EXPECT_EQ(entries_of(mesh, "Part 2", RegionKind::Cell), (V{1}));  // blank title
    EXPECT_EQ(mesh.Region(mesh.FindRegion("skin", RegionKind::Cell)).mDim, 2);

    EXPECT_EQ(entries_of(mesh, "fixed", RegionKind::Point), (V{0, 1, 2}));
    EXPECT_EQ(mesh.Region(mesh.FindRegion("fixed", RegionKind::Point)).mTag, 10);
    EXPECT_EQ(entries_of(mesh, "NODE set 20", RegionKind::Point), (V{0, 1, 2, 6, 7, 8}));
    EXPECT_EQ(entries_of(mesh, "SOLID set 30", RegionKind::Cell), (V{0, 1}));
    EXPECT_EQ(entries_of(mesh, "PART set 40", RegionKind::Cell), (V{0, 2}));
    EXPECT_EQ(mesh.Region(mesh.FindRegion("PART set 40", RegionKind::Cell)).mDim, -1);
    // Segment 1 is the hexahedron's bottom face (facet 4); segment 2 a triangular
    // face of the pyramid (cell 1), its facet 1.
    EXPECT_EQ(entries_of(mesh, "faces", RegionKind::Side), (V{0, 4, 1, 1}));
}

TEST(Lsdyna, SetsMayReferToElementsInAnIncludedFile) {
    const std::string dir = mt::temp_path("_lsd_include");
    std::filesystem::create_directories(dir + "/sub");
    {
        std::ofstream(dir + "/sub/nodes.k") << cube_nodes();
        std::ofstream(dir + "/sub/elems.k")
            << "*ELEMENT_SOLID\n"
            << ints({1, 1, 1, 2, 3, 4, 5, 6, 7, 8}) << "*PART\nblock\n"
            << ints({1, 1, 1}, 10);
        std::ofstream(dir + "/main.k")
            << "*KEYWORD\n*INCLUDE_PATH\nsub\n*INCLUDE\nnodes.k\n*INCLUDE\nsub/elems.k\n"
            << "*SET_SOLID_LIST_TITLE\neverything\n"
            << ints({1}, 10) << ints({1}, 10) << "*END\n";
    }
    const Mesh mesh = meshioplusplus::read_lsdyna(dir + "/main.k");
    EXPECT_EQ(mesh.NumPoints(), 9u);
    EXPECT_EQ(entries_of(mesh, "everything", RegionKind::Cell), (V{0}));
    EXPECT_EQ(entries_of(mesh, "block", RegionKind::Cell), (V{0}));
    std::filesystem::remove_all(dir);
}

TEST(Lsdyna, IncludeRecursionIsBounded) {
    const std::string dir = mt::temp_path("_lsd_loop");
    std::filesystem::create_directories(dir);
    std::ofstream(dir + "/loop.k") << "*INCLUDE\nloop.k\n";
    EXPECT_THROW(meshioplusplus::read_lsdyna(dir + "/loop.k"), ReadError);
    std::filesystem::remove_all(dir);
}

TEST(Lsdyna, FormatVariantsAndExtraCards) {
    const std::string path =
        write_deck(".k", cube_nodes() + "*ELEMENT_SHELL_THICKNESS\n" + ints({1, 1, 1, 2, 3, 4}) +
                             "     1.0     1.0     1.0     1.0\n" + ints({2, 1, 5, 6, 7, 8}) +
                             "     2.0     2.0     2.0     2.0\n*ELEMENT_BEAM\n" +
                             ints({3, 1, 1, 5, 9}) + "     1.0     1.0\n" + ints({4, 1, 2, 6, 9}) +
                             "*ELEMENT_DISCRETE\n" + ints({5, 1, 3, 7}) + "*ELEMENT_MASS\n" +
                             ipad(1, 8) + ipad(9, 8) + pad("1.5", 16) + ipad(1, 8) + "\n");
    const Mesh mesh = meshioplusplus::read_lsdyna(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 3u);
    EXPECT_EQ(std::string(mesh.Cells(0).Type()), "quad");
    EXPECT_EQ(mesh.Cells(0).NumCells(), 2u);
    EXPECT_EQ(std::string(mesh.Cells(1).Type()), "line");
    EXPECT_EQ(mesh.Cells(1).NumCells(), 3u);  // two beams and a discrete element
    EXPECT_EQ(std::string(mesh.Cells(2).Type()), "vertex");
}

TEST(Lsdyna, TwoLineTetra10) {
    std::string body = "*NODE\n";
    for (int i = 1; i <= 10; ++i)
        body += std_node(i, std::to_string(i).c_str(), std::to_string(i * i % 7).c_str(), "0.0");
    body += "*ELEMENT_SOLID\n" + ints({1, 1}) + ints({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    const Mesh mesh = meshioplusplus::read_lsdyna(write_deck(".k", body));
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(mesh.Cells(0).Type()), "tetra10");
    EXPECT_EQ(row_of(mesh, 0, 0), (V{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}));
}

TEST(Lsdyna, FortranRealSpellings) {
    const std::string body =
        "*NODE\n" + ipad(1, 8) + pad("1.5D+01", 16) + pad("2.5-1", 16) + pad("-3.0E0", 16) + "\n";
    const Mesh mesh = meshioplusplus::read_lsdyna(write_deck(".k", body));
    EXPECT_EQ(meshioplusplus::detail::read_double(mesh.Points(), 0), 15.0);
    EXPECT_EQ(meshioplusplus::detail::read_double(mesh.Points(), 1), 0.25);
    EXPECT_EQ(meshioplusplus::detail::read_double(mesh.Points(), 2), -3.0);
}

TEST(Lsdyna, ReadErrors) {
    using meshioplusplus::read_lsdyna;
    EXPECT_THROW(read_lsdyna(mt::temp_path("_missing.k")), ReadError);
    EXPECT_THROW(read_lsdyna(write_deck(
                     ".k", cube_nodes() + "*ELEMENT_SHELL\n" + ints({1, 1, 1, 2, 3, 77}))),
                 ReadError);  // undefined node
    EXPECT_THROW(read_lsdyna(write_deck(".k", cube_nodes() + std_node(3, "9.0", "9.0", "9.0"))),
                 ReadError);  // duplicate node id (the *NODE keyword carries on)
    EXPECT_THROW(
        read_lsdyna(write_deck(".k", cube_nodes() + "*ELEMENT_SHELL\n" + ints({1, 1, 1, 2, 3, 4}) +
                                         ints({1, 1, 1, 2, 3, 4}))),
        ReadError);  // duplicate element id within a family
    EXPECT_NO_THROW(
        read_lsdyna(write_deck(".k", cube_nodes() + "*ELEMENT_SHELL\n" + ints({1, 1, 1, 2, 3, 4}) +
                                         "*ELEMENT_BEAM\n" + ints({1, 1, 1, 2}))));
    EXPECT_THROW(read_lsdyna(write_deck(".k", "*NODE\n" + ipad(1, 8) + pad("abc", 16) + "\n")),
                 ReadError);
}

TEST(Lsdyna, WriterRejectsCellTypesWithoutACard) {
    EXPECT_THROW(meshioplusplus::write_lsdyna(mt::temp_path(".k"), mt::triangle6_mesh()),
                 WriteError);
}

TEST(Lsdyna, RegionsRoundTrip) {
    using namespace meshioplusplus;
    Mesh mesh = mt::hex_mesh();
    const std::size_t ncells = mesh.Cells(0).NumCells();
    auto entries = [](const std::vector<std::int64_t>& v, std::size_t stride) {
        NDArray a =
            NDArray::Uninit(DType::Int64, stride == 2 ? std::vector<std::size_t>{v.size() / 2, 2}
                                                      : std::vector<std::size_t>{v.size()});
        for (std::size_t k = 0; k < v.size(); ++k)
            a.As<std::int64_t>()[k] = v[k];
        return a;
    };
    std::vector<std::int64_t> all(ncells);
    for (std::size_t i = 0; i < ncells; ++i)
        all[i] = static_cast<std::int64_t>(i);
    mesh.AddRegion(Region("body", RegionKind::Cell, 3, 42, entries(all, 1)));
    mesh.AddRegion(Region("clamped", RegionKind::Point, -1, 5, entries({0, 1, 2}, 1)));
    mesh.AddRegion(Region("cells", RegionKind::Cell, -1, 6, entries({0}, 1)));
    mesh.AddRegion(Region("top", RegionKind::Side, -1, 7, entries({0, 5}, 2)));

    const std::string path = mt::temp_path(".k");
    write_lsdyna(path, mesh);
    const Mesh back = read_lsdyna(path);
    EXPECT_EQ(entries_of(back, "body", RegionKind::Cell), all);
    EXPECT_EQ(back.Region(back.FindRegion("body", RegionKind::Cell)).mTag, 42);
    EXPECT_EQ(back.Region(back.FindRegion("body", RegionKind::Cell)).mDim, 3);
    EXPECT_EQ(entries_of(back, "clamped", RegionKind::Point), (V{0, 1, 2}));
    EXPECT_EQ(entries_of(back, "cells", RegionKind::Cell), (V{0}));
    EXPECT_EQ(entries_of(back, "top", RegionKind::Side), (V{0, 5}));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(KeywordCard, WidthsPerMode) {
    using namespace meshioplusplus::detail;
    const CardField i8{'i', 8};
    EXPECT_EQ(card_field_width(i8, CardMode::Standard), 8);
    EXPECT_EQ(card_field_width(i8, CardMode::I10), 10);
    EXPECT_EQ(card_field_width(i8, CardMode::Long), 20);
    EXPECT_EQ(card_field_width(CardField{'r', 16}, CardMode::I10), 16);
    EXPECT_EQ(card_field_width(CardField{'i', 10}, CardMode::I10), 10);
}

TEST(KeywordCard, SplitAndFreeFormat) {
    using namespace meshioplusplus::detail;
    const std::vector<CardField> layout(4, CardField{'i', 8});
    const auto fixed = split_card("       1       2", layout, CardMode::Standard);
    ASSERT_EQ(fixed.size(), 4u);
    EXPECT_EQ(fixed[0], "1");
    EXPECT_EQ(fixed[1], "2");
    EXPECT_EQ(fixed[2], "");
    const auto free_form = split_card("1, 2,,4,5", layout, CardMode::Standard);
    ASSERT_EQ(free_form.size(), 5u);  // not truncated to the layout
    EXPECT_EQ(free_form[2], "");
    EXPECT_EQ(free_form[4], "5");
}

TEST(KeywordCard, ConvertsAndRejects) {
    using namespace meshioplusplus::detail;
    EXPECT_EQ(card_to_int("", ""), 0);
    EXPECT_EQ(card_to_int("-12", ""), -12);
    EXPECT_THROW(card_to_int("1.5", ""), meshioplusplus::ReadError);
    EXPECT_EQ(card_to_real("", ""), 0.0);
    EXPECT_EQ(card_to_real("1.5D+01", ""), 15.0);
    EXPECT_EQ(card_to_real("2.5-1", ""), 0.25);
    EXPECT_EQ(card_to_real(".5+2", ""), 50.0);
    EXPECT_THROW(card_to_real("1.5x", ""), meshioplusplus::ReadError);
}

TEST(KeywordCard, ErrorsNameTheFormat) {
    using namespace meshioplusplus::detail;
    try {
        card_to_int("x", " here");
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_EQ(std::string(e.what()), "LS-DYNA: invalid integer field 'x' here");
    }
    try {
        card_to_real("1.5x", " in a GRID card", "Nastran");
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_EQ(std::string(e.what()), "Nastran: invalid real field '1.5x' in a GRID card");
    }
    EXPECT_EQ(card_to_int("7", "", "Nastran"), 7);
    EXPECT_EQ(card_to_real("-2.45-16", "", "Nastran"), -2.45e-16);
}

TEST(KeywordCard, FormatReal16FitsAndRoundTrips) {
    using namespace meshioplusplus::detail;
    for (double x : {0.0, 1.0, -2.5, 0.1, 1.0 / 3.0, 123456.789, 1e-5, -1e-30, 1e100, 6.02e23}) {
        const std::string s = format_real16(x);
        EXPECT_LE(s.size(), 16u) << s;
        EXPECT_NEAR(card_to_real(s, ""), x, std::abs(x) * 2e-10) << s;
    }
    EXPECT_EQ(format_real16(0.1), "1.0e-01");
    EXPECT_EQ(card_to_real(format_real16(0.1), ""), 0.1);
    EXPECT_EQ(format_real16(0.0), "0.0");
}

TEST(KeywordCard, FormatRealShortIsPythonsRepr) {
    // The spellings Python's repr gives, which the Python twin reproduces.
    const std::pair<double, const char*> cases[] = {
        {0.1, "0.1"},
        {-106.644339063934, "-106.644339063934"},
        {1e16, "1e+16"},
        {1e15, "1000000000000000.0"},
        {1e-5, "1e-05"},
        {1e-4, "0.0001"},
        {2.5, "2.5"},
        {100.0, "100.0"},
        {-0.0001234, "-0.0001234"},
        {1.7976931348623157e308, "1.7976931348623157e+308"},
        {5e-324, "5e-324"},
        {0.0, "0.0"}};
    for (const auto& [value, text] : cases)
        EXPECT_EQ(meshioplusplus::detail::format_real_short(value), text) << text;
    // format_real_fit keeps a width; format_real16 is width 16.
    EXPECT_EQ(meshioplusplus::detail::format_real_fit(-106.644339063934, 20),
              "-1.0664433906393e+02");
    EXPECT_EQ(meshioplusplus::detail::format_real16(0.1),
              meshioplusplus::detail::format_real_fit(0.1, 16));
}
