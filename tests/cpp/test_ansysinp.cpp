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

// System includes
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/ansysinp.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::AnsysInfo;
using meshioplusplus::read_ansysinp;
using meshioplusplus::write_ansysinp;

namespace {

// Round-trip a mesh (no sets) through the C++ ansysInp writer/reader.
void roundtrip_plain(const meshioplusplus::Mesh& mesh, const std::string& suffix) {
    std::string path = mt::temp_path(suffix);
    AnsysInfo win, rout;
    write_ansysinp(path, mesh, win);
    meshioplusplus::Mesh out = read_ansysinp(path, rout);
    mt::expect_mesh_eq(mesh, out);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// Four nodes at the unit square's corners, then `Et` (ET/KEYOPT lines), one
// `(19i9)` element row and `Tail` (e.g. a CMBLOCK).
std::string deck(const std::string& rEt, const std::string& rRow, const std::string& rTail = "") {
    return "/PREP7\n" + rEt +
           "\nNBLOCK,6,SOLID\n(3i9,6e21.13e3)\n"
           "        1        0        0 0.0000000000000E+000 0.0000000000000E+000 "
           "0.0000000000000E+000\n"
           "        2        0        0 1.0000000000000E+000 0.0000000000000E+000 "
           "0.0000000000000E+000\n"
           "        3        0        0 1.0000000000000E+000 1.0000000000000E+000 "
           "0.0000000000000E+000\n"
           "        4        0        0 0.0000000000000E+000 1.0000000000000E+000 "
           "0.0000000000000E+000\n"
           "N,R5.3,LOC,       -1,\nEBLOCK,19,SOLID\n(19i9)\n" +
           rRow + "\n       -1\n" + rTail + "FINISH\n";
}

// An element row: mat, type, real, secnum, 4 zeros, the node count, 0, id 1, nodes.
std::string row(const std::vector<int>& rNodes) {
    char buf[16];
    std::string out;
    for (int v : {1, 1, 1, 1, 0, 0, 0, 0, static_cast<int>(rNodes.size()), 0, 1}) {
        std::snprintf(buf, sizeof(buf), "%9d", v);
        out += buf;
    }
    for (std::size_t k = 0; k < rNodes.size(); ++k) {
        if (k == 8)
            out += '\n';
        std::snprintf(buf, sizeof(buf), "%9d", rNodes[k]);
        out += buf;
    }
    return out;
}

meshioplusplus::Mesh read_text(const std::string& rText,
                               const meshioplusplus::ReadOptions& rOptions = {}) {
    const std::string path = mt::temp_path(".cdb");
    std::ofstream(path, std::ios::binary) << rText;
    AnsysInfo info;
    meshioplusplus::Mesh mesh = read_ansysinp(path, rOptions, info);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    return mesh;
}

}  // namespace

TEST(AnsysInp, TetraRoundtrip) {
    roundtrip_plain(mt::tet_mesh(), ".inp");
}
TEST(AnsysInp, HexRoundtrip) {
    roundtrip_plain(mt::hex_mesh(), ".inp");
}
TEST(AnsysInp, HybridRoundtrip) {
    roundtrip_plain(mt::tri_quad_mesh(), ".inp");
}

TEST(AnsysInp, SetsRoundtrip) {
    meshioplusplus::Mesh mesh = mt::tri_quad_mesh();  // blocks: triangle, quad, triangle
    std::string path = mt::temp_path(".inp");

    AnsysInfo win;
    win.mPointSets["CORNERS"] = {0, 1, 6};
    // mCellSets: one list per cell block (3 blocks). Select cell 0 of the first
    // triangle block and cell 0 of the quad block.
    win.mCellSets["SOME"] = {{0}, {0}, {}};

    write_ansysinp(path, mesh, win);

    AnsysInfo rout;
    meshioplusplus::Mesh out = read_ansysinp(path, rout);
    mt::expect_mesh_eq(mesh, out);

    ASSERT_TRUE(rout.mPointSets.count("CORNERS"));
    std::vector<std::int64_t> ps = rout.mPointSets["CORNERS"];
    std::sort(ps.begin(), ps.end());
    EXPECT_EQ(ps, (std::vector<std::int64_t>{0, 1, 6}));

    ASSERT_TRUE(rout.mCellSets.count("SOME"));
    // 2 element ids selected in total across the blocks.
    std::size_t total = 0;
    for (const auto& blk : rout.mCellSets["SOME"])
        total += blk.size();
    EXPECT_EQ(total, 2u);

    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(AnsysInp, UnknownTypeThrows) {
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}));
    m.AddCellBlock("polygon", mt::conn_from({{0, 1, 2}}));
    AnsysInfo info;
    std::string path = mt::temp_path(".inp");
    EXPECT_THROW(write_ansysinp(path, m, info), meshioplusplus::WriteError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(AnsysInp, EveryWritableTypeRoundTrips) {
    for (const auto& mesh : {mt::tri_mesh(), mt::quad_mesh(), mt::line_mesh(), mt::wedge_mesh(),
                             mt::triangle6_mesh(), mt::quad8_mesh(), mt::tet10_mesh(),
                             mt::hex20_mesh(), mt::wedge15_mesh(), mt::pyramid13_mesh()})
        roundtrip_plain(mesh, ".cdb");
}

TEST(AnsysInp, DegenerateBricksAndShells) {
    // SOLID185 with K == L and O == P is a wedge; with M == N == O == P a pyramid.
    auto wedge = read_text(deck("ET,1,185", row({1, 2, 3, 3, 1, 2, 3, 3})));
    ASSERT_EQ(wedge.NumCellBlocks(), 1u);
    EXPECT_EQ(wedge.Cells(0).Type(), "wedge");
    auto pyramid = read_text(deck("ET,1,SOLID185", row({1, 2, 3, 4, 1, 1, 1, 1})));
    EXPECT_EQ(pyramid.Cells(0).Type(), "pyramid");
    auto tri = read_text(deck("ET,1,181", row({1, 2, 3, 3})));
    EXPECT_EQ(tri.Cells(0).Type(), "triangle");
    // SURF152's fifth (orientation) node leaves a linear quad.
    auto surf = read_text(deck("ET,1,152", row({1, 2, 3, 4, 1})));
    EXPECT_EQ(surf.Cells(0).Type(), "quad");
}

TEST(AnsysInp, ElementTypeByKeyoptAbbreviation) {
    auto mesh = read_text(deck("ET,1,200\nKEYOP,1,1,6", row({1, 2, 3, 4})));
    EXPECT_EQ(mesh.Cells(0).Type(), "quad");
}

TEST(AnsysInp, MissingMidsidesArePlacedOnTheirEdges) {
    // PLANE183 with its midsides absent (three 0s, and a row cut short by one).
    auto mesh = read_text(deck("ET,1,183", row({1, 2, 3, 4, 0, 0, 0})));
    ASSERT_EQ(mesh.Cells(0).Type(), "quad8");
    EXPECT_EQ(mesh.NumPoints(), 8u);
    const auto& conn = mesh.Cells(0).Conn();
    const auto& pts = mesh.Points();
    const auto x = [&](std::int64_t p, int d) {
        return meshioplusplus::detail::read_double(pts, static_cast<std::size_t>(p) * 3 + d);
    };
    for (std::size_t k = 0; k < 4; ++k) {
        const std::int64_t mid = meshioplusplus::detail::read_int(conn, 4 + k);
        const std::int64_t a = meshioplusplus::detail::read_int(conn, k);
        const std::int64_t b = meshioplusplus::detail::read_int(conn, (k + 1) % 4);
        for (int d = 0; d < 3; ++d)
            EXPECT_DOUBLE_EQ(x(mid, d), 0.5 * (x(a, d) + x(b, d)));
    }
}

TEST(AnsysInp, UnknownElementTypeThrowsUnlessLenient) {
    const std::string text = deck("ET,1,999", row({1, 2, 3, 4}));
    EXPECT_THROW(read_text(text), meshioplusplus::ReadError);
    meshioplusplus::ReadOptions lenient;
    lenient.mLenient = true;
    EXPECT_EQ(read_text(text, lenient).NumCellBlocks(), 0u);
    EXPECT_THROW(read_text(deck("", row({1, 2, 3, 4}))), meshioplusplus::ReadError);
    EXPECT_THROW(read_text(deck("ET,1,181", row({1, 2, 3, 9}))), meshioplusplus::ReadError);
}

TEST(AnsysInp, ComponentsBecomeRegionsAndCellData) {
    // A short CMBLOCK (the count says 5) and a trailing `!` comment.
    auto mesh = read_text(deck("ET,1,181 ! shell", row({1, 2, 3, 4}),
                               "CMBLOCK,TOP  ,NODE,       5  ! nodes\n(8i10)\n"
                               "         3        -4\n"
                               "CMBLOCK,ALL,ELEM,1\n(8i10)\n         1\n"));
    const std::size_t top = mesh.FindRegion("TOP", meshioplusplus::RegionKind::Point);
    ASSERT_NE(top, meshioplusplus::Mesh::npos);
    const auto& reg = mesh.Region(top);
    EXPECT_EQ(std::vector<std::int64_t>(reg.Entries(), reg.Entries() + reg.NumEntries()),
              (std::vector<std::int64_t>{2, 3}));
    EXPECT_TRUE(mesh.HasRegion("ALL", meshioplusplus::RegionKind::Cell));
    for (const char* name : {"ansys:element", "ansys:type", "ansys:mat", "ansys:real"})
        EXPECT_TRUE(mesh.HasCellData(name)) << name;
    EXPECT_EQ(meshioplusplus::detail::read_int(mesh.CellData("ansys:element", 0), 0), 181);
}

TEST(AnsysInp, RegionsRoundTripAsComponents) {
    meshioplusplus::Mesh mesh = mt::hex20_mesh();
    mesh.AddRegion(meshioplusplus::Region("FIXED", meshioplusplus::RegionKind::Point, -1, -1,
                                          mt::conn_from({{0, 1, 2}})));
    mesh.AddRegion(meshioplusplus::Region("BODY", meshioplusplus::RegionKind::Cell, -1, -1,
                                          mt::conn_from({{0}})));
    const std::string path = mt::temp_path(".cdb");
    AnsysInfo none, info;
    write_ansysinp(path, mesh, none);
    const meshioplusplus::Mesh back = read_ansysinp(path, info);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    EXPECT_TRUE(back.HasRegion("FIXED", meshioplusplus::RegionKind::Point));
    EXPECT_TRUE(back.HasRegion("BODY", meshioplusplus::RegionKind::Cell));
    EXPECT_EQ(info.mPointSets["FIXED"], (std::vector<std::int64_t>{0, 1, 2}));
}
