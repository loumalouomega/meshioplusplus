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
#include <filesystem>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ansysinp.hpp"

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
