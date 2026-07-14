#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "mesh_fixtures.hpp"
#include "meshio/exceptions.hpp"
#include "meshio/formats/ansysinp.hpp"

using meshio::AnsysInfo;
using meshio::read_ansysinp;
using meshio::write_ansysinp;

namespace {

// Round-trip a mesh (no sets) through the C++ ansysInp writer/reader.
void roundtrip_plain(const meshio::Mesh& mesh, const std::string& suffix) {
    std::string path = mt::temp_path(suffix);
    AnsysInfo win, rout;
    write_ansysinp(path, mesh, win);
    meshio::Mesh out = read_ansysinp(path, rout);
    mt::expect_mesh_eq(mesh, out);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

}  // namespace

TEST(AnsysInp, TetraRoundtrip) { roundtrip_plain(mt::tet_mesh(), ".inp"); }
TEST(AnsysInp, HexRoundtrip) { roundtrip_plain(mt::hex_mesh(), ".inp"); }
TEST(AnsysInp, HybridRoundtrip) { roundtrip_plain(mt::tri_quad_mesh(), ".inp"); }

TEST(AnsysInp, SetsRoundtrip) {
    meshio::Mesh mesh = mt::tri_quad_mesh();  // blocks: triangle, quad, triangle
    std::string path = mt::temp_path(".inp");

    AnsysInfo win;
    win.point_sets["CORNERS"] = {0, 1, 6};
    // cell_sets: one list per cell block (3 blocks). Select cell 0 of the first
    // triangle block and cell 0 of the quad block.
    win.cell_sets["SOME"] = {{0}, {0}, {}};

    write_ansysinp(path, mesh, win);

    AnsysInfo rout;
    meshio::Mesh out = read_ansysinp(path, rout);
    mt::expect_mesh_eq(mesh, out);

    ASSERT_TRUE(rout.point_sets.count("CORNERS"));
    std::vector<std::int64_t> ps = rout.point_sets["CORNERS"];
    std::sort(ps.begin(), ps.end());
    EXPECT_EQ(ps, (std::vector<std::int64_t>{0, 1, 6}));

    ASSERT_TRUE(rout.cell_sets.count("SOME"));
    // 2 element ids selected in total across the blocks.
    std::size_t total = 0;
    for (const auto& blk : rout.cell_sets["SOME"]) total += blk.size();
    EXPECT_EQ(total, 2u);

    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(AnsysInp, UnknownTypeThrows) {
    meshio::Mesh m;
    m.points = mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}});
    m.cells.emplace_back("polygon", mt::conn_from({{0, 1, 2}}));
    AnsysInfo info;
    std::string path = mt::temp_path(".inp");
    EXPECT_THROW(write_ansysinp(path, m, info), meshio::WriteError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
