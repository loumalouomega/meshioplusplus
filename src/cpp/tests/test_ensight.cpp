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
#include <filesystem>
#include <fstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ensight.hpp"

namespace {

// EnSight writes a .case/.geo sibling pair; wrap mt::roundtrip so the
// sibling temp file is removed too.
void rt(const mt::Mesh& mesh, bool binary) {
    const double atol = binary ? 1e-6 : 1e-5;  // float32 vs %12.5e ASCII
    std::string sibling;
    mt::roundtrip(
        [&](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_ensight(p, m, binary);
            sibling = p.substr(0, p.size() - 5) + ".geo";
        },
        [](const std::string& p) { return meshioplusplus::read_ensight(p); }, mesh, ".case", atol);
    std::error_code ec;
    std::filesystem::remove(sibling, ec);
}

mt::Mesh wedge15_mesh() {
    return mt::make_mesh({{0, 0, 0},
                          {1, 0, 0},
                          {1, 1, 0},
                          {0, 0, 1},
                          {1, 0, 1},
                          {1, 1, 1},
                          {0.5, 0, 0},
                          {1, 0.5, 0},
                          {0.5, 0.5, 0},
                          {0.5, 0, 1},
                          {1, 0.5, 1},
                          {0.5, 0.5, 1},
                          {0, 0, 0.5},
                          {1, 0, 0.5},
                          {1, 1, 0.5}},
                         "wedge15", {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}});
}

}  // namespace

TEST(Ensight, AsciiLinear) {
    rt(mt::tri_mesh(), false);
    rt(mt::quad_mesh(), false);
    rt(mt::tet_mesh(), false);
    rt(mt::hex_mesh(), false);
    rt(mt::wedge_mesh(), false);
}
TEST(Ensight, AsciiQuadratic) {
    rt(mt::tet10_mesh(), false);
    rt(mt::hex20_mesh(), false);
    rt(wedge15_mesh(), false);
}
TEST(Ensight, AsciiHybrid) {
    rt(mt::tri_quad_mesh(), false);
}
TEST(Ensight, Ascii2D) {
    rt(mt::tri_mesh_2d(), false);
}
TEST(Ensight, BinaryLinear) {
    rt(mt::tri_mesh(), true);
    rt(mt::tet_mesh(), true);
    rt(mt::hex_mesh(), true);
    rt(mt::wedge_mesh(), true);
}
TEST(Ensight, BinaryQuadratic) {
    rt(mt::tet10_mesh(), true);
    rt(mt::hex20_mesh(), true);
    rt(wedge15_mesh(), true);
}

TEST(Ensight, AsciiKeywordsPresent) {
    std::string path = mt::temp_path(".case");
    meshioplusplus::write_ensight(path, mt::tet_mesh(), /*binary=*/false);
    const std::string geo_path = path.substr(0, path.size() - 5) + ".geo";
    std::ifstream in(geo_path, std::ios::binary);
    ASSERT_TRUE(in.good());
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("EnSight Gold Geometry File"), std::string::npos);
    EXPECT_NE(content.find("node id assign"), std::string::npos);
    EXPECT_NE(content.find("element id assign"), std::string::npos);
    EXPECT_NE(content.find("coordinates"), std::string::npos);
    EXPECT_NE(content.find("tetra4"), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(geo_path, ec);
}

TEST(Ensight, WriteRejectsUnknownType) {
    // No EnSight keyword for high-order Lagrange types -> WriteError.
    auto mesh =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}}, "line4", {{0, 1, 2, 3}});
    EXPECT_THROW(meshioplusplus::write_ensight(mt::temp_path(".case"), mesh, false),
                 meshioplusplus::WriteError);
}
