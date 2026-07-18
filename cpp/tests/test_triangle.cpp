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
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/triangle.hpp"

namespace {

// Triangle writes sibling files; wrap mt::roundtrip so the siblings are
// removed too.
void rt(const mt::Mesh& mesh, const std::string& rSuffix) {
    std::string stem;
    mt::roundtrip(
        [&](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_triangle(p, m);
            stem = p.substr(0, p.size() - rSuffix.size());
        },
        [](const std::string& p) { return meshioplusplus::read_triangle(p); }, mesh, rSuffix);
    std::error_code ec;
    std::filesystem::remove(stem + ".node", ec);
    std::filesystem::remove(stem + ".ele", ec);
    std::filesystem::remove(stem + ".poly", ec);
}

mt::Mesh triangle6_2d_mesh() {
    return mt::make_mesh({{0, 0}, {1, 0}, {1, 1}, {0.5, 0}, {1, 0.5}, {0.5, 0.5}}, "triangle6",
                         {{0, 1, 2, 3, 4, 5}});
}

mt::Mesh line_mesh_2d() {
    return mt::make_mesh({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, "line",
                         {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
}

}  // namespace

TEST(Triangle, NodeElePair) {
    rt(mt::tri_mesh_2d(), ".node");
    rt(mt::tri_mesh_2d(), ".ele");
}
TEST(Triangle, Triangle6) {
    rt(triangle6_2d_mesh(), ".node");
}
TEST(Triangle, Poly) {
    rt(line_mesh_2d(), ".poly");
}
TEST(Triangle, Rejects3DPoints) {
    EXPECT_THROW(meshioplusplus::write_triangle(mt::temp_path(".node"), mt::tri_mesh()),
                 meshioplusplus::WriteError);
}
