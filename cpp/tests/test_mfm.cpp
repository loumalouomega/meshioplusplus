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

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/mfm.hpp"

namespace {

void mfm_roundtrip(const mt::Mesh& mesh) {
    mt::roundtrip(
        [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_mfm(p, m, ".16e"); },
        [](const std::string& p) { return meshioplusplus::read_mfm(p); }, mesh, ".mfm");
}

}  // namespace

TEST(Mfm, Line) {
    mfm_roundtrip(mt::line_mesh());
}
TEST(Mfm, Triangle) {
    mfm_roundtrip(mt::tri_mesh());
}
TEST(Mfm, Triangle2D) {
    mfm_roundtrip(mt::tri_mesh_2d());
}
TEST(Mfm, Quad) {
    mfm_roundtrip(mt::quad_mesh());
}
TEST(Mfm, Tetra) {
    mfm_roundtrip(mt::tet_mesh());
}
TEST(Mfm, Hexahedron) {
    mfm_roundtrip(mt::hex_mesh());
}
TEST(Mfm, Wedge) {
    mfm_roundtrip(mt::wedge_mesh());
}

TEST(Mfm, RejectsMixedTypes) {
    std::string path = mt::temp_path(".mfm");
    EXPECT_THROW(meshioplusplus::write_mfm(path, mt::tri_quad_mesh(), ".16e"),
                 meshioplusplus::WriteError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
