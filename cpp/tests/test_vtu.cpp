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
#include "meshioplusplus/formats/vtu.hpp"

namespace {
void rt(const mt::Mesh& mesh, bool binary, bool zlib) {
    mt::roundtrip([=](const std::string& p,
                      const mt::Mesh& m) { meshioplusplus::write_vtu(p, m, binary, zlib); },
                  [](const std::string& p) { return meshioplusplus::read_vtu(p); }, mesh, ".vtu");
}
}  // namespace

TEST(Vtu, AsciiTri) {
    rt(mt::tri_mesh(), false, false);
}
TEST(Vtu, AsciiTetHexQuad) {
    rt(mt::tet_mesh(), false, false);
    rt(mt::hex_mesh(), false, false);
    rt(mt::quad_mesh(), false, false);
}
TEST(Vtu, AsciiHybrid) {
    rt(mt::tri_quad_mesh(), false, false);
}
TEST(Vtu, BinaryUncompressed) {
    rt(mt::tri_mesh(), true, false);
    rt(mt::tet_mesh(), true, false);
}
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
TEST(Vtu, BinaryZlib) {
    rt(mt::tri_mesh(), true, true);
    rt(mt::hex_mesh(), true, true);
}
#endif
