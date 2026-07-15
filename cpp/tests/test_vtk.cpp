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
#include "meshioplusplus/formats/vtk.hpp"

namespace {
void rt(const mt::Mesh& mesh, bool binary, bool v51) {
    mt::roundtrip([=](const std::string& p,
                      const mt::Mesh& m) { meshioplusplus::write_vtk(p, m, binary, v51); },
                  [](const std::string& p) { return meshioplusplus::read_vtk(p); }, mesh, ".vtk");
}
}  // namespace

TEST(Vtk, V51Ascii) {
    rt(mt::tri_mesh(), false, true);
    rt(mt::tet_mesh(), false, true);
    rt(mt::hex_mesh(), false, true);
}
TEST(Vtk, V51Binary) {
    rt(mt::tri_mesh(), true, true);
    rt(mt::tet_mesh(), true, true);
}
TEST(Vtk, V42Ascii) {
    rt(mt::tri_mesh(), false, false);
    rt(mt::tet_mesh(), false, false);
}
TEST(Vtk, V42Binary) {
    rt(mt::quad_mesh(), true, false);
}
TEST(Vtk, Hybrid) {
    rt(mt::tri_quad_mesh(), false, true);
}
