// VTK legacy round-trip tests (versions 4.2 and 5.1, ascii + big-endian binary).

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshio/formats/vtk.hpp"

namespace {
void rt(const mt::Mesh& mesh, bool binary, bool v51) {
    mt::roundtrip(
        [=](const std::string& p, const mt::Mesh& m) { meshio::write_vtk(p, m, binary, v51); },
        [](const std::string& p) { return meshio::read_vtk(p); }, mesh, ".vtk");
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
TEST(Vtk, V42Binary) { rt(mt::quad_mesh(), true, false); }
TEST(Vtk, Hybrid) { rt(mt::tri_quad_mesh(), false, true); }
