// Gmsh round-trip tests (formats 2.2 and 4.1, ascii + binary).

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshio/formats/gmsh.hpp"

namespace {
void rt22(const mt::Mesh& mesh, bool binary) {
    mt::roundtrip(
        [=](const std::string& p, const mt::Mesh& m) { meshio::write_gmsh22(p, m, binary); },
        [](const std::string& p) { return meshio::read_gmsh(p); }, mesh, ".msh");
}
void rt41(const mt::Mesh& mesh, bool binary) {
    mt::roundtrip(
        [=](const std::string& p, const mt::Mesh& m) { meshio::write_gmsh41(p, m, binary); },
        [](const std::string& p) { return meshio::read_gmsh(p); }, mesh, ".msh");
}
}  // namespace

TEST(Gmsh, V22Ascii) {
    rt22(mt::tri_mesh(), false);
    rt22(mt::tet_mesh(), false);
    rt22(mt::hex_mesh(), false);
}
TEST(Gmsh, V22Binary) {
    rt22(mt::tri_mesh(), true);
    rt22(mt::tet_mesh(), true);
}
TEST(Gmsh, V41Ascii) {
    rt41(mt::tri_mesh(), false);
    rt41(mt::tet_mesh(), false);
}
TEST(Gmsh, V41Binary) { rt41(mt::hex_mesh(), true); }
TEST(Gmsh, SecondOrder) {
    rt22(mt::tet10_mesh(), false);
    rt22(mt::hex20_mesh(), false);
}
