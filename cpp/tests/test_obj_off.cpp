// OBJ and OFF round-trip tests (surface formats: triangle/quad/polygon).

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshio/formats/obj_off.hpp"

TEST(Obj, TriQuadHybrid) {
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshio::write_obj(p, m); },
                  [](const std::string& p) { return meshio::read_obj(p); },
                  mt::tri_mesh(), ".obj");
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshio::write_obj(p, m); },
                  [](const std::string& p) { return meshio::read_obj(p); },
                  mt::quad_mesh(), ".obj");
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshio::write_obj(p, m); },
                  [](const std::string& p) { return meshio::read_obj(p); },
                  mt::tri_quad_mesh(), ".obj");
}

TEST(Off, Triangles) {
    // OFF in meshio is triangle-only.
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshio::write_off(p, m); },
                  [](const std::string& p) { return meshio::read_off(p); },
                  mt::tri_mesh(), ".off");
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshio::write_off(p, m); },
                  [](const std::string& p) { return meshio::read_off(p); },
                  mt::tri_mesh_2d(), ".off");
}
