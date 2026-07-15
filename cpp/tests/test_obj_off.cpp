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
#include "meshioplusplus/formats/obj_off.hpp"

TEST(Obj, TriQuadHybrid) {
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_obj(p, m); },
                  [](const std::string& p) { return meshioplusplus::read_obj(p); }, mt::tri_mesh(),
                  ".obj");
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_obj(p, m); },
                  [](const std::string& p) { return meshioplusplus::read_obj(p); }, mt::quad_mesh(),
                  ".obj");
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_obj(p, m); },
                  [](const std::string& p) { return meshioplusplus::read_obj(p); },
                  mt::tri_quad_mesh(), ".obj");
}

TEST(Off, Triangles) {
    // OFF in meshio is triangle-only.
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_off(p, m); },
                  [](const std::string& p) { return meshioplusplus::read_off(p); }, mt::tri_mesh(),
                  ".off");
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_off(p, m); },
                  [](const std::string& p) { return meshioplusplus::read_off(p); },
                  mt::tri_mesh_2d(), ".off");
}
