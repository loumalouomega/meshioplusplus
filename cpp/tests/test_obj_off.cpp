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
#include <cstdio>
#include <fstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
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
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_off(p, m); },
                  [](const std::string& p) { return meshioplusplus::read_off(p); }, mt::tri_mesh(),
                  ".off");
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_off(p, m); },
                  [](const std::string& p) { return meshioplusplus::read_off(p); },
                  mt::tri_mesh_2d(), ".off");
}

TEST(Off, QuadAndMixedFaces) {
    // https://github.com/loumalouomega/meshioplusplus/issues/35 — OFF files
    // with non-triangular faces used to be rejected outright.
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_off(p, m); },
                  [](const std::string& p) { return meshioplusplus::read_off(p); }, mt::quad_mesh(),
                  ".off");
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_off(p, m); },
                  [](const std::string& p) { return meshioplusplus::read_off(p); },
                  mt::tri_quad_mesh(), ".off");
}

TEST(Off, RejectsDegenerateFaceCount) {
    const std::string path = mt::temp_path(".off");
    {
        std::ofstream out(path);
        out << "OFF\n2 1 0\n0 0 0\n1 0 0\n2 0 1\n";
    }
    EXPECT_THROW(meshioplusplus::read_off(path), meshioplusplus::ReadError);
    std::remove(path.c_str());
}
