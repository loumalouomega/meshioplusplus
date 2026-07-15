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
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/formats/freefem.hpp"
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/formats/unv.hpp"

#define SIMPLE_RT(WRITER, READER, MESH, SUFFIX)                                  \
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { WRITER(p, m); }, \
                  [](const std::string& p) { return READER(p); }, MESH, SUFFIX)

TEST(FreeFem, Basic) {
    SIMPLE_RT(meshioplusplus::write_freefem, meshioplusplus::read_freefem, mt::tri_mesh_2d(),
              ".msh");
    SIMPLE_RT(meshioplusplus::write_freefem, meshioplusplus::read_freefem, mt::tri_mesh(), ".msh");
    SIMPLE_RT(meshioplusplus::write_freefem, meshioplusplus::read_freefem, mt::tet_mesh(), ".msh");
}

TEST(Mphtxt, Basic) {
    SIMPLE_RT(meshioplusplus::write_mphtxt, meshioplusplus::read_mphtxt, mt::tri_mesh(), ".mphtxt");
    SIMPLE_RT(meshioplusplus::write_mphtxt, meshioplusplus::read_mphtxt, mt::quad_mesh(),
              ".mphtxt");
    SIMPLE_RT(meshioplusplus::write_mphtxt, meshioplusplus::read_mphtxt, mt::hex_mesh(), ".mphtxt");
    SIMPLE_RT(meshioplusplus::write_mphtxt, meshioplusplus::read_mphtxt, mt::tri_quad_mesh(),
              ".mphtxt");
}

TEST(Flux, Basic) {
    SIMPLE_RT(meshioplusplus::write_flux, meshioplusplus::read_flux, mt::tri_mesh(), ".pf3");
    SIMPLE_RT(meshioplusplus::write_flux, meshioplusplus::read_flux, mt::tet_mesh(), ".pf3");
    SIMPLE_RT(meshioplusplus::write_flux, meshioplusplus::read_flux, mt::hex_mesh(), ".pf3");
}

TEST(Unv, LinearAndParabolic) {
    SIMPLE_RT(meshioplusplus::write_unv, meshioplusplus::read_unv, mt::tri_mesh(), ".unv");
    SIMPLE_RT(meshioplusplus::write_unv, meshioplusplus::read_unv, mt::tet_mesh(), ".unv");
    SIMPLE_RT(meshioplusplus::write_unv, meshioplusplus::read_unv, mt::hex_mesh(), ".unv");
    // parabolic elements exercise the Salome mid-node "sandwich" permutation
    SIMPLE_RT(meshioplusplus::write_unv, meshioplusplus::read_unv, mt::triangle6_mesh(), ".unv");
    SIMPLE_RT(meshioplusplus::write_unv, meshioplusplus::read_unv, mt::quad8_mesh(), ".unv");
    SIMPLE_RT(meshioplusplus::write_unv, meshioplusplus::read_unv, mt::tet10_mesh(), ".unv");
    SIMPLE_RT(meshioplusplus::write_unv, meshioplusplus::read_unv, mt::hex20_mesh(), ".unv");
}
