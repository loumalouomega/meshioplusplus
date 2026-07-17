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

// System includes
#include <filesystem>
#include <fstream>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ply.hpp"

namespace {

// PLY has no dedicated cpp/tests file elsewhere -- cpp/src/formats/ply.cpp is
// otherwise reached only through the Python shim. Exercise both the ASCII and
// binary reader/writer halves here.

TEST(Ply, AsciiRoundtrip) {
    auto w = [](const std::string& p, const mt::Mesh& m) {
        meshioplusplus::write_ply(p, m, /*binary=*/false);
    };
    auto r = [](const std::string& p) { return meshioplusplus::read_ply(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".ply");
    mt::roundtrip(w, r, mt::quad_mesh(), ".ply");
    mt::roundtrip(w, r, mt::tri_quad_mesh(), ".ply");
}

TEST(Ply, BinaryRoundtrip) {
    auto w = [](const std::string& p, const mt::Mesh& m) {
        meshioplusplus::write_ply(p, m, /*binary=*/true);
    };
    auto r = [](const std::string& p) { return meshioplusplus::read_ply(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".ply");
    mt::roundtrip(w, r, mt::quad_mesh(), ".ply");
    mt::roundtrip(w, r, mt::tri_quad_mesh(), ".ply");
}

TEST(Ply, ReadRejectsNonPly) {
    // A file that does not start with the "ply" magic must be rejected rather
    // than silently mis-parsed.
    std::string path = mt::temp_path(".ply");
    {
        std::ofstream f(path);
        f << "this is not a ply file\n";
    }
    EXPECT_THROW(meshioplusplus::read_ply(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

}  // namespace
