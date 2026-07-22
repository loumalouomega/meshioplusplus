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

// PLY has no dedicated tests/cpp file elsewhere -- src/cpp/src/formats/ply.cpp is
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

TEST(Ply, SkinHexMesh) {
    // Default skin=true: a hexahedron writes its 6 boundary quads with the
    // points compacted (all 8 cube corners are on the boundary).
    std::string path = mt::temp_path("_skin.ply");
    meshioplusplus::write_ply(path, mt::hex_mesh(), /*binary=*/true);
    mt::Mesh out = meshioplusplus::read_ply(path);
    EXPECT_EQ(out.NumPoints(), 8u);
    std::size_t quads = 0;
    for (const auto cb : out.CellRange())
        if (cb.Type() == "quad")
            quads += cb.NumCells();
    EXPECT_EQ(quads, 6u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Ply, SkinCompactsUnreferencedPoints) {
    // A tet mesh's skin references all 5 points; but a mesh with an interior
    // node (pyramid-decomposed cube, apex at the center) must not emit it.
    mt::Mesh volume;
    volume.AssignPoints(mt::points_from({{0, 0, 0},
                                         {1, 0, 0},
                                         {1, 1, 0},
                                         {0, 1, 0},
                                         {0, 0, 1},
                                         {1, 0, 1},
                                         {1, 1, 1},
                                         {0, 1, 1},
                                         {0.5, 0.5, 0.5}}));
    volume.AddCellBlock("pyramid", mt::conn_from({{0, 1, 2, 3, 8},
                                                  {7, 6, 5, 4, 8},
                                                  {0, 4, 5, 1, 8},
                                                  {2, 6, 7, 3, 8},
                                                  {0, 3, 7, 4, 8},
                                                  {1, 5, 6, 2, 8}}));
    std::string path = mt::temp_path("_skin_compact.ply");
    meshioplusplus::write_ply(path, volume, /*binary=*/false);
    mt::Mesh out = meshioplusplus::read_ply(path);
    EXPECT_EQ(out.NumPoints(), 8u);  // the center node is gone
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Ply, SkinFalseLegacyDropsVolumeCells) {
    // skin=false keeps the legacy behavior: volume cells are skipped and the
    // file carries the full vertex table with no faces.
    std::string path = mt::temp_path("_legacy.ply");
    meshioplusplus::write_ply(path, mt::tet_mesh(), /*binary=*/false, /*skin=*/false);
    mt::Mesh out = meshioplusplus::read_ply(path);
    EXPECT_EQ(out.NumPoints(), 5u);
    EXPECT_EQ(out.NumCellBlocks(), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
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
