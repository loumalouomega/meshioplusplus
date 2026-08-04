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

#include <fstream>
#include <iterator>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/operations/stats.hpp"

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

// --- VTK_POLYHEDRON (type 42) ------------------------------------------------

TEST(Vtu, PolyhedronRoundTrip) {
    // Both directions were refused before v9.19.0: the writer rejected any
    // polyhedron block, and the reader threw on merely SEEING the `faces`
    // array name -- before checking whether any cell was actually type 42.
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}));
    m.AddPolyhedronBlock(
        "polyhedron8",
        {{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}}});

    for (bool binary : {false, true}) {
        const std::string p = mt::temp_path(binary ? "_poly_b.vtu" : "_poly_a.vtu");
        meshioplusplus::write_vtu(p, m, binary, /*zlib=*/false);
        const meshioplusplus::Mesh back = meshioplusplus::read_vtu(p);
        ASSERT_EQ(back.NumCellBlocks(), 1u) << "binary=" << binary;
        const auto cb = back.Cells(0);
        EXPECT_TRUE(cb.IsPolyhedron());
        EXPECT_EQ(cb.NumFaces(0), 6u);
        EXPECT_NEAR(meshioplusplus::compute_stats(back).mUnsignedVolume, 1.0, 1e-12);
        std::error_code ec;
        std::filesystem::remove(p, ec);
    }
}

TEST(Vtu, PolyhedraMixWithOtherCellTypes) {
    // The Python reference forbids mixing (_vtu.py's "cannot mix polyhedral
    // cells with other cell types"), but the FORMAT does not: `faceoffsets`
    // carries -1 for a non-polyhedral cell, which IS the mixing mechanism. An
    // OpenFOAM mesh always mixes, so inheriting that restriction would defeat
    // the point.
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1},
                                    {2, 0, 0},
                                    {2, 1, 0},
                                    {2, 1, 1},
                                    {2, 0, 1}}));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}}));
    m.AddPolyhedronBlock("polyhedron8", {{{1, 2, 10, 5},
                                          {8, 11, 10, 9},
                                          {1, 8, 9, 2},
                                          {5, 10, 11, 4 + 4},
                                          {1, 5, 11, 8},
                                          {2, 9, 10, 2}}});

    const std::string p = mt::temp_path("_poly_mixed.vtu");
    meshioplusplus::write_vtu(p, m, /*binary=*/false, /*zlib=*/false);
    // Both arrays must be present, and faceoffsets must carry the -1 sentinel
    // for the hexahedron -- that is the whole mixing contract.
    std::ifstream in(p, std::ios::binary);
    ASSERT_TRUE(in.good());
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("\"faces\""), std::string::npos);
    EXPECT_NE(text.find("\"faceoffsets\""), std::string::npos);
    EXPECT_NE(text.find("-1"), std::string::npos);

    const meshioplusplus::Mesh back = meshioplusplus::read_vtu(p);
    bool saw_hex = false, saw_poly = false;
    for (const auto cb : back.CellRange()) {
        if (cb.IsPolyhedron())
            saw_poly = true;
        else if (std::string(cb.Type()) == "hexahedron")
            saw_hex = true;
    }
    EXPECT_TRUE(saw_hex) << "the hexahedron was lost in a mixed polyhedral file";
    EXPECT_TRUE(saw_poly) << "the polyhedron was lost in a mixed file";
    std::error_code ec;
    std::filesystem::remove(p, ec);
}
