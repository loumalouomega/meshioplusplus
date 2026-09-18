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

#include <chrono>
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

TEST(Vtu, ThreePolyhedraWithAnInterleavedCellReadCorrectly) {
    // Regression test for the O(ncells^2) `faceoffsets` rescan in
    // reconstruct_cells: with only one polyhedron (PolyhedronRoundTrip) or one
    // polyhedron preceded by one other cell (PolyhedraMixWithOtherCellTypes),
    // the inner rescan never actually runs more than once, so neither test
    // above can distinguish the correct fix (a running "last end offset"
    // carried across the WHOLE loop) from two plausible wrong ones: resetting
    // that running value at the start of each outer-loop run, or at the start
    // of each polyhedral run specifically. Both wrong forms are only exposed
    // when a polyhedral RUN is preceded by another polyhedral run separated by
    // an intervening non-polyhedral one -- exactly this mesh's block order:
    // polyhedron, polyhedron, hexahedron, polyhedron (four distinct blocks,
    // written and therefore read back as three vtk-type runs: [42,42], [12],
    // [42]). If the third polyhedron's face-stream start were computed wrong,
    // its decoded cell would come out non-manifold (wrong node count / a
    // volume of zero or a nonsensical value) rather than throwing, so the
    // volume check below is the one that actually catches it.
    auto unit_cube_faces = [](std::int64_t base) {
        return std::vector<std::vector<std::int64_t>>{
            {base + 0, base + 3, base + 2, base + 1}, {base + 4, base + 5, base + 6, base + 7},
            {base + 0, base + 1, base + 5, base + 4}, {base + 2, base + 3, base + 7, base + 6},
            {base + 0, base + 4, base + 7, base + 3}, {base + 1, base + 2, base + 6, base + 5}};
    };

    std::vector<std::vector<double>> pts;
    for (int block = 0; block < 4; ++block) {
        const double ox = 2.0 * block;  // disjoint, non-overlapping cubes
        for (const auto& c : {std::vector<double>{0, 0, 0},
                              {1, 0, 0},
                              {1, 1, 0},
                              {0, 1, 0},
                              {0, 0, 1},
                              {1, 0, 1},
                              {1, 1, 1},
                              {0, 1, 1}})
            pts.push_back({c[0] + ox, c[1], c[2]});
    }

    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock("polyhedron8", {unit_cube_faces(0)});
    m.AddPolyhedronBlock("polyhedron8", {unit_cube_faces(8)});
    m.AddCellBlock("hexahedron", mt::conn_from({{16, 17, 18, 19, 20, 21, 22, 23}}));
    m.AddPolyhedronBlock("polyhedron8", {unit_cube_faces(24)});

    const std::string p = mt::temp_path("_poly_interleaved.vtu");
    meshioplusplus::write_vtu(p, m, /*binary=*/false, /*zlib=*/false);
    const meshioplusplus::Mesh back = meshioplusplus::read_vtu(p);
    std::error_code ec;
    std::filesystem::remove(p, ec);

    // Two polyhedron<8> blocks (the two runs get bucketed together by unique
    // node count, since bucketing is per-run) plus one hexahedron block.
    std::size_t num_poly_cells = 0, num_hex_cells = 0;
    for (const auto cb : back.CellRange()) {
        if (cb.IsPolyhedron()) {
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                EXPECT_EQ(cb.NumFaces(r), 6u) << "cell " << num_poly_cells;
                ++num_poly_cells;
            }
        } else if (std::string(cb.Type()) == "hexahedron") {
            num_hex_cells += cb.NumCells();
        }
    }
    EXPECT_EQ(num_poly_cells, 3u);
    EXPECT_EQ(num_hex_cells, 1u);

    // Every cell -- all four blocks -- must be a unit cube; a wrong
    // begin-offset for the third polyhedron would slice into an unrelated
    // part of the faces stream (or throw as out-of-range), never quietly
    // reproduce the right geometry.
    const auto stats = meshioplusplus::compute_stats(back);
    EXPECT_NEAR(stats.mUnsignedVolume, 4.0, 1e-12);
}

TEST(Vtu, PolyhedronReconstructionIsNotQuadraticInCellCount) {
    // Direct regression test for the roadmap's own probe: the pre-fix
    // reconstruct_cells rescanned `faceoffsets` from index 0 for EVERY cell,
    // so total work was O(ncells^2) regardless of how much face data there
    // was -- a 100k-cell file cost ~5e9 trivial-loop iterations, on the order
    // of several seconds to "effectively hangs" depending on the machine. The
    // fix makes it O(ncells + |faces|); this asserts the read completes well
    // inside a generous budget that the O(n^2) version could not plausibly
    // meet. All 100k cells reuse the same 8 points/6 faces on purpose --
    // construction cost and memory are irrelevant here, only the read-side
    // cell-reconstruction loop is being timed.
    constexpr std::size_t kNumCells = 100000;
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}));
    std::vector<std::vector<std::vector<std::int64_t>>> cells(
        kNumCells,
        std::vector<std::vector<std::int64_t>>{
            {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}});
    m.AddPolyhedronBlock("polyhedron8", std::move(cells));

    const std::string p = mt::temp_path("_poly_100k.vtu");
    meshioplusplus::write_vtu(p, m, /*binary=*/false, /*zlib=*/false);

    const auto t0 = std::chrono::steady_clock::now();
    const meshioplusplus::Mesh back = meshioplusplus::read_vtu(p);
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0);

    std::error_code ec;
    std::filesystem::remove(p, ec);

    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).NumCells(), kNumCells);
    // A linear pass over 100k cells / 600k face entries is milliseconds; the
    // quadratic version was on the order of seconds to tens of seconds. 3s
    // is generous headroom for a loaded CI machine while still failing fast
    // if the O(n^2) behaviour ever comes back.
    EXPECT_LT(elapsed.count(), 3.0) << "reconstruct_cells took " << elapsed.count() << "s for "
                                    << kNumCells << " cells -- looks quadratic again";
}
