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

#include <algorithm>
#include <chrono>
#include <map>
#include <fstream>
#include <iterator>
#include <limits>

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

namespace {

/// Seconds to read back an ASCII VTU of @p NumCells identical unit-cube
/// polyhedra (best of two reads, which damps scheduler jitter). Every cell
/// reuses the same 8 points / 6 faces on purpose: construction cost and memory
/// are irrelevant here, only the read-side cell-reconstruction loop is timed.
double vtu_polyhedron_read_seconds(std::size_t NumCells) {
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}));
    std::vector<std::vector<std::vector<std::int64_t>>> cells(
        NumCells,
        std::vector<std::vector<std::int64_t>>{
            {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}});
    m.AddPolyhedronBlock("polyhedron8", std::move(cells));

    const std::string p = mt::temp_path("_poly_scaling.vtu");
    meshioplusplus::write_vtu(p, m, /*binary=*/false, /*zlib=*/false);

    double best = std::numeric_limits<double>::infinity();
    for (int rep = 0; rep < 2; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        const meshioplusplus::Mesh back = meshioplusplus::read_vtu(p);
        const double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        EXPECT_EQ(back.NumCellBlocks(), 1u);
        EXPECT_EQ(back.Cells(0).NumCells(), NumCells);
        best = std::min(best, elapsed);
    }

    std::error_code ec;
    std::filesystem::remove(p, ec);
    return best;
}

}  // namespace

TEST(Vtu, PolyhedraWithAlternatingNodeCountsKeepTheirCellData) {
    // A polyhedral run mixing node counts is bucketed into one block per count, so a
    // bucket's cells are not contiguous in the file: wedge(1), tet(2), wedge(3) must read
    // back as polyhedron6 -> {1, 3} and polyhedron4 -> {2}, not the file-order slices
    // {1, 2} and {3} that bucket-by-bucket slicing used to hand out.
    auto wedge = [](std::int64_t b) {
        return std::vector<std::vector<std::int64_t>>{{b, b + 2, b + 1},
                                                      {b + 3, b + 4, b + 5},
                                                      {b, b + 1, b + 4, b + 3},
                                                      {b + 1, b + 2, b + 5, b + 4},
                                                      {b + 2, b, b + 3, b + 5}};
    };
    auto tet = [](std::int64_t b) {
        return std::vector<std::vector<std::int64_t>>{
            {b, b + 2, b + 1}, {b, b + 1, b + 3}, {b + 1, b + 2, b + 3}, {b + 2, b, b + 3}};
    };
    std::vector<std::vector<double>> pts;
    for (int i = 0; i < 16; ++i)
        pts.push_back({0.1 * i, (i % 3) * 0.5, (i % 5) * 0.25});

    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock("polyhedron6", {wedge(0)});
    m.AddPolyhedronBlock("polyhedron4", {tet(6)});
    m.AddPolyhedronBlock("polyhedron6", {wedge(10)});
    std::vector<meshioplusplus::NDArray> tags;
    for (double v : {1.0, 2.0, 3.0}) {
        meshioplusplus::NDArray blk(meshioplusplus::DType::Float64, {std::size_t{1}});
        blk.As<double>()[0] = v;
        tags.push_back(std::move(blk));
    }
    m.AddCellData("tag", std::move(tags));

    const std::string p = mt::temp_path("_poly_alternating.vtu");
    meshioplusplus::write_vtu(p, m, /*binary=*/false, /*zlib=*/false);
    const meshioplusplus::Mesh back = meshioplusplus::read_vtu(p);
    std::error_code ec;
    std::filesystem::remove(p, ec);

    ASSERT_EQ(back.NumCellBlocks(), 2u);
    ASSERT_TRUE(back.HasCellData("tag"));
    std::map<std::string, std::vector<double>> tags_by_type;
    for (std::size_t b = 0; b < back.NumCellBlocks(); ++b) {
        const meshioplusplus::NDArray& d = back.CellData("tag", b);
        ASSERT_EQ(d.Size(), back.Cells(b).NumCells());
        tags_by_type[std::string(back.Cells(b).Type())] =
            std::vector<double>(d.As<double>(), d.As<double>() + d.Size());
    }
    EXPECT_EQ(tags_by_type["polyhedron6"], (std::vector<double>{1.0, 3.0}));
    EXPECT_EQ(tags_by_type["polyhedron4"], (std::vector<double>{2.0}));
}

TEST(Vtu, PolyhedronReconstructionIsNotQuadraticInCellCount) {
    // Direct regression test for the roadmap's own probe: the pre-fix
    // reconstruct_cells rescanned `faceoffsets` from index 0 for EVERY cell, so
    // total work was O(ncells^2) regardless of how much face data there was. The
    // fix makes it O(ncells + |faces|).
    //
    // This asserts how the time SCALES, not how long it takes. It used to assert
    // an absolute 3 s budget for 100k cells, which is a statement about the
    // machine as much as the code: under the instrumented coverage build the
    // linear read already used ~90% of it (4.7 s of ctest wall time on master),
    // so the test failed intermittently on unrelated changes. Growing the cell
    // count 5x costs about 5x when reconstruction is linear and about 25x when
    // it is quadratic, whatever the machine or instrumentation; the threshold
    // sits between the two (their geometric mean is ~11).
    constexpr std::size_t kSmall = 20000;
    constexpr std::size_t kLarge = 100000;
    constexpr double kMaxRatio = 12.0;

    const double t_small = vtu_polyhedron_read_seconds(kSmall);
    const double t_large = vtu_polyhedron_read_seconds(kLarge);
    // A floor keeps the ratio finite if a fast machine resolves the small read as ~0.
    const double ratio = t_large / std::max(t_small, 1e-3);

    GTEST_LOG_(INFO) << "read " << kSmall << " cells in " << t_small << "s, " << kLarge << " in "
                     << t_large << "s: x" << ratio << " for x" << kLarge / kSmall << " the cells";
    EXPECT_LT(ratio, kMaxRatio) << "reading " << kLarge / kSmall << "x the cells took " << ratio
                                << "x as long (" << t_small << "s -> " << t_large
                                << "s); linear is ~" << kLarge / kSmall
                                << "x, so reconstruct_cells looks quadratic again";
}
