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
#include <filesystem>
#include <fstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/operations/stats.hpp"

namespace {

// EnSight writes a .case/.geo sibling pair; wrap mt::roundtrip so the
// sibling temp file is removed too.
void rt(const mt::Mesh& mesh, bool binary) {
    const double atol = binary ? 1e-6 : 1e-5;  // float32 vs %12.5e ASCII
    std::string sibling;
    mt::roundtrip(
        [&](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_ensight(p, m, binary);
            sibling = p.substr(0, p.size() - 5) + ".geo";
        },
        [](const std::string& p) { return meshioplusplus::read_ensight(p); }, mesh, ".case", atol);
    std::error_code ec;
    std::filesystem::remove(sibling, ec);
}

mt::Mesh wedge15_mesh() {
    return mt::make_mesh({{0, 0, 0},
                          {1, 0, 0},
                          {1, 1, 0},
                          {0, 0, 1},
                          {1, 0, 1},
                          {1, 1, 1},
                          {0.5, 0, 0},
                          {1, 0.5, 0},
                          {0.5, 0.5, 0},
                          {0.5, 0, 1},
                          {1, 0.5, 1},
                          {0.5, 0.5, 1},
                          {0, 0, 0.5},
                          {1, 0, 0.5},
                          {1, 1, 0.5}},
                         "wedge15", {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}});
}

}  // namespace

TEST(Ensight, AsciiLinear) {
    rt(mt::tri_mesh(), false);
    rt(mt::quad_mesh(), false);
    rt(mt::tet_mesh(), false);
    rt(mt::hex_mesh(), false);
    rt(mt::wedge_mesh(), false);
}
TEST(Ensight, AsciiQuadratic) {
    rt(mt::tet10_mesh(), false);
    rt(mt::hex20_mesh(), false);
    rt(wedge15_mesh(), false);
}
TEST(Ensight, AsciiHybrid) {
    rt(mt::tri_quad_mesh(), false);
}
TEST(Ensight, Ascii2D) {
    rt(mt::tri_mesh_2d(), false);
}
TEST(Ensight, BinaryLinear) {
    rt(mt::tri_mesh(), true);
    rt(mt::tet_mesh(), true);
    rt(mt::hex_mesh(), true);
    rt(mt::wedge_mesh(), true);
}
TEST(Ensight, BinaryQuadratic) {
    rt(mt::tet10_mesh(), true);
    rt(mt::hex20_mesh(), true);
    rt(wedge15_mesh(), true);
}

TEST(Ensight, AsciiKeywordsPresent) {
    std::string path = mt::temp_path(".case");
    meshioplusplus::write_ensight(path, mt::tet_mesh(), /*binary=*/false);
    const std::string geo_path = path.substr(0, path.size() - 5) + ".geo";
    std::ifstream in(geo_path, std::ios::binary);
    ASSERT_TRUE(in.good());
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("EnSight Gold Geometry File"), std::string::npos);
    EXPECT_NE(content.find("node id assign"), std::string::npos);
    EXPECT_NE(content.find("element id assign"), std::string::npos);
    EXPECT_NE(content.find("coordinates"), std::string::npos);
    EXPECT_NE(content.find("tetra4"), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(geo_path, ec);
}

TEST(Ensight, WriteRejectsUnknownType) {
    // No EnSight keyword for high-order Lagrange types -> WriteError.
    auto mesh =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}}, "line4", {{0, 1, 2, 3}});
    EXPECT_THROW(meshioplusplus::write_ensight(mt::temp_path(".case"), mesh, false),
                 meshioplusplus::WriteError);
}

// --- ragged (nsided / nfaced) round trips ------------------------------------
//
// The reader has always parsed these; the writer refused them until v9.19.0.
// EnSight is the cheapest of the polyhedral writers precisely because its wire
// format is a direct CSR dump -- no orientation contract, no global face table
// -- so a round trip against the existing reader is a complete oracle.

TEST(Ensight, NsidedPolygonRoundTrip) {
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0.5, 0}}));
    m.AddPolygonBlock("polygon", {{0, 1, 2, 3}, {1, 4, 2}});  // a quad then a triangle

    for (bool binary : {false, true}) {
        const std::string path = mt::temp_path(binary ? "_nsided_b.case" : "_nsided_a.case");
        meshioplusplus::write_ensight(path, m, binary);
        const mt::Mesh back = meshioplusplus::read_ensight(path);
        ASSERT_EQ(back.NumCellBlocks(), 1u) << "binary=" << binary;
        const auto cb = back.Cells(0);
        EXPECT_TRUE(cb.IsRagged());
        ASSERT_EQ(cb.NumCells(), 2u);
        EXPECT_EQ(cb.RowSize(0), 4u);
        EXPECT_EQ(cb.RowSize(1), 3u);
        EXPECT_EQ(cb.Row(1)[0], 1);
        EXPECT_EQ(cb.Row(1)[1], 4);
        EXPECT_EQ(cb.Row(1)[2], 2);
        std::error_code ec;
        std::filesystem::remove(path, ec);
        std::filesystem::remove(path.substr(0, path.size() - 5) + ".geo", ec);
    }
}

TEST(Ensight, NfacedPolyhedronRoundTrip) {
    mt::Mesh m;
    m.AssignPoints(mt::points_from(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}));
    m.AddPolyhedronBlock(
        "polyhedron8",
        {{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}}});

    for (bool binary : {false, true}) {
        const std::string path = mt::temp_path(binary ? "_nfaced_b.case" : "_nfaced_a.case");
        meshioplusplus::write_ensight(path, m, binary);
        const mt::Mesh back = meshioplusplus::read_ensight(path);
        ASSERT_EQ(back.NumCellBlocks(), 1u) << "binary=" << binary;
        const auto cb = back.Cells(0);
        EXPECT_TRUE(cb.IsPolyhedron());
        ASSERT_EQ(cb.NumCells(), 1u);
        EXPECT_EQ(cb.NumFaces(0), 6u);
        // Geometry, not just arity: the cube must come back as a unit cube.
        EXPECT_NEAR(meshioplusplus::compute_stats(back).mUnsignedVolume, 1.0, 1e-12);
        std::error_code ec;
        std::filesystem::remove(path, ec);
        std::filesystem::remove(path.substr(0, path.size() - 5) + ".geo", ec);
    }
}

TEST(Ensight, TransientVariablesReadPerNodeAndPerElement) {
    // roadmap §1 tier B1: static geometry + transient VARIABLE files -- the
    // common case the plan calls out. One file per step (EnSight's own
    // convention, unlike MED/CGNS/Tecplot/Gmsh's one-file-many-steps), so
    // the fixture is two scalar-per-node files and a TIME/VARIABLE section
    // appended to the .case write_ensight already produces.
    mt::Mesh m = mt::tri_mesh();
    const std::string path = mt::temp_path("_transient.case");
    meshioplusplus::write_ensight(path, m, /*binary=*/false);
    const std::string geo = path.substr(0, path.size() - 5) + ".geo";
    const std::string dir = path.substr(0, path.find_last_of("/\\") + 1);

    // Append TIME + VARIABLE to the .case file.
    {
        std::ofstream cf(path, std::ios::app);
        cf << "TIME\n"
           << "time set:              1\n"
           << "number of steps:       2\n"
           << "filename start number: 0\n"
           << "filename increment:    1\n"
           << "time values:\n"
           << "0.0\n"
           << "2.5\n"
           << "VARIABLE\n"
           << "scalar per node:    1  pressure  pressure.****.scl\n"
           << "scalar per element: 1  density   density.****.escl\n";
    }

    // Two per-node scalar files, four points each (mt::tri_mesh()), and two
    // per-element scalar files, two triangles each (one cell block).
    for (int step = 0; step < 2; ++step) {
        char name[64];
        std::snprintf(name, sizeof(name), "pressure.%04d.scl", step);
        std::ofstream vf(dir + name);
        vf << "pressure\n";
        vf << "part\n";
        vf << "         1\n";
        vf << "coordinates\n";
        const double base = step == 0 ? 10.0 : 11.0;
        for (int i = 0; i < 4; ++i)
            vf << (base + i * 10.0) << "\n";

        std::snprintf(name, sizeof(name), "density.%04d.escl", step);
        std::ofstream ef(dir + name);
        ef << "density\n";
        ef << "part\n";
        ef << "         1\n";
        ef << "tria3\n";
        const double ebase = step == 0 ? 100.0 : 200.0;
        ef << ebase << "\n" << (ebase + 1.0) << "\n";
    }

    meshioplusplus::ReadOptions opts;
    const meshioplusplus::MeshMetadata meta = meshioplusplus::read_ensight_metadata(path, opts);
    ASSERT_EQ(meta.mTimeValues.size(), 2u);
    EXPECT_DOUBLE_EQ(meta.mTimeValues[0], 0.0);
    EXPECT_DOUBLE_EQ(meta.mTimeValues[1], 2.5);
    EXPECT_TRUE(meta.mFellBackToFullRead);

    meshioplusplus::ReadOptions first;
    first.mTimeStep = 0;
    const mt::Mesh out0 = meshioplusplus::read_ensight(path, first);
    ASSERT_TRUE(out0.HasPointData("pressure"));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out0.PointData("pressure"), 0), 10.0);
    ASSERT_TRUE(out0.HasCellData("density"));
    ASSERT_EQ(out0.CellDataNumBlocks("density"), 1u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out0.CellData("density", 0), 0), 100.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out0.CellData("density", 0), 1), 101.0);

    meshioplusplus::ReadOptions second;
    second.mTimeStep = 1;
    const mt::Mesh out1 = meshioplusplus::read_ensight(path, second);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out1.PointData("pressure"), 0), 11.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out1.CellData("density", 0), 0), 200.0);

    meshioplusplus::ReadOptions last;
    last.mTimeStep = -1;
    const mt::Mesh out_last = meshioplusplus::read_ensight(path, last);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out_last.PointData("pressure"), 0), 11.0);

    meshioplusplus::ReadOptions too_far;
    too_far.mTimeStep = 5;
    EXPECT_THROW(meshioplusplus::read_ensight(path, too_far), meshioplusplus::ReadError);

    EXPECT_TRUE(meshioplusplus::seq_format_may_have_steps("ensight"));
    EXPECT_EQ(meshioplusplus::sequence_num_steps(path, "ensight"), 2u);

    // The plain (no-options) overload delegates to ReadOptions{}, which wants
    // every array by default -- same as every other format -- so it reads
    // VARIABLE too, at the default (first) step.
    const mt::Mesh plain = meshioplusplus::read_ensight(path);
    ASSERT_TRUE(plain.HasPointData("pressure"));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(plain.PointData("pressure"), 0), 10.0);

    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(geo, ec);
    std::filesystem::remove(dir + "pressure.0000.scl", ec);
    std::filesystem::remove(dir + "pressure.0001.scl", ec);
    std::filesystem::remove(dir + "density.0000.escl", ec);
    std::filesystem::remove(dir + "density.0001.escl", ec);
}

TEST(Ensight, TensorVariablesReadSwapsSymmLastTwoComponents) {
    // roadmap §1.1: EnSight Gold's own file order for `tensor symm` is
    // `11 22 33 12 13 23` (xx yy zz xy xz yz) -- confirmed empirically
    // against ParaView's vtkEnSightGoldReader (see doc/formats/ensight.md).
    // meshio++'s own convention is `xx yy zz xy yz zx`, the same six values
    // with the last two swapped. `tensor asym` (9 comp, row-major) needs no
    // such swap -- also confirmed against pvpython.
    mt::Mesh m = mt::tri_mesh();
    const std::string path = mt::temp_path("_tensor.case");
    meshioplusplus::write_ensight(path, m, /*binary=*/false);
    const std::string geo = path.substr(0, path.size() - 5) + ".geo";
    const std::string dir = path.substr(0, path.find_last_of("/\\") + 1);

    {
        std::ofstream cf(path, std::ios::app);
        cf << "VARIABLE\n"
           << "tensor symm per node:    1  sig  sig.sig\n"
           << "tensor asym per element: 1  g    g.tasym\n";
    }

    // Per-node tensor symm, mt::tri_mesh() has 4 points; file order is
    // component-major (every comp1, then every comp2, ...). Node 0's file
    // components are 100..600 (xx=100, yy=200, zz=300, xy=400, xz=500,
    // yz=600); expected meshio order after the swap is
    // [100,200,300,400,600,500].
    {
        std::ofstream vf(dir + "sig.sig");
        vf << "per-node tensor\n";
        vf << "part\n";
        vf << "         1\n";
        vf << "coordinates\n";
        for (int comp = 1; comp <= 6; ++comp)
            for (int node = 0; node < 4; ++node)
                vf << (comp * 100 + node) << "\n";
    }
    // Per-element tensor asym, one triangle block, 2 cells: cell 0's file
    // components are 11..19, straight through with no reordering.
    {
        std::ofstream ef(dir + "g.tasym");
        ef << "per-element tensor\n";
        ef << "part\n";
        ef << "         1\n";
        ef << "tria3\n";
        for (int comp = 1; comp <= 9; ++comp)
            for (int cell = 0; cell < 2; ++cell)
                ef << (comp * 10 + cell) << "\n";
    }

    const mt::Mesh out = meshioplusplus::read_ensight(path);
    ASSERT_TRUE(out.HasPointData("sig"));
    const meshioplusplus::NDArray& sig = out.PointData("sig");
    ASSERT_EQ(sig.Shape().size(), 2u);
    ASSERT_EQ(sig.Shape()[1], 6u);
    const double expect_node0[6] = {100, 200, 300, 400, 600, 500};
    for (int c = 0; c < 6; ++c)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(sig, static_cast<std::size_t>(c)),
                         expect_node0[c]);

    ASSERT_TRUE(out.HasCellData("g"));
    ASSERT_EQ(out.CellDataNumBlocks("g"), 1u);
    const meshioplusplus::NDArray& g = out.CellData("g", 0);
    ASSERT_EQ(g.Shape()[1], 9u);
    for (int c = 0; c < 9; ++c)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(g, static_cast<std::size_t>(c)),
                         (c + 1) * 10);

    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(geo, ec);
    std::filesystem::remove(dir + "sig.sig", ec);
    std::filesystem::remove(dir + "g.tasym", ec);
}

void run_variable_write_round_trip(bool binary) {
    // roadmap §1.1: write_ensight now writes a VARIABLE section -- one
    // scalar, one vector (padded from 2 to 3), one tensor symm point_data;
    // one tensor asym cell_data; and one field_data scalar as a "constant
    // per case". mt::tri_mesh() is 4 points, one tria3 block of 2 cells.
    mt::Mesh m = mt::tri_mesh();
    m.AddPointData("temp", mt::data_array({1.0, 2.0, 3.0, 4.0}));
    m.AddPointData("vel2d", mt::data_array({1.0, 10.0, 2.0, 20.0, 3.0, 30.0, 4.0, 40.0}, 2));
    m.AddPointData("sig", mt::data_array({11.0, 21.0, 31.0, 41.0, 51.0, 61.0, 12.0, 22.0, 32.0,
                                          42.0, 52.0, 62.0, 13.0, 23.0, 33.0, 43.0, 53.0, 63.0,
                                          14.0, 24.0, 34.0, 44.0, 54.0, 64.0},
                                         6));
    std::vector<meshioplusplus::NDArray> eps;
    eps.push_back(mt::data_array(
        {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0,
         80.0, 90.0},
        9));
    m.AddCellData("eps", std::move(eps));
    m.AddFieldData("gravity", mt::data_array({9.81}));

    const std::string path = mt::temp_path(binary ? "_wv_bin.case" : "_wv_ascii.case");
    meshioplusplus::write_ensight(path, m, binary);
    const std::string geo = path.substr(0, path.size() - 5) + ".geo";
    const std::string dir = path.substr(0, path.find_last_of("/\\") + 1);

    const mt::Mesh out = meshioplusplus::read_ensight(path);
    const double atol = binary ? 1e-6 : 1e-5;

    ASSERT_TRUE(out.HasPointData("temp"));
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_NEAR(meshioplusplus::detail::read_double(out.PointData("temp"), i),
                    static_cast<double>(i + 1), atol);

    ASSERT_TRUE(out.HasPointData("vel2d"));
    const meshioplusplus::NDArray& vel = out.PointData("vel2d");
    ASSERT_EQ(vel.Shape()[1], 3u);  // padded from 2 to 3
    EXPECT_NEAR(meshioplusplus::detail::read_double(vel, 0), 1.0, atol);
    EXPECT_NEAR(meshioplusplus::detail::read_double(vel, 1), 10.0, atol);
    EXPECT_NEAR(meshioplusplus::detail::read_double(vel, 2), 0.0, atol);  // padded

    ASSERT_TRUE(out.HasPointData("sig"));
    const meshioplusplus::NDArray& sig = out.PointData("sig");
    ASSERT_EQ(sig.Shape()[1], 6u);
    for (std::size_t r = 0; r < 4; ++r)
        for (std::size_t c = 0; c < 6; ++c)
            EXPECT_NEAR(meshioplusplus::detail::read_double(sig, r * 6 + c),
                        meshioplusplus::detail::read_double(m.PointData("sig"), r * 6 + c), atol);

    ASSERT_TRUE(out.HasCellData("eps"));
    ASSERT_EQ(out.CellDataNumBlocks("eps"), 1u);
    const meshioplusplus::NDArray& eps_out = out.CellData("eps", 0);
    ASSERT_EQ(eps_out.Shape()[1], 9u);
    for (std::size_t r = 0; r < 2; ++r)
        for (std::size_t c = 0; c < 9; ++c)
            EXPECT_NEAR(meshioplusplus::detail::read_double(eps_out, r * 9 + c),
                        meshioplusplus::detail::read_double(m.CellData("eps", 0), r * 9 + c),
                        atol);

    ASSERT_TRUE(out.HasFieldData("gravity"));
    EXPECT_NEAR(meshioplusplus::detail::read_double(out.FieldData("gravity"), 0), 9.81, atol);

    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(geo, ec);
    for (const char* name : {"temp.scl", "vel2d.vec", "sig.tsym", "eps.etasym"})
        std::filesystem::remove(dir + name, ec);
}

TEST(Ensight, VariableWriteRoundTripsAscii) { run_variable_write_round_trip(false); }

TEST(Ensight, VariableWriteRoundTripsBinary) { run_variable_write_round_trip(true); }

TEST(Ensight, VariableWriteSkipsUnsupportedComponentCountsWithAWarning) {
    // A 5-component array has no EnSight representation (only 1, 2, 3, 6, 9
    // are); it is skipped, not an error, and the rest of the write succeeds.
    mt::Mesh m = mt::tri_mesh();
    m.AddPointData("weird", mt::data_array({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                            16, 17, 18, 19, 20},
                                           5));
    m.AddPointData("ok", mt::data_array({1.0, 2.0, 3.0, 4.0}));

    const std::string path = mt::temp_path("_skip.case");
    meshioplusplus::write_ensight(path, m, false);
    const std::string geo = path.substr(0, path.size() - 5) + ".geo";
    const std::string dir = path.substr(0, path.find_last_of("/\\") + 1);

    const mt::Mesh out = meshioplusplus::read_ensight(path);
    EXPECT_FALSE(out.HasPointData("weird"));
    EXPECT_TRUE(out.HasPointData("ok"));

    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(geo, ec);
    std::filesystem::remove(dir + "ok.scl", ec);
}

TEST(Ensight, FortranBinaryRoundTripsGeometryAndVariables) {
    // Fortran binary frames every record as a Fortran sequential unformatted
    // record under a "Fortran Binary" header; the reader strips the markers
    // (detail/fortran_records.hpp) for the geometry and its variable files.
    mt::Mesh m = mt::tri_mesh();
    m.AddPointData("temp", mt::data_array({1.0, 2.0, 3.0, 4.0}));
    const std::string path = mt::temp_path("_fortran.case");
    meshioplusplus::write_ensight(path, m, /*binary=*/true, /*fortran=*/true);
    const std::string geo = path.substr(0, path.size() - 5) + ".geo";
    std::string head(18, '\0');
    {
        std::ifstream in(geo, std::ios::binary);
        in.read(head.data(), 18);
    }
    const std::int32_t marker = 80;
    EXPECT_EQ(head.substr(0, 4), std::string(reinterpret_cast<const char*>(&marker), 4));
    EXPECT_EQ(head.substr(4, 14), "Fortran Binary");

    const mt::Mesh out = meshioplusplus::read_ensight(path);
    mt::expect_mesh_eq(m, out, 1e-6);
    ASSERT_TRUE(out.HasPointData("temp"));
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_NEAR(meshioplusplus::detail::read_double(out.PointData("temp"), i),
                    static_cast<double>(i + 1), 1e-6);
    EXPECT_THROW(meshioplusplus::write_ensight(path, m, /*binary=*/false, /*fortran=*/true),
                 meshioplusplus::WriteError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(geo, ec);
    std::filesystem::remove(path.substr(0, path.find_last_of("/\\") + 1) + "temp.scl", ec);
}

TEST(Ensight, BinaryVariableFilesStartWithTheirDescription) {
    // EnSight defines no format record for a variable file; until v16.17.0 a
    // "C Binary" one came first, and files with it still read.
    mt::Mesh m = mt::tri_mesh();
    m.AddPointData("temp", mt::data_array({1.0, 2.0, 3.0, 4.0}));
    const std::string path = mt::temp_path("_cbinvar.case");
    meshioplusplus::write_ensight(path, m, /*binary=*/true);
    const std::string dir = path.substr(0, path.find_last_of("/\\") + 1);
    const std::string var = dir + "temp.scl";
    std::string bytes;
    {
        std::ifstream in(var, std::ios::binary);
        bytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    ASSERT_GE(bytes.size(), 160u);
    EXPECT_EQ(bytes.substr(0, 8), "variable");
    EXPECT_EQ(bytes.substr(80, 4), "part");

    std::string legacy(80, '\0');
    legacy.replace(0, 8, "C Binary");
    {
        std::ofstream out(var, std::ios::binary);
        out << legacy << bytes;
    }
    const mt::Mesh back = meshioplusplus::read_ensight(path);
    ASSERT_TRUE(back.HasPointData("temp"));
    EXPECT_NEAR(meshioplusplus::detail::read_double(back.PointData("temp"), 3), 4.0, 1e-6);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(path.substr(0, path.size() - 5) + ".geo", ec);
    std::filesystem::remove(var, ec);
}
