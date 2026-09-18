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
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/formats/avsucd.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/netgen.hpp"
#include "meshioplusplus/formats/permas.hpp"
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/formats/ugrid.hpp"
#include "meshioplusplus/operations/sequence.hpp"

// Generic helper for `void write_X(path, mesh)` / `Mesh read_X(path)`.
#define SIMPLE_RT(WRITER, READER, MESH, SUFFIX, ATOL)                            \
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { WRITER(p, m); }, \
                  [](const std::string& p) { return READER(p); }, MESH, SUFFIX, ATOL)

TEST(Medit, Basic) {
    SIMPLE_RT(meshioplusplus::write_medit_ascii, meshioplusplus::read_medit_ascii, mt::tri_mesh(),
              ".mesh", 1e-12);
    SIMPLE_RT(meshioplusplus::write_medit_ascii, meshioplusplus::read_medit_ascii, mt::tet_mesh(),
              ".mesh", 1e-12);
    SIMPLE_RT(meshioplusplus::write_medit_ascii, meshioplusplus::read_medit_ascii, mt::hex_mesh(),
              ".mesh", 1e-12);
    SIMPLE_RT(meshioplusplus::write_medit_ascii, meshioplusplus::read_medit_ascii, mt::quad_mesh(),
              ".mesh", 1e-12);
}

TEST(Nastran, Basic) {
    // The C++ reader is sentinel-gated to meshio-written files, so a self
    // round-trip is the supported path.
    SIMPLE_RT(meshioplusplus::write_nastran, meshioplusplus::read_nastran, mt::tri_mesh(), ".bdf",
              1e-10);
    SIMPLE_RT(meshioplusplus::write_nastran, meshioplusplus::read_nastran, mt::tet_mesh(), ".bdf",
              1e-10);
    SIMPLE_RT(meshioplusplus::write_nastran, meshioplusplus::read_nastran, mt::hex_mesh(), ".bdf",
              1e-10);
}

TEST(Abaqus, Basic) {
    SIMPLE_RT(meshioplusplus::write_abaqus, meshioplusplus::read_abaqus, mt::tri_mesh(), ".inp",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_abaqus, meshioplusplus::read_abaqus, mt::tet_mesh(), ".inp",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_abaqus, meshioplusplus::read_abaqus, mt::hex_mesh(), ".inp",
              1e-12);
}

TEST(Avsucd, Basic) {
    SIMPLE_RT(meshioplusplus::write_avsucd, meshioplusplus::read_avsucd, mt::tri_mesh(), ".avs",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_avsucd, meshioplusplus::read_avsucd, mt::tet_mesh(), ".avs",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_avsucd, meshioplusplus::read_avsucd, mt::hex_mesh(), ".avs",
              1e-12);
}

TEST(Avsucd, WarnsWhenDemotingAnExtraIntegerCellDataArray) {
    // AVS-UCD has exactly one integer "material id" column; a second integer
    // cell_data array is not dropped (it still lands in the generic
    // cell-data section, demoted to real-valued) but that demotion used to
    // happen with no diagnostic at all.
    mt::Mesh m = mt::tri_mesh();
    const std::size_t nc = m.Cells(0).NumCells();
    meshioplusplus::NDArray a(meshioplusplus::DType::Int32, {nc});
    meshioplusplus::NDArray b(meshioplusplus::DType::Int32, {nc});
    for (std::size_t i = 0; i < nc; ++i) {
        a.As<std::int32_t>()[i] = static_cast<std::int32_t>(i);
        b.As<std::int32_t>()[i] = static_cast<std::int32_t>(2 * i);
    }
    m.AppendCellData("first_material", std::move(a));
    m.AppendCellData("second_material", std::move(b));

    const std::string p = mt::temp_path(".avs");
    testing::internal::CaptureStderr();
    meshioplusplus::write_avsucd(p, m);
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_NE(err.find("first_material"), std::string::npos) << err;
    EXPECT_NE(err.find("second_material"), std::string::npos) << err;

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Permas, Basic) {
    SIMPLE_RT(meshioplusplus::write_permas, meshioplusplus::read_permas, mt::tri_mesh(), ".post",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_permas, meshioplusplus::read_permas, mt::tri_quad_mesh(),
              ".post", 1e-12);
    SIMPLE_RT(meshioplusplus::write_permas, meshioplusplus::read_permas, mt::hex_mesh(), ".post",
              1e-12);
}

TEST(Tecplot, SingleType) {
    SIMPLE_RT(meshioplusplus::write_tecplot, meshioplusplus::read_tecplot, mt::tri_mesh(), ".dat",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_tecplot, meshioplusplus::read_tecplot, mt::quad_mesh(), ".dat",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_tecplot, meshioplusplus::read_tecplot, mt::tet_mesh(), ".dat",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_tecplot, meshioplusplus::read_tecplot, mt::hex_mesh(), ".dat",
              1e-12);
}

TEST(Ugrid, Basic) {
    SIMPLE_RT(meshioplusplus::write_ugrid, meshioplusplus::read_ugrid, mt::tri_mesh(), ".ugrid",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_ugrid, meshioplusplus::read_ugrid, mt::tet_mesh(), ".ugrid",
              1e-12);
    SIMPLE_RT(meshioplusplus::write_ugrid, meshioplusplus::read_ugrid, mt::hex_mesh(), ".ugrid",
              1e-12);
}

TEST(Ugrid, BinaryFlavours) {
    // The flavour (endianness, float/int width, Fortran record markers) is
    // decoded from the penultimate filename suffix, so round-trip through each
    // to exercise the byte-swap, width, and Fortran-record branches that the
    // ASCII ".ugrid" path never touches.
    for (const char* suffix :
         {".b8.ugrid", ".b4.ugrid", ".lb8.ugrid", ".lb4.ugrid", ".r8.ugrid", ".lr8.ugrid"}) {
        SIMPLE_RT(meshioplusplus::write_ugrid, meshioplusplus::read_ugrid, mt::tet_mesh(), suffix,
                  1e-12);
        SIMPLE_RT(meshioplusplus::write_ugrid, meshioplusplus::read_ugrid, mt::hex_mesh(), suffix,
                  1e-12);
    }
}

TEST(Ugrid, ReadRejectsTruncatedBinary) {
    // A binary-flavour file shorter than its fixed count header must raise
    // rather than read past EOF.
    std::string path = mt::temp_path(".lb8.ugrid");
    {
        std::ofstream f(path, std::ios::binary);
        const char partial[4] = {0, 0, 0, 0};
        f.write(partial, sizeof(partial));
    }
    EXPECT_THROW(meshioplusplus::read_ugrid(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Netgen, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) {
        meshioplusplus::write_netgen(p, m, ".16e");
    };
    auto r = [](const std::string& p) { return meshioplusplus::read_netgen(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".vol");
    mt::roundtrip(w, r, mt::tet_mesh(), ".vol");
    mt::roundtrip(w, r, mt::hex_mesh(), ".vol");
}

// --- Malformed-input reject paths ---
// These drive the readers' explicit `throw ReadError` branches that self
// round-trips never reach (medit.cpp:230, tecplot.cpp:185, nastran.cpp:265,
// netgen.cpp:233).

namespace {
std::string write_temp(const std::string& suffix, const std::string& contents) {
    std::string path = mt::temp_path(suffix);
    std::ofstream f(path);
    f << contents;
    return path;
}
}  // namespace

TEST(Medit, ReadRejectsMissingVertices) {
    std::string path = write_temp(".mesh", "MeshVersionFormatted 2\nDimension 3\n");
    EXPECT_THROW(meshioplusplus::read_medit_ascii(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tecplot, TransientReadSelectsOneZoneBySolutionTime) {
    // roadmap §1 tier B1: two FEBLOCK zones sharing STRANDID=1, distinct
    // SOLUTIONTIME. Hand-written text -- Tecplot ASCII needs no external
    // tool, unlike MED/CGNS/cgnslib.
    const std::string contents =
        "VARIABLES = \"X\" \"Y\" \"u\"\n"
        "ZONE N=3 E=1 DATAPACKING=BLOCK ZONETYPE=FETRIANGLE SOLUTIONTIME=0.0 STRANDID=1\n"
        "0.0 1.0 0.0\n"
        "0.0 0.0 1.0\n"
        "10.0 20.0 30.0\n"
        "1 2 3\n"
        "ZONE N=3 E=1 DATAPACKING=BLOCK ZONETYPE=FETRIANGLE SOLUTIONTIME=2.5 STRANDID=1\n"
        "0.0 1.0 0.0\n"
        "0.0 0.0 1.0\n"
        "11.0 21.0 31.0\n"
        "1 2 3\n";
    std::string path = write_temp(".dat", contents);

    meshioplusplus::ReadOptions opts;
    const meshioplusplus::MeshMetadata meta = meshioplusplus::read_tecplot_metadata(path, opts);
    ASSERT_EQ(meta.mTimeValues.size(), 2u);
    EXPECT_DOUBLE_EQ(meta.mTimeValues[0], 0.0);
    EXPECT_DOUBLE_EQ(meta.mTimeValues[1], 2.5);
    EXPECT_EQ(meta.mNumPoints, 3u);
    ASSERT_EQ(meta.mCellBlocks.size(), 1u);
    EXPECT_EQ(meta.mCellBlocks[0].mType, "triangle");
    EXPECT_EQ(meta.mCellBlocks[0].mNumCells, 1u);

    meshioplusplus::ReadOptions first;
    first.mTimeStep = 0;
    const mt::Mesh out0 = meshioplusplus::read_tecplot(path, first);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out0.PointData("u"), 0), 10.0);

    meshioplusplus::ReadOptions second;
    second.mTimeStep = 1;
    const mt::Mesh out1 = meshioplusplus::read_tecplot(path, second);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out1.PointData("u"), 0), 11.0);

    meshioplusplus::ReadOptions last;
    last.mTimeStep = -1;
    const mt::Mesh out_last = meshioplusplus::read_tecplot(path, last);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out_last.PointData("u"), 0), 11.0);

    meshioplusplus::ReadOptions too_far;
    too_far.mTimeStep = 5;
    EXPECT_THROW(meshioplusplus::read_tecplot(path, too_far), meshioplusplus::ReadError);

    EXPECT_TRUE(meshioplusplus::seq_format_may_have_steps("tecplot"));
    EXPECT_EQ(meshioplusplus::sequence_num_steps(path, "tecplot"), 2u);

    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tecplot, NonTransientMultiZoneWarnsAndReadsTheFirst) {
    // No SOLUTIONTIME anywhere: the "several static zones" case is not a
    // timeline (roadmap §7 owns concatenating them); only the first zone is
    // read, as before -- now with a warning rather than silent truncation.
    const std::string contents =
        "VARIABLES = \"X\" \"Y\"\n"
        "ZONE N=3 E=1 DATAPACKING=BLOCK ZONETYPE=FETRIANGLE\n"
        "0.0 1.0 0.0\n"
        "0.0 0.0 1.0\n"
        "1 2 3\n"
        "ZONE N=3 E=1 DATAPACKING=BLOCK ZONETYPE=FETRIANGLE\n"
        "2.0 3.0 2.0\n"
        "2.0 2.0 3.0\n"
        "1 2 3\n";
    std::string path = write_temp(".dat", contents);

    const mt::Mesh out = meshioplusplus::read_tecplot(path);
    EXPECT_EQ(out.NumPoints(), 3u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out.Points(), 0), 0.0);

    meshioplusplus::ReadOptions opts;
    const meshioplusplus::MeshMetadata meta = meshioplusplus::read_tecplot_metadata(path, opts);
    EXPECT_TRUE(meta.mTimeValues.empty());

    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tecplot, ReadRejectsMissingVariables) {
    std::string path = write_temp(".dat", "ZONE N=1 E=1\n");
    EXPECT_THROW(meshioplusplus::read_tecplot(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Nastran, ReadRejectsForeignFile) {
    // The C++ reader is sentinel-gated to meshio++-written files.
    std::string path = write_temp(".bdf", "$ some other tool's bulk data\nGRID,1,,0.,0.,0.\n");
    EXPECT_THROW(meshioplusplus::read_nastran(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Netgen, ReadRejectsInvalidFile) {
    std::string path = write_temp(".vol", "this is not a netgen mesh\n");
    EXPECT_THROW(meshioplusplus::read_netgen(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
