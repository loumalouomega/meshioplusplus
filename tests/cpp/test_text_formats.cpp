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
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/formats/avsucd.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/netgen.hpp"
#include "meshioplusplus/formats/permas.hpp"
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/formats/ugrid.hpp"

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
