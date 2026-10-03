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
#include <string>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtk.hpp"

namespace {
void rt(const mt::Mesh& mesh, bool binary, bool v51) {
    mt::roundtrip([=](const std::string& p,
                      const mt::Mesh& m) { meshioplusplus::write_vtk(p, m, binary, v51); },
                  [](const std::string& p) { return meshioplusplus::read_vtk(p); }, mesh, ".vtk");
}
}  // namespace

TEST(Vtk, V51Ascii) {
    rt(mt::tri_mesh(), false, true);
    rt(mt::tet_mesh(), false, true);
    rt(mt::hex_mesh(), false, true);
}
TEST(Vtk, V51Binary) {
    rt(mt::tri_mesh(), true, true);
    rt(mt::tet_mesh(), true, true);
}
TEST(Vtk, V42Ascii) {
    rt(mt::tri_mesh(), false, false);
    rt(mt::tet_mesh(), false, false);
}
TEST(Vtk, V42Binary) {
    rt(mt::quad_mesh(), true, false);
}
TEST(Vtk, Hybrid) {
    rt(mt::tri_quad_mesh(), false, true);
}

TEST(Vtk, ReadRejectsStructuredGridWithoutDimensions) {
    // Recognizing DATASET alone is not enough: geometry is required.
    std::string path = mt::temp_path(".vtk");
    {
        std::ofstream f(path);
        f << "# vtk DataFile Version 3.0\n"
          << "structured grid header\n"
          << "ASCII\n"
          << "DATASET STRUCTURED_GRID\n";
    }
    EXPECT_THROW(meshioplusplus::read_vtk(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Vtk, RefusesOutOfRangeFloatingTopology) {
    for (const char* dtype : {"float", "double"}) {
        for (const char* value : {"nan", "inf", "-inf", "9223372036854775808", "-1e30"}) {
            for (bool bad_offset : {false, true}) {
                SCOPED_TRACE(std::string(dtype) + " " + value +
                             (bad_offset ? " offsets" : " connectivity"));
                const std::string path = mt::temp_path(".vtk");
                {
                    std::ofstream f(path);
                    f << "# vtk DataFile Version 5.1\n"
                      << "invalid floating topology\nASCII\nDATASET UNSTRUCTURED_GRID\n"
                      << "POINTS 1 float\n0 0 0\nCELLS 2 1\nOFFSETS " << dtype << "\n0 "
                      << (bad_offset ? value : "1") << "\nCONNECTIVITY " << dtype << "\n"
                      << (bad_offset ? "0" : value) << "\nCELL_TYPES 1\n1\n";
                }
                try {
                    (void)meshioplusplus::read_vtk(path);
                    FAIL() << "An out-of-range topology value must be refused";
                } catch (const meshioplusplus::ReadError& e) {
                    EXPECT_NE(std::string(e.what()).find("integer field out of range"),
                              std::string::npos);
                }
                std::error_code ec;
                std::filesystem::remove(path, ec);
            }
        }
    }
}
