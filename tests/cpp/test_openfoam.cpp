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
#include <atomic>
#include <filesystem>
#include <fstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/formats/openfoam.hpp"

namespace fs = std::filesystem;

namespace {

// Write a minimal single-hex ASCII polyMesh case; return the case dir.
fs::path make_hex_case() {
    static std::atomic<unsigned> counter{0};
    fs::path base = fs::temp_directory_path() / ("meshio_of_" + std::to_string(counter++));
    fs::path poly = base / "constant" / "polyMesh";
    fs::create_directories(poly);

    auto hdr = [](const std::string& cls, const std::string& obj) {
        return "FoamFile\n{\n format ascii;\n class " + cls + ";\n object " + obj + ";\n}\n";
    };

    std::ofstream(poly / "points")
        << hdr("vectorField", "points")
        << "8\n(\n(0 0 0)\n(1 0 0)\n(1 1 0)\n(0 1 0)\n(0 0 1)\n(1 0 1)\n(1 1 "
           "1)\n(0 1 1)\n)\n";
    std::ofstream(poly / "faces")
        << hdr("faceList", "faces")
        << "6\n(\n4(0 3 2 1)\n4(4 5 6 7)\n4(0 1 5 4)\n4(2 3 7 6)\n4(1 2 6 "
           "5)\n4(0 4 7 3)\n)\n";
    std::ofstream(poly / "owner") << hdr("labelList", "owner") << "6\n(\n0\n0\n0\n0\n0\n0\n)\n";
    std::ofstream(poly / "boundary")
        << hdr("polyBoundaryMesh", "boundary")
        << "3\n(\nbottom { type wall; nFaces 1; startFace 0; }\ntop { type "
           "wall; nFaces 1; startFace 1; }\nsides { type wall; nFaces 4; "
           "startFace 2; }\n)\n";

    std::ofstream(base / "case.foam") << "";
    return base;
}

}  // namespace

TEST(OpenFoam, SingleHexAscii) {
    fs::path base = make_hex_case();
    meshioplusplus::OpenFoamInfo info;
    meshioplusplus::Mesh mesh = meshioplusplus::read_openfoam((base / "case.foam").string(), info);

    EXPECT_EQ(mesh.NumPoints(), 8u);
    bool has_hex = false;
    std::size_t nquad = 0;
    for (const auto cb : mesh.CellRange()) {
        if (cb.Type() == "hexahedron")
            has_hex = true;
        if (cb.Type() == "quad")
            nquad += cb.NumCells();
    }
    EXPECT_TRUE(has_hex);
    EXPECT_EQ(nquad, 6u);
    // 3 boundary patches -> 3 negative family tags.
    EXPECT_EQ(info.mCellTags.size(), 3u);

    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoam, ResolveViaCaseDir) {
    fs::path base = make_hex_case();
    meshioplusplus::OpenFoamInfo info;
    // Pass the case directory itself (not the .foam file).
    meshioplusplus::Mesh mesh = meshioplusplus::read_openfoam(base.string(), info);
    EXPECT_EQ(mesh.NumPoints(), 8u);
    std::error_code ec;
    fs::remove_all(base, ec);
}
