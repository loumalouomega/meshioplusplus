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
// Tests for content-based format sniffing (conservative signature matching).

// System includes
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/operations/sniff.hpp"

namespace {

std::string write_temp(const std::string& rContents, const char* pSuffix) {
    std::string path = std::string(std::tmpnam(nullptr)) + pSuffix;
    std::ofstream out(path, std::ios::binary);
    out.write(rContents.data(), static_cast<std::streamsize>(rContents.size()));
    out.close();
    return path;
}

}  // namespace

TEST(Sniff, RecognizesKnownSignatures) {
    struct Case {
        std::string contents;
        std::string expected;
    };
    const Case cases[] = {
        {"# vtk DataFile Version 3.0\n", "vtk"},
        {"$MeshFormat\n2.2 0 8\n", "gmsh"},
        {"ply\nformat ascii 1.0\n", "ply"},
        {"OFF\n8 6 0\n", "off"},
        {"solid mysolid\n facet normal 0 0 1\n", "stl"},
        {"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\">", "vtu"},
        {"<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\">", "vtp"},
        {"*Heading\n test\n*Node\n", "abaqus"},
        {"*KEYWORD\n*NODE\n", "lsdyna"},
        {"$ a comment\n*keyword long=y\n*NODE\n", "lsdyna"},
        {"% written by Salome\nTITRE\n mesh\nFINSF\n", "code_aster"},
        {"coor_3d NOM=INDEFINI\n N1 0 0 0\nFINSF\nFIN\n", "code_aster"},
        {"    1C\n    1UCALCULIX\n    2C\n", "frd"},
        {"    1C\r\n    2C                            20\r\n", "frd"},
        {"# Created by COMSOL\n0 1 \n1 # number of tags\n5 mesh1 \n", "mphtxt"},
        {std::string("\0\0\0\0\1\0\0\0\1\0\0\0\5\0\0\0m\0\0\0", 20), "mphbin"},
        {"<?xml version=\"1.0\"?>\n<febio_spec version=\"4.0\">", "febio"},
        {std::string("BEF\0\0\0\0\1", 8), "xplt"},
        {"25       0       0       1       0       0       0       0       0\ntitle\n", "patran"},
        {"26       0       0       1      33       7       1       2       0\r\n", "patran"},
        {"25 this is not a header card at all, just text\n", ""},
        {"   -1\n   100\n<NULL>\n8.2,\n   -1\n", "femap"},
        {"323152540\r\n   -1\r\n   100\r\n", "femap"},
        {"   -1\n   151\n", "unv"},
        {"MFEM mesh v1.0\n\ndimension\n2\n", "mfem"},
        {"MFEM NC mesh v1.0\n", "mfem"},
        {"libMesh-1.3.0\n2\t # number of elements\n", "libmesh"},
        {std::string("\0\0\0\x0dlibMesh-1.8.0\0\0\0", 20), "libmesh"},
        {"*I 19I 41921A6.23-1  A07-Nov-2A024     A16:50:01I 11I 18", "abaqus_fil"},
        {std::string("\0\x10\0\0\x09\0\0\0\0\0\0\0\x81\x07\0\0\0\0\0\0", 20), "abaqus_fil"},
        {std::string("\0\0\x10\0\0\0\0\0\0\0\0\x09\0\0\0\0\0\0\x07\x81", 20), "abaqus_fil"},
        {"#RADIOSS STARTER\n/BEGIN\nrun\n", "radioss"},
        {"# a comment\n$ another\n/BEGIN\nrun\n", "radioss"},
        {"title               pull\n$ a comment\nextended\nsizing 0 1 8\nend\n", "marc"},
        {"=beg=50100 (Analysis Title)\n          job1\n=end=\n", "marc_t19"},
        {"title               pull\nsizing 0 1 8\n", ""},  // no END or model option
        {"TITLE = \"a tecplot file\"\nVARIABLES = x y\nend\n", ""},
        {"*I am not a results file\n", ""},
        {"/NODE\n1 0 0 0\n", ""},
        {"25       0       0       0       0       0       0       0       0\n", ""},
        {"0 1 2 3\n", ""},
        {"    1C\n", ""},
        {"    1C\nsomething else\n", ""},
    };
    for (const Case& c : cases) {
        const std::string path = write_temp(c.contents, ".dat");
        EXPECT_EQ(meshioplusplus::sniff_format(path), c.expected) << "contents: " << c.contents;
        std::remove(path.c_str());
    }
}

TEST(Sniff, ReturnsEmptyOnAmbiguousOrUnknown) {
    // A generic HDF5 magic must not be claimed (med/h5m/cgns/hmf are ambiguous).
    const std::string hdf5 = write_temp(std::string("\x89HDF\r\n\x1a\n----", 12), ".dat");
    EXPECT_EQ(meshioplusplus::sniff_format(hdf5), "");
    std::remove(hdf5.c_str());

    const std::string garbage = write_temp("just some random text\n", ".dat");
    EXPECT_EQ(meshioplusplus::sniff_format(garbage), "");
    std::remove(garbage.c_str());

    EXPECT_EQ(meshioplusplus::sniff_format("/nonexistent/path/xyz.dat"), "");
}

namespace {

// A fresh, empty temporary directory.
std::filesystem::path make_temp_dir(const std::string& rName) {
    const std::filesystem::path dir = std::string(std::tmpnam(nullptr)) + "_" + rName;
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    return dir;
}

void touch(const std::filesystem::path& rPath) {
    std::filesystem::create_directories(rPath.parent_path());
    std::ofstream(rPath) << "x\n";
}

}  // namespace

TEST(Sniff, RecognizesDirectoryFormatsByTheirFiles) {
    namespace fs = std::filesystem;
    const fs::path elmer = make_temp_dir("elmer");
    touch(elmer / "mesh.header");
    EXPECT_EQ(meshioplusplus::sniff_format(elmer.string()), "elmer");
    // The header itself stands for its directory.
    EXPECT_EQ(meshioplusplus::sniff_format((elmer / "mesh.header").string()), "elmer");

    const fs::path parts = make_temp_dir("parts") / "partitioning.2";
    touch(parts / "part.1.header");
    EXPECT_EQ(meshioplusplus::sniff_format(parts.string()), "elmer");
    EXPECT_EQ(meshioplusplus::sniff_format(parts.parent_path().string()), "elmer");

    const fs::path foam = make_temp_dir("foam");
    touch(foam / "constant" / "polyMesh" / "owner");
    touch(foam / "constant" / "polyMesh" / "faces");
    EXPECT_EQ(meshioplusplus::sniff_format(foam.string()), "openfoam");
    EXPECT_EQ(meshioplusplus::sniff_format((foam / "constant" / "polyMesh").string()), "openfoam");

    const fs::path decomposed = make_temp_dir("decomposed");
    fs::create_directories(decomposed / "processor0" / "constant" / "polyMesh");
    EXPECT_EQ(meshioplusplus::sniff_format(decomposed.string()), "openfoam");

    // Both, or neither, is not a guess.
    const fs::path both = make_temp_dir("both");
    touch(both / "mesh.header");
    touch(both / "polyMesh" / "owner");
    touch(both / "polyMesh" / "faces");
    EXPECT_EQ(meshioplusplus::sniff_format(both.string()), "");
    const fs::path empty = make_temp_dir("empty");
    EXPECT_EQ(meshioplusplus::sniff_format(empty.string()), "");
    // A polyMesh without its faces file is not a case.
    const fs::path partial = make_temp_dir("partial");
    touch(partial / "polyMesh" / "owner");
    EXPECT_EQ(meshioplusplus::sniff_format(partial.string()), "");

    for (const fs::path& rDir :
         {elmer, parts.parent_path(), foam, decomposed, both, empty, partial})
        fs::remove_all(rDir);
}
