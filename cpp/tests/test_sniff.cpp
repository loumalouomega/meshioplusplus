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
