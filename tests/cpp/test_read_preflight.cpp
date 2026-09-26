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
/**
 * @file test_read_preflight.cpp
 * @brief Readers refuse a file they cannot read from its first bytes, before
 *        parsing the rest (roadmap §4, "A declined C++ read"): each input
 *        below has a refusing header and a body the full parse would reject
 *        with a different message, so the message shows which check ran.
 */

#include <cstdio>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"

namespace {

std::string rp_write(const std::string& rSuffix, const std::string& rText) {
    const std::string path = mt::temp_path(rSuffix);
    std::ofstream(path, std::ios::binary) << rText;
    return path;
}

template <class F>
std::string rp_message(F&& read) {
    try {
        read();
    } catch (const meshioplusplus::ReadError& e) {
        return e.what();
    }
    return "(no error)";
}

}  // namespace

TEST(ReadPreflight, VtkXmlRefusesLzmaFromTheStartTag) {
    // The start tag names lzma; the rest is not even well-formed XML, which a
    // full parse would report first.
    const std::string body =
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\"\n"
        "  version=\"1.0\" compressor='vtkLZMADataCompressor'>\n<<< broken";
    const std::string vtu = rp_write(".vtu", body);
    EXPECT_NE(rp_message([&] { meshioplusplus::read_vtu(vtu); }).find("lzma-compressed VTU"),
              std::string::npos);
    std::remove(vtu.c_str());

    std::string poly = body;
    poly.replace(poly.find("UnstructuredGrid"), 16, "PolyData");
    const std::string vtp = rp_write(".vtp", poly);
    EXPECT_NE(rp_message([&] { meshioplusplus::read_vtp(vtp); }).find("lzma-compressed VTP"),
              std::string::npos);
    std::remove(vtp.c_str());
}

TEST(ReadPreflight, VtkXmlLeavesAnotherFileTypeToTheFullParse) {
    // A PolyData start tag in a .vtu: the type check, not the pre-flight, answers.
    const std::string path =
        rp_write(".vtu",
                 "<VTKFile type=\"PolyData\" compressor=\"vtkLZMADataCompressor\">\n"
                 "</VTKFile>\n");
    EXPECT_NE(rp_message([&] { meshioplusplus::read_vtu(path); }).find("UnstructuredGrid"),
              std::string::npos);
    std::remove(path.c_str());
}

TEST(ReadPreflight, GmshRefusesPeriodicBeforeParsingNodes) {
    // $Nodes is malformed; $Periodic comes after it, as Gmsh writes it.
    const std::string path = rp_write(".msh",
                                      "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
                                      "$Nodes\nnot-a-count\n$EndNodes\n"
                                      "$Periodic  \n0\n$EndPeriodic\n");
    EXPECT_NE(rp_message([&] { meshioplusplus::read_gmsh(path); }).find("$Periodic"),
              std::string::npos);
    std::remove(path.c_str());
}

TEST(ReadPreflight, GmshIgnoresPeriodicInsideALine) {
    // "$Periodic" that is not a section header decides nothing.
    const std::string path = rp_write(".msh",
                                      "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
                                      "$PhysicalNames\n1\n2 1 \"x$Periodic\"\n$EndPhysicalNames\n"
                                      "$Nodes\n3\n1 0 0 0\n2 1 0 0\n3 0 1 0\n$EndNodes\n"
                                      "$Elements\n1\n1 2 0 1 2 3\n$EndElements\n");
    EXPECT_EQ(meshioplusplus::read_gmsh(path).NumPoints(), 3u);
    std::remove(path.c_str());
}

TEST(ReadPreflight, FluentRefusesAGmshFileAtItsFirstByte) {
    const std::string path =
        rp_write(".msh", "  \n$MeshFormat\n2.2 0 8\n$EndMeshFormat\n" + std::string(4096, 'x'));
    EXPECT_EQ(rp_message([&] { meshioplusplus::read_ansys(path); }),
              "Fluent: expected a section at byte 3");
    const std::string no_index = rp_write(".msh", "( x)\n");
    EXPECT_EQ(rp_message([&] { meshioplusplus::read_ansys(no_index); }),
              "Fluent: expected a section at byte 0");
    std::remove(path.c_str());
    std::remove(no_index.c_str());
}
