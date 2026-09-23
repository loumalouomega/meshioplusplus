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
 * @file test_femap.cpp
 * @brief Femap neutral reader/writer: a hand-written file with a tetra10 in
 *        Femap's slot layout, properties, groups, node lists, output sets and
 *        both output-vector encodings; the error paths; a round trip.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/femap.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;
using meshioplusplus::RegionKind;
namespace detail = meshioplusplus::detail;

std::string write_file(const std::string& rBody) {
    const std::string path = mt::temp_path(".neu");
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

// A unit tetra10 whose nodes sit in Femap's degenerate-brick slots: corners in
// slots 0, 1, 2 and 4 (the apex); mid-edges 0-1, 1-2, 2-0 in slots 8-10 and the
// apex edges in 12-14. Node ids step by 10. A MYSTRAN-style junk line precedes
// the first block; the file is 9.3, so a rigid element's node list is skipped.
std::string tet10_file() {
    const double xyz[10][3] = {{0, 0, 0},   {1, 0, 0},  {0, 1, 0},  {0, 0, 1},   {.5, 0, 0},
                               {.5, .5, 0}, {0, .5, 0}, {0, 0, .5}, {.5, 0, .5}, {0, .5, .5}};
    std::string out = "323152540\n   -1\n   100\n<NULL>\n9.3,\n   -1\n";
    out +=
        "   -1\n   402\n5,110,1,25,1,0,\nSOLID PART\n0,0,0,0,\n2,\n0,0,\n3,\n1.,2.,3.,\n"
        "0,\n0,\n   -1\n";
    out += "   -1\n   403\n";
    for (int k = 0; k < 10; ++k)
        out += std::to_string(10 * (k + 1)) + ",0,0,1,46,0,0,0,0,0,0," + std::to_string(xyz[k][0]) +
               "," + std::to_string(xyz[k][1]) + "," + std::to_string(xyz[k][2]) + ",0,1,\n";
    out += "   -1\n   -1\n   404\n";
    out += "3,124,5,26,10,1,0,0,0,0,0,0,0,\n";
    out += "10,20,30,0,40,0,0,0,50,60,\n70,0,80,90,100,0,0,0,0,0,\n";
    out += "0.,0.,0.,\n0.,0.,0.,\n0.,0.,0.,\n0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,\n";
    out += "4,124,5,29,13,1,0,0,0,0,0,0,0,\n10,0,0,0,0,0,0,0,0,0,\n0,0,0,0,0,0,0,0,0,0,\n";
    out += "0.,0.,0.,\n0.,0.,0.,\n0.,0.,0.,\n0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,\n";
    out += "20,0,1.,1,1,1,0,0,0,\n30,0,1.,1,1,1,0,0,0,\n-1,\n   -1\n";
    out += "   -1\n   408\n9,0,0,\nBASE\n0,0,0,\n0,0,0,0,0.,0.,\n0,0,\n";
    for (int k = 0; k < 6; ++k)
        out += "0,0,\n0.,0.,0.,\n0.,0.,0.,\n";
    out +=
        "91,\n17,\n10,30,10,1,\n-1,-1,-1,-1,\n-1,\n23,\n7,\n10,\n20,\n30,\n-1,\n8,\n3,\n-1,\n"
        "-1,\n   -1\n";
    // Two output sets; a nodal 451 vector in set 1 and a ranged 1051 one in set 2.
    out +=
        "   -1\n   450\n1,\nLOAD 1\n0,1,\n0.5,\n1,\na note\n2,\nLOAD 2\n0,1,\n1.5,\n0,\n"
        "   -1\n";
    out +=
        "   -1\n   451\n1,1,1,\nT\n0.,1.,1.,\n0,0,0,0,0,0,0,0,0,0,\n0,0,0,0,0,0,0,0,0,0,\n"
        "10,100,1,7,\n1,1,1,\n10,1.,\n100,10.,\n-1,0.,\n   -1\n";
    out +=
        "   -1\n  1051\n2,1,1,\nT\n0.,1.,1.,\n0,0,0,0,0,0,0,0,0,0,\n0,0,0,0,0,0,0,0,0,0,\n"
        "0,\n3,3,1,8,\n1,1,1,\n3,3,42.,\n-1,0.,\n   -1\n";
    return out;
}

double coord(const Mesh& rMesh, std::int64_t P, std::size_t C) {
    return detail::read_double(rMesh.Points(), static_cast<std::size_t>(P) * 3 + C);
}

}  // namespace

TEST(Femap, ReadsATetra10FromBrickSlots) {
    const Mesh mesh = meshioplusplus::read_femap(write_file(tet10_file()));
    ASSERT_EQ(mesh.NumPoints(), 10u);
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);  // the rigid element is skipped
    ASSERT_EQ(mesh.Cells(0).Type(), "tetra10");
    const auto& c = mesh.Cells(0).Conn();
    // meshio++ tetra10: corners, then 0-1, 1-2, 2-0, 0-3, 1-3, 2-3.
    const int edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    for (int e = 0; e < 6; ++e)
        for (std::size_t d = 0; d < 3; ++d)
            EXPECT_DOUBLE_EQ(coord(mesh, detail::read_int(c, 4 + e), d),
                             0.5 * (coord(mesh, detail::read_int(c, edges[e][0]), d) +
                                    coord(mesh, detail::read_int(c, edges[e][1]), d)))
                << "edge " << e;
    EXPECT_EQ(detail::read_int(mesh.CellData("femap:property", 0), 0), 5);
    EXPECT_EQ(detail::read_int(mesh.CellData("femap:type", 0), 0), 26);
}

TEST(Femap, PropertiesAndGroupsBecomeRegions) {
    const Mesh mesh = meshioplusplus::read_femap(write_file(tet10_file()));
    const std::size_t prop = mesh.FindRegion("SOLID PART", RegionKind::Cell);
    ASSERT_NE(prop, Mesh::npos);
    EXPECT_EQ(mesh.Region(prop).mTag, 5);
    const std::size_t nodes = mesh.FindRegion("BASE", RegionKind::Point);
    const std::size_t cells = mesh.FindRegion("BASE", RegionKind::Cell);
    ASSERT_NE(nodes, Mesh::npos);
    ASSERT_NE(cells, Mesh::npos);
    EXPECT_EQ(mesh.Region(nodes).NumEntries(), 3u);
    EXPECT_EQ(mesh.Region(cells).mTag, 9);
}

TEST(Femap, OutputSetsAreSteps) {
    const std::string path = write_file(tet10_file());
    EXPECT_EQ(meshioplusplus::femap_time_values(path), (std::vector<double>{0.5, 1.5}));
    const Mesh first = meshioplusplus::read_femap(path);
    EXPECT_DOUBLE_EQ(detail::read_double(first.FieldData("meshio:time"), 0), 0.5);
    const auto& t = first.PointData("T");
    EXPECT_DOUBLE_EQ(detail::read_double(t, 0), 1.0);
    EXPECT_DOUBLE_EQ(detail::read_double(t, 9), 10.0);
    EXPECT_TRUE(std::isnan(detail::read_double(t, 4)));
    ReadOptions last;
    last.mTimeStep = -1;
    const Mesh second = meshioplusplus::read_femap(path, last);
    EXPECT_DOUBLE_EQ(detail::read_double(second.CellData("T", 0), 0), 42.0);
    ReadOptions bad;
    bad.mTimeStep = 2;
    EXPECT_THROW(meshioplusplus::read_femap(path, bad), ReadError);
    const auto meta = meshioplusplus::read_femap_metadata(path, ReadOptions{});
    EXPECT_EQ(meta.mTimeValues.size(), 2u);
}

TEST(Femap, RoundTripsThroughItsOwnWriter) {
    const Mesh mesh = meshioplusplus::read_femap(write_file(tet10_file()));
    const std::string out = mt::temp_path(".neu");
    meshioplusplus::write_femap(out, mesh);
    const Mesh back = meshioplusplus::read_femap(out);
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    for (std::size_t k = 0; k < mesh.NumPoints() * 3; ++k)
        EXPECT_DOUBLE_EQ(detail::read_double(back.Points(), k),
                         detail::read_double(mesh.Points(), k));
    const auto& x = mesh.Cells(0).Conn();
    const auto& y = back.Cells(0).Conn();
    for (std::size_t k = 0; k < x.Size(); ++k)
        EXPECT_EQ(detail::read_int(x, k), detail::read_int(y, k));
    ASSERT_EQ(back.NumRegions(), mesh.NumRegions());
    for (std::size_t r = 0; r < mesh.NumRegions(); ++r) {
        EXPECT_EQ(back.Region(r).mName, mesh.Region(r).mName);
        EXPECT_EQ(back.Region(r).mTag, mesh.Region(r).mTag);
    }
}

TEST(Femap, ErrorsNameTheCulprit) {
    auto expect_error = [](const std::string& rBody, const std::string& rNeedle) {
        try {
            meshioplusplus::read_femap(write_file(rBody));
            ADD_FAILURE() << "no error for " << rNeedle;
        } catch (const ReadError& e) {
            EXPECT_NE(std::string(e.what()).find(rNeedle), std::string::npos) << e.what();
        }
    };
    expect_error("   -1\n   100\n<NULL>\n8.2,\n   -1\n", "no nodes");
    expect_error("   -1\n   403\n1,0,0,1,46,0,0,0,0,0,0,0.,0.,x,0,\n   -1\n", "bad coordinate");
    expect_error(
        "   -1\n   403\n1,0,0,1,46,0,0,0,0,0,0,0.,0.,0.,0,\n   -1\n   -1\n   404\n"
        "7,124,1,1,0,1,0,0,\n1,9,0,0,0,0,0,0,0,0,\n0,0,0,0,0,0,0,0,0,0,\n0.,0.,0.,\n"
        "0.,0.,0.,\n0.,0.,0.,\n0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,\n   -1\n",
        "undefined node 9");
    expect_error(
        "   -1\n   403\n1,0,0,1,46,0,0,0,0,0,0,0.,0.,0.,0,\n   -1\n   -1\n   404\n"
        "7,124,1,1,0,1,0,0,\n1,2,0,0,0,0,0,0,0,0,\n   -1\n",
        "ends inside");
}
