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
// Tests for the FEBio `.feb` reader (spec 2.5/3.0/4.0) and spec-4.0 writer.

// External includes
#include <gtest/gtest.h>

// System includes
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/febio.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/region.hpp"

using namespace meshioplusplus;
using namespace mt;

namespace {

std::string put(const std::string& rText) {
    const std::string path = temp_path(".feb");
    std::ofstream(path, std::ios::binary) << rText;
    return path;
}

std::string slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

const meshioplusplus::Region* find_region(const Mesh& rMesh, const std::string& rName,
                                          RegionKind Kind) {
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
        if (rMesh.Region(r).mName == rName && rMesh.Region(r).mKind == Kind)
            return &rMesh.Region(r);
    return nullptr;
}

// A unit hex27 in FEBio's numbering (FECore FEHex27): corners, bottom, top and
// vertical mid-edges, then face centres y-, x+, y+, x-, z-, z+ and the centre.
const double kHex27[27][3] = {
    {0, 0, 0},   {1, 0, 0},   {1, 1, 0},   {0, 1, 0},   {0, 0, 1},   {1, 0, 1},    {1, 1, 1},
    {0, 1, 1},   {.5, 0, 0},  {1, .5, 0},  {.5, 1, 0},  {0, .5, 0},  {.5, 0, 1},   {1, .5, 1},
    {.5, 1, 1},  {0, .5, 1},  {0, 0, .5},  {1, 0, .5},  {1, 1, .5},  {0, 1, .5},   {.5, 0, .5},
    {1, .5, .5}, {.5, 1, .5}, {0, .5, .5}, {.5, .5, 0}, {.5, .5, 1}, {.5, .5, .5},
};

std::string hex27_spec40() {
    std::ostringstream s;
    s << "<?xml version=\"1.0\"?>\n<febio_spec version=\"4.0\">\n<Module type=\"solid\"/>\n"
      << "<Material><material id=\"7\" name=\"steel\" type=\"isotropic elastic\"/></Material>\n"
      << "<Mesh>\n<Nodes name=\"all\">\n";
    for (int i = 0; i < 27; ++i)  // sparse ids: 10, 20, ...
        s << "<node id=\"" << 10 * (i + 1) << "\">" << kHex27[i][0] << ',' << kHex27[i][1] << ','
          << kHex27[i][2] << "</node>\n";
    s << "</Nodes>\n<Elements type=\"hex27\" name=\"brick\">\n<elem id=\"5\">";
    for (int i = 0; i < 27; ++i)
        s << (i ? "," : "") << 10 * (i + 1);
    s << "</elem>\n</Elements>\n"
      << "<NodeSet name=\"bottom\">10:40:10</NodeSet>\n"
      << "<ElementSet name=\"all_elems\">5</ElementSet>\n"
      << "<Surface name=\"base\"><quad4 id=\"1\">10,40,30,20</quad4></Surface>\n"
      << "<Surface name=\"loose\"><tri3 id=\"1\">10,20,70</tri3></Surface>\n"
      << "</Mesh>\n<MeshDomains><SolidDomain name=\"brick\" mat=\"steel\"/></MeshDomains>\n"
      << "<MeshData><NodeData name=\"t\" node_set=\"bottom\" data_type=\"scalar\">"
      << "<node lid=\"2\">5</node></NodeData></MeshData>\n</febio_spec>\n";
    return s.str();
}

}  // namespace

TEST(Febio, ReadsSpec40WithHex27InFebioOrder) {
    const std::string path = put(hex27_spec40());
    EXPECT_EQ(sniff_format(path), "febio");
    const Mesh mesh = read_febio(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);  // the brick, then the loose surface
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron27");
    EXPECT_EQ(mesh.Cells(1).Type(), "triangle");
    // meshio++ (VTK) hexahedron27 face centres: x-, x+, y-, y+, z-, z+.
    const double centres[6][3] = {{0, .5, .5}, {1, .5, .5}, {.5, 0, .5},
                                  {.5, 1, .5}, {.5, .5, 0}, {.5, .5, 1}};
    const NDArray& conn = mesh.Cells(0).Conn();
    for (int f = 0; f < 6; ++f)
        for (int d = 0; d < 3; ++d)
            EXPECT_DOUBLE_EQ(detail::read_double(
                                 mesh.Points(),
                                 static_cast<std::size_t>(detail::read_int(conn, 20 + f) * 3 + d)),
                             centres[f][d])
                << "face centre " << f;
    const auto* brick = find_region(mesh, "brick", RegionKind::Cell);
    ASSERT_NE(brick, nullptr);
    EXPECT_EQ(brick->mTag, 7);  // the domain's material id
    EXPECT_EQ(find_region(mesh, "bottom", RegionKind::Point)->NumEntries(), 4u);
    const auto* base = find_region(mesh, "base", RegionKind::Side);
    ASSERT_NE(base, nullptr);
    EXPECT_EQ(base->Entries()[1], 4);  // hexahedron face 4 = z-
    EXPECT_NE(find_region(mesh, "loose", RegionKind::Cell), nullptr);
    const NDArray& t = mesh.PointData("t");
    EXPECT_EQ(detail::read_double(t, 1), 5.0);  // lid 2 of `bottom` = node id 20
    EXPECT_TRUE(std::isnan(detail::read_double(t, 0)));
}

TEST(Febio, WritesAndReadsBack) {
    const Mesh mesh = read_febio(put(hex27_spec40()));
    const std::string out = temp_path(".feb");
    write_febio(out, mesh);
    const std::string text = slurp(out);
    EXPECT_NE(text.find("<febio_spec version=\"4.0\">"), std::string::npos);
    // `all_elems` and `brick` both cover the block exactly: the first in name
    // order names it, the other is written as an <ElementSet>.
    EXPECT_NE(text.find("<Elements type=\"hex27\" name=\"all_elems\">"), std::string::npos);
    EXPECT_NE(text.find("<SolidDomain name=\"all_elems\" mat=\"all_elems\"/>"),
              std::string::npos);
    EXPECT_NE(text.find("<ElementSet name=\"brick\">1</ElementSet>"), std::string::npos);
    // The loose triangle is not on a solid face: a shell, not a surface.
    EXPECT_NE(text.find("<ShellDomain name=\"loose\" mat=\"loose\"/>"), std::string::npos);
    const Mesh back = read_febio(out);
    ASSERT_EQ(back.NumCellBlocks(), 2u);
    EXPECT_EQ(cell_rows(back), cell_rows(mesh));
    expect_points_close(mesh, back, 0.0);
    for (const char* name : {"brick", "loose", "all_elems"})
        EXPECT_NE(find_region(back, name, RegionKind::Cell), nullptr) << name;
    EXPECT_NE(find_region(back, "base", RegionKind::Side), nullptr);
    EXPECT_NE(find_region(back, "bottom", RegionKind::Point), nullptr);
}

TEST(Febio, ReadsSpec25Geometry) {
    const std::string path =
        put("<febio_spec version=\"2.5\"><Module type=\"solid\"/><Geometry>"
            "<Nodes><node id=\"1\">0,0,0</node><node id=\"2\">1,0,0</node>"
            "<node id=\"3\">0,1,0</node><node id=\"4\">0,0,1</node></Nodes>"
            "<Elements type=\"tet4\" mat=\"3\" name=\"t\"><elem id=\"1\">1,2,3,4</elem></Elements>"
            "<NodeSet name=\"n\"><node id=\"1\"/><node id=\"4\"/></NodeSet>"
            "<ElementSet name=\"e\"><elem id=\"1\"/></ElementSet>"
            "</Geometry><MeshData><ElementData var=\"h\" elem_set=\"e\"><elem lid=\"1\">2</elem>"
            "</ElementData></MeshData></febio_spec>");
    const Mesh mesh = read_febio(path);
    EXPECT_EQ(find_region(mesh, "t", RegionKind::Cell)->mTag, 3);
    EXPECT_EQ(find_region(mesh, "n", RegionKind::Point)->NumEntries(), 2u);
    EXPECT_EQ(detail::read_double(mesh.CellData("h", 0), 0), 2.0);
}

TEST(Febio, Errors) {
    EXPECT_THROW(read_febio(put("<febio_spec version=\"1.2\"/>")), ReadError);
    EXPECT_THROW(read_febio(put("<notfebio/>")), ReadError);
    const std::string tet5 =
        put("<febio_spec version=\"4.0\"><Mesh><Nodes>"
            "<node id=\"1\">0,0,0</node><node id=\"2\">1,0,0</node><node id=\"3\">0,1,0</node>"
            "<node id=\"4\">0,0,1</node><node id=\"5\">.2,.2,.2</node></Nodes>"
            "<Elements type=\"tet5\" name=\"p\"><elem id=\"1\">1,2,3,4,5</elem></Elements>"
            "</Mesh></febio_spec>");
    EXPECT_THROW(read_febio(tet5), ReadError);
    ReadOptions lenient;
    lenient.mLenient = true;
    EXPECT_EQ(read_febio(tet5, lenient).Cells(0).Type(), "tetra");
    EXPECT_THROW(read_febio(put("<febio_spec version=\"4.0\"><Mesh><Elements type=\"tet4\" "
                                "name=\"p\"><elem id=\"1\">1,2,3,4</elem></Elements></Mesh>"
                                "</febio_spec>")),
                 ReadError);
    Mesh polygon;
    polygon.AssignPoints(points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 2, 0}}));
    polygon.AddCellBlock("polygon5", conn_from({{0, 1, 2, 3, 4}}));
    EXPECT_THROW(write_febio(temp_path(".feb"), polygon), WriteError);
}
