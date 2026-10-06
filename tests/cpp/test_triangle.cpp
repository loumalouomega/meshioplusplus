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
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/triangle.hpp"

namespace {

// Triangle writes sibling files; wrap mt::roundtrip so the siblings are
// removed too.
void rt(const mt::Mesh& mesh, const std::string& rSuffix) {
    std::string stem;
    mt::roundtrip(
        [&](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_triangle(p, m);
            stem = p.substr(0, p.size() - rSuffix.size());
        },
        [](const std::string& p) { return meshioplusplus::read_triangle(p); }, mesh, rSuffix);
    std::error_code ec;
    std::filesystem::remove(stem + ".node", ec);
    std::filesystem::remove(stem + ".ele", ec);
    std::filesystem::remove(stem + ".poly", ec);
}

mt::Mesh triangle6_2d_mesh() {
    return mt::make_mesh({{0, 0}, {1, 0}, {1, 1}, {0.5, 0}, {1, 0.5}, {0.5, 0.5}}, "triangle6",
                         {{0, 1, 2, 3, 4, 5}});
}

mt::Mesh line_mesh_2d() {
    return mt::make_mesh({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, "line",
                         {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
}

}  // namespace

TEST(Triangle, NodeElePair) {
    rt(mt::tri_mesh_2d(), ".node");
    rt(mt::tri_mesh_2d(), ".ele");
}
TEST(Triangle, Triangle6) {
    rt(triangle6_2d_mesh(), ".node");
}
TEST(Triangle, Poly) {
    rt(line_mesh_2d(), ".poly");
}
TEST(Triangle, Rejects3DPoints) {
    EXPECT_THROW(meshioplusplus::write_triangle(mt::temp_path(".node"), mt::tri_mesh()),
                 meshioplusplus::WriteError);
}

// ---- the text reader's tokenizer, pinned (roadmap 3.1.1.1) -----------------
//
// What the reader accepts, refuses and returns for hand-made text, so a change
// of its tokenizer can be checked against the reader it replaces. The
// odd-number cases record what the reader did when they were written.

namespace {

// Sibling `.node` / `.ele` / `.poly` files on disk, removed on scope exit.
class TriangleFiles {
public:
    TriangleFiles() : mStem(mt::temp_path("_tri")) {}
    ~TriangleFiles() {
        std::error_code ec;
        for (const char* ext : {".node", ".ele", ".poly"})
            std::filesystem::remove(mStem + ext, ec);
    }
    TriangleFiles& Node(const std::string& rText) { return Put(".node", rText); }
    TriangleFiles& Ele(const std::string& rText) { return Put(".ele", rText); }
    TriangleFiles& Poly(const std::string& rText) { return Put(".poly", rText); }
    std::string Path(const char* pExt) const { return mStem + pExt; }

private:
    TriangleFiles& Put(const char* pExt, const std::string& rText) {
        std::ofstream os(mStem + pExt, std::ios::binary);
        os << rText;
        return *this;
    }
    std::string mStem;
};

// The message of the `ReadError` reading @p pExt of @p rFiles raises, or "".
std::string tri_error(const TriangleFiles& rFiles, const char* pExt = ".node") {
    try {
        meshioplusplus::read_triangle(rFiles.Path(pExt));
    } catch (const meshioplusplus::ReadError& rExc) {
        return rExc.what();
    }
    return "";
}

const char kTriNode[] =
    "# vertices\n"
    "\n"
    "   # indented\n"
    "3 2 1 1\n"
    "1 0 0 0.5 7\n"
    "2 1 0 1.5 8\n"
    "3 0 1 2.5 9\n";
const char kTriEle[] =
    "# triangles\n"
    "1 3 1\n"
    "1 1 2 3 42\n";

std::string tri_replace_all(std::string Text, const std::string& rFrom, const std::string& rTo) {
    for (std::size_t pos = 0; (pos = Text.find(rFrom, pos)) != std::string::npos;
         pos += rTo.size())
        Text.replace(pos, rFrom.size(), rTo);
    return Text;
}

mt::Mesh tri_read(const std::string& rNode, const std::string& rEle) {
    TriangleFiles files;
    files.Node(rNode).Ele(rEle);
    return meshioplusplus::read_triangle(files.Path(".node"));
}

void expect_unit_triangle(const mt::Mesh& rMesh) {
    ASSERT_EQ(rMesh.NumPoints(), 3u);
    ASSERT_EQ(rMesh.PointDim(), 2u);
    const double* p = rMesh.Points().As<double>();
    const double expected[3][2] = {{0, 0}, {1, 0}, {0, 1}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t d = 0; d < 2; ++d)
            EXPECT_EQ(p[2 * i + d], expected[i][d]);
    const double* attr = rMesh.PointData("triangle:attr1").As<double>();
    const double* ref = rMesh.PointData("triangle:ref").As<double>();
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_EQ(attr[i], 0.5 + static_cast<double>(i));
        EXPECT_EQ(ref[i], 7.0 + static_cast<double>(i));
    }
    ASSERT_EQ(rMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(rMesh.Cells(0).Type(), "triangle");
    ASSERT_EQ(rMesh.Cells(0).NumCells(), 1u);
    const std::int64_t* c = rMesh.Cells(0).Conn().As<std::int64_t>();
    for (std::int64_t k = 0; k < 3; ++k)
        EXPECT_EQ(c[k], k);
    ASSERT_TRUE(rMesh.HasCellData("triangle:ref"));
    EXPECT_EQ(rMesh.CellData("triangle:ref", 0).As<double>()[0], 42.0);
}

// A zero-triangle file, so a test can read just its points.
mt::Mesh tri_points(const std::string& rBody, std::size_t Count) {
    return tri_read(std::to_string(Count) + " 2 0 0\n" + rBody, "0 3 0\n");
}

}  // namespace

TEST(TriangleText, ReadsAnAttributedTriangle) {
    expect_unit_triangle(tri_read(kTriNode, kTriEle));
}

TEST(TriangleText, AnEleFileNamesTheSamePair) {
    TriangleFiles files;
    files.Node(kTriNode).Ele(kTriEle);
    expect_unit_triangle(meshioplusplus::read_triangle(files.Path(".ele")));
}

TEST(TriangleText, CrlfLineEndings) {
    expect_unit_triangle(
        tri_read(tri_replace_all(kTriNode, "\n", "\r\n"), tri_replace_all(kTriEle, "\n", "\r\n")));
}

TEST(TriangleText, TabsFormFeedsAndVerticalTabsSeparate) {
    std::string node = tri_replace_all(kTriNode, " ", "\t");
    node = tri_replace_all(node, "1\t0\t0\t0.5", "\f1\v0\f0\v0.5");
    expect_unit_triangle(tri_read(node, tri_replace_all(kTriEle, " ", "\t")));
}

TEST(TriangleText, TokensRunAcrossLines) {
    std::string node;
    for (const char* tok : {"3", "2", "1", "1", "1", "0", "0", "0.5", "7", "2", "1", "0", "1.5",
                            "8", "3", "0", "1", "2.5", "9"})
        node += std::string(tok) + "\n";
    expect_unit_triangle(tri_read(node, "1\n3\n1\n1\n1\n2\n3\n42\n"));
}

TEST(TriangleText, ACommentRunsFromAnyHashToTheEndOfTheLine) {
    // Unlike TetGen's, an inline # starts a comment, wherever it is.
    expect_unit_triangle(
        tri_read(tri_replace_all(tri_replace_all(kTriNode, "3 2 1 1\n", "3 2 1 1 # header\n"),
                                 "2 1 0 1.5 8\n", "2 1 0 1.5 8# a vertex\n"),
                 tri_replace_all(kTriEle, "1 1 2 3 42\n", "1 1 2 3 42 #### done\n")));
    expect_unit_triangle(
        tri_read(tri_replace_all(kTriNode, "2 1 0 1.5 8\n", "2 1 0 # cut\n1.5 8\n"), kTriEle));
}

TEST(TriangleText, AFileOfOnlyCommentsEndsEarly) {
    TriangleFiles files;
    files.Node("# nothing\n\n   # at all\n");
    EXPECT_EQ(tri_error(files), "Triangle: unexpected end of file reading vertex count");
}

TEST(TriangleText, ALoneNodeFileIsAPointCloud) {
    TriangleFiles files;
    files.Node(kTriNode);
    const mt::Mesh m = meshioplusplus::read_triangle(files.Path(".node"));
    EXPECT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(m.NumCellBlocks(), 0u);
}

TEST(TriangleText, NumbersAreParsedLeniently) {
    const mt::Mesh m = tri_points("1 1e1 -2.5E-1\n"
                                  "2 1.5abc .5\n"
                                  "+3 -.5e1 12xyz\n",
                                  3);
    const double* p = m.Points().As<double>();
    const double expected[3][2] = {{10, -0.25}, {1.5, 0.5}, {-5, 12}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t d = 0; d < 2; ++d)
            EXPECT_EQ(p[2 * i + d], expected[i][d]) << i << ',' << d;
}

TEST(TriangleText, ANumberNeedsADigitToStart) {
    TriangleFiles files;
    files.Node("1 2 0 0\nabc 0 0\n");
    EXPECT_EQ(tri_error(files), "Triangle: expected an integer for vertex index");
    files.Node("1 2 0 0\n1 abc 0\n");
    EXPECT_EQ(tri_error(files), "Triangle: expected a number for x coordinate");
    files.Node("1 2 0 0\n1 0 -\n");
    EXPECT_EQ(tri_error(files), "Triangle: expected a number for y coordinate");
    files.Node("1 2 0 0\n+ 0 0\n");
    EXPECT_EQ(tri_error(files), "Triangle: expected an integer for vertex index");
}

TEST(TriangleText, SpecialDoubleTokens) {
    const mt::Mesh m = tri_points("1 nan inf\n"
                                  "2 -inf 1e999\n"
                                  "3 0x10 infinity\n"
                                  "4 -nan 1e-999\n",
                                  4);
    const double* p = m.Points().As<double>();
    // Recorded from the reader as it was: nan and inf words, an overflow to
    // infinity, a hexadecimal float, and an underflow to zero.
    EXPECT_TRUE(std::isnan(p[0]));
    EXPECT_EQ(p[1], std::numeric_limits<double>::infinity());
    EXPECT_EQ(p[2], -std::numeric_limits<double>::infinity());
    EXPECT_EQ(p[3], std::numeric_limits<double>::infinity());
    EXPECT_EQ(p[4], 16.0);
    EXPECT_EQ(p[5], std::numeric_limits<double>::infinity());
    EXPECT_TRUE(std::isnan(p[6]));
    EXPECT_TRUE(std::signbit(p[6]));
    EXPECT_EQ(p[7], 0.0);
}

TEST(TriangleText, ZeroBasedNumberingShiftsTheConnectivity) {
    const mt::Mesh m = tri_read("3 2 0 0\n0 0 0\n1 1 0\n2 0 1\n", "1 3 0\n0 0 1 2\n");
    ASSERT_EQ(m.Cells(0).NumCells(), 1u);
    const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
    for (std::int64_t k = 0; k < 3; ++k)
        EXPECT_EQ(c[k], k);
}

TEST(TriangleText, SixNodeTrianglesMakeATriangle6Block) {
    const mt::Mesh m = tri_read("6 2 0 0\n1 0 0\n2 1 0\n3 1 1\n4 .5 0\n5 1 .5\n6 .5 .5\n",
                                "1 6 0\n1 1 2 3 4 5 6\n");
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "triangle6");
}

TEST(TriangleText, RefusedInputsNameTheirReason) {
    TriangleFiles files;
    files.Node("3 3 0 0\n").Ele(kTriEle);
    EXPECT_EQ(tri_error(files), "Triangle: need 2D points");
    files.Node("-1 2 0 0\n");
    EXPECT_EQ(tri_error(files), "Triangle: malformed vertex header");
    files.Node("1 2 -1 0\n");
    EXPECT_EQ(tri_error(files), "Triangle: malformed vertex header");
    files.Node(tri_replace_all(kTriNode, "3 0 1 2.5 9\n", "5 0 1 2.5 9\n"));
    EXPECT_EQ(tri_error(files), "Triangle: vertices not numbered consecutively");
    files.Node(tri_replace_all(kTriNode, "3 0 1 2.5 9\n", "3 0 1 2.5\n"));
    EXPECT_EQ(tri_error(files), "Triangle: unexpected end of file reading boundary marker");

    files.Node(kTriNode).Ele(tri_replace_all(kTriEle, "1 3 1\n", "1 4 1\n"));
    EXPECT_EQ(tri_error(files), "Triangle: only 3- or 6-node triangles are supported");
    files.Ele(tri_replace_all(kTriEle, "1 1 2 3 42\n", "1 1 2 9 42\n"));
    EXPECT_EQ(tri_error(files), "Triangle: connectivity index out of range");
    files.Ele(tri_replace_all(kTriEle, "1 1 2 3 42\n", "1 1 2 99999999999999999999 42\n"));
    EXPECT_EQ(tri_error(files), "Triangle: connectivity index out of range");
    files.Ele(tri_replace_all(kTriEle, "1 1 2 3 42\n", "1 1 2 3\n"));
    EXPECT_EQ(tri_error(files), "Triangle: unexpected end of file reading triangle attribute");
}

TEST(TriangleText, TheReadErrorsOfAnUnreadableFile) {
    TriangleFiles files;
    files.Ele(kTriEle);
    EXPECT_NE(tri_error(files).find("Triangle: could not open file: "), std::string::npos);
    try {
        meshioplusplus::read_triangle(mt::temp_path(".vtk"));
        FAIL() << "expected a ReadError";
    } catch (const meshioplusplus::ReadError& rExc) {
        EXPECT_STREQ(rExc.what(), "Triangle: expected a .node, .ele, or .poly file");
    }
}

TEST(TriangleText, APolyFileHoldsInlineVerticesAndSegments) {
    TriangleFiles files;
    files.Poly(
        "# a square\n"
        "4 2 0 1\n"
        "1 0 0 5\n2 1 0 5\n3 1 1 5\n4 0 1 5\n"
        "4 1\n"
        "1 1 2 11\n2 2 3 12\n3 3 4 13\n4 4 1 14 # last\n"
        "0\n");
    const mt::Mesh m = meshioplusplus::read_triangle(files.Path(".poly"));
    EXPECT_EQ(m.NumPoints(), 4u);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "line");
    EXPECT_EQ(m.Cells(0).NumCells(), 4u);
    ASSERT_TRUE(m.HasCellData("triangle:ref"));
    EXPECT_EQ(m.CellData("triangle:ref", 0).As<std::int64_t>()[3], 14);
    EXPECT_EQ(m.PointData("triangle:ref").As<double>()[0], 5.0);
}

TEST(TriangleText, APolyFileMayTakeItsVerticesFromTheSiblingNodeFile) {
    TriangleFiles files;
    files.Node("3 2 0 0\n1 0 0\n2 1 0\n3 0 1\n");
    files.Poly("0 2 0 0\n2 0\n1 1 2\n2 2 3\n");
    const mt::Mesh m = meshioplusplus::read_triangle(files.Path(".poly"));
    EXPECT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(m.Cells(0).NumCells(), 2u);
}

TEST(TriangleText, APolyFileSkipsHolesAndRegions) {
    TriangleFiles files;
    files.Poly("3 2 0 0\n1 0 0\n2 1 0\n3 0 1\n"
               "1 0\n1 1 2\n"
               "1\n1 0.2 0.2\n"
               "1\n1 0.5 0.5 3 0.1\n");
    const mt::Mesh m = meshioplusplus::read_triangle(files.Path(".poly"));
    EXPECT_EQ(m.Cells(0).NumCells(), 1u);
}

TEST(TriangleText, APolyFileRefusesWhatItCannotRead) {
    TriangleFiles files;
    files.Poly("0 2 0 0\n");
    EXPECT_EQ(tri_error(files, ".poly"), "Triangle: .poly refers to a missing sibling .node file");
    const std::string verts = "3 2 0 0\n1 0 0\n2 1 0\n3 0 1\n";
    files.Poly(verts + "1 2\n1 1 2\n");
    EXPECT_EQ(tri_error(files, ".poly"), "Triangle: malformed segment header");
    files.Poly(verts + "1 0\n1 1 9\n");
    EXPECT_EQ(tri_error(files, ".poly"), "Triangle: segment endpoint out of range");
    files.Poly(verts + "2 0\n1 1 2\n");
    EXPECT_EQ(tri_error(files, ".poly"),
              "Triangle: unexpected end of file reading segment index");
}
