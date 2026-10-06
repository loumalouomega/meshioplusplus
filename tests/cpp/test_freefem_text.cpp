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
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/freefem.hpp"

// Pins what the FreeFEM `.msh` reader accepts, refuses and returns for
// hand-made text, so a change of its tokenizer (roadmap 3.1.1.1) can be checked
// against the reader it replaces. The odd-number cases record what the reader
// did when this file was written. Values are read through the dtype-neutral
// helpers, because a mesh backend may store points as float32 or float64.

namespace {

using meshioplusplus::detail::read_double;
using meshioplusplus::detail::read_int;

// A `.msh` file on disk, removed on scope exit.
class FemFile {
public:
    explicit FemFile(const std::string& rText) : mPath(mt::temp_path(".msh")) {
        std::ofstream os(mPath, std::ios::binary);
        os << rText;
    }
    ~FemFile() {
        std::error_code ec;
        std::filesystem::remove(mPath, ec);
    }
    const std::string& Path() const { return mPath; }

private:
    std::string mPath;
};

mt::Mesh fem_read(const std::string& rText) {
    const FemFile file(rText);
    return meshioplusplus::read_freefem(file.Path());
}

// The message of the `ReadError` a text raises, or "" when it reads.
std::string fem_error(const std::string& rText) {
    const FemFile file(rText);
    try {
        meshioplusplus::read_freefem(file.Path());
    } catch (const meshioplusplus::ReadError& rExc) {
        return rExc.what();
    }
    return "";
}

std::string fem_replace_all(std::string Text, const std::string& rFrom, const std::string& rTo) {
    for (std::size_t pos = 0; (pos = Text.find(rFrom, pos)) != std::string::npos;
         pos += rTo.size())
        Text.replace(pos, rFrom.size(), rTo);
    return Text;
}

// A right triangle with three boundary edges, vertex references 10..30, one
// triangle region (7) and edge references 1..3.
const char kFem2d[] =
    "3 1 3\n"
    "0 0 10\n"
    "1 0 20\n"
    "0 1 30\n"
    "1 2 3 7\n"
    "1 2 1\n"
    "2 3 2\n"
    "3 1 3\n";

// A unit tetrahedron with one boundary triangle.
const char kFem3d[] =
    "4 1 1\n"
    "0 0 0 1\n"
    "1 0 0 2\n"
    "0 1 0 3\n"
    "0 0 1 4\n"
    "1 2 3 4 9\n"
    "1 2 3 5\n";

std::int64_t cell_at(const mt::Mesh& rMesh, std::size_t Block, std::size_t Index) {
    return read_int(rMesh.Cells(Block).Conn(), Index);
}

void expect_2d(const mt::Mesh& rMesh) {
    ASSERT_EQ(rMesh.PointDim(), 2u);
    ASSERT_EQ(rMesh.NumPoints(), 3u);
    const double expected[3][2] = {{0, 0}, {1, 0}, {0, 1}};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t d = 0; d < 2; ++d)
            EXPECT_EQ(read_double(rMesh.Points(), 2 * i + d), expected[i][d]);
        EXPECT_EQ(read_int(rMesh.PointData("freefem:ref"), i), 10 * static_cast<std::int64_t>(i + 1));
    }
    ASSERT_EQ(rMesh.NumCellBlocks(), 2u);
    EXPECT_EQ(rMesh.Cells(0).Type(), "triangle");
    ASSERT_EQ(rMesh.Cells(0).NumCells(), 1u);
    for (std::size_t k = 0; k < 3; ++k)
        EXPECT_EQ(cell_at(rMesh, 0, k), static_cast<std::int64_t>(k));
    EXPECT_EQ(rMesh.Cells(1).Type(), "line");
    ASSERT_EQ(rMesh.Cells(1).NumCells(), 3u);
    const std::int64_t edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (std::size_t e = 0; e < 3; ++e)
        for (std::size_t k = 0; k < 2; ++k)
            EXPECT_EQ(cell_at(rMesh, 1, 2 * e + k), edges[e][k]);
    ASSERT_TRUE(rMesh.HasCellData("freefem:ref"));
    EXPECT_EQ(read_int(rMesh.CellData("freefem:ref", 0), 0), 7);
    for (std::size_t e = 0; e < 3; ++e)
        EXPECT_EQ(read_int(rMesh.CellData("freefem:ref", 1), e), static_cast<std::int64_t>(e + 1));
}

// A one-triangle-free file of `Count` 2-D vertices, to read just points.
mt::Mesh fem_points(const std::string& rBody, std::size_t Count) {
    return fem_read(std::to_string(Count) + " 0 0\n" + rBody);
}

}  // namespace

TEST(FreeFemText, ReadsATriangleWithItsBoundary) {
    expect_2d(fem_read(kFem2d));
}

TEST(FreeFemText, ReadsATetrahedronWithItsBoundaryTriangle) {
    const mt::Mesh m = fem_read(kFem3d);
    ASSERT_EQ(m.PointDim(), 3u);
    ASSERT_EQ(m.NumPoints(), 4u);
    EXPECT_EQ(read_double(m.Points(), 3 * 3 + 2), 1.0);
    EXPECT_EQ(read_int(m.PointData("freefem:ref"), 3), 4);
    ASSERT_EQ(m.NumCellBlocks(), 2u);
    EXPECT_EQ(m.Cells(0).Type(), "tetra");
    EXPECT_EQ(m.Cells(1).Type(), "triangle");
    for (std::size_t k = 0; k < 4; ++k)
        EXPECT_EQ(cell_at(m, 0, k), static_cast<std::int64_t>(k));
    EXPECT_EQ(read_int(m.CellData("freefem:ref", 0), 0), 9);
    EXPECT_EQ(read_int(m.CellData("freefem:ref", 1), 0), 5);
}

TEST(FreeFemText, CrlfLineEndings) {
    expect_2d(fem_read(fem_replace_all(kFem2d, "\n", "\r\n")));
}

TEST(FreeFemText, TabsFormFeedsAndVerticalTabsSeparateAndIndent) {
    std::string text = fem_replace_all(kFem2d, " ", "\t");
    text = fem_replace_all(text, "1\t0\t20", "\f1\v0\f20");
    expect_2d(fem_read(text));
}

TEST(FreeFemText, BlankAndWhitespaceOnlyLinesAreSkippedEverywhere) {
    std::string text = "\n  \t \n";                  // before the header
    text += fem_replace_all(kFem2d, "\n", "\n\n \n");  // between every row
    text += "\n\t\n";                                // after the last
    expect_2d(fem_read(text));
}

TEST(FreeFemText, ExtraTokensOnARowAreIgnored) {
    expect_2d(fem_read(fem_replace_all(fem_replace_all(kFem2d, "0 1 30\n", "0 1 30 99 98\n"),
                                       "1 2 3 7\n", "1 2 3 7 6 5\n")));
}

TEST(FreeFemText, ARowCannotBeSplitAcrossLines) {
    // One row per line: the second vertex's reference on its own line is a
    // row of one token.
    EXPECT_EQ(fem_error("2 0 0\n0 0 1\n1 0\n2\n"),
              "meshio++: FreeFem: expected 3 fields on a line, got 2");
    EXPECT_EQ(fem_error("2 0 0\n0 0 1\n1\n0 2\n"),
              "meshio++: FreeFem: expected 3 fields on a line, got 1");
}

TEST(FreeFemText, NumbersAreParsedLeniently) {
    const mt::Mesh m = fem_points("1e1 -2.5E-1 1\n"
                                  "1.5abc .5 2xyz\n"
                                  "+3 -.5e1 -4\n",
                                  3);
    const double expected[3][2] = {{10, -0.25}, {1.5, 0.5}, {3, -5}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t d = 0; d < 2; ++d)
            EXPECT_EQ(read_double(m.Points(), 2 * i + d), expected[i][d]) << i << ',' << d;
    EXPECT_EQ(read_int(m.PointData("freefem:ref"), 0), 1);
    EXPECT_EQ(read_int(m.PointData("freefem:ref"), 1), 2);
    EXPECT_EQ(read_int(m.PointData("freefem:ref"), 2), -4);
}

TEST(FreeFemText, ATokenWithoutANumberReadsAsZero) {
    const mt::Mesh m = fem_points("abc - x\n", 1);
    EXPECT_EQ(read_double(m.Points(), 0), 0.0);
    EXPECT_EQ(read_double(m.Points(), 1), 0.0);
    EXPECT_EQ(read_int(m.PointData("freefem:ref"), 0), 0);
}

TEST(FreeFemText, SpecialDoubleTokens) {
    const mt::Mesh m = fem_points("nan inf 1\n"
                                  "-inf 1e999 2\n"
                                  "0x10 infinity 3\n"
                                  "-nan 1e-999 4\n",
                                  4);
    // Recorded from the reader as it was: nan and inf words, an overflow to
    // infinity, a hexadecimal float, and an underflow to zero.
    const auto p = [&](std::size_t i) { return read_double(m.Points(), i); };
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_TRUE(std::isnan(p(0)));
    EXPECT_EQ(p(1), inf);
    EXPECT_EQ(p(2), -inf);
    EXPECT_EQ(p(3), inf);
    EXPECT_EQ(p(4), 16.0);
    EXPECT_EQ(p(5), inf);
    EXPECT_TRUE(std::isnan(p(6)));
    EXPECT_TRUE(std::signbit(p(6)));
    EXPECT_EQ(p(7), 0.0);
}

TEST(FreeFemText, AnOverflowingIntegerSaturates) {
    const mt::Mesh m = fem_read("3 1 0\n0 0 99999999999999999999\n1 0 -99999999999999999999\n"
                                "0 1 0\n1 2 99999999999999999999 99999999999999999999\n");
    // References saturate like strtoll; a node id saturates and then loses
    // one to the 1-based shift, and nothing checks it against the points.
    constexpr std::int64_t kMax = std::numeric_limits<std::int64_t>::max();
    EXPECT_EQ(read_int(m.PointData("freefem:ref"), 0), kMax);
    EXPECT_EQ(read_int(m.PointData("freefem:ref"), 1), std::numeric_limits<std::int64_t>::min());
    EXPECT_EQ(cell_at(m, 0, 2), kMax - 1);
    EXPECT_EQ(read_int(m.CellData("freefem:ref", 0), 0), kMax);
}

TEST(FreeFemText, NegativeElementCountsAreNoBlocks) {
    const mt::Mesh m = fem_read("2 -1 -2\n0 0 1\n1 0 2\n");
    EXPECT_EQ(m.NumPoints(), 2u);
    EXPECT_EQ(m.NumCellBlocks(), 0u);
    EXPECT_FALSE(m.HasCellData("freefem:ref"));
}

TEST(FreeFemText, ABlockWithNoElementsIsOmitted) {
    const mt::Mesh m = fem_read("3 1 0\n0 0 1\n1 0 2\n0 1 3\n1 2 3 7\n");
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "triangle");
}

TEST(FreeFemText, RefusedInputsNameTheirReason) {
    EXPECT_EQ(fem_error(""), "FreeFem: expected a 3-integer header");
    EXPECT_EQ(fem_error("\n \n"), "FreeFem: expected a 3-integer header");
    EXPECT_EQ(fem_error("3 1\n"), "FreeFem: expected a 3-integer header");
    EXPECT_EQ(fem_error("3 1 3 0\n"), "FreeFem: expected a 3-integer header");
    EXPECT_EQ(fem_error("0 0 0\n"), "FreeFem: missing vertices");
    EXPECT_EQ(fem_error("1 0 0\n5\n"), "FreeFem: bad vertex dimension");
    EXPECT_EQ(fem_error("1 0 0\n0 0 0 0 1\n"), "FreeFem: bad vertex dimension");
    EXPECT_EQ(fem_error("3 0 0\n0 0 1\n1 0 2\n"), "FreeFem: truncated vertices");
    EXPECT_EQ(fem_error("3 1 0\n0 0 1\n1 0 2\n0 1 3\n"), "FreeFem: truncated elements");
    EXPECT_EQ(fem_error("3 1 0\n0 0 1\n1 0 2\n0 1 3\n1 2 3\n"),
              "meshio++: FreeFem: expected 4 fields on a line, got 3");
}

TEST(FreeFemText, CountsAreBoundedByTheFile) {
    EXPECT_EQ(fem_error("-1 0 0\n0 0 1\n"), "meshio++: FreeFem: negative vertex count -1");
    const std::string big = fem_error("999999 0 0\n0 0 1\n");
    EXPECT_NE(big.find("meshio++: FreeFem: vertex count 999999 exceeds what the file holds ("),
              std::string::npos)
        << big;
    const std::string elems = fem_error("1 999999 0\n0 0 1\n");
    EXPECT_NE(elems.find("meshio++: FreeFem: element count 999999 exceeds what the file holds ("),
              std::string::npos)
        << elems;
}

TEST(FreeFemText, AMissingFileCannotBeOpened) {
    try {
        meshioplusplus::read_freefem(mt::temp_path(".msh"));
        FAIL() << "expected a ReadError";
    } catch (const meshioplusplus::ReadError& rExc) {
        EXPECT_NE(std::string(rExc.what()).find("Could not open file: "), std::string::npos);
    }
}
