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
#include "meshioplusplus/formats/ugrid.hpp"

// Pins what the ASCII UGRID reader accepts, refuses and returns for hand-made
// text, so a change of its tokenizer (roadmap 3.1.1.1) can be checked against
// the reader it replaces. An ASCII file is a stream of whitespace-separated
// tokens with no line structure: seven counts, the points, the surface
// connectivity (triangles then quads), their boundary tags, then the volume
// cells. The odd-number cases record what the reader did when this file was
// written. Values are read through the dtype-neutral helpers, because a mesh
// backend may store points as float32 or float64.

namespace {

using meshioplusplus::detail::read_double;
using meshioplusplus::detail::read_int;

// A `.ugrid` file on disk (no binary key in the name: ASCII), removed on exit.
class UgridFile {
public:
    explicit UgridFile(const std::string& rText) : mPath(mt::temp_path(".ugrid")) {
        std::ofstream os(mPath, std::ios::binary);
        os << rText;
    }
    ~UgridFile() {
        std::error_code ec;
        std::filesystem::remove(mPath, ec);
    }
    const std::string& Path() const { return mPath; }

private:
    std::string mPath;
};

mt::Mesh ug_read(const std::string& rText) {
    const UgridFile file(rText);
    return meshioplusplus::read_ugrid(file.Path());
}

std::string ug_error(const std::string& rText) {
    const UgridFile file(rText);
    try {
        meshioplusplus::read_ugrid(file.Path());
    } catch (const meshioplusplus::ReadError& rExc) {
        return rExc.what();
    }
    return "";
}

std::string ug_replace_all(std::string Text, const std::string& rFrom, const std::string& rTo) {
    for (std::size_t pos = 0; (pos = Text.find(rFrom, pos)) != std::string::npos;
         pos += rTo.size())
        Text.replace(pos, rFrom.size(), rTo);
    return Text;
}

// One tetrahedron with one boundary triangle tagged 5.
const char kUgrid[] =
    "4 1 0 1 0 0 0\n"
    "0 0 0\n"
    "1 0 0\n"
    "0 1 0\n"
    "0 0 1\n"
    "1 2 3\n"
    "5\n"
    "1 2 3 4\n";

std::int64_t ug_cell(const mt::Mesh& rMesh, std::size_t Block, std::size_t Index) {
    return read_int(rMesh.Cells(Block).Conn(), Index);
}

void expect_tet(const mt::Mesh& rMesh) {
    ASSERT_EQ(rMesh.NumPoints(), 4u);
    const double expected[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t d = 0; d < 3; ++d)
            EXPECT_EQ(read_double(rMesh.Points(), 3 * i + d), expected[i][d]);
    ASSERT_EQ(rMesh.NumCellBlocks(), 2u);
    EXPECT_EQ(rMesh.Cells(0).Type(), "triangle");
    EXPECT_EQ(rMesh.Cells(1).Type(), "tetra");
    for (std::size_t k = 0; k < 3; ++k)
        EXPECT_EQ(ug_cell(rMesh, 0, k), static_cast<std::int64_t>(k));
    for (std::size_t k = 0; k < 4; ++k)
        EXPECT_EQ(ug_cell(rMesh, 1, k), static_cast<std::int64_t>(k));
    ASSERT_TRUE(rMesh.HasCellData("ugrid:ref"));
    EXPECT_EQ(read_int(rMesh.CellData("ugrid:ref", 0), 0), 5);
    EXPECT_EQ(read_int(rMesh.CellData("ugrid:ref", 1), 0), 0);
}

// `Count` points and no cells, to read just coordinates.
mt::Mesh ug_points(const std::string& rBody, std::size_t Count) {
    return ug_read(std::to_string(Count) + " 0 0 0 0 0 0\n" + rBody);
}

}  // namespace

TEST(UgridText, ReadsATetrahedronWithABoundaryTriangle) {
    expect_tet(ug_read(kUgrid));
}

TEST(UgridText, ReadsAQuadWithItsTag) {
    const mt::Mesh m = ug_read("4 0 1 0 0 0 0\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n1 2 3 4\n8\n");
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "quad");
    for (std::size_t k = 0; k < 4; ++k)
        EXPECT_EQ(ug_cell(m, 0, k), static_cast<std::int64_t>(k));
    EXPECT_EQ(read_int(m.CellData("ugrid:ref", 0), 0), 8);
}

TEST(UgridText, CrlfLineEndings) {
    expect_tet(ug_read(ug_replace_all(kUgrid, "\n", "\r\n")));
}

TEST(UgridText, TabsFormFeedsAndVerticalTabsSeparate) {
    std::string text = ug_replace_all(kUgrid, " ", "\t");
    text = ug_replace_all(text, "1\t0\t0\n", "\f1\v0\f0\n");
    expect_tet(ug_read(text));
}

TEST(UgridText, TheLayoutHasNoLineStructure) {
    expect_tet(ug_read(ug_replace_all(kUgrid, "\n", " ")));
    expect_tet(ug_read(ug_replace_all(ug_replace_all(kUgrid, " ", "\n"), "\n\n", "\n")));
    expect_tet(ug_read("  \n\n" + ug_replace_all(kUgrid, "\n", "\n\n  \t")));
}

TEST(UgridText, NumbersAreParsedLeniently) {
    const mt::Mesh m = ug_points("1e1 -2.5E-1 1.5abc\n.5 +3 -.5e1\n", 2);
    const double expected[6] = {10, -0.25, 1.5, 0.5, 3, -5};
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_EQ(read_double(m.Points(), i), expected[i]) << i;
}

TEST(UgridText, ATokenWithoutANumberReadsAsZero) {
    const mt::Mesh m = ug_points("abc - x\n", 1);
    for (std::size_t i = 0; i < 3; ++i)
        EXPECT_EQ(read_double(m.Points(), i), 0.0);
}

TEST(UgridText, SpecialDoubleTokens) {
    const mt::Mesh m = ug_points("nan inf -inf\n1e999 0x10 infinity\n-nan 1e-999 2\n", 3);
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
    EXPECT_EQ(p(8), 2.0);
}

TEST(UgridText, ConnectivityIsCheckedAgainstThePoints) {
    const std::string what = "UGRID: an element references a node outside the point table";
    const std::string head = "4 1 0 1 0 0 0\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n";
    EXPECT_EQ(ug_error(head + "1 2 3\n5\n1 2 3 4\n"), "");
    EXPECT_EQ(ug_error(head + "1 2 5\n5\n1 2 3 4\n"), what);                       // past the end
    EXPECT_EQ(ug_error(head + "0 2 3\n5\n1 2 3 4\n"), what);                       // 1-based
    EXPECT_EQ(ug_error(head + "1 2 -3\n5\n1 2 3 4\n"), what);                      // negative
    EXPECT_EQ(ug_error(head + "1 2 99999999999999999999\n5\n1 2 3 4\n"), what);    // saturated
    EXPECT_EQ(ug_error(head + "1 2 3\n5\n1 2 3 99999999999999999999\n"), what);    // in a tet
    EXPECT_EQ(ug_error(head + "1 2 3\n5\n1 2 3 -99999999999999999999\n"), what);
}

TEST(UgridText, ATagTakesWhatStrtollTakes) {
    const mt::Mesh m = ug_read(ug_replace_all(kUgrid, "\n5\n", "\n99999999999999999999\n"));
    EXPECT_EQ(read_int(m.CellData("ugrid:ref", 0), 0), std::numeric_limits<std::int64_t>::max());
    const mt::Mesh n = ug_read(ug_replace_all(kUgrid, "\n5\n", "\n7xyz\n"));
    EXPECT_EQ(read_int(n.CellData("ugrid:ref", 0), 0), 7);
}

TEST(UgridText, AShortStreamEndsTheFile) {
    const std::string eof = "UGRID: unexpected end of file";
    EXPECT_EQ(ug_error(""), eof);
    EXPECT_EQ(ug_error(" \n\t \r\n"), eof);
    EXPECT_EQ(ug_error("4 1 0 1 0 0\n"), eof);                 // six counts
    EXPECT_EQ(ug_error("1 0 0 0 0 0 0\n0 0\n"), eof);          // two coordinates
    EXPECT_EQ(ug_error("1 1 0 0 0 0 0\n0 0 0\n1 1\n"), eof);   // two connectivity ids
    EXPECT_EQ(ug_error(std::string(kUgrid, sizeof kUgrid - 1 - 4)), eof);  // ids missing
}

TEST(UgridText, AnEmptyMeshHasNoPointsAndNoCells) {
    const mt::Mesh m = ug_read("0 0 0 0 0 0 0\n");
    EXPECT_EQ(m.NumPoints(), 0u);
    EXPECT_EQ(m.NumCellBlocks(), 0u);
}

TEST(UgridText, CountsAreBoundedByTheFile) {
    EXPECT_EQ(ug_error("-1 0 0 0 0 0 0\n"), "meshio++: UGRID: negative header count -1");
    EXPECT_EQ(ug_error("0 0 0 0 0 0 -7\n"), "meshio++: UGRID: negative header count -7");
    const std::string big = ug_error("999999 0 0 0 0 0 0\n0 0 0\n");
    EXPECT_NE(big.find("meshio++: UGRID: header count 999999 exceeds what the file holds ("),
              std::string::npos)
        << big;
}

TEST(UgridText, AMissingFileCannotBeOpened) {
    try {
        meshioplusplus::read_ugrid(mt::temp_path(".ugrid"));
        FAIL() << "expected a ReadError";
    } catch (const meshioplusplus::ReadError& rExc) {
        EXPECT_NE(std::string(rExc.what()).find("Could not open file: "), std::string::npos);
    }
}
