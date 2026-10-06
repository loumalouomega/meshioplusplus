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
#include <stdexcept>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/obj_off.hpp"

// Pins what the OBJ reader accepts, refuses and returns for hand-made text, so
// a change of its tokenizer (roadmap 3.1.1.1) can be checked against the reader
// it replaces. The odd-number cases record what the reader did when this file
// was written, not what the format asks for: face indices are never range
// checked or made relative, a face item is its leading integer, and a face
// with a word where an index belongs throws the C++ library's
// `std::invalid_argument`. Values are read through the dtype-neutral helpers,
// because a mesh backend may store points as float32 or float64.

namespace {

using meshioplusplus::detail::read_double;
using meshioplusplus::detail::read_int;

// An `.obj` file on disk, removed on exit.
class ObjFile {
public:
    explicit ObjFile(const std::string& rText) : mPath(mt::temp_path(".obj")) {
        std::ofstream os(mPath, std::ios::binary);
        os << rText;
    }
    ~ObjFile() {
        std::error_code ec;
        std::filesystem::remove(mPath, ec);
    }
    const std::string& Path() const { return mPath; }

private:
    std::string mPath;
};

mt::Mesh obj_read(const std::string& rText) {
    const ObjFile file(rText);
    return meshioplusplus::read_obj(file.Path());
}

// The `ReadError` message of a file, or "" when it reads.
std::string obj_error(const std::string& rText) {
    const ObjFile file(rText);
    try {
        meshioplusplus::read_obj(file.Path());
    } catch (const meshioplusplus::ReadError& rExc) {
        return rExc.what();
    }
    return "";
}

double obj_point(const mt::Mesh& rMesh, std::size_t Row, std::size_t Col) {
    return read_double(rMesh.Points(), 3 * Row + Col);
}

std::int64_t obj_cell(const mt::Mesh& rMesh, std::size_t Block, std::size_t Index) {
    return read_int(rMesh.Cells(Block).Conn(), Index);
}

std::int64_t obj_group(const mt::Mesh& rMesh, std::size_t Block, std::size_t Index) {
    return read_int(rMesh.CellData("obj:group_ids", Block), Index);
}

// The first `Count` points of a file made of `v` lines only.
mt::Mesh obj_points(const std::string& rBody) {
    return obj_read(rBody);
}

}  // namespace

TEST(ObjText, ReadsATriangle) {
    const mt::Mesh m = obj_read("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n");
    ASSERT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(obj_point(m, 1, 0), 1.0);
    EXPECT_EQ(obj_point(m, 2, 1), 1.0);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "triangle");
    EXPECT_EQ(m.Cells(0).NumCells(), 1u);
    for (std::size_t k = 0; k < 3; ++k)
        EXPECT_EQ(obj_cell(m, 0, k), static_cast<std::int64_t>(k));
    ASSERT_TRUE(m.HasCellData("obj:group_ids"));
    EXPECT_EQ(obj_group(m, 0, 0), -1);
    EXPECT_FALSE(m.HasPointData("obj:vn"));
    EXPECT_FALSE(m.HasPointData("obj:vt"));
}

TEST(ObjText, AnEmptyFileIsAnEmptyMesh) {
    const mt::Mesh m = obj_read("");
    EXPECT_EQ(m.NumPoints(), 0u);
    EXPECT_EQ(m.NumCellBlocks(), 0u);
    EXPECT_FALSE(m.HasCellData("obj:group_ids"));
}

TEST(ObjText, CommentsAndBlankLinesAreSkipped) {
    const mt::Mesh m = obj_read(
        "# a comment\n"
        "\n"
        "   \t  \n"
        "   # an indented comment\n"
        "v 0 0 0\n"
        "v 1 0 0\n"
        "v 0 1 0\n"
        "f 1 2 3\n");
    EXPECT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(m.Cells(0).NumCells(), 1u);
}

TEST(ObjText, ATrailingCommentOnADataLineAddsNothing) {
    const mt::Mesh m = obj_read("v 1 2 3 # corner\nv 4 5 6\nv 7 8 9\nf 1 2 3 # face\n");
    ASSERT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(obj_point(m, 0, 2), 3.0);
    EXPECT_EQ(m.Cells(0).Type(), "triangle");
}

TEST(ObjText, CrlfTabsAndFormFeedsSeparateTokens) {
    const mt::Mesh m = obj_read(
        "v\t0\t0\t0\r\n"
        "v 1\f0\v0\r\n"
        "  v 0 1 0  \r\n"
        "f\t1 2\t3\r\n");
    ASSERT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(obj_point(m, 1, 0), 1.0);
    EXPECT_EQ(obj_point(m, 2, 1), 1.0);
    EXPECT_EQ(m.Cells(0).Type(), "triangle");
    EXPECT_EQ(obj_cell(m, 0, 2), 2);
}

TEST(ObjText, ALastLineWithoutANewlineIsRead) {
    const mt::Mesh m = obj_read("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3");
    EXPECT_EQ(m.Cells(0).NumCells(), 1u);
}

TEST(ObjText, AShortVertexLeavesTheRestZero) {
    const mt::Mesh m = obj_read("v 5\nv 5 6\nv 5 6 7\n");
    ASSERT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(obj_point(m, 0, 0), 5.0);
    EXPECT_EQ(obj_point(m, 0, 1), 0.0);
    EXPECT_EQ(obj_point(m, 0, 2), 0.0);
    EXPECT_EQ(obj_point(m, 1, 1), 6.0);
    EXPECT_EQ(obj_point(m, 1, 2), 0.0);
    EXPECT_EQ(obj_point(m, 2, 2), 7.0);
}

TEST(ObjText, ABareVertexTagIsTheOrigin) {
    const mt::Mesh m = obj_read("v\nv 1 1 1\n");
    ASSERT_EQ(m.NumPoints(), 2u);
    for (std::size_t c = 0; c < 3; ++c)
        EXPECT_EQ(obj_point(m, 0, c), 0.0);
}

TEST(ObjText, ExtraVertexTokensAreIgnored) {
    // A fourth coordinate (`w`) or a vertex colour follows the position.
    const mt::Mesh m = obj_read("v 1 2 3 1.0\nv 4 5 6 0.5 0.5 0.5\n");
    ASSERT_EQ(m.NumPoints(), 2u);
    EXPECT_EQ(obj_point(m, 0, 2), 3.0);
    EXPECT_EQ(obj_point(m, 1, 0), 4.0);
    EXPECT_EQ(obj_point(m, 1, 2), 6.0);
}

TEST(ObjText, ExponentsSignsAndASeparateDecimalPoint) {
    const mt::Mesh m = obj_points("v 1e2 -2.5E-1 +3\nv .5 5. -0\n");
    EXPECT_EQ(obj_point(m, 0, 0), 100.0);
    EXPECT_EQ(obj_point(m, 0, 1), -0.25);
    EXPECT_EQ(obj_point(m, 0, 2), 3.0);
    EXPECT_EQ(obj_point(m, 1, 0), 0.5);
    EXPECT_EQ(obj_point(m, 1, 1), 5.0);
    EXPECT_TRUE(std::signbit(obj_point(m, 1, 2)));
}

TEST(ObjText, NonFiniteOverflowAndUnderflowValues) {
    const mt::Mesh m = obj_points("v inf -inf nan\nv 1e999 -1e999 1e-999\n");
    EXPECT_TRUE(std::isinf(obj_point(m, 0, 0)));
    EXPECT_GT(obj_point(m, 0, 0), 0.0);
    EXPECT_TRUE(std::isinf(obj_point(m, 0, 1)));
    EXPECT_LT(obj_point(m, 0, 1), 0.0);
    EXPECT_TRUE(std::isnan(obj_point(m, 0, 2)));
    EXPECT_TRUE(std::isinf(obj_point(m, 1, 0)));
    EXPECT_TRUE(std::isinf(obj_point(m, 1, 1)));
    EXPECT_EQ(obj_point(m, 1, 2), 0.0);
}

TEST(ObjText, AHexadecimalFloatIsParsed) {
    const mt::Mesh m = obj_points("v 0x1p3 0x.8 0\n");
    EXPECT_EQ(obj_point(m, 0, 0), 8.0);
    EXPECT_EQ(obj_point(m, 0, 1), 0.5);
}

TEST(ObjText, ABadCoordinateStopsTheVertexAndKeepsTheRestZero) {
    // A stream extraction that fails sets its value to zero and ends the line.
    const mt::Mesh m = obj_read("v 1 x 3\nv x 2 3\nv 4 5 y\n");
    ASSERT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(obj_point(m, 0, 0), 1.0);
    EXPECT_EQ(obj_point(m, 0, 1), 0.0);
    EXPECT_EQ(obj_point(m, 0, 2), 0.0);
    for (std::size_t c = 0; c < 3; ++c)
        EXPECT_EQ(obj_point(m, 1, c), 0.0);
    EXPECT_EQ(obj_point(m, 2, 0), 4.0);
    EXPECT_EQ(obj_point(m, 2, 1), 5.0);
    EXPECT_EQ(obj_point(m, 2, 2), 0.0);
}

TEST(ObjText, ATrailingWordEndsANumberedVertexLine) {
    const mt::Mesh m = obj_read("v 1.5abc 2 3\n");
    ASSERT_EQ(m.NumPoints(), 1u);
    EXPECT_EQ(obj_point(m, 0, 0), 1.5);
    EXPECT_EQ(obj_point(m, 0, 1), 0.0);
    EXPECT_EQ(obj_point(m, 0, 2), 0.0);
}

TEST(ObjText, OnlyTheExactTagsCount) {
    // `v1`, `vp`, `o`, `s`, `usemtl`, `mtllib`, `l` and `p` carry nothing.
    const mt::Mesh m = obj_read(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\n"
        "v1 9 9 9\nvp 0.5 0.5\no object\ns off\nusemtl red\nmtllib a.mtl\nl 1 2\np 1\n"
        "f 1 2 3\n");
    EXPECT_EQ(m.NumPoints(), 3u);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).NumCells(), 1u);
}

TEST(ObjText, NormalsAndTextureCoordinatesBecomePointData) {
    const mt::Mesh m = obj_read(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\n"
        "vn 0 0 1\nvn 0 0 1\nvn 0 0 -1\n"
        "vt 0 0\nvt 1 0\nvt 0 1\n"
        "f 1/1/1 2/2/2 3/3/3\n");
    ASSERT_TRUE(m.HasPointData("obj:vn"));
    ASSERT_TRUE(m.HasPointData("obj:vt"));
    const auto& vn = m.PointData("obj:vn");
    const auto& vt = m.PointData("obj:vt");
    ASSERT_EQ(vn.Shape().size(), 2u);
    EXPECT_EQ(vn.Shape()[0], 3u);
    EXPECT_EQ(vn.Shape()[1], 3u);
    EXPECT_EQ(vt.Shape()[0], 3u);
    EXPECT_EQ(vt.Shape()[1], 2u);
    EXPECT_EQ(read_double(vn, 8), -1.0);
    EXPECT_EQ(read_double(vt, 2), 1.0);
    EXPECT_EQ(read_double(vt, 5), 1.0);
}

TEST(ObjText, AttributeLinesAreCountedByTheirOwnTokens) {
    // `vt` may carry one, two or three components; `vn` any number.
    const mt::Mesh one = obj_read("vt 0.5\nvt 0.25\n");
    EXPECT_EQ(one.PointData("obj:vt").Shape()[1], 1u);
    const mt::Mesh three = obj_read("vt 0 1 2\nvn 1 2 3 4\n");
    EXPECT_EQ(three.PointData("obj:vt").Shape()[1], 3u);
    EXPECT_EQ(three.PointData("obj:vn").Shape()[1], 4u);
}

TEST(ObjText, AnAttributeLineStopsAtTheFirstNonNumber) {
    const mt::Mesh m = obj_read("vn 1 2 x 4\nvn 5 6\n");
    ASSERT_TRUE(m.HasPointData("obj:vn"));
    EXPECT_EQ(m.PointData("obj:vn").Shape()[0], 2u);
    EXPECT_EQ(m.PointData("obj:vn").Shape()[1], 2u);
    EXPECT_EQ(read_double(m.PointData("obj:vn"), 1), 2.0);
}

TEST(ObjText, RowsOfOneAttributeWithDifferentLengthsAreRefused) {
    EXPECT_EQ(obj_error("vt 0 0\nvt 1 0 0\n"),
              "OBJ: rows of one attribute with different lengths");
    EXPECT_EQ(obj_error("vn 0 0 1\nvn 0 1\n"),
              "OBJ: rows of one attribute with different lengths");
    // One attribute at a time: a `vt` of one width beside a `vn` of another is fine.
    EXPECT_EQ(obj_error("vt 0 0\nvn 0 0 1\n"), "");
}

TEST(ObjText, FaceItemsAreTheirLeadingInteger) {
    const mt::Mesh m = obj_read(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\n"
        "f 1/1/1 2//2 3/3\n"
        "f 1 2 3\n");
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    ASSERT_EQ(m.Cells(0).NumCells(), 2u);
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t k = 0; k < 3; ++k)
            EXPECT_EQ(obj_cell(m, 0, 3 * c + k), static_cast<std::int64_t>(k));
}

TEST(ObjText, FaceIndicesAreNotRangeCheckedOrMadeRelative) {
    const mt::Mesh m = obj_read("v 0 0 0\nf 1 2 99\nf -1 -2 -3\nf 0 +1 07\n");
    ASSERT_EQ(m.Cells(0).NumCells(), 3u);
    EXPECT_EQ(obj_cell(m, 0, 2), 98);
    EXPECT_EQ(obj_cell(m, 0, 3), -2);
    EXPECT_EQ(obj_cell(m, 0, 4), -3);
    EXPECT_EQ(obj_cell(m, 0, 5), -4);
    EXPECT_EQ(obj_cell(m, 0, 6), -1);
    EXPECT_EQ(obj_cell(m, 0, 7), 0);
    EXPECT_EQ(obj_cell(m, 0, 8), 6);
}

TEST(ObjText, AFaceItemKeepsOnlyTheIntegerBeforeAnythingElse) {
    // `stoll` stops at the first character that is not part of the number.
    const mt::Mesh m = obj_read("f 1x 2.9 3e5\n");
    ASSERT_EQ(m.Cells(0).NumCells(), 1u);
    EXPECT_EQ(obj_cell(m, 0, 0), 0);
    EXPECT_EQ(obj_cell(m, 0, 1), 1);
    EXPECT_EQ(obj_cell(m, 0, 2), 2);
}

TEST(ObjText, AFaceItemWithoutADigitThrowsTheLibraryException) {
    EXPECT_THROW(obj_read("f a b c\n"), std::invalid_argument);
    EXPECT_THROW(obj_read("f 1 2 /3\n"), std::invalid_argument);
    EXPECT_THROW(obj_read("f 1 - 3\n"), std::invalid_argument);
}

TEST(ObjText, AFaceIndexOutOfInt64RangeThrowsTheLibraryException) {
    EXPECT_THROW(obj_read("f 1 2 99999999999999999999\n"), std::out_of_range);
}

TEST(ObjText, FaceSizesSplitBlocksAndNameTheirTypes) {
    const mt::Mesh m = obj_read(
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 2 0 0\n"
        "f 1 2 3\nf 1 3 4\n"
        "f 1 2 3 4\n"
        "f 1 2 3 4 5\nf 1 2 3 4 5\n"
        "f 1 2 3\n");
    ASSERT_EQ(m.NumCellBlocks(), 4u);
    EXPECT_EQ(m.Cells(0).Type(), "triangle");
    EXPECT_EQ(m.Cells(0).NumCells(), 2u);
    EXPECT_EQ(m.Cells(1).Type(), "quad");
    EXPECT_EQ(m.Cells(1).NumCells(), 1u);
    EXPECT_EQ(m.Cells(2).Type(), "polygon");
    EXPECT_EQ(m.Cells(2).NumCells(), 2u);
    EXPECT_EQ(obj_cell(m, 2, 9), 4);
    EXPECT_EQ(m.Cells(3).Type(), "triangle");
    EXPECT_EQ(m.Cells(3).NumCells(), 1u);
}

TEST(ObjText, GroupsStartBlocksAndCountFromZero) {
    const mt::Mesh m = obj_read(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\n"
        "f 1 2 3\n"
        "g first\n"
        "f 1 2 3\nf 3 2 1\n"
        "g second\n"
        "f 1 2 3\n"
        "g trailing\n");
    ASSERT_EQ(m.NumCellBlocks(), 3u);
    EXPECT_EQ(m.Cells(0).NumCells(), 1u);
    EXPECT_EQ(m.Cells(1).NumCells(), 2u);
    EXPECT_EQ(m.Cells(2).NumCells(), 1u);
    EXPECT_EQ(obj_group(m, 0, 0), -1);
    EXPECT_EQ(obj_group(m, 1, 0), 0);
    EXPECT_EQ(obj_group(m, 1, 1), 0);
    EXPECT_EQ(obj_group(m, 2, 0), 1);
}

TEST(ObjText, AGroupWithNoFacesIsDropped) {
    const mt::Mesh m = obj_read(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\n"
        "g a\ng b\ng c\nf 1 2 3\n");
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(obj_group(m, 0, 0), 2);
}

TEST(ObjText, AGroupLineWithNoNameStillStartsAGroup) {
    const mt::Mesh m = obj_read("v 0 0 0\nv 1 0 0\nv 0 1 0\ng\nf 1 2 3\n");
    EXPECT_EQ(obj_group(m, 0, 0), 0);
}

TEST(ObjText, AFaceWithNoItemsIsAnEmptyPolygon) {
    const mt::Mesh m = obj_read("v 0 0 0\nf\nf\n");
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "polygon");
    EXPECT_EQ(m.Cells(0).NumCells(), 2u);
}

TEST(ObjText, ALongTokenIsParsed) {
    const std::string digits(80, '1');
    const mt::Mesh m = obj_read("v 0." + digits + " 0 0\nf 1 2 3\n");
    EXPECT_NEAR(obj_point(m, 0, 0), 0.1111111111111111, 1e-15);
}

TEST(ObjText, ManyLinesKeepTheirOrder) {
    std::string text;
    constexpr std::size_t n = 2000;
    for (std::size_t i = 0; i < n; ++i)
        text += "v " + std::to_string(i) + " " + std::to_string(2 * i) + " 0.5\n";
    for (std::size_t i = 0; i + 2 < n; ++i)
        text += "f " + std::to_string(i + 1) + "/1 " + std::to_string(i + 2) + "/1 " +
                std::to_string(i + 3) + "/1\n";
    const mt::Mesh m = obj_read(text);
    ASSERT_EQ(m.NumPoints(), n);
    EXPECT_EQ(obj_point(m, n - 1, 1), static_cast<double>(2 * (n - 1)));
    ASSERT_EQ(m.Cells(0).NumCells(), n - 2);
    EXPECT_EQ(obj_cell(m, 0, 3 * (n - 3) + 2), static_cast<std::int64_t>(n - 1));
}

TEST(ObjText, AMissingFileIsAReadError) {
    EXPECT_THROW(meshioplusplus::read_obj("/nonexistent/dir/mesh.obj"),
                 meshioplusplus::ReadError);
}
