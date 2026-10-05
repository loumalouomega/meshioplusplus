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
#include "meshioplusplus/formats/tetgen.hpp"

// Pins what the TetGen text reader accepts, refuses and returns for hand-made
// `.node` / `.ele` text, so a change of its tokenizer (roadmap 3.1.1.1) can be
// checked against the reader it replaces. The odd-number cases record what the
// reader did when this file was written.

namespace {

// A `.node` / `.ele` pair on disk, removed on scope exit. `Ele` is absent for
// a lone `.node` file.
class TetgenFiles {
public:
    TetgenFiles(const std::string& rNode, const std::optional<std::string>& rEle)
        : mStem(mt::temp_path("_tg")) {
        write(mStem + ".node", rNode);
        if (rEle)
            write(mStem + ".ele", *rEle);
    }
    ~TetgenFiles() {
        std::error_code ec;
        std::filesystem::remove(mStem + ".node", ec);
        std::filesystem::remove(mStem + ".ele", ec);
    }
    std::string Node() const { return mStem + ".node"; }
    std::string Ele() const { return mStem + ".ele"; }

private:
    static void write(const std::string& rPath, const std::string& rText) {
        std::ofstream os(rPath, std::ios::binary);
        os << rText;
    }
    std::string mStem;
};

mt::Mesh read_pair(const std::string& rNode, const std::string& rEle) {
    const TetgenFiles files(rNode, rEle);
    return meshioplusplus::read_tetgen(files.Node());
}

// The message of the `ReadError` a pair raises, or "" when it reads.
std::string read_error(const std::string& rNode, const std::optional<std::string>& rEle) {
    const TetgenFiles files(rNode, rEle);
    try {
        meshioplusplus::read_tetgen(files.Node());
    } catch (const meshioplusplus::ReadError& rExc) {
        return rExc.what();
    }
    return "";
}

// One unit tetrahedron with a point attribute and marker and a region number.
const char kNode[] =
    "# a comment\n"
    "\n"
    "   # an indented comment\n"
    "4 3 1 1\n"
    "1 0 0 0 0.5 7\n"
    "2 1 0 0 1.5 8\n"
    "3 0 1 0 2.5 9\n"
    "4 0 0 1 3.5 10\n";
const char kEle[] =
    "# elements\n"
    "1 4 1\n"
    "1 1 2 3 4 42\n";

void expect_unit_tet(const mt::Mesh& rMesh) {
    ASSERT_EQ(rMesh.NumPoints(), 4u);
    const double* p = rMesh.Points().As<double>();
    const double expected[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t d = 0; d < 3; ++d)
            EXPECT_EQ(p[3 * i + d], expected[i][d]);
    const double* attr = rMesh.PointData("tetgen:attr1").As<double>();
    const double* ref = rMesh.PointData("tetgen:ref").As<double>();
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_EQ(attr[i], 0.5 + static_cast<double>(i));
        EXPECT_EQ(ref[i], 7.0 + static_cast<double>(i));
    }
    ASSERT_EQ(rMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(rMesh.Cells(0).Type(), "tetra");
    ASSERT_EQ(rMesh.Cells(0).NumCells(), 1u);
    const std::int64_t* c = rMesh.Cells(0).Conn().As<std::int64_t>();
    for (std::int64_t k = 0; k < 4; ++k)
        EXPECT_EQ(c[k], k);
    ASSERT_TRUE(rMesh.HasCellData("tetgen:ref"));
    EXPECT_EQ(rMesh.CellData("tetgen:ref", 0).As<std::int64_t>()[0], 42);
}

std::string replace_all(std::string Text, const std::string& rFrom, const std::string& rTo) {
    for (std::size_t pos = 0; (pos = Text.find(rFrom, pos)) != std::string::npos;
         pos += rTo.size())
        Text.replace(pos, rFrom.size(), rTo);
    return Text;
}

// One zero-cell tetrahedral file, so a test can read just its points.
mt::Mesh read_points(const std::string& rNodeBody, std::size_t Count) {
    return read_pair(std::to_string(Count) + " 3 0 0\n" + rNodeBody, "0 4 0\n");
}

}  // namespace

TEST(Tetgen, ReadsAnAttributedTetrahedron) {
    expect_unit_tet(read_pair(kNode, kEle));
}

TEST(Tetgen, AnEleFileNamesTheSamePair) {
    const TetgenFiles files(kNode, kEle);
    expect_unit_tet(meshioplusplus::read_tetgen(files.Ele()));
}

TEST(Tetgen, CrlfLineEndings) {
    expect_unit_tet(read_pair(replace_all(kNode, "\n", "\r\n"), replace_all(kEle, "\n", "\r\n")));
}

TEST(Tetgen, TabsFormFeedsAndVerticalTabsSeparateAndIndent) {
    std::string node = replace_all(kNode, " ", "\t");
    node = replace_all(node, "1\t0\t0\t0\t0.5", "\f1\v0\f0\v0\f0.5");
    expect_unit_tet(read_pair(node, replace_all(kEle, " ", "\t")));
}

TEST(Tetgen, DataTokensRunAcrossLines) {
    // Every token on its own line: the data are one stream after the header.
    std::string node = "4 3 1 1\n";
    for (const char* tok : {"1", "0", "0", "0", "0.5", "7", "2", "1", "0", "0", "1.5", "8", "3", "0",
                            "1", "0", "2.5", "9", "4", "0", "0", "1", "3.5", "10"})
        node += std::string(tok) + "\n";
    expect_unit_tet(read_pair(node, "1 4 1\n1\n1\n2\n3\n4\n42\n"));
}

TEST(Tetgen, AnInlineCommentOnTheHeaderLineIsTolerated) {
    // Header tokens past the fourth are ignored.
    expect_unit_tet(read_pair(replace_all(kNode, "4 3 1 1\n", "4 3 1 1 # points\n"),
                              replace_all(kEle, "1 4 1\n", "1 4 1 # tets\n")));
}

TEST(Tetgen, AnInlineCommentAfterDataIsNotAComment) {
    // Only a line that *starts* with # is skipped; its tokens here are data.
    EXPECT_EQ(read_error(replace_all(kNode, "1 0 0 0 0.5 7\n", "1 0 0 0 0.5 7 # first\n"), kEle),
              "TetGen: .node data size mismatch");
}

TEST(Tetgen, NumbersAreParsedLeniently) {
    const mt::Mesh m = read_points("1 1e1 -2.5E-1 +3\n"
                                   "2 1.5abc .5 -.5e1\n"
                                   "3 abc 0 12\n",
                                   3);
    const double* p = m.Points().As<double>();
    const double expected[3][3] = {{10, -0.25, 3}, {1.5, 0.5, -5}, {0, 0, 12}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t d = 0; d < 3; ++d)
            EXPECT_EQ(p[3 * i + d], expected[i][d]) << i << ',' << d;
}

TEST(Tetgen, ANodeIndexTakesItsIntegerPart) {
    // The index is a number cast to an integer: 1.9 numbers the first node 1.
    const mt::Mesh m = read_points("1.9 0 0 0\n2 1 0 0\n", 2);
    EXPECT_EQ(m.NumPoints(), 2u);
}

TEST(Tetgen, ZeroBasedNodeNumberingShiftsTheConnectivity) {
    const mt::Mesh m = read_pair("4 3 0 0\n0 0 0 0\n1 1 0 0\n2 0 1 0\n3 0 0 1\n", "1 4 0\n0 0 1 2 3\n");
    ASSERT_EQ(m.Cells(0).NumCells(), 1u);
    const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
    for (std::int64_t k = 0; k < 4; ++k)
        EXPECT_EQ(c[k], k);
}

TEST(Tetgen, AnOverflowingIntegerSaturates) {
    const mt::Mesh m = read_pair("4 3 0 0\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n",
                                 "1 4 0\n1 1 2 3 99999999999999999999\n");
    const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
    EXPECT_EQ(c[3], std::numeric_limits<std::int64_t>::max() - 1);
}

TEST(Tetgen, SpecialDoubleTokens) {
    const mt::Mesh m = read_points("1 nan inf 0\n"
                                   "2 -inf 1e999 0x10\n"
                                   "3 infinity -nan 1e-999\n",
                                   3);
    const double* p = m.Points().As<double>();
    // Recorded from the reader as it was: nan and inf words, an overflow to
    // infinity, a hexadecimal float, and an underflow to zero.
    EXPECT_TRUE(std::isnan(p[0]));
    EXPECT_EQ(p[1], std::numeric_limits<double>::infinity());
    EXPECT_EQ(p[2], 0.0);
    EXPECT_EQ(p[3], -std::numeric_limits<double>::infinity());
    EXPECT_EQ(p[4], std::numeric_limits<double>::infinity());
    EXPECT_EQ(p[5], 16.0);
    EXPECT_EQ(p[6], std::numeric_limits<double>::infinity());
    EXPECT_TRUE(std::isnan(p[7]));
    EXPECT_TRUE(std::signbit(p[7]));
    EXPECT_EQ(p[8], 0.0);
}

TEST(Tetgen, AnEmptyMeshHasNoPointsAndNoCells) {
    const mt::Mesh m = read_pair("0 3 0 0\n", "0 4 0\n");
    EXPECT_EQ(m.NumPoints(), 0u);
}

TEST(Tetgen, RefusedInputsNameTheirReason) {
    EXPECT_EQ(read_error("4 3 1\n", kEle), "TetGen: malformed .node header");
    EXPECT_EQ(read_error("4 2 0 0\n", kEle), "TetGen: need 3D points");
    EXPECT_EQ(read_error(replace_all(kNode, "4 0 0 1 3.5 10\n", ""), kEle),
              "TetGen: .node data size mismatch");
    EXPECT_EQ(read_error(kNode, replace_all(kEle, "1 4 1\n", "1 4\n")),
              "TetGen: malformed .ele header");
    EXPECT_EQ(read_error(kNode, replace_all(kEle, "1 4 1\n", "1 10 1\n")),
              "TetGen: only 4-node tetrahedra supported");
    EXPECT_EQ(read_error(kNode, replace_all(kEle, "1 1 2 3 4 42\n", "1 1 2 3 4\n")),
              "TetGen: .ele data size mismatch");
    EXPECT_EQ(read_error(replace_all(kNode, "3 0 1 0 2.5 9\n", "5 0 1 0 2.5 9\n"), kEle),
              "TetGen: nodes not numbered consecutively");
}

TEST(Tetgen, AFileOfOnlyCommentsHasNoHeader) {
    const std::string what = read_error("# nothing\n\n   \n# at all\n", kEle);
    EXPECT_NE(what.find("TetGen: missing header line in "), std::string::npos) << what;
    EXPECT_NE(read_error("", kEle).find("TetGen: missing header line in "), std::string::npos);
}

TEST(Tetgen, AMissingElementFileCannotBeOpened) {
    const std::string what = read_error(kNode, std::nullopt);
    EXPECT_NE(what.find("Could not open file: "), std::string::npos) << what;
    EXPECT_NE(what.find(".ele"), std::string::npos) << what;
}

TEST(Tetgen, OnlyNodeAndEleSuffixesAreRead) {
    try {
        meshioplusplus::read_tetgen(mt::temp_path(".vtk"));
        FAIL() << "expected a ReadError";
    } catch (const meshioplusplus::ReadError& rExc) {
        EXPECT_STREQ(rExc.what(), "TetGen: expected a .node or .ele file");
    }
}
