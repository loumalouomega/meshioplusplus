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
 * @file test_frd.cpp
 * @brief CalculiX `.frd` reader and ASCII writer: node permutations, steps, short and
 *        long formats, derived invariants, the write/read round trip and the error paths.
 *
 * The files are built here from fixed-column records rather than read from
 * tests/python/meshes/frd/: this suite has no test-data path. The Python suite runs the
 * same reader over real `ccx` output.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <ios>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/frd.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::FrdReadOptions;
using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;

std::string ipad(std::int64_t value, std::size_t width) {
    std::string s = std::to_string(value);
    return std::string(width > s.size() ? width - s.size() : 0, ' ') + s;
}

std::string e12(double value) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%12.5E", value);
    return buf;
}

/// Builds `.frd` text in the short (I5) or long (I10) layout.
class FrdText {
public:
    explicit FrdText(bool longFormat = true) : mLong(longFormat) {}

    FrdText& Nodes(const std::vector<std::array<double, 3>>& rNodes) {
        mText += "    1C\n    1UTEST\n    2C" + ipad(static_cast<std::int64_t>(rNodes.size()), 30) +
                 std::string(37, ' ') + (mLong ? "1" : "0") + "\n";
        for (std::size_t i = 0; i < rNodes.size(); ++i)
            mText += " -1" + Id(static_cast<std::int64_t>(i) + 1) + e12(rNodes[i][0]) +
                     e12(rNodes[i][1]) + e12(rNodes[i][2]) + "\n";
        mText += " -3\n";
        return *this;
    }

    FrdText& Elements(const std::vector<std::pair<int, std::vector<std::int64_t>>>& rElements) {
        mText += "    3C" + ipad(static_cast<std::int64_t>(rElements.size()), 30) +
                 std::string(37, ' ') + (mLong ? "1" : "0") + "\n";
        std::int64_t eid = 1;
        for (const auto& [type, nodes] : rElements) {
            mText += " -1" + Id(eid++) + ipad(type, 5) + ipad(0, 5) + ipad(1, 5) + "\n";
            const std::size_t per_line = mLong ? 10 : 15;
            for (std::size_t k = 0; k < nodes.size(); k += per_line) {
                mText += " -2";
                for (std::size_t j = k; j < std::min(nodes.size(), k + per_line); ++j)
                    mText += Id(nodes[j]);
                mText += "\n";
            }
        }
        mText += " -3\n";
        return *this;
    }

    /// One `100C` header and result block. `Rows` are (node id, values); values wrap at six.
    FrdText& Result(double value, int analysis, int stepNumber, const std::string& rName,
                    const std::vector<std::string>& rComponents,
                    const std::vector<std::pair<std::int64_t, std::vector<double>>>& rRows,
                    bool withCalculated = false) {
        mText += "    1PSTEP" + std::string(25, ' ') + "1           1           1          \n";
        char head[64];
        std::snprintf(head, sizeof(head), "%12.9f", value);
        mText += "  100CL " + ipad(100 + stepNumber, 4) + head +
                 ipad(static_cast<std::int64_t>(rRows.size()), 12) + std::string(20, ' ') +
                 ipad(analysis, 2) + ipad(stepNumber, 5) + std::string(10, ' ') +
                 ipad(mLong ? 1 : 0, 2) + "\n";
        const std::size_t declared = rComponents.size() + (withCalculated ? 1 : 0);
        mText += " -4  " + Left(rName) + ipad(static_cast<std::int64_t>(declared), 5) + ipad(1, 5) +
                 "\n";
        for (std::size_t c = 0; c < rComponents.size(); ++c)
            mText += " -5  " + Left(rComponents[c]) + ipad(1, 5) + ipad(1, 5) +
                     ipad(static_cast<std::int64_t>(c) + 1, 5) + ipad(0, 5) + ipad(0, 5) + "\n";
        if (withCalculated)
            mText += " -5  " + Left("ALL") + ipad(1, 5) + ipad(2, 5) + ipad(0, 5) + ipad(0, 5) +
                     ipad(1, 5) + "ALL\n";
        for (const auto& [node, values] : rRows) {
            mText += " -1" + Id(node);
            for (std::size_t k = 0; k < values.size(); ++k) {
                if (k > 0 && k % 6 == 0)
                    mText += "\n -2" + std::string(mLong ? 10 : 5, ' ');
                mText += e12(values[k]);
            }
            mText += "\n";
        }
        mText += " -3\n";
        return *this;
    }

    std::string Finish() const { return mText + "  9999\n"; }

private:
    bool mLong;
    std::string mText;

    std::string Id(std::int64_t id) const { return ipad(id, mLong ? 10 : 5); }
    static std::string Left(const std::string& rText) {
        return rText + std::string(rText.size() < 8 ? 8 - rText.size() : 0, ' ');
    }
};

/// Appends the raw little-endian bytes of @p Value to @p rOut.
template <class T>
void append_raw(std::string& rOut, T Value) {
    rOut.append(reinterpret_cast<const char*>(&Value), sizeof(Value));
}

/// Builds the binary `.frd` layout (`*NODE OUTPUT`/`*ELEMENT OUTPUT`) by hand,
/// mirroring FrdText's column layout exactly (see doc/formats/frd.md and the real
/// ccx 2.23 output this was verified against). No test-data path: everything is
/// built here, like FrdText.
class FrdBinaryText {
public:
    FrdBinaryText& Nodes(const std::vector<std::array<double, 3>>& rNodes,
                         std::size_t RealBytes = 8) {
        mText += "    2C" + ipad(static_cast<std::int64_t>(rNodes.size()), 30) +
                 std::string(37, ' ') + (RealBytes == 8 ? "3" : "2") + "\n";
        for (std::size_t i = 0; i < rNodes.size(); ++i) {
            append_raw(mText, static_cast<std::int32_t>(i + 1));
            for (int k = 0; k < 3; ++k)
                AppendReal(rNodes[i][static_cast<std::size_t>(k)], RealBytes);
        }
        return *this;
    }

    FrdBinaryText& Elements(const std::vector<std::pair<int, std::vector<std::int64_t>>>& rElements) {
        mText += "    3C" + ipad(static_cast<std::int64_t>(rElements.size()), 30) +
                 std::string(37, ' ') + "2\n";
        std::int32_t eid = 1;
        for (const auto& [type, nodes] : rElements) {
            append_raw(mText, eid++);
            append_raw(mText, static_cast<std::int32_t>(type));
            append_raw(mText, static_cast<std::int32_t>(0));  // group
            append_raw(mText, static_cast<std::int32_t>(1));  // material
            for (std::int64_t n : nodes)
                append_raw(mText, static_cast<std::int32_t>(n));
        }
        return *this;
    }

    /// One `100C` header and result block, values in a single `-4` block (no `-2`
    /// continuation needed: the binary layout has no line wrapping at all).
    FrdBinaryText& Result(double value, int analysis, int stepNumber, const std::string& rName,
                          const std::vector<std::string>& rComponents,
                          const std::vector<std::pair<std::int64_t, std::vector<double>>>& rRows,
                          std::size_t RealBytes = 4, bool withCalculated = false) {
        mText += "    1PSTEP" + std::string(25, ' ') + "1           1           1          \n";
        char head[64];
        std::snprintf(head, sizeof(head), "%12.9f", value);
        mText += "  100CL " + ipad(100 + stepNumber, 4) + head +
                 ipad(static_cast<std::int64_t>(rRows.size()), 12) + std::string(20, ' ') +
                 ipad(analysis, 2) + ipad(stepNumber, 5) + std::string(10, ' ') +
                 ipad(RealBytes == 8 ? 3 : 2, 2) + "\n";
        const std::size_t declared = rComponents.size() + (withCalculated ? 1 : 0);
        mText += " -4  " + Left(rName) + ipad(static_cast<std::int64_t>(declared), 5) + ipad(1, 5) +
                 "\n";
        for (std::size_t c = 0; c < rComponents.size(); ++c)
            mText += " -5  " + Left(rComponents[c]) + ipad(1, 5) + ipad(1, 5) +
                     ipad(static_cast<std::int64_t>(c) + 1, 5) + ipad(0, 5) + ipad(0, 5) + "\n";
        if (withCalculated)
            mText += " -5  " + Left("ALL") + ipad(1, 5) + ipad(2, 5) + ipad(0, 5) + ipad(0, 5) +
                     ipad(1, 5) + "ALL\n";
        for (const auto& [node, values] : rRows) {
            append_raw(mText, static_cast<std::int32_t>(node));
            for (double v : values)
                AppendReal(v, RealBytes);
        }
        return *this;
    }

    std::string Finish() const { return mText + " 9999\n"; }

private:
    std::string mText;

    void AppendReal(double Value, std::size_t RealBytes) {
        if (RealBytes == 8)
            append_raw(mText, Value);
        else
            append_raw(mText, static_cast<float>(Value));
    }
    static std::string Left(const std::string& rText) {
        return rText + std::string(rText.size() < 8 ? 8 - rText.size() : 0, ' ');
    }
};

std::string text_file(const std::string& rText) {
    const std::string path = mt::temp_path(".frd");
    std::ofstream out(path, std::ios::binary);
    out << rText;
    return path;
}

struct Temp {
    std::string mPath;
    explicit Temp(const std::string& rText) : mPath(text_file(rText)) {}
    ~Temp() { std::remove(mPath.c_str()); }
};

std::vector<std::int64_t> row_of(const Mesh& rMesh, std::size_t block, std::size_t row) {
    const auto cb = rMesh.Cells(block);
    const std::size_t k = cb.NodesPerCell();
    std::vector<std::int64_t> out;
    for (std::size_t c = 0; c < k; ++c)
        out.push_back(meshioplusplus::detail::read_int(cb.Conn(), row * k + c));
    return out;
}

double at(const meshioplusplus::NDArray& rArray, std::size_t i) {
    return meshioplusplus::detail::read_double(rArray, i);
}

using Nodes = std::vector<std::array<double, 3>>;
using V = std::vector<std::int64_t>;

// Reference element nodes, shifted so no two are the same point.
Nodes numbered_nodes(std::size_t n) {
    Nodes out;
    for (std::size_t i = 0; i < n; ++i)
        out.push_back({static_cast<double>(i), 0.5 * static_cast<double>(i * i), -1.0});
    return out;
}

V iota(std::int64_t n) {
    V out;
    for (std::int64_t i = 1; i <= n; ++i)
        out.push_back(i);
    return out;
}

}  // namespace

TEST(FrdRead, Hexahedron20MidEdgeGroupsAreSwapped) {
    // file order: 12 mid-edge nodes of the end faces, then the 4 verticals, then the top ring
    Temp t(FrdText().Nodes(numbered_nodes(20)).Elements({{4, iota(20)}}).Finish());
    const Mesh mesh = meshioplusplus::read_frd(t.mPath);
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron20");
    EXPECT_EQ(row_of(mesh, 0, 0),
              (V{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15}));
}

TEST(FrdRead, Wedge15AndLine3Permutations) {
    Temp t(
        FrdText().Nodes(numbered_nodes(18)).Elements({{5, iota(15)}, {12, {16, 17, 18}}}).Finish());
    const Mesh mesh = meshioplusplus::read_frd(t.mPath);
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    EXPECT_EQ(mesh.Cells(0).Type(), "wedge15");
    EXPECT_EQ(row_of(mesh, 0, 0), (V{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11}));
    EXPECT_EQ(mesh.Cells(1).Type(), "line3");
    EXPECT_EQ(row_of(mesh, 1, 0), (V{15, 17, 16}));  // end, mid, end -> end, end, mid
}

TEST(FrdRead, EveryElementTypeMapsToItsCellType) {
    const std::vector<std::pair<int, std::string>> types = {
        {1, "hexahedron"}, {2, "wedge"},   {3, "tetra"},    {4, "hexahedron20"},
        {5, "wedge15"},    {6, "tetra10"}, {7, "triangle"}, {8, "triangle6"},
        {9, "quad"},       {10, "quad8"},  {11, "line"},    {12, "line3"}};
    const std::vector<std::size_t> counts = {8, 6, 4, 20, 15, 10, 3, 6, 4, 8, 2, 3};
    std::vector<std::pair<int, V>> elements;
    for (std::size_t i = 0; i < types.size(); ++i)
        elements.push_back({types[i].first, iota(static_cast<std::int64_t>(counts[i]))});
    Temp t(FrdText().Nodes(numbered_nodes(20)).Elements(elements).Finish());
    const Mesh mesh = meshioplusplus::read_frd(t.mPath);
    ASSERT_EQ(mesh.NumCellBlocks(), types.size());
    for (std::size_t i = 0; i < types.size(); ++i)
        EXPECT_EQ(mesh.Cells(i).Type(), types[i].second) << types[i].first;
}

TEST(FrdRead, ConsecutiveElementsOfOneTypeShareABlock) {
    Temp t(FrdText()
               .Nodes(numbered_nodes(6))
               .Elements({{3, {1, 2, 3, 4}}, {3, {2, 3, 4, 5}}, {7, {1, 2, 3}}, {3, {3, 4, 5, 6}}})
               .Finish());
    const Mesh mesh = meshioplusplus::read_frd(t.mPath);
    ASSERT_EQ(mesh.NumCellBlocks(), 3u);
    EXPECT_EQ(mesh.Cells(0).NumCells(), 2u);
    EXPECT_EQ(mesh.Cells(1).Type(), "triangle");
    EXPECT_EQ(mesh.Cells(2).NumCells(), 1u);
}

TEST(FrdRead, ShortAndLongFormatsAgree) {
    const V hex = iota(8);
    const std::vector<std::pair<std::int64_t, std::vector<double>>> disp = {
        {1, {0.1, -0.2, 0.3}}, {2, {-1.5e-7, 2.0, -3.0}},
        {3, {0, 0, 0}},        {4, {1, 2, 3}},
        {5, {4, 5, 6}},        {6, {7, 8, 9}},
        {7, {1, 1, 1}},        {8, {-1, -1, -1}}};
    Mesh meshes[2];
    for (int i = 0; i < 2; ++i) {
        Temp t(FrdText(i == 1)
                   .Nodes(numbered_nodes(8))
                   .Elements({{1, hex}})
                   .Result(0.5, 0, 1, "DISP", {"D1", "D2", "D3"}, disp, true)
                   .Finish());
        meshes[i] = meshioplusplus::read_frd(t.mPath);
    }
    ASSERT_TRUE(meshes[0].HasPointData("DISP"));
    for (std::size_t k = 0; k < 24; ++k)
        EXPECT_EQ(at(meshes[0].PointData("DISP"), k), at(meshes[1].PointData("DISP"), k));
    EXPECT_EQ(meshes[0].PointData("DISP").Shape(), (std::vector<std::size_t>{8, 3}));
    EXPECT_DOUBLE_EQ(at(meshes[0].PointData("DISP"), 3), -1.5e-7);
}

TEST(FrdRead, BinaryLayoutMatchesTheAsciiRenditionOfTheSameRun) {
    const V hex = iota(8);
    const std::vector<std::pair<std::int64_t, std::vector<double>>> disp = {
        {1, {0.1, -0.2, 0.3}}, {2, {-1.5e-7, 2.0, -3.0}},
        {3, {0, 0, 0}},        {4, {1, 2, 3}},
        {5, {4, 5, 6}},        {6, {7, 8, 9}},
        {7, {1, 1, 1}},        {8, {-1, -1, -1}}};
    const std::vector<std::pair<std::int64_t, std::vector<double>>> stress = {
        {1, {1, 2, 3, 4, 5, 6}}, {2, {6, 5, 4, 3, 2, 1}}, {3, {0, 0, 0, 0, 0, 0}},
        {4, {1, 1, 1, 1, 1, 1}}, {5, {2, 3, 4, 5, 6, 7}}, {6, {-1, -2, -3, -4, -5, -6}},
        {7, {10, 20, 30, 40, 50, 60}}, {8, {-10, -20, -30, -40, -50, -60}}};

    Temp ascii_t(FrdText()
                     .Nodes(numbered_nodes(8))
                     .Elements({{1, hex}})
                     .Result(0.5, 0, 1, "DISP", {"D1", "D2", "D3"}, disp, true)
                     .Result(0.5, 0, 1, "STRESS", {"SXX", "SYY", "SZZ", "SXY", "SYZ", "SZX"},
                             stress)
                     .Finish());
    Temp bin_t(FrdBinaryText()
                   .Nodes(numbered_nodes(8), /*RealBytes=*/8)
                   .Elements({{1, hex}})
                   .Result(0.5, 0, 1, "DISP", {"D1", "D2", "D3"}, disp, /*RealBytes=*/4, true)
                   .Result(0.5, 0, 1, "STRESS", {"SXX", "SYY", "SZZ", "SXY", "SYZ", "SZX"}, stress,
                           /*RealBytes=*/4)
                   .Finish());

    const Mesh ascii_mesh = meshioplusplus::read_frd(ascii_t.mPath);
    const Mesh bin_mesh = meshioplusplus::read_frd(bin_t.mPath);
    mt::expect_same_geometry(ascii_mesh, bin_mesh);
    for (const char* name : {"DISP", "STRESS"}) {
        const auto& a = ascii_mesh.PointData(name);
        const auto& b = bin_mesh.PointData(name);
        ASSERT_EQ(a.Shape(), b.Shape());
        for (std::size_t k = 0; k < a.Size(); ++k)
            EXPECT_NEAR(at(a, k), at(b, k), 1e-4) << name << " index " << k;
    }
    // ASCII derives group/material from its own -1 header line, binary from its own
    // packed header quad -- confirm both paths agree.
    EXPECT_EQ(meshioplusplus::detail::read_int(ascii_mesh.CellData("frd:material", 0), 0),
             meshioplusplus::detail::read_int(bin_mesh.CellData("frd:material", 0), 0));
}

TEST(FrdRead, BinaryDerivedInvariantsMatchAscii) {
    const V hex = iota(8);
    const std::vector<std::pair<std::int64_t, std::vector<double>>> stress = {
        {1, {1, 2, 3, 0.5, 0.6, 0.7}}, {2, {6, 5, 4, 0.1, 0.2, 0.3}},
        {3, {0, 0, 0, 0, 0, 0}},       {4, {1, 1, 1, 1, 1, 1}},
        {5, {2, 3, 4, 5, 6, 7}},       {6, {-1, -2, -3, -4, -5, -6}},
        {7, {10, 20, 30, 40, 50, 60}}, {8, {-10, -20, -30, -40, -50, -60}}};
    Temp t(FrdBinaryText()
               .Nodes(numbered_nodes(8), 8)
               .Elements({{1, hex}})
               .Result(0.5, 0, 1, "STRESS", {"SXX", "SYY", "SZZ", "SXY", "SYZ", "SZX"}, stress, 4)
               .Finish());
    FrdReadOptions opts;
    opts.mDerived = true;
    const Mesh mesh = meshioplusplus::read_frd(t.mPath, ReadOptions{}, opts);
    ASSERT_TRUE(mesh.HasPointData("STRESS_mises"));
    ASSERT_TRUE(mesh.HasPointData("STRESS_principal"));
    // Node 1: xx=1 yy=2 zz=3 xy=0.5 yz=0.6 zx=0.7.
    const double expect =
        std::sqrt(0.5 * ((1 - 2) * (1 - 2) + (2 - 3) * (2 - 3) + (3 - 1) * (3 - 1) +
                         6.0 * (0.5 * 0.5 + 0.6 * 0.6 + 0.7 * 0.7)));
    EXPECT_NEAR(at(mesh.PointData("STRESS_mises"), 0), expect, 1e-4);
}

TEST(FrdRead, BinaryClaimOnAsciiTextFailsCleanly) {
    // A hand-edited ASCII file with the 2C flag flipped to a binary one: the ASCII
    // bytes get misread as binary records and must fail on a structural check
    // (an invalid element type or an out-of-bounds record), never crash.
    std::string text = FrdText().Nodes(numbered_nodes(4)).Elements({{3, iota(4)}}).Finish();
    const std::size_t at2c = text.find("    2C");
    text[text.find('\n', at2c) - 1] = '2';
    Temp t(text);
    EXPECT_THROW(meshioplusplus::read_frd(t.mPath), ReadError);
}

TEST(FrdRead, GluedNegativeValuesAreSlicedByColumn) {
    std::string text =
        FrdText()
            .Nodes(numbered_nodes(4))
            .Elements({{3, iota(4)}})
            .Result(
                1.0, 0, 1, "DISP", {"D1", "D2", "D3"},
                {{1, {-1.0, -2.0, -3.0}}, {2, {-1.18144e-6, 0, 0}}, {3, {0, 0, 0}}, {4, {0, 0, 0}}})
            .Finish();
    EXPECT_NE(text.find("-1.00000E+00-2.00000E+00-3.00000E+00"), std::string::npos);
    Temp t(text);
    const Mesh mesh = meshioplusplus::read_frd(t.mPath);
    EXPECT_DOUBLE_EQ(at(mesh.PointData("DISP"), 1), -2.0);
    EXPECT_DOUBLE_EQ(at(mesh.PointData("DISP"), 3), -1.18144e-6);
}

TEST(FrdRead, ContinuationLinesJoinIntoOneRow) {
    std::vector<double> wide;
    for (int i = 0; i < 8; ++i)
        wide.push_back(10.0 + i);
    Temp t(FrdText(false)
               .Nodes(numbered_nodes(4))
               .Elements({{3, iota(4)}})
               .Result(1.0, 0, 1, "WIDE", {"A", "B", "C", "D", "E", "F", "G", "H"},
                       {{1, wide}, {2, wide}, {3, wide}, {4, wide}})
               .Finish());
    const Mesh mesh = meshioplusplus::read_frd(t.mPath);
    EXPECT_EQ(mesh.PointData("WIDE").Shape(), (std::vector<std::size_t>{4, 8}));
    EXPECT_DOUBLE_EQ(at(mesh.PointData("WIDE"), 7), 17.0);
    EXPECT_DOUBLE_EQ(at(mesh.PointData("WIDE"), 8), 10.0);
}

TEST(FrdRead, StepsAreIncrementsAndNodesAResultOmitsAreNaN) {
    FrdText text;
    text.Nodes(numbered_nodes(4)).Elements({{3, iota(4)}});
    text.Result(0.25, 1, 1, "NDTEMP", {"T"}, {{1, {10}}, {2, {20}}, {3, {30}}, {4, {40}}});
    text.Result(0.75, 1, 2, "NDTEMP", {"T"}, {{1, {11}}, {3, {33}}});
    text.Result(0.75, 1, 2, "ERROR", {"STR(%)"}, {{1, {1}}, {2, {2}}, {3, {3}}, {4, {4}}});
    Temp t(text.Finish());

    const Mesh first = meshioplusplus::read_frd(t.mPath);
    EXPECT_DOUBLE_EQ(at(first.FieldData("meshio:time"), 0), 0.25);
    EXPECT_EQ(meshioplusplus::detail::read_int(first.FieldData("frd:step"), 0), 1);
    EXPECT_EQ(meshioplusplus::detail::read_int(first.FieldData("frd:analysis"), 0), 1);
    EXPECT_EQ(first.PointData("NDTEMP").Shape(), (std::vector<std::size_t>{4}));

    ReadOptions last;
    last.mTimeStep = -1;
    const Mesh second = meshioplusplus::read_frd(t.mPath, last);
    EXPECT_DOUBLE_EQ(at(second.FieldData("meshio:time"), 0), 0.75);
    EXPECT_DOUBLE_EQ(at(second.PointData("NDTEMP"), 2), 33.0);
    EXPECT_TRUE(std::isnan(at(second.PointData("NDTEMP"), 1)));
    EXPECT_TRUE(second.HasPointData("ERROR"));
    EXPECT_FALSE(first.HasPointData("ERROR"));

    const auto meta = meshioplusplus::read_frd_metadata(t.mPath);
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.25, 0.75}));
    EXPECT_EQ(meta.mFormat, "frd");

    ReadOptions bad;
    bad.mTimeStep = 2;
    EXPECT_THROW(meshioplusplus::read_frd(t.mPath, bad), ReadError);
    bad.mTimeStep = -3;
    EXPECT_THROW(meshioplusplus::read_frd(t.mPath, bad), ReadError);
}

TEST(FrdRead, FramesSharingAValueButNotAnIdAreDistinctSteps) {
    // degenerate frequency modes: same frequency, different frame id
    FrdText text;
    text.Nodes(numbered_nodes(4)).Elements({{3, iota(4)}});
    for (int mode = 1; mode <= 2; ++mode)
        text.Result(5.0, 2, mode, "DISP", {"D1", "D2", "D3"},
                    {{1, {double(mode), 0, 0}}, {2, {0, 0, 0}}, {3, {0, 0, 0}}, {4, {0, 0, 0}}});
    Temp t(text.Finish());
    EXPECT_EQ(meshioplusplus::read_frd_metadata(t.mPath).mTimeValues.size(), 2u);
    ReadOptions second;
    second.mTimeStep = 1;
    EXPECT_DOUBLE_EQ(at(meshioplusplus::read_frd(t.mPath, second).PointData("DISP"), 0), 2.0);
}

TEST(FrdRead, SelectiveRead) {
    FrdText text;
    text.Nodes(numbered_nodes(4)).Elements({{3, iota(4)}});
    text.Result(1.0, 0, 1, "DISP", {"D1", "D2", "D3"},
                {{1, {1, 2, 3}}, {2, {1, 2, 3}}, {3, {1, 2, 3}}, {4, {1, 2, 3}}});
    text.Result(1.0, 0, 1, "NDTEMP", {"T"}, {{1, {1}}, {2, {1}}, {3, {1}}, {4, {1}}});
    Temp t(text.Finish());
    ReadOptions only;
    only.mDataArrays = std::vector<std::string>{"NDTEMP"};
    const Mesh narrowed = meshioplusplus::read_frd(t.mPath, only);
    EXPECT_TRUE(narrowed.HasPointData("NDTEMP"));
    EXPECT_FALSE(narrowed.HasPointData("DISP"));
    ReadOptions points;
    points.mPointsOnly = true;
    const Mesh geometry = meshioplusplus::read_frd(t.mPath, points);
    EXPECT_EQ(geometry.NumPointData(), 0u);
    EXPECT_TRUE(geometry.HasFieldData("meshio:time"));
}

TEST(FrdRead, DerivedInvariants) {
    // uniaxial xx = 100, then a pure shear xy = 50: mises 100 and 50*sqrt(3)
    const std::vector<std::pair<std::int64_t, std::vector<double>>> rows = {
        {1, {100, 0, 0, 0, 0, 0}},
        {2, {0, 0, 0, 50, 0, 0}},
        {3, {10, 20, 30, 0, 0, 0}},
        {4, {1, 2, 3, 4, 5, 6}}};
    Temp t(FrdText()
               .Nodes(numbered_nodes(4))
               .Elements({{3, iota(4)}})
               .Result(1.0, 0, 1, "STRESS", {"SXX", "SYY", "SZZ", "SXY", "SYZ", "SZX"}, rows)
               .Result(1.0, 0, 1, "DISP", {"D1", "D2", "D3"},
                       {{1, {1, 2, 3}}, {2, {1, 2, 3}}, {3, {1, 2, 3}}, {4, {1, 2, 3}}})
               .Finish());
    FrdReadOptions frd;
    frd.mDerived = true;
    const Mesh mesh = meshioplusplus::read_frd(t.mPath, ReadOptions{}, frd);
    ASSERT_TRUE(mesh.HasPointData("STRESS_mises"));
    ASSERT_TRUE(mesh.HasPointData("STRESS_principal"));
    EXPECT_FALSE(mesh.HasPointData("DISP_mises"));
    EXPECT_NEAR(at(mesh.PointData("STRESS_mises"), 0), 100.0, 1e-12);
    EXPECT_NEAR(at(mesh.PointData("STRESS_mises"), 1), 50.0 * std::sqrt(3.0), 1e-12);
    const auto& principal = mesh.PointData("STRESS_principal");
    EXPECT_EQ(principal.Shape(), (std::vector<std::size_t>{4, 3}));
    EXPECT_NEAR(at(principal, 0), 0.0, 1e-12);
    EXPECT_NEAR(at(principal, 2), 100.0, 1e-12);
    EXPECT_NEAR(at(principal, 3), -50.0, 1e-12);  // shear: -50, 0, 50
    EXPECT_NEAR(at(principal, 5), 50.0, 1e-12);
    for (std::size_t r = 0; r < 4; ++r) {
        EXPECT_LE(at(principal, r * 3), at(principal, r * 3 + 1));
        EXPECT_LE(at(principal, r * 3 + 1), at(principal, r * 3 + 2));
    }
    // off by default
    EXPECT_FALSE(meshioplusplus::read_frd(t.mPath).HasPointData("STRESS_mises"));
}

TEST(FrdRead, GroupAndMaterialAreCellData) {
    Temp t(FrdText().Nodes(numbered_nodes(4)).Elements({{3, iota(4)}}).Finish());
    const Mesh mesh = meshioplusplus::read_frd(t.mPath);
    EXPECT_EQ(meshioplusplus::detail::read_int(mesh.CellData("frd:group", 0), 0), 0);
    EXPECT_EQ(meshioplusplus::detail::read_int(mesh.CellData("frd:material", 0), 0), 1);
}

TEST(FrdRead, ErrorsAreReadErrors) {
    const std::string good = FrdText().Nodes(numbered_nodes(4)).Elements({{3, iota(4)}}).Finish();
    // A binary flag on an otherwise-ASCII node header is covered by its own
    // BinaryClaimOnAsciiTextFailsCleanly test.
    {
        Temp t("hello\nworld\n");
        EXPECT_THROW(meshioplusplus::read_frd(t.mPath), ReadError);
    }
    EXPECT_THROW(meshioplusplus::read_frd("/no/such/file.frd"), ReadError);
    {
        Temp t(FrdText().Nodes(numbered_nodes(4)).Elements({{3, {1, 2, 3, 44}}}).Finish());
        EXPECT_THROW(meshioplusplus::read_frd(t.mPath), ReadError);
    }
    {
        Temp t(FrdText().Nodes(numbered_nodes(4)).Elements({{3, {1, 2, 3}}}).Finish());
        EXPECT_THROW(meshioplusplus::read_frd(t.mPath), ReadError);  // tetra needs four
    }
    {
        Temp t(FrdText()
                   .Nodes(numbered_nodes(4))
                   .Elements({{3, iota(4)}})
                   .Result(1.0, 0, 1, "NDTEMP", {"T"}, {{9, {1.0}}})
                   .Finish());
        EXPECT_THROW(meshioplusplus::read_frd(t.mPath), ReadError);  // result on node 9
    }
}

TEST(FrdRegistry, ReadsAndWritesWithAnExtensionAndMetadata) {
    EXPECT_EQ(meshioplusplus::resolve_format("results.frd", ""), "frd");
    EXPECT_EQ(meshioplusplus::registry_readers().count("frd"), 1u);
    EXPECT_EQ(meshioplusplus::registry_writers().count("frd"), 1u);
    EXPECT_EQ(meshioplusplus::registry_metadata_readers().count("frd"), 1u);
    EXPECT_TRUE(meshioplusplus::seq_format_may_have_steps("frd"));
}

// --- Writing -----------------------------------------------------------------------

namespace {

using meshioplusplus::DType;
using meshioplusplus::NDArray;
using meshioplusplus::WriteError;

NDArray doubles(std::initializer_list<std::size_t> shape, const std::vector<double>& values) {
    NDArray out(DType::Float64, std::vector<std::size_t>(shape));
    std::copy(values.begin(), values.end(), out.As<double>());
    return out;
}

NDArray ints(std::initializer_list<std::size_t> shape, const std::vector<std::int64_t>& values) {
    NDArray out(DType::Int64, std::vector<std::size_t>(shape));
    std::copy(values.begin(), values.end(), out.As<std::int64_t>());
    return out;
}

/// A tetrahedron, a quadratic wedge and a quadratic line (their node orders differ
/// from the file's), with results and a step.
Mesh sample_mesh() {
    Mesh mesh;
    std::vector<double> xyz;
    for (std::size_t i = 0; i < 20; ++i) {
        xyz.push_back(static_cast<double>(i));
        xyz.push_back(0.5 * static_cast<double>(i * i));
        xyz.push_back(-1.25 + 0.001 * static_cast<double>(i));
    }
    mesh.AssignPoints(doubles({20, 3}, xyz));
    mesh.AddCellBlock("tetra", ints({1, 4}, {0, 1, 2, 3}));
    mesh.AddCellBlock("wedge15", ints({1, 15}, {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}));
    mesh.AddCellBlock("line3", ints({1, 3}, {2, 19, 3}));
    std::vector<double> disp, stress, temp;
    for (std::size_t i = 0; i < 20; ++i) {
        for (std::size_t k = 0; k < 3; ++k)
            disp.push_back(1.0e-3 * static_cast<double>(i + k));
        for (std::size_t k = 0; k < 6; ++k)
            stress.push_back(1.5e6 * static_cast<double>(i) - 2.5e5 * static_cast<double>(k));
        temp.push_back(20.0 + 0.25 * static_cast<double>(i));
    }
    mesh.AddPointData("DISP", doubles({20, 3}, disp));
    mesh.AddPointData("STRESS", doubles({20, 6}, stress));
    mesh.AddPointData("NDTEMP", doubles({20}, temp));
    mesh.AddCellData("frd:group", {ints({1}, {2}), ints({1}, {3}), ints({1}, {4})});
    mesh.AddCellData("frd:material", {ints({1}, {1}), ints({1}, {5}), ints({1}, {6})});
    mesh.AddFieldData(meshioplusplus::kSequenceTimeKey, doubles({1}, {0.125}));
    mesh.AddFieldData("frd:step", ints({1}, {4}));
    mesh.AddFieldData("frd:analysis", ints({1}, {1}));
    return mesh;
}

struct TempFile {
    std::string mPath = mt::temp_path(".frd");
    ~TempFile() { std::remove(mPath.c_str()); }
};

std::string slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

void expect_close(const NDArray& rA, const NDArray& rB, const std::string& rWhat) {
    ASSERT_EQ(rA.Size(), rB.Size()) << rWhat;
    for (std::size_t i = 0; i < rA.Size(); ++i) {
        const double a = at(rA, i);
        const double b = at(rB, i);
        EXPECT_NEAR(a, b, 1e-5 * std::fabs(b) + 1e-12) << rWhat << " [" << i << "]";
    }
}

}  // namespace

TEST(FrdWrite, RoundTripKeepsMeshDataAndStep) {
    const Mesh mesh = sample_mesh();
    for (const bool long_ids : {true, false}) {
        TempFile file;
        meshioplusplus::write_frd(file.mPath, mesh, long_ids);
        const Mesh back = meshioplusplus::read_frd(file.mPath);
        ASSERT_EQ(back.NumCellBlocks(), 3u);
        EXPECT_EQ(back.Cells(0).Type(), "tetra");
        EXPECT_EQ(back.Cells(1).Type(), "wedge15");
        EXPECT_EQ(back.Cells(2).Type(), "line3");
        EXPECT_EQ(row_of(back, 1, 0), row_of(mesh, 1, 0));  // the node permutation round-trips
        EXPECT_EQ(row_of(back, 2, 0), (V{2, 19, 3}));
        expect_close(back.Points(), mesh.Points(), "points");
        for (const char* name : {"DISP", "STRESS", "NDTEMP"})
            expect_close(back.PointData(name), mesh.PointData(name), name);
        EXPECT_EQ(at(back.CellData("frd:group", 1), 0), 3.0);
        EXPECT_EQ(at(back.CellData("frd:material", 2), 0), 6.0);
        EXPECT_EQ(at(back.FieldData(meshioplusplus::kSequenceTimeKey), 0), 0.125);
        EXPECT_EQ(at(back.FieldData("frd:step"), 0), 4.0);
        EXPECT_EQ(at(back.FieldData("frd:analysis"), 0), 1.0);
    }
}

TEST(FrdWrite, MatchesTheHandBuiltFileTheReaderTestsUse) {
    // The scratch builder is the oracle: the library writer must read back the same.
    Nodes nodes = numbered_nodes(20);
    V conn;
    for (std::int64_t i = 1; i <= 15; ++i)
        conn.push_back(i);
    std::vector<std::pair<std::int64_t, std::vector<double>>> rows;
    for (std::int64_t i = 1; i <= 20; ++i)
        rows.push_back({i, {0.5 * static_cast<double>(i), -2.0, 1.0e3}});
    Temp oracle(FrdText()
                    .Nodes(nodes)
                    .Elements({{5, conn}})
                    .Result(0.5, 0, 1, "DISP", {"D1", "D2", "D3"}, rows, true)
                    .Finish());
    const Mesh expected = meshioplusplus::read_frd(oracle.mPath);

    TempFile file;
    meshioplusplus::write_frd(file.mPath, expected);
    const Mesh actual = meshioplusplus::read_frd(file.mPath);
    ASSERT_EQ(actual.NumCellBlocks(), 1u);
    EXPECT_EQ(row_of(actual, 0, 0), row_of(expected, 0, 0));
    expect_close(actual.Points(), expected.Points(), "points");
    expect_close(actual.PointData("DISP"), expected.PointData("DISP"), "DISP");
}

TEST(FrdWrite, TheLayoutIsTheOneCcxWrites) {
    TempFile file;
    meshioplusplus::write_frd(file.mPath, sample_mesh());
    const std::string text = slurp(file.mPath);
    EXPECT_EQ(text.rfind("    1C\n    1UWritten by meshio++", 0), 0u);
    EXPECT_NE(text.find("    2C                            20" + std::string(37, ' ') + "1\n"),
              std::string::npos);
    EXPECT_NE(text.find(" -4  DISP        4    1\n"), std::string::npos);  // DISP + calculated ALL
    EXPECT_NE(text.find(" -5  ALL         1    2    0    0    1ALL\n"), std::string::npos);
    EXPECT_NE(text.find(" -4  STRESS      6    1\n -5  SXX         1    4    1    1\n"),
              std::string::npos);
    EXPECT_EQ(text.substr(text.size() - 11), " -3\n  9999\n");
}

TEST(FrdWrite, TheShortFormUsesI5Ids) {
    TempFile file;
    meshioplusplus::write_frd(file.mPath, sample_mesh(), false);
    const std::string text = slurp(file.mPath);
    EXPECT_NE(text.find(std::string(37, ' ') + "0\n -1    1 "), std::string::npos);
}

TEST(FrdWrite, DropsWhatTheFormatCannotHold) {
    Mesh mesh = sample_mesh();
    mesh.AddCellBlock("pyramid", ints({1, 5}, {0, 1, 2, 3, 4}));
    mesh.AddCellData("frd:group", {ints({1}, {2}), ints({1}, {3}), ints({1}, {4}), ints({1}, {5})});
    mesh.AddCellData("pressure", {doubles({1}, {1}), doubles({1}, {2}), doubles({1}, {3}),
                                  doubles({1}, {4})});
    mesh.AddPointData("TOOLONGNAME", doubles({20}, std::vector<double>(20, 1.0)));
    mesh.AddPointData("per_cell", doubles({4}, {1, 2, 3, 4}));
    mesh.AddFieldData("note", ints({1}, {1}));
    TempFile file;
    meshioplusplus::write_frd(file.mPath, mesh);
    const Mesh back = meshioplusplus::read_frd(file.mPath);
    EXPECT_EQ(back.NumCellBlocks(), 3u);  // the pyramid is gone
    EXPECT_FALSE(back.HasCellData("pressure"));
    EXPECT_FALSE(back.HasFieldData("note"));
    EXPECT_FALSE(back.HasPointData("TOOLONGNAME"));
    EXPECT_FALSE(back.HasPointData("per_cell"));
    EXPECT_TRUE(back.HasPointData("DISP"));
}

TEST(FrdWrite, LeavesOutNodesWithANonFiniteValue) {
    Mesh mesh = sample_mesh();
    std::vector<double> temp(20, 1.0);
    temp[3] = std::nan("");
    mesh.AddPointData("NDTEMP", doubles({20}, temp));
    TempFile file;
    meshioplusplus::write_frd(file.mPath, mesh);
    const Mesh back = meshioplusplus::read_frd(file.mPath);
    EXPECT_TRUE(std::isnan(at(back.PointData("NDTEMP"), 3)));
    EXPECT_EQ(at(back.PointData("NDTEMP"), 4), 1.0);
}

TEST(FrdWrite, RefusesWhatItCannotRepresentAndLeavesNoFile) {
    {
        Mesh mesh = sample_mesh();
        std::vector<double> temp(20, 0.0);
        temp[0] = -1.0e120;
        mesh.AddPointData("NDTEMP", doubles({20}, temp));
        TempFile file;
        EXPECT_THROW(meshioplusplus::write_frd(file.mPath, mesh), WriteError);
        std::ifstream in(file.mPath);
        EXPECT_FALSE(in.good());
        temp[0] = -1.0e-120;  // tiny values flush to zero instead
        mesh.AddPointData("NDTEMP", doubles({20}, temp));
        meshioplusplus::write_frd(file.mPath, mesh);
        EXPECT_EQ(at(meshioplusplus::read_frd(file.mPath).PointData("NDTEMP"), 0), 0.0);
    }
    {
        Mesh mesh;
        mesh.AssignPoints(doubles({1, 3}, {INFINITY, 0.0, 0.0}));
        TempFile file;
        EXPECT_THROW(meshioplusplus::write_frd(file.mPath, mesh), WriteError);
    }
    {
        Mesh mesh;
        mesh.AssignPoints(doubles({3, 3}, {0, 0, 0, 1, 0, 0, 0, 1, 0}));
        mesh.AddCellBlock("triangle", ints({1, 3}, {0, 1, 5}));
        TempFile file;
        EXPECT_THROW(meshioplusplus::write_frd(file.mPath, mesh), WriteError);
    }
    {
        Mesh mesh = sample_mesh();
        mesh.AddCellData("frd:group", {ints({1}, {100000}), ints({1}, {3}), ints({1}, {4})});
        TempFile file;
        EXPECT_THROW(meshioplusplus::write_frd(file.mPath, mesh), WriteError);
    }
    {
        Mesh mesh;
        mesh.AssignPoints(NDArray(DType::Float64, {std::size_t{100000}, std::size_t{3}}));
        TempFile file;
        EXPECT_THROW(meshioplusplus::write_frd(file.mPath, mesh, false), WriteError);
        meshioplusplus::write_frd(file.mPath, mesh);  // the long form holds it
        EXPECT_EQ(meshioplusplus::read_frd(file.mPath).NumPoints(), 100000u);
    }
}

TEST(FrdWrite, TheRegistryWritesIt) {
    const auto& writers = meshioplusplus::registry_writers();
    ASSERT_EQ(writers.count("frd"), 1u);
    TempFile file;
    writers.at("frd")(file.mPath, sample_mesh());
    EXPECT_EQ(meshioplusplus::read_frd(file.mPath).NumCellBlocks(), 3u);
}

// --- The series writer ----------------------------------------------------------------

namespace {

/// `sample_mesh()` with its step fields set.
Mesh step_mesh(double time, std::int64_t step, double scale) {
    Mesh mesh = sample_mesh();
    std::vector<double> temp;
    for (std::size_t i = 0; i < 20; ++i)
        temp.push_back(scale * (20.0 + 0.25 * static_cast<double>(i)));
    mesh.AddPointData("NDTEMP", doubles({20}, temp));
    mesh.AddFieldData(meshioplusplus::kSequenceTimeKey, doubles({1}, {time}));
    mesh.AddFieldData("frd:step", ints({1}, {step}));
    return mesh;
}

}  // namespace

TEST(FrdSeries, EveryStepReadsBackByItsIndex) {
    for (const bool long_ids : {true, false}) {
        TempFile file;
        {
            meshioplusplus::FrdSeriesWriter writer(file.mPath, long_ids);
            for (int k = 0; k < 3; ++k)
                writer.Write(0.5 * (k + 1), step_mesh(0.0, k + 1, 1.0 + k));
            EXPECT_EQ(writer.NumSteps(), 3u);
            writer.Finalize();
            writer.Finalize();  // idempotent
        }
        EXPECT_EQ(meshioplusplus::read_frd_metadata(file.mPath).mTimeValues,
                  (std::vector<double>{0.5, 1.0, 1.5}));
        for (int k = 0; k < 3; ++k) {
            meshioplusplus::ReadOptions opts;
            opts.mTimeStep = k;
            const Mesh back = meshioplusplus::read_frd(file.mPath, opts);
            EXPECT_EQ(back.NumCellBlocks(), 3u);
            EXPECT_EQ(at(back.FieldData(meshioplusplus::kSequenceTimeKey), 0), 0.5 * (k + 1));
            EXPECT_EQ(at(back.FieldData("frd:step"), 0), static_cast<double>(k + 1));
            EXPECT_NEAR(at(back.PointData("NDTEMP"), 4), (1.0 + k) * 21.0, 1e-4 * (1.0 + k));
        }
    }
}

TEST(FrdSeries, ASeriesWithOneStepIsTheSingleStepFile) {
    const Mesh mesh = step_mesh(0.25, 1, 1.0);
    TempFile series, single;
    {
        meshioplusplus::FrdSeriesWriter writer(series.mPath);
        writer.Write(0.25, mesh);
        writer.Finalize();
    }
    meshioplusplus::write_frd(single.mPath, mesh);
    EXPECT_EQ(slurp(series.mPath), slurp(single.mPath));
}

TEST(FrdSeries, ARepeatedStepNumberGetsTheNextFreeOne) {
    TempFile file;
    {
        meshioplusplus::FrdSeriesWriter writer(file.mPath);
        for (int k = 0; k < 3; ++k)
            writer.Write(1.0 * k, step_mesh(0.0, 5, 1.0));
        writer.Finalize();
    }
    std::vector<double> steps;
    for (int k = 0; k < 3; ++k) {
        meshioplusplus::ReadOptions opts;
        opts.mTimeStep = k;
        steps.push_back(at(meshioplusplus::read_frd(file.mPath, opts).FieldData("frd:step"), 0));
    }
    EXPECT_EQ(steps, (std::vector<double>{5, 1, 2}));
}

TEST(FrdSeries, RefusesOtherCellsAndAnEmptySeries) {
    TempFile file;
    meshioplusplus::FrdSeriesWriter writer(file.mPath);
    EXPECT_THROW(writer.Finalize(), WriteError);  // nothing written yet
    writer.Write(0.0, step_mesh(0.0, 1, 1.0));
    Mesh other;
    other.AssignPoints(doubles({20, 3}, std::vector<double>(60, 0.0)));
    other.AddCellBlock("tetra", ints({1, 4}, {0, 1, 2, 4}));
    EXPECT_THROW(writer.Write(1.0, other), WriteError);
    EXPECT_EQ(writer.NumSteps(), 1u);  // the refused step left the file intact
    writer.Finalize();
    EXPECT_EQ(meshioplusplus::read_frd_metadata(file.mPath).mTimeValues.size(), 1u);
}

TEST(FrdSeries, TheSequenceEngineWritesItAndSaysSo) {
    std::string why;
    EXPECT_TRUE(meshioplusplus::sequence_write_supports_time("frd", why));
    EXPECT_TRUE(why.empty());
}
