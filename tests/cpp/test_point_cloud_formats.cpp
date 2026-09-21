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
 * @file test_point_cloud_formats.cpp
 * @brief The PCD (PCL v0.7) and headerless-XYZ point-cloud readers/writers.
 *
 * The files are written by hand rather than read from tests/python/meshes/: this suite
 * has no test-data path, and a format's own writer is no oracle for its reader. The PCL
 * corpus (`binary_compressed` included) is exercised from the Python suite.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/pcd.hpp"
#include "meshioplusplus/formats/xyz.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/write_options.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::NDArray;
using meshioplusplus::PcdData;
using meshioplusplus::PcdReadOptions;
using meshioplusplus::XyzReadOptions;

std::string write_text(const std::string& rContents, const std::string& rSuffix) {
    const std::string path = mt::temp_path(rSuffix);
    std::ofstream out(path, std::ios::binary);
    out.write(rContents.data(), static_cast<std::streamsize>(rContents.size()));
    return path;
}

std::string slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::vector<double> values(const NDArray& rA) {
    std::vector<double> out;
    for (std::size_t i = 0; i < rA.Size(); ++i)
        out.push_back(meshioplusplus::detail::read_double(rA, i));
    return out;
}

void expect_values(const NDArray& rA, const std::vector<double>& rExpected) {
    const std::vector<double> got = values(rA);
    ASSERT_EQ(got.size(), rExpected.size());
    for (std::size_t i = 0; i < got.size(); ++i)
        EXPECT_DOUBLE_EQ(got[i], rExpected[i]) << "element " << i;
}

std::size_t num_blocks(const mt::Mesh& rMesh, const std::string& rType) {
    std::size_t n = 0;
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() == rType)
            n += cb.NumCells();
    return n;
}

// Five points with normals, uint8 colours, a scalar and an integer label.
mt::Mesh cloud() {
    mt::Mesh m;
    m.AssignPoints(
        mt::points_from({{0.1, 0.2, 0.3}, {1, 2, 3}, {-1.5, 0.25, 8}, {1e-3, 2e5, -7}, {0, 0, 0}}));
    m.AddCellBlock("vertex", mt::conn_from({{0}, {1}, {2}, {3}, {4}}));
    m.AddPointData("normals",
                   mt::data_array({0, 0, 1, 1, 0, 0, 0, 1, 0, 0.6, 0.8, 0, 0, 0, -1}, 3));
    NDArray rgb(DType::UInt8, {5, 3});
    const std::uint8_t colours[15] = {255, 0, 0, 0, 255, 0, 0, 0, 255, 1, 2, 3, 128, 64, 32};
    std::memcpy(rgb.Data(), colours, sizeof(colours));
    m.AddPointData("rgb", std::move(rgb));
    m.AddPointData("intensity", mt::data_array({0.5, 1.5, 2.5, 3.5, 4.5}));
    m.AddPointData("label", mt::int_data_array({1, 2, 3, 4, 5}));
    return m;
}

const char* const kHeaderXyz =
    "# .PCD v0.7 - Point Cloud Data file format\nVERSION 0.7\nFIELDS x y z\nSIZE 4 4 4\n"
    "TYPE F F F\nCOUNT 1 1 1\n";

}  // namespace

// ---------------------------------------------------------------------------
// PCD

TEST(PcdIo, RoundTripsEveryEncoding) {
    const mt::Mesh in = cloud();
    for (PcdData mode : {PcdData::Ascii, PcdData::Binary, PcdData::BinaryCompressed}) {
        const std::string path = mt::temp_path("_pcd_rt.pcd");
        meshioplusplus::write_pcd(path, in, mode, /*float64_points=*/true);
        const mt::Mesh out = meshioplusplus::read_pcd(path);
        expect_values(out.Points(), values(in.Points()));
        expect_values(out.PointData("normals"), values(in.PointData("normals")));
        expect_values(out.PointData("rgb"), values(in.PointData("rgb")));
        expect_values(out.PointData("intensity"), values(in.PointData("intensity")));
        expect_values(out.PointData("label"), values(in.PointData("label")));
        ASSERT_EQ(out.PointData("rgb").Shape().size(), 2u);
        EXPECT_EQ(out.PointData("rgb").Shape()[1], 3u);
        EXPECT_EQ(num_blocks(out, "vertex"), 5u);
        std::remove(path.c_str());
    }
}

TEST(PcdIo, WritesFloat32ByDefaultAndFloat64OnRequest) {
    const std::string a = mt::temp_path("_pcd_f32.pcd"), b = mt::temp_path("_pcd_f64.pcd");
    meshioplusplus::write_pcd(a, cloud(), PcdData::Binary);
    meshioplusplus::write_pcd(b, cloud(), PcdData::Binary, true);
    EXPECT_NE(slurp(a).find("SIZE 4 4 4"), std::string::npos);
    EXPECT_NE(slurp(b).find("SIZE 8 8 8"), std::string::npos);
    EXPECT_EQ(meshioplusplus::read_pcd(a).Points().Dtype(), DType::Float32);
    std::remove(a.c_str());
    std::remove(b.c_str());
}

TEST(PcdIo, ColoursAreBitCastFromAFloatSlot) {
    // 0x00FF8040 held in a float32: converting the value would give 0..255 garbage.
    const std::uint32_t packed = 0x00FF8040u;
    float slot;
    std::memcpy(&slot, &packed, 4);
    char text[64];
    std::snprintf(text, sizeof(text), "%.9g", static_cast<double>(slot));
    const std::string path = write_text(
        "VERSION 0.7\nFIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH 1\n"
        "HEIGHT 1\nPOINTS 1\nDATA ascii\n0 0 0 " +
            std::string(text) + "\n",
        ".pcd");
    const mt::Mesh m = meshioplusplus::read_pcd(path);
    expect_values(m.PointData("rgb"), {255, 128, 64});
    std::remove(path.c_str());
}

TEST(PcdIo, RgbaIsAUintAndRgbWritesAsAFloatSlot) {
    const std::string path = write_text(
        "VERSION 0.7\nFIELDS x y z rgba\nSIZE 4 4 4 4\nTYPE F F F U\nCOUNT 1 1 1 1\nWIDTH 1\n"
        "HEIGHT 1\nPOINTS 1\nDATA ascii\n0 0 0 2164228160\n",  // 0x80FF8040
        ".pcd");
    expect_values(meshioplusplus::read_pcd(path).PointData("rgba"), {255, 128, 64, 128});
    std::remove(path.c_str());

    const std::string out = mt::temp_path("_pcd_rgb.pcd");
    meshioplusplus::write_pcd(out, cloud(), PcdData::Ascii);
    const std::string text = slurp(out);
    EXPECT_NE(text.find("FIELDS x y z intensity label normal_x normal_y normal_z rgb"),
              std::string::npos);
    EXPECT_NE(text.find("TYPE F F F F I F F F F"), std::string::npos);
    std::remove(out.c_str());
}

TEST(PcdIo, MultiCountPaddingAndIntegerTypes) {
    const std::string path = write_text(
        "VERSION 0.7\nFIELDS x y z hist _ ring delta\nSIZE 4 4 4 4 1 2 4\n"
        "TYPE F F F F U U I\nCOUNT 1 1 1 3 1 1 1\nWIDTH 2\nHEIGHT 1\nPOINTS 2\nDATA ascii\n"
        "1 2 3 10 11 12 0 5 100\n4 5 6 13 14 15 0 6 -200\n",
        ".pcd");
    const mt::Mesh m = meshioplusplus::read_pcd(path);
    expect_values(m.PointData("hist"), {10, 11, 12, 13, 14, 15});
    EXPECT_EQ(m.PointData("hist").Shape()[1], 3u);
    expect_values(m.PointData("ring"), {5, 6});
    expect_values(m.PointData("delta"), {100, -200});
    EXPECT_FALSE(m.HasPointData("_"));
    std::remove(path.c_str());
}

TEST(PcdIo, OrganisedCloudsAreKeptOrDropped) {
    const std::string path = write_text(
        std::string(kHeaderXyz) +
            "WIDTH 3\nHEIGHT 2\nPOINTS 6\nDATA ascii\n0 0 1\n1 0 1\nnan nan nan\n1 1 1\n0 1 1\n"
            "nan nan nan\n",
        ".pcd");
    const mt::Mesh kept = meshioplusplus::read_pcd(path);
    EXPECT_EQ(kept.NumPoints(), 6u);
    expect_values(kept.FieldData("pcd:width"), {3});
    expect_values(kept.FieldData("pcd:height"), {2});

    PcdReadOptions options;
    options.mDropInvalid = true;
    const mt::Mesh dropped = meshioplusplus::read_pcd(path, options);
    EXPECT_EQ(dropped.NumPoints(), 4u);
    EXPECT_FALSE(dropped.HasFieldData("pcd:width"));
    EXPECT_EQ(num_blocks(dropped, "vertex"), 4u);

    // and the organisation survives a write
    const std::string out = mt::temp_path("_pcd_org.pcd");
    meshioplusplus::write_pcd(out, kept, PcdData::Binary);
    EXPECT_NE(slurp(out).find("WIDTH 3\nHEIGHT 2\n"), std::string::npos);
    std::remove(path.c_str());
    std::remove(out.c_str());
}

TEST(PcdIo, ViewpointIsRecordedNotApplied) {
    const std::string path =
        write_text(std::string(kHeaderXyz) +
                       "WIDTH 1\nHEIGHT 1\nVIEWPOINT 1 2 3 0 1 0 0\nPOINTS 1\nDATA ascii\n1 2 3\n",
                   ".pcd");
    const mt::Mesh m = meshioplusplus::read_pcd(path);
    expect_values(m.Points(), {1, 2, 3});
    expect_values(m.FieldData("pcd:viewpoint"), {1, 2, 3, 0, 1, 0, 0});
    const std::string out = mt::temp_path("_pcd_vp.pcd");
    meshioplusplus::write_pcd(out, m, PcdData::Ascii);
    EXPECT_NE(slurp(out).find("VIEWPOINT 1 2 3 0 1 0 0\n"), std::string::npos);
    std::remove(path.c_str());
    std::remove(out.c_str());
}

TEST(PcdIo, DecodesAHandEncodedLzfStream) {
    // x = 1 1, y = 2 2, z = 3 3 as float32, struct-of-arrays: each 4-byte literal is
    // followed by a back reference (length 4, distance 4) -- written from the stream
    // format by hand, not produced by this library's own encoder.
    const unsigned char stream[] = {3, 0x00, 0x00, 0x80, 0x3F, 0x40, 0x03,  //
                                    3, 0x00, 0x00, 0x00, 0x40, 0x40, 0x03,  //
                                    3, 0x00, 0x00, 0x40, 0x40, 0x40, 0x03};
    std::string file =
        std::string(kHeaderXyz) + "WIDTH 2\nHEIGHT 1\nPOINTS 2\nDATA binary_compressed\n";
    const std::uint32_t sizes[2] = {sizeof(stream), 24};
    file.append(reinterpret_cast<const char*>(sizes), 8);
    file.append(reinterpret_cast<const char*>(stream), sizeof(stream));
    const std::string path = write_text(file, ".pcd");
    expect_values(meshioplusplus::read_pcd(path).Points(), {1, 2, 3, 1, 2, 3});
    std::remove(path.c_str());
}

TEST(PcdIo, CorruptCompressedDataThrows) {
    const std::string good = mt::temp_path("_pcd_good.pcd");
    meshioplusplus::write_pcd(good, cloud(), PcdData::BinaryCompressed);
    const std::string bytes = slurp(good);
    const std::size_t at = bytes.find("DATA binary_compressed\n") + 23;
    for (const std::string& blob :
         {bytes.substr(0, at + 3), bytes.substr(0, at + 13),
          bytes.substr(0, at + 8) + std::string(bytes.size() - at - 8, '\xff')}) {
        const std::string bad = write_text(blob, ".pcd");
        EXPECT_THROW(meshioplusplus::read_pcd(bad), meshioplusplus::ReadError);
        std::remove(bad.c_str());
    }
    std::remove(good.c_str());
}

TEST(PcdIo, MalformedFilesThrowReadError) {
    const std::vector<std::string> bad = {
        "VERSION 0.7\nFIELDS x y z\n",  // no DATA
        std::string(kHeaderXyz) + "WIDTH 1\nHEIGHT 1\nPOINTS 1\nDATA zip\n1 2 3\n",
        "VERSION 0.7\nFIELDS x y z\nSIZE 4 4\nTYPE F F F\nCOUNT 1 1 1\nWIDTH 1\nPOINTS 1\nDATA "
        "ascii\n",
        "VERSION 0.7\nFIELDS x y z\nSIZE 4 4 4\nTYPE F F X\nCOUNT 1 1 1\nWIDTH 1\nPOINTS 1\nDATA "
        "ascii\n",
        "VERSION 0.7\nFIELDS x y w\nSIZE 4 4 4\nTYPE F F F\nCOUNT 1 1 1\nWIDTH 1\nPOINTS 1\nDATA "
        "ascii\n1 2 3\n",
        std::string(kHeaderXyz) + "WIDTH 1\nHEIGHT 1\nPOINTS 1\nDATA ascii\n1 2\n",
        std::string(kHeaderXyz) + "WIDTH 1\nHEIGHT 1\nPOINTS 1\nDATA ascii\n1 2 abc\n",
        std::string(kHeaderXyz) + "WIDTH 1\nHEIGHT 1\nPOINTS 1\nDATA binary\n\x01\x02",
        std::string(kHeaderXyz) + "WIDTH 1\nHEIGHT 1\nVIEWPOINT 1 2\nPOINTS 1\nDATA ascii\n1 2 3\n",
    };
    for (const std::string& text : bad) {
        const std::string path = write_text(text, ".pcd");
        EXPECT_THROW(meshioplusplus::read_pcd(path), meshioplusplus::ReadError) << text;
        std::remove(path.c_str());
    }
    EXPECT_THROW(meshioplusplus::read_pcd("/nonexistent/none.pcd"), meshioplusplus::ReadError);
}

TEST(PcdIo, WriterDropsWhatPcdCannotHold) {
    mt::Mesh m = mt::tri_mesh();  // a surface mesh: the cells go, the points stay
    const std::string path = mt::temp_path("_pcd_tri.pcd");
    meshioplusplus::write_pcd(path, m, PcdData::Binary);
    const mt::Mesh out = meshioplusplus::read_pcd(path);
    EXPECT_EQ(out.NumPoints(), m.NumPoints());
    EXPECT_EQ(num_blocks(out, "triangle"), 0u);
    EXPECT_EQ(num_blocks(out, "vertex"), m.NumPoints());
    std::remove(path.c_str());
}

TEST(PcdIo, EmptyCloudRoundTrips) {
    mt::Mesh m;
    m.AssignPoints(NDArray(DType::Float64, {0, 3}));
    for (PcdData mode : {PcdData::Ascii, PcdData::Binary, PcdData::BinaryCompressed}) {
        const std::string path = mt::temp_path("_pcd_empty.pcd");
        meshioplusplus::write_pcd(path, m, mode);
        EXPECT_EQ(meshioplusplus::read_pcd(path).NumPoints(), 0u);
        std::remove(path.c_str());
    }
}

TEST(PcdIo, IsSniffedAndRegistered) {
    const std::string path = mt::temp_path("_pcd_sniff.cloud");
    meshioplusplus::write_pcd(path, cloud(), PcdData::Binary);
    EXPECT_EQ(meshioplusplus::sniff_format(path), "pcd");
    std::remove(path.c_str());
    const std::string bare = write_text(std::string("VERSION .7\nFIELDS x y z\n"), ".cloud");
    EXPECT_EQ(meshioplusplus::sniff_format(bare), "pcd");
    std::remove(bare.c_str());
    EXPECT_EQ(meshioplusplus::resolve_format("a.pcd", ""), "pcd");
    EXPECT_TRUE(meshioplusplus::registry_readers().count("pcd"));
    EXPECT_TRUE(meshioplusplus::registry_writers().count("pcd"));
}

TEST(PcdIo, WriteOptionsSelectTheEncoding) {
    meshioplusplus::WriteOptions ascii;
    ascii.mEncoding = meshioplusplus::WriteEncoding::Ascii;
    const std::string a = mt::temp_path("_pcd_ex_a.pcd");
    meshioplusplus::registry_write_ex(a, cloud(), "pcd", ascii);
    EXPECT_NE(slurp(a).find("DATA ascii"), std::string::npos);
    meshioplusplus::WriteOptions binary;
    binary.mEncoding = meshioplusplus::WriteEncoding::Binary;
    meshioplusplus::registry_write_ex(a, cloud(), "pcd", binary);
    EXPECT_NE(slurp(a).find("DATA binary\n"), std::string::npos);
    std::remove(a.c_str());

    meshioplusplus::WriteOptions codec;
    codec.mCodecSet = true;
    std::string why;
    EXPECT_FALSE(meshioplusplus::registry_write_supports("pcd", codec, why));
}

// ---------------------------------------------------------------------------
// XYZ

TEST(XyzIo, RoundTripsExactly) {
    const mt::Mesh in = cloud();
    const std::string path = mt::temp_path("_xyz_rt.xyz");
    meshioplusplus::write_xyz(path, in);
    const mt::Mesh out = meshioplusplus::read_xyz(path);
    expect_values(out.Points(), values(in.Points()));
    expect_values(out.PointData("normals"), values(in.PointData("normals")));
    expect_values(out.PointData("rgb"), values(in.PointData("rgb")));
    expect_values(out.PointData("intensity"), values(in.PointData("intensity")));
    expect_values(out.PointData("label"), values(in.PointData("label")));
    EXPECT_EQ(num_blocks(out, "vertex"), 5u);
    std::remove(path.c_str());
}

TEST(XyzIo, HeaderRowsAndFloatFormat) {
    const std::string path = mt::temp_path("_xyz_w.xyz");
    meshioplusplus::write_xyz(path, cloud());
    std::istringstream text(slurp(path));
    std::string line, header;
    while (std::getline(text, line))
        if (line.rfind("# x y z", 0) == 0)
            header = line;
    EXPECT_EQ(header, "# x y z intensity label nx ny nz r g b");
    meshioplusplus::write_xyz(path, cloud(), ".3f");
    EXPECT_NE(slurp(path).find("0.100 0.200 0.300 0.500 1 0.000 0.000 1.000 255 0 0"),
              std::string::npos);
    EXPECT_THROW(meshioplusplus::write_xyz(path, cloud(), "%s%n"), meshioplusplus::WriteError);
    std::remove(path.c_str());
}

TEST(XyzIo, ColumnsAreResolvedByCountExtensionAndRange) {
    struct Case {
        const char* mSuffix;
        const char* mText;
        const char* mKey;  // the point data expected ("" = none)
        std::size_t mWidth;
    };
    const Case cases[] = {
        {".xyz", "1 2 3\n4 5 6\n", "", 0},
        {".xyz", "1 2 3 0.5\n4 5 6 0.7\n", "scalar", 1},
        {".pts", "2\n1 2 3 9\n4 5 6 8\n", "intensity", 1},
        {".xyzn", "1 2 3 0 0 1\n4 5 6 1 0 0\n", "normals", 3},
        {".xyzrgb", "1 2 3 255 0 0\n4 5 6 0 255 0\n", "rgb", 3},
        {".xyz", "1 2 3 0 0 1\n4 5 6 0 1 0\n", "normals", 3},  // unit vectors beat bytes
        {".xyz", "1 2 3 10 20 30\n4 5 6 40 50 60\n", "rgb", 3},
        {".xyz", "1 2 3 0.5 0.25 1\n4 5 6 0 0.1 0.2\n", "rgb", 3},
    };
    for (const Case& c : cases) {
        const std::string path = write_text(c.mText, c.mSuffix);
        const mt::Mesh m = meshioplusplus::read_xyz(path);
        EXPECT_EQ(m.NumPoints(), 2u) << c.mText;
        if (std::string(c.mKey).empty())
            EXPECT_EQ(m.NumPointData(), 0u) << c.mText;
        else {
            ASSERT_TRUE(m.HasPointData(c.mKey)) << c.mText;
            const auto& shape = m.PointData(c.mKey).Shape();
            EXPECT_EQ(shape.size() == 2 ? shape[1] : 1u, c.mWidth) << c.mText;
        }
        std::remove(path.c_str());
    }
    const std::string ambiguous = write_text("1 2 3 5.5 6.5 7.5\n4 5 6 4.5 5.5 6.5\n", ".xyz");
    EXPECT_THROW(meshioplusplus::read_xyz(ambiguous), meshioplusplus::ReadError);
    std::remove(ambiguous.c_str());
    const std::string five = write_text("1 2 3 4 5\n", ".xyz");
    EXPECT_THROW(meshioplusplus::read_xyz(five), meshioplusplus::ReadError);
    std::remove(five.c_str());
}

TEST(XyzIo, UnitFloatColoursAreScaled) {
    const std::string path = write_text("0 0 0 0.5 0.25 1\n0 0 0 0 0.1 0.2\n", ".xyz");
    expect_values(meshioplusplus::read_xyz(path).PointData("rgb"), {128, 64, 255, 0, 26, 51});
    std::remove(path.c_str());
}

TEST(XyzIo, DelimitersCommentsAndNumberForms) {
    for (const char* text : {"1 2 3\n4 5 6\n", "1,2,3\n4,5,6\n", "1;2;3\n4;5;6\n",
                             "1 , 2 , 3\n4, 5, 6\n", "  1\t2   3  \r\n4 5 6\r\n",
                             "# c\n\n1 2 3\n// again\n4 5 6\n# tail\n", "1,2,3,\n4,5,6,\n"}) {
        const std::string path = write_text(text, ".xyz");
        expect_values(meshioplusplus::read_xyz(path).Points(), {1, 2, 3, 4, 5, 6});
        std::remove(path.c_str());
    }
    XyzReadOptions options;
    options.mDelimiter = "|";
    const std::string bar = write_text("1|2|3\n4|5|6\n", ".xyz");
    expect_values(meshioplusplus::read_xyz(bar, options).Points(), {1, 2, 3, 4, 5, 6});
    std::remove(bar.c_str());

    const std::string forms = write_text("1e3 -2.5E-2 +3\nnan inf -inf\n", ".xyz");
    const std::vector<double> p = values(meshioplusplus::read_xyz(forms).Points());
    EXPECT_DOUBLE_EQ(p[0], 1000);
    EXPECT_DOUBLE_EQ(p[1], -0.025);
    EXPECT_TRUE(std::isnan(p[3]) && std::isinf(p[4]) && p[5] < 0);
    std::remove(forms.c_str());
}

TEST(XyzIo, ExplicitAndHeaderColumns) {
    XyzReadOptions options;
    options.mColumns = {"x", "y", "z", "_", "temperature", "u", "v"};
    const std::string path = write_text("1 2 3 4 5 0.1 0.2\n6 7 8 9 10 0.3 0.4\n", ".xyz");
    const mt::Mesh m = meshioplusplus::read_xyz(path, options);
    expect_values(m.PointData("temperature"), {5, 10});
    expect_values(m.PointData("u"), {0.1, 0.3});
    EXPECT_FALSE(m.HasPointData("_"));
    options.mColumns = {"x", "y", "y", "a", "b", "c", "d"};
    EXPECT_THROW(meshioplusplus::read_xyz(path, options), meshioplusplus::ReadError);
    std::remove(path.c_str());

    const std::string header =
        write_text("# made by hand\n# x y z Temperature nx ny nz\n1 2 3 9 0 0 1\n", ".xyz");
    const mt::Mesh h = meshioplusplus::read_xyz(header);
    expect_values(h.PointData("Temperature"), {9});
    expect_values(h.PointData("normals"), {0, 0, 1});
    std::remove(header.c_str());

    const std::string cc = write_text("//X,Y,Z,R,G,B,A\n1,2,3,255,0,10,128\n", ".asc");
    expect_values(meshioplusplus::read_xyz(cc).PointData("rgba"), {255, 0, 10, 128});
    std::remove(cc.c_str());
}

TEST(XyzIo, PtsCountIsValidated) {
    const std::string ok = write_text("2\n1 2 3\n4 5 6\n", ".pts");
    EXPECT_EQ(meshioplusplus::read_xyz(ok).NumPoints(), 2u);
    const std::string bad = write_text("3\n1 2 3\n4 5 6\n", ".pts");
    EXPECT_THROW(meshioplusplus::read_xyz(bad), meshioplusplus::ReadError);
    const std::string notpts = write_text("2\n1 2 3\n4 5 6\n", ".xyz");
    EXPECT_THROW(meshioplusplus::read_xyz(notpts), meshioplusplus::ReadError);
    for (const auto* p : {&ok, &bad, &notpts})
        std::remove(p->c_str());
}

TEST(XyzIo, ChemistryXyzIsRefusedByName) {
    const std::string water = write_text(
        "3\nwater molecule\nO  0.000  0.000  0.000\nH  0.757  0.586  0.000\nH -0.757  0.586  "
        "0.000\n",
        ".xyz");
    try {
        meshioplusplus::read_xyz(water);
        FAIL() << "chemistry XYZ was accepted";
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("chemistry"), std::string::npos);
    }
    std::remove(water.c_str());
    // a leading integer row is a count, not an atom count, when a number follows
    const std::string pts = write_text("2\n1 2 3\nnan 5 6\n", ".pts");
    EXPECT_EQ(meshioplusplus::read_xyz(pts).NumPoints(), 2u);
    std::remove(pts.c_str());
}

TEST(XyzIo, MalformedRowsNameTheLine) {
    const std::string ragged = write_text("1 2 3\n4 5 6\n7 8\n", ".xyz");
    try {
        meshioplusplus::read_xyz(ragged);
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("line 3"), std::string::npos);
    }
    const std::string text = write_text("1 2 3\n4 x 6\n", ".xyz");
    try {
        meshioplusplus::read_xyz(text);
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("not numeric"), std::string::npos);
    }
    std::remove(ragged.c_str());
    std::remove(text.c_str());
}

TEST(XyzIo, EmptyFilesAreEmptyClouds) {
    for (const char* text : {"", "\n\n", "# nothing here\n"}) {
        const std::string path = write_text(text, ".xyz");
        const mt::Mesh m = meshioplusplus::read_xyz(path);
        EXPECT_EQ(m.NumPoints(), 0u);
        std::remove(path.c_str());
    }
}

TEST(XyzIo, RegistryKnowsTheAliases) {
    for (const char* ext : {"a.xyz", "a.xyzn", "a.xyzrgb", "a.asc", "a.pts", "a.txt"})
        EXPECT_EQ(meshioplusplus::resolve_format(ext, ""), "xyz") << ext;
    EXPECT_TRUE(meshioplusplus::registry_readers().count("xyz"));
    meshioplusplus::WriteOptions ascii;
    ascii.mEncoding = meshioplusplus::WriteEncoding::Ascii;
    std::string why;
    EXPECT_FALSE(meshioplusplus::registry_write_supports("xyz", ascii, why));  // no variants
    meshioplusplus::WriteOptions ff;
    ff.mFloatFormat = ".3f";
    EXPECT_TRUE(meshioplusplus::registry_write_supports("xyz", ff, why));
    const std::string path = mt::temp_path("_xyz_ex.xyz");
    meshioplusplus::registry_write_ex(path, cloud(), "xyz", ff);
    EXPECT_NE(slurp(path).find("0.100 0.200 0.300"), std::string::npos);
    std::remove(path.c_str());
}
