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

// External includes
#include <gtest/gtest.h>

// System includes
#include <filesystem>
#include <fstream>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/pcd.hpp"

namespace {

// A 4-point cloud with normals, byte colours and one scalar, as the pcd
// readers and subsample_points produce it: one `vertex` block.
meshioplusplus::Mesh pcd_cloud() {
    using meshioplusplus::DType;
    using meshioplusplus::NDArray;
    meshioplusplus::Mesh m;
    NDArray pts(DType::Float32, {4, 3});
    const float xyz[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            pts.As<float>()[3 * i + j] = xyz[i][j];
    m.AssignPoints(std::move(pts));
    NDArray conn(DType::Int64, {4, 1});
    for (std::size_t i = 0; i < 4; ++i)
        conn.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    m.AddCellBlock("vertex", std::move(conn));
    NDArray normals(DType::Float32, {4, 3});
    for (std::size_t i = 0; i < 4; ++i) {
        normals.As<float>()[3 * i + 0] = 0.0f;
        normals.As<float>()[3 * i + 1] = 0.0f;
        normals.As<float>()[3 * i + 2] = 1.0f;
    }
    m.AddPointData("normals", std::move(normals));
    NDArray rgb(DType::UInt8, {4, 3});
    rgb.As<std::uint8_t>()[0] = 255;
    rgb.As<std::uint8_t>()[1] = 0;
    rgb.As<std::uint8_t>()[2] = 0;
    for (std::size_t i = 3; i < 12; ++i)
        rgb.As<std::uint8_t>()[i] = static_cast<std::uint8_t>(i);
    m.AddPointData("rgb", std::move(rgb));
    NDArray intensity(DType::Float32, {4});
    for (std::size_t i = 0; i < 4; ++i)
        intensity.As<float>()[i] = 0.25f * static_cast<float>(i);
    m.AddPointData("intensity", std::move(intensity));
    return m;
}

void pcd_write_text(const std::string& path, const std::string& body) {
    std::ofstream out(path, std::ios::binary);
    out << body;
}

}  // namespace

TEST(Pcd, EveryDataModeRoundTrips) {
    for (const meshioplusplus::PcdData data :
         {meshioplusplus::PcdData::Ascii, meshioplusplus::PcdData::Binary,
          meshioplusplus::PcdData::BinaryCompressed}) {
        const std::string path = mt::temp_path(".pcd");
        meshioplusplus::write_pcd(path, pcd_cloud(), data);
        const meshioplusplus::Mesh back = meshioplusplus::read_pcd(path);
        EXPECT_EQ(back.NumPoints(), 4u);
        ASSERT_EQ(back.NumCellBlocks(), 1u);
        EXPECT_EQ(back.Cells(0).Type(), "vertex");
        EXPECT_TRUE(back.HasPointData("normals"));
        EXPECT_TRUE(back.HasPointData("rgb"));
        EXPECT_TRUE(back.HasPointData("intensity"));
        EXPECT_EQ(back.PointData("rgb").Dtype(), meshioplusplus::DType::UInt8);
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
}

TEST(Pcd, FloatRgbSlotIsUnpackedByBitCast) {
    // PCL's `rgb` is a uint32 0x00RRGGBB living in a float32 slot: 0x00FF0000
    // as a float reads back as red, which a value-cast would turn to black.
    const std::string path = mt::temp_path(".pcd");
    pcd_write_text(path,
                   "# .PCD v0.7 - Point Cloud Data file format\n"
                   "VERSION 0.7\n"
                   "FIELDS x y z rgb\n"
                   "SIZE 4 4 4 4\n"
                   "TYPE F F F F\n"
                   "COUNT 1 1 1 1\n"
                   "WIDTH 1\n"
                   "HEIGHT 1\n"
                   "VIEWPOINT 0 0 0 1 0 0 0\n"
                   "POINTS 1\n"
                   "DATA ascii\n"
                   "0 0 0 2.34180515e-38\n");
    const meshioplusplus::Mesh back = meshioplusplus::read_pcd(path);
    ASSERT_TRUE(back.HasPointData("rgb"));
    const meshioplusplus::NDArray& rgb = back.PointData("rgb");
    EXPECT_EQ(rgb.As<std::uint8_t>()[0], 255);
    EXPECT_EQ(rgb.As<std::uint8_t>()[1], 0);
    EXPECT_EQ(rgb.As<std::uint8_t>()[2], 0);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Pcd, OrganisedCloudsAreKeptWholeUnlessDropped) {
    const std::string path = mt::temp_path(".pcd");
    pcd_write_text(path,
                   "VERSION 0.7\n"
                   "FIELDS x y z\n"
                   "SIZE 4 4 4\n"
                   "TYPE F F F\n"
                   "COUNT 1 1 1\n"
                   "WIDTH 3\n"
                   "HEIGHT 2\n"
                   "VIEWPOINT 0 0 0 1 0 0 0\n"
                   "POINTS 6\n"
                   "DATA ascii\n"
                   "0 0 1\n1 0 1\nnan nan nan\n1 1 1\n0 1 1\nnan nan nan\n");
    const meshioplusplus::Mesh kept = meshioplusplus::read_pcd(path);
    EXPECT_EQ(kept.NumPoints(), 6u);
    EXPECT_TRUE(kept.HasFieldData("pcd:width"));
    meshioplusplus::PcdReadOptions drop;
    drop.mDropInvalid = true;
    const meshioplusplus::Mesh clean = meshioplusplus::read_pcd(path, drop);
    EXPECT_EQ(clean.NumPoints(), 4u);
    EXPECT_FALSE(clean.HasFieldData("pcd:width"));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Pcd, ViewpointIsRecordedNeverApplied) {
    const std::string path = mt::temp_path(".pcd");
    pcd_write_text(path,
                   "VERSION 0.7\n"
                   "FIELDS x y z\n"
                   "SIZE 4 4 4\n"
                   "TYPE F F F\n"
                   "COUNT 1 1 1\n"
                   "WIDTH 1\n"
                   "HEIGHT 1\n"
                   "VIEWPOINT 1 2 3 0 1 0 0\n"
                   "POINTS 1\n"
                   "DATA ascii\n"
                   "1 2 3\n");
    const meshioplusplus::Mesh back = meshioplusplus::read_pcd(path);
    EXPECT_EQ(back.Points().As<float>()[0], 1.0f);
    EXPECT_TRUE(back.HasFieldData("pcd:viewpoint"));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Pcd, MalformedHeadersAreRefused) {
    const std::string path = mt::temp_path(".pcd");
    pcd_write_text(path, "VERSION 0.7\nPOINTS 1\nDATA ascii\n1 2 3\n");
    EXPECT_THROW(meshioplusplus::read_pcd(path), meshioplusplus::ReadError);
    pcd_write_text(path,
                   "VERSION 0.7\n"
                   "FIELDS x y z\n"
                   "SIZE 4 4 4\n"
                   "TYPE F F F\n"
                   "COUNT 1 1 1\n"
                   "WIDTH 1\n"
                   "HEIGHT 1\n"
                   "VIEWPOINT 0 0 0 1 0 0 0\n"
                   "POINTS 1\n"
                   "DATA zip\n");
    EXPECT_THROW(meshioplusplus::read_pcd(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
