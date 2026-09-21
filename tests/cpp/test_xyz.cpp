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
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/xyz.hpp"

namespace {

// A 4-point cloud with normals, byte colours and one scalar, as the xyz
// readers and subsample_points produce it: one `vertex` block.
meshioplusplus::Mesh xyz_cloud() {
    using meshioplusplus::DType;
    using meshioplusplus::NDArray;
    meshioplusplus::Mesh m;
    NDArray pts(DType::Float64, {4, 3});
    const double xyz[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            pts.As<double>()[3 * i + j] = xyz[i][j];
    m.AssignPoints(std::move(pts));
    NDArray conn(DType::Int64, {4, 1});
    for (std::size_t i = 0; i < 4; ++i)
        conn.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    m.AddCellBlock("vertex", std::move(conn));
    NDArray normals(DType::Float64, {4, 3});
    for (std::size_t i = 0; i < 4; ++i) {
        normals.As<double>()[3 * i + 0] = 0.0;
        normals.As<double>()[3 * i + 1] = 0.0;
        normals.As<double>()[3 * i + 2] = 1.0;
    }
    m.AddPointData("normals", std::move(normals));
    NDArray rgb(meshioplusplus::DType::UInt8, {4, 3});
    for (std::size_t i = 0; i < 12; ++i)
        rgb.As<std::uint8_t>()[i] = static_cast<std::uint8_t>(10 * i);
    m.AddPointData("rgb", std::move(rgb));
    NDArray scalar(DType::Float64, {4});
    for (std::size_t i = 0; i < 4; ++i)
        scalar.As<double>()[i] = 0.5 * static_cast<double>(i);
    m.AddPointData("temperature", std::move(scalar));
    return m;
}

void xyz_write_text(const std::string& path, const std::string& body) {
    std::ofstream out(path, std::ios::binary);
    out << body;
}

}  // namespace

TEST(Xyz, RoundTripsPointsNormalsColoursAndScalars) {
    const std::string path = mt::temp_path(".xyz");
    meshioplusplus::write_xyz(path, xyz_cloud());
    const meshioplusplus::Mesh back = meshioplusplus::read_xyz(path);
    EXPECT_EQ(back.NumPoints(), 4u);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "vertex");
    EXPECT_TRUE(back.HasPointData("normals"));
    EXPECT_TRUE(back.HasPointData("rgb"));
    EXPECT_TRUE(back.HasPointData("temperature"));
#if defined(MESHIOPLUSPLUS_MESH_BACKEND_MESHIO)
    EXPECT_EQ(back.PointData("rgb").Dtype(), meshioplusplus::DType::UInt8);
#else
    // NATIVE and KRATOS widen every integer array to Int64 on ingest: assert the kind.
    EXPECT_FALSE(meshioplusplus::detail::is_float_dtype(back.PointData("rgb").Dtype()));
#endif
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Xyz, SixColumnsResolveByRangeAndByExtension) {
    const std::string normals = mt::temp_path("_n.xyz");
    xyz_write_text(normals, "0 0 0 0 0 1\n");
    EXPECT_TRUE(meshioplusplus::read_xyz(normals).HasPointData("normals"));
    const std::string colours = mt::temp_path("_c.xyz");
    xyz_write_text(colours, "0 0 0 255 0 0\n");
    EXPECT_TRUE(meshioplusplus::read_xyz(colours).HasPointData("rgb"));
    const std::string alias = mt::temp_path("_a.xyzn");
    xyz_write_text(alias, "0 0 0 0.5 0.5 0.5\n");
    EXPECT_TRUE(meshioplusplus::read_xyz(alias).HasPointData("normals"));
    std::error_code ec;
    std::filesystem::remove(normals, ec);
    std::filesystem::remove(colours, ec);
    std::filesystem::remove(alias, ec);
}

TEST(Xyz, ExplicitColumnsAndDelimiter) {
    const std::string path = mt::temp_path(".xyz");
    xyz_write_text(path, "1,2,3,9,0.5\n");
    meshioplusplus::XyzReadOptions opts;
    opts.mColumns = {"x", "y", "z", "_", "p"};
    opts.mDelimiter = ",";
    const meshioplusplus::Mesh back = meshioplusplus::read_xyz(path, opts);
    EXPECT_FALSE(back.HasPointData("_"));
    EXPECT_TRUE(back.HasPointData("p"));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Xyz, ChemistryXyzIsRefusedByName) {
    const std::string path = mt::temp_path(".xyz");
    xyz_write_text(path, "2\ncomment\nC 0 0 0\nH 1 0 0\n");
    EXPECT_THROW(meshioplusplus::read_xyz(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Xyz, PtsCountIsValidated) {
    const std::string path = mt::temp_path(".pts");
    xyz_write_text(path, "3\n0 0 0\n");
    EXPECT_THROW(meshioplusplus::read_xyz(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Xyz, RaggedRowsAreRefused) {
    const std::string path = mt::temp_path(".xyz");
    xyz_write_text(path, "1 2 3\n1 2\n");
    EXPECT_THROW(meshioplusplus::read_xyz(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
