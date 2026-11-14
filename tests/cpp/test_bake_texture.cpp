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
// bake_texture and the texture round trips through glTF and OBJ (roadmap 2.4.8).

#include <cstdio>
#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/gltf.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/operations/bake_texture.hpp"
#include "meshioplusplus/region.hpp"

namespace {

meshioplusplus::Mesh bake_test_triangle() {
    using namespace meshioplusplus;
    Mesh m;
    NDArray pts(DType::Float64, {3, 3});
    const double xyz[] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    std::copy(std::begin(xyz), std::end(xyz), pts.As<double>());
    m.AssignPoints(std::move(pts));
    NDArray tri(DType::Int64, {1, 3});
    for (std::int64_t k = 0; k < 3; ++k)
        tri.As<std::int64_t>()[k] = k;
    m.AddCellBlock("triangle", std::move(tri));
    m.AddPointData("uv", mt::data_array({0, 0, 1, 0, 0, 1}, 2));
    m.AddPointData("h", mt::data_array({0, 1, 2}));
    return m;
}

}  // namespace

TEST(BakeTexture, DrawsTheTriangleAndCarriesThePng) {
    meshioplusplus::BakeTextureOptions o;
    o.mArray = "h";
    o.mWidth = 16;
    o.mHeight = 16;
    const auto r = meshioplusplus::bake_texture(bake_test_triangle(), o);
    EXPECT_EQ(r.mImage.Shape(), (std::vector<std::size_t>{16, 16, 4}));
    EXPECT_GT(r.mCoveredPixels, 100);
    EXPECT_LT(r.mCoveredPixels, 256);
    EXPECT_DOUBLE_EQ(r.mVMin, 0.0);
    EXPECT_DOUBLE_EQ(r.mVMax, 2.0);
    ASSERT_TRUE(r.mMesh.HasFieldData("texture:png"));
    o.mArray = "missing";
    EXPECT_THROW(meshioplusplus::bake_texture(bake_test_triangle(), o), std::invalid_argument);
}

TEST(BakeTexture, TheTextureRoundTripsThroughGltfAndObj) {
    meshioplusplus::BakeTextureOptions o;
    o.mArray = "h";
    o.mWidth = 8;
    o.mHeight = 8;
    auto mesh = meshioplusplus::bake_texture(bake_test_triangle(), o).mMesh;
    const auto& png = mesh.FieldData("texture:png");
    const std::string glb = mt::temp_path("_baked.glb");
    meshioplusplus::write_gltf(glb, mesh);
    const auto back = meshioplusplus::read_gltf(glb);
    ASSERT_TRUE(back.HasFieldData("texture:png"));
    EXPECT_EQ(back.FieldData("texture:png").Size(), png.Size());
    const std::filesystem::path dir =
        std::filesystem::path(mt::temp_path("_baked_obj")).replace_extension("");
    std::filesystem::create_directories(dir);
    const std::string obj = (dir / "part.obj").string();
    meshioplusplus::write_obj(obj, mesh);
    EXPECT_TRUE(std::filesystem::exists(dir / "part.mtl"));
    EXPECT_TRUE(std::filesystem::exists(dir / "part.png"));
    const auto back_obj = meshioplusplus::read_obj(obj);
    ASSERT_TRUE(back_obj.HasFieldData("texture:png"));
    EXPECT_EQ(back_obj.FieldData("texture:png").Size(), png.Size());
    std::remove(glb.c_str());
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}
