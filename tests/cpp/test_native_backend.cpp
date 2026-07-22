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
 * @file test_native_backend.cpp
 * @brief NATIVE-backend-specific tests: dtype canonicalization, the
 * move-not-copy ingest fast path, CSR ragged storage invariants, and the
 * fast-consumer surface (`GlobalConnectivity`, `ConnSpan`, ...).
 *
 * The whole file compiles away under other backends (the CMake glob picks it
 * up unconditionally; the `#ifdef` keeps non-NATIVE builds clean).
 */

#ifdef MESHIOPLUSPLUS_MESH_BACKEND_NATIVE

// System includes
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/mesh.hpp"

using meshioplusplus::CellType;
using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

TEST(NativeBackend, BackendName) {
    EXPECT_STREQ(meshioplusplus::mesh_backend_name(), "native");
}

TEST(NativeBackend, CanonicalizesDtypesWithinKind) {
    Mesh m;
    // Float32 points -> Float64.
    NDArray pts32 = NDArray::Uninit(DType::Float32, {2, 3});
    for (int i = 0; i < 6; ++i)
        pts32.As<float>()[i] = 0.5f * static_cast<float>(i);
    m.AssignPoints(std::move(pts32));
    EXPECT_EQ(m.Points().Dtype(), DType::Float64);
    EXPECT_DOUBLE_EQ(m.Points().As<double>()[5], 2.5);

    // Int32 connectivity -> Int64.
    NDArray conn32 = NDArray::Uninit(DType::Int32, {1, 3});
    for (int i = 0; i < 3; ++i)
        conn32.As<std::int32_t>()[i] = i;
    m.AddCellBlock("triangle", std::move(conn32));
    EXPECT_EQ(m.Cells(0).Conn().Dtype(), DType::Int64);
    EXPECT_EQ(m.Cells(0).Conn().As<std::int64_t>()[2], 2);

    // Integer data stays integer kind (Int64), float data becomes Float64.
    NDArray tag = NDArray::Uninit(DType::UInt16, {1});
    tag.As<std::uint16_t>()[0] = 7;
    m.AddCellData("tag", {std::move(tag)});
    EXPECT_EQ(m.CellData("tag", 0).Dtype(), DType::Int64);
    NDArray temp = NDArray::Uninit(DType::Float32, {2});
    temp.As<float>()[0] = 1.5f;
    temp.As<float>()[1] = 2.5f;
    m.AddPointData("temp", std::move(temp));
    EXPECT_EQ(m.PointData("temp").Dtype(), DType::Float64);
    EXPECT_DOUBLE_EQ(m.PointData("temp").As<double>()[1], 2.5);
}

TEST(NativeBackend, CanonicalOwningArraysAreMovedNotCopied) {
    NDArray pts = mt::points_from({{0, 0, 0}, {1, 0, 0}});
    const std::byte* p_before = pts.Data();
    Mesh m;
    m.AssignPoints(std::move(pts));
    EXPECT_EQ(m.Points().Data(), p_before);  // same buffer: moved, not copied

    NDArray conn = mt::conn_from({{0, 1}});
    const std::byte* c_before = conn.Data();
    m.AddCellBlock("line", std::move(conn));
    EXPECT_EQ(m.Cells(0).Conn().Data(), c_before);
}

TEST(NativeBackend, ViewsAreOwnedOnIngest) {
    std::vector<double> buf = {0, 0, 0, 1, 0, 0};
    NDArray view =
        NDArray::MakeView(DType::Float64, {2, 3}, reinterpret_cast<std::byte*>(buf.data()));
    Mesh m;
    m.AssignPoints(std::move(view));
    EXPECT_FALSE(m.Points().IsView());
    buf[3] = 99.0;  // mutating the source must not affect the mesh
    EXPECT_DOUBLE_EQ(m.Points().As<double>()[3], 1.0);
}

TEST(NativeBackend, PolygonCsrLayout) {
    Mesh m;
    m.AssignPoints(NDArray(DType::Float64, {6, 3}));
    m.AddPolygonBlock("polygon", {{0, 1, 2}, {1, 3, 4, 2, 5}});
    const auto& r_block = m.Blocks()[0];
    EXPECT_EQ(r_block.mType, CellType::Polygon);
    EXPECT_EQ(r_block.mRowOffsets, (std::vector<std::int64_t>{0, 3, 8}));
    EXPECT_EQ(r_block.mFlat, (std::vector<std::int64_t>{0, 1, 2, 1, 3, 4, 2, 5}));
    EXPECT_TRUE(r_block.mFaceOffsets.empty());
}

TEST(NativeBackend, PolyhedronCsrLayout) {
    Mesh m;
    m.AssignPoints(NDArray(DType::Float64, {8, 3}));
    m.AddPolyhedronBlock("polyhedron4", {
                                            {{0, 1, 2}, {0, 1, 3}},
                                            {{4, 5, 6}, {4, 5, 7}, {5, 6, 7}},
                                        });
    const auto& r_block = m.Blocks()[0];
    EXPECT_EQ(r_block.mFaceOffsets, (std::vector<std::int64_t>{0, 2, 5}));
    EXPECT_EQ(r_block.mRowOffsets, (std::vector<std::int64_t>{0, 3, 6, 9, 12, 15}));
    ASSERT_EQ(m.Cells(0).NumCells(), 2u);
    EXPECT_EQ(m.Cells(0).NumFaces(1), 3u);
    const auto face = m.Cells(0).Face(1, 2);
    ASSERT_EQ(face.second, 3u);
    EXPECT_EQ(face.first[0], 5);
}

TEST(NativeBackend, FastConsumerSurface) {
    Mesh m = mt::tri_quad_mesh();
    EXPECT_EQ(m.BlockType(0), CellType::Triangle);
    EXPECT_EQ(m.BlockType(1), CellType::Quad);
    EXPECT_EQ(m.ConnSpan(1).size(), 4u);
    EXPECT_EQ(m.ConnSpan(1)[2], 4);
    EXPECT_DOUBLE_EQ(m.PointsData()[3], 1.0);  // point 1, x
}

TEST(NativeBackend, GlobalConnectivityCsr) {
    Mesh m = mt::tri_quad_mesh();  // 2 tri + 1 quad + 1 tri = 4 cells
    const auto& r_csr = m.GlobalConnectivity();
    ASSERT_EQ(r_csr.mTypes.size(), 4u);
    EXPECT_EQ(r_csr.mOffsets, (std::vector<std::int64_t>{0, 3, 6, 10, 13}));
    EXPECT_EQ(r_csr.mTypes[2], CellType::Quad);
    // Cell 2 (the quad) spans mConn[6..10).
    EXPECT_EQ(r_csr.mConn[6], 1);
    EXPECT_EQ(r_csr.mConn[9], 5);

    // The cache is invalidated when a block is added.
    m.AddCellBlock("line", mt::conn_from({{0, 1}}));
    EXPECT_EQ(m.GlobalConnectivity().mTypes.size(), 5u);
}

TEST(NativeBackend, CustomTypeNameIsPreserved) {
    Mesh m;
    m.AssignPoints(NDArray(DType::Float64, {8, 3}));
    m.AddPolyhedronBlock("polyhedron12", {{{0, 1, 2}}});
    EXPECT_EQ(m.Cells(0).Type(), "polyhedron12");
    EXPECT_EQ(m.BlockType(0), CellType::Custom);
}

#endif  // MESHIOPLUSPLUS_MESH_BACKEND_NATIVE
