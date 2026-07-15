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
#include <cstdint>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

using meshioplusplus::CellBlock;
using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

TEST(NDArray, OwningConstructionAndAccess) {
    NDArray a(DType::Float64, {2, 3});
    EXPECT_EQ(a.ndim(), 2u);
    EXPECT_EQ(a.size(), 6u);
    EXPECT_EQ(a.nbytes(), 6u * 8u);
    EXPECT_FALSE(a.is_view());
    // zero-initialised
    for (std::size_t i = 0; i < a.size(); ++i) EXPECT_EQ(a.as<double>()[i], 0.0);
    a.as<double>()[5] = 3.5;
    EXPECT_EQ(a.as<double>()[5], 3.5);
}

TEST(NDArray, IntegerDtypeSizes) {
    EXPECT_EQ(NDArray(DType::Int32, {4}).nbytes(), 16u);
    EXPECT_EQ(NDArray(DType::Int64, {4}).nbytes(), 32u);
    EXPECT_EQ(NDArray(DType::UInt8, {4}).nbytes(), 4u);
    EXPECT_EQ(meshioplusplus::dtype_size(DType::Float32), 4u);
}

TEST(NDArray, Reshape) {
    NDArray a(DType::Int32, {2, 3});
    a.reshape({3, 2});
    EXPECT_EQ(a.shape()[0], 3u);
    EXPECT_EQ(a.shape()[1], 2u);
    // inconsistent reshape is ignored
    a.reshape({5, 5});
    EXPECT_EQ(a.size(), 6u);
}

TEST(NDArray, ViewBecomesOwnedCopy) {
    std::vector<std::int64_t> buf = {1, 2, 3, 4};
    NDArray v = NDArray::make_view(DType::Int64, {4},
                                   reinterpret_cast<std::byte*>(buf.data()));
    EXPECT_TRUE(v.is_view());
    EXPECT_EQ(v.as<std::int64_t>()[2], 3);
    v.make_owned();
    EXPECT_FALSE(v.is_view());
    buf[2] = 99;  // mutating the original no longer affects the owned copy
    EXPECT_EQ(v.as<std::int64_t>()[2], 3);
}

TEST(Mesh, CountsAndCellBlock) {
    Mesh m;
    m.points = NDArray(DType::Float64, {5, 3});
    EXPECT_EQ(m.num_points(), 5u);
    m.cells.emplace_back("tetra", NDArray(DType::Int64, {2, 4}));
    EXPECT_EQ(m.cells.size(), 1u);
    EXPECT_EQ(m.cells[0].type, "tetra");
    EXPECT_EQ(m.cells[0].num_cells(), 2u);
    Mesh empty;
    EXPECT_EQ(empty.num_points(), 0u);
}

TEST(Mesh, RaggedCellBlocks) {
    // 1-level ragged (jagged polygon): cells with varying node counts.
    CellBlock poly;
    poly.type = "polygon";
    poly.polygon_rows = {{0, 1, 2}, {1, 2, 3, 4, 5}};
    EXPECT_TRUE(poly.is_ragged());
    EXPECT_EQ(poly.num_cells(), 2u);
    EXPECT_EQ(poly.polygon_rows[1].size(), 5u);

    // 2-level ragged (polyhedron): each cell is a list of faces.
    CellBlock poh;
    poh.type = "polyhedron4";
    poh.polyhedron_rows = {
        {{1, 2, 5}, {1, 2, 7}, {1, 5, 7}, {2, 5, 7}},
        {{2, 5, 6}, {2, 6, 7}, {2, 5, 7}, {5, 6, 7}},
    };
    EXPECT_TRUE(poh.is_ragged());
    EXPECT_EQ(poh.num_cells(), 2u);
    EXPECT_EQ(poh.polyhedron_rows[0].size(), 4u);  // 4 faces

    // A rectangular block is not ragged.
    CellBlock tri("triangle", NDArray(DType::Int64, {3, 3}));
    EXPECT_FALSE(tri.is_ragged());
    EXPECT_EQ(tri.num_cells(), 3u);
}

TEST(Exceptions, ReadWriteErrorMessages) {
    try {
        throw meshioplusplus::ReadError("boom-read");
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_STREQ(e.what(), "boom-read");
    }
    try {
        throw meshioplusplus::WriteError("boom-write");
    } catch (const std::runtime_error& e) {  // both derive from runtime_error
        EXPECT_STREQ(e.what(), "boom-write");
    }
    EXPECT_THROW(throw meshioplusplus::ReadError("x"), meshioplusplus::ReadError);
    EXPECT_THROW(throw meshioplusplus::WriteError("x"), meshioplusplus::WriteError);
}
