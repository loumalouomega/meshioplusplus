// Tests for the core value types: NDArray, Mesh/CellBlock, exceptions.

#include <gtest/gtest.h>

#include <cstdint>

#include "meshio/exceptions.hpp"
#include "meshio/mesh.hpp"
#include "meshio/ndarray.hpp"

using meshio::CellBlock;
using meshio::DType;
using meshio::Mesh;
using meshio::NDArray;

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
    EXPECT_EQ(meshio::dtype_size(DType::Float32), 4u);
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

TEST(Exceptions, ReadWriteErrorMessages) {
    try {
        throw meshio::ReadError("boom-read");
    } catch (const meshio::ReadError& e) {
        EXPECT_STREQ(e.what(), "boom-read");
    }
    try {
        throw meshio::WriteError("boom-write");
    } catch (const std::runtime_error& e) {  // both derive from runtime_error
        EXPECT_STREQ(e.what(), "boom-write");
    }
    EXPECT_THROW(throw meshio::ReadError("x"), meshio::ReadError);
    EXPECT_THROW(throw meshio::WriteError("x"), meshio::WriteError);
}
