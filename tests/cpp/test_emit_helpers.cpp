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

// The private emission helpers the operations share: data_plumbing.hpp's
// append_rows (convert_cells, refine, tet_poly) and decimate_shared.hpp's
// decim_emit_blocks (decimate, decimate_volume).

#include <cstddef>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

#include "../../src/cpp/src/detail/data_plumbing.hpp"
#include "../../src/cpp/src/detail/decimate_shared.hpp"

using meshioplusplus::CellType;
using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

namespace {

NDArray emit_int64_block(std::size_t Rows, std::size_t Cols, std::int64_t First) {
    NDArray a = NDArray::Uninit(DType::Int64, {Rows, Cols});
    for (std::size_t i = 0; i < Rows * Cols; ++i)
        a.As<std::int64_t>()[i] = First + static_cast<std::int64_t>(i);
    return a;
}

}  // namespace

TEST(AppendRows, KeepsTheOriginalBytesAndStoresTheNewRowsInTheDtype) {
    NDArray a = NDArray::Uninit(DType::Float32, {2, 3});
    for (std::size_t i = 0; i < 6; ++i)
        a.As<float>()[i] = 0.5f * static_cast<float>(i);
    const NDArray b = meshioplusplus::detail::append_rows(a, 2, 2, 3, [](double* pNew) {
        for (std::size_t i = 0; i < 6; ++i)
            pNew[i] = 10.0 + static_cast<double>(i) + 0.25;
    });
    ASSERT_EQ(b.Dtype(), DType::Float32);
    ASSERT_EQ(b.Shape(), (std::vector<std::size_t>{4, 3}));
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_EQ(b.As<float>()[i], 0.5f * static_cast<float>(i));
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_EQ(b.As<float>()[6 + i], static_cast<float>(10.0 + static_cast<double>(i) + 0.25));
}

TEST(AppendRows, OneDimensionalArrayAndNoNewRows) {
    NDArray a = NDArray::Uninit(DType::Int32, {3});
    for (std::size_t i = 0; i < 3; ++i)
        a.As<std::int32_t>()[i] = static_cast<std::int32_t>(7 * i);
    const NDArray same = meshioplusplus::detail::append_rows(a, 3, 0, 1, [](double*) {});
    ASSERT_EQ(same.Shape(), (std::vector<std::size_t>{3}));
    for (std::size_t i = 0; i < 3; ++i)
        EXPECT_EQ(same.As<std::int32_t>()[i], static_cast<std::int32_t>(7 * i));
    const NDArray grown =
        meshioplusplus::detail::append_rows(a, 3, 1, 1, [](double* pNew) { pNew[0] = -4.0; });
    ASSERT_EQ(grown.Shape(), (std::vector<std::size_t>{4}));
    EXPECT_EQ(grown.As<std::int32_t>()[3], -4);
}

TEST(DecimEmitBlocks, IdentityParentsKeepTheAliveCellsInOrder) {
    // Two tet blocks of 2 and 1 cells over 6 points; cell 1 of block 0 dies.
    Mesh in;
    in.AssignPoints(NDArray::Uninit(DType::Float64, {6, 3}));
    in.AddCellBlock("tetra", emit_int64_block(2, 4, 0));
    in.AddCellBlock("tetra", emit_int64_block(1, 4, 2));
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
    const std::vector<std::uint8_t> alive = {1, 0, 1};
    const std::vector<std::int64_t> bases = {0, 2};
    const std::vector<std::int64_t> remap = {5, 4, 3, 2, 1, 0};
    Mesh out;
    std::vector<NDArray> maps;
    meshioplusplus::detail::decim_emit_blocks<4>(in, in, nullptr, alive, conn, bases, remap,
                                                 CellType::Tetra, out, maps);
    ASSERT_EQ(out.NumCellBlocks(), 2u);
    ASSERT_EQ(maps.size(), 2u);
    const auto blocks = out.CellRange();
    auto it = blocks.begin();
    const NDArray& b0 = (*it).Conn();
    ASSERT_EQ(b0.Shape(), (std::vector<std::size_t>{1, 4}));
    EXPECT_EQ(std::vector<std::int64_t>(b0.As<std::int64_t>(), b0.As<std::int64_t>() + 4),
              (std::vector<std::int64_t>{5, 4, 3, 2}));
    ++it;
    const NDArray& b1 = (*it).Conn();
    EXPECT_EQ(std::vector<std::int64_t>(b1.As<std::int64_t>(), b1.As<std::int64_t>() + 4),
              (std::vector<std::int64_t>{3, 2, 1, 0}));
    EXPECT_EQ(std::vector<std::int64_t>(maps[0].As<std::int64_t>(),
                                        maps[0].As<std::int64_t>() + maps[0].Size()),
              (std::vector<std::int64_t>{0, -1}));
    EXPECT_EQ(maps[1].As<std::int64_t>()[0], 0);
}

TEST(DecimEmitBlocks, ParentArrayMapsEachInputCellToItsFirstSurvivingChild) {
    // One input quad split into two triangles (parent 0 twice) plus one
    // triangle (parent 1); the first child dies, so the quad maps to row 0,
    // the second child's row, and the triangle to row 1.
    Mesh in;
    in.AssignPoints(NDArray::Uninit(DType::Float64, {5, 3}));
    in.AddCellBlock("quad", emit_int64_block(2, 4, 0));
    Mesh simp;
    simp.AssignPoints(NDArray::Uninit(DType::Float64, {5, 3}));
    simp.AddCellBlock("triangle", emit_int64_block(3, 3, 0));
    NDArray parent = NDArray::Uninit(DType::Int64, {3});
    parent.As<std::int64_t>()[0] = 0;
    parent.As<std::int64_t>()[1] = 0;
    parent.As<std::int64_t>()[2] = 1;
    std::vector<NDArray> parent_blocks;
    parent_blocks.push_back(std::move(parent));
    simp.AddCellData("convert:parent_cell", std::move(parent_blocks));
    const std::vector<std::int64_t> conn = {0, 1, 2, 0, 2, 3, 1, 3, 4};
    const std::vector<std::uint8_t> alive = {0, 1, 1};
    const std::vector<std::int64_t> bases = {0};
    const std::vector<std::int64_t> remap = {0, 1, 2, 3, 4};
    Mesh out;
    std::vector<NDArray> maps;
    meshioplusplus::detail::decim_emit_blocks<3>(in, simp, "convert:parent_cell", alive, conn,
                                                 bases, remap, CellType::Triangle, out, maps);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    ASSERT_EQ(maps.size(), 1u);
    EXPECT_EQ(std::vector<std::int64_t>(maps[0].As<std::int64_t>(),
                                        maps[0].As<std::int64_t>() + maps[0].Size()),
              (std::vector<std::int64_t>{0, 1}));
    const NDArray& tris = (*out.CellRange().begin()).Conn();
    ASSERT_EQ(tris.Shape(), (std::vector<std::size_t>{2, 3}));
    EXPECT_EQ(std::vector<std::int64_t>(tris.As<std::int64_t>(), tris.As<std::int64_t>() + 6),
              (std::vector<std::int64_t>{0, 2, 3, 1, 3, 4}));
}
