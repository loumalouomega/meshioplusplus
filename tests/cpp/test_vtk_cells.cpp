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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// detail/vtk_cells: `slice_rows` on arrays with no rows to copy.
//
// A cell-data array of shape (n, 0) owns no buffer, and `memcpy` with a null
// pointer is undefined even for zero bytes. A mutated VTKHDF file reached this
// through `vtkhdf_append_cell_data` (the fuzz campaign's finding under UBSan's
// nonnull-attribute check), so the test is meant to run in the sanitized job.

// External includes
#include <gtest/gtest.h>

// System includes
#include <cstddef>
#include <vector>

// Project includes
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace {

namespace mio = meshioplusplus;

TEST(VtkCells, SliceRowsOfAZeroWidthArrayCopiesNothing) {
    const mio::NDArray zero_width = mio::NDArray::Uninit(mio::DType::Float64, {4, 0});
    const mio::NDArray out = mio::detail::slice_rows(zero_width, 1, 3);
    ASSERT_EQ(out.Shape().size(), 2u);
    EXPECT_EQ(out.Shape()[0], 2u);
    EXPECT_EQ(out.Shape()[1], 0u);
    EXPECT_EQ(out.Size(), 0u);
}

TEST(VtkCells, SliceRowsOfAnEmptyRangeCopiesNothing) {
    const mio::NDArray a = mio::NDArray::Uninit(mio::DType::Float64, {4, 2});
    const mio::NDArray out = mio::detail::slice_rows(a, 2, 2);
    ASSERT_EQ(out.Shape().size(), 2u);
    EXPECT_EQ(out.Shape()[0], 0u);
    EXPECT_EQ(out.Shape()[1], 2u);
}

TEST(VtkCells, SliceRowsStillCopiesTheRowsOfAnOrdinaryArray) {
    mio::NDArray a = mio::NDArray::Uninit(mio::DType::Float64, {4, 2});
    for (std::size_t i = 0; i < 8; ++i)
        a.As<double>()[i] = static_cast<double>(i);
    const mio::NDArray out = mio::detail::slice_rows(a, 1, 3);
    ASSERT_EQ(out.Shape().size(), 2u);
    ASSERT_EQ(out.Shape()[0], 2u);
    ASSERT_EQ(out.Shape()[1], 2u);
    const std::vector<double> want = {2.0, 3.0, 4.0, 5.0};
    for (std::size_t i = 0; i < want.size(); ++i)
        EXPECT_EQ(out.As<double>()[i], want[i]);
}

}  // namespace
