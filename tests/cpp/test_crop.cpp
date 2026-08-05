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
// Tests for the crop operation (bounding box / half-space / cell_data predicate).

// System includes
#include <cmath>
#include <limits>
#include <stdexcept>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/crop.hpp"

namespace {

using meshioplusplus::crop_bbox;
using meshioplusplus::crop_halfspace;
using meshioplusplus::crop_predicate;
using meshioplusplus::CropMode;
using meshioplusplus::CropResult;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::RefineCompare;

// 2x1 grid of quads (points x in [0,2], y in [0,1]).
Mesh quad_grid() {
    return mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, "quad",
                         {{0, 1, 4, 3}, {1, 2, 5, 4}});
}

std::size_t count_used(const Mesh& m) {
    std::size_t used = 0;
    std::vector<char> seen(m.NumPoints(), 0);
    for (const auto cb : m.CellRange()) {
        const auto& conn = cb.Conn();
        for (std::size_t i = 0; i < conn.Size(); ++i)
            seen[static_cast<std::size_t>(conn.As<std::int64_t>()[i])] = 1;
    }
    for (char c : seen)
        used += c ? 1 : 0;
    return used;
}

TEST(Crop, BBoxModeAll) {
    Mesh m = quad_grid();
    double lo[3] = {-0.1, -0.1, -0.1}, hi[3] = {1.1, 1.1, 1.1};
    CropResult r = crop_bbox(m, lo, hi, CropMode::All);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 1u);           // only the fully-inside quad
    EXPECT_EQ(count_used(r.mMesh), r.mMesh.NumPoints());  // no unused points
}

TEST(Crop, BBoxModeAnyKeepsStraddling) {
    Mesh m = quad_grid();
    double lo[3] = {-0.1, -0.1, -0.1}, hi[3] = {1.1, 1.1, 1.1};
    CropResult r = crop_bbox(m, lo, hi, CropMode::Any);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 2u);
}

TEST(Crop, HalfSpace) {
    Mesh m = quad_grid();
    double point[3] = {1, 0, 0}, normal[3] = {1, 0, 0};  // keep x >= 1
    CropResult r = crop_halfspace(m, point, normal, CropMode::All);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 1u);
}

// A quad grid whose cells carry a scalar `t` = 0, 1, ... and, when asked, a NaN.
Mesh tagged_quad_grid(bool WithNan) {
    Mesh m = mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}, {3, 1, 0}},
        "quad", {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}});
    NDArray t = NDArray::Uninit(meshioplusplus::DType::Float64, {3u});
    t.As<double>()[0] = 0.0;
    t.As<double>()[1] = WithNan ? std::numeric_limits<double>::quiet_NaN() : 1.0;
    t.As<double>()[2] = 2.0;
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(t));
    m.AddCellData("t", std::move(blocks));
    return m;
}

TEST(Crop, PredicateKeepsTheMatchingCells) {
    const CropResult r = crop_predicate(tagged_quad_grid(false), "t", RefineCompare::Less, 1.5);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 2u);
    // Points the dropped cell alone referenced are pruned, exactly as the
    // spatial crops do -- all three funnel through the same subset builder.
    EXPECT_EQ(count_used(r.mMesh), r.mMesh.NumPoints());
    EXPECT_EQ(r.mMesh.NumPoints(), 6u);
    // The predicate array itself rides through, subsetted.
    ASSERT_TRUE(r.mMesh.HasCellData("t"));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(r.mMesh.CellData("t", 0), 0), 0.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(r.mMesh.CellData("t", 0), 1), 1.0);
}

TEST(Crop, PredicateHonoursEveryComparison) {
    const auto n = [](RefineCompare op, double v) {
        return crop_predicate(tagged_quad_grid(false), "t", op, v).mMesh.Cells(0).NumCells();
    };
    EXPECT_EQ(n(RefineCompare::Less, 1.0), 1u);
    EXPECT_EQ(n(RefineCompare::LessEqual, 1.0), 2u);
    EXPECT_EQ(n(RefineCompare::Greater, 1.0), 1u);
    EXPECT_EQ(n(RefineCompare::GreaterEqual, 1.0), 2u);
    EXPECT_EQ(n(RefineCompare::Equal, 1.0), 1u);
    EXPECT_EQ(n(RefineCompare::NotEqual, 1.0), 2u);
}

// The rule `refine` already states and this operation inherits by sharing the
// evaluator: a non-finite value never matches, whatever the comparison. It has
// to hold for `!=` too, which is the case a naive implementation gets wrong --
// `NaN != 1.0` is true in IEEE.
TEST(Crop, ANonFiniteCellValueNeverMatches) {
    for (RefineCompare op :
         {RefineCompare::Less, RefineCompare::LessEqual, RefineCompare::Greater,
          RefineCompare::GreaterEqual, RefineCompare::Equal, RefineCompare::NotEqual}) {
        const CropResult r = crop_predicate(tagged_quad_grid(true), "t", op, 1.0);
        for (std::size_t c = 0; c < r.mMesh.Cells(0).NumCells(); ++c)
            EXPECT_TRUE(
                std::isfinite(meshioplusplus::detail::read_double(r.mMesh.CellData("t", 0), c)))
                << "the NaN cell survived comparison " << static_cast<int>(op);
    }
}

TEST(Crop, PredicateRejectsWhatIsNotAScalarCellArray) {
    Mesh m = tagged_quad_grid(false);
    // An unknown name, with the available ones listed.
    EXPECT_THROW(crop_predicate(m, "nope", RefineCompare::Less, 1.0), std::invalid_argument);
    // point_data is refused BY NAME rather than averaged onto the cells.
    NDArray pd = NDArray(meshioplusplus::DType::Float64, {m.NumPoints()});
    m.AddPointData("p", std::move(pd));
    try {
        crop_predicate(m, "p", RefineCompare::Less, 1.0);
        FAIL() << "expected a throw";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("point_data"), std::string::npos);
    }
    // A vector array has no single value per cell to compare.
    NDArray v = NDArray(meshioplusplus::DType::Float64, {3u, 3u});
    std::vector<NDArray> vb;
    vb.push_back(std::move(v));
    m.AddCellData("v", std::move(vb));
    EXPECT_THROW(crop_predicate(m, "v", RefineCompare::Less, 1.0), std::invalid_argument);
}

TEST(Crop, PredicateRecordsIdsLikeTheOthers) {
    const CropResult r =
        crop_predicate(tagged_quad_grid(false), "t", RefineCompare::Less, 1.5, true);
    EXPECT_TRUE(r.mMesh.HasPointData("crop:original_point_id"));
    EXPECT_TRUE(r.mMesh.HasCellData("crop:original_cell_id"));
    EXPECT_EQ(meshioplusplus::detail::read_int(r.mMesh.CellData("crop:original_cell_id", 0), 1), 1);
}

TEST(Crop, RecordIds) {
    Mesh m = quad_grid();
    double lo[3] = {-0.1, -0.1, -0.1}, hi[3] = {1.1, 1.1, 1.1};
    CropResult r = crop_bbox(m, lo, hi, CropMode::All, /*record_ids=*/true);
    EXPECT_TRUE(r.mMesh.HasPointData("crop:original_point_id"));
    EXPECT_TRUE(r.mMesh.HasCellData("crop:original_cell_id"));
}

}  // namespace
