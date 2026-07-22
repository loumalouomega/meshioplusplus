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
// Tests for the crop operation (bounding box / half-space).

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/crop.hpp"

namespace {

using meshioplusplus::crop_bbox;
using meshioplusplus::crop_halfspace;
using meshioplusplus::CropMode;
using meshioplusplus::CropResult;
using meshioplusplus::Mesh;

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

TEST(Crop, RecordIds) {
    Mesh m = quad_grid();
    double lo[3] = {-0.1, -0.1, -0.1}, hi[3] = {1.1, 1.1, 1.1};
    CropResult r = crop_bbox(m, lo, hi, CropMode::All, /*record_ids=*/true);
    EXPECT_TRUE(r.mMesh.HasPointData("crop:original_point_id"));
    EXPECT_TRUE(r.mMesh.HasCellData("crop:original_cell_id"));
}

}  // namespace
