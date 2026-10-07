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
// Tests for the affine transform operation.

// System includes
#include <cmath>
#include <cstdint>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/transform.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::transform;

Mesh unit_cube() {
    return mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
}

double px(const Mesh& m, std::size_t node, std::size_t comp) {
    return meshioplusplus::detail::read_double(m.Points(), node * m.PointDim() + comp);
}

TEST(Transform, TranslateMovesPointsOnly) {
    Mesh m = unit_cube();
    Mesh out = transform(m, meshioplusplus::transform_translation(10, 20, 30));
    EXPECT_NEAR(px(out, 1, 0), 11.0, 1e-12);
    EXPECT_NEAR(px(out, 1, 1), 20.0, 1e-12);
    // connectivity unchanged
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).NumCells(), 1u);
}

TEST(Transform, Rotate90AboutZ) {
    Mesh m = unit_cube();
    Mesh out = transform(m, meshioplusplus::transform_rotation(0, 0, 1, M_PI / 2.0));
    // (1,0,0) -> (0,1,0)
    EXPECT_NEAR(px(out, 1, 0), 0.0, 1e-12);
    EXPECT_NEAR(px(out, 1, 1), 1.0, 1e-12);
}

TEST(Transform, UnitScale) {
    Mesh m = unit_cube();
    Mesh out = transform(m, meshioplusplus::transform_units(0.001));
    EXPECT_NEAR(px(out, 6, 0), 0.001, 1e-15);
    EXPECT_NEAR(px(out, 6, 2), 0.001, 1e-15);
}

TEST(Transform, MatrixEqualsTranslate) {
    Mesh m = unit_cube();
    double mat[16] = {1, 0, 0, 5, 0, 1, 0, 6, 0, 0, 1, 7, 0, 0, 0, 1};
    Mesh out = transform(m, meshioplusplus::transform_from_matrix(mat));
    EXPECT_NEAR(px(out, 0, 0), 5.0, 1e-12);
    EXPECT_NEAR(px(out, 0, 1), 6.0, 1e-12);
    EXPECT_NEAR(px(out, 0, 2), 7.0, 1e-12);
}

double cd(const Mesh& m, const char* name, std::size_t comp) {
    return meshioplusplus::detail::read_double(m.CellData(name, 0), comp);
}

// A vector living on a CELL rotates exactly as one living on a point does.
// Before v10.33.0 `rotate_vector_data` reached point data only, and cell data
// rode through in the old frame -- silently, since nothing about the array
// says which frame it is in.
TEST(Transform, RotateVectorDataReachesCellDataToo) {
    Mesh m = unit_cube();
    m.AddCellData("vel", {mt::data_array({1.0, 0.0, 0.0}, 3)});
    Mesh out = transform(m, meshioplusplus::transform_rotation(0, 0, 1, M_PI / 2.0), true);
    // (1,0,0) -> (0,1,0)
    EXPECT_NEAR(cd(out, "vel", 0), 0.0, 1e-12);
    EXPECT_NEAR(cd(out, "vel", 1), 1.0, 1e-12);
    EXPECT_NEAR(cd(out, "vel", 2), 0.0, 1e-12);
}

TEST(Transform, RotateVectorDataRotatesACellTensor) {
    Mesh m = unit_cube();
    // diag(1, 2, 3) under a quarter turn about z becomes diag(2, 1, 3).
    m.AddCellData("stress", {mt::data_array({1, 0, 0, 0, 2, 0, 0, 0, 3}, 9)});
    Mesh out = transform(m, meshioplusplus::transform_rotation(0, 0, 1, M_PI / 2.0), true);
    EXPECT_NEAR(cd(out, "stress", 0), 2.0, 1e-12);
    EXPECT_NEAR(cd(out, "stress", 4), 1.0, 1e-12);
    EXPECT_NEAR(cd(out, "stress", 8), 3.0, 1e-12);
}

TEST(Transform, CellDataIsUntouchedWithoutTheFlag) {
    Mesh m = unit_cube();
    m.AddCellData("vel", {mt::data_array({1.0, 0.0, 0.0}, 3)});
    Mesh out = transform(m, meshioplusplus::transform_rotation(0, 0, 1, M_PI / 2.0));
    EXPECT_NEAR(cd(out, "vel", 0), 1.0, 1e-12);
    EXPECT_NEAR(cd(out, "vel", 1), 0.0, 1e-12);
}

TEST(Transform, IntegerCellDataIsNeverRotated) {
    Mesh m = unit_cube();
    m.AddCellData("tag", {mt::int_data_array({7, 8, 9})});
    Mesh out = transform(m, meshioplusplus::transform_rotation(0, 0, 1, M_PI / 2.0), true);
    EXPECT_EQ(cd(out, "tag", 0), 7.0);
    EXPECT_EQ(cd(out, "tag", 1), 8.0);
}

TEST(Transform, RaggedBlocksAreCopiedUnchanged) {
    // Polygon and polyhedron blocks pass through the per-block copy unchanged
    // while the points move; the Python binding relies on this since it stopped
    // refusing ragged meshes.
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {0.5, 0.5, 1}},
                           "vertex", {{0}});
    m.AddPolygonBlock("polygon", {{0, 1, 4}, {0, 1, 2, 3}});
    m.AddPolyhedronBlock("polyhedron5",
                         {{{0, 3, 2, 1}, {0, 1, 5}, {1, 2, 5}, {2, 3, 5}, {3, 0, 5}}});
    Mesh out = transform(m, meshioplusplus::transform_translation(0, 0, 1));
    EXPECT_NEAR(px(out, 5, 2), 2.0, 1e-12);
    ASSERT_EQ(out.NumCellBlocks(), 3u);

    const Mesh::CellView poly = out.Cells(1);
    ASSERT_TRUE(poly.IsRagged());
    const std::vector<std::vector<std::int64_t>> want_rows = {{0, 1, 4}, {0, 1, 2, 3}};
    ASSERT_EQ(poly.NumCells(), want_rows.size());
    for (std::size_t c = 0; c < want_rows.size(); ++c) {
        ASSERT_EQ(poly.RowSize(c), want_rows[c].size()) << c;
        for (std::size_t k = 0; k < want_rows[c].size(); ++k)
            EXPECT_EQ(poly.Row(c)[k], want_rows[c][k]) << c << "," << k;
    }

    const Mesh::CellView hedra = out.Cells(2);
    ASSERT_TRUE(hedra.IsPolyhedron());
    const std::vector<std::vector<std::int64_t>> want_faces = {
        {0, 3, 2, 1}, {0, 1, 5}, {1, 2, 5}, {2, 3, 5}, {3, 0, 5}};
    ASSERT_EQ(hedra.NumCells(), 1u);
    ASSERT_EQ(hedra.NumFaces(0), want_faces.size());
    for (std::size_t f = 0; f < want_faces.size(); ++f) {
        const auto face = hedra.Face(0, f);
        ASSERT_EQ(face.second, want_faces[f].size()) << f;
        for (std::size_t k = 0; k < face.second; ++k)
            EXPECT_EQ(face.first[k], want_faces[f][k]) << f << "," << k;
    }
}

}  // namespace
