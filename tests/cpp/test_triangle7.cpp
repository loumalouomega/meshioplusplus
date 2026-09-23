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
 * @file test_triangle7.cpp
 * @brief `triangle7` (v16.0.0): a `triangle6` with a centre node, as Code_Aster
 *        (`TRIA7`), MED (`TR7`) and VTK (34) define it. The cell-type tables,
 *        the edge table, the corner count and the linear base.
 */

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/types.hpp"

using namespace meshioplusplus;

namespace {

Mesh triangle7_mesh() {
    return mt::make_mesh({{0, 0, 0},
                          {1, 0, 0},
                          {0, 1, 0},
                          {0.5, 0, 0},
                          {0.5, 0.5, 0},
                          {0, 0.5, 0},
                          {1.0 / 3.0, 1.0 / 3.0, 0}},
                         "triangle7", {{0, 1, 2, 3, 4, 5, 6}});
}

}  // namespace

TEST(Triangle7, IsAKnownSevenNodeSurfaceCell) {
    const CellType t = cell_type_from_name("triangle7");
    ASSERT_EQ(t, CellType::Triangle7);
    EXPECT_EQ(cell_type_name(t), "triangle7");
    EXPECT_EQ(cell_type_num_nodes(t), 7);
    EXPECT_EQ(cell_type_dimension(t), 2);
    EXPECT_EQ(detail::cell_corner_count(t), 3);
}

TEST(Triangle7, HasTheEdgesOfATriangle6) {
    const auto& t7 = detail::cell_edges(CellType::Triangle7);
    const auto& t6 = detail::cell_edges(CellType::Triangle6);
    ASSERT_EQ(t7.size(), 3u);
    for (std::size_t i = 0; i < t7.size(); ++i) {
        EXPECT_EQ(t7[i].mEdgeType, t6[i].mEdgeType);
        EXPECT_EQ(t7[i].mNodes, t6[i].mNodes);
    }
}

TEST(Triangle7, LinearizesToATriangle) {
    const ConvertCellsResult r = convert_cells(triangle7_mesh());
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "triangle");
    EXPECT_EQ(r.mMesh.NumPoints(), 3u);
}

TEST(Triangle7, SurfaceIsItsThreeQuadraticEdges) {
    const Mesh s = extract_surface(triangle7_mesh());
    ASSERT_EQ(s.NumCellBlocks(), 1u);
    EXPECT_EQ(s.Cells(0).Type(), "line3");
    EXPECT_EQ(s.Cells(0).NumCells(), 3u);
}
