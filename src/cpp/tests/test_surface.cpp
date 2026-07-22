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
// Tests for the surface/boundary extraction operation and its edge tables.

// System includes
#include <array>
#include <cstdint>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/surface.hpp"

namespace {

using meshioplusplus::extract_surface;
using meshioplusplus::Mesh;

std::size_t count_type(const Mesh& rMesh, const std::string& rType) {
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() == rType)
            return cb.NumCells();
    return 0;
}

}  // namespace

TEST(Surface, SingleHexHasSixQuads) {
    Mesh s = extract_surface(mt::hex_mesh());
    EXPECT_EQ(count_type(s, "quad"), 6u);
    EXPECT_EQ(s.NumPoints(), 8u);  // all points preserved
}

TEST(Surface, TwoHexStackDropsInternalFace) {
    // A 2x stack: 12 exterior quads (6+6-2 shared) => 10.
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1},
                                    {0, 0, 2},
                                    {1, 0, 2},
                                    {1, 1, 2},
                                    {0, 1, 2}}));
    m.AddCellBlock("hexahedron",
                   mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}, {4, 5, 6, 7, 8, 9, 10, 11}}));
    Mesh s = extract_surface(m, /*recordParentIds=*/true);
    EXPECT_EQ(count_type(s, "quad"), 10u);
    ASSERT_TRUE(s.HasCellData("surface:parent_cell"));
    // Every parent id is 0 or 1 (two input cells).
    const meshioplusplus::NDArray& par = s.CellData("surface:parent_cell", 0);
    for (std::size_t i = 0; i < meshioplusplus::detail::rows(par); ++i) {
        const std::int64_t p = meshioplusplus::detail::read_int(par, i);
        EXPECT_TRUE(p == 0 || p == 1);
    }
}

TEST(Surface, TetraHasFourTriangles) {
    Mesh s = extract_surface(
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra", {{0, 1, 2, 3}}));
    EXPECT_EQ(count_type(s, "triangle"), 4u);
    EXPECT_EQ(s.NumPoints(), 4u);
}

TEST(Surface, TwoDPatchGivesBoundaryEdges) {
    // Two triangles forming a unit square -> 4 boundary line segments.
    Mesh m = mt::make_mesh({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, "triangle", {{0, 1, 2}, {0, 2, 3}});
    Mesh s = extract_surface(m);
    EXPECT_EQ(count_type(s, "line"), 4u);
}

TEST(Surface, PointsPreserved) {
    Mesh src = mt::hex_mesh();
    Mesh s = extract_surface(src);
    // The cube's 8 corner coordinates survive unchanged.
    mt::expect_points_close(src, s, 1e-12);
}

TEST(Surface, NoSupportedCellsThrows) {
    EXPECT_THROW(extract_surface(mt::line_mesh()), std::invalid_argument);
}

TEST(CellEdges, MidNodesAreMidpoints) {
    // On a reference triangle6, each line3 mid node is the exact midpoint of
    // its two corner nodes (mirrors the cell_faces winding invariant).
    const std::array<std::array<double, 2>, 6> tri6 = {{
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},  // corners
        {0.5, 0.0},
        {0.5, 0.5},
        {0.0, 0.5},  // mids: (0,1), (1,2), (2,0)
    }};
    for (const meshioplusplus::detail::CellEdgeDef& e :
         meshioplusplus::detail::cell_edges(meshioplusplus::CellType::Triangle6)) {
        const auto& a = tri6[e.mNodes[0]];
        const auto& b = tri6[e.mNodes[1]];
        const auto& mid = tri6[e.mNodes[2]];
        EXPECT_NEAR(mid[0], 0.5 * (a[0] + b[0]), 1e-12);
        EXPECT_NEAR(mid[1], 0.5 * (a[1] + b[1]), 1e-12);
    }
}
