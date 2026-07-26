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
 * @file test_mesh_api.cpp
 * @brief Backend-agnostic tests of the uniform format-facing mesh API
 * (`mesh_api.hpp`).
 *
 * Every assertion here must hold for ALL mesh backends (MESHIO, NATIVE,
 * KRATOS) — this file is the executable form of the API contract, run under
 * each backend leg of CI. Backend-specific behavior belongs in
 * `test_native_backend.cpp` / `test_kratos_backend.cpp`.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

namespace {

NDArray int32_array(const std::vector<std::int32_t>& rVals) {
    NDArray a = NDArray::Uninit(DType::Int32, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int32_t>()[i] = rVals[i];
    return a;
}

NDArray float64_array(const std::vector<double>& rVals) {
    NDArray a = NDArray::Uninit(DType::Float64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<double>()[i] = rVals[i];
    return a;
}

}  // namespace

TEST(MeshApi, BackendNameIsReported) {
    const std::string name = meshioplusplus::mesh_backend_name();
    EXPECT_TRUE(name == "meshio" || name == "native" || name == "kratos");
}

TEST(MeshApi, PointsRoundTrip) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}));
    EXPECT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(m.PointDim(), 3u);
    const NDArray& pts = m.Points();
    EXPECT_EQ(meshioplusplus::detail::read_double(pts, 3), 1.0);
    Mesh empty;
    EXPECT_EQ(empty.NumPoints(), 0u);
    EXPECT_EQ(empty.PointDim(), 0u);
}

TEST(MeshApi, CellBlocksAndRange) {
    Mesh m = mt::tri_quad_mesh();
    ASSERT_EQ(m.NumCellBlocks(), 3u);
    EXPECT_EQ(m.Cells(0).Type(), "triangle");
    EXPECT_EQ(m.Cells(1).Type(), "quad");
    EXPECT_EQ(m.Cells(2).Type(), "triangle");
    EXPECT_EQ(m.Cells(1).NumCells(), 1u);
    EXPECT_EQ(m.Cells(1).NodesPerCell(), 4u);
    EXPECT_FALSE(m.Cells(1).IsRagged());

    // Range iteration visits blocks in insertion order.
    std::vector<std::string> types;
    std::size_t total = 0;
    for (const auto cb : m.CellRange()) {
        types.push_back(cb.Type());
        total += cb.NumCells();
    }
    EXPECT_EQ(types, (std::vector<std::string>{"triangle", "quad", "triangle"}));
    EXPECT_EQ(total, 4u);

    // Connectivity is readable through the view regardless of backend.
    const NDArray& conn = m.Cells(1).Conn();
    EXPECT_EQ(meshioplusplus::detail::read_int(conn, 2), 4);
}

TEST(MeshApi, RaggedBlocks) {
    Mesh m;
    m.AssignPoints(NDArray(DType::Float64, {8, 3}));
    m.AddPolygonBlock("polygon", {{0, 1, 2}, {1, 2, 3, 4, 5}});
    m.AddPolyhedronBlock("polyhedron", {{{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}}});

    const auto poly = m.Cells(0);
    EXPECT_TRUE(poly.IsRagged());
    EXPECT_FALSE(poly.IsPolyhedron());
    ASSERT_EQ(poly.NumCells(), 2u);
    ASSERT_EQ(poly.RowSize(0), 3u);
    ASSERT_EQ(poly.RowSize(1), 5u);
    EXPECT_EQ(poly.Row(1)[4], 5);

    const auto poh = m.Cells(1);
    EXPECT_TRUE(poh.IsRagged());
    EXPECT_TRUE(poh.IsPolyhedron());
    ASSERT_EQ(poh.NumCells(), 1u);
    ASSERT_EQ(poh.NumFaces(0), 4u);
    const auto face = poh.Face(0, 3);
    ASSERT_EQ(face.second, 3u);
    EXPECT_EQ(face.first[1], 2);
}

TEST(MeshApi, DataNamesStaySortedAfterLaterInserts) {
    // The name lists are memoized (v9.0.0). Reading them BEFORE each insert is
    // the point: a cache that failed to invalidate would keep returning the
    // earlier answer, silently changing on-disk field order rather than
    // crashing. Every backend must pass this identically.
    Mesh m = mt::tri_mesh();
    m.AddPointData("zeta", float64_array({1, 2, 3, 4}));
    EXPECT_EQ(m.PointDataNames(), (std::vector<std::string>{"zeta"}));
    m.AddPointData("alpha", float64_array({4, 3, 2, 1}));
    EXPECT_EQ(m.PointDataNames(), (std::vector<std::string>{"alpha", "zeta"}));
    // Overwriting an existing name must NOT disturb the order.
    m.AddPointData("zeta", float64_array({9, 9, 9, 9}));
    EXPECT_EQ(m.PointDataNames(), (std::vector<std::string>{"alpha", "zeta"}));

    m.AddFieldData("mu", float64_array({7}));
    EXPECT_EQ(m.FieldDataNames(), (std::vector<std::string>{"mu"}));
    m.AddFieldData("beta", float64_array({8}));
    EXPECT_EQ(m.FieldDataNames(), (std::vector<std::string>{"beta", "mu"}));

    m.AddCellData("tags", {int32_array({1, 2})});
    EXPECT_EQ(m.CellDataNames(), (std::vector<std::string>{"tags"}));
    // AppendCellData creates a new list too, so it must invalidate as well.
    m.AppendCellData("area", float64_array({0.5, 0.5}));
    EXPECT_EQ(m.CellDataNames(), (std::vector<std::string>{"area", "tags"}));
}

TEST(MeshApi, DataNamesAreSorted) {
    Mesh m = mt::tri_mesh();
    m.AddPointData("zeta", float64_array({1, 2, 3, 4}));
    m.AddPointData("alpha", float64_array({4, 3, 2, 1}));
    m.AddFieldData("mu", float64_array({7}));
    m.AddFieldData("beta", float64_array({8}));
    m.AddCellData("tags", {int32_array({1, 2})});
    m.AddCellData("area", {float64_array({0.5, 0.5})});

    EXPECT_EQ(m.PointDataNames(), (std::vector<std::string>{"alpha", "zeta"}));
    EXPECT_EQ(m.FieldDataNames(), (std::vector<std::string>{"beta", "mu"}));
    EXPECT_EQ(m.CellDataNames(), (std::vector<std::string>{"area", "tags"}));
    EXPECT_EQ(m.NumPointData(), 2u);
    EXPECT_EQ(m.NumFieldData(), 2u);
    EXPECT_EQ(m.NumCellData(), 2u);
    EXPECT_TRUE(m.HasPointData("alpha"));
    EXPECT_FALSE(m.HasPointData("nope"));
    EXPECT_TRUE(m.HasFieldData("mu"));
    EXPECT_TRUE(m.HasCellData("tags"));
    EXPECT_EQ(m.CellDataNumBlocks("tags"), 1u);

    EXPECT_EQ(meshioplusplus::detail::read_double(m.PointData("zeta"), 1), 2.0);
    EXPECT_EQ(meshioplusplus::detail::read_double(m.FieldData("beta"), 0), 8.0);
    EXPECT_EQ(meshioplusplus::detail::read_int(m.CellData("tags", 0), 1), 2);
}

TEST(MeshApi, IntegerKindSurvivesIngest) {
    // The "first integer cell_data array is the tag" convention (su2, medit,
    // ...) requires that integer-kind data never silently becomes float —
    // regardless of how a backend canonicalizes dtypes internally.
    Mesh m = mt::tri_mesh();
    m.AddCellData("tag", {int32_array({7, 9})});
    m.AddPointData("temp", float64_array({0, 1, 2, 3}));
    EXPECT_FALSE(meshioplusplus::detail::is_float_dtype(m.CellData("tag", 0).Dtype()));
    EXPECT_TRUE(meshioplusplus::detail::is_float_dtype(m.PointData("temp").Dtype()));
    EXPECT_EQ(meshioplusplus::detail::read_int(m.CellData("tag", 0), 0), 7);
}

TEST(MeshApi, AppendCellDataBuildsPerBlockLists) {
    Mesh m = mt::tri_quad_mesh();
    m.AppendCellData("ref", int32_array({1, 1}));
    m.AppendCellData("ref", int32_array({2}));
    m.AppendCellData("ref", int32_array({3}));
    ASSERT_EQ(m.CellDataNumBlocks("ref"), 3u);
    EXPECT_EQ(meshioplusplus::detail::read_int(m.CellData("ref", 1), 0), 2);
    EXPECT_EQ(meshioplusplus::detail::read_int(m.CellData("ref", 2), 0), 3);
}

TEST(MeshApi, AddReplacesExistingName) {
    Mesh m = mt::tri_mesh();
    m.AddPointData("v", float64_array({1, 1, 1, 1}));
    m.AddPointData("v", float64_array({2, 2, 2, 2}));
    EXPECT_EQ(m.NumPointData(), 1u);
    EXPECT_EQ(meshioplusplus::detail::read_double(m.PointData("v"), 0), 2.0);
}
