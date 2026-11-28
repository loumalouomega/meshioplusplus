#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "meshioplusplus/detail/crs.hpp"
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/coordinates.hpp"
#include "meshioplusplus/operations/mesh_boolean.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/units.hpp"
#include "meshioplusplus/providers.hpp"

using namespace meshioplusplus;

namespace {

Mesh make_unit_tet() {
    Mesh m;
    NDArray pts = NDArray::Uninit(DType::Float64, {4, 3});
    auto* p = pts.As<double>();
    p[0] = 0.0;
    p[1] = 0.0;
    p[2] = 0.0;
    p[3] = 1000.0;
    p[4] = 0.0;
    p[5] = 0.0;
    p[6] = 0.0;
    p[7] = 1000.0;
    p[8] = 0.0;
    p[9] = 0.0;
    p[10] = 0.0;
    p[11] = 1000.0;
    m.AssignPoints(std::move(pts));
    NDArray conn = NDArray::Uninit(DType::Int64, {1, 4});
    auto* c = conn.As<std::int64_t>();
    c[0] = 0;
    c[1] = 1;
    c[2] = 2;
    c[3] = 3;
    m.AddCellBlock("tetra", std::move(conn));
    return m;
}

Mesh make_quad_grid_3x3() {
    Mesh m;
    NDArray pts = NDArray::Uninit(DType::Float64, {16, 3});
    auto* p = pts.As<double>();
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            int idx = j * 4 + i;
            p[idx * 3 + 0] = static_cast<double>(i);
            p[idx * 3 + 1] = static_cast<double>(j);
            p[idx * 3 + 2] = 0.0;
        }
    }
    m.AssignPoints(std::move(pts));
    NDArray conn = NDArray::Uninit(DType::Int64, {9, 4});
    auto* c = conn.As<std::int64_t>();
    int q = 0;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            c[q * 4 + 0] = j * 4 + i;
            c[q * 4 + 1] = j * 4 + (i + 1);
            c[q * 4 + 2] = (j + 1) * 4 + (i + 1);
            c[q * 4 + 3] = (j + 1) * 4 + i;
            ++q;
        }
    }
    m.AddCellBlock("quad", std::move(conn));
    return m;
}

Mesh make_closed_tet_surface() {
    Mesh m;
    NDArray pts = NDArray::Uninit(DType::Float64, {4, 3});
    auto* p = pts.As<double>();
    p[0] = 0.0;
    p[1] = 0.0;
    p[2] = 0.0;
    p[3] = 1.0;
    p[4] = 0.0;
    p[5] = 0.0;
    p[6] = 0.0;
    p[7] = 1.0;
    p[8] = 0.0;
    p[9] = 0.0;
    p[10] = 0.0;
    p[11] = 1.0;
    m.AssignPoints(std::move(pts));
    NDArray conn = NDArray::Uninit(DType::Int64, {4, 3});
    auto* c = conn.As<std::int64_t>();
    std::int64_t tris[4][3] = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
    for (int i = 0; i < 4; ++i) {
        c[i * 3 + 0] = tris[i][0];
        c[i * 3 + 1] = tris[i][1];
        c[i * 3 + 2] = tris[i][2];
    }
    m.AddCellBlock("triangle", std::move(conn));
    return m;
}

BooleanProvider make_mock_boolean_provider() {
    BooleanProvider p;
    p.mName = "mock_bool";
    p.mVersion = "1.0";
    p.mBoolean = [](const BooleanRequest& rReq) {
        BooleanReply reply;
        if (rReq.mOp == 1) {
            reply.mPoints = rReq.mPointsB;
            reply.mTriangles = rReq.mTrianglesB;
        } else {
            reply.mPoints = rReq.mPointsA;
            reply.mTriangles = rReq.mTrianglesA;
        }
        return reply;
    };
    return p;
}

}  // namespace

TEST(Wave3Ops, ConvertUnitsRoundTripAndStrict) {
    Mesh m = make_unit_tet();
    detail::set_mesh_length_unit(m, 1e-3);  // mm

    NDArray stress = NDArray::Uninit(DType::Float64, {4});
    for (std::size_t i = 0; i < 4; ++i)
        stress.As<double>()[i] = 250.0;  // 250 MPa
    m.AddPointData("stress", std::move(stress));
    set_field_unit(m, "stress", field_unit_from_name("MPa"));

    FieldUnit read_back;
    ASSERT_TRUE(get_field_unit(m, "stress", read_back));
    EXPECT_DOUBLE_EQ(read_back.mScale, 1e6);

    ConvertUnitsOptions to_si;
    to_si.mSystem = "SI";
    to_si.mStrict = true;
    auto si_res = convert_units(m, to_si);
    EXPECT_EQ(si_res.mNumArraysConverted, 1);
    EXPECT_NEAR(si_res.mMesh.Points().As<double>()[3], 1.0, 1e-12);
    EXPECT_NEAR(si_res.mMesh.PointData("stress").As<double>()[0], 250e6, 1e-3);

    ConvertUnitsOptions back_mmts;
    back_mmts.mSystem = "mm-t-s";
    back_mmts.mStrict = true;
    auto back = convert_units(si_res.mMesh, back_mmts);
    EXPECT_NEAR(back.mMesh.Points().As<double>()[3], 1000.0, 1e-9);
    EXPECT_NEAR(back.mMesh.PointData("stress").As<double>()[0], 250.0, 1e-9);

    // Undeclared array: strict throws, non-strict warns and keeps values.
    NDArray unknown = NDArray::Uninit(DType::Float64, {4});
    unknown.As<double>()[0] = 42.0;
    m.AddPointData("mystery", std::move(unknown));
    EXPECT_THROW(convert_units(m, to_si), std::invalid_argument);

    to_si.mStrict = false;
    auto relaxed = convert_units(m, to_si);
    EXPECT_FALSE(relaxed.mWarnings.empty());
    EXPECT_DOUBLE_EQ(relaxed.mMesh.PointData("mystery").As<double>()[0], 42.0);
}

TEST(Wave3Ops, GeoreferencingCrs) {
    Mesh m = make_unit_tet();
    detail::set_mesh_crs(m, "EPSG:4326");
    EXPECT_EQ(detail::mesh_crs(m), "EPSG:4326");

    CoordinateOptions opts;
    opts.mSystem = CoordinateSystem::Crs;
    opts.mTargetCrs = "EPSG:4326";
    auto same = convert_coordinates(m, opts);
    EXPECT_EQ(detail::mesh_crs(same.mMesh), "EPSG:4326");

    opts.mTargetCrs = "EPSG:3857";
    EXPECT_THROW(convert_coordinates(m, opts), Unsupported);
}

TEST(Wave3Ops, MeshBooleanProviderDispatch) {
    Mesh a = make_closed_tet_surface();
    Mesh b = make_closed_tet_surface();

    EXPECT_FALSE(has_boolean_provider());
    EXPECT_THROW(mesh_boolean(a, b), Unsupported);

    register_boolean_provider(make_mock_boolean_provider());
    EXPECT_TRUE(has_boolean_provider());

    MeshBooleanOptions opts;
    opts.mOp = BooleanOp::Union;
    auto res = mesh_boolean(a, b, opts);
    EXPECT_EQ(res.mBackend, "mock_bool");
    EXPECT_EQ(res.mNumTriangles, 4);
    EXPECT_EQ(res.mMesh.Cells(0).NumCells(), 4u);

    unregister_boolean_provider();
    EXPECT_FALSE(has_boolean_provider());
}

TEST(Wave3Ops, RefineMixedGreenQuad) {
    Mesh m = make_quad_grid_3x3();
    RefineOptions opts;
    opts.mClosure = RefineClosure::MixedGreen;
    opts.mCells = {4};  // center quad only
    auto res = refine(m, opts);

    // Center quad -> 4 quads; 4 corner quads stay untouched -> 8 quads total.
    // 4 edge-adjacent neighbor quads each have 1 split edge (5 perimeter segments) -> centroid fan
    // produces 5 triangles each = 20 triangles and +4 face centers (16 + 4 + 1 + 4 = 25 points).
    std::size_t n_quads = 0;
    std::size_t n_tris = 0;
    for (std::size_t bi = 0; bi < res.mMesh.NumCellBlocks(); ++bi) {
        const auto& cb = res.mMesh.Cells(bi);
        if (cb.Type() == "quad")
            n_quads += cb.NumCells();
        else if (cb.Type() == "triangle")
            n_tris += cb.NumCells();
    }
    EXPECT_EQ(n_quads, 8u);
    EXPECT_EQ(n_tris, 20u);
    EXPECT_EQ(res.mMesh.NumPoints(), 16u + 4u + 1u + 4u);
}
