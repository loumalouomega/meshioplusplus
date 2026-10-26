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
// boundary_layer (roadmap 6.9): prisms grown from a surface along its most
// normal normals, cut short where they would collide, never written folded.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

// Project includes
#include "brep_solids.hpp"
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "../../src/cpp/src/detail/layer_sweep.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/boundary_layer.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/quality.hpp"

namespace {

namespace mio = meshioplusplus;
namespace mdt = meshioplusplus::detail;

constexpr double kPi = 3.14159265358979323846;

std::uint64_t digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

mio::Mesh skin(const mio::BrepInfo& rB, double H, mio::MeshCadQuads Q = mio::MeshCadQuads::None) {
    mio::MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    mio::MeshCadSurfaceExtras ex;
    ex.mQuads = Q;
    return mio::mesh_cad_surface(rB, o, mio::MeshCadSizeSources{}, ex).mMesh;
}

mio::BoundaryLayerOptions layers(double H, int N, double G = 1.2) {
    mio::BoundaryLayerOptions o;
    o.mFirstHeight = H;
    o.mNumLayers = N;
    o.mGrowth = G;
    return o;
}

/// The volume a block's closed faces enclose (divergence theorem); quads
/// by their two triangles from corner 0.
double enclosed(const mio::Mesh& rM, std::size_t Block) {
    const double* x = rM.Points().As<double>();
    const auto cells = rM.Cells(Block);
    const std::size_t k = cells.NodesPerCell();
    double v = 0.0;
    for (std::size_t c = 0; c < cells.NumCells(); ++c) {
        std::array<std::int64_t, 4> p{};
        for (std::size_t i = 0; i < k; ++i)
            p[i] = mio::detail::read_int(cells.Conn(), k * c + i);
        const auto tet = [&](std::int64_t A, std::int64_t B, std::int64_t C) {
            const double* a = x + 3 * A;
            const double* b = x + 3 * B;
            const double* d = x + 3 * C;
            v += (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
                  a[2] * (b[0] * d[1] - b[1] * d[0])) /
                 6.0;
        };
        tet(p[0], p[1], p[2]);
        if (k == 4)
            tet(p[0], p[2], p[3]);
    }
    return v;
}

/// compute_quality's sum of the volume cells' volumes and their least
/// scaled Jacobian.
std::pair<double, double> cells_volume(const mio::Mesh& rM) {
    const mio::QualityReport q = mio::compute_quality(rM);
    double vol = 0.0, sj = 1.0;
    for (const auto& [name, r_blocks] : q.mCellArrays) {
        if (name != "quality:volume" && name != "quality:scaled_jacobian")
            continue;
        for (std::size_t b = 0; b < r_blocks.size(); ++b) {
            const std::string type(rM.Cells(b).Type());
            if (type != "wedge" && type != "hexahedron" && type != "pyramid" && type != "tetra")
                continue;
            const double* v = r_blocks[b].As<double>();
            for (std::size_t c = 0; c < rM.Cells(b).NumCells(); ++c) {
                if (name == "quality:volume")
                    vol += v[c];
                else
                    sj = std::min(sj, v[c]);
            }
        }
    }
    return {vol, sj};
}

}  // namespace

// The unit cube grown outward: at its edges and corners the most normal
// normal ((1, 1, 1) / sqrt 3 at a corner) and its stretch keep every face's
// layers their thickness, so the cap lies on the planes at the total
// thickness off the faces; the cells fill the shell between wall and cap.
TEST(BoundaryLayer, TheCubesCapKeepsItsThickness) {
    const mio::Mesh s = skin(brep_solids::box(), 0.25);
    const auto o = layers(0.01, 4, 1.3);
    const double total = 0.01 * (1.0 + 1.3 + 1.3 * 1.3 + 1.3 * 1.3 * 1.3);
    const auto r = mio::boundary_layer(s, o);
    EXPECT_EQ(r.mNumCollapsed, 0);
    EXPECT_EQ(r.mNumShrunk, 0);
    EXPECT_TRUE(r.mCapClosed);
    EXPECT_GT(r.mMinScaledJacobian, 0.1);
    EXPECT_NEAR(r.mMinThickness, total, 1e-12);
    // Every cap point lies on the cube grown by the total thickness.
    const double* x = r.mMesh.Points().As<double>();
    std::size_t cap = 0, wall = 0;
    for (std::size_t b = 0; b < r.mMesh.NumCellBlocks(); ++b) {
        const auto& r_l = r.mMesh.CellData(mio::kBoundaryLayerName, b);
        if (mio::detail::read_int(r_l, 0) == o.mNumLayers + 1)
            cap = b;
        if (mio::detail::read_int(r_l, 0) == 0)
            wall = b;
    }
    const auto cells = r.mMesh.Cells(cap);
    for (std::size_t c = 0; c < cells.NumCells(); ++c)
        for (std::size_t i = 0; i < 3; ++i) {
            const auto p = static_cast<std::size_t>(mio::detail::read_int(cells.Conn(), 3 * c + i));
            double far = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
                far = std::max({far, -x[3 * p + k], x[3 * p + k] - 1.0});
            EXPECT_NEAR(far, total, 1e-12) << p;
        }
    const auto [vol, sj] = cells_volume(r.mMesh);
    EXPECT_NEAR(vol, enclosed(r.mMesh, cap) - enclosed(r.mMesh, wall), 1e-12);
    EXPECT_GT(sj, 0.1);
}

// The sphere grown outward and inward: closed caps, the cells filling the
// shell exactly, every cell valid, and the result independent of threads.
TEST(BoundaryLayer, TheSphereGrowsBothWays) {
    const mio::Mesh s = skin(brep_solids::sphere(1.0), 0.2);
    for (const bool reverse : {false, true}) {
        auto o = layers(0.02, 3);
        o.mReverse = reverse;
        const auto r = mio::boundary_layer(s, o);
        std::printf("[6.9] sphere %s: %lld wedges, least scaled Jacobian %.3f, thickness %.4f\n",
                    reverse ? "inward" : "outward", static_cast<long long>(r.mNumWedges),
                    r.mMinScaledJacobian, r.mMinThickness);
        EXPECT_EQ(r.mNumCollapsed, 0) << reverse;
        EXPECT_EQ(r.mNumShrunk, 0) << reverse;
        EXPECT_TRUE(r.mCapClosed);
        EXPECT_GT(r.mMinScaledJacobian, 0.1);
        const auto [vol, sj] = cells_volume(r.mMesh);
        EXPECT_GT(sj, 0.1);
        // Thickness 0.02 (1 + 1.2 + 1.44): the shell's volume.
        const double t = 0.02 * 3.64;
        const double shell = reverse ? 4.0 / 3.0 * kPi * (1.0 - std::pow(1.0 - t, 3))
                                     : 4.0 / 3.0 * kPi * (std::pow(1.0 + t, 3) - 1.0);
        EXPECT_NEAR(vol / shell, 1.0, 0.03) << reverse;
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
        const std::uint64_t ref = digest(r.mMesh);
        const int saved = omp_get_max_threads();
        for (int threads : {1, 3}) {
            omp_set_num_threads(threads);
            EXPECT_EQ(digest(mio::boundary_layer(s, o).mMesh), ref) << threads;
        }
        omp_set_num_threads(saved);
#endif
    }
}

// A slab 0.05 thick grown inward 0.1: the layers growing toward each other
// stop at 0.45 of the gap, valid, the cap not crossing the middle.
TEST(BoundaryLayer, LayersStopBeforeTheyMeet) {
    const mio::Mesh s = skin(brep_solids::box(1.0, 1.0, 0.05), 0.1);
    auto o = layers(0.02, 4, 1.0);
    o.mReverse = true;
    const auto r = mio::boundary_layer(s, o);
    std::printf("[6.9] slab: %lld of the points thinner, thickness %.4f to %.4f\n",
                static_cast<long long>(r.mNumShrunk), r.mMinThickness, r.mMaxThickness);
    EXPECT_GT(r.mNumShrunk, 0);
    EXPECT_LE(r.mMaxThickness, 0.45 * 0.05 + 1e-12);
    EXPECT_GT(r.mMinScaledJacobian, 0.0);
    const auto [vol, sj] = cells_volume(r.mMesh);
    EXPECT_LT(vol, 0.05);  // inside the slab
    EXPECT_GT(sj, 0.0);
}

// An all-quad cube grows hexahedra; the face data is carried to the cells.
TEST(BoundaryLayer, QuadsGrowHexahedra) {
    const mio::Mesh s = skin(brep_solids::box(), 0.25, mio::MeshCadQuads::All);
    const auto r = mio::boundary_layer(s, layers(0.02, 2));
    EXPECT_GT(r.mNumHexahedra, 0);
    EXPECT_EQ(r.mNumWedges, 0);
    EXPECT_TRUE(r.mCapClosed);
    EXPECT_GT(r.mMinScaledJacobian, 0.1);
    ASSERT_TRUE(r.mMesh.HasCellData(mio::kMeshCadFaceName));
    for (std::size_t b = 0; b < r.mMesh.NumCellBlocks(); ++b) {
        const auto& r_f = r.mMesh.CellData(mio::kMeshCadFaceName, b);
        for (std::size_t c = 0; c < r.mMesh.Cells(b).NumCells(); ++c) {
            const auto f = mio::detail::read_int(r_f, c);
            EXPECT_TRUE(f >= 0 && f < 6) << b << " " << c;
        }
    }
}

// The sweep alone: a point whose offsets are zero collapses its cells into
// a pyramid (one corner) and a tetrahedron (two), both positive, and the
// wall and cap faces come out split like their column.
TEST(BoundaryLayer, TheSweepCollapsesAPointIntoPyramidsAndTetrahedra) {
    // Two triangles of the unit square, normals +z; point 0 collapsed.
    const std::vector<double> x = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    mdt::LayerSweepInput in;
    in.pPoints = x.data();
    in.mNumPoints = 4;
    in.mFaces = {{0, 1, 2, -1}, {0, 2, 3, -1}};
    in.mLayers = 2;
    in.mDirections = {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1};
    in.mOffsets = {0, 0, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    auto r = mdt::layer_sweep(in);
    EXPECT_EQ(r.mWedges.size(), 0u);
    EXPECT_EQ(r.mPyramids.size(), 4u);
    for (const auto& r_c : r.mPyramids)
        EXPECT_GT(mdt::layer_cell_scaled_jacobian(r.mPoints.data(), r_c, 5), 0.0);
    // Collapsing point 1 too leaves a tetrahedron per layer of the first
    // triangle.
    in.mOffsets = {0, 0, 0, 0, 0.1, 0.2, 0.1, 0.2};
    r = mdt::layer_sweep(in);
    EXPECT_EQ(r.mTetrahedra.size(), 2u);
    for (const auto& r_c : r.mTetrahedra)
        EXPECT_GT(mdt::layer_cell_scaled_jacobian(r.mPoints.data(), r_c, 4), 0.0);
    // A quad with a collapsed corner is split into two triangles' columns.
    in.mFaces = {{0, 1, 2, 3}};
    in.mOffsets = {0, 0, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    r = mdt::layer_sweep(in);
    EXPECT_EQ(r.mHexahedra.size(), 0u);
    EXPECT_EQ(r.mPyramids.size(), 4u);
    EXPECT_EQ(r.mWallTriangles.size(), 2u);
    EXPECT_EQ(r.mCapTriangles.size(), 2u);
    // Nothing collapsed: a hexahedron per layer, a right one.
    in.mOffsets = {0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2};
    r = mdt::layer_sweep(in);
    ASSERT_EQ(r.mHexahedra.size(), 2u);
    EXPECT_NEAR(mdt::layer_cell_scaled_jacobian(r.mPoints.data(), r.mHexahedra[0], 8), 1.0, 1e-12);
}

// Two triangles folded back onto each other across their shared edge: no
// direction sees both from the edge's points, which stay on the surface;
// the rest grows, and every cell written is valid.
TEST(BoundaryLayer, APointNoDirectionSeesStaysPut) {
    mio::Mesh m;
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {4, 3});
    const double coords[12] = {0, 0, 0, 1, 0, 0, 0.5, 1, 0, 0.5, 1, 0.05};
    std::copy(coords, coords + 12, p.As<double>());
    m.AssignPoints(std::move(p));
    mio::NDArray c = mio::NDArray::Uninit(mio::DType::Int64, {2, 3});
    const std::int64_t tris[6] = {0, 1, 2, 1, 0, 3};
    std::copy(tris, tris + 6, c.As<std::int64_t>());
    m.AddCellBlock("triangle", std::move(c));
    const auto r = mio::boundary_layer(m, layers(0.01, 2));
    EXPECT_EQ(r.mNumNotVisible, 2);
    EXPECT_GE(r.mNumCollapsed, 2);
    EXPECT_EQ(r.mNumWedges, 0);
    for (std::size_t b = 0; b < r.mMesh.NumCellBlocks(); ++b) {
        const std::string type(r.mMesh.Cells(b).Type());
        EXPECT_TRUE(type == "tetra" || type == "pyramid" || type == "triangle") << type;
    }
    if (r.mNumCells > 0)
        EXPECT_GT(r.mMinScaledJacobian, 0.0);
    EXPECT_FALSE(r.mWarnings.empty());
}

TEST(BoundaryLayer, RefusalsAreByName) {
    const mio::Mesh s = skin(brep_solids::box(), 0.5);
    EXPECT_THROW((void)mio::boundary_layer(s, layers(0.0, 3)), std::invalid_argument);
    EXPECT_THROW((void)mio::boundary_layer(s, layers(0.01, 0)), std::invalid_argument);
    EXPECT_THROW((void)mio::boundary_layer(s, layers(0.01, 3, 0.5)), std::invalid_argument);
    auto o = layers(0.01, 3);
    o.mCollisionFraction = 0.0;
    EXPECT_THROW((void)mio::boundary_layer(s, o), std::invalid_argument);
    // A tetrahedral mesh is not a surface.
    mio::Mesh tet;
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {4, 3});
    const double coords[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    std::copy(coords, coords + 12, p.As<double>());
    tet.AssignPoints(std::move(p));
    mio::NDArray c = mio::NDArray::Uninit(mio::DType::Int64, {1, 4});
    for (std::int64_t i = 0; i < 4; ++i)
        c.As<std::int64_t>()[i] = i;
    tet.AddCellBlock("tetra", std::move(c));
    EXPECT_THROW((void)mio::boundary_layer(tet, layers(0.01, 3)), std::invalid_argument);
}
