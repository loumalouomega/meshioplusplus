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

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/stats.hpp"

namespace {

using meshioplusplus::compute_stats;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::smooth;
using meshioplusplus::SmoothMethod;
using meshioplusplus::SmoothOptions;
using meshioplusplus::SmoothResult;
using meshioplusplus::StatsReport;

// --- fixtures ---------------------------------------------------------------

// An n x n regular grid of unit quads in the z = 0 plane.
Mesh quad_grid(int n) {
    std::vector<std::vector<double>> pts;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            pts.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});
    auto pid = [n](int i, int j) { return static_cast<std::int64_t>(i * n + j); };
    std::vector<std::vector<std::int64_t>> cells;
    for (int i = 0; i < n - 1; ++i)
        for (int j = 0; j < n - 1; ++j)
            cells.push_back({pid(i, j), pid(i + 1, j), pid(i + 1, j + 1), pid(i, j + 1)});
    return mt::make_mesh(std::move(pts), "quad", std::move(cells));
}

// An n x n x n block of unit hexahedra.
Mesh hex_block(int n) {
    std::vector<std::vector<double>> pts;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            for (int k = 0; k <= n; ++k)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    auto pid = [n](int i, int j, int k) {
        return static_cast<std::int64_t>((i * (n + 1) + j) * (n + 1) + k);
    };
    std::vector<std::vector<std::int64_t>> cells;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                cells.push_back({pid(i, j, k), pid(i + 1, j, k), pid(i + 1, j + 1, k),
                                 pid(i, j + 1, k), pid(i, j, k + 1), pid(i + 1, j, k + 1),
                                 pid(i + 1, j + 1, k + 1), pid(i, j + 1, k + 1)});
    return mt::make_mesh(std::move(pts), "hexahedron", std::move(cells));
}

double point_at(const Mesh& rMesh, std::size_t Node, std::size_t Comp) {
    return meshioplusplus::detail::read_double(rMesh.Points(), Node * rMesh.PointDim() + Comp);
}

// A cheap deterministic jitter, so the fixtures need no RNG.
void jitter(Mesh& rMesh, const std::vector<std::size_t>& rNodes, double Scale) {
    NDArray pts = rMesh.Points();
    pts.MakeOwned();
    double* p = pts.As<double>();
    const std::size_t dim = rMesh.PointDim();
    for (std::size_t k = 0; k < rNodes.size(); ++k)
        for (std::size_t d = 0; d < dim; ++d)
            p[rNodes[k] * dim + d] +=
                Scale * std::sin(static_cast<double>(7 * k + 13 * d + 1) * 1.7);
    rMesh.AssignPoints(std::move(pts));
}

std::vector<std::size_t> grid_interior(int n) {
    std::vector<std::size_t> out;
    for (int i = 1; i < n - 1; ++i)
        for (int j = 1; j < n - 1; ++j)
            out.push_back(static_cast<std::size_t>(i * n + j));
    return out;
}

SmoothOptions opts(SmoothMethod method, int iterations, bool guard = true) {
    SmoothOptions o;
    o.mMethod = method;
    o.mIterations = iterations;
    o.mGuardInversion = guard;
    return o;
}

// --- the core promise -------------------------------------------------------

TEST(Smooth, DisplacedInteriorNodeReturnsTowardItsCentroid) {
    Mesh mesh = quad_grid(9);
    const std::size_t node = 4 * 9 + 4;
    NDArray pts = mesh.Points();
    pts.MakeOwned();
    pts.As<double>()[node * 3 + 0] += 0.45;
    pts.As<double>()[node * 3 + 1] -= 0.35;
    mesh.AssignPoints(std::move(pts));

    const double before =
        std::hypot(point_at(mesh, node, 0) - 4.0, point_at(mesh, node, 1) - 4.0);
    const SmoothResult r = smooth(mesh, opts(SmoothMethod::Laplacian, 8));
    const double after =
        std::hypot(point_at(r.mMesh, node, 0) - 4.0, point_at(r.mMesh, node, 1) - 4.0);

    EXPECT_GT(before, 0.5);
    EXPECT_LT(after, before / 10.0);
}

TEST(Smooth, TaubinDoesNotShrinkWhereLaplacianDoes) {
    Mesh base = quad_grid(9);
    jitter(base, grid_interior(9), 0.18);

    SmoothOptions lap = opts(SmoothMethod::Laplacian, 40, /*guard=*/false);
    lap.mFixBoundary = false;
    lap.mPreserveFeatures = false;
    SmoothOptions tau = opts(SmoothMethod::Taubin, 40, /*guard=*/false);
    tau.mFixBoundary = false;
    tau.mPreserveFeatures = false;

    const StatsReport s0 = compute_stats(base);
    const StatsReport sl = compute_stats(smooth(base, lap).mMesh);
    const StatsReport st = compute_stats(smooth(base, tau).mMesh);

    // Laplacian collapses the footprint; Taubin holds it.
    EXPECT_LT(sl.mExtent[0], 0.6 * s0.mExtent[0]);
    EXPECT_GT(st.mExtent[0], 0.95 * s0.mExtent[0]);
}

TEST(Smooth, FixBoundaryPinsBoundaryNodesExactly) {
    Mesh mesh = quad_grid(9);
    jitter(mesh, grid_interior(9), 0.15);
    const SmoothResult r = smooth(mesh, opts(SmoothMethod::Taubin, 12));

    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j) {
            if (i != 0 && i != 8 && j != 0 && j != 8)
                continue;
            const std::size_t node = static_cast<std::size_t>(i * 9 + j);
            // Bit-identical: pinned means pinned.
            EXPECT_EQ(point_at(mesh, node, 0), point_at(r.mMesh, node, 0));
            EXPECT_EQ(point_at(mesh, node, 1), point_at(r.mMesh, node, 1));
        }
}

TEST(Smooth, GeometryOnlyTopologyAndDataUnchanged) {
    Mesh mesh = mt::data_mesh();
    const SmoothResult r = smooth(mesh, opts(SmoothMethod::Taubin, 5));
    // Cells, cell_data, field_data and point_data values all pass through; only
    // the coordinates may differ.
    mt::expect_same_geometry(mesh, r.mMesh);
}

TEST(Smooth, PointCountAndBlockStructurePreserved) {
    Mesh mesh = hex_block(3);
    jitter(mesh, {31, 32, 33}, 0.2);
    const SmoothResult r = smooth(mesh, opts(SmoothMethod::Taubin, 6));

    EXPECT_EQ(r.mMesh.NumPoints(), mesh.NumPoints());
    EXPECT_EQ(r.mMesh.NumCellBlocks(), mesh.NumCellBlocks());
    EXPECT_EQ(mt::cell_rows(r.mMesh), mt::cell_rows(mesh));
}

// --- the inversion guard ----------------------------------------------------

TEST(Smooth, GuardPreventsNewInversions) {
    Mesh mesh = hex_block(4);
    std::vector<std::size_t> interior;
    for (std::size_t i = 0; i < mesh.NumPoints(); ++i) {
        const double x = point_at(mesh, i, 0);
        const double y = point_at(mesh, i, 1);
        const double z = point_at(mesh, i, 2);
        if (x > 0 && x < 4 && y > 0 && y < 4 && z > 0 && z < 4)
            interior.push_back(i);
    }
    jitter(mesh, interior, 0.75);

    const SmoothResult r = smooth(mesh, opts(SmoothMethod::Taubin, 15, /*guard=*/true));
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
}

TEST(Smooth, GuardDoesNotLockInPreExistingInversions) {
    // The guard is "do no harm", not "preserve the sign": a cell that arrives
    // inverted must still be allowed to be repaired. An earlier version rejected
    // any sign change and so pinned every pre-existing inversion permanently.
    // Drag one interior node clean through the plane of its neighbours, which
    // turns every hexahedron touching it inside out.
    Mesh mesh = hex_block(3);
    const std::size_t tangled = (1 * 4 + 1) * 4 + 1;  // the node at (1, 1, 1)
    NDArray pts = mesh.Points();
    pts.MakeOwned();
    pts.As<double>()[tangled * 3 + 0] += 3.5;
    mesh.AssignPoints(std::move(pts));
    ASSERT_GT(compute_stats(mesh).mNumInverted, 0) << "fixture must start tangled";

    const SmoothResult r = smooth(mesh, opts(SmoothMethod::Taubin, 30, /*guard=*/true));
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
}

TEST(Smooth, SmoothingNeverIncreasesTheInvertedCount) {
    // Guards against smooth_facefan_volume and quality.cpp's
    // quality_facefan_volume (a deliberate duplicate) drifting apart in sign.
    for (double scale : {0.1, 0.3, 0.5}) {
        Mesh mesh = hex_block(3);
        std::vector<std::size_t> all;
        for (std::size_t i = 0; i < mesh.NumPoints(); ++i)
            all.push_back(i);
        jitter(mesh, all, scale);
        const std::int64_t before = compute_stats(mesh).mNumInverted;
        const std::int64_t after =
            compute_stats(smooth(mesh, opts(SmoothMethod::Taubin, 10)).mMesh).mNumInverted;
        EXPECT_LE(after, before) << "scale = " << scale;
    }
}

TEST(Smooth, QualityImproves) {
    Mesh mesh = hex_block(4);
    std::vector<std::size_t> interior;
    for (std::size_t i = 0; i < mesh.NumPoints(); ++i) {
        const double x = point_at(mesh, i, 0);
        const double y = point_at(mesh, i, 1);
        const double z = point_at(mesh, i, 2);
        if (x > 0 && x < 4 && y > 0 && y < 4 && z > 0 && z < 4)
            interior.push_back(i);
    }
    jitter(mesh, interior, 0.35);

    auto min_sj = [](const Mesh& rM) {
        for (const auto& kv : meshioplusplus::compute_quality(rM).mMetrics)
            if (kv.first == "quality:scaled_jacobian")
                return kv.second.mMin;
        return 0.0;
    };
    const double before = min_sj(mesh);
    const double after = min_sj(smooth(mesh, opts(SmoothMethod::Taubin, 20)).mMesh);
    EXPECT_GT(after, before);
}

// --- pinning ----------------------------------------------------------------

TEST(Smooth, FrozenMaskPinsExtraNodes) {
    Mesh mesh = quad_grid(9);
    jitter(mesh, grid_interior(9), 0.15);
    SmoothOptions o = opts(SmoothMethod::Taubin, 10);
    o.mFrozen.assign(mesh.NumPoints(), 0);
    const std::size_t pinned = 4 * 9 + 4;
    o.mFrozen[pinned] = 1;

    const SmoothResult r = smooth(mesh, o);
    EXPECT_EQ(point_at(mesh, pinned, 0), point_at(r.mMesh, pinned, 0));
    EXPECT_EQ(point_at(mesh, pinned, 1), point_at(r.mMesh, pinned, 1));
}

TEST(Smooth, MissizedFrozenMaskThrows) {
    Mesh mesh = quad_grid(5);
    SmoothOptions o;
    o.mFrozen.assign(3, 0);
    EXPECT_THROW(smooth(mesh, o), std::invalid_argument);
}

TEST(Smooth, HigherOrderCellsArePinnedNotDistorted) {
    // triangle6 has no cell_refine_edges row, so its nodes hold still rather
    // than being smoothed toward a guessed neighbourhood.
    Mesh mesh = mt::triangle6_mesh();
    SmoothOptions o = opts(SmoothMethod::Taubin, 10, /*guard=*/false);
    const SmoothResult r = smooth(mesh, o);
    for (std::size_t i = 0; i < mesh.NumPoints(); ++i)
        for (std::size_t d = 0; d < mesh.PointDim(); ++d)
            EXPECT_EQ(point_at(mesh, i, d), point_at(r.mMesh, i, d));
    EXPECT_EQ(r.mNumNodesMoved, 0);
}

// --- validation -------------------------------------------------------------

TEST(Smooth, InvalidOptionsThrow) {
    Mesh mesh = quad_grid(5);
    SmoothOptions bad_lambda;
    bad_lambda.mLambda = 1.5;
    EXPECT_THROW(smooth(mesh, bad_lambda), std::invalid_argument);

    SmoothOptions bad_mu;
    bad_mu.mMethod = SmoothMethod::Taubin;
    bad_mu.mLambda = 0.4;
    bad_mu.mMu = -0.2;  // must be < -lambda
    EXPECT_THROW(smooth(mesh, bad_mu), std::invalid_argument);

    EXPECT_THROW(meshioplusplus::smooth_method_from_name("bogus"), std::invalid_argument);
    EXPECT_EQ(meshioplusplus::smooth_method_from_name("taubin"), SmoothMethod::Taubin);
    EXPECT_EQ(meshioplusplus::smooth_method_from_name("laplacian"), SmoothMethod::Laplacian);
}

TEST(Smooth, ZeroIterationsIsANoOp) {
    Mesh mesh = quad_grid(6);
    jitter(mesh, grid_interior(6), 0.2);
    const SmoothResult r = smooth(mesh, opts(SmoothMethod::Taubin, 0));
    for (std::size_t i = 0; i < mesh.NumPoints(); ++i)
        for (std::size_t d = 0; d < mesh.PointDim(); ++d)
            EXPECT_EQ(point_at(mesh, i, d), point_at(r.mMesh, i, d));
    EXPECT_EQ(r.mNumNodesMoved, 0);
    EXPECT_EQ(r.mMaxDisplacement, 0.0);
}

TEST(Smooth, EmptyMeshIsHandled) {
    Mesh mesh;
    const SmoothResult r = smooth(mesh, opts(SmoothMethod::Taubin, 3));
    EXPECT_EQ(r.mMesh.NumPoints(), 0u);
    EXPECT_EQ(r.mNumNodesMoved, 0);
}

// --- determinism ------------------------------------------------------------

TEST(Smooth, ResultIsStableAcrossRepeatedRuns) {
    Mesh mesh = hex_block(3);
    std::vector<std::size_t> all;
    for (std::size_t i = 0; i < mesh.NumPoints(); ++i)
        all.push_back(i);
    jitter(mesh, all, 0.12);

    const SmoothResult first = smooth(mesh, opts(SmoothMethod::Taubin, 7));
    for (int run = 0; run < 5; ++run) {
        const SmoothResult again = smooth(mesh, opts(SmoothMethod::Taubin, 7));
        ASSERT_EQ(again.mMesh.NumPoints(), first.mMesh.NumPoints());
        for (std::size_t i = 0; i < first.mMesh.NumPoints(); ++i)
            for (std::size_t d = 0; d < 3; ++d)
                // Byte-identical, not merely close -- the Jacobi update plus the
                // ascending-neighbour sum order is what guarantees this across
                // backends and thread counts.
                ASSERT_EQ(point_at(first.mMesh, i, d), point_at(again.mMesh, i, d))
                    << "run " << run << ", node " << i << ", component " << d;
        ASSERT_EQ(again.mNumSkippedInversion, first.mNumSkippedInversion);
    }
}

}  // namespace
