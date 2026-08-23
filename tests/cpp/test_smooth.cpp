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
#include <array>
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

// An n x n x n block of unit tetrahedra: the same cube6_mesh() 6-tet
// decomposition (shared main diagonal 0-6, the standard FEM split) applied
// per hex cell of hex_block(n)'s own point lattice -- reuses that lattice
// exactly, so neighbouring cells already share vertices with no separate
// welding pass.
Mesh tet_block(int n) {
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
            for (int k = 0; k < n; ++k) {
                const std::int64_t h[8] = {
                    pid(i, j, k),         pid(i + 1, j, k),         pid(i + 1, j + 1, k),
                    pid(i, j + 1, k),     pid(i, j, k + 1),         pid(i + 1, j, k + 1),
                    pid(i + 1, j + 1, k + 1), pid(i, j + 1, k + 1)};
                cells.push_back({h[0], h[1], h[2], h[6]});
                cells.push_back({h[0], h[2], h[3], h[6]});
                cells.push_back({h[0], h[3], h[7], h[6]});
                cells.push_back({h[0], h[7], h[4], h[6]});
                cells.push_back({h[0], h[4], h[5], h[6]});
                cells.push_back({h[0], h[5], h[1], h[6]});
            }
    return mt::make_mesh(std::move(pts), "tetra", std::move(cells));
}

// Node ids of tet_block(n) strictly away from every boundary face (1 <= i,j,k
// <= n-1), the tetrahedral analogue of grid_interior/hex_block's own interior
// selection idiom.
std::vector<std::size_t> tet_block_interior(int n) {
    std::vector<std::size_t> out;
    auto pid = [n](int i, int j, int k) {
        return static_cast<std::size_t>((i * (n + 1) + j) * (n + 1) + k);
    };
    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j)
            for (int k = 1; k < n; ++k)
                out.push_back(pid(i, j, k));
    return out;
}

// A regular octahedron (apexes at unit distance along each axis) split into 8
// tets sharing its own center O. O is strictly interior -- every one of the
// octahedron's 8 outer triangular faces uses only the 6 apex points, never O
// -- so under the default mFixBoundary it is the mesh's only free vertex.
// Winding is deliberately not normalised per tet (nothing here exercises the
// inversion guard); tets use whichever apex order falls out of the loop.
Mesh octahedron_mesh() {
    std::vector<std::vector<double>> pts = {
        {0, 0, 0},    // 0: center O
        {1, 0, 0},   // 1: +x
        {-1, 0, 0},  // 2: -x
        {0, 1, 0},   // 3: +y
        {0, -1, 0},  // 4: -y
        {0, 0, 1},   // 5: +z
        {0, 0, -1},  // 6: -z
    };
    const int xs[2] = {1, 2};
    const int ys[2] = {3, 4};
    const int zs[2] = {5, 6};
    std::vector<std::vector<std::int64_t>> cells;
    for (int xi : xs)
        for (int yi : ys)
            for (int zi : zs)
                cells.push_back({0, xi, yi, zi});
    return mt::make_mesh(std::move(pts), "tetra", std::move(cells));
}

// Independent reference circumcenter of a tetrahedron, the standard explicit
// cross-product closed form -- a genuinely different algebraic derivation
// from smooth.cpp's own (private, hence unreachable from a test TU anyway)
// linear-system Cramer's-rule solve, so the two cannot share a transcription
// bug. Used only to compute expected values below.
std::array<double, 3> ref_circumcenter(const std::array<double, 3>& p0,
                                       const std::array<double, 3>& p1,
                                       const std::array<double, 3>& p2,
                                       const std::array<double, 3>& p3) {
    auto sub = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };
    auto cross = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return std::array<double, 3>{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                                     a[0] * b[1] - a[1] * b[0]};
    };
    auto dot = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    const std::array<double, 3> d1 = sub(p1, p0), d2 = sub(p2, p0), d3 = sub(p3, p0);
    const double n1 = dot(d1, d1), n2 = dot(d2, d2), n3 = dot(d3, d3);
    const std::array<double, 3> c23 = cross(d2, d3), c31 = cross(d3, d1), c12 = cross(d1, d2);
    const double denom = 2.0 * dot(d1, c23);
    return {p0[0] + (n1 * c23[0] + n2 * c31[0] + n3 * c12[0]) / denom,
            p0[1] + (n1 * c23[1] + n2 * c31[1] + n3 * c12[1]) / denom,
            p0[2] + (n1 * c23[2] + n2 * c31[2] + n3 * c12[2]) / denom};
}

// The worst (smallest) per-cell min-dihedral / (largest) per-cell max-dihedral
// across the whole mesh -- the two summary fields `Smooth.*Odt*` tests use to
// judge overall element quality.
double worst_min_dihedral(const Mesh& rM) {
    for (const auto& kv : meshioplusplus::compute_quality(rM).mMetrics)
        if (kv.first == "quality:min_dihedral")
            return kv.second.mMin;
    return std::nan("");
}
double worst_max_dihedral(const Mesh& rM) {
    for (const auto& kv : meshioplusplus::compute_quality(rM).mMetrics)
        if (kv.first == "quality:max_dihedral")
            return kv.second.mMax;
    return std::nan("");
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

    const double before = std::hypot(point_at(mesh, node, 0) - 4.0, point_at(mesh, node, 1) - 4.0);
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
    //
    // The displacement was 3.5 before v9.16.0, when the signed volume fanned
    // each face about its own first node. That decomposition reported those
    // cells as inverted; the corner-average fan (detail/polyhedron.hpp) does
    // not, and neither answer is "wrong" -- a self-intersecting hexahedron has
    // no decomposition-independent signed volume. 5.0 is inverted under the
    // canonical measure, so the fixture states its premise honestly rather than
    // relying on an artefact of one triangulation.
    Mesh mesh = hex_block(3);
    const std::size_t tangled = (1 * 4 + 1) * 4 + 1;  // the node at (1, 1, 1)
    NDArray pts = mesh.Points();
    pts.MakeOwned();
    pts.As<double>()[tangled * 3 + 0] += 5.0;
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
    EXPECT_EQ(meshioplusplus::smooth_method_from_name("odt"), SmoothMethod::Odt);
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

// --- ODT ---------------------------------------------------------------------

TEST(Smooth, OdtMovesEachVertexOfASingleTetTowardItsCircumcenter) {
    // A fully asymmetric, general-position tetrahedron. mFixBoundary is off --
    // every vertex of a lone tet is trivially its own boundary, so leaving the
    // default on would freeze the whole mesh -- and so is the inversion guard,
    // isolating the formula itself from the guard's own logic. Each of the 4
    // vertices sees only this one tet, so its ODT target is exactly that tet's
    // own circumcenter, computed here via ref_circumcenter -- a genuinely
    // independent formula from smooth.cpp's own (private, and so unreachable
    // from a test TU regardless) linear-system solve. This is the strong,
    // discriminating oracle: for a generic asymmetric tetrahedron the
    // circumcenter and the centroid are different points, so a bug that
    // silently computed the wrong one (e.g. a Laplacian-style centroid) would
    // fail this test with a wrong numeric target, unlike a symmetry-only check.
    const std::array<double, 3> p0 = {0.0, 0.0, 0.0};
    const std::array<double, 3> p1 = {3.0, 0.0, 0.0};
    const std::array<double, 3> p2 = {0.0, 2.0, 0.3};
    const std::array<double, 3> p3 = {0.5, 0.7, 4.0};
    Mesh mesh = mt::make_mesh({{p0[0], p0[1], p0[2]},
                               {p1[0], p1[1], p1[2]},
                               {p2[0], p2[1], p2[2]},
                               {p3[0], p3[1], p3[2]}},
                              "tetra", {{0, 1, 2, 3}});

    SmoothOptions o;
    o.mMethod = SmoothMethod::Odt;
    o.mIterations = 1;
    o.mLambda = 0.6;
    o.mFixBoundary = false;
    o.mGuardInversion = false;

    const SmoothResult r = smooth(mesh, o);
    const std::array<double, 3> cc = ref_circumcenter(p0, p1, p2, p3);
    const std::array<double, 3> pts[4] = {p0, p1, p2, p3};
    for (int i = 0; i < 4; ++i)
        for (int d = 0; d < 3; ++d) {
            const double expected = pts[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)] +
                                    o.mLambda *
                                        (cc[static_cast<std::size_t>(d)] -
                                         pts[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)]);
            EXPECT_NEAR(
                point_at(r.mMesh, static_cast<std::size_t>(i), static_cast<std::size_t>(d)),
                expected, 1e-9)
                << "vertex " << i << ", component " << d;
        }
}

TEST(Smooth, OdtInteriorVertexOfASymmetricOctahedronStaysAtTheCenter) {
    // A weaker, WEAK-BY-ITSELF oracle, documented as such: the octahedral
    // symmetry group forces ANY equivariant target function's average over
    // the 8 surrounding tets to land back on the center (e.g. an
    // accidentally-centroid-computing implementation would ALSO pass this
    // test), so it does not by itself pin the circumcenter formula -- that is
    // OdtMovesEachVertexOfASingleTetTowardItsCircumcenter's job. What this
    // test DOES pin, and the other cannot: the node -> incident-tet CSR
    // (`incidence`/`cells`, built and shared with the guard) correctly
    // aggregates MULTIPLE incident tets with the right volume weighting, and
    // boundary pinning correctly leaves exactly one node free.
    Mesh mesh = octahedron_mesh();
    SmoothOptions o;
    o.mMethod = SmoothMethod::Odt;
    o.mIterations = 3;
    o.mGuardInversion = false;  // winding not normalised in this fixture
    const SmoothResult r = smooth(mesh, o);
    for (std::size_t d = 0; d < 3; ++d)
        EXPECT_NEAR(point_at(r.mMesh, 0, d), 0.0, 1e-12);
    // The 6 apex points are all boundary and must not move at all.
    for (std::size_t i = 1; i <= 6; ++i)
        for (std::size_t d = 0; d < 3; ++d)
            EXPECT_EQ(point_at(r.mMesh, i, d), point_at(mesh, i, d));
}

TEST(Smooth, OdtImprovesDihedralAnglesMoreThanTaubinOnAJitteredTetBlock) {
    Mesh mesh = tet_block(3);
    jitter(mesh, tet_block_interior(3), 0.2);

    const double before_min = worst_min_dihedral(mesh);
    const double before_max = worst_max_dihedral(mesh);

    const Mesh taubin_out = smooth(mesh, opts(SmoothMethod::Taubin, 15)).mMesh;
    SmoothOptions odt_opts;
    odt_opts.mMethod = SmoothMethod::Odt;
    odt_opts.mIterations = 15;
    const Mesh odt_out = smooth(mesh, odt_opts).mMesh;

    const double taubin_min = worst_min_dihedral(taubin_out);
    const double taubin_max = worst_max_dihedral(taubin_out);
    const double odt_min = worst_min_dihedral(odt_out);
    const double odt_max = worst_max_dihedral(odt_out);

    EXPECT_GT(odt_min, before_min);
    EXPECT_LT(odt_max, before_max);
    EXPECT_GT(odt_min, taubin_min) << "odt=" << odt_min << " taubin=" << taubin_min;
    EXPECT_LT(odt_max, taubin_max) << "odt=" << odt_max << " taubin=" << taubin_max;
}

TEST(Smooth, OdtGuardPreventsNewInversions) {
    Mesh mesh = tet_block(4);
    jitter(mesh, tet_block_interior(4), 0.75);

    SmoothOptions o;
    o.mMethod = SmoothMethod::Odt;
    o.mIterations = 15;
    o.mGuardInversion = true;
    const SmoothResult r = smooth(mesh, o);
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
}

TEST(Smooth, OdtIsDeterministicAcrossRepeatedRuns) {
    Mesh mesh = tet_block(3);
    jitter(mesh, tet_block_interior(3), 0.15);

    SmoothOptions o;
    o.mMethod = SmoothMethod::Odt;
    o.mIterations = 5;

    const SmoothResult first = smooth(mesh, o);
    for (int run = 0; run < 5; ++run) {
        const SmoothResult again = smooth(mesh, o);
        ASSERT_EQ(again.mMesh.NumPoints(), first.mMesh.NumPoints());
        for (std::size_t i = 0; i < first.mMesh.NumPoints(); ++i)
            for (std::size_t d = 0; d < 3; ++d)
                ASSERT_EQ(point_at(first.mMesh, i, d), point_at(again.mMesh, i, d))
                    << "run " << run << ", node " << i << ", component " << d;
        ASSERT_EQ(again.mNumSkippedInversion, first.mNumSkippedInversion);
    }
}

TEST(Smooth, OdtSilentlyIgnoresMu) {
    Mesh mesh = tet_block(2);
    SmoothOptions o;
    o.mMethod = SmoothMethod::Odt;
    o.mIterations = 3;
    o.mMu = 12.5;  // nonsensical for taubin, but odt never reads mMu at all
    EXPECT_NO_THROW(smooth(mesh, o));
}

TEST(Smooth, OdtRejectsBlocksOutsideTheTetOnlyScope) {
    SmoothOptions odt;
    odt.mMethod = SmoothMethod::Odt;

    EXPECT_THROW(smooth(hex_block(2), odt), std::invalid_argument);  // 3D, not linear tetra

    Mesh mixed = tet_block(2);
    mixed.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));  // non-3D alongside its tets
    EXPECT_THROW(smooth(mixed, odt), std::invalid_argument);

    Mesh empty;
    EXPECT_THROW(smooth(empty, odt), std::invalid_argument);  // no tetra cell block at all
}

}  // namespace
