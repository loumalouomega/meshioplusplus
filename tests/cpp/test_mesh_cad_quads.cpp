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
// mesh_cad_surface's quadrilaterals (roadmap 6.8): quad-dominant meshes by
// greedy recombination, all-quad ones by recombining a mesh of twice the size
// and splitting it through its edges' midpoints.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <cstdio>
#include <map>
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
#include "../../src/cpp/src/detail/curved_quad.hpp"
#include "../../src/cpp/src/detail/quad_recombine.hpp"
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/quality.hpp"

namespace {

namespace mio = meshioplusplus;

constexpr double kPi = 3.14159265358979323846;

std::uint64_t digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

mio::MeshCadSurfaceOptions sized(double H, double MinAngle = 0.0) {
    mio::MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    o.mMinAngleDeg = MinAngle;
    return o;
}

struct Quads {
    mio::MeshCadSurfaceResult mResult;
    mio::MeshCadSurfaceExtrasReport mReport;
};

Quads run(const mio::BrepInfo& rB, const mio::MeshCadSurfaceOptions& rO, mio::MeshCadQuads Q) {
    mio::MeshCadSurfaceExtras ex;
    ex.mQuads = Q;
    Quads out;
    out.mResult =
        mio::mesh_cad_surface(rB, rO, mio::MeshCadSizeSources{}, ex, nullptr, &out.mReport);
    return out;
}

/// The block of type `rType`, or -1.
std::int64_t block(const mio::Mesh& rM, const std::string& rType) {
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b)
        if (rM.Cells(b).Type() == rType)
            return static_cast<std::int64_t>(b);
    return -1;
}

/// Every edge of the triangle and quad blocks used twice, once each way.
bool watertight(const mio::Mesh& rM) {
    std::map<std::pair<std::int64_t, std::int64_t>, int> uses;
    for (const std::string type : {"triangle", "quad"}) {
        const std::int64_t b = block(rM, type);
        if (b < 0)
            continue;
        const auto cells = rM.Cells(static_cast<std::size_t>(b));
        const std::size_t n = type == "triangle" ? 3 : 4;
        for (std::size_t c = 0; c < cells.NumCells(); ++c)
            for (std::size_t k = 0; k < n; ++k)
                ++uses[{mio::detail::read_int(cells.Conn(), n * c + k),
                        mio::detail::read_int(cells.Conn(), n * c + (k + 1) % n)}];
    }
    for (const auto& [e, n] : uses)
        if (n != 1 || uses.count({e.second, e.first}) == 0)
            return false;
    return true;
}

/// The surface area of the triangle and quad blocks (each quad by its two
/// triangles from corner 0).
double area(const mio::Mesh& rM) {
    const double* x = rM.Points().As<double>();
    const auto tri = [&](std::int64_t A, std::int64_t B, std::int64_t C) {
        const double* a = x + 3 * A;
        const double* b = x + 3 * B;
        const double* c = x + 3 * C;
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        return 0.5 * std::hypot(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                u[0] * v[1] - u[1] * v[0]);
    };
    double s = 0.0;
    for (const std::string type : {"triangle", "quad"}) {
        const std::int64_t b = block(rM, type);
        if (b < 0)
            continue;
        const auto cells = rM.Cells(static_cast<std::size_t>(b));
        const std::size_t n = type == "triangle" ? 3 : 4;
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            std::array<std::int64_t, 4> p{};
            for (std::size_t k = 0; k < n; ++k)
                p[k] = mio::detail::read_int(cells.Conn(), n * c + k);
            s += tri(p[0], p[1], p[2]);
            if (n == 4)
                s += tri(p[0], p[2], p[3]);
        }
    }
    return s;
}

/// compute_quality's least scaled Jacobian over the mesh.
double min_scaled_jacobian(const mio::Mesh& rM) {
    const mio::QualityReport q = mio::compute_quality(rM);
    for (const auto& [name, r_s] : q.mMetrics)
        if (name == "quality:scaled_jacobian")
            return r_s.mMin;
    return std::numeric_limits<double>::quiet_NaN();
}

// Roadmap 6.8.3's probes: the box and the capped cylinder all-quad. No
// triangle is left, every quad is convex (its scaled Jacobian, the report's
// and compute_quality's alike, well above zero), the surface is closed, every
// point lies on the solid, each quad keeps its CAD face, the edges are linear
// lines, and the mesh does not depend on the threads.
TEST(MeshCadQuads, TheBoxAndTheCylinderAreAllQuad) {
    struct S {
        const char* mName;
        mio::BrepInfo mB;
        double mArea;
        double mAreaTol;
        std::function<bool(const double*)> mOn;
    };
    const std::vector<S> shapes = {
        {"box", brep_solids::box(), 6.0, 1e-12,
         [](const double* pX) {
             bool on = false;
             for (int k = 0; k < 3; ++k)
                 on = on || std::fabs(pX[k]) < 1e-12 || std::fabs(pX[k] - 1.0) < 1e-12;
             return on;
         }},
        {"cylinder", brep_solids::capped_cylinder(1.0, 2.0), 6.0 * kPi, 5e-3, [](const double* pX) {
             return std::fabs(std::hypot(pX[0], pX[1]) - 1.0) < 1e-9 || std::fabs(pX[2]) < 1e-12 ||
                    std::fabs(pX[2] - 2.0) < 1e-12;
         }}};
    for (const auto& r_s : shapes) {
        const auto r = run(r_s.mB, sized(0.2, 25.0), mio::MeshCadQuads::All);
        const mio::Mesh& m = r.mResult.mMesh;
        ASSERT_EQ(r.mResult.mNumFailed, 0) << r_s.mName;
        EXPECT_EQ(block(m, "triangle"), -1) << r_s.mName;
        ASSERT_GE(block(m, "quad"), 0) << r_s.mName;
        EXPECT_GE(block(m, "line"), 0) << r_s.mName;
        EXPECT_EQ(block(m, "line3"), -1) << r_s.mName;
        EXPECT_EQ(r.mReport.mNumLeftoverTriangles, 0);
        EXPECT_EQ(r.mResult.mNumTriangles, 0);
        EXPECT_EQ(static_cast<std::size_t>(r.mReport.mNumQuads),
                  m.Cells(static_cast<std::size_t>(block(m, "quad"))).NumCells());
        EXPECT_EQ(r.mResult.mNumOpenEdges, 0) << r_s.mName;
        EXPECT_EQ(r.mResult.mNumNonManifoldEdges, 0) << r_s.mName;
        EXPECT_EQ(r.mResult.mNumMisorientedEdges, 0) << r_s.mName;
        EXPECT_TRUE(watertight(m)) << r_s.mName;
        const double sj = min_scaled_jacobian(m);
        std::printf(
            "[6.8] %s all-quad at h 0.2: %lld quads, least scaled Jacobian %.3f, area "
            "error %.2e\n",
            r_s.mName, static_cast<long long>(r.mReport.mNumQuads), sj,
            std::fabs(area(m) - r_s.mArea) / r_s.mArea);
        EXPECT_NEAR(sj, r.mReport.mMinQuadScaledJacobian, 1e-12) << r_s.mName;
        EXPECT_GT(sj, 0.3) << r_s.mName;
        EXPECT_LT(std::fabs(area(m) - r_s.mArea) / r_s.mArea, r_s.mAreaTol) << r_s.mName;
        const double* x = m.Points().As<double>();
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            EXPECT_TRUE(r_s.mOn(x + 3 * p)) << r_s.mName << " " << p;
        // Each quad keeps its face; every face has quads.
        const auto qb = static_cast<std::size_t>(block(m, "quad"));
        const auto& face = m.CellData(mio::kMeshCadFaceName, qb);
        std::vector<std::int64_t> per(r_s.mB.mFaces.size(), 0);
        for (std::size_t c = 0; c < m.Cells(qb).NumCells(); ++c)
            ++per[static_cast<std::size_t>(mio::detail::read_int(face, c))];
        for (const std::int64_t n : per)
            EXPECT_GT(n, 0) << r_s.mName;
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
        const std::uint64_t ref = digest(m);
        const int saved = omp_get_max_threads();
        for (int threads : {1, 2, 4}) {
            omp_set_num_threads(threads);
            EXPECT_EQ(digest(run(r_s.mB, sized(0.2, 25.0), mio::MeshCadQuads::All).mResult.mMesh),
                      ref)
                << threads;
        }
        omp_set_num_threads(saved);
#endif
    }
}

// Quad-dominant: the triangles that pair into good quads are recombined, the
// rest stay beside them; the surface stays closed, the quads convex, and
// with quads off the mesh is the plain one byte for byte.
TEST(MeshCadQuads, TheSphereIsQuadDominant) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    const auto plain = mio::mesh_cad_surface(sphere, sized(0.1, 25.0));
    const auto r = run(sphere, sized(0.1, 25.0), mio::MeshCadQuads::Dominant);
    const mio::Mesh& m = r.mResult.mMesh;
    ASSERT_EQ(r.mResult.mNumFailed, 0);
    ASSERT_GE(block(m, "quad"), 0);
    ASSERT_GE(block(m, "triangle"), 0);
    EXPECT_EQ(static_cast<std::size_t>(r.mReport.mNumLeftoverTriangles),
              m.Cells(static_cast<std::size_t>(block(m, "triangle"))).NumCells());
    EXPECT_EQ(2 * r.mReport.mNumQuads + r.mReport.mNumLeftoverTriangles, plain.mNumTriangles);
    EXPECT_EQ(r.mResult.mNumPoints, plain.mNumPoints);
    EXPECT_TRUE(watertight(m));
    EXPECT_EQ(r.mResult.mNumOpenEdges, 0);
    const double share = double(2 * r.mReport.mNumQuads) / double(plain.mNumTriangles);
    std::printf(
        "[6.8] sphere quad-dominant at h 0.1: %lld quads and %lld triangles (%.0f%% of "
        "the triangles recombined), least scaled Jacobian %.3f\n",
        static_cast<long long>(r.mReport.mNumQuads),
        static_cast<long long>(r.mReport.mNumLeftoverTriangles), 100.0 * share,
        r.mReport.mMinQuadScaledJacobian);
    EXPECT_GT(share, 0.75);
    EXPECT_GT(r.mReport.mMinQuadScaledJacobian, 0.3);
    EXPECT_NEAR(area(m), area(plain.mMesh), 1e-9);
    // Off: byte for byte.
    const auto none = run(sphere, sized(0.1, 25.0), mio::MeshCadQuads::None);
    EXPECT_EQ(digest(none.mResult.mMesh), digest(plain.mMesh));
    EXPECT_EQ(none.mReport.mNumQuads, 0);
    // A sized run recombines its second meshing only.
    mio::MeshCadSurfaceOptions o = sized(0.2, 25.0);
    o.mCurvatureDeviation = 0.01;
    const auto c = run(sphere, o, mio::MeshCadQuads::Dominant);
    EXPECT_EQ(c.mResult.mNumFailed, 0);
    EXPECT_GT(c.mReport.mNumQuads, 0);
    EXPECT_TRUE(watertight(c.mResult.mMesh));
}

// Roadmap 6.8.1: pairing by a maximum matching (Blossom-quad, the default)
// leaves under half the triangles greedy pairing does on the sphere at
// h = 0.1 (224 against 542), and on a face of more triangles than one patch
// (h = 0.04, the patches and the band) 30% fewer at least (1,382 against
// 2,426); the result does not depend on threads.
TEST(MeshCadQuads, AMatchingLeavesFewerTriangles) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    for (const double h : {0.1, 0.04}) {
        mio::MeshCadSurfaceExtras ex;
        ex.mQuads = mio::MeshCadQuads::Dominant;
        mio::MeshCadSurfaceExtrasReport matched, greedy;
        const auto rm = mio::mesh_cad_surface(sphere, sized(h, 25.0), mio::MeshCadSizeSources{}, ex,
                                              nullptr, &matched);
        ex.mQuadPairs = mio::MeshCadQuadPairs::Greedy;
        const auto rg = mio::mesh_cad_surface(sphere, sized(h, 25.0), mio::MeshCadSizeSources{}, ex,
                                              nullptr, &greedy);
        std::printf(
            "[6.8.1] sphere at h %.2f: %lld triangles left by the matching, %lld by "
            "greedy pairing; least scaled Jacobian %.3f against %.3f\n",
            h, static_cast<long long>(matched.mNumLeftoverTriangles),
            static_cast<long long>(greedy.mNumLeftoverTriangles), matched.mMinQuadScaledJacobian,
            greedy.mMinQuadScaledJacobian);
        EXPECT_EQ(2 * matched.mNumQuads + matched.mNumLeftoverTriangles,
                  2 * greedy.mNumQuads + greedy.mNumLeftoverTriangles);
        // Exact on one patch: under half; in patches: at least 30% fewer.
        if (h > 0.05)
            EXPECT_LT(2 * matched.mNumLeftoverTriangles, greedy.mNumLeftoverTriangles) << h;
        else
            EXPECT_LT(10 * matched.mNumLeftoverTriangles, 7 * greedy.mNumLeftoverTriangles) << h;
        EXPECT_GT(matched.mMinQuadScaledJacobian, 0.3);
        EXPECT_TRUE(watertight(rm.mMesh));
        EXPECT_TRUE(watertight(rg.mMesh));
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
        ex.mQuadPairs = mio::MeshCadQuadPairs::Matching;
        const std::uint64_t ref = digest(rm.mMesh);
        const int saved = omp_get_max_threads();
        for (int threads : {1, 4}) {
            omp_set_num_threads(threads);
            EXPECT_EQ(
                digest(mio::mesh_cad_surface(sphere, sized(h, 25.0), mio::MeshCadSizeSources{}, ex)
                           .mMesh),
                ref)
                << threads;
        }
        omp_set_num_threads(saved);
#endif
    }
}

TEST(MeshCadQuads, RefusalsAreByName) {
    const mio::BrepInfo box = brep_solids::box();
    mio::MeshCadSurfaceExtras ex;
    ex.mQuads = mio::MeshCadQuads::All;
    mio::MeshCadSurfaceOptions o = sized(0.3);
    // All-quad meshes at order 2, whose curved elements cannot cross a
    // suppressed corner.
    mio::MeshCadSizeSources ve;
    ve.mVirtualEdges = true;
    EXPECT_THROW((void)mio::mesh_cad_surface(box, o, ve, ex), std::invalid_argument);
    ex.mQuads = static_cast<mio::MeshCadQuads>(7);
    EXPECT_THROW((void)mio::mesh_cad_surface(box, o, mio::MeshCadSizeSources{}, ex),
                 std::invalid_argument);
    ex.mQuads = mio::MeshCadQuads::Dominant;
    ex.mQuadPairs = static_cast<mio::MeshCadQuadPairs>(5);
    EXPECT_THROW((void)mio::mesh_cad_surface(box, o, mio::MeshCadSizeSources{}, ex),
                 std::invalid_argument);
    // What 6.8.3 opened: every order, untangling and virtual faces.
    ex.mQuadPairs = mio::MeshCadQuadPairs::Matching;
    for (const auto q : {mio::MeshCadQuads::Dominant, mio::MeshCadQuads::All})
        for (const int order : {1, 2, 3}) {
            ex.mQuads = q;
            o.mOrder = order;
            ex.mUntangle = order > 1;
            ex.mVirtualFaces = true;
            EXPECT_NO_THROW((void)mio::mesh_cad_surface(box, o, mio::MeshCadSizeSources{}, ex))
                << order;
        }
}

// Roadmap 6.8.3: quads="all" asks a size field for twice its size. The
// capped cylinder sized by its curvature is all-quad, closed, finer than
// the uniform mesh, and deterministic.
TEST(MeshCadQuads, AllQuadTakesSizeSources) {
    const mio::BrepInfo cyl = brep_solids::capped_cylinder(1.0, 2.0);
    mio::MeshCadSurfaceOptions o = sized(0.2);
    o.mCurvatureDeviation = 0.002;
    const auto u = run(cyl, sized(0.2), mio::MeshCadQuads::All);
    const auto r = run(cyl, o, mio::MeshCadQuads::All);
    std::printf("[6.8.3] the cylinder all-quad: %lld quads uniform, %lld sized\n",
                static_cast<long long>(u.mReport.mNumQuads),
                static_cast<long long>(r.mReport.mNumQuads));
    EXPECT_EQ(r.mResult.mNumFailed, 0);
    EXPECT_EQ(r.mReport.mNumLeftoverTriangles, 0);
    EXPECT_EQ(block(r.mResult.mMesh, "triangle"), -1);
    EXPECT_GT(r.mReport.mNumQuads, 2 * u.mReport.mNumQuads);
    EXPECT_TRUE(watertight(r.mResult.mMesh));
    EXPECT_GT(r.mReport.mMinQuadScaledJacobian, 0.2);
    EXPECT_EQ(digest(run(cyl, o, mio::MeshCadQuads::All).mResult.mMesh), digest(r.mResult.mMesh));
}

/// The area of the quad9 block, by 5 x 5 Gauss points on each quad.
double quad9_area(const mio::Mesh& rM) {
    const std::int64_t b = block(rM, "quad9");
    if (b < 0)
        return 0.0;
    const double* x = rM.Points().As<double>();
    const auto cells = rM.Cells(static_cast<std::size_t>(b));
    const double g[5] = {0.0469100770306680, 0.2307653449471585, 0.5, 0.7692346550528415,
                         0.9530899229693320};
    const double w[5] = {0.1184634425280945, 0.2393143352496832, 0.2844444444444444,
                         0.2393143352496832, 0.1184634425280945};
    double s = 0.0;
    std::array<mio::detail::CurvedVec3, 9> n{};
    for (std::size_t c = 0; c < cells.NumCells(); ++c) {
        for (std::size_t k = 0; k < 9; ++k) {
            const auto id =
                static_cast<std::size_t>(mio::detail::read_int(cells.Conn(), 9 * c + k));
            n[k] = {x[3 * id], x[3 * id + 1], x[3 * id + 2]};
        }
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j) {
                // |x_xi x x_eta| by the Jacobian against its own direction.
                const double e = 1e-7;
                const auto p0 = mio::detail::curved_quad_point(n.data(), g[i], g[j]);
                const auto pu = mio::detail::curved_quad_point(n.data(), g[i] + e, g[j]);
                const auto pv = mio::detail::curved_quad_point(n.data(), g[i], g[j] + e);
                const double a[3] = {(pu[0] - p0[0]) / e, (pu[1] - p0[1]) / e, (pu[2] - p0[2]) / e};
                const double v[3] = {(pv[0] - p0[0]) / e, (pv[1] - p0[1]) / e, (pv[2] - p0[2]) / e};
                const mio::detail::CurvedVec3 cr{a[1] * v[2] - a[2] * v[1],
                                                 a[2] * v[0] - a[0] * v[2],
                                                 a[0] * v[1] - a[1] * v[0]};
                const double len = std::sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
                s += w[i] * w[j] *
                     (len > 0.0 ? mio::detail::curved_quad_jacobian(
                                      n.data(), {cr[0] / len, cr[1] / len, cr[2] / len}, g[i], g[j])
                                : 0.0);
            }
    }
    return s;
}

// Roadmap 6.8.3: quads="all" with order 2 writes quad9. On the sphere every
// node lies on it, every quad is proved valid, the surface is closed, the
// edges are line3, and the area converges at order 4 where the linear
// quads' does at order 2.
TEST(MeshCadQuads, CurvedQuadsLieOnTheSphereAndConverge) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    std::array<double, 2> err{}, err_lin{};
    for (std::size_t i = 0; i < 2; ++i) {
        const double h = i == 0 ? 0.4 : 0.2;
        mio::MeshCadSurfaceOptions o = sized(h);
        o.mOrder = 2;
        const auto r = run(sphere, o, mio::MeshCadQuads::All);
        const mio::Mesh& m = r.mResult.mMesh;
        ASSERT_EQ(r.mResult.mNumFailed, 0);
        EXPECT_EQ(r.mResult.mNumCurvedStraightened, 0);
        EXPECT_EQ(r.mResult.mNumOpenEdges, 0);
        const std::int64_t b = block(m, "quad9");
        ASSERT_GE(b, 0);
        EXPECT_EQ(block(m, "quad"), -1);
        EXPECT_EQ(block(m, "triangle"), -1);
        const double* x = m.Points().As<double>();
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            EXPECT_NEAR(std::hypot(x[3 * p], x[3 * p + 1], x[3 * p + 2]), 1.0, 1e-12) << p;
        const auto cells = m.Cells(static_cast<std::size_t>(b));
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            std::array<mio::detail::CurvedVec3, 9> n{};
            for (std::size_t k = 0; k < 9; ++k) {
                const auto id =
                    static_cast<std::size_t>(mio::detail::read_int(cells.Conn(), 9 * c + k));
                n[k] = {x[3 * id], x[3 * id + 1], x[3 * id + 2]};
            }
            // The outward normal at the centre.
            EXPECT_TRUE(mio::detail::curved_quad_check(n.data(), n[8]).mValid) << c;
        }
        err[i] = std::fabs(quad9_area(m) - 4.0 * kPi);
        err_lin[i] = std::fabs(area(run(sphere, sized(h), mio::MeshCadQuads::All).mResult.mMesh) -
                               4.0 * kPi);
        if (i == 1) {
            EXPECT_EQ(digest(run(sphere, o, mio::MeshCadQuads::All).mResult.mMesh), digest(m));
        }
    }
    const double rate = std::log2(err[0] / err[1]);
    const double rate_lin = std::log2(err_lin[0] / err_lin[1]);
    std::printf(
        "[6.8.3] the sphere's area: quad9 %.2e -> %.2e (order %.2f), quad %.2e -> %.2e "
        "(order %.2f)\n",
        err[0], err[1], rate, err_lin[0], err_lin[1], rate_lin);
    EXPECT_GT(rate, 3.5);
    EXPECT_GT(rate_lin, 1.8);
    EXPECT_LT(err[1], 1e-4);
}

// The quad9 validity check alone: the unit square is valid; its centre
// pulled past a corner folds it at the nodes; an edge node pushed across
// the opposite one folds it between them though every node's Jacobian
// stays positive.
TEST(MeshCadQuads, TheQuad9CheckFindsFolds) {
    using V = mio::detail::CurvedVec3;
    std::array<V, 9> sq{V{0, 0, 0},   V{1, 0, 0},   V{1, 1, 0},   V{0, 1, 0},    V{0.5, 0, 0},
                        V{1, 0.5, 0}, V{0.5, 1, 0}, V{0, 0.5, 0}, V{0.5, 0.5, 0}};
    const V up{0, 0, 1};
    EXPECT_TRUE(mio::detail::curved_quad_check(sq.data(), up).mValid);
    EXPECT_NEAR(mio::detail::curved_quad_jacobian(sq.data(), up, 0.3, 0.7), 1.0, 1e-14);
    auto bad = sq;
    bad[8] = {1.4, 1.4, 0};
    EXPECT_FALSE(mio::detail::curved_quad_check(bad.data(), up).mValid);
    // A bent edge: the bottom's midpoint raised to 0.6 of the height keeps
    // the nodal Jacobians positive but folds the element near it.
    bad = sq;
    bad[4] = {0.5, 0.6, 0};
    const auto c = mio::detail::curved_quad_check(bad.data(), up);
    EXPECT_FALSE(c.mValid);
    EXPECT_TRUE(c.mDecided);
    // Mildly bent, it stays valid.
    bad[4] = {0.5, 0.1, 0};
    EXPECT_TRUE(mio::detail::curved_quad_check(bad.data(), up).mValid);
    // Against the opposite normal every element is invalid.
    EXPECT_FALSE(mio::detail::curved_quad_check(sq.data(), V{0, 0, -1}).mValid);
}

/// The area of a quad9 or quad16 block (`Order` 2 or 3), by 6 x 6 Gauss
/// points on each quad, from its tangents' cross product.
double curved_quad_area(const mio::Mesh& rM, int Order) {
    const std::int64_t b = block(rM, Order == 2 ? "quad9" : "quad16");
    if (b < 0)
        return 0.0;
    const std::size_t nn = mio::detail::curved_quad_num_nodes(Order);
    const double* x = rM.Points().As<double>();
    const auto cells = rM.Cells(static_cast<std::size_t>(b));
    const double g[6] = {0.0337652428984240, 0.1693953067668677, 0.3806904069584015,
                         0.6193095930415985, 0.8306046932331323, 0.9662347571015760};
    const double w[6] = {0.0856622461895852, 0.1803807865240693, 0.2339569672863455,
                         0.2339569672863455, 0.1803807865240693, 0.0856622461895852};
    double s = 0.0;
    std::array<mio::detail::CurvedVec3, 16> n{};
    for (std::size_t c = 0; c < cells.NumCells(); ++c) {
        for (std::size_t k = 0; k < nn; ++k) {
            const auto id =
                static_cast<std::size_t>(mio::detail::read_int(cells.Conn(), nn * c + k));
            n[k] = {x[3 * id], x[3 * id + 1], x[3 * id + 2]};
        }
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                const double e = 1e-7;
                const auto p0 = mio::detail::curved_quad_point(Order, n.data(), g[i], g[j]);
                const auto pu = mio::detail::curved_quad_point(Order, n.data(), g[i] + e, g[j]);
                const auto pv = mio::detail::curved_quad_point(Order, n.data(), g[i], g[j] + e);
                const double a[3] = {(pu[0] - p0[0]) / e, (pu[1] - p0[1]) / e, (pu[2] - p0[2]) / e};
                const double v[3] = {(pv[0] - p0[0]) / e, (pv[1] - p0[1]) / e, (pv[2] - p0[2]) / e};
                const mio::detail::CurvedVec3 cr{a[1] * v[2] - a[2] * v[1],
                                                 a[2] * v[0] - a[0] * v[2],
                                                 a[0] * v[1] - a[1] * v[0]};
                const double len = std::sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
                s += w[i] * w[j] *
                     (len > 0.0 ? mio::detail::curved_quad_jacobian(
                                      Order, n.data(), {cr[0] / len, cr[1] / len, cr[2] / len},
                                      g[i], g[j])
                                : 0.0);
            }
    }
    return s;
}

/// Points no cell uses.
std::size_t unused_points(const mio::Mesh& rM) {
    std::vector<char> used(rM.NumPoints(), 0);
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
        const auto cells = rM.Cells(b);
        const std::size_t nn = cells.NodesPerCell();
        for (std::size_t i = 0; i < cells.NumCells() * nn; ++i) {
            const std::int64_t v = mio::detail::read_int(cells.Conn(), i);
            if (v < 0 || static_cast<std::size_t>(v) >= used.size()) {
                ADD_FAILURE() << "block " << b << " names point " << v;
                continue;
            }
            used[static_cast<std::size_t>(v)] = 1;
        }
    }
    return static_cast<std::size_t>(std::count(used.begin(), used.end(), 0));
}

// Roadmap 6.8.3: quads="all" with order 3 writes quad16 whose edge and inside
// nodes sit at the thirds of the parameters of one chart (a cubic's accuracy
// rests on one lattice under the surface's map, 6.6.2), placed by projection
// where that chart leaves its patch or folds an element. On the sphere every
// node lies on it, every quad16 is proved valid, the edges are line4, no
// point is unused, and the area beats quad9's at each size, converging at
// order 4 or better.
TEST(MeshCadQuads, CubicQuadsLieOnTheSphereAndBeatQuad9) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    std::array<double, 2> err{}, err9{};
    for (std::size_t i = 0; i < 2; ++i) {
        const double h = i == 0 ? 0.4 : 0.2;
        mio::MeshCadSurfaceOptions o = sized(h);
        o.mOrder = 3;
        const auto r = run(sphere, o, mio::MeshCadQuads::All);
        const mio::Mesh& m = r.mResult.mMesh;
        ASSERT_EQ(r.mResult.mNumFailed, 0);
        EXPECT_EQ(r.mResult.mNumCurvedStraightened, 0);
        EXPECT_EQ(r.mResult.mNumOpenEdges, 0);
        EXPECT_EQ(unused_points(m), 0u);
        const std::int64_t b = block(m, "quad16");
        ASSERT_GE(b, 0);
        EXPECT_GE(block(m, "line4"), 0);
        const double* x = m.Points().As<double>();
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            EXPECT_NEAR(std::hypot(x[3 * p], x[3 * p + 1], x[3 * p + 2]), 1.0, 1e-12) << p;
        const auto cells = m.Cells(static_cast<std::size_t>(b));
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            std::array<mio::detail::CurvedVec3, 16> n{};
            mio::detail::CurvedVec3 mid{0.0, 0.0, 0.0};
            for (std::size_t k = 0; k < 16; ++k) {
                const auto id =
                    static_cast<std::size_t>(mio::detail::read_int(cells.Conn(), 16 * c + k));
                n[k] = {x[3 * id], x[3 * id + 1], x[3 * id + 2]};
                for (std::size_t d = 0; d < 3; ++d)
                    mid[d] += n[k][d] / 16.0;
            }
            EXPECT_TRUE(mio::detail::curved_quad_check(3, n.data(), mid).mValid) << c;
        }
        err[i] = std::fabs(curved_quad_area(m, 3) - 4.0 * kPi);
        o.mOrder = 2;
        err9[i] = std::fabs(
            curved_quad_area(run(sphere, o, mio::MeshCadQuads::All).mResult.mMesh, 2) - 4.0 * kPi);
        EXPECT_LT(err[i], err9[i]) << h;
    }
    const double rate = std::log2(err[0] / err[1]);
    std::printf("[6.8.3] the sphere's area: quad16 %.2e -> %.2e (order %.2f), quad9 %.2e -> %.2e\n",
                err[0], err[1], rate, err9[0], err9[1]);
    EXPECT_GT(rate, 3.5);
}

// The quad16 check alone: the unit square on its lattice is valid with a
// Jacobian of 1; an inside node pulled past a corner folds it; against the
// opposite normal it is invalid.
TEST(MeshCadQuads, TheQuad16CheckFindsFolds) {
    using V = mio::detail::CurvedVec3;
    std::array<V, 16> sq{};
    for (std::size_t k = 0; k < 16; ++k) {
        const auto p = mio::detail::curved_quad_node_position(3, k);
        sq[k] = {p[0], p[1], 0.0};
    }
    const V up{0, 0, 1};
    EXPECT_TRUE(mio::detail::curved_quad_check(3, sq.data(), up).mValid);
    EXPECT_NEAR(mio::detail::curved_quad_jacobian(3, sq.data(), up, 0.3, 0.7), 1.0, 1e-13);
    EXPECT_NEAR(mio::detail::curved_quad_scaled_bound(3, sq.data(), up), 1.0, 1e-12);
    auto bad = sq;
    bad[12] = {1.3, 1.3, 0};
    EXPECT_FALSE(mio::detail::curved_quad_check(3, bad.data(), up).mValid);
    EXPECT_LT(mio::detail::curved_quad_scaled_bound(3, bad.data(), up), 0.0);
    EXPECT_FALSE(mio::detail::curved_quad_check(3, sq.data(), V{0, 0, -1}).mValid);
    // The order-2 entry points are quad9's.
    std::array<V, 9> sq9{};
    for (std::size_t k = 0; k < 9; ++k) {
        const auto p = mio::detail::curved_quad_node_position(2, k);
        sq9[k] = {p[0], p[1], 0.0};
    }
    EXPECT_TRUE(mio::detail::curved_quad_check(sq9.data(), up).mValid);
    EXPECT_EQ(mio::detail::curved_quad_num_nodes(2), 9u);
    EXPECT_EQ(mio::detail::curved_quad_num_nodes(3), 16u);
}

// Roadmap 6.8.3: quad-dominant meshes at order 2 and 3 keep their triangles'
// curved edges: triangle6 beside quad9, triangle10 beside quad16, every node
// on the sphere, every element proved valid or straightened, the surface
// closed, and no point unused (a pair's diagonal nodes are dropped). The
// all-quad split drops the pairs' diagonal midpoints too, at every order.
TEST(MeshCadQuads, CurvedQuadDominantAndNoUnusedPoints) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    for (const int order : {2, 3}) {
        mio::MeshCadSurfaceOptions o = sized(0.3);
        o.mOrder = order;
        const auto r = run(sphere, o, mio::MeshCadQuads::Dominant);
        const mio::Mesh& m = r.mResult.mMesh;
        ASSERT_EQ(r.mResult.mNumFailed, 0);
        EXPECT_EQ(r.mResult.mNumOpenEdges, 0);
        EXPECT_GE(block(m, order == 2 ? "quad9" : "quad16"), 0);
        EXPECT_GE(block(m, order == 2 ? "triangle6" : "triangle10"), 0);
        EXPECT_EQ(unused_points(m), 0u) << order;
        EXPECT_GT(r.mReport.mNumQuads, 0);
        EXPECT_EQ(digest(run(sphere, o, mio::MeshCadQuads::Dominant).mResult.mMesh), digest(m));
    }
    for (const int order : {1, 2, 3}) {
        mio::MeshCadSurfaceOptions o = sized(0.3);
        o.mOrder = order;
        EXPECT_EQ(unused_points(run(sphere, o, mio::MeshCadQuads::All).mResult.mMesh), 0u) << order;
    }
}

// Roadmap 6.6.3 with quads: untangling moves the quads' and the triangles'
// own nodes inside each face. On the sphere's quad-dominant mesh at h = 0.3
// it rescues every element straightening alone would.
TEST(MeshCadQuads, UntanglingTakesQuads) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    for (const int order : {2, 3}) {
        mio::MeshCadSurfaceOptions o = sized(0.3);
        o.mOrder = order;
        const auto plain = run(sphere, o, mio::MeshCadQuads::Dominant);
        mio::MeshCadSurfaceExtras ex;
        ex.mQuads = mio::MeshCadQuads::Dominant;
        ex.mUntangle = true;
        mio::MeshCadSurfaceExtrasReport rep;
        const auto r =
            mio::mesh_cad_surface(sphere, o, mio::MeshCadSizeSources{}, ex, nullptr, &rep);
        std::printf(
            "[6.6.3] sphere quad-dominant order %d: %lld straightened plain, %lld untangled\n",
            order, static_cast<long long>(plain.mResult.mNumCurvedStraightened),
            static_cast<long long>(rep.mNumUntangled));
        EXPECT_LE(r.mNumCurvedStraightened, plain.mResult.mNumCurvedStraightened);
        EXPECT_EQ(rep.mNumUntangled,
                  plain.mResult.mNumCurvedStraightened - r.mNumCurvedStraightened);
    }
}

// The kernel's side curves (6.8.3): a pair whose quad would have both sides
// at a corner on one curve is no candidate, since a curved element is a
// straight angle there.
TEST(MeshCadQuads, NoCurvedQuadHasTwoSidesOnOneCurve) {
    const std::vector<double> x = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    const std::vector<std::array<std::int64_t, 3>> square = {{0, 1, 2}, {0, 2, 3}};
    // (0, 1, 2) and (0, 2, 3) pair into the quad (0, 1, 2, 3), whose sides
    // 1-2 (the first's side 1) and 2-3 (the second's side 1) meet at corner 2.
    std::vector<std::array<std::int64_t, 3>> curves = {{-1, 7, -1}, {-1, 7, -1}};
    EXPECT_TRUE(mio::detail::quad_recombine(square, {0, 0}, x.data(), 0.3,
                                            mio::detail::QuadPairing::Greedy, &curves)
                    .mQuads.empty());
    curves = {{-1, 7, -1}, {-1, 8, -1}};
    EXPECT_EQ(mio::detail::quad_recombine(square, {0, 0}, x.data(), 0.3,
                                          mio::detail::QuadPairing::Greedy, &curves)
                  .mQuads.size(),
              1u);
}

// The recombination kernel alone: two right triangles of a unit square make
// one quad whose quality is 1; a pair folded across a crease of 90 degrees
// scores 0 and stays apart.
TEST(MeshCadQuads, TheKernelMergesTheBestPairs) {
    const std::vector<double> x = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1};
    const std::vector<std::array<std::int64_t, 3>> square = {{0, 1, 2}, {0, 2, 3}};
    const auto r = mio::detail::quad_recombine(square, {0, 0}, x.data(), 0.3);
    ASSERT_EQ(r.mQuads.size(), 1u);
    EXPECT_TRUE(r.mLeft.empty());
    EXPECT_NEAR(mio::detail::quad_recombine_quality(x.data(), r.mQuads[0]), 1.0, 1e-12);
    EXPECT_NEAR(mio::detail::quad_scaled_jacobian(x.data(), r.mQuads[0]), 1.0, 1e-12);
    // Different groups never merge.
    EXPECT_TRUE(mio::detail::quad_recombine(square, {0, 1}, x.data(), 0.3).mQuads.empty());
    // Folded: (0, 2, 3) against (0, 3, 4) up the y = 1 wall.
    const std::vector<std::array<std::int64_t, 3>> fold = {{0, 2, 3}, {3, 2, 4}};
    EXPECT_TRUE(mio::detail::quad_recombine(fold, {0, 0}, x.data(), 0.3).mQuads.empty());
}

}  // namespace
