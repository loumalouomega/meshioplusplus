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
// Surface reconstruction from an oriented point cloud (roadmap 7.5.2).

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "thread_sweep.hpp"

#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/reconstruct_surface.hpp"

using namespace meshioplusplus;

namespace {

/// The Fibonacci sphere of @p N points and radius @p R, with radial normals
/// (or none).
Mesh rs_sphere(std::int64_t N, double R, bool Normals = true) {
    const Mesh pts = sphere_points(N, R);
    Mesh m;
    NDArray p(DType::Float64, {pts.NumPoints(), std::size_t{3}});
    NDArray n(DType::Float64, {pts.NumPoints(), std::size_t{3}});
    const double* src = pts.Points().As<double>();
    for (std::size_t i = 0; i < pts.NumPoints() * 3; ++i) {
        p.As<double>()[i] = src[i];
        n.As<double>()[i] = src[i] / R;
    }
    m.AssignPoints(std::move(p));
    if (Normals)
        m.AddPointData("normals", std::move(n));
    return m;
}

/// The edges used by exactly one triangle.
std::size_t rs_open_edges(const Mesh& rMesh) {
    std::map<std::pair<std::int64_t, std::int64_t>, int> uses;
    const auto* c = rMesh.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < rMesh.Cells(0).NumCells(); ++t)
        for (std::size_t e = 0; e < 3; ++e) {
            const std::int64_t a = c[t * 3 + e], b = c[t * 3 + (e + 1) % 3];
            ++uses[{std::min(a, b), std::max(a, b)}];
        }
    std::size_t open = 0;
    for (const auto& [edge, count] : uses)
        open += count == 1 ? 1 : 0;
    return open;
}

/// Six times the signed volume the triangles enclose.
double rs_volume6(const Mesh& rMesh) {
    const auto* c = rMesh.Cells(0).Conn().As<std::int64_t>();
    const double* p = rMesh.Points().As<double>();
    double v = 0.0;
    for (std::size_t t = 0; t < rMesh.Cells(0).NumCells(); ++t) {
        const double* a = p + 3 * c[t * 3];
        const double* b = p + 3 * c[t * 3 + 1];
        const double* d = p + 3 * c[t * 3 + 2];
        v += a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
             a[2] * (b[0] * d[1] - b[1] * d[0]);
    }
    return v;
}

}  // namespace

TEST(ReconstructSurface, SphereWithinOneSpacing) {
    // The roadmap's probe: a cloud sampled on a sphere reconstructs within one
    // sample spacing in Hausdorff distance, closed and outward.
    const ReconstructSurfaceResult r = reconstruct_surface(rs_sphere(4000, 1.0));
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "triangle");
    EXPECT_EQ(r.mNumSamples, 4000);
    EXPECT_FALSE(r.mEstimatedNormals);
    EXPECT_NEAR(r.mSpacing, std::sqrt(4.0 * M_PI / 4000.0), 0.02);
    EXPECT_DOUBLE_EQ(r.mCellSize, r.mSpacing);
    EXPECT_EQ(r.mNumComponents, 1);
    EXPECT_EQ(rs_open_edges(r.mMesh), 0u);
    EXPECT_GT(rs_volume6(r.mMesh), 0.0);  // wound outward
    const HausdorffResult h = hausdorff_distance(r.mMesh, icosa_sphere(1.0, 5));
    EXPECT_LT(h.mDistance, r.mSpacing);
    EXPECT_LT(h.mDistance, 0.25 * r.mSpacing);  // in fact far closer
}

TEST(ReconstructSurface, RandomSamplesClose) {
    // Uniformly random samples clump: their median nearest-neighbour distance
    // is about half the gaps between them. The spacing is taken to the fourth
    // nearest sample so the trim and the band bridge those gaps.
    std::mt19937_64 rng(11);
    std::normal_distribution<double> g;
    const std::size_t n = 6000;
    Mesh m;
    NDArray p(DType::Float64, {n, std::size_t{3}});
    for (std::size_t i = 0; i < n; ++i) {
        double v[3] = {g(rng), g(rng), g(rng)};
        const double r = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        for (std::size_t a = 0; a < 3; ++a)
            p.As<double>()[i * 3 + a] = v[a] / r;
    }
    NDArray nrm = p;
    m.AssignPoints(std::move(p));
    m.AddPointData("normals", std::move(nrm));
    const ReconstructSurfaceResult r = reconstruct_surface(m);
    EXPECT_EQ(r.mNumComponents, 1);
    EXPECT_EQ(rs_open_edges(r.mMesh), 0u);
    EXPECT_LT(hausdorff_distance(r.mMesh, icosa_sphere(1.0, 5)).mDistance, 0.1 * r.mSpacing);
}

TEST(ReconstructSurface, OneNeighbourIsHoppesField) {
    ReconstructSurfaceOptions o;
    o.mK = 1;
    const ReconstructSurfaceResult r = reconstruct_surface(rs_sphere(3000, 2.0), o);
    EXPECT_EQ(r.mNumComponents, 1);
    EXPECT_EQ(rs_open_edges(r.mMesh), 0u);
    EXPECT_LT(hausdorff_distance(r.mMesh, icosa_sphere(2.0, 5)).mDistance, r.mSpacing);
}

TEST(ReconstructSurface, EstimatesMissingNormals) {
    const ReconstructSurfaceResult r = reconstruct_surface(rs_sphere(3000, 1.0, false));
    EXPECT_TRUE(r.mEstimatedNormals);
    EXPECT_EQ(r.mNumComponents, 1);
    EXPECT_GT(rs_volume6(r.mMesh), 0.0);
    EXPECT_LT(hausdorff_distance(r.mMesh, icosa_sphere(1.0, 5)).mDistance, r.mSpacing);
    ReconstructSurfaceOptions o;
    o.mEstimateNormals = false;
    EXPECT_THROW(reconstruct_surface(rs_sphere(300, 1.0, false), o), std::invalid_argument);
}

TEST(ReconstructSurface, OpenScanStaysOpen) {
    // The upper hemisphere: the rim is trimmed, not capped.
    const Mesh full = rs_sphere(6000, 1.0);
    const double* p = full.Points().As<double>();
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < full.NumPoints(); ++i)
        if (p[i * 3 + 2] > 0.0)
            keep.push_back(i);
    Mesh half;
    NDArray hp(DType::Float64, {keep.size(), std::size_t{3}});
    NDArray hn(DType::Float64, {keep.size(), std::size_t{3}});
    for (std::size_t j = 0; j < keep.size(); ++j)
        for (std::size_t a = 0; a < 3; ++a) {
            hp.As<double>()[j * 3 + a] = p[keep[j] * 3 + a];
            hn.As<double>()[j * 3 + a] = p[keep[j] * 3 + a];
        }
    half.AssignPoints(std::move(hp));
    half.AddPointData("normals", std::move(hn));
    const ReconstructSurfaceResult r = reconstruct_surface(half);
    EXPECT_EQ(r.mNumComponents, 1);
    EXPECT_GT(r.mNumUndefined, 0);
    EXPECT_GT(rs_open_edges(r.mMesh), 0u);
    // Nothing reaches far below the rim: at most the trim distance and a cell.
    const double* q = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        EXPECT_GT(q[i * 3 + 2], -(r.mMaxDistance + r.mCellSize));
        const double rad = std::sqrt(q[i * 3] * q[i * 3] + q[i * 3 + 1] * q[i * 3 + 1] +
                                     q[i * 3 + 2] * q[i * 3 + 2]);
        EXPECT_NEAR(rad, 1.0, r.mSpacing);
    }
}

TEST(ReconstructSurface, LatticeChoicesAndUnit) {
    Mesh m = rs_sphere(2000, 1.0);
    detail::set_mesh_length_unit(m, 1e-3);
    ReconstructSurfaceOptions o;
    o.mResolution = 20;
    const ReconstructSurfaceResult a = reconstruct_surface(m, o);
    EXPECT_NEAR(a.mCellSize, 2.0 / 20.0, 1e-3);
    EXPECT_EQ(a.mNumComponents, 1);
    EXPECT_DOUBLE_EQ(detail::mesh_length_unit(a.mMesh), 1e-3);
    o.mResolution = 0;
    o.mCellSize = 0.05;
    o.mMaxDistance = 0.3;
    const ReconstructSurfaceResult b = reconstruct_surface(m, o);
    EXPECT_DOUBLE_EQ(b.mCellSize, 0.05);
    EXPECT_DOUBLE_EQ(b.mMaxDistance, 0.3);
    EXPECT_GT(b.mNumTriangles, a.mNumTriangles);
    EXPECT_EQ(rs_open_edges(b.mMesh), 0u);
}

TEST(ReconstructSurface, SameOnEveryThreadCount) {
    const Mesh m = rs_sphere(2500, 1.0);
    std::vector<double> ref_p;
    std::vector<std::int64_t> ref_c;
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    const int saved = mio_test::max_threads();
    for (int threads : {1, 2, 4}) {
        mio_test::set_threads(threads);
#else
    for (int threads : {1, 2}) {
        (void)threads;
#endif
        const ReconstructSurfaceResult r = reconstruct_surface(m);
        const double* p = r.mMesh.Points().As<double>();
        const auto* c = r.mMesh.Cells(0).Conn().As<std::int64_t>();
        std::vector<double> vp(p, p + r.mMesh.NumPoints() * 3);
        std::vector<std::int64_t> vc(c, c + r.mMesh.Cells(0).NumCells() * 3);
        if (ref_p.empty()) {
            ref_p = vp;
            ref_c = vc;
        }
        EXPECT_EQ(vp, ref_p);
        EXPECT_EQ(vc, ref_c);
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(saved);
#endif
}

TEST(ReconstructSurface, Refusals) {
    EXPECT_THROW(reconstruct_surface(Mesh{}), std::invalid_argument);
    const Mesh m = rs_sphere(200, 1.0);
    ReconstructSurfaceOptions o;
    o.mK = 0;
    EXPECT_THROW(reconstruct_surface(m, o), std::invalid_argument);
    o = {};
    o.mCellSize = 0.1;
    o.mResolution = 10;
    EXPECT_THROW(reconstruct_surface(m, o), std::invalid_argument);
    o = {};
    o.mMaxDistance = -1.0;
    EXPECT_THROW(reconstruct_surface(m, o), std::invalid_argument);
    o = {};
    o.mCellSize = 1e-3;
    o.mMaxCells = 1000;
    EXPECT_THROW(reconstruct_surface(m, o), std::invalid_argument);
    o = {};
    o.mNormalsName = "missing";
    o.mEstimateNormals = false;
    EXPECT_THROW(reconstruct_surface(m, o), std::invalid_argument);
}
