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
// detail/tri_bvh: the fast winding number against the exact sum (β = 0) on a
// sphere, a torus, an open hemisphere, two overlapping spheres, an inverted
// sphere and a self-intersecting soup; its value inside and outside; the
// solid angles of a closed surface adding up to 4π; ray hits against a brute
// force; and bit-identical results on every thread count.

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

// External includes
#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

// Project includes
#include "analytic_surfaces.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "../../src/cpp/src/detail/tri_bvh.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::detail::Vec3;

// Seeded points in a box, from raw mt19937 output.
std::vector<Vec3> tb_points(std::size_t N, double Lo, double Hi, std::uint32_t Seed) {
    std::mt19937 rng(Seed);
    std::vector<Vec3> p(N);
    for (auto& q : p)
        for (double& c : q)
            c = Lo + (Hi - Lo) * (static_cast<double>(rng()) / 4294967295.0);
    return p;
}

// The same surface with every triangle's winding reversed.
mio::Mesh tb_flipped(const mio::Mesh& rM) {
    std::vector<std::vector<double>> pts;
    const double* p = rM.Points().As<double>();
    for (std::size_t i = 0; i < rM.NumPoints(); ++i)
        pts.push_back({p[3 * i], p[3 * i + 1], p[3 * i + 2]});
    std::vector<std::vector<std::int64_t>> tris;
    const auto block = rM.Cells(0);
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < block.NumCells(); ++t)
        tris.push_back({c[3 * t], c[3 * t + 2], c[3 * t + 1]});
    return mt::make_mesh(std::move(pts), "triangle", std::move(tris));
}

// The sphere minus every `Every`-th triangle.
mio::Mesh tb_holed(const mio::Mesh& rM, std::size_t Every) {
    std::vector<std::vector<double>> pts;
    const double* p = rM.Points().As<double>();
    for (std::size_t i = 0; i < rM.NumPoints(); ++i)
        pts.push_back({p[3 * i], p[3 * i + 1], p[3 * i + 2]});
    std::vector<std::vector<std::int64_t>> tris;
    const auto block = rM.Cells(0);
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < block.NumCells(); ++t)
        if (t % Every != 0)
            tris.push_back({c[3 * t], c[3 * t + 1], c[3 * t + 2]});
    return mt::make_mesh(std::move(pts), "triangle", std::move(tris));
}

void tb_expect_fast_matches_exact(const mio::Mesh& rSkin, double Tol) {
    const auto soup = mio::detail::build_triangle_soup(rSkin, "");
    const auto tree = mio::detail::build_tri_bvh(soup);
    for (const Vec3& q : tb_points(300, -2.0, 2.0, 7u)) {
        const double exact = mio::detail::fast_winding(tree, q, 0.0);
        EXPECT_NEAR(mio::detail::fast_winding(tree, q, 2.0), exact, Tol);
        EXPECT_NEAR(mio::detail::fast_winding(tree, q, 1e9), exact, 1e-12);
    }
}

}  // namespace

TEST(TriBvh, SolidAnglesOfAClosedSurfaceAddUpToFourPi) {
    const mio::Mesh sphere = mio::icosa_sphere(1.0, 2);
    const auto soup = mio::detail::build_triangle_soup(sphere, "");
    for (const Vec3& q : {Vec3{0.0, 0.0, 0.0}, Vec3{0.3, -0.2, 0.5}}) {
        double sum = 0.0;
        for (std::size_t t = 0; t < soup.NumTriangles(); ++t)
            sum += mio::detail::tri_solid_angle(soup.mCorners[3 * t], soup.mCorners[3 * t + 1],
                                                soup.mCorners[3 * t + 2], q);
        EXPECT_NEAR(sum, 4.0 * 3.14159265358979323846, 1e-9);
    }
}

TEST(TriBvh, WindingIsOneInsideAndZeroOutside) {
    const mio::Mesh sphere = mio::icosa_sphere(1.0, 3);
    const auto soup = mio::detail::build_triangle_soup(sphere, "");
    const auto tree = mio::detail::build_tri_bvh(soup);
    for (const Vec3& q : tb_points(400, -2.0, 2.0, 3u)) {
        const double r = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
        const double w = mio::detail::fast_winding(tree, q, 2.0);
        if (r < 0.95)
            EXPECT_NEAR(w, 1.0, 1e-2);
        else if (r > 1.05)
            EXPECT_NEAR(w, 0.0, 1e-2);
    }
}

TEST(TriBvh, FastWindingMatchesTheExactSum) {
    tb_expect_fast_matches_exact(mio::icosa_sphere(1.2, 3), 5e-3);
    tb_expect_fast_matches_exact(analytic::torus(1.0, 0.4, 48, 24), 5e-3);
    tb_expect_fast_matches_exact(tb_holed(mio::icosa_sphere(1.0, 3), 7), 5e-3);
}

TEST(TriBvh, HolesOverlapsAndInversionDegradeGracefully) {
    const mio::Mesh a = mio::icosa_sphere(1.0, 3);
    const mio::Mesh b = mio::icosa_sphere(1.0, 3, {{0.8, 0.0, 0.0}});
    mio::MergeOptions mo;
    const mio::Mesh both = mio::merge({&a, &b}, mo).mMesh;
    const auto soup = mio::detail::build_triangle_soup(both, "");
    const auto tree = mio::detail::build_tri_bvh(soup);
    EXPECT_NEAR(mio::detail::fast_winding(tree, {0.4, 0.0, 0.0}, 2.0), 2.0, 1e-2);
    EXPECT_NEAR(mio::detail::fast_winding(tree, {-0.6, 0.0, 0.0}, 2.0), 1.0, 1e-2);
    // A self-intersecting soup of two overlapping spheres still matches the sum.
    tb_expect_fast_matches_exact(both, 5e-3);
    // Inverted: -1 inside.
    const auto flipped = mio::detail::build_triangle_soup(tb_flipped(a), "");
    const auto ftree = mio::detail::build_tri_bvh(flipped);
    EXPECT_NEAR(mio::detail::fast_winding(ftree, {0.1, 0.2, 0.0}, 2.0), -1.0, 1e-2);
    // A holed sphere: still well inside at the centre, outside far away.
    const auto holed = mio::detail::build_triangle_soup(tb_holed(a, 11), "");
    const auto htree = mio::detail::build_tri_bvh(holed);
    EXPECT_GT(mio::detail::fast_winding(htree, {0.0, 0.0, 0.0}, 2.0), 0.8);
    EXPECT_LT(mio::detail::fast_winding(htree, {3.0, 0.0, 0.0}, 2.0), 0.05);
}

TEST(TriBvh, RaysHitWhatABruteForceHits) {
    const mio::Mesh skin = analytic::torus(1.0, 0.4, 32, 16);
    const auto soup = mio::detail::build_triangle_soup(skin, "");
    const auto tree = mio::detail::build_tri_bvh(soup);
    const auto origins = tb_points(200, -1.5, 1.5, 11u);
    const auto dirs = tb_points(200, -1.0, 1.0, 13u);
    for (std::size_t i = 0; i < origins.size(); ++i) {
        const Vec3& o = origins[i];
        const Vec3& d = dirs[i];
        const auto hit = mio::detail::tri_bvh_ray(tree, o, d, 1e-9, 1e30);
        double best = 1e30;
        std::int64_t best_t = -1;
        for (std::size_t t = 0; t < soup.NumTriangles(); ++t) {
            const Vec3* c = &soup.mCorners[3 * t];
            const Vec3 e1 = mio::detail::vec3_sub(c[1], c[0]),
                       e2 = mio::detail::vec3_sub(c[2], c[0]);
            const Vec3 p = mio::detail::vec3_cross(d, e2);
            const double det = mio::detail::vec3_dot(e1, p);
            if (std::abs(det) < 1e-300)
                continue;
            const Vec3 s = mio::detail::vec3_sub(o, c[0]);
            const double u = mio::detail::vec3_dot(s, p) / det;
            const Vec3 q = mio::detail::vec3_cross(s, e1);
            const double v = mio::detail::vec3_dot(d, q) / det;
            const double th = mio::detail::vec3_dot(e2, q) / det;
            if (u >= -1e-10 && v >= -1e-10 && u + v <= 1 + 1e-10 && th > 1e-9 && th < best) {
                best = th;
                best_t = static_cast<std::int64_t>(t);
            }
        }
        EXPECT_EQ(hit.mTriangle, best_t);
        if (best_t >= 0)
            EXPECT_NEAR(hit.mT, best, 1e-9 * (1.0 + best));
    }
}

TEST(TriBvh, WindingIsTheSameOnEveryThreadCount) {
    const mio::Mesh skin = analytic::torus(1.0, 0.4, 48, 24);
    const auto soup = mio::detail::build_triangle_soup(skin, "");
    const auto tree = mio::detail::build_tri_bvh(soup);
    const auto pts = tb_points(5000, -1.6, 1.6, 17u);
    std::vector<double> ref;
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const auto w = mio::detail::fast_windings(tree, pts, 2.0);
        if (ref.empty())
            ref = w;
        EXPECT_EQ(w, ref);
    }
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
#endif
}
