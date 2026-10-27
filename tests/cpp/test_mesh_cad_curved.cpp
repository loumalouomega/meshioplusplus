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

// Curved surface meshes from mesh_cad_surface (roadmap 6.6): quadratic and
// cubic triangles whose nodes lie on the exact surface and edges.

// System includes
#include <array>
#include <cmath>
#include <cstdint>
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
#include "../../src/cpp/src/detail/curved_triangle.hpp"
#include "../../src/cpp/src/operations/mesh_cad_testing.hpp"
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"

namespace {

namespace mio = meshioplusplus;
namespace mdt = meshioplusplus::detail;

constexpr double kPi = 3.14159265358979323846;

mio::MeshCadSurfaceOptions curved(double H, int Order, double MinAngle = 0.0) {
    mio::MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    o.mOrder = Order;
    o.mMinAngleDeg = MinAngle;
    return o;
}

std::uint64_t digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

/// The elements of the mesh's first block, their nodes as points.
std::vector<std::vector<mdt::CurvedVec3>> elements(const mio::Mesh& rM, int Order) {
    const auto cells = rM.Cells(0);
    const std::size_t nn = static_cast<std::size_t>(mdt::curved_triangle_num_nodes(Order));
    const double* x = rM.Points().As<double>();
    std::vector<std::vector<mdt::CurvedVec3>> out(cells.NumCells());
    for (std::size_t c = 0; c < cells.NumCells(); ++c)
        for (std::size_t i = 0; i < nn; ++i) {
            const auto id =
                static_cast<std::size_t>(mio::detail::read_int(cells.Conn(), nn * c + i));
            out[c].push_back({x[3 * id], x[3 * id + 1], x[3 * id + 2]});
        }
    return out;
}

/// Area and enclosed volume of the curved mesh, by a collapsed 8 x 8
/// Gauss-Legendre rule on each element (the area element's components are
/// J against the three axes).
std::array<double, 2> area_volume(const mio::Mesh& rM, int Order) {
    static const double kX[8] = {0.0198550717512319, 0.1016667612931866, 0.2372337950418355,
                                 0.4082826787521751, 0.5917173212478249, 0.7627662049581645,
                                 0.8983332387068134, 0.9801449282487681};
    static const double kW[8] = {0.0506142681451881, 0.1111905172266872, 0.1568533229389436,
                                 0.1813418916891810, 0.1813418916891810, 0.1568533229389436,
                                 0.1111905172266872, 0.0506142681451881};
    double area = 0.0, volume = 0.0;
    for (const auto& r_e : elements(rM, Order))
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j) {
                const double u = kX[i], v = kX[j];
                const double xi = u, eta = v * (1.0 - u);
                const double w = kW[i] * kW[j] * (1.0 - u);
                mdt::CurvedVec3 n{};
                for (std::size_t k = 0; k < 3; ++k) {
                    mdt::CurvedVec3 axis{0.0, 0.0, 0.0};
                    axis[k] = 1.0;
                    n[k] = mdt::curved_triangle_jacobian(Order, r_e.data(), axis, xi, eta);
                }
                const auto x = mdt::curved_triangle_point(Order, r_e.data(), xi, eta);
                area += w * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                volume += w * (x[0] * n[0] + x[1] * n[1] + x[2] * n[2]) / 3.0;
            }
    return {area, volume};
}

TEST(MeshCadCurved, EveryNodeLiesOnTheSphere) {
    // Interior edge nodes are the surface at parameter points, boundary ones
    // the edges' arc-length points, the centres the surface at the
    // parameter centroid: on the unit sphere all of them, to rounding.
    for (const int order : {2, 3}) {
        const auto r = mio::mesh_cad_surface(brep_solids::sphere(1.0), curved(0.3, order));
        ASSERT_EQ(r.mNumFailed, 0);
        EXPECT_EQ(r.mMesh.Cells(0).Type(), order == 2 ? "triangle6" : "triangle10");
        EXPECT_EQ(r.mMesh.Cells(1).Type(), order == 2 ? "line3" : "line4");
        const double* x = r.mMesh.Points().As<double>();
        for (std::size_t p = 0; p < r.mMesh.NumPoints(); ++p)
            EXPECT_NEAR(std::hypot(x[3 * p], x[3 * p + 1], x[3 * p + 2]), 1.0, 1e-12) << p;
        EXPECT_EQ(r.mNumOpenEdges, 0);
        EXPECT_EQ(r.mNumMisorientedEdges, 0);
        EXPECT_EQ(r.mNumCurvedStraightened, 0);
    }
}

TEST(MeshCadCurved, NodesAreSharedAndCornersAreTheLinearMesh) {
    // The corners are the linear mesh's points and triangles; every edge
    // node of a closed solid serves two elements (the line block's too on
    // an edge), every centre one.
    const mio::BrepInfo cyl = brep_solids::capped_cylinder(1.0, 2.0);
    const auto lin = mio::mesh_cad_surface(cyl, curved(0.3, 1));
    for (const int order : {2, 3}) {
        const auto r = mio::mesh_cad_surface(cyl, curved(0.3, order));
        ASSERT_EQ(r.mNumFailed, 0);
        const auto tri = r.mMesh.Cells(0);
        const auto lt = lin.mMesh.Cells(0);
        ASSERT_EQ(tri.NumCells(), lt.NumCells());
        const std::size_t nn = static_cast<std::size_t>(mdt::curved_triangle_num_nodes(order));
        for (std::size_t c = 0; c < tri.NumCells(); ++c)
            for (std::size_t k = 0; k < 3; ++k)
                EXPECT_EQ(mio::detail::read_int(tri.Conn(), nn * c + k),
                          mio::detail::read_int(lt.Conn(), 3 * c + k));
        const double* a = lin.mMesh.Points().As<double>();
        const double* b = r.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < 3 * lin.mMesh.NumPoints(); ++i)
            EXPECT_EQ(a[i], b[i]);
        std::map<std::int64_t, int> uses;
        for (std::size_t c = 0; c < tri.NumCells(); ++c)
            for (std::size_t k = 3; k < nn; ++k)
                ++uses[mio::detail::read_int(tri.Conn(), nn * c + k)];
        for (const auto& [id, n] : uses) {
            const bool centre =
                static_cast<std::size_t>(id) >= lin.mMesh.NumPoints() && order == 3 && n == 1;
            EXPECT_TRUE(n == 2 || centre) << id << " used " << n;
        }
        // Every node is on the cylinder: the band, or a cap.
        const double* x = r.mMesh.Points().As<double>();
        for (std::size_t p = 0; p < r.mMesh.NumPoints(); ++p) {
            const double rad = std::hypot(x[3 * p], x[3 * p + 1]);
            const double z = x[3 * p + 2];
            const bool band = std::fabs(rad - 1.0) < 1e-12 && z > -1e-12 && z < 2.0 + 1e-12;
            const bool cap =
                (std::fabs(z) < 1e-12 || std::fabs(z - 2.0) < 1e-12) && rad < 1.0 + 1e-12;
            EXPECT_TRUE(band || cap) << p;
        }
    }
}

TEST(MeshCadCurved, AreaAndVolumeConvergeFasterThanLinear) {
    // Roadmap 6.6's probe: on the cylinder and the sphere, the error of the
    // area and the enclosed volume against the exact ones, by order, as h
    // halves from 0.4 to 0.1. Measured: about 2 for linear triangles, 3.9
    // to 4.3 for quadratic and 4.0 to 4.3 for cubic, the cubic errors 6 to
    // 20 times smaller than the quadratic (its interior nodes at parameter
    // fractions, roadmap 6.6.2, are what keeps it near 4). Pinned below
    // with margin: at least 1.8, 3.5 and 3.5, and cubic below quadratic.
    const double cyl_area = 2.0 * kPi * 2.0 + 2.0 * kPi, cyl_volume = 2.0 * kPi;
    const double sph_area = 4.0 * kPi, sph_volume = 4.0 / 3.0 * kPi;
    std::map<std::pair<int, int>, std::array<double, 4>> errors;  // (order, step)
    for (const int order : {1, 2, 3}) {
        int step = 0;
        double prev[4] = {0, 0, 0, 0};
        for (const double h : {0.4, 0.2, 0.1}) {
            const auto c = area_volume(
                mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), curved(h, order))
                    .mMesh,
                order);
            const auto s = area_volume(
                mio::mesh_cad_surface(brep_solids::sphere(1.0), curved(h, order)).mMesh, order);
            const double err[4] = {std::fabs(c[0] - cyl_area), std::fabs(c[1] - cyl_volume),
                                   std::fabs(s[0] - sph_area), std::fabs(s[1] - sph_volume)};
            std::printf(
                "[curved] order %d h %.2f  cylinder area %.3e volume %.3e  sphere area "
                "%.3e volume %.3e",
                order, h, err[0], err[1], err[2], err[3]);
            if (prev[0] > 0.0)
                std::printf("  orders %.2f %.2f %.2f %.2f", std::log2(prev[0] / err[0]),
                            std::log2(prev[1] / err[1]), std::log2(prev[2] / err[2]),
                            std::log2(prev[3] / err[3]));
            std::printf("\n");
            for (int k = 0; k < 4; ++k) {
                if (prev[k] > 0.0)
                    EXPECT_GT(std::log2(prev[k] / err[k]), order == 1 ? 1.8 : 3.5)
                        << "order " << order << ", h " << h << ", measure " << k;
                prev[k] = err[k];
            }
            errors[{order, step++}] = {err[0], err[1], err[2], err[3]};
        }
    }
    for (int step = 0; step < 3; ++step)
        for (int k = 0; k < 4; ++k) {
            EXPECT_LT((errors[{2, step}][k]), (errors[{1, step}][k])) << step << " " << k;
            EXPECT_LT((errors[{3, step}][k]), (errors[{2, step}][k])) << step << " " << k;
        }
}

TEST(MeshCadCurved, AnElementTheCurvedBoundaryFoldsIsKeptStraightSided) {
    // At h = 1.75 the cylinder's faces hold a few large triangles, and some
    // cubic elements fold (their Jacobian not positive everywhere against
    // their face's normal): they are reported, warned and made
    // straight-sided, and every element written is valid.
    const auto r = mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), curved(1.75, 3));
    ASSERT_EQ(r.mNumFailed, 0);
    EXPECT_GT(r.mNumCurvedStraightened, 0);
    bool warned = false;
    for (const std::string& r_w : r.mWarnings)
        warned = warned || r_w.find("kept straight-sided") != std::string::npos;
    EXPECT_TRUE(warned);
    std::int64_t straight = 0;
    for (const auto& r_e : elements(r.mMesh, 3)) {
        // The written element against its corners' own normal.
        const mdt::CurvedVec3 a{r_e[1][0] - r_e[0][0], r_e[1][1] - r_e[0][1],
                                r_e[1][2] - r_e[0][2]};
        const mdt::CurvedVec3 b{r_e[2][0] - r_e[0][0], r_e[2][1] - r_e[0][1],
                                r_e[2][2] - r_e[0][2]};
        const mdt::CurvedVec3 n{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                                a[0] * b[1] - a[1] * b[0]};
        EXPECT_TRUE(mdt::curved_triangle_check(3, r_e.data(), n).mValid);
        // A straight-sided element: its edge nodes on its chords.
        bool on_chords = true;
        for (std::size_t k = 0; k < 3; ++k)
            for (int j = 1; j < 3; ++j) {
                const auto& p = r_e[3 + 2 * k + static_cast<std::size_t>(j - 1)];
                for (std::size_t c = 0; c < 3; ++c) {
                    const double want =
                        r_e[k][c] + (double(j) / 3.0) * (r_e[(k + 1) % 3][c] - r_e[k][c]);
                    on_chords = on_chords && std::fabs(p[c] - want) < 1e-12;
                }
            }
        straight += on_chords ? 1 : 0;
    }
    EXPECT_EQ(straight, r.mNumCurvedStraightened);
    // At h = 0.3 nothing folds.
    EXPECT_EQ(mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), curved(0.3, 3))
                  .mNumCurvedStraightened,
              0);
}

TEST(MeshCadCurved, UntanglingKeepsFoldedElementsCurved) {
    // Roadmap 6.6.3 (Toulorge et al. 2013): the holed disc, 0.5 thick, at
    // h = 0.6. Across the thin wall between the hole and the rim some
    // curved elements fold. With `untangle` the caps' own nodes move on
    // their planes first: none is straightened, every element written is
    // valid, every node lies on the solid's faces, and the area is closer.
    const mio::BrepInfo disc = brep_solids::thick_disc_with_hole();
    const auto on_solid = [](const double* pX) {
        return std::fabs(pX[2]) < 1e-12 || std::fabs(pX[2] - 0.5) < 1e-12 ||
               std::fabs(std::hypot(pX[0], pX[1]) - 2.0) < 1e-9 ||
               std::fabs(std::hypot(pX[0] - 0.98, pX[1]) - 1.0) < 1e-9;
    };
    // The surface area: two annuli and the two bands.
    const double exact = 2.0 * (kPi * 4.0 - kPi) + 0.5 * 2.0 * kPi * (2.0 + 1.0);
    mio::MeshCadSurfaceExtras ex;
    ex.mUntangle = true;
    for (const int order : {2, 3}) {
        const auto plain = mio::mesh_cad_surface(disc, curved(0.6, order));
        mio::MeshCadSurfaceExtrasReport x;
        const auto r = mio::mesh_cad_surface(disc, curved(0.6, order), mio::MeshCadSizeSources{},
                                             ex, nullptr, &x);
        ASSERT_EQ(r.mNumFailed, 0);
        const double e_plain = std::fabs(area_volume(plain.mMesh, order)[0] - exact);
        const double e = std::fabs(area_volume(r.mMesh, order)[0] - exact);
        std::printf(
            "[6.6.3] order %d: straightened %lld without untangling, %lld with (%lld "
            "untangled); area error %.3e against %.3e\n",
            order, static_cast<long long>(plain.mNumCurvedStraightened),
            static_cast<long long>(r.mNumCurvedStraightened),
            static_cast<long long>(x.mNumUntangled), e_plain, e);
        EXPECT_GT(plain.mNumCurvedStraightened, 0) << order;
        EXPECT_EQ(r.mNumCurvedStraightened, 0) << order;
        EXPECT_EQ(x.mNumUntangled, plain.mNumCurvedStraightened) << order;
        EXPECT_LT(e, e_plain) << order;
        for (const auto& r_e : elements(r.mMesh, order)) {
            const mdt::CurvedVec3 a{r_e[1][0] - r_e[0][0], r_e[1][1] - r_e[0][1],
                                    r_e[1][2] - r_e[0][2]};
            const mdt::CurvedVec3 b{r_e[2][0] - r_e[0][0], r_e[2][1] - r_e[0][1],
                                    r_e[2][2] - r_e[0][2]};
            const mdt::CurvedVec3 n{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                                    a[0] * b[1] - a[1] * b[0]};
            EXPECT_TRUE(mdt::curved_triangle_check(order, r_e.data(), n).mValid);
        }
        const double* p = r.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
            EXPECT_TRUE(on_solid(p + 3 * i)) << i;
            // Deterministic across threads.
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
        const std::uint64_t ref = digest(r.mMesh);
        const int saved = omp_get_max_threads();
        for (int threads : {1, 2, 4}) {
            omp_set_num_threads(threads);
            EXPECT_EQ(digest(mio::mesh_cad_surface(disc, curved(0.6, order),
                                                   mio::MeshCadSizeSources{}, ex)
                                 .mMesh),
                      ref)
                << threads;
        }
        omp_set_num_threads(saved);
#endif
    }
    // A fold at a corner where an element's two sides run along one edge
    // (the cylinder's caps at h = 1.75, above: each cap is two triangles
    // with every corner on its rim) is a straight angle there, where the
    // Jacobian vanishes whatever the nodes do: with untangling the face is
    // meshed again with an interior grid, so every element has a corner
    // inside it, and none is straightened.
    const mio::BrepInfo cyl = brep_solids::capped_cylinder(1.0, 2.0);
    mio::MeshCadSurfaceExtrasReport xc;
    const auto rc =
        mio::mesh_cad_surface(cyl, curved(1.75, 3), mio::MeshCadSizeSources{}, ex, nullptr, &xc);
    const auto pc = mio::mesh_cad_surface(cyl, curved(1.75, 3));
    std::printf(
        "[6.6.3] capped cylinder h 1.75 order 3: %lld straightened, %lld with untangling "
        "(%lld triangles against %lld)\n",
        static_cast<long long>(pc.mNumCurvedStraightened),
        static_cast<long long>(rc.mNumCurvedStraightened), static_cast<long long>(rc.mNumTriangles),
        static_cast<long long>(pc.mNumTriangles));
    EXPECT_GT(pc.mNumCurvedStraightened, 0);
    EXPECT_EQ(rc.mNumCurvedStraightened, 0);
    EXPECT_GT(rc.mNumTriangles, pc.mNumTriangles);
    EXPECT_EQ(rc.mNumOpenEdges, 0);
    // Where nothing folds the mesh is the plain one, byte for byte.
    EXPECT_EQ(
        digest(mio::mesh_cad_surface(cyl, curved(0.3, 3), mio::MeshCadSizeSources{}, ex).mMesh),
        digest(mio::mesh_cad_surface(cyl, curved(0.3, 3)).mMesh));
    // Linear elements have nothing to untangle: refused.
    EXPECT_THROW((void)mio::mesh_cad_surface(cyl, curved(0.3, 1), mio::MeshCadSizeSources{}, ex),
                 std::invalid_argument);
}

TEST(MeshCadCurved, TheMeshDoesNotDependOnTheFaceOrderOrThreads) {
    const mio::BrepInfo cyl = brep_solids::capped_cylinder(1.0, 2.0);
    for (const int order : {2, 3}) {
        const std::uint64_t ref =
            digest(mio::mesh_cad_surface(cyl, curved(0.2, order, 30.0)).mMesh);
        for (const std::uint64_t seed : {1u, 7u, 12345u}) {
            mio::detail::mesh_cad_set_face_order_seed(seed);
            EXPECT_EQ(digest(mio::mesh_cad_surface(cyl, curved(0.2, order, 30.0)).mMesh), ref)
                << seed;
        }
        mio::detail::mesh_cad_set_face_order_seed(0);
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
        const int saved = omp_get_max_threads();
        for (int threads : {1, 2, 4}) {
            omp_set_num_threads(threads);
            EXPECT_EQ(digest(mio::mesh_cad_surface(cyl, curved(0.2, order, 30.0)).mMesh), ref)
                << threads;
        }
        omp_set_num_threads(saved);
#endif
    }
}

}  // namespace
