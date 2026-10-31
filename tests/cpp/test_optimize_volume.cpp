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
// Tests for predicate-free ODT remeshing (operations/optimize_volume.cpp). The
// two analytic flip fixtures below were found by an independent numpy search
// for a 5-point configuration whose edge/face bistellar re-triangulation
// strictly raises the worst scaled Jacobian -- so a wrong winding, a wrong
// convexity test, or a wrong acceptance comparison shows up as the flip NOT
// firing (num_*_flips == 0) or as an inverted output cell, not merely a
// cosmetic difference. The "remeshing not smoothing" oracle is the load-bearing
// one: on a fixture with no free interior vertex, relocation alone cannot move
// the worst quality at all, so any improvement there is proof the connectivity
// changed.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/surface.hpp"

namespace {

using meshioplusplus::compute_stats;
using meshioplusplus::extract_surface;
using meshioplusplus::Mesh;
using meshioplusplus::optimize_volume;
using meshioplusplus::OptimizeVolumeObjective;
using meshioplusplus::OptimizeVolumeOptions;
using meshioplusplus::OptimizeVolumeResult;

// The two tets (a,b,c,d) and (a,b,c,e) share face abc; the 2-3 flip to three
// tets around edge de raises the worst scaled Jacobian from ~0.066 to ~0.362.
Mesh two_tet_23_fixture() {
    return mt::make_mesh({{0.0, 0.0, 0.0},
                          {1.0, 0.0, 0.0},
                          {0.3, 0.9, 0.0},
                          {0.2252, 0.2808, 0.6977},
                          {0.4659, 0.3149, -0.0367}},
                         "tetra", {{0, 1, 2, 3}, {0, 1, 2, 4}});
}

// Three tets around edge (3,4); the 3-2 flip to two caps over triangle (0,1,2)
// raises the worst scaled Jacobian from ~0.299 to ~0.904.
Mesh three_tet_32_fixture() {
    return mt::make_mesh({{0.0, 0.0, 0.0},
                          {1.0, 0.0, 0.0},
                          {0.3, 0.9, 0.0},
                          {0.45, 0.35, 0.9},
                          {0.45, 0.35, -0.9}},
                         "tetra", {{3, 4, 0, 1}, {3, 4, 1, 2}, {3, 4, 2, 0}});
}

// The sorted set of boundary triangle facets (as sorted node-id triples). Point
// indices are invariant across optimize_volume, so this is directly comparable
// before/after.
std::set<std::array<std::int64_t, 3>> boundary_facets(const Mesh& rMesh) {
    std::set<std::array<std::int64_t, 3>> out;
    const Mesh surf = extract_surface(rMesh);
    for (const auto cb : surf.CellRange()) {
        if (cb.Type() != std::string("triangle"))
            continue;
        const auto& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::array<std::int64_t, 3> f{meshioplusplus::detail::read_int(conn, c * 3),
                                          meshioplusplus::detail::read_int(conn, c * 3 + 1),
                                          meshioplusplus::detail::read_int(conn, c * 3 + 2)};
            std::sort(f.begin(), f.end());
            out.insert(f);
        }
    }
    return out;
}

double signed_vol6(const Mesh& rMesh, std::int64_t a, std::int64_t b, std::int64_t c,
                   std::int64_t d) {
    const auto& p = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    auto rd = [&](std::int64_t v, std::size_t k) {
        return meshioplusplus::detail::read_double(p, static_cast<std::size_t>(v) * dim + k);
    };
    const double e1[3] = {rd(b, 0) - rd(a, 0), rd(b, 1) - rd(a, 1), rd(b, 2) - rd(a, 2)};
    const double e2[3] = {rd(c, 0) - rd(a, 0), rd(c, 1) - rd(a, 1), rd(c, 2) - rd(a, 2)};
    const double e3[3] = {rd(d, 0) - rd(a, 0), rd(d, 1) - rd(a, 1), rd(d, 2) - rd(a, 2)};
    const double cx = e2[1] * e3[2] - e2[2] * e3[1];
    const double cy = e2[2] * e3[0] - e2[0] * e3[2];
    const double cz = e2[0] * e3[1] - e2[1] * e3[0];
    return e1[0] * cx + e1[1] * cy + e1[2] * cz;
}

// A regular tetrahedron split by an interior vertex (index 4, placed off the
// true centroid so it is genuinely sub-optimal) into four sub-tets. The four
// corners are on the boundary; vertex 4 is the one free interior node, so this
// exercises the relocation half while staying a guaranteed-valid mesh.
Mesh centroid_split_fixture() {
    return mt::make_mesh({{0.0, 0.0, 0.0},
                          {1.0, 0.0, 0.0},
                          {0.5, 0.866, 0.0},
                          {0.5, 0.289, 0.816},
                          {0.55, 0.25, 0.30}},  // interior, jittered off centroid
                         "tetra", {{4, 1, 2, 3}, {0, 4, 2, 3}, {0, 1, 4, 3}, {0, 1, 2, 4}});
}

TEST(OptimizeVolume, TwentyThreeFlipImprovesWorstQuality) {
    OptimizeVolumeOptions o;
    o.mObjective = OptimizeVolumeObjective::Energy;  // the scaled-Jacobian oracle
    o.mRelocate = false;                             // isolate the flip
    const OptimizeVolumeResult r = optimize_volume(two_tet_23_fixture(), o);
    EXPECT_EQ(r.mNum23Flips, 1);
    EXPECT_EQ(r.mNum32Flips, 0);
    EXPECT_EQ(r.mNumTets, 3);
    EXPECT_GT(r.mMinQualityAfter, r.mMinQualityBefore + 0.2);
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
}

// The winding oracle named in optimize_volume.cpp: on a convex union the three
// new tets are all positively oriented. A wrong (d,e) winding would either
// invert them (caught here) or fail the convexity test and not flip at all.
TEST(OptimizeVolume, TwentyThreeFlipWindingIsPositiveOnAConvexUnion) {
    OptimizeVolumeOptions o;
    o.mObjective = OptimizeVolumeObjective::Energy;
    o.mRelocate = false;
    const OptimizeVolumeResult r = optimize_volume(two_tet_23_fixture(), o);
    ASSERT_EQ(r.mNumTets, 3);
    const auto cv = r.mMesh.Cells(0);
    const auto& conn = cv.Conn();
    for (std::size_t c = 0; c < 3; ++c)
        EXPECT_GT(signed_vol6(r.mMesh, meshioplusplus::detail::read_int(conn, c * 4),
                              meshioplusplus::detail::read_int(conn, c * 4 + 1),
                              meshioplusplus::detail::read_int(conn, c * 4 + 2),
                              meshioplusplus::detail::read_int(conn, c * 4 + 3)),
                  0.0);
}

TEST(OptimizeVolume, ThirtyTwoFlipImprovesWorstQuality) {
    OptimizeVolumeOptions o;
    o.mObjective = OptimizeVolumeObjective::Energy;
    o.mRelocate = false;
    const OptimizeVolumeResult r = optimize_volume(three_tet_32_fixture(), o);
    EXPECT_EQ(r.mNum32Flips, 1);
    EXPECT_EQ(r.mNum23Flips, 0);
    EXPECT_EQ(r.mNumTets, 2);
    EXPECT_GT(r.mMinQualityAfter, r.mMinQualityBefore + 0.2);
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
}

// The load-bearing "remeshing, not smoothing" oracle: every vertex of the 2-tet
// fixture is on its boundary, so relocation with the default boundary pinning
// cannot move a single point -- yet quality still improves, which is only
// possible by changing connectivity.
TEST(OptimizeVolume, FlipsHelpWhereSmoothingCannot) {
    OptimizeVolumeOptions relocate_only;
    relocate_only.mObjective = OptimizeVolumeObjective::Energy;
    relocate_only.mFlip = false;
    const OptimizeVolumeResult sm = optimize_volume(two_tet_23_fixture(), relocate_only);
    EXPECT_EQ(sm.mNumVerticesMoved, 0);
    EXPECT_NEAR(sm.mMinQualityAfter, sm.mMinQualityBefore, 1e-12);

    OptimizeVolumeOptions full;
    full.mObjective = OptimizeVolumeObjective::Energy;
    const OptimizeVolumeResult r = optimize_volume(two_tet_23_fixture(), full);
    EXPECT_GT(r.mNumFlips, 0);
    EXPECT_GT(r.mMinQualityAfter, sm.mMinQualityAfter + 0.2);
}

TEST(OptimizeVolume, MonotoneQualityAndNoInversionOnACluster) {
    OptimizeVolumeOptions o;
    o.mObjective = OptimizeVolumeObjective::Energy;
    const OptimizeVolumeResult r = optimize_volume(centroid_split_fixture(), o);
    EXPECT_GE(r.mMinQualityAfter, r.mMinQualityBefore - 1e-12);
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
}

TEST(OptimizeVolume, BoundaryIsInvariant) {
    const Mesh in = two_tet_23_fixture();
    const auto before = boundary_facets(in);
    const OptimizeVolumeResult r = optimize_volume(in);
    const auto after = boundary_facets(r.mMesh);
    EXPECT_EQ(before, after);
}

TEST(OptimizeVolume, PointSetIsInvariant) {
    const Mesh in = centroid_split_fixture();
    const OptimizeVolumeResult r = optimize_volume(in);
    EXPECT_EQ(r.mMesh.NumPoints(), in.NumPoints());
}

TEST(OptimizeVolume, NoOpWhenBothHalvesDisabled) {
    OptimizeVolumeOptions o;
    o.mRelocate = false;
    o.mFlip = false;
    const OptimizeVolumeResult r = optimize_volume(two_tet_23_fixture(), o);
    EXPECT_EQ(r.mNumFlips, 0);
    EXPECT_EQ(r.mNumVerticesMoved, 0);
    EXPECT_EQ(r.mNumTets, 2);
}

// A Kuhn-split cube of n^3 hexes (6 tets each) whose points are moved by up
// to 0.35 of a lattice step -- boundary points within their face planes -- by
// an integer LCG scaled by a power of two, so the input is bit-identical on
// every platform (no libm). Enough of its tets are poor that both flip passes
// and the relocation fire over several sweeps.
Mesh jittered_cube_fixture(std::size_t n) {
    const std::size_t np = n + 1;
    std::vector<std::array<double, 3>> pts(np * np * np);
    std::uint32_t state = 2463534242u;
    for (std::size_t k = 0; k < np; ++k)
        for (std::size_t j = 0; j < np; ++j)
            for (std::size_t i = 0; i < np; ++i) {
                const std::size_t ijk[3] = {i, j, k};
                auto& p = pts[(k * np + j) * np + i];
                for (int d = 0; d < 3; ++d) {
                    state = state * 1664525u + 1013904223u;
                    // (state >> 16) in [0, 65536): an offset in [-0.35, 0.35) steps.
                    const double off = 0.35 * (static_cast<double>(state >> 16) / 32768.0 - 1.0);
                    const bool on_face = ijk[d] == 0 || ijk[d] == n;
                    p[d] = (static_cast<double>(ijk[d]) + (on_face ? 0.0 : off)) /
                           static_cast<double>(n);
                }
            }
    static const int tets[6][4][3] = {
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {1, 1, 1}},
        {{0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 1}},
        {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}}, {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {1, 1, 1}},
    };
    std::vector<std::vector<std::int64_t>> cells;
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                for (const auto& t : tets) {
                    std::vector<std::int64_t> c;
                    for (const auto& v : t)
                        c.push_back(static_cast<std::int64_t>(((k + v[2]) * np + (j + v[1])) * np +
                                                              (i + v[0])));
                    cells.push_back(std::move(c));
                }
    std::vector<std::vector<double>> coords;
    for (const auto& p : pts)
        coords.push_back({p[0], p[1], p[2]});
    return mt::make_mesh(coords, "tetra", cells);
}

std::uint64_t optimize_volume_digest(const OptimizeVolumeResult& rR) {
    meshioplusplus::bench::MeshDigest d;
    d.Of(rR.mMesh);
    d.U64(static_cast<std::uint64_t>(rR.mNum23Flips));
    d.U64(static_cast<std::uint64_t>(rR.mNum32Flips));
    d.U64(static_cast<std::uint64_t>(rR.mNumVerticesMoved));
    return d.Value();
}

// Byte-identical results from repeated runs: the gate for the sweep's
// relocation context and the sort-based flip tables (roadmap §3).
TEST(OptimizeVolume, ResultIsStableAcrossRepeatedRuns) {
    const Mesh in = jittered_cube_fixture(5);
    for (const auto objective : {OptimizeVolumeObjective::Energy, OptimizeVolumeObjective::Worst}) {
        OptimizeVolumeOptions o;
        o.mObjective = objective;
        const OptimizeVolumeResult first = optimize_volume(in, o);
        ASSERT_GT(first.mNum23Flips, 0);
        ASSERT_GT(first.mNumVerticesMoved, 0);
        const std::uint64_t digest = optimize_volume_digest(first);
        for (int run = 0; run < 3; ++run)
            EXPECT_EQ(optimize_volume_digest(optimize_volume(in, o)), digest) << "run " << run;
    }
}

// The result the implementation produced before its rewrite, pinned: a change
// to the sweep's data structures must not change a single byte. x86-64
// GCC/Clang only -- elsewhere a compiler may contract a*b+c into an FMA, which
// legitimately changes the last bits.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
TEST(OptimizeVolume, ResultMatchesTheGoldenDigest) {
    // The energy objective, byte for byte what the default was up to v18.0.0.
    OptimizeVolumeOptions o;
    o.mObjective = OptimizeVolumeObjective::Energy;
    const OptimizeVolumeResult r = optimize_volume(jittered_cube_fixture(5), o);
    EXPECT_EQ(optimize_volume_digest(r), 0x6c0938c3f829d96bull);
    EXPECT_EQ(r.mNum23Flips, 3);
    EXPECT_EQ(r.mNum32Flips, 2);
    EXPECT_EQ(r.mNumVerticesMoved, 62);
}
#endif

// --- the worst-cell objective (the default since ABI 37) ---------------------

// The least dihedral angle of a mesh's tets, in degrees.
double least_dihedral(const Mesh& rMesh) {
    double m = 180.0;
    const auto cv = rMesh.Cells(0);
    const auto& conn = cv.Conn();
    const auto& p = rMesh.Points();
    auto pt = [&](std::int64_t v) {
        return std::array<double, 3>{meshioplusplus::detail::read_double(p, v * 3),
                                     meshioplusplus::detail::read_double(p, v * 3 + 1),
                                     meshioplusplus::detail::read_double(p, v * 3 + 2)};
    };
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
    static const int kE[6][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                 {1, 2, 0, 3}, {1, 3, 0, 2}, {2, 3, 0, 1}};
    for (std::size_t c = 0; c < cv.NumCells(); ++c) {
        std::array<std::array<double, 3>, 4> q;
        for (int k = 0; k < 4; ++k)
            q[static_cast<std::size_t>(k)] =
                pt(meshioplusplus::detail::read_int(conn, c * 4 + static_cast<std::size_t>(k)));
        for (const auto& e : kE) {
            const auto ab = sub(q[e[1]], q[e[0]]);
            const auto n1 = cross(ab, sub(q[e[2]], q[e[0]]));
            const auto n2 = cross(ab, sub(q[e[3]], q[e[0]]));
            const double cs = dot(n1, n2) / std::sqrt(dot(n1, n1) * dot(n2, n2));
            m = std::min(m, std::acos(std::clamp(cs, -1.0, 1.0)) * 57.29577951308232);
        }
    }
    return m;
}

// The k tets around the interior edge (u, v) = (k, k + 1), whose ring r_0 ..
// r_{k-1} has a triangulation with a much better least angle than the ring's
// tets (found by a search over perturbed rings; the comment gives both).
Mesh edge_ring_fixture(int k) {
    std::vector<std::vector<double>> pts;
    if (k == 4)  // around 25.0 degrees, best triangulation 41.8
        pts = {{0.71, 0.03, -0.04},   {0.39, 0.82, 0.07},  {-0.42, 0.32, 0.05},
               {-0.42, -1.08, -0.02}, {0.12, -0.06, 1.31}, {0.02, -0.12, -0.98}};
    else if (k == 5)  // 11.6 -> 31.3
        pts = {{1.28, 0.65, -0.11},  {-0.84, 0.21, 0.08}, {-0.47, -1.16, 0.06},
               {0.88, -0.96, -0.09}, {0.34, -0.36, 0.1},  {-0.06, -0.16, 1.22},
               {-0.02, -0.16, -0.57}};
    else if (k == 6)  // 8.7 -> 27.0
        pts = {{-0.37, 0.72, -0.08},  {-1.45, 0.47, 0.12}, {-0.85, 0.09, 0.1},
               {-0.16, -0.39, -0.05}, {0.34, -1.12, 0.09}, {1.48, -0.35, -0.05},
               {-0.12, -0.03, 0.98},  {-0.04, -0.17, -0.8}};
    else  // k == 7: 5.1 -> 25.6
        pts = {{0.33, 1.03, 0.11},  {-0.29, 0.3, 0.03},   {-1.02, 0.88, -0.02},
               {-1.04, 0.32, 0.11}, {-0.82, -0.4, -0.02}, {-0.22, -0.4, -0.14},
               {0.17, -0.44, 0.05}, {0.01, -0.12, 0.83},  {-0.04, 0.11, -0.99}};
    std::vector<std::vector<std::int64_t>> tets;
    for (int i = 0; i < k; ++i)
        tets.push_back({k, k + 1, i, (i + 1) % k});
    return mt::make_mesh(pts, "tetra", tets);
}

TEST(OptimizeVolume, EdgeRemovalReplacesARingByItsBestTriangulation) {
    const double best[] = {0.0, 0.0, 0.0, 0.0, 41.78, 31.33, 26.98, 25.61};
    for (int k = 4; k <= 7; ++k) {
        OptimizeVolumeOptions o;
        o.mRelocate = false;  // the ring's vertices are all on its boundary anyway
        const Mesh in = edge_ring_fixture(k);
        const OptimizeVolumeResult r = optimize_volume(in, o);
        EXPECT_GE(r.mNumEdgeRemovals, 1) << "k = " << k;
        EXPECT_EQ(r.mNumTets, 2 * (k - 2)) << "k = " << k;
        EXPECT_GT(r.mMinDihedralAfter, best[k] - 0.01) << "k = " << k;
        EXPECT_NEAR(r.mMinDihedralAfter, least_dihedral(r.mMesh), 1e-9);
        EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
        EXPECT_EQ(boundary_facets(in), boundary_facets(r.mMesh)) << "k = " << k;
        // A ring larger than mMaxRingSize is left alone.
        o.mMaxRingSize = 3;
        const OptimizeVolumeResult small = optimize_volume(in, o);
        EXPECT_EQ(small.mNumEdgeRemovals, 0) << "k = " << k;
    }
}

// The inverse: three faces (0,1,2), (0,2,3), (0,3,4) of a pentagon ring
// sandwiched between apexes 5 and 6. No 2-3 flip of one face and no edge
// removal of one diagonal raises the least angle of the tets it touches, but
// the five tets around the new edge (5, 6) do: 12.9 -> 28.1 degrees (found by
// a search over perturbed rings).
TEST(OptimizeVolume, MultiFaceRemovalInsertsTheEdgeBetweenTwoApexes) {
    const Mesh in = mt::make_mesh(
        {{0.38, 0.5, -0.08},
         {-1.46, 0.23, -0.08},
         {-0.96, -0.72, -0.03},
         {-0.15, -0.61, 0.09},
         {1.29, -0.55, -0.09},
         {0.25, 0.07, 0.27},
         {0.04, -0.02, -0.76}},
        "tetra",
        {{0, 1, 2, 5}, {0, 2, 1, 6}, {0, 2, 3, 5}, {0, 3, 2, 6}, {0, 3, 4, 5}, {0, 4, 3, 6}});
    OptimizeVolumeOptions o;
    o.mRelocate = false;
    const OptimizeVolumeResult r = optimize_volume(in, o);
    EXPECT_EQ(r.mNumMultiFaceRemovals, 1);
    EXPECT_EQ(r.mNum23Flips, 0);
    EXPECT_EQ(r.mNumEdgeRemovals, 0);
    EXPECT_EQ(r.mNumTets, 5);
    EXPECT_NEAR(r.mMinDihedralBefore, 12.95, 0.01);
    EXPECT_NEAR(r.mMinDihedralAfter, 28.13, 0.01);
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
    EXPECT_EQ(boundary_facets(in), boundary_facets(r.mMesh));
    // Capped below the disc's ring, it is left alone.
    o.mMaxRingSize = 4;
    EXPECT_EQ(optimize_volume(in, o).mNumMultiFaceRemovals, 0);
}

// The probe of roadmap 3.2.1 on a mesh the core can build without Netgen:
// the least dihedral angle never drops, here it rises, and the boundary holds.
TEST(OptimizeVolume, WorstObjectiveNeverLowersTheLeastDihedralAngle) {
    for (std::size_t n : {3u, 5u}) {
        const Mesh in = jittered_cube_fixture(n);
        const OptimizeVolumeResult r = optimize_volume(in);
        EXPECT_NEAR(r.mMinDihedralBefore, least_dihedral(in), 1e-9);
        EXPECT_NEAR(r.mMinDihedralAfter, least_dihedral(r.mMesh), 1e-9);
        EXPECT_GT(r.mMinDihedralAfter, r.mMinDihedralBefore) << "n = " << n;
        EXPECT_GT(r.mNumWorstSmoothed, 0);
        EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
        EXPECT_EQ(boundary_facets(in), boundary_facets(r.mMesh));
        // Each half alone keeps the guarantee too.
        for (const bool relocate : {true, false}) {
            OptimizeVolumeOptions half;
            half.mRelocate = relocate;
            half.mFlip = !relocate;
            const OptimizeVolumeResult h = optimize_volume(in, half);
            EXPECT_GE(h.mMinDihedralAfter, h.mMinDihedralBefore) << "relocate = " << relocate;
        }
    }
}

// Where the energy objective lowers the least angle, the worst one does not.
TEST(OptimizeVolume, EnergyMayLowerTheLeastAngleWorstDoesNot) {
    const Mesh in = jittered_cube_fixture(5);
    OptimizeVolumeOptions energy;
    energy.mObjective = OptimizeVolumeObjective::Energy;
    const OptimizeVolumeResult e = optimize_volume(in, energy);
    const OptimizeVolumeResult w = optimize_volume(in);
    EXPECT_GE(w.mMinDihedralAfter, w.mMinDihedralBefore);
    EXPECT_GE(w.mMinDihedralAfter, e.mMinDihedralAfter);
    EXPECT_EQ(e.mNumEdgeRemovals + e.mNumMultiFaceRemovals + e.mNumWorstSmoothed, 0);
}

TEST(OptimizeVolume, ObjectiveNamesAndBadRingSize) {
    EXPECT_EQ(meshioplusplus::optimize_volume_objective_from_name("worst"),
              OptimizeVolumeObjective::Worst);
    EXPECT_EQ(meshioplusplus::optimize_volume_objective_from_name("energy"),
              OptimizeVolumeObjective::Energy);
    EXPECT_THROW(meshioplusplus::optimize_volume_objective_from_name("best"),
                 std::invalid_argument);
    OptimizeVolumeOptions o;
    o.mMaxRingSize = 2;
    EXPECT_THROW(optimize_volume(two_tet_23_fixture(), o), std::invalid_argument);
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
// The worst objective's result, pinned like the energy one above.
TEST(OptimizeVolume, WorstResultMatchesTheGoldenDigest) {
    const OptimizeVolumeResult r = optimize_volume(jittered_cube_fixture(5));
    EXPECT_EQ(optimize_volume_digest(r), 0x3ebd8d80ee6ad35full);
}
#endif

TEST(OptimizeVolume, RejectsNonTetByName) {
    EXPECT_THROW(optimize_volume(mt::hex_mesh()), std::invalid_argument);
    EXPECT_THROW(optimize_volume(mt::tri_mesh()), std::invalid_argument);
    EXPECT_THROW(optimize_volume(mt::tet10_mesh()), std::invalid_argument);
}

}  // namespace
