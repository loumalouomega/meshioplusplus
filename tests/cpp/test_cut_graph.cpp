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
// detail/cut_graph and the parametrisation ladder (fit_cad's freeform modes,
// v17.12.0): topology, cuts into discs, injectivity, distortion, determinism.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

#include "analytic_surfaces.hpp"

#include "../../src/cpp/src/detail/cut_graph.hpp"
#include "../../src/cpp/src/detail/param_disc.hpp"

using namespace meshioplusplus;
using detail::PatchTri;
using detail::PatchVec;

namespace {

struct Patch {
    std::vector<PatchVec> mPts;
    std::vector<PatchTri> mTris;
};

// A triangle mesh as a patch, keeping the triangles `rKeep` accepts and only
// the points they use (renumbered in order of first use).
template <class Keep>
Patch to_patch(const Mesh& rMesh, Keep rKeep) {
    std::vector<PatchVec> all;
    const NDArray& x = rMesh.Points();
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
        all.push_back(
            {x.As<double>()[i * 3], x.As<double>()[i * 3 + 1], x.As<double>()[i * 3 + 2]});
    Patch p;
    std::vector<std::size_t> local(all.size(), SIZE_MAX);
    const NDArray& c = rMesh.Cells(0).Conn();
    for (std::size_t t = 0; t < rMesh.Cells(0).NumCells(); ++t) {
        const PatchTri tri{static_cast<std::size_t>(c.As<std::int64_t>()[t * 3]),
                           static_cast<std::size_t>(c.As<std::int64_t>()[t * 3 + 1]),
                           static_cast<std::size_t>(c.As<std::int64_t>()[t * 3 + 2])};
        PatchVec m{0.0, 0.0, 0.0};
        for (std::size_t v : tri)
            for (int j = 0; j < 3; ++j)
                m[static_cast<std::size_t>(j)] += all[v][static_cast<std::size_t>(j)] / 3.0;
        if (!rKeep(m))
            continue;
        PatchTri out{};
        for (int k = 0; k < 3; ++k) {
            std::size_t& l = local[tri[static_cast<std::size_t>(k)]];
            if (l == SIZE_MAX) {
                l = p.mPts.size();
                p.mPts.push_back(all[tri[static_cast<std::size_t>(k)]]);
            }
            out[static_cast<std::size_t>(k)] = l;
        }
        p.mTris.push_back(out);
    }
    return p;
}

Patch all_of(const Mesh& rMesh) {
    return to_patch(rMesh, [](const PatchVec&) { return true; });
}

// The unit cube's four sides (no top, no bottom): an annulus.
Patch cube_band(int n) {
    return to_patch(analytic::triangle_box(n),
                    [](const PatchVec& rM) { return rM[2] > 1e-9 && rM[2] < 1.0 - 1e-9; });
}

std::vector<double> sizes_of(const std::vector<std::int64_t>& rLabel, std::size_t N) {
    std::vector<double> s(N, 0.0);
    for (std::int64_t l : rLabel)
        s[static_cast<std::size_t>(l)] += 1.0;
    return s;
}

}  // namespace

TEST(CutGraph, TopologyOfDiscAnnulusSphereAndTorus) {
    const Patch top =
        to_patch(analytic::triangle_box(4), [](const PatchVec& rM) { return rM[2] > 1.0 - 1e-9; });
    const auto disc = detail::patch_topology(top.mTris);
    EXPECT_TRUE(disc.IsDisc());
    EXPECT_EQ(disc.mGenus, 0);

    const auto band = detail::patch_topology(cube_band(4).mTris);
    EXPECT_EQ(band.mNumBoundaryLoops, 2u);
    EXPECT_EQ(band.mEuler, 0);
    EXPECT_EQ(band.mGenus, 0);
    EXPECT_FALSE(band.IsDisc());

    const auto sphere = detail::patch_topology(all_of(analytic::cc_cage(3, 2)).mTris);
    EXPECT_EQ(sphere.mNumBoundaryLoops, 0u);
    EXPECT_EQ(sphere.mEuler, 2);
    EXPECT_EQ(sphere.mGenus, 0);

    const auto torus = detail::patch_topology(all_of(analytic::organic_torus(24, 12)).mTris);
    EXPECT_EQ(torus.mEuler, 0);
    EXPECT_EQ(torus.mGenus, 1);

    const auto two = detail::patch_topology(all_of(analytic::double_torus(24)).mTris);
    EXPECT_TRUE(two.mManifold);
    EXPECT_EQ(two.mGenus, 2);
}

TEST(CutGraph, AnnulusSplitsIntoTwoDiscs) {
    const Patch band = cube_band(6);
    std::size_t n = 0;
    const std::vector<std::int64_t> label =
        detail::patch_split_to_discs(band.mPts, band.mTris, 8, n);
    ASSERT_EQ(n, 2u);
    ASSERT_EQ(label.size(), band.mTris.size());
    const std::vector<double> s = sizes_of(label, n);
    // Two cuts half way round each other: comparable halves.
    EXPECT_GT(std::min(s[0], s[1]) / std::max(s[0], s[1]), 0.3);
    // Deterministic.
    std::size_t n2 = 0;
    EXPECT_EQ(detail::patch_split_to_discs(band.mPts, band.mTris, 8, n2), label);
}

TEST(CutGraph, ThreeHolesSplitIntoTwoDiscs) {
    // The four sides of the cube with a hole in the front one.
    const Patch p = to_patch(analytic::triangle_box(6), [](const PatchVec& rM) {
        if (rM[2] < 1e-9 || rM[2] > 1.0 - 1e-9)
            return false;
        const bool hole =
            rM[1] < 1e-9 && rM[0] > 0.34 && rM[0] < 0.66 && rM[2] > 0.34 && rM[2] < 0.66;
        return !hole;
    });
    ASSERT_EQ(detail::patch_topology(p.mTris).mNumBoundaryLoops, 3u);
    std::size_t n = 0;
    const std::vector<std::int64_t> label = detail::patch_split_to_discs(p.mPts, p.mTris, 4, n);
    EXPECT_EQ(n, 2u);
    EXPECT_EQ(label.size(), p.mTris.size());
}

TEST(CutGraph, ClosedSphereSplitsIntoTwoCaps) {
    const Patch p = all_of(analytic::cc_cage(3, 3));
    std::size_t n = 0;
    const std::vector<std::int64_t> label = detail::patch_split_to_discs(p.mPts, p.mTris, 16, n);
    ASSERT_EQ(n, 2u);
    const std::vector<double> s = sizes_of(label, n);
    EXPECT_GT(std::min(s[0], s[1]) / std::max(s[0], s[1]), 0.5);
}

TEST(CutGraph, DiscAndTorusAreNotCut) {
    const Patch top =
        to_patch(analytic::triangle_box(4), [](const PatchVec& rM) { return rM[2] > 1.0 - 1e-9; });
    std::size_t n = 7;
    EXPECT_TRUE(detail::patch_split_to_discs(top.mPts, top.mTris, 4, n).empty());
    EXPECT_EQ(n, 0u);
    const Patch torus = all_of(analytic::organic_torus(24, 12));
    EXPECT_TRUE(detail::patch_split_to_discs(torus.mPts, torus.mTris, 4, n).empty());
    EXPECT_EQ(n, 0u);
}

TEST(CutGraph, ShortestPathIsShortest) {
    const Patch top =
        to_patch(analytic::triangle_box(4), [](const PatchVec& rM) { return rM[2] > 1.0 - 1e-9; });
    // From a corner to the opposite corner of the top: along the a-c diagonals
    // where they run that way, else a staircase; never longer than two sides.
    std::size_t c0 = 0, c1 = 0;
    for (std::size_t i = 0; i < top.mPts.size(); ++i) {
        if (top.mPts[i] == PatchVec{0.0, 0.0, 1.0})
            c0 = i;
        if (top.mPts[i] == PatchVec{1.0, 1.0, 1.0})
            c1 = i;
    }
    const std::vector<std::size_t> path =
        detail::patch_shortest_path(top.mPts, top.mTris, {c0}, {c1}, {});
    ASSERT_GE(path.size(), 2u);
    EXPECT_EQ(path.front(), c0);
    EXPECT_EQ(path.back(), c1);
    double len = 0.0;
    for (std::size_t k = 0; k + 1 < path.size(); ++k) {
        const PatchVec& a = top.mPts[path[k]];
        const PatchVec& b = top.mPts[path[k + 1]];
        len += std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                         (a[2] - b[2]) * (a[2] - b[2]));
    }
    EXPECT_LE(len, 2.0 + 1e-12);
    EXPECT_GE(len, std::sqrt(2.0) - 1e-12);
}

// Every rung the ladder takes is injective (no flipped triangle).
TEST(ParamLadder, InjectiveOnTheFixtures) {
    for (const Mesh& m : {analytic::monkey_saddle_box(12), analytic::gaussian_bump_box(12),
                          analytic::s_band_box(12, true), analytic::hood_box(12)}) {
        // The top: the triangles on none of the box's five planes.
        const Patch p = to_patch(m, [](const PatchVec& rC) {
            return rC[0] > 1e-12 && rC[0] < 1.0 - 1e-12 && rC[1] > 1e-12 && rC[1] < 1.0 - 1e-12 &&
                   rC[2] > 0.5;
        });
        for (detail::ParamLadder which :
             {detail::ParamLadder::Auto, detail::ParamLadder::Projection, detail::ParamLadder::Lscm,
              detail::ParamLadder::Harmonic}) {
            const detail::DiscParam par = detail::param_ladder(p.mPts, p.mTris, which);
            ASSERT_TRUE(par.mOk) << static_cast<int>(which);
            EXPECT_EQ(detail::param_flips(p.mTris, par.mU, par.mV), 0u);
            EXPECT_EQ(par.mMethod != 0, true);
        }
    }
    // A cap that folds over its plane takes a later rung.
    const Patch cap =
        to_patch(analytic::cc_cage(3, 3), [](const PatchVec& rC) { return rC[2] > -0.2; });
    const detail::DiscParam par =
        detail::param_ladder(cap.mPts, cap.mTris, detail::ParamLadder::Auto);
    ASSERT_TRUE(par.mOk);
    EXPECT_EQ(detail::param_flips(cap.mTris, par.mU, par.mV), 0u);
}

// Least-squares conformal maps keep angles on a smooth cap: the mean ratio
// of each triangle's Jacobian singular values stays near 1.
TEST(ParamLadder, LscmConformalDistortionIsBounded) {
    const Patch cap =
        to_patch(analytic::cc_cage(3, 3), [](const PatchVec& rC) { return rC[2] > 0.0; });
    const detail::DiscParam par =
        detail::param_ladder(cap.mPts, cap.mTris, detail::ParamLadder::Lscm);
    ASSERT_TRUE(par.mOk);
    double sum = 0.0;
    for (const PatchTri& t : cap.mTris) {
        // The 3-D triangle in its own frame, and its (u, v) image.
        const PatchVec& a = cap.mPts[t[0]];
        const PatchVec& b = cap.mPts[t[1]];
        const PatchVec& c = cap.mPts[t[2]];
        double e1[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        double e2[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double l1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
        const double x2 = (e2[0] * e1[0] + e2[1] * e1[1] + e2[2] * e1[2]) / l1;
        const double y2 =
            std::sqrt(std::max(0.0, e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2] - x2 * x2));
        const double du1 = par.mU[t[1]] - par.mU[t[0]], dv1 = par.mV[t[1]] - par.mV[t[0]];
        const double du2 = par.mU[t[2]] - par.mU[t[0]], dv2 = par.mV[t[2]] - par.mV[t[0]];
        // J maps (l1, 0) -> (du1, dv1) and (x2, y2) -> (du2, dv2).
        const double j00 = du1 / l1, j10 = dv1 / l1;
        const double j01 = (du2 - j00 * x2) / y2, j11 = (dv2 - j10 * x2) / y2;
        const double f = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
        const double det = std::fabs(j00 * j11 - j01 * j10);
        const double disc = std::sqrt(std::max(0.0, f * f - 4.0 * det * det));
        const double s1 = std::sqrt(0.5 * (f + disc)),
                     s2 = std::sqrt(std::max(1e-300, 0.5 * (f - disc)));
        sum += s1 / s2;
    }
    EXPECT_LT(sum / static_cast<double>(cap.mTris.size()), 1.25);
}

// The same parameters at any thread count.
TEST(ParamLadder, DeterministicAcrossThreadCounts) {
    const Patch cap =
        to_patch(analytic::cc_cage(5, 3), [](const PatchVec& rC) { return rC[2] > 0.1; });
    std::vector<double> ref_u, ref_v;
    for (int threads : {1, 2, 8}) {
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
        omp_set_num_threads(threads);
#endif
        for (detail::ParamLadder which :
             {detail::ParamLadder::Lscm, detail::ParamLadder::Harmonic}) {
            const detail::DiscParam par = detail::param_ladder(cap.mPts, cap.mTris, which);
            ASSERT_TRUE(par.mOk);
            if (which == detail::ParamLadder::Lscm) {
                if (ref_u.empty()) {
                    ref_u = par.mU;
                    ref_v = par.mV;
                } else {
                    EXPECT_EQ(par.mU, ref_u) << threads;
                    EXPECT_EQ(par.mV, ref_v) << threads;
                }
            }
        }
    }
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(omp_get_num_procs());
#endif
}
