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
// detail/quad_layout (fit_cad's network mode, v17.12.0): valid four-sided
// layouts on discs, on closed skins of genus 0, 1 and 2, along features, and
// deterministic.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

#include "analytic_surfaces.hpp"

#include "../../src/cpp/src/detail/quad_layout.hpp"

using namespace meshioplusplus;
using namespace analytic;
using detail::PatchTri;
using detail::PatchVec;

namespace {

struct QlPatch {
    std::vector<PatchVec> mPts;
    std::vector<PatchTri> mTris;
};

QlPatch ql_patch(const Mesh& rMesh) {
    QlPatch p;
    const double* x = rMesh.Points().As<double>();
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
        p.mPts.push_back({x[i * 3], x[i * 3 + 1], x[i * 3 + 2]});
    const std::int64_t* c = rMesh.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < rMesh.Cells(0).NumCells(); ++t)
        p.mTris.push_back({static_cast<std::size_t>(c[t * 3]),
                           static_cast<std::size_t>(c[t * 3 + 1]),
                           static_cast<std::size_t>(c[t * 3 + 2])});
    return p;
}

// The triangles of a skin facing up (a heightfield box's top).
QlPatch ql_top(const QlPatch& rAll) {
    QlPatch p;
    p.mPts = rAll.mPts;
    for (const PatchTri& t : rAll.mTris) {
        const PatchVec& a = rAll.mPts[t[0]];
        const PatchVec& b = rAll.mPts[t[1]];
        const PatchVec& c = rAll.mPts[t[2]];
        const double nx = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]);
        const double ny = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]);
        const double nz = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if (nz > 0.3 * std::sqrt(nx * nx + ny * ny + nz * nz))
            p.mTris.push_back(t);
    }
    return p;
}

double ql_diag(const QlPatch& rP) {
    PatchVec lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const PatchTri& t : rP.mTris)
        for (std::size_t v : t)
            for (std::size_t j = 0; j < 3; ++j) {
                lo[j] = std::min(lo[j], rP.mPts[v][j]);
                hi[j] = std::max(hi[j], rP.mPts[v][j]);
            }
    return std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                     (hi[2] - lo[2]) * (hi[2] - lo[2]));
}

void expect_valid(const detail::QuadLayout& rL, const QlPatch& rP) {
    ASSERT_TRUE(rL.mOk) << "failure " << static_cast<int>(rL.mFailure);
    EXPECT_EQ(detail::quad_layout_check(rL, rP.mTris), detail::QuadLayoutFailure::None);
    EXPECT_FALSE(rL.mPatches.empty());
    for (std::int64_t p : rL.mPatchOfTri)
        EXPECT_GE(p, 0);
    for (int v : rL.mValence) {
        EXPECT_GE(v, 2);
        EXPECT_NE(v, 4);
    }
}

}  // namespace

TEST(QuadLayout, SaddleTopIsOneQuad) {
    const QlPatch top = ql_top(ql_patch(monkey_saddle_box(12)));
    const detail::QuadLayout l = detail::quad_layout(top.mPts, top.mTris, {});
    expect_valid(l, top);
    EXPECT_EQ(l.mPatches.size(), 1U);
    EXPECT_TRUE(l.mExtraordinary.empty());
    // The square's four corners are the top's.
    for (std::size_t j = 0; j < 4; ++j) {
        const PatchVec& c = top.mPts[l.mPatches[0].mBoundary[l.mPatches[0].mCornerPos[j]]];
        EXPECT_TRUE((c[0] == 0.0 || c[0] == 1.0) && (c[1] == 0.0 || c[1] == 1.0));
    }
}

TEST(QuadLayout, SaddleTopSplitsToTheTargetSize) {
    const QlPatch top = ql_top(ql_patch(monkey_saddle_box(24)));
    detail::QuadLayoutOptions o;
    o.mTargetSize = 0.3;
    const detail::QuadLayout l = detail::quad_layout(top.mPts, top.mTris, o);
    expect_valid(l, top);
    EXPECT_GE(l.mPatches.size(), 4U);
    EXPECT_GT(l.mNumInteriorArcs, 0U);
}

TEST(QuadLayout, FeaturesMakeTheCubeSixQuads) {
    const QlPatch box = ql_patch(triangle_box(6));
    const detail::QuadLayout l = detail::quad_layout(box.mPts, box.mTris, {});
    expect_valid(l, box);
    EXPECT_EQ(l.mPatches.size(), 6U);
    EXPECT_EQ(l.mExtraordinary.size(), 8U);
    for (int v : l.mValence)
        EXPECT_EQ(v, 3);
    EXPECT_EQ(l.mChains.size(), 12U);
    EXPECT_EQ(l.mNumTJunctions, 0U);
}

TEST(QuadLayout, ClosedSkinsOfEveryGenus) {
    struct Case {
        const char* mName;
        Mesh mMesh;
    };
    std::vector<Case> cases;
    cases.push_back({"cc_cage3", cc_cage(3, 3)});
    cases.push_back({"cc_cage5", cc_cage(5, 3)});
    cases.push_back({"organic_torus", organic_torus(32, 16)});
    cases.push_back({"sphere_with_handle", sphere_with_handle(24)});
    cases.push_back({"double_torus", double_torus(32)});
    for (const Case& c : cases) {
        SCOPED_TRACE(c.mName);
        const QlPatch p = ql_patch(c.mMesh);
        detail::QuadLayoutOptions o;
        o.mTargetSize = 0.25 * ql_diag(p);
        const detail::QuadLayout l = detail::quad_layout(p.mPts, p.mTris, o);
        expect_valid(l, p);
        std::size_t multi = 0;
        for (const auto& sides : l.mSides)
            for (const auto& sd : sides)
                multi += sd.mSegments.size() > 1 ? 1 : 0;
        std::printf(
            "  %-20s %zu patches, %zu chains, %zu EVs, %zu T, %zu separatrix splits, "
            "%zu multi-segment sides\n",
            c.mName, l.mPatches.size(), l.mChains.size(), l.mExtraordinary.size(), l.mNumTJunctions,
            l.mNumSeparatrixSplits, multi);
    }
}

TEST(QuadLayout, HoodAndBumpTops) {
    for (const Mesh& m : {hood_box(16, true), gaussian_bump_box(16), s_band_box(16)}) {
        const QlPatch top = ql_top(ql_patch(m));
        detail::QuadLayoutOptions o;
        o.mTargetSize = 0.4;
        const detail::QuadLayout l = detail::quad_layout(top.mPts, top.mTris, o);
        expect_valid(l, top);
    }
}

TEST(QuadLayout, DeterministicAcrossThreads) {
    const QlPatch p = ql_patch(double_torus(24));
    detail::QuadLayoutOptions o;
    o.mTargetSize = 0.25 * ql_diag(p);
    std::vector<std::int64_t> ref;
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 8}) {
        omp_set_num_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const detail::QuadLayout l = detail::quad_layout(p.mPts, p.mTris, o);
        if (ref.empty())
            ref = l.mPatchOfTri;
        else
            EXPECT_EQ(l.mPatchOfTri, ref);
    }
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
#endif
}
