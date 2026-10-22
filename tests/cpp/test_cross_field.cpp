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
// detail/cross_field (fit_cad's network mode, v17.12.0): Poincare-Hopf on
// closed skins of genus 0, 1 and 2, boundary alignment, determinism.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

#include "analytic_surfaces.hpp"

#include "../../src/cpp/src/detail/cross_field.hpp"

using namespace meshioplusplus;
using namespace analytic;
using detail::PatchTri;
using detail::PatchVec;

namespace {

struct CfPatch {
    std::vector<PatchVec> mPts;
    std::vector<PatchTri> mTris;
};

CfPatch cf_patch(const Mesh& rMesh) {
    CfPatch p;
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

int index_sum(const detail::CrossField& rF) {
    int s = 0;
    for (int i : rF.mIndex)
        s += i;
    return s;
}

}  // namespace

// Poincare-Hopf: the indices of a closed skin sum to its Euler characteristic
// (in quarters, 4 chi), whatever the field.
TEST(CrossField, IndicesSumToTheEulerCharacteristic) {
    struct Case {
        const char* mName;
        Mesh mMesh;
        int mChi;
    };
    std::vector<Case> cases;
    cases.push_back({"cc_cage3", cc_cage(3, 3), 2});
    cases.push_back({"cc_cage5", cc_cage(5, 3), 2});
    cases.push_back({"organic_torus", organic_torus(32, 16), 0});
    cases.push_back({"sphere_with_handle", sphere_with_handle(24), 0});
    cases.push_back({"double_torus", double_torus(32), -2});
    for (const Case& c : cases) {
        SCOPED_TRACE(c.mName);
        const CfPatch p = cf_patch(c.mMesh);
        const detail::CrossField f = detail::cross_field(p.mPts, p.mTris, {});
        ASSERT_TRUE(f.mOk);
        EXPECT_EQ(index_sum(f), 4 * c.mChi);
    }
}

// The cube's limit surface keeps eight valence-3 singularities, one per
// corner, and nothing else.
TEST(CrossField, CubeCageHasEightValenceThreeSingularities) {
    const CfPatch p = cf_patch(cc_cage(3, 3));
    const detail::CrossField f = detail::cross_field(p.mPts, p.mTris, {});
    ASSERT_TRUE(f.mOk);
    ASSERT_EQ(f.mSingular.size(), 8U);
    for (std::size_t v : f.mSingular) {
        EXPECT_EQ(f.mIndex[v], 1);
        EXPECT_EQ(detail::cross_valence(f, v), 3);
    }
}

// The sharp cube with its edges as features: the field follows the faces'
// grids, so the only singularities are the eight corners.
TEST(CrossField, FeatureEdgesAlignTheField) {
    const Mesh box = triangle_box(6);
    const CfPatch p = cf_patch(box);
    detail::CrossFieldOptions o;
    // The cube's edges: both ends on two of the planes x, y, z in {0, 1}.
    const auto on = [&](std::size_t v) {
        int k = 0;
        for (int j = 0; j < 3; ++j) {
            const double c = p.mPts[v][static_cast<std::size_t>(j)];
            k += (std::abs(c) < 1e-12 || std::abs(c - 1.0) < 1e-12) ? 1 : 0;
        }
        return k;
    };
    for (const PatchTri& t : p.mTris)
        for (int k = 0; k < 3; ++k) {
            const std::size_t a = t[static_cast<std::size_t>(k)];
            const std::size_t b = t[static_cast<std::size_t>((k + 1) % 3)];
            if (a < b && on(a) >= 2 && on(b) >= 2) {
                // Both on an edge line, and the same one: they share two planes.
                int shared = 0;
                for (int j = 0; j < 3; ++j) {
                    const double ca = p.mPts[a][static_cast<std::size_t>(j)];
                    const double cb = p.mPts[b][static_cast<std::size_t>(j)];
                    const bool pa = std::abs(ca) < 1e-12 || std::abs(ca - 1.0) < 1e-12;
                    shared += (pa && ca == cb) ? 1 : 0;
                }
                if (shared >= 2)
                    o.mFeatureEdges.push_back({a, b});
            }
        }
    ASSERT_EQ(o.mFeatureEdges.size(), 12U * 6U);
    const detail::CrossField f = detail::cross_field(p.mPts, p.mTris, o);
    ASSERT_TRUE(f.mOk);
    EXPECT_EQ(f.mSingular.size(), 8U);
    EXPECT_EQ(index_sum(f), 8);
    // Every triangle's cross holds the face's grid directions.
    for (std::size_t t = 0; t < p.mTris.size(); ++t) {
        const PatchVec d = detail::cross_branch(f, t, 0);
        const double m = std::max({std::abs(d[0]), std::abs(d[1]), std::abs(d[2])});
        EXPECT_NEAR(m, 1.0, 1e-6) << t;
    }
}

// A disc patch with its boundary constrained: a square's field is the grid,
// with no interior singularity.
TEST(CrossField, ConstrainedDiscHasNoInteriorSingularity) {
    const Mesh top = monkey_saddle_box(12);
    CfPatch all = cf_patch(top);
    // The saddle's top: the triangles facing up (the sides face sideways).
    CfPatch p;
    std::vector<std::size_t> local(all.mPts.size(), SIZE_MAX);
    for (const PatchTri& t : all.mTris) {
        const PatchVec& a = all.mPts[t[0]];
        const PatchVec& b = all.mPts[t[1]];
        const PatchVec& c = all.mPts[t[2]];
        const double nx = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]);
        const double ny = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]);
        const double nz = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if (!(nz > 0.3 * std::sqrt(nx * nx + ny * ny + nz * nz)))
            continue;
        PatchTri o{};
        for (int k = 0; k < 3; ++k) {
            std::size_t& l = local[t[static_cast<std::size_t>(k)]];
            if (l == SIZE_MAX) {
                l = p.mPts.size();
                p.mPts.push_back(all.mPts[t[static_cast<std::size_t>(k)]]);
            }
            o[static_cast<std::size_t>(k)] = l;
        }
        p.mTris.push_back(o);
    }
    ASSERT_GT(p.mTris.size(), 100U);
    const detail::CrossField f = detail::cross_field(p.mPts, p.mTris, {});
    ASSERT_TRUE(f.mOk);
    EXPECT_TRUE(f.mSingular.empty());
}

TEST(CrossField, DeterministicAcrossThreads) {
    const CfPatch p = cf_patch(double_torus(24));
    std::vector<double> ref;
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 8}) {
        omp_set_num_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const detail::CrossField f = detail::cross_field(p.mPts, p.mTris, {});
        if (ref.empty())
            ref = f.mTheta;
        else
            EXPECT_EQ(f.mTheta, ref);
    }
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
#endif
}
