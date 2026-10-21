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
// The OpenCASCADE-free half of fit_cad's freeform faces: flattening a patch,
// the fold test, and splitting a patch no flattening parametrises.

#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/patch_flatten.hpp"

using namespace meshioplusplus::detail;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Patch {
    std::vector<PatchVec> mPts;
    std::vector<PatchTri> mTris;
};

// A grid of (I + 1) x (J + 1) points through `rAt(u, v)`, u and v in [0, 1],
// two triangles per cell. With `Wrap`, u is periodic (the last column is the
// first): a closed band.
template <class F>
Patch grid(int I, int J, bool Wrap, F&& rAt) {
    Patch p;
    const int ni = Wrap ? I : I + 1;
    for (int j = 0; j <= J; ++j)
        for (int i = 0; i < ni; ++i)
            p.mPts.push_back(rAt(static_cast<double>(i) / I, static_cast<double>(j) / J));
    auto id = [&](int i, int j) { return static_cast<std::size_t>(j * ni + (i % ni)); };
    for (int j = 0; j < J; ++j)
        for (int i = 0; i < I; ++i) {
            p.mTris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            p.mTris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
    return p;
}

// A band of revolution whose radius waves along z: no primitive, and closed,
// so neither its plane nor an unrolled cylinder parametrises it.
Patch wavy_band(int N = 48, int M = 24) {
    return grid(N, M, true, [](double U, double V) {
        const double a = 2.0 * kPi * U, z = 2.0 * V;
        const double r = 1.0 + 0.15 * std::sin(2.0 * kPi * z);
        return PatchVec{r * std::cos(a), r * std::sin(a), z};
    });
}

bool flattenable(const Patch& rP, const std::vector<std::int64_t>& rPiece, std::int64_t K) {
    std::vector<std::size_t> sub;
    for (std::size_t t = 0; t < rPiece.size(); ++t)
        if (rPiece[t] == K)
            sub.push_back(t);
    std::set<std::size_t> verts;
    for (std::size_t t : sub)
        for (std::size_t v : rP.mTris[t])
            verts.insert(v);
    std::vector<std::size_t> local(rP.mPts.size(), 0);
    std::vector<PatchVec> pts;
    for (std::size_t v : verts) {
        local[v] = pts.size();
        pts.push_back(rP.mPts[v]);
    }
    std::vector<PatchTri> tris;
    for (std::size_t t : sub)
        tris.push_back({local[rP.mTris[t][0]], local[rP.mTris[t][1]], local[rP.mTris[t][2]]});
    return patch_flattenable(pts, tris);
}

}  // namespace

TEST(PatchFlatten, APlaneAndAnOpenCylinderBandFlatten) {
    const Patch plane =
        grid(8, 8, false, [](double U, double V) { return PatchVec{U, V, 0.1 * U * V}; });
    EXPECT_TRUE(patch_flattenable(plane.mPts, plane.mTris));
    // Three quarters of a tall cylinder: its best-fit plane is vertical and
    // folds it; the unrolled cylinder does not.
    const Patch arc = grid(24, 6, false, [](double U, double V) {
        const double a = 1.5 * kPi * U;
        return PatchVec{std::cos(a), std::sin(a), 4.0 * V};
    });
    std::vector<PatchVec> normals;
    std::vector<double> areas;
    patch_normals(arc.mPts, arc.mTris, normals, areas);
    PatchFlat flat;
    ASSERT_TRUE(patch_flatten(0, arc.mPts, normals, areas, flat));
    EXPECT_TRUE(patch_flat_folds(flat, arc.mTris));
    ASSERT_TRUE(patch_flatten(1, arc.mPts, normals, areas, flat));
    EXPECT_FALSE(patch_flat_folds(flat, arc.mTris));
    EXPECT_NEAR(flat.mRadius, 1.0, 1e-9);
    // The map back reproduces the points.
    for (std::size_t i = 0; i < arc.mPts.size(); ++i) {
        const PatchVec q = flat.At(flat.mS[i], flat.mT[i], flat.mH[i]);
        for (int a = 0; a < 3; ++a)
            EXPECT_NEAR(q[a], arc.mPts[i][a], 1e-9);
    }
}

TEST(PatchFlatten, AClosedWavyBandFoldsOverBoth) {
    const Patch band = wavy_band();
    EXPECT_FALSE(patch_flattenable(band.mPts, band.mTris));
}

TEST(PatchFlatten, SplittingLeavesPiecesThatFlatten) {
    const Patch band = wavy_band();
    std::size_t n = 0;
    const std::vector<std::int64_t> piece = patch_split_until_flat(band.mPts, band.mTris, {}, n);
    ASSERT_GE(n, 2u);
    std::size_t covered = 0;
    for (std::int64_t k : piece)
        covered += k >= 0;
    EXPECT_EQ(covered, piece.size());  // every triangle lands in a piece
    for (std::size_t k = 0; k < n; ++k)
        EXPECT_TRUE(flattenable(band, piece, static_cast<std::int64_t>(k))) << "piece " << k;
    // Pieces are numbered in the order of their lowest triangle.
    std::int64_t seen = -1;
    for (std::int64_t k : piece)
        if (k > seen) {
            EXPECT_EQ(k, seen + 1);
            seen = k;
        }
}

TEST(PatchFlatten, SplittingIsDeterministicAndRespectsItsLimits) {
    const Patch band = wavy_band();
    std::size_t a = 0, b = 0;
    EXPECT_EQ(patch_split_until_flat(band.mPts, band.mTris, {}, a),
              patch_split_until_flat(band.mPts, band.mTris, {}, b));
    EXPECT_EQ(a, b);
    // No halving allowed: the band stays facets.
    PatchSplitOptions none;
    none.mMaxDepth = 0;
    std::size_t n = 0;
    for (std::int64_t k : patch_split_until_flat(band.mPts, band.mTris, none, n))
        EXPECT_EQ(k, -1);
    EXPECT_EQ(n, 0u);
    // Pieces too small to keep stay facets too.
    PatchSplitOptions huge;
    huge.mMinTriangles = band.mTris.size() + 1;
    patch_split_until_flat(band.mPts, band.mTris, huge, n);
    EXPECT_EQ(n, 0u);
}

TEST(PatchFlatten, BisectingHalvesAPatchAcrossItsLongestExtent) {
    // A 40 x 10 sheet, long in x: the cut is across x, so each half is 20 x 10.
    const Patch sheet = grid(40, 10, false, [](double U, double V) {
        return PatchVec{40.0 * U, 10.0 * V, 0.2 * std::sin(6.0 * U)};
    });
    std::size_t n = 0;
    const std::vector<std::int64_t> piece = patch_bisect(sheet.mPts, sheet.mTris, 16, n);
    ASSERT_EQ(n, 2u);
    ASSERT_EQ(piece.size(), sheet.mTris.size());
    std::size_t count[2] = {0, 0};
    double max_x[2] = {-1e300, -1e300}, min_x[2] = {1e300, 1e300};
    for (std::size_t t = 0; t < piece.size(); ++t) {
        ASSERT_GE(piece[t], 0);
        ++count[piece[t]];
        for (std::size_t v : sheet.mTris[t]) {
            max_x[piece[t]] = std::max(max_x[piece[t]], sheet.mPts[v][0]);
            min_x[piece[t]] = std::min(min_x[piece[t]], sheet.mPts[v][0]);
        }
    }
    EXPECT_EQ(count[0] + count[1], piece.size());
    EXPECT_NEAR(static_cast<double>(count[0]), static_cast<double>(count[1]), 20.0);
    // Pieces are numbered by their lowest triangle, so piece 0 holds triangle 0 (x = 0).
    EXPECT_EQ(piece[0], 0);
    EXPECT_LE(max_x[0], 21.0);
    EXPECT_GE(min_x[1], 19.0);
}

TEST(PatchFlatten, BisectingIsDeterministicAndRespectsItsMinimum) {
    const Patch band = wavy_band();
    std::size_t a = 0, b = 0;
    EXPECT_EQ(patch_bisect(band.mPts, band.mTris, 16, a), patch_bisect(band.mPts, band.mTris, 16, b));
    EXPECT_EQ(a, b);
    EXPECT_GE(a, 2u);
    // A patch under twice the minimum is not halved.
    std::size_t n = 7;
    EXPECT_TRUE(patch_bisect(band.mPts, band.mTris, band.mTris.size(), n).empty());
    EXPECT_EQ(n, 0u);
}
