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
// DXF R2000 drawings (roadmap 7.6.4).
// The outline of a part along a direction (roadmap 7.6.3).

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/fit_curves.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/silhouette.hpp"
#include "meshioplusplus/operations/slice.hpp"

using namespace meshioplusplus;

namespace {

constexpr double kSilTestPi = 3.14159265358979323846;

// Two unit squares overlapping in a quarter: crossing edges, one outline.
Mesh sil_two_squares() {
    return mt::make_mesh({{0, 0, 0},
                          {1, 0, 0},
                          {1, 1, 0},
                          {0, 1, 0},
                          {0.5, 0.5, 0.2},
                          {1.5, 0.5, 0.2},
                          {1.5, 1.5, 0.2},
                          {0.5, 1.5, 0.2}},
                         "quad", {{0, 1, 2, 3}, {4, 5, 6, 7}});
}

// Each edge of a triangle surface is shared by exactly two triangles, once
// each way: closed and consistently wound.
bool sil_closed(const Mesh& rMesh) {
    std::map<std::array<std::int64_t, 2>, int> uses;
    const auto& c = rMesh.Cells(0);
    const std::int64_t* p = c.Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < c.NumCells(); ++t)
        for (std::size_t k = 0; k < 3; ++k)
            ++uses[{p[3 * t + k], p[3 * t + (k + 1) % 3]}];
    for (const auto& [e, n] : uses)
        if (n != 1 || !uses.count({e[1], e[0]}))
            return false;
    return true;
}

}  // namespace

TEST(Silhouette, ASphereIsOneCircleOfItsRadius) {
    for (const std::array<double, 3> n :
         {std::array<double, 3>{0, 0, 1}, {1, 0, 0}, {1, 1, 1}, {0.3, -0.2, 0.9}}) {
        SilhouetteOptions o;
        o.mNormal = n;
        const SilhouetteResult r = silhouette(icosa_sphere(2.0, 3), o);
        EXPECT_EQ(r.mNumLoops, 1);
        EXPECT_EQ(r.mNumHoles, 0);
        const double* p = r.mMesh.Points().As<double>();
        const double nn = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
            const double h = (p[3 * i] * n[0] + p[3 * i + 1] * n[1] + p[3 * i + 2] * n[2]) / nn;
            const double rr = std::sqrt(p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 1] +
                                        p[3 * i + 2] * p[3 * i + 2] - h * h);
            EXPECT_LE(rr, 2.0 + 1e-12);
            EXPECT_GE(rr, 2.0 * 0.98);                 // within the tessellation's sag
            EXPECT_NEAR(h, -2.0 * 0.99, 2.0 * 0.011);  // on the plane at the lowest point
        }
        EXPECT_NEAR(r.mArea, kSilTestPi * 4.0, 0.15);
    }
}

TEST(Silhouette, APlateWithAHoleHasItsHole) {
    const SilhouetteResult r = silhouette(annulus(0.5, 1.0, 24, 3));
    EXPECT_EQ(r.mNumOuter, 1);
    EXPECT_EQ(r.mNumHoles, 1);
    EXPECT_TRUE(r.mMesh.HasCellData(kSilhouetteHoleName));
    const SilhouetteResult s = silhouette(sil_two_squares());
    EXPECT_EQ(s.mNumLoops, 1);
    EXPECT_NEAR(s.mArea, 1.75, 1e-12);
    EXPECT_NEAR(s.mPerimeter, 6.0, 1e-12);
}

TEST(Silhouette, ClearanceOffsetsExactly) {
    const Mesh cube = box({{2, 2, 2}});  // a unit cube of hexahedra: its skin
    SilhouetteOptions o;
    o.mClearance = 0.1;
    o.mRoundAngleDeg = 80.0;  // the square's corners turn by 90: mitred
    EXPECT_NEAR(silhouette(cube, o).mArea, 1.44, 1e-12);
    o.mRoundAngleDeg = 180.0;  // rounded, the arcs' sides tangent to the circle
    const double rounded = silhouette(cube, o).mArea;
    EXPECT_GT(rounded, 1.0 + 0.4 + kSilTestPi * 0.01);
    EXPECT_LT(rounded, 1.0 + 0.4 + kSilTestPi * 0.01 * 1.01);
    o.mClearance = -0.1;  // a press fit: the convex corners stay sharp
    EXPECT_NEAR(silhouette(cube, o).mArea, 0.64, 1e-12);
    o.mClearance = -0.6;
    const SilhouetteResult gone = silhouette(cube, o);
    EXPECT_EQ(gone.mNumLoops, 0);
    EXPECT_FALSE(gone.mWarnings.empty());
}

TEST(Silhouette, TheOutlineFeedsFitCurves) {
    SilhouetteOptions o;
    o.mNormal = {0, 0, 1};
    const SilhouetteResult r =
        silhouette(annulus(0.5, 1.0, 64, 2, 2 * kSilTestPi, {0, 0, 0}, "triangle"), o);
    FitCurvesOptions fo;
    fo.mTolerance = 1e-3;
    const FitCurvesResult f = fit_curves(r.mMesh, fo);
    EXPECT_EQ(f.mNumCircles, 2);
}

TEST(Silhouette, ExtrudeMakesAClosedPrism) {
    SilhouetteOptions o;
    o.mExtrude = true;
    o.mNormal = {0, 0, 1};
    const SilhouetteResult r = silhouette(icosa_sphere(1.0, 2), o);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_TRUE(sil_closed(r.mMesh));
    EXPECT_NEAR(r.mMinHeight, -1.0, 1e-12);
    EXPECT_NEAR(r.mMaxHeight, 1.0, 1e-12);
}

TEST(Silhouette, IsDeterministicAndRefusesBadInput) {
    SilhouetteOptions o;
    o.mNormal = {0.2, 0.5, 1.0};
    o.mClearance = 0.05;
    const SilhouetteResult a = silhouette(icosa_sphere(1.0, 2), o);
    const SilhouetteResult b = silhouette(icosa_sphere(1.0, 2), o);
    ASSERT_EQ(a.mMesh.NumPoints(), b.mMesh.NumPoints());
    const double* pa = a.mMesh.Points().As<double>();
    const double* pb = b.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < 3 * a.mMesh.NumPoints(); ++i)
        EXPECT_EQ(pa[i], pb[i]);
    o.mNormal = {0, 0, 0};
    EXPECT_THROW(silhouette(icosa_sphere(), o), std::invalid_argument);
    EXPECT_THROW(silhouette(mt::line_mesh()), std::invalid_argument);
    SilhouetteOptions bad;
    bad.mArcSegments = 2;
    EXPECT_THROW(silhouette(icosa_sphere(), bad), std::invalid_argument);
}
