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
// Lines, arcs and circles from chained sections (roadmap 7.6.2).

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/fit_curves.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/slice.hpp"

using namespace meshioplusplus;

namespace {

constexpr double kFcTestPi = 3.14159265358979323846;

// The side walls of a prism over a closed outline, z in [0, 1], as triangles.
Mesh fc_prism(const std::vector<std::array<double, 2>>& rOutline) {
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> tris;
    const auto n = static_cast<std::int64_t>(rOutline.size());
    for (const auto& r_p : rOutline)
        pts.push_back({r_p[0], r_p[1], 0.0});
    for (const auto& r_p : rOutline)
        pts.push_back({r_p[0], r_p[1], 1.0});
    for (std::int64_t i = 0; i < n; ++i) {
        const std::int64_t j = (i + 1) % n;
        tris.push_back({i, j, j + n});
        tris.push_back({i, j + n, i + n});
    }
    return mt::make_mesh(pts, "triangle", tris);
}

// A slot: two semicircles of radius R about (+-A, 0), K segments each.
std::vector<std::array<double, 2>> fc_stadium(double A, double R, int K) {
    std::vector<std::array<double, 2>> out;
    for (int k = 0; k <= K; ++k) {
        const double t = -kFcTestPi / 2 + kFcTestPi * k / K;
        out.push_back({A + R * std::cos(t), R * std::sin(t)});
    }
    for (int k = 0; k <= K; ++k) {
        const double t = kFcTestPi / 2 + kFcTestPi * k / K;
        out.push_back({-A + R * std::cos(t), R * std::sin(t)});
    }
    return out;
}

Mesh fc_section(const Mesh& rMesh) {
    SliceOptions o;
    o.mOrigin = {0.0, 0.0, 0.5};
    o.mChain = true;
    return slice(rMesh, o);
}

}  // namespace

TEST(FitCurves, ASlotGivesTwoLinesAndTwoArcsWithExactRadii) {
    const Mesh slot = fc_prism(fc_stadium(2.0, 0.75, 24));
    FitCurvesOptions o;
    o.mTolerance = 1e-6;
    const FitCurvesResult r = fit_curves(fc_section(slot), o);
    EXPECT_EQ(r.mNumLines, 2);
    EXPECT_EQ(r.mNumArcs, 2);
    EXPECT_EQ(r.mNumCircles, 0);
    for (const SketchPiece& p : r.mPieces) {
        if (p.mType == SketchPieceType::Arc) {
            EXPECT_NEAR(p.mRadius, 0.75, 1e-6);
            EXPECT_NEAR(std::fabs(p.mSweepDeg), 180.0, 1e-4);
            EXPECT_NEAR(std::fabs(p.mCentre[0]), 2.0, 1e-6);
        } else {
            // The straight sides run from x = -2 to 2, corners on the tangents.
            EXPECT_NEAR(std::hypot(p.mEnd[0] - p.mStart[0], p.mEnd[1] - p.mStart[1]), 4.0, 1e-6);
        }
    }
    ASSERT_EQ(r.mFrames.size(), 1u);
    EXPECT_NEAR(std::fabs(r.mFrames[0][10]), 1.0, 1e-12);  // the normal is z
    EXPECT_TRUE(r.mMesh.HasFieldData(kSketchParamsName));
    EXPECT_TRUE(r.mMesh.HasCellData(kSketchPieceName));
    // Every segment belongs to a piece.
    const NDArray& piece = r.mMesh.CellData(kSketchPieceName, 0);
    for (std::size_t k = 0; k < piece.Size(); ++k)
        EXPECT_GE(piece.As<std::int64_t>()[k], 0);
}

TEST(FitCurves, APlateWithAHoleIsFourLinesAndACircle) {
    std::vector<std::array<double, 2>> square{{0, 0},   {1, 0},   {2, 0},   {3, 0},
                                              {3, 1},   {3, 2},   {3, 2.5}, {2, 2.5},
                                              {1, 2.5}, {0, 2.5}, {0, 1}};
    std::vector<std::array<double, 2>> hole;
    for (int k = 0; k < 32; ++k) {
        const double t = -2 * kFcTestPi * k / 32;
        hole.push_back({1.5 + 0.5 * std::cos(t), 1.25 + 0.5 * std::sin(t)});
    }
    const Mesh outer = fc_prism(square), inner = fc_prism(hole);
    const Mesh plate = merge({&outer, &inner}).mMesh;
    FitCurvesOptions o;
    o.mTolerance = 1e-6;
    const FitCurvesResult r = fit_curves(fc_section(plate), o);
    EXPECT_EQ(r.mNumLines, 4);
    EXPECT_EQ(r.mNumCircles, 1);
    EXPECT_EQ(r.mNumArcs, 0);
    for (const SketchPiece& p : r.mPieces)
        if (p.mType == SketchPieceType::Circle)
            EXPECT_NEAR(p.mRadius, 0.5, 1e-9);
}

TEST(FitCurves, SnappingAndRefusals) {
    // A rectangle turned by 0.3 degrees snaps back onto the axes.
    const double t = 0.3 * kFcTestPi / 180.0;
    std::vector<std::array<double, 2>> rect;
    for (const auto& r_p : std::vector<std::array<double, 2>>{
             {0, 0}, {2, 0}, {4, 0}, {4, 1}, {4, 2}, {2, 2}, {0, 2}, {0, 1}})
        rect.push_back({r_p[0] * std::cos(t) - r_p[1] * std::sin(t),
                        r_p[0] * std::sin(t) + r_p[1] * std::cos(t)});
    FitCurvesOptions o;
    o.mTolerance = 1e-6;
    o.mSnapAngleDeg = 1.0;
    o.mNormal = {0.0, 0.0, 1.0};
    const FitCurvesResult r = fit_curves(fc_section(fc_prism(rect)), o);
    ASSERT_EQ(r.mNumLines, 4);
    for (const SketchPiece& p : r.mPieces) {
        const double dx = p.mEnd[0] - p.mStart[0], dy = p.mEnd[1] - p.mStart[1];
        EXPECT_LT(std::min(std::fabs(dx), std::fabs(dy)), 1e-12);
    }
    EXPECT_THROW(fit_curves(mt::tri_mesh()), std::invalid_argument);
    EXPECT_THROW(fit_curves(mt::line_mesh()), std::invalid_argument);  // not chained
    o.mSnapAngleDeg = 30.0;
    EXPECT_THROW(fit_curves(fc_section(fc_prism(rect)), o), std::invalid_argument);
}

TEST(FitCurves, ARectangleIsFourLinesNotACircle) {
    // Its four corners lie on one circle; a sampled arc steps at most 30
    // degrees a point, so the outline is four lines.
    const Mesh box = fc_prism({{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}});
    FitCurvesOptions o;
    o.mTolerance = 1e-6;
    const FitCurvesResult r = fit_curves(fc_section(box), o);
    EXPECT_EQ(r.mNumLines, 4);
    EXPECT_EQ(r.mNumCircles, 0);
    EXPECT_EQ(r.mNumArcs, 0);
}
