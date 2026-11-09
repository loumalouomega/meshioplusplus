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
// detail/dense_lp (roadmap 7.1.4) and the minimum-zone fits on it.

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/dense_lp.hpp"
#include "../../src/cpp/src/detail/zone_fit.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/segment.hpp"

using namespace meshioplusplus;
using detail::DenseLpStatus;

TEST(DenseLp, TextbookProblems) {
    // max 3x + 5y s.t. x <= 4, 2y <= 12, 3x + 2y <= 18 (x, y free; the
    // optimum is (2, 6), 36).
    const auto r = detail::dense_lp_inequality({-3.0, -5.0}, {1, 0, 0, 2, 3, 2, -1, 0, 0, -1},
                                               {4, 12, 18, 0, 0});
    ASSERT_EQ(r.mStatus, DenseLpStatus::Optimal);
    EXPECT_NEAR(r.mX[0], 2.0, 1e-12);
    EXPECT_NEAR(r.mX[1], 6.0, 1e-12);
    EXPECT_NEAR(r.mObjective, -36.0, 1e-12);
    // min x + y s.t. x + y = 1, x - y = 0.5, y >= 0 (equality form).
    const auto e = detail::dense_lp_equality({1.0, 1.0}, {1, 1, 1, -1}, {1.0, 0.5});
    ASSERT_EQ(e.mStatus, DenseLpStatus::Optimal);
    EXPECT_NEAR(e.mX[0], 0.75, 1e-12);
    EXPECT_NEAR(e.mX[1], 0.25, 1e-12);
    // Infeasible and unbounded.
    EXPECT_EQ(detail::dense_lp_inequality({1.0}, {1, -1}, {-1, -1}).mStatus,
              DenseLpStatus::Infeasible);  // x <= -1 and x >= 1
    EXPECT_EQ(detail::dense_lp_inequality({-1.0}, {-1}, {0}).mStatus, DenseLpStatus::Unbounded);
}

TEST(ZoneFit, TheFlatnessOfTwoParallelPlanesIsTheirGap) {
    // Points between z = 0 and z = 0.1, tilted, with points on both planes.
    std::vector<double> pts;
    const double c = std::cos(0.3), s = std::sin(0.3);
    for (int i = 0; i < 20; ++i)
        for (int j = 0; j < 20; ++j) {
            const double x = i * 0.1, y = j * 0.1;
            const double z = (i + j) % 3 == 0
                                 ? 0.0
                                 : ((i * 7 + j * 3) % 5 == 0 ? 0.1 : 0.05 * ((i * j) % 7) / 7.0);
            pts.insert(pts.end(), {x, c * y - s * z, s * y + c * z});
        }
    const double start[10] = {1, 1, 0.05, 0.0, -s + 0.02, c, 0, 0, 0, 0};
    const auto r =
        detail::zone_fit(1, start, pts.data(), pts.size() / 3, detail::ZoneCriterion::MinimumZone);
    ASSERT_TRUE(r.mOk) << r.mWhy;
    EXPECT_NEAR(r.mValue, 0.1, 1e-12);
    EXPECT_NEAR(std::fabs(r.mParams[4]), s, 1e-12);
    EXPECT_NEAR(std::fabs(r.mParams[5]), c, 1e-12);
}

TEST(ZoneFit, EnvelopesOfACircleOfPoints) {
    // Points on a cylinder of radius 1 about z, two bumps of 0.01 outwards
    // and one dent of 0.02 inwards.
    std::vector<double> pts;
    for (int k = 0; k < 360; k += 5)
        for (int h = 0; h < 5; ++h) {
            const double t = k * 3.14159265358979323846 / 180.0;
            const double r = k == 0 ? 1.01 : (k == 180 ? 1.01 : (k == 90 ? 0.98 : 1.0));
            pts.insert(pts.end(), {r * std::cos(t), r * std::sin(t), 0.5 * h});
        }
    const double start[10] = {0.003, -0.002, 0, 0.01, 0.0, 1.0, 1.0, 0, 0, 0};
    const auto mc = detail::zone_fit(3, start, pts.data(), pts.size() / 3,
                                     detail::ZoneCriterion::MinCircumscribed);
    ASSERT_TRUE(mc.mOk) << mc.mWhy;
    EXPECT_NEAR(mc.mValue, 1.01, 1e-9);  // the two bumps are a diameter apart
    const auto mi =
        detail::zone_fit(3, start, pts.data(), pts.size() / 3, detail::ZoneCriterion::MaxInscribed);
    ASSERT_TRUE(mi.mOk) << mi.mWhy;
    EXPECT_GT(mi.mValue, 0.98);  // the axis moves away from the dent
    EXPECT_LT(mi.mValue, 1.0);
    const auto mz =
        detail::zone_fit(3, start, pts.data(), pts.size() / 3, detail::ZoneCriterion::MinimumZone);
    ASSERT_TRUE(mz.mOk) << mz.mWhy;
    EXPECT_LE(mz.mValue, 0.03 + 1e-12);  // never wider than about the true axis
    EXPECT_GT(mz.mValue, 0.0);
    EXPECT_FALSE(detail::zone_fit(1, start, pts.data(), pts.size() / 3,
                                  detail::ZoneCriterion::MinCircumscribed)
                     .mOk);
}

TEST(ZoneFit, FitPrimitivesTakesTheCriterion) {
    SegmentOptions so;
    so.mMethod = SegmentMethod::Primitive;
    const Mesh seg = segment(icosa_sphere(2.0, 3), so).mMesh;
    FitPrimitivesOptions o;
    const FitPrimitivesResult ls = fit_primitives(seg, o);
    o.mCriterion = FitCriterion::MinCircumscribed;
    const FitPrimitivesResult mc = fit_primitives(seg, o);
    o.mCriterion = FitCriterion::MaxInscribed;
    const FitPrimitivesResult mi = fit_primitives(seg, o);
    const auto radius = [](const FitPrimitivesResult& rR) {
        return rR.mMesh.FieldData("fit:patch_params").As<double>()[6];
    };
    // The vertices lie on the sphere of radius 2, the facets inside it.
    EXPECT_NEAR(radius(mc), 2.0, 1e-9);
    EXPECT_LE(radius(mi), radius(ls) + 1e-12);
    EXPECT_GE(radius(mc), radius(ls) - 1e-12);
    EXPECT_EQ(fit_criterion_from_name("minimum_zone"), FitCriterion::MinimumZone);
    EXPECT_THROW(fit_criterion_from_name("chebyshev"), std::invalid_argument);
}
