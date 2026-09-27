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

// check_quality / parse_quality_thresholds (roadmap §5 quality gate, v16.24.0).

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/quality_gate.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// Two triangles: a right isosceles one (min angle 45) and a sliver.
Mesh two_triangles() {
    return mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {3, 0, 0}, {5, 0, 0}, {4, 0.1, 0}},
                         "triangle", {{0, 1, 2}, {3, 4, 5}});
}

}  // namespace

TEST(QualityGate, ParsesTheSpecificationText) {
    const auto t = parse_quality_thresholds(
        "scaled_jacobian >= 0.2; quality:aspect_ratio<=5 @ 1%\n# comment, ignored\nskewness <= "
        "0.5 @ 0.25");
    ASSERT_EQ(t.size(), 3u);
    EXPECT_EQ(t[0].mMetric, "quality:scaled_jacobian");
    EXPECT_EQ(t[0].mMin, 0.2);
    EXPECT_TRUE(std::isnan(t[0].mMax));
    EXPECT_EQ(t[1].mMetric, "quality:aspect_ratio");
    EXPECT_EQ(t[1].mMax, 5.0);
    EXPECT_DOUBLE_EQ(t[1].mMaxFraction, 0.01);
    EXPECT_DOUBLE_EQ(t[2].mMaxFraction, 0.25);
    EXPECT_TRUE(parse_quality_thresholds("  ;\n# only a comment\n").empty());
}

TEST(QualityGate, RejectsMalformedClauses) {
    EXPECT_THROW(parse_quality_thresholds("bogus >= 1"), std::invalid_argument);
    EXPECT_THROW(parse_quality_thresholds("min_angle > 1"), std::invalid_argument);
    EXPECT_THROW(parse_quality_thresholds("min_angle >= x"), std::invalid_argument);
    EXPECT_THROW(parse_quality_thresholds("min_angle >= 1 <= 2"), std::invalid_argument);
    QualityGateOptions o;
    o.mThresholds.push_back({"min_angle", std::nan(""), std::nan(""), 0.0});
    EXPECT_THROW(check_quality(two_triangles(), o), std::invalid_argument);
    o.mThresholds[0].mMin = 1.0;
    o.mThresholds[0].mMaxFraction = 2.0;
    EXPECT_THROW(check_quality(two_triangles(), o), std::invalid_argument);
}

TEST(QualityGate, CountsViolationsPerCell) {
    QualityGateOptions o;
    o.mThresholds = parse_quality_thresholds("min_angle >= 30");
    const QualityGateResult r = check_quality(two_triangles(), o);
    ASSERT_EQ(r.mChecks.size(), 3u);  // the threshold, then inverted and degenerate
    const QualityCheck& c = r.mChecks[0];
    EXPECT_EQ(c.mName, "min_angle >= 30");
    EXPECT_EQ(c.mEvaluated, 2);
    EXPECT_EQ(c.mViolations, 1);
    EXPECT_EQ(c.mWorstCell, 1);
    EXPECT_LT(c.mWorst, 10.0);
    EXPECT_FALSE(c.mPassed);
    EXPECT_FALSE(r.mPassed);
    EXPECT_NE(quality_gate_summary(r).find("FAIL"), std::string::npos);
}

TEST(QualityGate, AnAllowedFractionPasses) {
    QualityGateOptions o;
    o.mThresholds = parse_quality_thresholds("min_angle >= 30 @ 50%");
    const QualityGateResult r = check_quality(two_triangles(), o);
    EXPECT_TRUE(r.mChecks[0].mPassed);
    EXPECT_DOUBLE_EQ(r.mChecks[0].mFraction, 0.5);
    EXPECT_TRUE(r.mPassed);
}

TEST(QualityGate, CountLimitsCanBeDisabled) {
    QualityGateOptions o;
    o.mMaxInverted = -1;
    o.mMaxDegenerate = -1;
    const QualityGateResult r = check_quality(two_triangles(), o);
    EXPECT_TRUE(r.mChecks.empty());
    EXPECT_TRUE(r.mPassed);
}

TEST(QualityGate, TheMetricListMatchesComputeQuality) {
    // parse_quality_thresholds validates against its own copy of the metric
    // names; every name compute_quality reports must parse.
    const QualityReport rep = compute_quality(two_triangles());
    for (const auto& kv : rep.mMetrics)
        EXPECT_NO_THROW(parse_quality_thresholds(kv.first + " >= 0")) << kv.first;
    EXPECT_EQ(rep.mMetrics.size(), 11u);
}
