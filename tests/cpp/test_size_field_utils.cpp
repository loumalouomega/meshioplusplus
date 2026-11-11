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

// System includes
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/size_field_utils.hpp"

namespace {

using meshioplusplus::detail::SizeFieldUtils;

TEST(SizeFieldUtils, MetricLengthClosedForms) {
    // uniform sizes: Euclidean over h
    EXPECT_DOUBLE_EQ(SizeFieldUtils::EdgeLength({0, 0, 0}, {2, 0, 0}, 0.5, 0.5), 4.0);
    // geometric sizes: |ab| (hb - ha) / (ha hb ln(hb / ha)), symmetric
    const double l = SizeFieldUtils::MetricLength(1.0, 0.1, 0.4);
    EXPECT_NEAR(l, 0.3 / (0.04 * std::log(4.0)), 1e-12);
    EXPECT_NEAR(SizeFieldUtils::MetricLength(1.0, 0.4, 0.1), l, 1e-12);
    // the metric midpoint halves the metric length, at size sqrt-free r^t
    const double t = SizeFieldUtils::MidpointParameter(0.1, 0.4);
    EXPECT_GT(t, 0.0);
    EXPECT_LT(t, 0.5);  // nearer the small end
    const double hm = SizeFieldUtils::Interpolate(0.1, 0.4, t);
    EXPECT_NEAR(SizeFieldUtils::MetricLength(t, 0.1, hm), 0.5 * l, 1e-9);
    EXPECT_NEAR(SizeFieldUtils::MetricLength(1.0 - t, hm, 0.4), 0.5 * l, 1e-9);
    EXPECT_DOUBLE_EQ(SizeFieldUtils::MidpointParameter(0.3, 0.3), 0.5);
}

TEST(SizeFieldUtils, GradationBoundsTheGrowthAlongEveryEdge) {
    // a line of points 0.1 apart, one small size at the start
    const std::size_t n = 50;
    std::vector<double> xyz(n * 3, 0.0);
    std::vector<std::array<std::int64_t, 2>> edges;
    for (std::size_t i = 0; i < n; ++i) {
        xyz[i * 3] = 0.1 * static_cast<double>(i);
        if (i + 1 < n)
            edges.push_back({static_cast<std::int64_t>(i), static_cast<std::int64_t>(i + 1)});
    }
    std::vector<double> h(n, 10.0);
    h[0] = 0.01;
    const double hgrad = 1.3;
    const std::int64_t lowered = SizeFieldUtils::Grade(xyz, edges, hgrad, h);
    EXPECT_EQ(lowered, static_cast<std::int64_t>(n - 1));
    for (const auto& e : edges) {
        const double d = 0.1;
        EXPECT_LE(h[static_cast<std::size_t>(e[1])],
                  h[static_cast<std::size_t>(e[0])] + (hgrad - 1) * d + 1e-12);
        EXPECT_LE(h[static_cast<std::size_t>(e[0])],
                  h[static_cast<std::size_t>(e[1])] + (hgrad - 1) * d + 1e-12);
    }
    EXPECT_NEAR(h[n - 1], 0.01 + 0.3 * 0.1 * static_cast<double>(n - 1), 1e-9);
}

}  // namespace
