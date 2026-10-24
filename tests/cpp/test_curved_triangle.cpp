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

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/curved_triangle.hpp"

namespace {

namespace mdt = meshioplusplus::detail;
using V = mdt::CurvedVec3;

constexpr double kThird = 1.0 / 3.0;

/// The reference positions of the nodes, in the element's layout.
std::vector<V> reference(int Order) {
    if (Order == 2)
        return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.5, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0}};
    return {{0, 0, 0},
            {1, 0, 0},
            {0, 1, 0},
            {kThird, 0, 0},
            {2 * kThird, 0, 0},
            {2 * kThird, kThird, 0},
            {kThird, 2 * kThird, 0},
            {0, 2 * kThird, 0},
            {0, kThird, 0},
            {kThird, kThird, 0}};
}

constexpr V kUp{0.0, 0.0, 1.0};

/// The smallest J over a fine lattice, a sampled bound from above.
double sampled_min(int Order, const std::vector<V>& rN) {
    double m = 1e300;
    const int n = 90;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; i + j <= n; ++j)
            m = std::min(m, mdt::curved_triangle_jacobian(Order, rN.data(), kUp, double(i) / n,
                                                          double(j) / n));
    return m;
}

TEST(CurvedTriangle, TheShapeFunctionsInterpolateTheirNodes) {
    // On the reference element each node maps to itself, and the identity
    // element has J = 1 everywhere.
    for (const int order : {2, 3}) {
        const std::vector<V> ref = reference(order);
        for (const V& r_x : ref) {
            const V x = mdt::curved_triangle_point(order, ref.data(), r_x[0], r_x[1]);
            for (int k = 0; k < 3; ++k)
                EXPECT_NEAR(x[k], r_x[k], 1e-14) << order;
        }
        for (const double xi : {0.0, 0.2, 0.7})
            for (const double eta : {0.0, 0.1, 0.3})
                if (xi + eta <= 1.0)
                    EXPECT_NEAR(mdt::curved_triangle_jacobian(order, ref.data(), kUp, xi, eta), 1.0,
                                1e-13);
        const auto c = mdt::curved_triangle_check(order, ref.data(), kUp);
        EXPECT_TRUE(c.mValid);
        EXPECT_TRUE(c.mDecided);
        // The other side: against the opposite normal, invalid.
        EXPECT_FALSE(mdt::curved_triangle_check(order, ref.data(), {0.0, 0.0, -1.0}).mValid);
    }
}

TEST(CurvedTriangle, TheBernsteinFormReproducesAPolynomial) {
    // p = 1 + xi - 2 eta^2 (degree 2) and xi^4 + eta (degree 4): the
    // coefficients sum back to p through the Bernstein basis.
    for (const int d : {2, 4}) {
        std::vector<std::array<int, 2>> lat;
        for (int i = 0; i <= d; ++i)
            for (int j = 0; i + j <= d; ++j)
                lat.push_back({i, j});
        const auto p = [d](double x, double y) {
            return d == 2 ? 1.0 + x - 2.0 * y * y : x * x * x * x + y;
        };
        std::vector<double> vals, coef(lat.size());
        for (const auto& r_l : lat)
            vals.push_back(p(double(r_l[0]) / d, double(r_l[1]) / d));
        mdt::curved_triangle_bernstein(d, vals.data(), coef.data());
        for (const auto& [x, y] : {std::array<double, 2>{0.1, 0.2}, {0.5, 0.25}, {0.0, 0.9}}) {
            double s = 0.0;
            for (std::size_t b = 0; b < lat.size(); ++b) {
                const int i = lat[b][0], j = lat[b][1], k = d - i - j;
                double f = 1.0;
                for (int m = 2; m <= d; ++m)
                    f *= m;
                for (const int m : {i, j, k})
                    for (int q = 2; q <= m; ++q)
                        f /= q;
                s += coef[b] * f * std::pow(x, i) * std::pow(y, j) * std::pow(1.0 - x - y, k);
            }
            EXPECT_NEAR(s, p(x, y), 1e-12) << d;
        }
    }
}

TEST(CurvedTriangle, AnElementPositiveAtItsNodesThatFoldsInsideIsInvalid) {
    // Roadmap 6.6.1: the nodal check is necessary, not sufficient. These
    // elements have J > 0 at every node and J < 0 inside. (With only one
    // mid-edge node moved, a P2 triangle positive at its nodes was found
    // valid at every position of a 121 x 121 grid; three moved nodes fold.)
    std::vector<V> p2 = reference(2);
    p2[3] = {0.0, -0.07, 0.0};
    p2[4] = {0.8, 0.47, 0.0};
    p2[5] = {-0.2, 0.0, 0.0};
    std::vector<V> p3 = reference(3);
    const double moved[7][2] = {{0.82, -0.18}, {0.66, 0.4},  {0.71, 0.59}, {0.44, 0.37},
                                {-0.17, 0.72}, {-0.22, 0.6}, {0.45, -0.04}};
    for (int i = 0; i < 7; ++i)
        p3[static_cast<std::size_t>(3 + i)] = {moved[i][0], moved[i][1], 0.0};
    for (const auto& [order, r_n] : {std::pair<int, const std::vector<V>*>{2, &p2}, {3, &p3}}) {
        const auto c = mdt::curved_triangle_check(order, r_n->data(), kUp);
        EXPECT_GT(c.mMinNodal, 0.0) << order;               // every node passes
        EXPECT_LT(sampled_min(order, *r_n), 0.0) << order;  // yet it folds
        EXPECT_FALSE(c.mValid) << order;
        EXPECT_TRUE(c.mDecided) << order;
    }
}

TEST(CurvedTriangle, AMidEdgeNodePushedAcrossTheOppositeEdgeIsInvalid) {
    // The roadmap's own probe: the mid-edge node of (0, 1) pushed past the
    // opposite corner. That already makes a nodal J non-positive.
    std::vector<V> p2 = reference(2);
    p2[3] = {0.0, 1.4, 0.0};
    const auto c = mdt::curved_triangle_check(2, p2.data(), kUp);
    EXPECT_LE(c.mMinNodal, 0.0);
    EXPECT_FALSE(c.mValid);
}

TEST(CurvedTriangle, ACurvedValidElementIsProvedValid) {
    // A P2 and a P3 triangle bent onto a cylinder of radius 2 around the x
    // axis (z = 2 - sqrt(4 - y^2)), its normal the straight element's: valid,
    // decided without reaching the depth limit.
    for (const int order : {2, 3}) {
        std::vector<V> n = reference(order);
        for (V& r_x : n)
            r_x[2] = 2.0 - std::sqrt(4.0 - r_x[1] * r_x[1]);
        const auto c = mdt::curved_triangle_check(order, n.data(), kUp, 2);
        EXPECT_TRUE(c.mValid) << order;
        EXPECT_TRUE(c.mDecided) << order;
        EXPECT_GT(sampled_min(order, n), 0.0) << order;
    }
}

}  // namespace
