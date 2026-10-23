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
// detail/predicates: the exact sign of orient2d, incircle, orient3d and
// insphere. The reference is the same determinant in 128-bit integer
// arithmetic, on integer (or scaled dyadic) coordinates whose products the
// doubles cannot hold exactly, so the floating-point filter is uncertain and
// the exact stage decides; degenerate configurations (collinear, cocircular,
// coplanar, cospherical) must come back exactly zero. Kettner et al.'s
// "classroom examples" grid (Computational Geometry 40, 2008) is included, on
// which the naive determinant is wrong.

// System includes
#include <array>
#include <cmath>
#include <cstdint>
#include <random>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/predicates.hpp"

namespace {

namespace mio = meshioplusplus;

int sign_of(double V) {
    return (V > 0.0) - (V < 0.0);
}

#ifdef __SIZEOF_INT128__
using I128 = __int128;

int sign_of(I128 V) {
    return (V > 0) - (V < 0);
}

I128 ii(double V) {
    return static_cast<I128>(V);
}

I128 exact_orient2d(const double* pA, const double* pB, const double* pC) {
    return (ii(pA[0]) - ii(pC[0])) * (ii(pB[1]) - ii(pC[1])) -
           (ii(pA[1]) - ii(pC[1])) * (ii(pB[0]) - ii(pC[0]));
}

I128 exact_incircle(const double* pA, const double* pB, const double* pC, const double* pD) {
    const I128 adx = ii(pA[0]) - ii(pD[0]), ady = ii(pA[1]) - ii(pD[1]);
    const I128 bdx = ii(pB[0]) - ii(pD[0]), bdy = ii(pB[1]) - ii(pD[1]);
    const I128 cdx = ii(pC[0]) - ii(pD[0]), cdy = ii(pC[1]) - ii(pD[1]);
    return (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) +
           (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy) +
           (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
}

I128 det3(const I128 (&rM)[3][3]) {
    return rM[0][0] * (rM[1][1] * rM[2][2] - rM[1][2] * rM[2][1]) -
           rM[0][1] * (rM[1][0] * rM[2][2] - rM[1][2] * rM[2][0]) +
           rM[0][2] * (rM[1][0] * rM[2][1] - rM[1][1] * rM[2][0]);
}

I128 exact_orient3d(const double* pA, const double* pB, const double* pC, const double* pD) {
    I128 m[3][3];
    const double* p[3] = {pA, pB, pC};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            m[i][k] = ii(p[i][k]) - ii(pD[k]);
    return det3(m);
}

// insphere(a, b, c, d, e) = det[[a-e, |a-e|^2], ..., [d-e, |d-e|^2]], which is
// positive for e inside when orient3d(a, b, c, d) > 0.
I128 exact_insphere(const double* pA, const double* pB, const double* pC, const double* pD,
                    const double* pE) {
    I128 r[4][4];
    const double* p[4] = {pA, pB, pC, pD};
    for (int i = 0; i < 4; ++i) {
        I128 lift = 0;
        for (int k = 0; k < 3; ++k) {
            r[i][k] = ii(p[i][k]) - ii(pE[k]);
            lift += r[i][k] * r[i][k];
        }
        r[i][3] = lift;
    }
    I128 det = 0;
    for (int col = 0; col < 4; ++col) {
        I128 m[3][3];
        for (int i = 1; i < 4; ++i)
            for (int k = 0, c = 0; k < 4; ++k)
                if (k != col)
                    m[i - 1][c++] = r[i][k];
        const I128 term = r[0][col] * det3(m);
        det += (col % 2 == 0) ? term : -term;
    }
    return det;
}
#endif

TEST(Predicates, TwoProductIsExact) {
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> u(-1e10, 1e10);
    for (int i = 0; i < 100000; ++i) {
        const double a = u(rng);
        const double b = u(rng);
        double x = 0.0, y = 0.0;
        mio::detail::predicates_two_product(a, b, x, y);
        EXPECT_EQ(x, a * b);
        EXPECT_EQ(y, std::fma(a, b, -x));
    }
}

TEST(Predicates, Conventions) {
    const double a[3] = {0.0, 0.0, 0.0};
    const double b[3] = {1.0, 0.0, 0.0};
    const double c[3] = {0.0, 1.0, 0.0};
    const double in2[2] = {0.25, 0.25};
    const double out2[2] = {2.0, 2.0};
    const double below[3] = {0.0, 0.0, -1.0};
    const double centre[3] = {0.2, 0.2, -0.2};
    const double far[3] = {5.0, 5.0, 5.0};
    EXPECT_GT(mio::detail::orient2d(a, b, c), 0.0);
    EXPECT_LT(mio::detail::orient2d(a, c, b), 0.0);
    EXPECT_GT(mio::detail::incircle(a, b, c, in2), 0.0);
    EXPECT_LT(mio::detail::incircle(a, b, c, out2), 0.0);
    EXPECT_GT(mio::detail::orient3d(a, b, c, below), 0.0);
    EXPECT_GT(mio::detail::insphere(a, b, c, below, centre), 0.0);
    EXPECT_LT(mio::detail::insphere(a, b, c, below, far), 0.0);
}

#ifdef __SIZEOF_INT128__
// Kettner, Mehlhorn, Pion, Schirra and Yap: p = (0.5 + i u, 0.5 + j u),
// u = 2^-53, against q = (12, 12) and r = (24, 24). Scaled by 2^53 the
// coordinates are integers, so the 128-bit determinant is the truth.
TEST(Predicates, Orient2dClassroomGrid) {
    const double u = std::ldexp(1.0, -53);
    const double q[2] = {12.0, 12.0};
    const double r[2] = {24.0, 24.0};
    const double s = std::ldexp(1.0, 53);
    const double qs[2] = {q[0] * s, q[1] * s};
    const double rs[2] = {r[0] * s, r[1] * s};
    int naive_wrong = 0;
    for (int i = 0; i < 256; ++i) {
        for (int j = 0; j < 256; ++j) {
            const double p[2] = {0.5 + i * u, 0.5 + j * u};
            const double ps[2] = {p[0] * s, p[1] * s};
            const int truth = sign_of(exact_orient2d(ps, qs, rs));
            ASSERT_EQ(sign_of(mio::detail::orient2d(p, q, r)), truth) << i << " " << j;
            const double naive = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
            naive_wrong += sign_of(-naive) != truth;  // orient2d(p, q, r) = -naive
        }
    }
    EXPECT_GT(naive_wrong, 0) << "the grid no longer exercises the exact stage";
}

TEST(Predicates, Orient2dRandomAndCollinear) {
    std::mt19937_64 rng(11);
    std::uniform_int_distribution<std::int64_t> big(-(std::int64_t{1} << 40),
                                                    std::int64_t{1} << 40);
    std::uniform_int_distribution<std::int64_t> step(-(1 << 18), 1 << 18);
    for (int it = 0; it < 20000; ++it) {
        const double a[2] = {double(big(rng)), double(big(rng))};
        const double d[2] = {double(step(rng)), double(step(rng))};
        const std::int64_t k = step(rng) % 64;
        // c on the line through a and b exactly, then nudged by -1, 0 or +1.
        const double b[2] = {a[0] + d[0], a[1] + d[1]};
        const double nudge = double(it % 3 - 1);
        const double c[2] = {a[0] + k * d[0] + nudge, a[1] + k * d[1]};
        ASSERT_EQ(sign_of(mio::detail::orient2d(a, b, c)), sign_of(exact_orient2d(a, b, c)));
        const double e[2] = {double(big(rng)), double(big(rng))};
        ASSERT_EQ(sign_of(mio::detail::orient2d(a, b, e)), sign_of(exact_orient2d(a, b, e)));
    }
}

TEST(Predicates, IncircleRandomAndCocircular) {
    // Integer points on the circle of radius 5^4 = 625 (from Gaussian
    // integers), far from the origin so the products are inexact.
    const std::array<std::array<int, 2>, 8> circle = {{{625, 0},
                                                       {0, 625},
                                                       {-375, 500},
                                                       {500, -375},
                                                       {-527, -336},
                                                       {336, 527},
                                                       {-600, 175},
                                                       {175, -600}}};
    const double off[2] = {double(std::int64_t{1} << 29), -double(std::int64_t{1} << 28)};
    std::mt19937_64 rng(13);
    std::uniform_int_distribution<int> pick(0, 7);
    for (int it = 0; it < 5000; ++it) {
        int i = pick(rng), j = pick(rng), k = pick(rng), l = pick(rng);
        if (i == j || j == k || i == k)
            continue;
        const double a[2] = {off[0] + circle[i][0], off[1] + circle[i][1]};
        const double b[2] = {off[0] + circle[j][0], off[1] + circle[j][1]};
        const double c[2] = {off[0] + circle[k][0], off[1] + circle[k][1]};
        const double nudge = double(it % 3 - 1);
        const double d[2] = {off[0] + circle[l][0] + nudge, off[1] + circle[l][1]};
        const I128 truth = exact_incircle(a, b, c, d);
        if (nudge == 0.0) {
            ASSERT_EQ(truth, 0);
        }
        ASSERT_EQ(sign_of(mio::detail::incircle(a, b, c, d)), sign_of(truth));
    }
    std::uniform_int_distribution<std::int64_t> big(-(std::int64_t{1} << 30),
                                                    std::int64_t{1} << 30);
    for (int it = 0; it < 20000; ++it) {
        double p[4][2];
        for (auto& r_p : p)
            r_p[0] = double(big(rng)), r_p[1] = double(big(rng));
        ASSERT_EQ(sign_of(mio::detail::incircle(p[0], p[1], p[2], p[3])),
                  sign_of(exact_incircle(p[0], p[1], p[2], p[3])));
    }
}

TEST(Predicates, Orient3dRandomAndCoplanar) {
    std::mt19937_64 rng(17);
    std::uniform_int_distribution<std::int64_t> big(-(std::int64_t{1} << 38),
                                                    std::int64_t{1} << 38);
    std::uniform_int_distribution<std::int64_t> small(-(1 << 12), 1 << 12);
    for (int it = 0; it < 20000; ++it) {
        double a[3], u[3], v[3];
        for (int k = 0; k < 3; ++k)
            a[k] = double(big(rng)), u[k] = double(small(rng)), v[k] = double(small(rng));
        const double s = double(small(rng) % 32), t = double(small(rng) % 32);
        double b[3], c[3], d[3];
        for (int k = 0; k < 3; ++k) {
            b[k] = a[k] + u[k];
            c[k] = a[k] + v[k];
            d[k] = a[k] + s * u[k] + t * v[k];
        }
        d[2] += double(it % 3 - 1);
        ASSERT_EQ(sign_of(mio::detail::orient3d(a, b, c, d)), sign_of(exact_orient3d(a, b, c, d)));
        double e[3];
        for (double& r_e : e)
            r_e = double(big(rng));
        ASSERT_EQ(sign_of(mio::detail::orient3d(a, b, c, e)), sign_of(exact_orient3d(a, b, c, e)));
    }
}

TEST(Predicates, InsphereRandomAndCospherical) {
    // Integer points on the sphere of radius 9 (1^2 + 4^2 + 8^2 = 81 and
    // permutations, 0^2 + 0^2 + 9^2), offset so the products are inexact.
    const std::array<std::array<int, 3>, 8> sphere = {{{1, 4, 8},
                                                       {-8, 1, 4},
                                                       {4, -8, -1},
                                                       {0, 0, 9},
                                                       {0, -9, 0},
                                                       {9, 0, 0},
                                                       {-4, -1, -8},
                                                       {8, -4, 1}}};
    const double off[3] = {double(1 << 21), -double(1 << 20), double(1 << 19)};
    std::mt19937_64 rng(19);
    std::uniform_int_distribution<int> pick(0, 7);
    int zeros = 0;
    for (int it = 0; it < 4000; ++it) {
        int idx[5];
        for (int& r_i : idx)
            r_i = pick(rng);
        bool distinct = true;
        for (int x = 0; x < 4; ++x)
            for (int y = x + 1; y < 4; ++y)
                distinct = distinct && idx[x] != idx[y];
        if (!distinct)
            continue;
        double p[5][3];
        for (int x = 0; x < 5; ++x)
            for (int k = 0; k < 3; ++k)
                p[x][k] = off[k] + sphere[idx[x]][k];
        p[4][0] += double(it % 3 - 1);
        const I128 truth = exact_insphere(p[0], p[1], p[2], p[3], p[4]);
        zeros += truth == 0;
        ASSERT_EQ(sign_of(mio::detail::insphere(p[0], p[1], p[2], p[3], p[4])), sign_of(truth));
    }
    EXPECT_GT(zeros, 0);
    std::uniform_int_distribution<std::int64_t> big(-(1 << 22), 1 << 22);
    for (int it = 0; it < 20000; ++it) {
        double p[5][3];
        for (auto& r_p : p)
            for (double& r_c : r_p)
                r_c = double(big(rng));
        ASSERT_EQ(sign_of(mio::detail::insphere(p[0], p[1], p[2], p[3], p[4])),
                  sign_of(exact_insphere(p[0], p[1], p[2], p[3], p[4])));
    }
}
#endif

TEST(Predicates, PermutationSigns) {
    std::mt19937_64 rng(23);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (int it = 0; it < 20000; ++it) {
        double p[5][3];
        for (auto& r_p : p)
            for (double& r_c : r_p)
                r_c = u(rng);
        // Make some inputs exactly degenerate: p[2] the midpoint of p[0], p[1].
        if (it % 4 == 0)
            for (int k = 0; k < 3; ++k)
                p[2][k] = 0.5 * (p[0][k] + p[1][k]);
        const int o = sign_of(mio::detail::orient2d(p[0], p[1], p[2]));
        EXPECT_EQ(sign_of(mio::detail::orient2d(p[1], p[2], p[0])), o);
        EXPECT_EQ(sign_of(mio::detail::orient2d(p[1], p[0], p[2])), -o);
        const int ic = sign_of(mio::detail::incircle(p[0], p[1], p[2], p[3]));
        EXPECT_EQ(sign_of(mio::detail::incircle(p[1], p[2], p[0], p[3])), ic);
        EXPECT_EQ(sign_of(mio::detail::incircle(p[1], p[0], p[2], p[3])), -ic);
        const int o3 = sign_of(mio::detail::orient3d(p[0], p[1], p[2], p[3]));
        EXPECT_EQ(sign_of(mio::detail::orient3d(p[1], p[0], p[2], p[3])), -o3);
        EXPECT_EQ(sign_of(mio::detail::orient3d(p[0], p[1], p[3], p[2])), -o3);
        const int is = sign_of(mio::detail::insphere(p[0], p[1], p[2], p[3], p[4]));
        EXPECT_EQ(sign_of(mio::detail::insphere(p[1], p[0], p[2], p[3], p[4])), -is);
    }
}

}  // namespace
