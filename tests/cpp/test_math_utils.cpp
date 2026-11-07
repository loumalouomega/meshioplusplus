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
// Unit tests for detail/math_utils.hpp -- MathUtils, the vector and 3x3
// helpers that replaced hand-written copies across the operations. Each member
// must compute what the copy computed, bit for bit, so the tests compare it
// against the copy's expression as well as against known values.

// System includes
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/math_utils.hpp"

using meshioplusplus::detail::MathUtils;
using meshioplusplus::detail::Vec3;

namespace {

std::uint64_t bits(double X) {
    return std::bit_cast<std::uint64_t>(X);
}

void expect_same_bits(const Vec3& rA, const Vec3& rB) {
    for (std::size_t k = 0; k < 3; ++k)
        EXPECT_EQ(bits(rA[k]), bits(rB[k])) << "component " << k;
}

// Deterministic, awkward inputs: mixed magnitudes and signs.
std::vector<Vec3> sample_vectors() {
    std::mt19937_64 rng(20261107);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<Vec3> out;
    for (int i = 0; i < 64; ++i) {
        const double s = std::pow(10.0, static_cast<double>(i % 9) - 4.0);
        out.push_back({u(rng) * s, u(rng) * s * 3.0, u(rng) / s});
    }
    out.push_back({0.0, 0.0, 0.0});
    out.push_back({-0.0, 0.0, -0.0});
    out.push_back({1e-300, -1e-300, 1e-300});
    return out;
}

}  // namespace

TEST(MathUtils, VectorBasicsKnownValues) {
    const Vec3 a{1.0, 2.0, 3.0}, b{4.0, -5.0, 6.0};
    EXPECT_EQ(MathUtils::Sub(a, b), (Vec3{-3.0, 7.0, -3.0}));
    EXPECT_EQ(MathUtils::Add(a, b), (Vec3{5.0, -3.0, 9.0}));
    EXPECT_EQ(MathUtils::Scale(a, 2.0), (Vec3{2.0, 4.0, 6.0}));
    EXPECT_DOUBLE_EQ(MathUtils::Dot(a, b), 12.0);
    EXPECT_EQ(MathUtils::Cross(a, b), (Vec3{27.0, 6.0, -13.0}));
    EXPECT_DOUBLE_EQ(MathUtils::NormSq(a), 14.0);
    EXPECT_DOUBLE_EQ(MathUtils::Norm(Vec3{3.0, 4.0, 12.0}), 13.0);
    EXPECT_DOUBLE_EQ(MathUtils::Distance(Vec3{1.0, 1.0, 1.0}, Vec3{4.0, 5.0, 1.0}), 5.0);
}

TEST(MathUtils, PointerOverloadsMatchVec3OverloadsBitForBit) {
    const std::vector<Vec3> v = sample_vectors();
    for (std::size_t i = 0; i + 1 < v.size(); ++i) {
        const Vec3& a = v[i];
        const Vec3& b = v[i + 1];
        expect_same_bits(MathUtils::Sub(a.data(), b.data()), MathUtils::Sub(a, b));
        expect_same_bits(MathUtils::Cross(a.data(), b.data()), MathUtils::Cross(a, b));
        Vec3 out{};
        MathUtils::Cross(a.data(), b.data(), out.data());
        expect_same_bits(out, MathUtils::Cross(a, b));
        EXPECT_EQ(bits(MathUtils::Dot(a.data(), b.data())), bits(MathUtils::Dot(a, b)));
        EXPECT_EQ(bits(MathUtils::Dot(a.data(), b)), bits(MathUtils::Dot(a, b)));
        EXPECT_EQ(bits(MathUtils::Distance(a.data(), b.data())), bits(MathUtils::Distance(a, b)));
    }
}

TEST(MathUtils, CrossAndDotKeepTheHandWrittenOrder) {
    // The copies these replaced wrote the components out; the shared ones
    // must round the same way.
    const std::vector<Vec3> v = sample_vectors();
    for (std::size_t i = 0; i + 1 < v.size(); ++i) {
        const Vec3& u = v[i];
        const Vec3& w = v[i + 1];
        const Vec3 copy{u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                        u[0] * w[1] - u[1] * w[0]};
        expect_same_bits(MathUtils::Cross(u, w), copy);
        EXPECT_EQ(bits(MathUtils::Dot(u, w)), bits(u[0] * w[0] + u[1] * w[1] + u[2] * w[2]));
    }
}

TEST(MathUtils, UnitByDivisionIsTheDivisionForm) {
    for (const Vec3& a : sample_vectors()) {
        const double l = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        const Vec3 copy = l > 0.0 ? Vec3{a[0] / l, a[1] / l, a[2] / l} : Vec3{0.0, 0.0, 0.0};
        expect_same_bits(MathUtils::UnitByDivision(a), copy);
    }
    const Vec3 u = MathUtils::UnitByDivision(Vec3{0.0, 3.0, 4.0});
    EXPECT_DOUBLE_EQ(u[1], 0.6);
    EXPECT_DOUBLE_EQ(u[2], 0.8);
    // A zero or NaN length gives the zero vector, never NaN.
    EXPECT_EQ(MathUtils::UnitByDivision(Vec3{0.0, 0.0, 0.0}), (Vec3{0.0, 0.0, 0.0}));
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(MathUtils::UnitByDivision(Vec3{nan, 1.0, 0.0}), (Vec3{0.0, 0.0, 0.0}));
}

TEST(MathUtils, AngleDeg) {
    EXPECT_NEAR(MathUtils::AngleDeg(Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}), 90.0, 1e-12);
    EXPECT_NEAR(MathUtils::AngleDeg(Vec3{1.0, 0.0, 0.0}, Vec3{1.0, 1.0, 0.0}), 45.0, 1e-12);
    EXPECT_NEAR(MathUtils::AngleDeg(Vec3{1.0, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0}), 180.0, 1e-12);
    EXPECT_EQ(MathUtils::AngleDeg(Vec3{2.0, 0.0, 0.0}, Vec3{5.0, 0.0, 0.0}), 0.0);
    // atan2 stays accurate where acos of the normalised dot does not.
    const double tiny = MathUtils::AngleDeg(Vec3{1.0, 0.0, 0.0}, Vec3{1.0, 1e-9, 0.0});
    EXPECT_NEAR(tiny, 1e-9 * 180.0 / MathUtils::kPi, 1e-20);
    for (std::size_t i = 0; i + 1 < sample_vectors().size(); ++i) {
        const Vec3 a = sample_vectors()[i], b = sample_vectors()[i + 1];
        const Vec3 c{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]};
        const double copy = std::atan2(std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]),
                                       a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) *
                            (180.0 / 3.14159265358979323846);
        EXPECT_EQ(bits(MathUtils::AngleDeg(a, b)), bits(copy));
    }
}

TEST(MathUtils, SolveSymmetric3KnownSystem) {
    // M = [[4, 1, 2], [1, 3, 0], [2, 0, 5]], x = (1, -2, 3) => b = M x.
    const double b[3] = {4.0 - 2.0 + 6.0, 1.0 - 6.0, 2.0 + 15.0};
    double x[3] = {0.0, 0.0, 0.0};
    ASSERT_TRUE(MathUtils::SolveSymmetric3(4.0, 1.0, 2.0, 3.0, 0.0, 5.0, b, x));
    EXPECT_NEAR(x[0], 1.0, 1e-14);
    EXPECT_NEAR(x[1], -2.0, 1e-14);
    EXPECT_NEAR(x[2], 3.0, 1e-14);
}

TEST(MathUtils, SolveSymmetric3RefusesSingularAndLeavesTheOutputAlone) {
    const double b[3] = {1.0, 2.0, 3.0};
    double x[3] = {7.0, 8.0, 9.0};
    // Rank 2: the third row is the sum of the first two.
    EXPECT_FALSE(MathUtils::SolveSymmetric3(1.0, 0.0, 1.0, 1.0, 1.0, 2.0, b, x));
    EXPECT_EQ(x[0], 7.0);
    EXPECT_EQ(x[1], 8.0);
    EXPECT_EQ(x[2], 9.0);
    // All zero, and NaN entries: never solved.
    EXPECT_FALSE(MathUtils::SolveSymmetric3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, b, x));
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(MathUtils::SolveSymmetric3(nan, 0.0, 0.0, 1.0, 0.0, 1.0, b, x));
    // Scale-relative: a well-conditioned system at a tiny scale still solves
    // (while its determinant, 1e-150, is still a normal double).
    double y[3];
    EXPECT_TRUE(MathUtils::SolveSymmetric3(1e-50, 0.0, 0.0, 1e-50, 0.0, 1e-50, b, y));
    EXPECT_DOUBLE_EQ(y[2], 3e50);
}

TEST(MathUtils, SolveSymmetric3MatchesTheDecimateAndGradientCopies) {
    // gradient's copy, verbatim (row-major M, no negation).
    const auto copy = [](const double* pM, const double* pB, double* pOut) {
        const double c00 = pM[4] * pM[8] - pM[5] * pM[5];
        const double c01 = pM[2] * pM[5] - pM[1] * pM[8];
        const double c02 = pM[1] * pM[5] - pM[2] * pM[4];
        const double det = pM[0] * c00 + pM[1] * c01 + pM[2] * c02;
        double scale = std::abs(pM[0]);
        scale = std::max(scale, std::abs(pM[1]));
        scale = std::max(scale, std::abs(pM[2]));
        scale = std::max(scale, std::abs(pM[4]));
        scale = std::max(scale, std::abs(pM[5]));
        scale = std::max(scale, std::abs(pM[8]));
        if (!(std::abs(det) > 1e-12 * (scale * scale * scale)))
            return false;
        const double c11 = pM[0] * pM[8] - pM[2] * pM[2];
        const double c12 = pM[1] * pM[2] - pM[0] * pM[5];
        const double c22 = pM[0] * pM[4] - pM[1] * pM[1];
        const double inv = 1.0 / det;
        pOut[0] = (c00 * pB[0] + c01 * pB[1] + c02 * pB[2]) * inv;
        pOut[1] = (c01 * pB[0] + c11 * pB[1] + c12 * pB[2]) * inv;
        pOut[2] = (c02 * pB[0] + c12 * pB[1] + c22 * pB[2]) * inv;
        return true;
    };
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> u(-2.0, 2.0);
    for (int t = 0; t < 500; ++t) {
        double m[9];
        for (int r = 0; r < 3; ++r)
            for (int c = r; c < 3; ++c)
                m[r * 3 + c] = m[c * 3 + r] = u(rng) * (t % 7 == 0 ? 1e-150 : 1.0);
        const double b[3] = {u(rng), u(rng), t % 5 == 0 ? 0.0 : u(rng)};
        double want[3] = {0, 0, 0}, got[3] = {0, 0, 0};
        const bool ok_want = copy(m, b, want);
        const bool ok_got = MathUtils::SolveSymmetric3(m[0], m[1], m[2], m[4], m[5], m[8], b, got);
        ASSERT_EQ(ok_want, ok_got);
        for (int k = 0; k < 3; ++k)
            EXPECT_EQ(bits(want[k]), bits(got[k]));
        // decimate's form: -(adj . b) * inv, the kernel's result negated.
        if (ok_got)
            for (int k = 0; k < 3; ++k) {
                const double neg = -got[k];
                EXPECT_EQ(bits(neg), bits(-(want[k])));
            }
    }
}
