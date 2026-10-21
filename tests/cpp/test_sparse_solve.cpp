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
// The OpenCASCADE-free linear algebra of the skin-to-CAD fits: compressed-row
// assembly, conjugate gradients, BiCGSTAB, dense Cholesky/LDLT and the skyline
// LDLT of banded least-squares normal equations.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/dense_chol.hpp"
#include "../../src/cpp/src/detail/sparse_solve.hpp"

using namespace meshioplusplus::detail;

namespace {

// The 5-point Laplacian (Dirichlet) on an N x N grid: symmetric positive definite.
SparseCsr laplacian(std::size_t N) {
    std::vector<SparseTriplet> t;
    const auto id = [N](std::size_t I, std::size_t J) { return I + N * J; };
    for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < N; ++i) {
            t.push_back({id(i, j), id(i, j), 4.0});
            if (i > 0)
                t.push_back({id(i, j), id(i - 1, j), -1.0});
            if (i + 1 < N)
                t.push_back({id(i, j), id(i + 1, j), -1.0});
            if (j > 0)
                t.push_back({id(i, j), id(i, j - 1), -1.0});
            if (j + 1 < N)
                t.push_back({id(i, j), id(i, j + 1), -1.0});
        }
    return csr_from_triplets(N * N, N * N, t);
}

// A deterministic pseudo-random sequence (splitmix64) in [-1, 1).
struct Rng {
    std::uint64_t mState;
    double Next() {
        std::uint64_t z = (mState += 0x9e3779b97f4a7c15ull);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        z ^= z >> 31;
        return static_cast<double>(z >> 11) * (2.0 / 9007199254740992.0) - 1.0;
    }
};

std::vector<double> residual(const SparseCsr& rA, const std::vector<double>& rX,
                             const std::vector<double>& rB) {
    std::vector<double> ax(rB.size());
    csr_multiply(rA, rX.data(), ax.data());
    for (std::size_t i = 0; i < ax.size(); ++i)
        ax[i] = rB[i] - ax[i];
    return ax;
}

double norm(const std::vector<double>& rV) {
    double s = 0.0;
    for (double v : rV)
        s += v * v;
    return std::sqrt(s);
}

}  // namespace

TEST(SparseSolve, CsrFromTripletsSumsDuplicatesInOrder) {
    const SparseCsr a = csr_from_triplets(
        3, 3, {{2, 0, 1.0}, {0, 1, 2.0}, {2, 0, 3.0}, {1, 1, 5.0}, {0, 0, 7.0}, {0, 1, -2.0}});
    ASSERT_EQ(a.mStart, (std::vector<std::size_t>{0, 2, 3, 4}));
    EXPECT_EQ(a.mCol, (std::vector<std::size_t>{0, 1, 1, 0}));
    EXPECT_EQ(a.mValue, (std::vector<double>{7.0, 0.0, 5.0, 4.0}));
    EXPECT_EQ(csr_diagonal(a), (std::vector<double>{7.0, 5.0, 0.0}));
    EXPECT_THROW(csr_from_triplets(2, 2, {{2, 0, 1.0}}), std::invalid_argument);
}

TEST(SparseSolve, PcgSolvesLaplacian) {
    const SparseCsr a = laplacian(40);
    Rng rng{1};
    std::vector<double> b(a.mRows);
    for (double& v : b)
        v = rng.Next();
    std::vector<double> x(a.mRows, 0.0);
    const SparseSolveReport rep = csr_pcg(a, b.data(), x.data(), 1e-12);
    EXPECT_TRUE(rep.mConverged);
    EXPECT_GT(rep.mIterations, 0u);
    EXPECT_LT(norm(residual(a, x, b)) / norm(b), 1e-10);
    // A zero right-hand side gives zero.
    std::vector<double> z(a.mRows, 0.0), xz(a.mRows, 1.0);
    EXPECT_TRUE(csr_pcg(a, z.data(), xz.data()).mConverged);
    EXPECT_EQ(norm(xz), 0.0);
}

TEST(SparseSolve, BicgstabSolvesNonSymmetric) {
    // Convection-diffusion: the Laplacian plus an upwind first difference.
    const std::size_t n = 30;
    std::vector<SparseTriplet> t;
    const auto id = [n](std::size_t I, std::size_t J) { return I + n * J; };
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i) {
            t.push_back({id(i, j), id(i, j), 4.5});
            if (i > 0)
                t.push_back({id(i, j), id(i - 1, j), -1.5});
            if (i + 1 < n)
                t.push_back({id(i, j), id(i + 1, j), -1.0});
            if (j > 0)
                t.push_back({id(i, j), id(i, j - 1), -1.0});
            if (j + 1 < n)
                t.push_back({id(i, j), id(i, j + 1), -1.0});
        }
    const SparseCsr a = csr_from_triplets(n * n, n * n, t);
    Rng rng{7};
    std::vector<double> b(a.mRows);
    for (double& v : b)
        v = rng.Next();
    std::vector<double> x(a.mRows, 0.0);
    const SparseSolveReport rep = csr_bicgstab(a, b.data(), x.data(), 1e-11);
    EXPECT_TRUE(rep.mConverged);
    EXPECT_LT(norm(residual(a, x, b)) / norm(b), 1e-9);
}

TEST(SparseSolve, ResultIndependentOfThreadCount) {
    // The dot product sums fixed chunks: it equals the chunked serial sum bit
    // for bit, whatever the scheduling, and so does a whole solve run twice.
    Rng rng{3};
    const std::size_t n = 5 * kSparseChunk + 17;
    std::vector<double> a(n), b(n);
    for (std::size_t i = 0; i < n; ++i) {
        a[i] = rng.Next();
        b[i] = rng.Next() * 1e3;
    }
    double expect = 0.0;
    for (std::size_t c = 0; c * kSparseChunk < n; ++c) {
        double s = 0.0;
        for (std::size_t i = c * kSparseChunk; i < std::min(n, (c + 1) * kSparseChunk); ++i)
            s += a[i] * b[i];
        expect += s;
    }
    EXPECT_EQ(sparse_dot(a.data(), b.data(), n), expect);

    const SparseCsr lap = laplacian(90);
    std::vector<double> rhs(lap.mRows);
    for (double& v : rhs)
        v = rng.Next();
    std::vector<double> x1(lap.mRows, 0.0), x2(lap.mRows, 0.0);
    csr_pcg(lap, rhs.data(), x1.data(), 1e-10);
    csr_pcg(lap, rhs.data(), x2.data(), 1e-10);
    EXPECT_EQ(x1, x2);
}

TEST(SparseSolve, DenseCholeskyAndLdlt) {
    // A = M Mᵀ + I, SPD.
    const std::size_t n = 6;
    Rng rng{11};
    std::vector<double> m(n * n), a(n * n, 0.0);
    for (double& v : m)
        v = rng.Next();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t k = 0; k < n; ++k)
                a[i * n + j] += m[i * n + k] * m[j * n + k];
            if (i == j)
                a[i * n + j] += 1.0;
        }
    std::vector<double> b(n);
    for (double& v : b)
        v = rng.Next();
    const auto check = [&](const std::vector<double>& rX, double Shift) {
        for (std::size_t i = 0; i < n; ++i) {
            double s = Shift * rX[i];
            for (std::size_t j = 0; j < n; ++j)
                s += a[i * n + j] * rX[j];
            EXPECT_NEAR(s, b[i], 1e-12);
        }
    };
    std::vector<double> l = a, x = b;
    ASSERT_TRUE(chol_factor(l, n));
    chol_solve(l, n, x.data());
    check(x, 0.0);
    std::vector<double> ld = a, y = b;
    ASSERT_TRUE(ldlt_factor(ld, n, 0.5));
    ldlt_solve(ld, n, y.data());
    check(y, 0.5);
    // Not positive definite: Cholesky refuses.
    std::vector<double> bad = {1.0, 2.0, 2.0, 1.0};
    EXPECT_FALSE(chol_factor(bad, 2));
}

TEST(SparseSolve, SkylineLdltMatchesDense) {
    // A banded SPD matrix with a ragged envelope.
    const std::size_t n = 12;
    const std::vector<std::size_t> first = {0, 0, 1, 0, 2, 4, 3, 6, 6, 5, 8, 9};
    Rng rng{5};
    SkylineMatrix s = skyline_zero(first);
    std::vector<double> dense(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = first[i]; j < i; ++j) {
            const double v = 0.3 * rng.Next();
            s.At(i, j) = v;
            dense[i * n + j] = dense[j * n + i] = v;
        }
    for (std::size_t i = 0; i < n; ++i) {
        s.At(i, i) = 4.0;
        dense[i * n + i] = 4.0;
    }
    EXPECT_TRUE(s.Contains(4, 2));
    EXPECT_FALSE(s.Contains(4, 1));
    std::vector<double> b(n);
    for (double& v : b)
        v = rng.Next();
    std::vector<double> xs = b, xd = b;
    ASSERT_TRUE(skyline_ldlt(s, 0.25));
    skyline_solve(s, xs.data());
    for (std::size_t i = 0; i < n; ++i)
        dense[i * n + i] += 0.25;
    ASSERT_TRUE(chol_factor(dense, n));
    chol_solve(dense, n, xd.data());
    for (std::size_t i = 0; i < n; ++i)
        EXPECT_NEAR(xs[i], xd[i], 1e-13);
}
