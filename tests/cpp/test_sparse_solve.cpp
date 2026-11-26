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
//                   SPDX-License-Identifier: MIT
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

TEST(SparseSolve, PcgWorkspaceAndBlockMatchScalarBitForBit) {
    const SparseCsr a = laplacian(25);
    const std::size_t n = a.mRows;
    Rng rng{19};
    constexpr std::size_t kCols = 5;
    std::vector<std::vector<double>> rhs(kCols, std::vector<double>(n));
    for (std::size_t c = 0; c < kCols; ++c)
        for (double& v : rhs[c])
            v = (c == 2) ? 0.0 : rng.Next();  // include a zero RHS column

    std::vector<std::vector<double>> x_scalar(kCols, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> x_ws(kCols, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> x_blk(kCols, std::vector<double>(n, 0.0));
    CsrPcgWorkspace ws;
    ws.Prepare(a);
    for (std::size_t c = 0; c < kCols; ++c) {
        const SparseSolveReport rs = csr_pcg(a, rhs[c].data(), x_scalar[c].data(), 1e-11);
        const SparseSolveReport rw = csr_pcg(a, rhs[c].data(), x_ws[c].data(), ws, 1e-11);
        EXPECT_TRUE(rs.mConverged);
        EXPECT_EQ(rw.mIterations, rs.mIterations);
        EXPECT_EQ(x_ws[c], x_scalar[c]);
    }
    std::vector<const double*> bp(kCols);
    std::vector<double*> xp(kCols);
    for (std::size_t c = 0; c < kCols; ++c) {
        bp[c] = rhs[c].data();
        xp[c] = x_blk[c].data();
    }
    const std::vector<SparseSolveReport> rb =
        csr_pcg_block(a, kCols, bp.data(), xp.data(), ws, 1e-11);
    ASSERT_EQ(rb.size(), kCols);
    for (std::size_t c = 0; c < kCols; ++c) {
        EXPECT_TRUE(rb[c].mConverged);
        EXPECT_EQ(x_blk[c], x_scalar[c]);
    }
}

TEST(SparseSolve, Ic0PcgBeatsJacobiAndFallsBackOnIndefinite) {
    const SparseCsr a = laplacian(40);
    const std::size_t n = a.mRows;
    Rng rng{23};
    std::vector<double> b0(n), b1(n), b2(n);
    for (std::size_t i = 0; i < n; ++i) {
        b0[i] = rng.Next();
        b1[i] = rng.Next();
        b2[i] = rng.Next();
    }
    std::vector<double> x_jac(n, 0.0), x_ic0(n, 0.0);
    const SparseSolveReport r_jac = csr_pcg(a, b0.data(), x_jac.data(), 1e-11);
    CsrIc0 precond;
    ASSERT_TRUE(precond.Build(a));
    EXPECT_TRUE(precond.mValid);
    CsrPcgWorkspace ws;
    const SparseSolveReport r_ic0 = csr_ic0_pcg(a, precond, b0.data(), x_ic0.data(), ws, 1e-11);
    EXPECT_TRUE(r_ic0.mConverged);
    EXPECT_LT(r_ic0.mIterations, r_jac.mIterations);
    EXPECT_LT(norm(residual(a, x_ic0, b0)) / norm(b0), 1e-10);

    // Block IC(0)-PCG is bit-identical to single-column IC(0)-PCG.
    std::vector<double> x0(n, 0.0), x1(n, 0.0), x2(n, 0.0);
    const double* bp[3] = {b0.data(), b1.data(), b2.data()};
    double* xp[3] = {x0.data(), x1.data(), x2.data()};
    const std::vector<SparseSolveReport> rb = csr_ic0_pcg_block(a, precond, 3, bp, xp, ws, 1e-11);
    ASSERT_EQ(rb.size(), 3u);
    for (const SparseSolveReport& r : rb)
        EXPECT_TRUE(r.mConverged);
    EXPECT_EQ(x0, x_ic0);
    EXPECT_LT(norm(residual(a, x1, b1)) / norm(b1), 1e-10);
    EXPECT_LT(norm(residual(a, x2, b2)) / norm(b2), 1e-10);

    // Non-SPD matrix: IC(0) build declines cleanly and Apply falls back to Jacobi.
    const SparseCsr bad =
        csr_from_triplets(2, 2, {{0, 0, -2.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, -2.0}});
    CsrIc0 bad_ic;
    EXPECT_FALSE(bad_ic.Build(bad));
    EXPECT_FALSE(bad_ic.mValid);
    double r_in[2] = {4.0, -6.0}, z_out[2] = {0.0, 0.0};
    bad_ic.Apply(r_in, z_out);
    EXPECT_DOUBLE_EQ(z_out[0], -2.0);
    EXPECT_DOUBLE_EQ(z_out[1], 3.0);
}

TEST(SparseSolve, LdltDirectBlockDisconnectedAndGuards) {
    // Two disconnected 2-D Laplacian blocks combined into one matrix.
    const SparseCsr g1 = laplacian(20);
    const SparseCsr g2 = laplacian(15);
    const std::size_t n1 = g1.mRows, n2 = g2.mRows, n = n1 + n2;
    std::vector<SparseTriplet> trip;
    for (std::size_t r = 0; r < n1; ++r)
        for (std::size_t k = g1.mStart[r]; k < g1.mStart[r + 1]; ++k)
            trip.push_back({r, g1.mCol[k], g1.mValue[k]});
    for (std::size_t r = 0; r < n2; ++r)
        for (std::size_t k = g2.mStart[r]; k < g2.mStart[r + 1]; ++k)
            trip.push_back({n1 + r, n1 + g2.mCol[k], g2.mValue[k]});
    const SparseCsr a = csr_from_triplets(n, n, std::move(trip));

    // Fill-budget guard declines without throwing.
    CsrLdlt tight;
    EXPECT_FALSE(tight.Factor(a, 10));
    EXPECT_FALSE(tight.mValid);

    // Non-SPD matrix declines cleanly.
    const SparseCsr indef =
        csr_from_triplets(2, 2, {{0, 0, 1.0}, {0, 1, 3.0}, {1, 0, 3.0}, {1, 1, 1.0}});
    CsrLdlt bad;
    EXPECT_FALSE(bad.Factor(indef));
    EXPECT_FALSE(bad.mValid);

    CsrLdlt ldlt;
    ASSERT_TRUE(ldlt.Factor(a));
    ASSERT_TRUE(ldlt.mValid);
    EXPECT_GT(ldlt.mLi.size(), 0u);

    // SolveBlock (7 RHS: one full 4-wide chunk + 3 remainder) is bit-identical to Solve.
    Rng rng{31};
    constexpr std::size_t kCols = 7;
    std::vector<std::vector<double>> rhs(kCols, std::vector<double>(n));
    std::vector<std::vector<double>> x_single(kCols, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> x_block(kCols, std::vector<double>(n, 0.0));
    std::vector<const double*> bp(kCols);
    std::vector<double*> xp(kCols);
    for (std::size_t c = 0; c < kCols; ++c) {
        for (double& v : rhs[c])
            v = rng.Next();
        ldlt.Solve(rhs[c].data(), x_single[c].data());
        EXPECT_LT(norm(residual(a, x_single[c], rhs[c])) / norm(rhs[c]), 1e-12);
        bp[c] = rhs[c].data();
        xp[c] = x_block[c].data();
    }
    ldlt.SolveBlock(kCols, bp.data(), xp.data());
    for (std::size_t c = 0; c < kCols; ++c)
        EXPECT_EQ(x_block[c], x_single[c]);
}

TEST(SparseSolve, LdltScalesTo100kSurfaceGrid) {
    // 320 x 320 = 102,400 unknowns: Nested Dissection keeps fill O(n log n),
    // well below the 20M default fill budget, and solves to machine precision.
    const SparseCsr a = laplacian(320);
    CsrLdlt ldlt;
    ASSERT_TRUE(ldlt.Factor(a));
    EXPECT_LT(ldlt.mLi.size(), std::size_t{4000000});
    Rng rng{41};
    std::vector<double> b(a.mRows), x(a.mRows, 0.0);
    for (double& v : b)
        v = rng.Next();
    ldlt.Solve(b.data(), x.data());
    EXPECT_LT(norm(residual(a, x, b)) / norm(b), 1e-11);
}
