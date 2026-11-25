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
#pragma once

/**
 * @file detail/sparse_solve.hpp
 * @brief Sparse linear algebra for the mesh operations and skin-to-CAD fits:
 * a compressed-row matrix assembled from triplets, its product with a vector,
 * Jacobi- and IC(0)-preconditioned conjugate gradients (with workspace reuse
 * and multi-RHS block solves), BiCGSTAB for non-symmetric systems, and a
 * zero-dependency RCM-ordered sparse direct LDL^T factorisation for repeated
 * SPD solves on a fixed matrix (roadmap 5.24).
 *
 * A **core-private** header (the `lsq_fit.hpp` precedent): no installed header
 * names it, and it adds nothing to the API or the ABI. Every result is the
 * same on every backend and at every thread count: the rows of a product are
 * independent, every dot product sums fixed chunks of `kSparseChunk` entries
 * in parallel and then the chunk sums serially in index order
 * (`parallel_reduce`), IC(0) and LDL^T factorisations run in deterministic order,
 * and multi-RHS solves run independent columns without cross-column FP drift.
 * Iteration caps are fixed, there is no randomness, and triplets with the same
 * row and column are summed in the order they were given.
 */

// System includes
#include <cstddef>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// The entries every chunk of a deterministic dot product sums.
inline constexpr std::size_t kSparseChunk = 4096;

/// One `(row, column, value)` entry of a matrix being assembled.
struct SparseTriplet {
    std::size_t mRow = 0;
    std::size_t mCol = 0;
    double mValue = 0.0;
};

/// A compressed-row matrix: row `r`'s entries are `mCol` / `mValue` in
/// `[mStart[r], mStart[r + 1])`, columns ascending and unique.
struct SparseCsr {
    std::size_t mRows = 0;
    std::size_t mCols = 0;
    std::vector<std::size_t> mStart;
    std::vector<std::size_t> mCol;
    std::vector<double> mValue;

    /// The stored entries.
    std::size_t NumNonzeros() const { return mValue.size(); }
};

/**
 * @brief Assemble a compressed-row matrix.
 * @param Rows, Cols the shape; every triplet must lie inside it.
 * @param rTriplets the entries, in any order; duplicates are summed in the
 *        order given (a stable sort by row, then column).
 * @throws std::invalid_argument on a triplet outside the shape.
 */
SparseCsr csr_from_triplets(std::size_t Rows, std::size_t Cols,
                            std::vector<SparseTriplet> rTriplets);

/// `pY = A pX`, rows in parallel (each row serial, so deterministic).
void csr_multiply(const SparseCsr& rA, const double* pX, double* pY);

/// The diagonal of a square matrix (zero where none is stored).
std::vector<double> csr_diagonal(const SparseCsr& rA);

/// `sum(pA[i] * pB[i])` over fixed chunks (see the file comment).
double sparse_dot(const double* pA, const double* pB, std::size_t N);

/// How an iterative solve ended.
struct SparseSolveReport {
    bool mConverged = false;
    std::size_t mIterations = 0;
    /// The final `||b - A x|| / ||b||` (0 when `b` is zero).
    double mRelativeResidual = 0.0;
};

/// Reusable Jacobi preconditioner and scratch vectors for repeated `csr_pcg`
/// solves on one fixed matrix (roadmap 5.24.1).
struct CsrPcgWorkspace {
    std::size_t mRows = 0;
    std::vector<double> mMinv;
    std::vector<double> mR;
    std::vector<double> mZ;
    std::vector<double> mP;
    std::vector<double> mAp;

    /// Cache the inverse diagonal and size the scratch buffers for @p rA.
    void Prepare(const SparseCsr& rA);
};

/**
 * @brief Solve a symmetric positive definite system by conjugate gradients
 * with a Jacobi (diagonal) preconditioner.
 * @param rA square, symmetric positive definite.
 * @param pB the right-hand side, `rA.mRows` entries.
 * @param pX the start on entry (zeros is fine), the solution on exit.
 * @param RelativeTolerance stop once `||b - A x|| <= RelativeTolerance ||b||`.
 * @param MaxIterations the cap; 0 takes `2 * rA.mRows + 10`.
 */
SparseSolveReport csr_pcg(const SparseCsr& rA, const double* pB, double* pX,
                          double RelativeTolerance = 1e-12, std::size_t MaxIterations = 0);

/// `csr_pcg` reusing @p rWorkspace (bit-identical to the allocating overload).
SparseSolveReport csr_pcg(const SparseCsr& rA, const double* pB, double* pX,
                          CsrPcgWorkspace& rWorkspace, double RelativeTolerance = 1e-12,
                          std::size_t MaxIterations = 0);

/// Multi-RHS Jacobi-PCG on one fixed matrix @p rA, streaming @p rA once per
/// iteration across active columns (or running independent columns in parallel
/// when `rA.mRows` is small). Bit-identical per column to `csr_pcg`.
std::vector<SparseSolveReport> csr_pcg_block(const SparseCsr& rA, std::size_t NumRhs,
                                             const double* const* ppB, double* const* ppX,
                                             CsrPcgWorkspace& rWorkspace,
                                             double RelativeTolerance = 1e-12,
                                             std::size_t MaxIterations = 0);

/// Zero-fill incomplete Cholesky `IC(0)` preconditioner for a symmetric
/// positive definite CSR matrix (roadmap 5.24.2). Deterministic, no external
/// dependency; retries with diagonal shift and falls back to Jacobi if the
/// matrix is not SPD.
struct CsrIc0 {
    std::size_t mRows = 0;
    std::vector<std::size_t> mLStart;
    std::vector<std::size_t> mLCol;
    std::vector<double> mLVal;
    std::vector<double> mDiagInv;
    std::vector<double> mJacobiInv;
    bool mValid = false;

    /// Build the IC(0) factor for @p rA. Returns true when IC(0) succeeded;
    /// when false, `Apply` falls back to `mJacobiInv`.
    bool Build(const SparseCsr& rA);

    /// Apply `M^-1 r` into `pZ` (`pZ` must not alias `pR`).
    void Apply(const double* pR, double* pZ) const;
};

/// Conjugate gradients preconditioned by @p rPrecond, reusing @p rWorkspace.
SparseSolveReport csr_ic0_pcg(const SparseCsr& rA, const CsrIc0& rPrecond, const double* pB,
                              double* pX, CsrPcgWorkspace& rWorkspace,
                              double RelativeTolerance = 1e-12, std::size_t MaxIterations = 0);

/// Convenience overload: builds `CsrIc0` once and solves.
SparseSolveReport csr_ic0_pcg(const SparseCsr& rA, const double* pB, double* pX,
                              double RelativeTolerance = 1e-12, std::size_t MaxIterations = 0);

/// Multi-RHS IC(0)-preconditioned conjugate gradients on one fixed matrix @p rA.
std::vector<SparseSolveReport> csr_ic0_pcg_block(const SparseCsr& rA, const CsrIc0& rPrecond,
                                                 std::size_t NumRhs, const double* const* ppB,
                                                 double* const* ppX, CsrPcgWorkspace& rWorkspace,
                                                 double RelativeTolerance = 1e-12,
                                                 std::size_t MaxIterations = 0);

/// Cached sparse direct `P A P^T = L D L^T` factorisation with RCM fill-
/// reducing ordering and Davis's up-looking sparse LDL^T (roadmap 5.24.2).
/// Zero external dependency and deterministic across backends and thread counts.
struct CsrLdlt {
    std::size_t mRows = 0;
    std::vector<std::size_t> mPerm;     ///< permuted index -> original index
    std::vector<std::size_t> mInvPerm;  ///< original index -> permuted index
    std::vector<std::size_t> mLp;       ///< column pointers of unit lower L (size n + 1)
    std::vector<std::size_t> mLi;       ///< row indices of L (strictly > col)
    std::vector<double> mLx;            ///< entries of L
    std::vector<double> mDInv;          ///< inverse diagonal entries `1 / D[j]`
    bool mValid = false;

    /// Factor SPD @p rA once. Returns false if symbolic `nnz(L)` exceeds
    /// @p MaxLNonzeros or a non-positive pivot is encountered.
    bool Factor(const SparseCsr& rA, std::size_t MaxLNonzeros = 4000000);

    /// Solve `A x = b` using the cached factorisation and caller scratch @p rWork.
    void Solve(const double* pB, double* pX, std::vector<double>& rWork) const;

    /// Solve `A x = b` using the cached factorisation.
    void Solve(const double* pB, double* pX) const;

    /// Solve `NumRhs` systems independently (parallel across RHS columns).
    void SolveBlock(std::size_t NumRhs, const double* const* ppB, double* const* ppX) const;
};

/**
 * @brief Solve a general square system by Jacobi-preconditioned BiCGSTAB (van
 * der Vorst, 1992). Arguments as `csr_pcg`; a breakdown (a vanishing inner
 * product) ends the solve unconverged with the best iterate so far.
 */
SparseSolveReport csr_bicgstab(const SparseCsr& rA, const double* pB, double* pX,
                               double RelativeTolerance = 1e-12, std::size_t MaxIterations = 0);

}  // namespace detail
}  // namespace meshioplusplus
