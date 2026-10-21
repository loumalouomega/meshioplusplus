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
#pragma once

/**
 * @file detail/sparse_solve.hpp
 * @brief Sparse linear algebra for the skin-to-CAD fits: a compressed-row
 * matrix assembled from triplets, its product with a vector, and Jacobi-
 * preconditioned conjugate gradients (symmetric positive definite systems,
 * such as least-squares conformal maps) and BiCGSTAB (the non-symmetric
 * mean-value parametrisation).
 *
 * A **core-private** header (the `lsq_fit.hpp` precedent): no installed header
 * names it, and it adds nothing to the API or the ABI. Every result is the
 * same on every backend and at every thread count: the rows of a product are
 * independent, and every dot product sums fixed chunks of `kSparseChunk`
 * entries in parallel and then the chunk sums serially in index order
 * (`parallel_reduce`), so no floating-point sum depends on the scheduling.
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

/**
 * @brief Solve a general square system by Jacobi-preconditioned BiCGSTAB (van
 * der Vorst, 1992). Arguments as `csr_pcg`; a breakdown (a vanishing inner
 * product) ends the solve unconverged with the best iterate so far.
 */
SparseSolveReport csr_bicgstab(const SparseCsr& rA, const double* pB, double* pX,
                               double RelativeTolerance = 1e-12, std::size_t MaxIterations = 0);

}  // namespace detail
}  // namespace meshioplusplus
