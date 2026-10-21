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
 * @file detail/dense_chol.hpp
 * @brief Symmetric factorisations for the skin-to-CAD fits: a dense Cholesky
 * and LDLᵀ (with a diagonal shift, for Levenberg-Marquardt damping), and a
 * skyline (envelope) LDLᵀ for the banded normal equations of a tensor-product
 * B-spline least-squares fit, whose `n` poles couple only within about
 * `(p + 1) * nu` of each other.
 *
 * A **core-private** header (the `lsq_fit.hpp` precedent): nothing installed
 * names it. Every routine is serial with a fixed operation order, so the
 * result is the same on every backend and at every thread count. No pivoting:
 * the matrices these factor are symmetric positive (semi)definite by
 * construction, and a non-positive pivot is reported, not repaired.
 */

// System includes
#include <cstddef>
#include <vector>

namespace meshioplusplus {
namespace detail {

/**
 * @brief Factor a dense symmetric positive definite `N x N` matrix as `L Lᵀ`.
 * @param rA row-major; only the lower triangle is read, and it is
 *        overwritten by `L` (the strict upper triangle is left as it was).
 * @return false when a pivot is not positive or not finite.
 */
bool chol_factor(std::vector<double>& rA, std::size_t N);

/// Solve `L Lᵀ x = b` in place (`pB` holds `b`, then `x`).
void chol_solve(const std::vector<double>& rL, std::size_t N, double* pB);

/**
 * @brief Factor `A + Shift I` (dense, symmetric) as `L D Lᵀ`, unit `L` in the
 * strict lower triangle and `D` on the diagonal of `rA`.
 * @return false when a pivot is zero or not finite (an indefinite matrix with
 *         non-zero pivots factors: the caller decides what that means).
 */
bool ldlt_factor(std::vector<double>& rA, std::size_t N, double Shift = 0.0);

/// Solve `L D Lᵀ x = b` in place.
void ldlt_solve(const std::vector<double>& rLd, std::size_t N, double* pB);

/**
 * @brief A symmetric matrix stored by rows of its lower triangle, row `i`
 * from column `mFirst[i]` to the diagonal (`mFirst[i] <= i`): the envelope.
 * Fill-in of an LDLᵀ factorisation stays inside the envelope, so the factor
 * overwrites it in place.
 */
struct SkylineMatrix {
    std::size_t mN = 0;
    std::vector<std::size_t> mFirst;
    /// Row `i`'s entries start at `mStart[i]`; the diagonal is at `mStart[i + 1] - 1`.
    std::vector<std::size_t> mStart;
    std::vector<double> mValues;

    /// Whether `(i, j)` (either triangle) lies inside the envelope.
    bool Contains(std::size_t I, std::size_t J) const {
        if (I < J) {
            const std::size_t t = I;
            I = J;
            J = t;
        }
        return J >= mFirst[I];
    }
    /// The entry `(i, j)` (either triangle), which must lie inside the envelope.
    double& At(std::size_t I, std::size_t J) {
        if (I < J) {
            const std::size_t t = I;
            I = J;
            J = t;
        }
        return mValues[mStart[I] + (J - mFirst[I])];
    }
    double At(std::size_t I, std::size_t J) const {
        if (I < J) {
            const std::size_t t = I;
            I = J;
            J = t;
        }
        return mValues[mStart[I] + (J - mFirst[I])];
    }
};

/**
 * @brief A zero matrix with the given envelope.
 * @param rFirst per row, its first stored column (`<=` the row); clamped so
 *        the envelope is monotone (`mFirst[i] >= mFirst[i - 1]` is not
 *        required, but a row never starts after its diagonal).
 */
SkylineMatrix skyline_zero(const std::vector<std::size_t>& rFirst);

/**
 * @brief Factor a skyline matrix plus `Shift I` as `L D Lᵀ` in place: unit `L`
 * below the diagonal, `D` on it.
 * @return false when a pivot is not positive or not finite (the matrices this
 *         factors are positive definite; a regularised least-squares system
 *         with too little data shows up here).
 */
bool skyline_ldlt(SkylineMatrix& rA, double Shift = 0.0);

/// Solve `L D Lᵀ x = b` in place with the factor `skyline_ldlt` left.
void skyline_solve(const SkylineMatrix& rLd, double* pB);

}  // namespace detail
}  // namespace meshioplusplus
