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
 * @file detail/lsq_fit.hpp
 * @brief The small dense linear algebra `fit_primitives` needs: a symmetric
 * 3x3 eigen-decomposition with vectors, a pivoted solver for a few unknowns,
 * an orthonormal frame around an axis, and a Levenberg-Marquardt loop.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI. Every routine is
 * serial and written so the numpy twin (`_fit_primitives.py`) can follow it
 * step for step: fixed sweep and iteration caps, no randomness, ties broken by
 * index.
 */

// System includes
#include <cstddef>
#include <functional>
#include <vector>

namespace meshioplusplus {
namespace detail {

/**
 * @brief Eigen-decomposition of a symmetric 3x3 matrix by cyclic Jacobi sweeps.
 * @param pA the matrix, row-major (9 entries); only read.
 * @param pW the eigenvalues, ascending (ties keep their index order).
 * @param pV the eigenvectors, one per row, matching `pW`; each is signed so its
 *        largest-magnitude component (the first on a tie) is positive.
 */
void lsq_sym3_eigen(const double* pA, double* pW, double* pV);

/**
 * @brief Solve the `K x K` system `A x = b` by Gaussian elimination with
 * partial pivoting.
 * @param rA row-major, overwritten.
 * @param rB the right-hand side on entry, the solution on success.
 * @return false when a pivot is zero or not finite.
 */
bool lsq_solve(std::vector<double>& rA, std::vector<double>& rB, std::size_t K);

/// An orthonormal `(e1, e2)` completing the unit vector `pE3` to a right-handed frame.
void lsq_frame(const double* pE3, double* pE1, double* pE2);

/// Fills `N` residuals for the parameter vector.
using LsqResidual = std::function<void(const std::vector<double>& rTheta, double* pOut)>;

/**
 * @brief Minimise the sum of squared residuals by Levenberg-Marquardt with a
 * forward-difference Jacobian.
 * @param rTheta the start on entry, the minimiser on exit.
 * @param rTypical per parameter, its typical magnitude (sets the difference step).
 * @param N the number of residuals.
 * @param MaxIterations the cap on accepted steps.
 */
void lsq_levenberg_marquardt(std::vector<double>& rTheta, const std::vector<double>& rTypical,
                             std::size_t N, const LsqResidual& rResidual, int MaxIterations);

/**
 * @brief Kasa's algebraic circle through the points `(rX[i], rY[i])`: `x^2 +
 * y^2 + D x + E y + F = 0` by its 3x3 normal equations, solved by Cramer's
 * rule, accumulated in point order (the arithmetic `fit_cad`'s cylinder
 * unroll has always used).
 * @return false when the system is singular or the radius is not real.
 */
bool lsq_circle_kasa(const std::vector<double>& rX, const std::vector<double>& rY, double& rCx,
                     double& rCy, double& rR);

/**
 * @brief The geometric circle: Gauss-Newton on the residuals `|p - c| - r`
 * from the start `(rCx, rCy, rR)` (Kasa's, say), at most `MaxIterations`
 * steps, stopping when a step moves the centre and radius by less than
 * `1e-15` of the radius.
 * @return false when a step is singular; the start is then left unchanged.
 */
bool lsq_circle_geometric(const std::vector<double>& rX, const std::vector<double>& rY, double& rCx,
                          double& rCy, double& rR, int MaxIterations = 50);

}  // namespace detail
}  // namespace meshioplusplus
