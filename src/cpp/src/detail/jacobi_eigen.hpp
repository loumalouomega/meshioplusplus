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
 * @file detail/jacobi_eigen.hpp
 * @brief Cyclic Jacobi eigen-decomposition of a small dense symmetric matrix,
 * shared by `eigenmodes` (its Rayleigh-Ritz step) and `registration` (the
 * quaternion and 6 x 6 point-to-plane systems).
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): it lives beside
 * the `.cpp` files, no installed header names it, and it adds nothing to the
 * API or the ABI. Each caller passes its own convergence tolerance, so the
 * arithmetic -- and every digest built on it -- is unchanged.
 */

// System includes
#include <cmath>
#include <cstddef>

namespace meshioplusplus {
namespace detail {

/**
 * @brief Diagonalises the symmetric N x N row-major matrix @p pA in place by
 * cyclic Jacobi sweeps (at most 100).
 *
 * Sweeps stop once the squared off-diagonal mass is zero or at most
 * `RelTol` times the squared diagonal mass.
 * @param pA the matrix (N * N doubles); overwritten by the rotated matrix.
 * @param N the matrix order.
 * @param pW receives the N eigenvalues (the final diagonal, unsorted).
 * @param pV receives the eigenvectors as the columns of an N x N row-major
 *        matrix (N * N doubles; need not be initialised).
 * @param RelTol the relative convergence tolerance.
 */
inline void jacobi_eigen_sym(double* pA, std::size_t N, double* pW, double* pV, double RelTol) {
    for (std::size_t i = 0; i < N * N; ++i)
        pV[i] = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        pV[i * N + i] = 1.0;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0, diag = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            diag += pA[i * N + i] * pA[i * N + i];
            for (std::size_t j = i + 1; j < N; ++j)
                off += pA[i * N + j] * pA[i * N + j];
        }
        if (off == 0.0 || off <= RelTol * diag)
            break;
        for (std::size_t p = 0; p < N; ++p)
            for (std::size_t q = p + 1; q < N; ++q) {
                const double apq = pA[p * N + q];
                if (apq == 0.0)
                    continue;
                const double theta = (pA[q * N + q] - pA[p * N + p]) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (std::size_t k = 0; k < N; ++k) {  // A <- A J
                    const double akp = pA[k * N + p], akq = pA[k * N + q];
                    pA[k * N + p] = c * akp - s * akq;
                    pA[k * N + q] = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < N; ++k) {  // A <- J^T A
                    const double apk = pA[p * N + k], aqk = pA[q * N + k];
                    pA[p * N + k] = c * apk - s * aqk;
                    pA[q * N + k] = s * apk + c * aqk;
                }
                for (std::size_t k = 0; k < N; ++k) {
                    const double vkp = pV[k * N + p], vkq = pV[k * N + q];
                    pV[k * N + p] = c * vkp - s * vkq;
                    pV[k * N + q] = s * vkp + c * vkq;
                }
            }
    }
    for (std::size_t i = 0; i < N; ++i)
        pW[i] = pA[i * N + i];
}

}  // namespace detail
}  // namespace meshioplusplus
