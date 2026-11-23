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
 * @file detail/sym3_eigen.hpp
 * @brief Invariants of a symmetric 3x3 tensor stored as six components in the
 * `xx yy zz xy yz zx` order used throughout this codebase (see
 * `doc/mesh_data_model.md`).
 *
 * This was originally private to the CalculiX `.frd` reader's `derived=`
 * option; it now backs `operations/tensor_invariants.hpp` as well, so both
 * callers share one eigensolver instead of two copies drifting apart.
 */

// System includes
#include <cstddef>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief The von Mises equivalent of a symmetric tensor.
 * @param pT six components, `xx yy zz xy yz zx`.
 * @return `sqrt(0.5*((xx-yy)^2+(yy-zz)^2+(zz-xx)^2+6*(xy^2+yz^2+zx^2)))`.
 */
MESHIOPLUSPLUS_API double sym3_mises(const double* pT);

/**
 * @brief Eigenvalues of a symmetric tensor, ascending, by cyclic Jacobi.
 * @param pT six components, `xx yy zz xy yz zx`.
 * @param pOut the three eigenvalues, ascending.
 */
MESHIOPLUSPLUS_API void sym3_principal(const double* pT, double* pOut);

/**
 * @brief Eigenvalues and eigenvectors of a symmetric tensor, by the same cyclic
 * Jacobi as `sym3_principal` (whose values this returns, bit for bit).
 * @param pT six components, `xx yy zz xy yz zx`.
 * @param pValues the three eigenvalues, ascending.
 * @param pVectors the row-major 3x3 whose column `k` (entries `k`, `k + 3`,
 *        `k + 6`) is the unit eigenvector of `pValues[k]`. Each vector's sign is
 *        fixed: its first component of magnitude above `1e-8` is positive. A
 *        repeated eigenvalue has no unique eigenvectors; the pair returned is
 *        one orthonormal basis of its eigenspace.
 */
MESHIOPLUSPLUS_API void sym3_eigen(const double* pT, double* pValues, double* pVectors);

}  // namespace detail
}  // namespace meshioplusplus
