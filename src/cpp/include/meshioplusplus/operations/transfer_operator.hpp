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
 * @file operations/transfer_operator.hpp
 * @brief Sparse CSR inter-mesh prolongation (`P`: coarse -> fine) and
 * restriction (`R = P^T`: fine -> coarse) transfer operators (roadmap 15.5).
 *
 * Uses exact nested stencils from `refine:entity` when available, or
 * barycentric interpolation via `detail::PointLocator` / surface projection
 * for arbitrary non-nested mesh pairs. Every located row sums to 1 (partition
 * of unity) with deterministic ascending column indices.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/mesh_operator.hpp"

namespace meshioplusplus {

/// How `transfer_operator` builds the coarse-to-fine interpolation weights.
enum class TransferOperatorMode : std::int32_t {
    /// Use `refine:entity` nested stencils when compatible, else barycentric location.
    Auto = 0,
    /// Require `refine:entity` nested stencils (throws if absent/incompatible).
    Nested = 1,
    /// Always use geometric barycentric location (`PointLocator` / surface projection).
    Barycentric = 2,
};

/// Parse `"auto"`, `"nested"`, or `"barycentric"`.
MESHIOPLUSPLUS_API TransferOperatorMode transfer_operator_mode_from_name(const std::string& rName);

/// Options for `transfer_operator`.
struct TransferOperatorOptions {
    /// Weight-construction strategy (`Auto`, `Nested`, or `Barycentric`).
    TransferOperatorMode mMode = TransferOperatorMode::Auto;
    /// If true, return the transpose restriction operator `R = P^T` (`n_coarse x n_fine`).
    bool mTranspose = false;
    /// Geometric containment tolerance for barycentric location.
    double mTolerance = 1e-12;
    /// If true, fall back to nearest coarse vertex (weight 1) for fine points outside `coarse`.
    bool mFallbackNearest = false;
};

/**
 * @brief Build the sparse CSR prolongation matrix `P` (`n_fine x n_coarse`), or
 * restriction `R = P^T` (`n_coarse x n_fine`) when `rOptions.mTranspose` is true.
 *
 * Returns `MeshOperatorResult` where `mNumDegenerate` counts fine points that
 * could not be located (0 when all points are located or `mFallbackNearest` is true).
 */
MESHIOPLUSPLUS_API MeshOperatorResult transfer_operator(const Mesh& rCoarse, const Mesh& rFine,
                                                        const TransferOperatorOptions& rOptions = {});

}  // namespace meshioplusplus
