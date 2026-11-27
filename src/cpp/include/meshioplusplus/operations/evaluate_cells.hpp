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
 * @file operations/evaluate_cells.hpp
 * @brief High-order cell evaluation (`evaluate_cells`) and exact
 * Lagrange <-> Bezier change of basis (`change_basis`) (roadmap 15.3.1, 15.3.2).
 */

// System includes
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// Options for `evaluate_cells`.
struct EvaluateCellsOptions {
    /// Global block-major cell ids to evaluate (`N` entries).
    std::vector<std::int64_t> mCells;
    /// Reference coordinates in `[0, 1]` (flat `3 * N` entries).
    std::vector<double> mReferencePoints;
    /// Point-data field names to evaluate (empty = all numeric point_data).
    std::vector<std::string> mFields;
    /// If true, also compute physical gradients `nabla_x f` via the inverse Jacobian.
    bool mDerivatives = false;
};

/// Result of `evaluate_cells`.
struct EvaluateCellsResult {
    /// Evaluated physical coordinates, Float64 `(N, 3)`.
    NDArray mPoints;
    /// Evaluated field values per requested field name (`(N,)` or `(N, C)`).
    std::vector<std::pair<std::string, NDArray>> mValues;
    /// Evaluated physical gradients per requested field name (`(N, 3)` or `(N, C, 3)`), if requested.
    std::vector<std::pair<std::string, NDArray>> mGradients;
    /// Number of queries successfully evaluated.
    std::int64_t mNumEvaluated = 0;
};

/// Target polynomial basis for `change_basis`.
enum class BasisTarget : std::int32_t {
    Bezier = 0,
    Lagrange = 1,
};

/// Parse `"bezier"` or `"lagrange"`.
MESHIOPLUSPLUS_API BasisTarget basis_target_from_name(const std::string& rName);

/// Options for `change_basis`.
struct ChangeBasisOptions {
    /// Target basis (`Bezier` or `Lagrange`).
    BasisTarget mTarget = BasisTarget::Bezier;
    /// If true, transform numeric `point_data` arrays alongside coordinates.
    bool mTransformFields = true;
};

/// Result of `change_basis`.
struct ChangeBasisResult {
    /// Converted mesh (`VtkBezier*` with `HigherOrderDegrees` or `VtkLagrange*`).
    Mesh mMesh;
    /// Number of cells transformed.
    std::int64_t mNumConverted = 0;
};

/**
 * @brief Evaluate physical positions, point fields, and optional physical
 * gradients at arbitrary reference points inside mesh cells.
 */
MESHIOPLUSPLUS_API EvaluateCellsResult evaluate_cells(const Mesh& rMesh,
                                                      const EvaluateCellsOptions& rOptions);

/**
 * @brief Convert between Lagrange and Bezier cell representations via exact
 * Bernstein-Lagrange transformation matrices.
 */
MESHIOPLUSPLUS_API ChangeBasisResult change_basis(const Mesh& rMesh,
                                                  const ChangeBasisOptions& rOptions = {});

}  // namespace meshioplusplus
