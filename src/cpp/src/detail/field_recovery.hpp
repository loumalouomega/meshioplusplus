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
 * @file detail/field_recovery.hpp
 * @brief Array-level internal entry points for cell gradient and point/cell
 * averaging (roadmap 5.27).
 *
 * Core-private (not installed; no ABI impact). Allows composite operations
 * such as `estimate_error` (and point-located `gradient`) to run the recovery
 * pipeline directly on arrays without cloning the whole `Mesh` at every stage
 * and without recomputing per-cell measures.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/gradient.hpp"

namespace meshioplusplus {
namespace detail {

struct GradientArrayResult {
    std::vector<NDArray> mBlocks;
    std::int64_t mNumSkipped = 0;
    std::int64_t mNumFallback = 0;
};

/**
 * @brief Compute per-block cell gradient arrays without cloning @p rMesh.
 *
 * Runs the same validation and arithmetic as `gradient()` with
 * `mLocation == DataLocation::Cell` (skipping the output-name collision check
 * so internal callers need not fabricate a dummy key).
 */
GradientArrayResult gradient_array(const Mesh& rMesh, const GradientOptions& rOptions);

/**
 * @brief Average one per-block cell-data array to points without cloning @p rMesh.
 *
 * When @p pRawMeasures is non-null and `rOpts.weight == CellPointWeight::Measure`,
 * populates `*pRawMeasures` with the raw `detail::cell_measure(...)` values
 * before non-positive / non-finite fallback clamping so callers such as
 * `estimate_error` can reuse them bit-for-bit.
 */
NDArray cell_to_point_array(const Mesh& rMesh, const std::vector<NDArray>& rBlocks,
                            const std::string& rName, const DataAverageOptions& rOpts,
                            std::vector<std::vector<double>>* pRawMeasures = nullptr);

/**
 * @brief Average one point-data array to per-block cell arrays without cloning @p rMesh.
 */
std::vector<NDArray> point_to_cell_array(const Mesh& rMesh, const NDArray& rSrc,
                                         const std::string& rName, const DataAverageOptions& rOpts);

}  // namespace detail
}  // namespace meshioplusplus
