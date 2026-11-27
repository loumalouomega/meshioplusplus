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
 * @file operations/ip_ops.hpp
 * @brief Integration-point data operations: `ip_to_points`, `ip_to_cells`,
 *        and `points_to_ip` (roadmap 8.3).
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Options for `ip_to_points`.
struct IpToPointsOptions {
    /// Specific integration-point `cell_data` array names to extrapolate.
    /// Empty means all `cell_data` arrays that carry an integration-point layout.
    std::vector<std::string> mArrays;
    /// When true, drop the converted `cell_data` arrays and their `meshio:ip:` metadata.
    bool mDropSource = false;
};

/// Options for `ip_to_cells`.
struct IpToCellsOptions {
    /// Specific integration-point `cell_data` array names to average per cell.
    /// Empty means all `cell_data` arrays that carry an integration-point layout.
    std::vector<std::string> mArrays;
    /// When true, drop the `meshio:ip:` metadata after replacing the cell array.
    bool mDropMetadata = true;
};

/// Options for `points_to_ip`.
struct PointsToIpOptions {
    /// Specific `point_data` array names to sample at integration points.
    /// Empty means all numeric `point_data` arrays.
    std::vector<std::string> mArrays;
    /// Quadrature rule id (`0` = element-nodal, `>= 1` = Gauss exactness order).
    std::int64_t mRuleId = 2;
    /// When true, drop the source `point_data` arrays after sampling.
    bool mDropSource = false;
};

/// Result returned by `ip_to_points`, `ip_to_cells`, and `points_to_ip`.
struct IpOpResult {
    Mesh mMesh;
    std::vector<std::string> mConverted;
    std::vector<std::string> mWarnings;
};

/**
 * @brief Extrapolate integration-point `cell_data` arrays to `point_data` using
 *        element shape functions (least-squares when `n_ip != n_nodes`) and
 *        nodal averaging across incident cells.
 */
MESHIOPLUSPLUS_API IpOpResult ip_to_points(const Mesh& rMesh,
                                           const IpToPointsOptions& rOptions = {});

/**
 * @brief Average integration-point `cell_data` arrays to one value per cell
 *        using quadrature weights when known (or uniform weights otherwise).
 */
MESHIOPLUSPLUS_API IpOpResult ip_to_cells(const Mesh& rMesh, const IpToCellsOptions& rOptions = {});

/**
 * @brief Evaluate `point_data` arrays at element integration points and attach
 *        them as `cell_data` with `field_data["meshio:ip:<name>"]`.
 */
MESHIOPLUSPLUS_API IpOpResult points_to_ip(const Mesh& rMesh,
                                           const PointsToIpOptions& rOptions = {});

}  // namespace meshioplusplus
