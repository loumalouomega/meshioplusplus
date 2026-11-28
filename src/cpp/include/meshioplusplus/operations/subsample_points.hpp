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
 * @file operations/subsample_points.hpp
 * @brief Count-driven point-cloud subsampling (roadmap 9.9.2).
 *
 * `subsample_points(mesh, options)` selects exactly `mCount` points from
 * `mesh` (optionally restricted to an axis-aligned box `mBounds`) in
 * selection order using one of three methods:
 *
 * - `Farthest`: exact greedy farthest-point sampling (`detail/farthest_points.hpp`),
 *   deterministic given `mStart`.
 * - `Grid`: isotropic lattice refinement picking the nearest-to-centre
 *   representative of each occupied cell, followed by farthest-point sampling
 *   over those representatives.
 * - `Random`: uniform random draw without replacement seeded by `mSeed`.
 *
 * The output is a point cloud carrying one `vertex` block, gathered
 * `point_data`, copied `field_data`, and remapped `Point` regions (empty
 * groups kept). Non-vertex cells, `cell_data`, and `Cell`/`Side` regions are
 * dropped with warnings. Setting `mRecordIds = true` attaches
 * `budget:original_point_id`.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The point_data array attached when `mRecordIds` is `true`.
inline constexpr const char* kBudgetOriginalPointIdName = "budget:original_point_id";

/// Selection method for `subsample_points`.
enum class SubsampleMethod : std::int32_t {
    Farthest = 0,
    Grid = 1,
    Random = 2,
};

/// Options for `subsample_points`.
struct SubsamplePointsOptions {
    /// Number of points to select (>= 1, <= number of candidates).
    std::int64_t mCount = 0;
    /// Selection method (`Farthest`, `Grid`, or `Random`).
    SubsampleMethod mMethod = SubsampleMethod::Farthest;
    /// Seed for `Random` and for `mStart < 0`.
    std::int64_t mSeed = 0;
    /// Index of the first selected point for `Farthest`/`Grid`; `< 0` draws from `mSeed`.
    std::int64_t mStart = 0;
    /// When `true`, restrict candidates to `[xlo, ylo, zlo, xhi, yhi, zhi]`.
    bool mHasBounds = false;
    std::array<double, 6> mBounds{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    /// Lattice cell ceiling for `Grid` (default matches Python `DEFAULT_MAX_CELLS`).
    std::int64_t mMaxCells = 16000000;
    /// Attach `budget:original_point_id` to `point_data`.
    bool mRecordIds = false;
};

/// Output of `subsample_points`.
struct SubsamplePointsResult {
    Mesh mMesh;
    /// Selected input point indices in selection order.
    std::vector<std::int64_t> mPointMap;
    std::vector<std::string> mWarnings;
};

/**
 * @brief Select a fixed count of points and return them as a point cloud.
 * @throws std::invalid_argument on invalid count, method, bounds, start index,
 *         or too few distinct positions / lattice cells.
 */
MESHIOPLUSPLUS_API SubsamplePointsResult subsample_points(const Mesh& rMesh,
                                                          const SubsamplePointsOptions& rOptions);

}  // namespace meshioplusplus
