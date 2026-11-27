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
 * @file operations/fvm_quality.hpp
 * @brief Finite-volume (CFD) interior face-pair mesh quality diagnostics
 * (roadmap 9.11): non-orthogonality (degrees), FVM face skewness, adjacent
 * cell volume growth ratio, and face pyramid-height / concavity checks across
 * every interior owner-neighbour face pair (including polyhedra).
 *
 * Worst-face values are projected onto owner and neighbour cells as
 * `fvm:non_orthogonality`, `fvm:face_skewness`, `fvm:volume_ratio`, and
 * `fvm:min_pyramid_height`, which `check_quality` (`quality_gate.hpp`) can
 * gate on directly (e.g. `"fvm:non_orthogonality <= 70"`).
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
#include "meshioplusplus/operations/quality.hpp"

namespace meshioplusplus {

inline constexpr const char* kFvmNonOrthogonalityName = "fvm:non_orthogonality";
inline constexpr const char* kFvmFaceSkewnessName = "fvm:face_skewness";
inline constexpr const char* kFvmVolumeRatioName = "fvm:volume_ratio";
inline constexpr const char* kFvmMinPyramidHeightName = "fvm:min_pyramid_height";

/// Summary report of `compute_fvm_quality`.
struct FvmQualityReport {
    /// Number of 3-D volume cells evaluated.
    std::int64_t mNumCells = 0;
    /// Total faces found (interior + boundary).
    std::int64_t mNumFaces = 0;
    /// Number of interior owner-neighbour face pairs evaluated.
    std::int64_t mNumInteriorFaces = 0;
    /// Number of interior faces with non-orthogonality > 70 deg.
    std::int64_t mNumSevereNonOrtho = 0;
    /// Number of interior faces with negative owner or neighbour pyramid height.
    std::int64_t mNumNegativePyramids = 0;
    /// Per-interior-face summary statistics keyed by `fvm:*` metric name.
    std::vector<std::pair<std::string, QualityMetricSummary>> mMetrics;
    /// Worst-face values projected onto cells, per block, keyed by `fvm:*` name.
    std::vector<std::pair<std::string, std::vector<NDArray>>> mCellArrays;
};

/**
 * @brief Compute finite-volume face-pair quality metrics on a 3-D volume mesh.
 */
MESHIOPLUSPLUS_API FvmQualityReport compute_fvm_quality(const Mesh& rMesh);

/**
 * @brief Return a copy of @p rMesh with the `fvm:*` cell_data arrays attached.
 */
MESHIOPLUSPLUS_API Mesh attach_fvm_quality(const Mesh& rMesh);

}  // namespace meshioplusplus
