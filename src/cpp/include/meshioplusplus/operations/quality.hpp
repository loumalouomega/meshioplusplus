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
 * @file operations/quality.hpp
 * @brief Per-cell mesh quality metrics + inverted-element detection, part of
 * the dependency-free mesh-operations layer.
 *
 * `compute_quality` scores every cell of every block on a fixed set of
 * geometric metrics (area/volume, scaled Jacobian, aspect ratio, skewness,
 * interior/dihedral angles, warpage) evaluated on the cell's *corner* nodes
 * (quadratic variants reduce to their linear parent), and aggregates a
 * per-metric summary plus counts of inverted (negative signed volume) and
 * degenerate (near-zero) cells. `attach_quality` returns a copy of the mesh
 * with every metric attached as `cell_data` (name `"quality:<metric>"`), ready
 * to write out and visualise. All metrics use only standard C++ math and the
 * uniform mesh API, so they run under every mesh backend.
 */

// System includes
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Aggregate statistics for one quality metric over all cells that have
 * a finite value for it.
 */
struct QualityMetricSummary {
    /// Number of histogram bins spanning `[mMin, mMax]`.
    static constexpr int K = 10;
    double mMin = std::numeric_limits<double>::quiet_NaN();
    double mMax = std::numeric_limits<double>::quiet_NaN();
    double mMean = std::numeric_limits<double>::quiet_NaN();
    /// Cells contributing a finite value.
    std::int64_t mCount = 0;
    /// Counts per bin over `[mMin, mMax]`; all zero when `mCount == 0`.
    std::array<std::int64_t, K> mHistogram{};
};

/**
 * @brief The result of `compute_quality`: per-metric summaries, per-cell metric
 * arrays (ready to attach as `cell_data`), and inverted/degenerate counts.
 */
struct QualityReport {
    /// Total cells scored (across all blocks, including unsupported ones).
    std::int64_t mNumCells = 0;
    /// Cells with negative signed volume/area.
    std::int64_t mNumInverted = 0;
    /// Cells that are degenerate (near-zero volume/edge/Jacobian).
    std::int64_t mNumDegenerate = 0;
    /// Metric name -> summary, in the fixed metric order (see `quality.cpp`).
    std::vector<std::pair<std::string, QualityMetricSummary>> mMetrics;
    /// Metric name -> one `Float64` array of shape `(num_cells_in_block, 1)`
    /// per cell block (NaN where the metric does not apply), in the same order.
    std::vector<std::pair<std::string, std::vector<NDArray>>> mCellArrays;
};

/**
 * @brief Compute per-cell quality metrics and their aggregate report.
 * @param rMesh the mesh to evaluate
 * @return the filled `QualityReport`
 */
MESHIOPLUSPLUS_API QualityReport compute_quality(const Mesh& rMesh);

/**
 * @brief A copy of @p rMesh with every quality metric attached as `cell_data`.
 *
 * Each metric becomes a `cell_data` entry named `"quality:<metric>"` with one
 * array per cell block (NaN where the metric does not apply to that block's
 * cell type).
 * @param rMesh the mesh to score
 * @return a copy carrying the metric arrays
 */
MESHIOPLUSPLUS_API Mesh attach_quality(const Mesh& rMesh);

}  // namespace meshioplusplus
