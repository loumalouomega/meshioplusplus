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
 * @file operations/quality_gate.hpp
 * @brief Pass/fail thresholds over `compute_quality`'s metrics: the check a CI
 * job over meshes scripts.
 *
 * `check_quality(mesh, {thresholds})` scores every cell with `compute_quality`
 * and tests each threshold against the **per-cell** values (the report's
 * histograms span each metric's own data range, so they cannot answer "how
 * many cells fall below 0.2"). A threshold bounds one metric from below
 * (`>=`), above (`<=`) or both, and may allow a fraction of the cells it
 * applies to to violate it. Cells where the metric does not apply (NaN) are
 * not evaluated. Inverted and degenerate cells are gated by count.
 *
 * ### The specification text
 *
 * Every surface that cannot pass a structure (the CLIs, the C ABI, the
 * pipeline) spells thresholds as `parse_quality_thresholds` reads them:
 * clauses separated by `;`, `,` or a newline (`#` starts a comment to the end
 * of its line, so a gate file is one clause per line), each `METRIC >= VALUE` or
 * `METRIC <= VALUE`, optionally followed by `@ FRACTION` -- a number in
 * `[0, 1]` or a percentage (`@1%`). `METRIC` is a `compute_quality` metric
 * name, with or without its `quality:` prefix:
 *
 *     scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%
 */

// System includes
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/quality.hpp"

namespace meshioplusplus {

/// One bound on one metric.
struct QualityThreshold {
    /// A `compute_quality` metric, with or without the `quality:` prefix.
    std::string mMetric;
    /// The smallest acceptable value, or NaN for none.
    double mMin = std::numeric_limits<double>::quiet_NaN();
    /// The largest acceptable value, or NaN for none.
    double mMax = std::numeric_limits<double>::quiet_NaN();
    /// The fraction of evaluated cells allowed to violate the bound, in `[0, 1]`.
    double mMaxFraction = 0.0;
};

/// What `check_quality` checks.
struct QualityGateOptions {
    std::vector<QualityThreshold> mThresholds;
    /// The most inverted cells allowed; negative disables the check.
    std::int64_t mMaxInverted = 0;
    /// The most degenerate cells allowed; negative disables the check.
    std::int64_t mMaxDegenerate = 0;
};

/// The outcome of one check.
struct QualityCheck {
    /// A readable spelling: `"scaled_jacobian >= 0.2"`, `"inverted <= 0"`.
    std::string mName;
    /// The full metric name (`"quality:scaled_jacobian"`), or `"inverted"` /
    /// `"degenerate"` for the count checks.
    std::string mMetric;
    double mMin = std::numeric_limits<double>::quiet_NaN();
    double mMax = std::numeric_limits<double>::quiet_NaN();
    double mMaxFraction = 0.0;
    /// Cells the check applied to (finite metric value), and how many violate it.
    std::int64_t mEvaluated = 0;
    std::int64_t mViolations = 0;
    /// `mViolations / mEvaluated`, 0 when nothing was evaluated.
    double mFraction = 0.0;
    /// The evaluated value closest to (or furthest beyond) the bound, and its
    /// global (block-major) cell; NaN / -1 when nothing was evaluated.
    double mWorst = std::numeric_limits<double>::quiet_NaN();
    std::int64_t mWorstCell = -1;
    bool mPassed = true;
};

/// Every check, and whether all passed.
struct QualityGateResult {
    bool mPassed = true;
    std::vector<QualityCheck> mChecks;
    /// The `compute_quality` report the checks were evaluated on.
    QualityReport mReport;
};

/**
 * @brief Parse the threshold specification text (see the file comment).
 * @throws std::invalid_argument on a malformed clause or an unknown metric.
 */
MESHIOPLUSPLUS_API std::vector<QualityThreshold> parse_quality_thresholds(const std::string& rText);

/**
 * @brief Score @p rMesh and test every threshold and count limit.
 * @throws std::invalid_argument on an unknown metric, a threshold with neither
 *         bound, or a fraction outside `[0, 1]`.
 */
MESHIOPLUSPLUS_API QualityGateResult check_quality(const Mesh& rMesh,
                                                   const QualityGateOptions& rOptions = {});

/// A multi-line, human-readable summary of @p rResult (one line per check).
MESHIOPLUSPLUS_API std::string quality_gate_summary(const QualityGateResult& rResult);

}  // namespace meshioplusplus
