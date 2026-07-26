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
 * @file operations/data_condition.hpp
 * @brief Value conditioning for data arrays: clamp, normalize, standardize.
 *
 * - `Clamp` maps each value to `min(max(x, lo), hi)`.
 * - `Normalize` maps the array's own `[min, max]` linearly onto `[lo, hi]`
 *   (`[0, 1]` by default).
 * - `Standardize` maps to zero mean and unit standard deviation (population
 *   standard deviation, divisor N).
 *
 * `ConditionScope::Component` (the default) treats each trailing component
 * independently, so each column of a vector array is conditioned on its own
 * statistics. `ConditionScope::Magnitude` instead computes statistics over each
 * row's Euclidean magnitude and rescales whole rows, preserving direction.
 *
 * For a `cell_data` array the statistics are computed **jointly over all cell
 * blocks** and one transform is then applied to every block — the only
 * consistent meaning for normalizing an array that happens to be split across
 * blocks.
 *
 * Degenerate reductions are handled rather than producing NaN: `Normalize` on a
 * constant array fills `lo` and warns, and `Standardize` with zero standard
 * deviation fills `0` and warns.
 *
 * Non-finite values follow the policy documented in
 * `operations/data_common.hpp` — always excluded from the statistics, with
 * `NanPolicy` deciding what reaches the output.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format — it is not
 * in the format registry.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {

/// Which conditioning transform to apply.
enum class ConditionMode {
    Clamp,        ///< `x -> min(max(x, lo), hi)`.
    Normalize,    ///< Affine map from the array's `[min, max]` onto `[lo, hi]`.
    Standardize,  ///< `(x - mean) / stddev`.
};

/// Whether components are conditioned independently or by row magnitude.
enum class ConditionScope {
    Component,  ///< Each trailing component independently (default).
    Magnitude,  ///< Statistics over row magnitude; rows rescaled, direction kept.
};

/// Options for `data_condition`.
struct DataConditionOptions {
    /// Which data map the arrays live in.
    DataLocation location = DataLocation::Point;
    /// Names to condition; empty means every array at the location.
    std::vector<std::string> names;
    /// The transform to apply.
    ConditionMode mode = ConditionMode::Clamp;
    /// Component-wise or by row magnitude.
    ConditionScope scope = ConditionScope::Component;
    /// Clamp lower bound, or Normalize target lower bound.
    double lo = 0.0;
    /// Clamp upper bound, or Normalize target upper bound.
    double hi = 1.0;
    /// What reaches the output for non-finite values.
    NanPolicy nan_policy = NanPolicy::Ignore;
    /// Replacement used when `nan_policy` is `Replace`.
    double nan_replacement = 0.0;
    /// Empty replaces the array in place; otherwise the result is stored as
    /// `name + suffix` and the original is left alone.
    std::string suffix;
    /// `Clamp` only: keep the input dtype instead of producing `Float64`.
    /// `Normalize` and `Standardize` always produce `Float64`.
    bool preserve_dtype = true;
};

/**
 * @brief Conditions the values of the selected data arrays.
 * @param rMesh the source mesh (unmodified).
 * @param rOpts which arrays, which transform, and how to store the result.
 * @return a new mesh with the conditioned arrays; geometry is untouched.
 * @throws std::invalid_argument on an unknown name, on a `cell_data` array
 *         whose block count disagrees with the mesh, on `lo > hi`, or under
 *         `NanPolicy::Fail`.
 */
MESHIOPLUSPLUS_API Mesh data_condition(const Mesh& rMesh, const DataConditionOptions& rOpts);

/**
 * @brief Parses a conditioning-mode name.
 * @param rName one of `clamp`, `normalize`/`rescale`, `standardize`/`zscore`.
 * @return the matching `ConditionMode`.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API ConditionMode condition_mode_from_name(const std::string& rName);

/**
 * @brief Parses a conditioning-scope name.
 * @param rName one of `component`, `magnitude`.
 * @return the matching `ConditionScope`.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API ConditionScope condition_scope_from_name(const std::string& rName);

}  // namespace meshioplusplus
