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
 * @file operations/data_common.hpp
 * @brief Shared vocabulary for the *data* operations — the ones that act on a
 * mesh's `point_data` / `cell_data` / `field_data` arrays rather than on its
 * geometry.
 *
 * `DataLocation` and `NanPolicy` appear in the public signature of every data
 * operation, so they live here rather than in `detail/`. The `*_from_name`
 * parsers exist because both CLIs and the WASM binding pass these enums across
 * their boundary as plain strings.
 *
 * ### The non-finite (NaN / inf) policy, documented once
 *
 * Every data operation follows the same rule:
 *
 * - Non-finite values are **always excluded from every reduction** — min, max,
 *   mean, standard deviation, and the accumulators used by point/cell
 *   averaging. An array that is entirely non-finite therefore reduces to NaN.
 * - What reaches the *output* is chosen by `NanPolicy`: `Ignore` (the default)
 *   passes the value through unchanged, `Replace` substitutes the caller's
 *   replacement value, and `Fail` throws.
 * - `data_info` never throws — it *counts* non-finite values instead.
 * - Division by zero inside `data_calc` produces an IEEE infinity or NaN and is
 *   deliberately **not** an error.
 *
 * Everything here is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. These are operations, not file formats — they are
 * not in the format registry.
 */

// System includes
#include <cstddef>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Which of a mesh's three data maps an array lives in.
enum class DataLocation {
    Point,  ///< `point_data` — one row per mesh point.
    Cell,   ///< `cell_data` — one array *per cell block*, one row per cell.
    Field,  ///< `field_data` — mesh-global, no entity correspondence.
};

/// What happens to non-finite values on the way to the output (see the file
/// comment: they are excluded from reductions regardless of this setting).
enum class NanPolicy {
    Ignore,   ///< Pass NaN/inf through to the output unchanged (default).
    Replace,  ///< Substitute the caller's replacement value.
    Fail,     ///< Throw `std::invalid_argument` naming the array and index.
};

/**
 * @brief Parses a location name.
 * @param rName one of `point`/`point_data`/`node`, `cell`/`cell_data`/`element`,
 *        `field`/`field_data` (case-sensitive).
 * @return the matching `DataLocation`.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API DataLocation data_location_from_name(const std::string& rName);

/**
 * @brief The canonical name of a location, as used in diagnostics.
 * @param Location the location to name.
 * @return `"point_data"`, `"cell_data"` or `"field_data"`.
 */
MESHIOPLUSPLUS_API const char* data_location_name(DataLocation Location);

/**
 * @brief Parses a NaN-policy name.
 * @param rName one of `ignore`, `replace`, `fail`.
 * @return the matching `NanPolicy`.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API NanPolicy nan_policy_from_name(const std::string& rName);

/**
 * @brief The names present at @p Location, sorted.
 *
 * Goes through the uniform API's `*DataNames()` accessors, which guarantee
 * sorted order on every backend.
 * @param rMesh the mesh to inspect.
 * @param Location which data map to list.
 * @return the sorted array names.
 */
MESHIOPLUSPLUS_API std::vector<std::string> data_names(const Mesh& rMesh, DataLocation Location);

/**
 * @brief Whether @p rName exists at @p Location in @p rMesh.
 * @param rMesh the mesh to inspect.
 * @param Location which data map to look in.
 * @param rName the array name.
 * @return `true` if present.
 */
MESHIOPLUSPLUS_API bool data_has(const Mesh& rMesh, DataLocation Location, const std::string& rName);

/**
 * @brief Builds the standard "unknown key" diagnostic, listing what *is* there.
 *
 * Produces e.g. `meshio++: no point_data array named 'foo' (available: T, p, u)`
 * or, when the location is empty, `... (the mesh has no point_data)`.
 * @param rMesh the mesh whose available keys are listed.
 * @param Location the location that was searched.
 * @param rName the key that was not found.
 * @return the formatted message.
 */
MESHIOPLUSPLUS_API std::string data_unknown_key_message(const Mesh& rMesh, DataLocation Location,
                                     const std::string& rName);

/**
 * @brief Number of components of an array — the product of its trailing
 * dimensions, with a 1-D (or shapeless) array counting as 1 component.
 * @param rArray the array to measure.
 * @return the component count (at least 1).
 */
MESHIOPLUSPLUS_API std::size_t data_num_components(const NDArray& rArray);

}  // namespace meshioplusplus
