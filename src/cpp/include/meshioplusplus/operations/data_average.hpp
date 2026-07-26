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
 * @file operations/data_average.hpp
 * @brief Moves data between the point and cell locations by averaging, without
 * touching the geometry.
 *
 * `point_data_to_cell_data` gives each cell the mean of the values at its own
 * nodes. `cell_data_to_point_data` gives each point the mean of the values on
 * the cells incident to it, optionally weighted by each cell's |measure| (area
 * for a 2D cell, volume for a 3D one) so that large cells count for more.
 *
 * Both work component-wise, so scalar, vector and tensor arrays are all handled
 * the same way. Because a mean is not an integer, the output is **always**
 * `Float64` regardless of the input dtype.
 *
 * Non-finite values follow the policy documented in
 * `operations/data_common.hpp`: they never contribute to an average, and a
 * point or cell with no finite contribution at all yields NaN.
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

/// Weighting for the cell_data -> point_data direction.
enum class CellPointWeight {
    Uniform,  ///< Every incident cell contributes equally.
    Measure,  ///< Weight by |cell measure| (length / area / volume by dimension).
};

/// Options for both averaging directions.
struct DataAverageOptions {
    /// Names to convert; empty means every array at the source location.
    std::vector<std::string> names;
    /// Weighting, cell -> point only; ignored in the other direction.
    CellPointWeight weight = CellPointWeight::Uniform;
    /// Output name is `prefix + input name` (empty keeps the same name).
    std::string prefix;
    /// Output name is `input name + suffix` (empty keeps the same name).
    std::string suffix;
    /// When false, an output name that already exists at the target location
    /// is an error rather than being overwritten.
    bool overwrite = true;
    /// What reaches the output for non-finite results.
    NanPolicy nan_policy = NanPolicy::Ignore;
    /// Replacement used when `nan_policy` is `Replace`.
    double nan_replacement = 0.0;
};

/**
 * @brief Averages `point_data` onto the cells: each cell's value is the mean
 * over its own nodes.
 *
 * Ragged polygon blocks average over the row's nodes; polyhedron blocks average
 * over the *distinct* nodes across all their faces.
 * @param rMesh the source mesh (unmodified).
 * @param rOpts which arrays to convert and how to name the results.
 * @return a new mesh carrying the produced `cell_data` alongside everything the
 *         input already had.
 * @throws std::invalid_argument on an unknown name, on an output-name collision
 *         when `overwrite` is false, or under `NanPolicy::Fail`.
 */
MESHIOPLUSPLUS_API Mesh point_data_to_cell_data(const Mesh& rMesh, const DataAverageOptions& rOpts = {});

/**
 * @brief Averages `cell_data` onto the points: each point's value is the
 * (optionally measure-weighted) mean over the cells incident to it.
 *
 * A point touched by no cell — or by no cell with a finite value — gets NaN.
 * Cells whose measure cannot be computed (ragged and polyhedron blocks) fall
 * back to a unit weight under `CellPointWeight::Measure`, with one warning per
 * call.
 * @param rMesh the source mesh (unmodified).
 * @param rOpts which arrays to convert, the weighting, and result naming.
 * @return a new mesh carrying the produced `point_data` alongside everything
 *         the input already had.
 * @throws std::invalid_argument on an unknown name, on a `cell_data` array
 *         whose block count disagrees with the mesh, on an output-name
 *         collision when `overwrite` is false, or under `NanPolicy::Fail`.
 */
MESHIOPLUSPLUS_API Mesh cell_data_to_point_data(const Mesh& rMesh, const DataAverageOptions& rOpts = {});

/**
 * @brief Parses a weighting name.
 * @param rName one of `uniform`/`simple`, or `measure`/`area`/`volume`.
 * @return the matching `CellPointWeight`.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API CellPointWeight cell_point_weight_from_name(const std::string& rName);

}  // namespace meshioplusplus
