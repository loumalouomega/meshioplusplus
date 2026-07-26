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
 * @file operations/data_info.hpp
 * @brief Read-only per-array summary of everything a mesh's data maps carry.
 *
 * `data_info` is the *data* view that complements the topological `info` verb
 * and the geometric `compute_stats`: for every `point_data`, `cell_data` and
 * `field_data` array it reports the location, dtype, shape and component count,
 * the number of entries, min / max / mean (whole-array and per component), and
 * the counts of NaN and infinite values. The mesh is not modified.
 *
 * Unlike the other data operations this one **never throws** on non-finite
 * values — it counts them, which is the whole point of the report.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format — it is not
 * in the format registry.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {

/// Read-only summary of one data array (see `data_info`).
struct DataArrayInfo {
    DataLocation mLocation = DataLocation::Point;  ///< Which data map it lives in.
    std::string mName;                             ///< The array name.
    DType mDtype = DType::Float64;                 ///< Its dtype, as stored.
    /// Shape as stored. For `cell_data` this is block 0's shape; see `mNumBlocks`.
    std::vector<std::size_t> mShape;
    /// `cell_data`: the number of cell blocks. 1 for the other locations.
    std::int64_t mNumBlocks = 1;
    /// Rows: points, cells summed over all blocks, or the field length.
    std::int64_t mNumEntries = 0;
    /// Product of the trailing dimensions (1 for a scalar array).
    std::int64_t mNumComponents = 1;
    /// `mNumEntries * mNumComponents`.
    std::int64_t mNumValues = 0;
    double mMin = 0.0;                      ///< Smallest finite value, or NaN when there are none.
    double mMax = 0.0;                      ///< Largest finite value, or NaN when there are none.
    double mMean = 0.0;                     ///< Mean of the finite values, or NaN.
    std::vector<double> mMinPerComponent;   ///< Size `mNumComponents`.
    std::vector<double> mMaxPerComponent;   ///< Size `mNumComponents`.
    std::vector<double> mMeanPerComponent;  ///< Size `mNumComponents`.
    std::int64_t mNumNan = 0;               ///< Count of NaN values.
    std::int64_t mNumInf = 0;               ///< Count of +/-inf values.
    std::int64_t mNumFinite = 0;            ///< Count of finite values.
    /// `cell_data` whose blocks disagree in component count. Reported rather
    /// than thrown, since `data_info` is a diagnostic.
    bool mInconsistentBlocks = false;
};

/// The full read-only report (see `data_info`).
struct DataInfoReport {
    /// `point_data`, then `cell_data`, then `field_data`; each group in the
    /// sorted-name order the uniform API guarantees.
    std::vector<DataArrayInfo> mArrays;
    std::int64_t mNumPointData = 0;  ///< Number of `point_data` arrays.
    std::int64_t mNumCellData = 0;   ///< Number of `cell_data` arrays.
    std::int64_t mNumFieldData = 0;  ///< Number of `field_data` arrays.
};

/**
 * @brief Summarizes every data array a mesh carries (read-only).
 * @param rMesh the mesh to inspect (unmodified).
 * @return the populated report.
 */
MESHIOPLUSPLUS_API DataInfoReport data_info(const Mesh& rMesh);

/**
 * @brief Summarizes one named data array.
 * @param rMesh the mesh to inspect (unmodified).
 * @param Location which data map the array lives in.
 * @param rName the array name.
 * @return the populated summary.
 * @throws std::invalid_argument, listing the available keys, if absent.
 */
MESHIOPLUSPLUS_API DataArrayInfo data_array_info(const Mesh& rMesh, DataLocation Location, const std::string& rName);

}  // namespace meshioplusplus
