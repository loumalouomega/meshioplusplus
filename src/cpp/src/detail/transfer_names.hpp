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
 * @file detail/transfer_names.hpp
 * @brief Plumbing shared by the data-transfer operations (`interpolate`,
 * `conservative_interpolate`, `merge`): the truncating element store, the
 * unknown-array diagnostic and the `on_conflict` name policy.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): it lives beside
 * the `.cpp` files, no installed header names it, and it adds nothing to the
 * API or the ABI. The operation name is a parameter, so every diagnostic stays
 * byte-identical to the per-operation copies it replaces.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {
namespace detail {

/// Store `v` at flat element index `flat` of `rOut`, casting to `rOut`'s dtype
/// with a plain `static_cast` (truncating, unlike `value_io.hpp`'s rounding,
/// saturating `write_double`).
inline void store_double_truncating(NDArray& rOut, std::size_t flat, double v) {
    switch (rOut.Dtype()) {
        case DType::Float32:
            rOut.As<float>()[flat] = static_cast<float>(v);
            return;
        case DType::Float64:
            rOut.As<double>()[flat] = v;
            return;
        case DType::Int8:
            rOut.As<std::int8_t>()[flat] = static_cast<std::int8_t>(v);
            return;
        case DType::Int16:
            rOut.As<std::int16_t>()[flat] = static_cast<std::int16_t>(v);
            return;
        case DType::Int32:
            rOut.As<std::int32_t>()[flat] = static_cast<std::int32_t>(v);
            return;
        case DType::Int64:
            rOut.As<std::int64_t>()[flat] = static_cast<std::int64_t>(v);
            return;
        case DType::UInt8:
            rOut.As<std::uint8_t>()[flat] = static_cast<std::uint8_t>(v);
            return;
        case DType::UInt16:
            rOut.As<std::uint16_t>()[flat] = static_cast<std::uint16_t>(v);
            return;
        case DType::UInt32:
            rOut.As<std::uint32_t>()[flat] = static_cast<std::uint32_t>(v);
            return;
        case DType::UInt64:
            rOut.As<std::uint64_t>()[flat] = static_cast<std::uint64_t>(v);
            return;
    }
}

/// The standard unknown-key diagnostic of operation @p pOpName, covering both
/// source locations.
inline std::string transfer_unknown_array_message(const char* pOpName, const Mesh& rSource,
                                                  const std::string& rName) {
    auto join = [](const std::vector<std::string>& rNames) -> std::string {
        if (rNames.empty())
            return "none";
        std::string s;
        for (std::size_t i = 0; i < rNames.size(); ++i) {
            if (i > 0)
                s += ", ";
            s += rNames[i];
        }
        return s;
    };
    return std::string("meshio++: ") + pOpName +
           ": no source point_data or cell_data array named '" + rName +
           "' (available point_data: " + join(data_names(rSource, DataLocation::Point)) +
           "; cell_data: " + join(data_names(rSource, DataLocation::Cell)) + ")";
}

/// The output name for a transferred array under the conflict policy, resolved
/// against the *original* target before any sampling happens. `TConflict` is
/// the operation's `Error` / `Overwrite` / `Suffix` enum.
template <class TConflict>
std::string transfer_resolve_name(const char* pOpName, const Mesh& rTarget, bool isPointData,
                                  const std::string& rName, TConflict Conflict) {
    const bool exists = isPointData ? rTarget.HasPointData(rName) : rTarget.HasCellData(rName);
    if (!exists)
        return rName;
    switch (Conflict) {
        case TConflict::Overwrite:
            return rName;
        case TConflict::Suffix: {
            const std::string suffixed = rName + "_interp";
            const bool taken =
                isPointData ? rTarget.HasPointData(suffixed) : rTarget.HasCellData(suffixed);
            if (taken)
                throw std::invalid_argument(std::string("meshio++: ") + pOpName + ": array '" +
                                            rName + "' already exists on the target and so does '" +
                                            suffixed + "' (on_conflict='suffix')");
            return suffixed;
        }
        case TConflict::Error:
        default:
            throw std::invalid_argument(std::string("meshio++: ") + pOpName + ": array '" + rName +
                                        "' already exists on the target (on_conflict='error'; "
                                        "pass 'overwrite' or 'suffix')");
    }
}

/// Parses an `on_conflict` name ('error', 'overwrite' or 'suffix') into
/// operation @p pOpName's `TConflict` enum.
template <class TConflict>
TConflict transfer_conflict_from_name(const char* pOpName, const std::string& rName) {
    if (rName == "error")
        return TConflict::Error;
    if (rName == "overwrite")
        return TConflict::Overwrite;
    if (rName == "suffix")
        return TConflict::Suffix;
    throw std::invalid_argument(std::string("meshio++: ") + pOpName + ": unknown on_conflict '" +
                                rName + "' (expected 'error', 'overwrite' or 'suffix')");
}

}  // namespace detail
}  // namespace meshioplusplus
