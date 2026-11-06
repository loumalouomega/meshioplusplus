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
 * @file detail/point_field_check.hpp
 * @brief The argument check `gradient` and `hessian` share: a named
 * point_data field, and a result location on the mesh.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <stdexcept>
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {
namespace detail {

/// Throw `std::invalid_argument` (prefixed with `pPrefix`) unless
/// `rArrayName` names point_data of `rMesh` and `Location` is not Field. A
/// cell_data field of that name is called out: it is piecewise constant and
/// has no derivative.
inline void require_point_field(const Mesh& rMesh, const std::string& rArrayName,
                                DataLocation Location, const char* pPrefix) {
    if (rArrayName.empty())
        throw std::invalid_argument(std::string(pPrefix) + "an array name is required");
    if (Location == DataLocation::Field)
        throw std::invalid_argument(std::string(pPrefix) +
                                    "field_data has no location on the mesh; the result must be "
                                    "'point' or 'cell'");
    if (!rMesh.HasPointData(rArrayName)) {
        // A cell_data field is piecewise constant and has no derivative. Name
        // the fix, exactly as isosurface does for the same reason.
        if (rMesh.HasCellData(rArrayName))
            throw std::invalid_argument(
                std::string(pPrefix) + "'" + rArrayName +
                "' is cell_data, which is piecewise constant and has no derivative; convert it "
                "first with cell_data_to_point_data (CLI: meshioplusplus data to-point)");
        throw std::invalid_argument(
            std::string(pPrefix) +
            data_unknown_key_message(rMesh, DataLocation::Point, rArrayName));
    }
}

}  // namespace detail
}  // namespace meshioplusplus
