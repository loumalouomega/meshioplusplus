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
 * @file detail/vtk_write_cell.hpp
 * @brief The VTK cell type and node order a VTK writer writes a block as.
 *
 * A **core-private** header (the `vtk_hex_order.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * `meshio_to_vtk_type` and `meshio_to_vtk_order` (`vtk_common.hpp`) cover the
 * types VTK has an id for. The fixed-order types beyond those (`triangle10`,
 * `quad16`, `tetra20`, `wedge40`, `hexahedron64`, ...) are written as VTK
 * Lagrange cells, their Gmsh node order permuted into VTK's
 * (`lagrange::fixed_lagrange`, `lagrange::gmsh_to_vtk_order`), and read back
 * as the `VTK_LAGRANGE_*` type: the twin of `_vtk_common.FIXED_LAGRANGE`. The
 * VTU, legacy VTK, VTKHDF and VTX writers ask here.
 */

// System includes
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "../formats/lagrange_common.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {
namespace detail {

/// The VTK cell type id `rType` is written as, or none.
inline std::optional<std::int64_t> vtk_type_for_write(const std::string& rType) {
    const auto& types = meshio_to_vtk_type();
    if (const auto it = types.find(rType); it != types.end())
        return it->second;
    lagrange::FixedLagrange fixed;
    if (lagrange::fixed_lagrange(rType, fixed))
        return fixed.mVtkType;
    return std::nullopt;
}

/// The node permutation `rType` is written with: `order[j]` is the meshio++
/// node at VTK position `j`; empty for none.
inline std::vector<int> vtk_order_for_write(const std::string& rType) {
    lagrange::FixedLagrange fixed;
    if (meshio_to_vtk_type().count(rType) == 0 && lagrange::fixed_lagrange(rType, fixed)) {
        std::vector<int> order = lagrange::gmsh_to_vtk_order(fixed.mShape, fixed.mOrder);
        for (std::size_t j = 0; j < order.size(); ++j)
            if (order[j] != static_cast<int>(j))
                return order;
        return {};
    }
    return meshio_to_vtk_order(rType);
}

}  // namespace detail
}  // namespace meshioplusplus
