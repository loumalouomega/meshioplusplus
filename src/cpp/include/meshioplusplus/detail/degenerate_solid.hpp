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
 * @file detail/degenerate_solid.hpp
 * @brief Collapse an 8-node brick with repeated nodes into the solid it is.
 *
 * Keyword decks that have only a brick card (LS-DYNA `*ELEMENT_SOLID`, Radioss
 * `/BRICK`) write tetrahedra, pyramids and wedges as bricks whose trailing
 * nodes repeat. `collapse_brick` recognises the patterns both solvers document,
 * most degenerate first, and returns the cell in meshio++'s node order. The
 * Python twin is `_collapse_solid` in `lsdyna/_lsdyna.py`.
 */

// System includes
#include <array>
#include <cstdint>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {
namespace detail {

/// A brick's cell after collapsing repeated nodes.
struct CollapsedBrick {
    const char* mType;                 ///< "tetra", "pyramid", "wedge" or "hexahedron"
    std::vector<std::int64_t> mNodes;  ///< in meshio++ order
};

/**
 * @brief The cell an 8-node brick with repeated nodes stands for.
 *
 * Patterns, checked in this order: `1 2 3 4 4 4 4 4` and `1 2 3 3 4 4 4 4`
 * (tetra), `1 2 3 4 5 5 5 5` (pyramid), `1 2 3 3 4 5 6 6` and `1 2 3 4 5 5 6 6`
 * (wedge), then any one side edge collapsed in both the bottom and the top face,
 * such as Radioss's `1 2 3 1 5 6 7 5` (wedge). Anything else stays a hexahedron.
 * @param rNodes The brick's eight node ids in file order.
 */
MESHIOPLUSPLUS_API CollapsedBrick collapse_brick(const std::array<std::int64_t, 8>& rNodes);

/**
 * @brief The 8 brick nodes a `tetra`, `pyramid` or `wedge` is written as, the
 * inverse of `collapse_brick` (`1 2 3 4 4 4 4 4`, `1 2 3 4 5 5 5 5`, and a wedge
 * `1 3 6 4 2 2 5 5`); any other type's first eight nodes as they are. The
 * LS-DYNA and Radioss writers share it; the Python twin is `_expand_solid` in
 * `lsdyna/_lsdyna.py`. Since v16.17.0.
 * @param Type The meshio++ cell type.
 * @param pRow The cell's nodes, in meshio++ order.
 */
MESHIOPLUSPLUS_API std::array<std::int64_t, 8> expand_brick(std::string_view Type,
                                                            const std::int64_t* pRow);

}  // namespace detail
}  // namespace meshioplusplus
