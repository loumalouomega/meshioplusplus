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
 * @file detail/node_order.hpp
 * @brief One registry of node-ordering permutations between file formats and
 * meshio++'s (VTK) cell layouts, keyed by (format, cell type).
 *
 * Every entry holds both directions as *gather* tables, so a caller never has
 * to know whether a table is its own inverse:
 *
 * - `mToMeshio[k]`: the file slot that meshio++ node `k` comes from on read,
 *   `meshio[k] = file[mToMeshio[k]]`;
 * - `mFromMeshio[j]`: the meshio++ node that file slot `j` comes from on write,
 *   `file[j] = meshio[mFromMeshio[j]]`.
 *
 * A type with no entry uses the identity. The formats covered, and where each
 * table was pinned, are listed in `doc/node_ordering.md`.
 */

// System includes
#include <string_view>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {
namespace detail {

/// Both directions of one (format, cell type) node permutation.
struct NodeOrder {
    std::vector<int> mToMeshio;
    std::vector<int> mFromMeshio;
};

/**
 * @brief The node permutation between a format and meshio++ for one cell type.
 * @param Format The format key (`"med"`, `"code_aster"`, `"frd"`, `"unv"`).
 * @param CellType The meshio++ cell type name (e.g. `"hexahedron20"`).
 * @return The permutation, or `nullptr` when the two orders are identical (or
 *         the format is unknown).
 */
MESHIOPLUSPLUS_API const NodeOrder* node_order(std::string_view Format, std::string_view CellType);

/**
 * @brief Every (format, cell type) key in the registry, for self-tests.
 * @return Pairs of format key and cell type name, sorted.
 */
MESHIOPLUSPLUS_API std::vector<std::pair<std::string_view, std::string_view>> node_order_keys();

}  // namespace detail
}  // namespace meshioplusplus
