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
 * @file detail/cell_components.hpp
 * @brief Node-connected components of a mesh's cells: two cells are in one
 * component when a chain of cells, each sharing a node with the next, joins
 * them.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI. `split`'s
 * `by="component"` pieces and `select_region`'s component criterion share it,
 * so the k-th piece of one is the k-th component of the other.
 *
 * Components are numbered in the order of each one's lowest global
 * (block-major) cell. The union sweep is order-dependent and serial; the label
 * resolution is a scan, so the numbering is the same on every backend and
 * thread count.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// The component of every cell, and how many components there are.
struct CellComponents {
    /// One label per global (block-major) cell, in `[0, mCount)`.
    std::vector<std::int64_t> mLabel;
    /// Number of components.
    std::int64_t mCount = 0;
};

/// Label the node-connected components of @p rMesh's cells.
CellComponents cell_components(const Mesh& rMesh);

/// Label the edge-connected components (cells sharing a facet, over
/// `build_facet_dual_graph` with polygon rings), numbered alike: in the order
/// of each one's lowest global cell. Serial (roadmap 7.8.1).
CellComponents cell_components_by_edge(const Mesh& rMesh);

}  // namespace detail
}  // namespace meshioplusplus
