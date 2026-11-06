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

#pragma once

/**
 * @file detail/cell_centroids.hpp
 * @brief Per-cell centroids of a whole mesh, shared by `partition` and
 * `interpolate` (doc/geometry_kernels.md).
 */

// System includes
#include <cstddef>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus::detail {

/**
 * @brief Per-cell centroids in global (block-major) cell order, always three
 * doubles per cell with unused axes zero: the mean of the cell's nodes, where a
 * polyhedron counts each distinct node of its faces once.
 *
 * The accumulation order is fixed (it is pinned by `partition`'s and
 * `interpolate`'s numpy parity tests, and by `_partition.py`'s `_centroids`).
 * @param rMesh The mesh.
 * @param TotalCells The mesh's total cell count over every block.
 */
std::vector<double> cell_centroids(const Mesh& rMesh, std::size_t TotalCells);

}  // namespace meshioplusplus::detail
