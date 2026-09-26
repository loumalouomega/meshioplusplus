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
 * @file operations/neighbors.hpp
 * @brief Exact radius and k-nearest neighbour search over a point cloud.
 *
 * The search behind `meshioplusplus.proximity_graph` (roadmap §3): the Python
 * layer keeps the input handling and the graph assembly (sorting,
 * symmetrising, deduplicating -- the Non-goals keep graph construction in
 * Python) and hands the core the pair search, which dominated it (200k points
 * took 6.4 s for a radius graph and 23 s for k = 16 in numpy on one core).
 *
 * The answer is exactly numpy's: a pair's squared distance is `dx*dx + dy*dy
 * + dz*dz` in that order, a radius is inclusive, k-nearest ties go to the
 * lower neighbour index, and a periodic displacement is reduced to its
 * minimum image with round-half-to-even (`np.round`). The lattice is only an
 * accelerator: any cell size gives the same pairs.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// Which neighbours `neighbor_pairs` returns.
enum class NeighborMethod : std::uint8_t {
    Radius,    ///< every pair within `mRadius` (inclusive), once, as (i, j) with i < j
    KNearest,  ///< every point's `mK` nearest, as directed (point, neighbour) pairs
};

/// Options for `neighbor_pairs`.
struct NeighborOptions {
    NeighborMethod mMethod = NeighborMethod::Radius;
    /// The cutoff for `Radius`; must be positive.
    double mRadius = 0.0;
    /// The neighbour count for `KNearest`; clamped to `N - 1`.
    std::int64_t mK = 0;
    /// The periodic box, one side per coordinate, or empty for none. The
    /// points must already lie in `[0, side)` on every axis (the caller wraps
    /// them), and a radius may not exceed half the smallest side.
    std::vector<double> mBox;
    /// The lattice cell side, or 0 for the automatic choice. It changes how
    /// many candidates are examined, never the answer.
    double mCellSize = 0.0;
};

/// The pairs `neighbor_pairs` found, as two Int64 arrays of equal length.
struct NeighborPairs {
    NDArray mSource;
    NDArray mTarget;
};

/**
 * @brief The radius or k-nearest neighbour pairs of a point cloud.
 * @param rPoints Float64 `(N, d)` coordinates, `d` in 1..3, all finite.
 * @param rOptions The method and its parameter.
 * @return For `Radius`, each pair once with `source < target`; for
 *         `KNearest`, `k` pairs per point in ascending (distance, index)
 *         order. Pairs come grouped by source, ascending.
 * @throws std::invalid_argument on a bad shape, a non-finite coordinate, a
 *         non-positive radius, a bad box, or a radius over half the box.
 */
MESHIOPLUSPLUS_API NeighborPairs neighbor_pairs(const NDArray& rPoints,
                                                const NeighborOptions& rOptions);

}  // namespace meshioplusplus
