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
 * @file detail/bucket_build.hpp
 * @brief Build a frozen `BucketTable` from boxes of grid keys, in parallel.
 *
 * A **core-private** header (the `detail/slot_runs.hpp` precedent). The table
 * is the one serial `SpatialGrid::InsertBox(rLo[i], rHi[i], i)` calls for
 * `i = 0 .. n-1` would leave -- every bucket's ids ascending, the buckets in
 * the order those inserts first met them, the occupied box covering every
 * box's two corners -- in compressed-sparse-row form, with no per-cell node
 * or vector. A point is the box whose corners are its key. The distance
 * query builds its triangle grid with it, and the interpolation family its
 * point and simplex grids (roadmap 5.4.1).
 */

// System includes
#include <vector>

// Project includes
#include "meshioplusplus/detail/bucket_table.hpp"

namespace meshioplusplus {
namespace detail {

/// The cell key of a (z-padded, 3-component) coordinate: `SpatialGrid::KeyOf`.
inline GridKey bucket_key(const double* pCoords, double Cell) {
    return GridKey{grid_quantize(pCoords[0], Cell), grid_quantize(pCoords[1], Cell),
                   grid_quantize(pCoords[2], Cell)};
}

/**
 * @brief The bucket table of `rLo.size()` key boxes (`rLo[i]` to `rHi[i]`,
 * inclusive), id `i` for box `i`. Empty input gives an empty table.
 */
BucketTable build_bucket_table(double Cell, const std::vector<GridKey>& rLo,
                               const std::vector<GridKey>& rHi);

}  // namespace detail
}  // namespace meshioplusplus
