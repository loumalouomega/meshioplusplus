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
 * @file detail/tri_pairs.hpp
 * @brief The triangle pairs of one or two soups whose bounding boxes
 * overlap, found over a bucket grid: the candidate search of
 * `find_intersections`.
 *
 * The grid is `check_conformity`'s (a `SpatialGrid` over the boxes, a cell
 * as large as the triangles' mean extent and never finer than 1/256 of the
 * diagonal), so a large triangle among small ones spans a bounded number of
 * buckets. Each triangle's candidates are collected, sorted and made unique
 * on their own, so which bucket a pair is met in is unobservable and the
 * candidates are the same on every backend and thread count.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

// Project includes
#include "meshioplusplus/detail/surface_distance.hpp"

namespace meshioplusplus::detail {

/// Called once per triangle `I` of the first soup with its candidates `rJs`,
/// ascending: the triangles of the second soup (of the first, above `I`, when
/// there is no second) whose boxes meet `I`'s, touching included.
using TriPairVisitor = std::function<void(std::size_t I, const std::vector<std::int64_t>& rJs)>;

/// Visit the box-overlapping pairs of @p rA with itself (@p pB null), each
/// once at its lower triangle, or of @p rA against @p *pB. The visits run in
/// parallel: @p rVisit writes to per-`I` storage only.
void soup_box_pairs(const TriangleSoup& rA, const TriangleSoup* pB, const TriPairVisitor& rVisit);

}  // namespace meshioplusplus::detail
