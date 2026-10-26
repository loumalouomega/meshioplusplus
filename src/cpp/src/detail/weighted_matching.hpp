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
 * @file detail/weighted_matching.hpp
 * @brief A maximum-weight matching of a general graph (roadmap 6.8.1).
 *
 * Edmonds' blossom algorithm in its primal-dual form (J. Edmonds, *Paths,
 * trees, and flowers*, 1965, and *Maximum matching and a polyhedron with
 * 0,1-vertices*, 1965), in the O(n^3) organisation of Z. Galil, *Efficient
 * algorithms for finding maximum matching in graphs*, ACM Computing Surveys
 * 18 (1986): stages, each growing alternating trees from the free vertices
 * over the tight edges, shrinking an odd cycle into a blossom, expanding a
 * blossom whose dual reaches zero, and changing the duals by the least of
 * the four slacks until an augmenting path appears. The weights are
 * integers and are doubled on entry, so every dual and slack stays an
 * integer and the result is exact; vertices, edges and stacks are visited
 * in index order, so it is deterministic.
 *
 * With `MaxCardinality` the matching is the heaviest among those with the
 * most edges (the dual variables are not allowed to stop the search).
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <cstdint>
#include <vector>

namespace meshioplusplus::detail {

struct MatchingEdge {
    std::int64_t mA = 0;
    std::int64_t mB = 0;
    std::int64_t mWeight = 0;  ///< |weight| below 2^60
};

/// Each vertex's mate in a maximum-weight matching, or -1. Self-loops are
/// ignored; parallel edges are allowed.
std::vector<std::int64_t> max_weight_matching(std::int64_t NumVertices,
                                              const std::vector<MatchingEdge>& rEdges,
                                              bool MaxCardinality);

}  // namespace meshioplusplus::detail
