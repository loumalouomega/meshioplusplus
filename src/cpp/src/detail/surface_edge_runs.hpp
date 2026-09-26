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
 * @file detail/surface_edge_runs.hpp
 * @brief A soup's (triangle, corner) edge records grouped by undirected edge,
 * computed once for the two tables built from them.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * `build_surface_edges` (the watertightness report) and `build_distance_query`
 * (the edge normals) group the same keys -- corner `3 t + i`'s edge from
 * corner `i` to corner `i + 1`, as (low, high) -- the same way. An operation
 * that needs both (`sample_distance`, `distance_to_surface`, `shrinkwrap`,
 * `remesh_volume`) groups them once and hands the runs to both, with results
 * identical to the one-argument forms. Roadmap §4, "Build once".
 */

// System includes
#include <vector>

// Project includes
#include "meshioplusplus/detail/surface_distance.hpp"
#include "slot_runs.hpp"

namespace meshioplusplus {
namespace detail {

/// Every corner's undirected edge key and the runs of equal keys, in
/// ascending key order (bucketed by the lower endpoint).
struct SurfaceEdgeRuns {
    std::vector<SurfaceEdgeKey> mKeys;  ///< per corner `3 t + i`
    SlotRuns mRuns;
};

/// Group @p rSoup's corner edges.
SurfaceEdgeRuns surface_edge_runs(const TriangleSoup& rSoup);

/// `build_surface_edges(rSoup)` from runs `surface_edge_runs(rSoup)` built.
SurfaceEdgeMap build_surface_edges_from_runs(const TriangleSoup& rSoup,
                                             const SurfaceEdgeRuns& rRuns);

/// `build_distance_query(rSoup, rOptions)` from runs `surface_edge_runs(rSoup)`
/// built.
DistanceQuery build_distance_query_from_runs(const TriangleSoup& rSoup,
                                             const SurfaceDistanceOptions& rOptions,
                                             const SurfaceEdgeRuns& rRuns);

}  // namespace detail
}  // namespace meshioplusplus
