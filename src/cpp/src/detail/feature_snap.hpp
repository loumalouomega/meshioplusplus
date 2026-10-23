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
 * @file detail/feature_snap.hpp
 * @brief The sharp features of a skin as corners and chains, and snapping a
 * volume mesh's boundary onto them, for `octree_mesh`.
 *
 * **The feature graph.** A skin's sharp and non-manifold edges form a graph; a
 * vertex where it branches, ends, or turns by more than the corner angle is a
 * *corner*, and the paths between corners are *chains* (a closed loop with no
 * corner is cut at its lowest-numbered vertex and half way round).
 *
 * **Snapping.** Each corner, in ascending order, moves the nearest boundary
 * vertex within reach (1.5 local sizes) that can take it onto itself, exactly,
 * and pins it. Each chain whose ends are snapped is then traced on the
 * boundary's edges by a shortest path that stays within a tube around the
 * polyline (Dijkstra, ties by vertex id, edge cost growing with the distance
 * from the polyline), and the path's inner vertices are moved onto the
 * polyline and pinned. Every move goes through one guard: no tetrahedron
 * around the vertex may lose its volume or drop below a two-degree dihedral
 * angle it did not already have; a refused chain vertex tries a half and a
 * quarter of the move, a refused corner the next candidate. A corner or chain
 * that cannot be recovered is counted, never silently dropped.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstdint>
#include <functional>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

/// The corners and chains of a skin's feature graph.
struct FeatureGraph {
    /// Corner points (skin point ids), ascending.
    std::vector<std::int64_t> mCorners;
    /// Chains of skin point ids, each from a corner to a corner.
    std::vector<std::vector<std::int64_t>> mChains;
};

/// The feature graph of the edges @p rEdges (pairs of skin point ids) over
/// the skin points @p rPoints.
FeatureGraph feature_graph(const std::vector<Vec3>& rPoints,
                           const std::vector<std::array<std::int64_t, 2>>& rEdges,
                           double CornerAngleDeg);

/// What `feature_snap` recovered.
struct FeatureSnapReport {
    std::int64_t mNumCorners = 0;
    std::int64_t mNumCornersSnapped = 0;
    std::int64_t mNumChains = 0;
    std::int64_t mNumChainsMissed = 0;
    std::int64_t mNumMoved = 0;  ///< boundary vertices moved
};

/**
 * @brief Snap the boundary of the tetrahedra @p rTets (faces of the tetrahedra
 * with `rInside[t]` that no other such tetrahedron shares) onto the features.
 * @param rPoints the mesh points, moved in place.
 * @param rPinned per point; snapped points are set to 1 and pinned points are
 *        never moved.
 * @param rFeature per point, set to 1 on a corner and 2 on a chain.
 * @param rLocalSize the local mesh size at a point in space.
 */
FeatureSnapReport feature_snap(std::vector<Vec3>& rPoints,
                               const std::vector<std::array<std::int64_t, 4>>& rTets,
                               const std::vector<std::uint8_t>& rInside, const FeatureGraph& rGraph,
                               const std::vector<Vec3>& rSkinPoints,
                               const std::function<double(const Vec3&)>& rLocalSize,
                               double MinVolume, std::vector<std::uint8_t>& rPinned,
                               std::vector<std::uint8_t>& rFeature);

}  // namespace detail
}  // namespace meshioplusplus
