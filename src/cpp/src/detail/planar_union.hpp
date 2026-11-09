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
 * @file detail/planar_union.hpp
 * @brief The region a soup of planar triangles covers, minus the region
 * another covers, as triangles and boundary loops (roadmap 7.6.3: the
 * silhouette of a surface, and its offset).
 *
 *  - **Candidate edges.** An edge two triangles of one soup share, lying on
 *    either side of it, is inside that soup's region; every other edge (on
 *    one triangle, on two on the same side, on three or more) may be on its
 *    boundary. For the projection of a surface these are its boundary edges
 *    and the edges where the surface turns from facing the plane to facing
 *    away: the contour generators. Triangles of zero area are left out.
 *  - **Crossings.** The candidate edges of both soups are split where they
 *    meet (`cdt2d_find_crossings`, exact): at the endpoint one touches the
 *    other at, at the endpoints of a collinear overlap, and at the crossing
 *    point of a proper crossing (computed in double from the exact
 *    orientations), until none is left.
 *  - **Coverage.** The split edges are the constraints of a planar CDT
 *    (`detail/cdt2d.hpp`). Coverage is constant across an unconstrained edge,
 *    so each component of triangles joined by unconstrained edges is decided
 *    once: outside the hull when it reaches the super triangle, else by its
 *    largest triangle's centroid, tested exactly (`orient2d`) against the
 *    triangles of each soup found through a uniform bucket grid. A component
 *    is covered when some triangle of the first soup holds the point and none
 *    of the second does.
 *  - **Loops.** The covered triangles' edges whose other side is not covered,
 *    walked with the covered side on the left, turning round each vertex
 *    through the covered triangles only (so two loops touching at a point
 *    stay two loops): outer boundaries counter-clockwise, holes clockwise.
 *    Each loop starts at its lowest point id, and loops are ordered by it.
 *
 * Points closer than the merge tolerance are merged into the first (in
 * input order, the first soup's before the second's). Sequential, and the
 * output depends on the input alone, so it is deterministic.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// Planar triangles over shared points, in any orientation.
struct PlanarSoup {
    std::vector<std::array<double, 2>> mPoints;
    std::vector<std::array<std::int64_t, 3>> mTriangles;
};

/// A covered region: triangles (counter-clockwise) and boundary loops.
struct PlanarRegion {
    bool mOk = true;
    std::string mWhy;  ///< why not, when not
    std::vector<std::array<double, 2>> mPoints;
    std::vector<std::array<std::int64_t, 3>> mTriangles;
    /// Closed polylines, the last point not repeated, the region on the left.
    std::vector<std::vector<std::int64_t>> mLoops;
};

/**
 * @brief The region `rAdd` covers minus the region `rSubtract` covers.
 * @param MergeTolerance points closer than this are one point (0: only
 *        points at exactly the same place).
 */
PlanarRegion planar_region(const PlanarSoup& rAdd, const PlanarSoup& rSubtract,
                           double MergeTolerance);

/// Twice the signed area of a closed loop (positive counter-clockwise).
double planar_loop_twice_area(const std::vector<std::array<double, 2>>& rPoints,
                              const std::vector<std::int64_t>& rLoop);

}  // namespace detail
}  // namespace meshioplusplus
