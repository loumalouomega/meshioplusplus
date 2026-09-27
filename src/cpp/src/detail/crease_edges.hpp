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
 * @file detail/crease_edges.hpp
 * @brief The one crease test: which edges of a polygonal surface are sharp,
 * open, non-manifold or wound inconsistently.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * `feature_edges` reports these edges; `decimate`, `decimate_volume` and
 * `smooth` pin their endpoints. Before v16.23.0 each of the three carried its
 * own per-*vertex* test -- any two incident faces further apart than the
 * feature angle, adjacent or not -- which pinned every vertex of a coarsely
 * tessellated sphere. The test here is per *edge*: only the two faces that
 * share an edge are compared, so a smooth-but-coarse patch stays free.
 *
 * ### The rule
 *
 * The faces are polygons (a triangle is a 3-gon) given as rings of point ids;
 * each ring edge `(v_k, v_{k+1})` is one *use* of the undirected edge. Per
 * edge, with `u` uses:
 *
 *  - `u == 1`: **boundary**.
 *  - `u >= 3`: **non-manifold**.
 *  - `u == 2`, both uses from one face: not an edge of the surface (a ring
 *    that revisits an edge); skipped.
 *  - `u == 2` otherwise: the two faces' unit normals are compared. When the
 *    faces walk the edge the same way they disagree about "out" -- the pair is
 *    **inconsistent** -- and the second normal is negated first, so the angle
 *    is the dihedral the faces would have after a reorientation. The edge is
 *    **sharp** when `dot < cos(FeatureAngleDeg)`. A face with a zero normal
 *    (degenerate) gives no angle: the edge is never sharp.
 *
 * Collapsed ring edges (`v_k == v_{k+1}`) are not uses.
 *
 * ### Determinism
 *
 * Uses are sorted by a total order (`parallel_sort`), and every edge is
 * classified from its own run alone, so the result -- one record per reported
 * edge, ascending in `(lo, hi)` -- is the same on every backend and thread
 * count.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// One reported edge of a surface.
struct CreaseEdge {
    std::int64_t mLo = 0;  ///< the smaller endpoint id
    std::int64_t mHi = 0;  ///< the larger endpoint id
    /// How many face-ring edges reference it: 1 boundary, 2 a pair, 3+ non-manifold.
    std::int32_t mUses = 0;
    /// A pair whose two faces walk the edge the same way.
    bool mInconsistent = false;
    /// A pair whose (orientation-corrected) dihedral exceeds the feature angle.
    bool mSharp = false;
    /// That dihedral in degrees, in `[0, 180]`; NaN unless the edge is a pair of
    /// two non-degenerate faces.
    double mAngleDeg = 0.0;

    bool IsBoundary() const { return mUses == 1; }
    bool IsNonManifold() const { return mUses >= 3; }
    /// What a pinning operation freezes: a sharp or non-manifold edge. Open
    /// edges are the caller's separate "preserve boundary" decision.
    bool IsCrease() const { return mSharp || mUses >= 3; }
};

/**
 * @brief Classify the edges of a polygonal surface.
 *
 * @param rStart CSR offsets into @p rNodes, one more than the number of faces.
 * @param rNodes the face rings, concatenated.
 * @param rUnitNormal three entries per face: its unit normal, or `{0, 0, 0}`
 *        for a face with no direction.
 * @param FeatureAngleDeg the largest dihedral still treated as smooth.
 * @return every edge that is boundary, non-manifold, inconsistent or sharp,
 *         ascending in `(lo, hi)`. Smooth consistent pairs are omitted.
 */
std::vector<CreaseEdge> crease_edges(const std::vector<std::int64_t>& rStart,
                                     const std::vector<std::int64_t>& rNodes,
                                     const std::vector<double>& rUnitNormal,
                                     double FeatureAngleDeg);

/// `crease_edges` over a triangle list (three corners per face).
std::vector<CreaseEdge> crease_edges_triangles(const std::vector<std::int64_t>& rCorners,
                                               const std::vector<double>& rUnitNormal,
                                               double FeatureAngleDeg);

/// Set `rPinned[v] = 1` for both endpoints of every `IsCrease()` edge.
void pin_crease_endpoints(const std::vector<CreaseEdge>& rEdges,
                          std::vector<std::uint8_t>& rPinned);

/**
 * @brief Newell's unit normal of a ring of points, `{0, 0, 0}` when it has no
 * direction. For a triangle this is the direction of `cross(b - a, c - a)`.
 * @param pXyz three coordinates per point id.
 */
void ring_unit_normal(const double* pXyz, const std::int64_t* pRing, std::size_t Size,
                      double* pOut);

}  // namespace detail
}  // namespace meshioplusplus
