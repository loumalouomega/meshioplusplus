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
 * @file detail/mesh_boolean.hpp
 * @brief The surface-boolean spike (roadmap 7.8): union, intersection and
 * difference of two closed triangle surfaces through an exact arrangement.
 * Not a public operation; `doc/spikes.md` records what it measured and the
 * decision it led to.
 *
 * The method:
 *
 *  - **Pairs.** Triangle pairs of A and B whose boxes overlap
 *    (`soup_box_pairs`), decided exactly (`tri_tri_intersect`).
 *  - **Exact points.** Every arrangement vertex is an input vertex or the
 *    point where an edge of one surface crosses the plane of a triangle of
 *    the other, held as exact homogeneous coordinates `(x, y, z, w)` in
 *    expansion arithmetic (`detail/expansion.hpp`): with `d_p = n . (p - a)`
 *    for the plane's normal `n`, the edge `pq` meets it at
 *    `(d_p q - d_q p) / (d_p - d_q)`. As long as neither surface intersects
 *    itself no other kind of point arises: two intersection segments inside
 *    one triangle meet only where the other surface's triangles meet, on its
 *    edges. Each pair's segment is the set of these points lying in both
 *    triangles, decided by exact `orient2d` on the homogeneous coordinates.
 *  - **Retriangulation.** Each cut triangle is retriangulated in its own
 *    plane (the dominant axis dropped) with its points and segments: the
 *    points inserted by exact point location, the segments recovered by edge
 *    flips. Points on a shared edge are inserted into both triangles beside
 *    it, so the result stays welded.
 *  - **Classification.** The pieces of each surface are grouped into patches
 *    bounded by the intersection curves; each patch is inside or outside the
 *    other surface by the parity of an exact ray cast from one sample point
 *    (re-cast along another direction when the ray grazes an edge).
 *  - **Selection.** Union keeps the outside patches of both; intersection the
 *    inside patches of both; difference A's outside and B's inside, turned.
 *
 * Coplanar contact (two triangles in one plane) is not handled: the call
 * fails naming it. Output points are rounded to doubles; the topology is
 * exact.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

enum class BooleanOp { Union, Intersection, Difference };

struct BooleanReport {
    std::int64_t mNumCandidatePairs = 0;  ///< triangle pairs whose boxes overlap
    std::int64_t mNumIntersectingPairs = 0;
    std::int64_t mNumCoplanarPairs = 0;
    std::int64_t mNumImplicitPoints = 0;  ///< edge-plane points made
    std::int64_t mNumSegments = 0;        ///< intersection segments
    std::int64_t mNumCutTriangles = 0;
    std::int64_t mNumSubTriangles = 0;
    std::int64_t mNumPatchesA = 0;
    std::int64_t mNumPatchesB = 0;
    std::int64_t mNumRayRetries = 0;
    std::int64_t mNumLabelConflicts = 0;  ///< neighbours across a curve on the same side
    std::int64_t mNumOutputTriangles = 0;
    /// Output pieces that rounding the points to doubles left collinear or
    /// turned over (every piece is positive exactly).
    std::int64_t mNumRoundingInverted = 0;
    double mSecondsPairs = 0.0;
    double mSecondsArrangement = 0.0;
    double mSecondsTriangulate = 0.0;
    double mSecondsClassify = 0.0;
};

struct BooleanResult {
    bool mOk = false;
    /// Why it failed (empty when `mOk`): "coplanar contact", "not closed", ...
    std::string mFailure;
    /// Points (rounded) and one `triangle` block; cell data `boolean:source`
    /// (0 from A, 1 from B).
    Mesh mMesh;
    BooleanReport mReport;
};

/// The boolean of two closed, manifold, consistently oriented triangle
/// surfaces (see the file comment). Never throws for geometry: a case it
/// cannot do comes back with `mOk` false and `mFailure`.
BooleanResult mesh_boolean_spike(const Mesh& rA, const Mesh& rB, BooleanOp Op);

/// The volume a closed triangle surface encloses (divergence theorem about
/// its centroid), the spike's oracle.
double closed_surface_volume(const Mesh& rSurface);

}  // namespace detail
}  // namespace meshioplusplus
