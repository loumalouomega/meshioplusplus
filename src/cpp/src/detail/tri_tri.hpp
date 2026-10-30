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
 * @file detail/tri_tri.hpp
 * @brief Whether and where two triangles in 3-D meet, decided exactly (roadmap
 * 6.9.1; shared with `find_intersections`).
 *
 * Every decision is a sign of Shewchuk's exact `orient3d` or `orient2d`
 * (`predicates.hpp`), so the answer is the exact one for the triangles'
 * double coordinates, in the spirit of P. Guigue and O. Devillers, "Fast and
 * robust triangle-triangle overlap test using orientation predicates",
 * Journal of Graphics Tools 8(1), 2003:
 *
 *  - **Apart.** When every corner of one triangle lies strictly on one side
 *    of the other's plane, they do not meet.
 *  - **Not in one plane.** Two triangles whose planes differ meet exactly
 *    when an edge of one meets the other (their intersection is a segment
 *    whose ends lie on edges). A segment `ab` meets the triangle `pqr` when
 *    `a` and `b` are not strictly on one side of its plane and the three
 *    orientations of `ab` against the triangle's edges, `orient3d(a, b, p,
 *    q)` and its turns, have no two strictly opposite signs.
 *  - **In one plane.** Projected on the plane of the two coordinates the
 *    normal leaves (dropping one coordinate keeps them exact), they meet
 *    exactly when two edges cross or touch, or a corner of one lies inside
 *    or on the other, by `orient2d`.
 *
 * Touching counts as meeting: a caller that wants proper crossings only
 * leaves out the pairs that share a corner. A degenerate triangle (its
 * corners on a line) is treated as its segments.
 *
 * `tri_tri_intersect_unshared` asks the question a self-intersection check
 * needs: whether two triangles of one surface meet anywhere other than at the
 * corners they share. `tri_tri_segment` constructs where two triangles meet;
 * its decision is the exact one above, its end points are computed in double.
 *
 * A **core-private** header: nothing installed names it.
 */

namespace meshioplusplus::detail {

/// What `tri_tri_segment` found.
enum class TriTriContact : int {
    None = 0,      ///< the triangles do not meet
    Segment = 1,   ///< they meet along a segment `(mP, mQ)`
    Point = 2,     ///< they meet at one point, `mP == mQ`
    Coplanar = 3,  ///< they lie in one plane and overlap; no segment is given
};

/// Where two triangles meet.
struct TriTriSegment {
    TriTriContact mKind = TriTriContact::None;
    double mP[3] = {0.0, 0.0, 0.0};
    double mQ[3] = {0.0, 0.0, 0.0};
};

/// Whether the triangles `(pA, pB, pC)` and `(pD, pE, pF)` (three coordinates
/// each) meet, touching included.
bool tri_tri_intersect(const double* pA, const double* pB, const double* pC, const double* pD,
                       const double* pE, const double* pF);

/// Whether the segment `(pA, pB)` meets the triangle `(pP, pQ, pR)`.
bool segment_tri_intersect(const double* pA, const double* pB, const double* pP, const double* pQ,
                           const double* pR);

/// Whether the triangles meet anywhere other than at the corners they share,
/// a corner being shared when its coordinates are bitwise equal in both.
/// Two triangles on one edge meet beyond it only when they lie in one plane
/// and fold onto each other; two on one corner meet beyond it when an edge
/// of one crosses the other or runs into it from the corner. Two identical
/// triangles always meet. Exact, like `tri_tri_intersect`.
bool tri_tri_intersect_unshared(const double* pA, const double* pB, const double* pC,
                                const double* pD, const double* pE, const double* pF);

/// Where the triangles meet: exact decision (`tri_tri_intersect`), end points
/// in double, where the edges that cross each plane cross it, clipped to the
/// line both planes hold.
TriTriSegment tri_tri_segment(const double* pA, const double* pB, const double* pC,
                              const double* pD, const double* pE, const double* pF);

}  // namespace meshioplusplus::detail
