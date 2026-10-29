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
 * @brief Whether two triangles in 3-D meet, decided exactly (roadmap 6.9.1;
 * the kernel `find_intersections`, roadmap 3.1.2, will share).
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
 * A **core-private** header: nothing installed names it.
 */

namespace meshioplusplus::detail {

/// Whether the triangles `(pA, pB, pC)` and `(pD, pE, pF)` (three coordinates
/// each) meet, touching included.
bool tri_tri_intersect(const double* pA, const double* pB, const double* pC, const double* pD,
                       const double* pE, const double* pF);

/// Whether the segment `(pA, pB)` meets the triangle `(pP, pQ, pR)`.
bool segment_tri_intersect(const double* pA, const double* pB, const double* pP, const double* pQ,
                           const double* pR);

}  // namespace meshioplusplus::detail
