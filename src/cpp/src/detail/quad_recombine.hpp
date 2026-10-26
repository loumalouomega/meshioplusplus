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
 * @file detail/quad_recombine.hpp
 * @brief Greedy recombination of a surface triangulation into a
 * quad-dominant mesh (roadmap 6.8).
 *
 * Each interior edge whose two triangles make a convex quadrilateral is a
 * candidate, scored by the quad's quality: the least over its corners of
 * `1 - |theta - 90| / 90` (1 for a rectangle, 0 at a corner of 0 or 180
 * degrees), times the cosine of the dihedral angle between the two
 * triangles, so a pair folded across a crease scores low. Candidates are
 * merged best first (ties to the lower edge key) while both triangles are
 * free and the score reaches the floor, the greedy order of H. Borouchaki
 * and P. J. Frey, "Adaptive triangular-quadrilateral mesh generation", Int.
 * J. Numer. Meth. Engng. 41 (1998); C. K. Lee and S. H. Lo, "A new scheme for
 * the generation of a graded quadrilateral mesh", Comput. Struct. 52 (1994),
 * merge front by front instead. A minimum-cost perfect matching (Blossom-quad,
 * Remacle et al. 2012) would leave fewer triangles; its libraries are
 * excluded by licence (roadmap 6.8.1). Deterministic: the order depends on
 * the input alone.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

namespace meshioplusplus::detail {

/// The quality of quad `rQ` (point ids into `pX`, three coordinates each,
/// counter-clockwise about its normal) in [0, 1]: 0 when it is not convex.
double quad_recombine_quality(const double* pX, const std::array<std::int64_t, 4>& rQ);

/// Verdict's scaled Jacobian of a quad against its normal (the measure
/// `compute_quality` reports): the least corner `(e1 x e2) . n / (|e1| |e2|)`.
double quad_scaled_jacobian(const double* pX, const std::array<std::int64_t, 4>& rQ);

struct QuadRecombineResult {
    /// Each quad, and the two triangles (indices into the input) it merges:
    /// `(a, d, b, c)` from `(a, b, c)` and `(b, a, d)`, so it runs as they do.
    std::vector<std::array<std::int64_t, 4>> mQuads;
    std::vector<std::array<std::size_t, 2>> mPairs;
    /// The triangles left, as indices into the input, ascending.
    std::vector<std::size_t> mLeft;
};

/**
 * @brief Recombine consistently oriented triangles `rTris` (point ids) whose
 * pairs may merge only within one `rGroup` (one value per triangle, e.g. the
 * CAD face it lies in); quads come out in the order of their first triangle.
 */
QuadRecombineResult quad_recombine(const std::vector<std::array<std::int64_t, 3>>& rTris,
                                   const std::vector<std::int64_t>& rGroup, const double* pX,
                                   double MinQuality);

}  // namespace meshioplusplus::detail
