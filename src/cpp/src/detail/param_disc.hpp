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
 * @file detail/param_disc.hpp
 * @brief Parametrise a disc-shaped triangle patch over the plane, for the
 * general freeform B-spline faces of `fit_cad`: least-squares conformal maps
 * (Lévy, Petitjean, Ray and Maillot, ACM TOG 2002) with a free boundary, and
 * mean-value coordinates (Floater, CAGD 2003) with the boundary on a circle,
 * which is bijective for any disc. A patch whose boundary is not one loop
 * (an annulus, a closed band, a sphere) is not a disc and is not
 * parametrised here: `fit_cad` splits it first.
 *
 * A **core-private** header (the `patch_flatten.hpp` precedent): nothing
 * installed names it. Deterministic: the linear systems are solved by the
 * thread-count independent `sparse_solve.hpp` and every loop runs in index
 * order.
 */

// System includes
#include <cstddef>
#include <vector>

// Project includes (private, not installed)
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief The boundary loops of a triangle patch, each a vertex cycle in the
 * direction the triangles run along it, in the order of their lowest vertex.
 * @return empty when the patch has no boundary; `rOk` false when a vertex has
 *         two boundary edges leaving it (a pinched or non-manifold boundary).
 */
std::vector<std::vector<std::size_t>> param_boundary_loops(const std::vector<PatchTri>& rTris,
                                                           bool& rOk);

/// Whether the patch is a topological disc: one boundary loop, Euler
/// characteristic `V - E + F = 1` over the vertices the triangles use.
bool param_is_disc(const std::vector<PatchTri>& rTris);

/// A parametrisation: per vertex `(u, v)`, and the triangles it flips.
struct DiscParam {
    bool mOk = false;
    std::vector<double> mU;
    std::vector<double> mV;
    /// Triangles whose orientation in the plane differs from the majority's
    /// (or of no area): a parametrisation is usable only with none.
    std::size_t mNumFlips = 0;
};

/**
 * @brief Least-squares conformal map, its two most distant boundary vertices
 * pinned at their 3-D distance apart (so lengths are roughly kept), solved by
 * conjugate gradients on the normal equations. Oriented so the triangles keep
 * the orientation they have in 3-D.
 */
DiscParam param_lscm(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris);

/**
 * @brief Mean-value coordinates with the boundary on a circle by arc length
 * (radius: the loop's length over 2π), solved by BiCGSTAB (the weights are not
 * symmetric). Bijective for a disc with positive weights.
 */
DiscParam param_mean_value(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris);

/// The triangles whose signed area in `(u, v)` has the minority sign or is zero.
std::size_t param_flips(const std::vector<PatchTri>& rTris, const std::vector<double>& rU,
                        const std::vector<double>& rV);

/**
 * @brief The parametrisation `fit_cad` uses: least-squares conformal maps,
 * then mean-value coordinates if those flip a triangle; `mOk` false for a
 * patch that is not a disc or that both flip.
 */
DiscParam param_disc(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris);

}  // namespace detail
}  // namespace meshioplusplus
