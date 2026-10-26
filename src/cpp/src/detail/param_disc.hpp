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
#include <array>
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
    /// Which map `param_disc` or `param_ladder` settled on: 0 none, 1
    /// least-squares conformal, 2 mean-value, 3 the best-fit plane, 4 the
    /// unrolled cylinder, 5 the square, 6 the annulus; and whether the
    /// patch was a disc at all.
    int mMethod = 0;
    bool mDisc = false;
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

/**
 * @brief Mean-value coordinates with the boundary on the unit square: the
 * loop `rLoop` (in the triangles' orientation) split at its four corners
 * (positions in it, ascending), side `j` laid by chord length along the
 * square's edge from corner `j` -- `(0, 0), (1, 0), (1, 1), (0, 1)` -- to
 * the next (v17.12.0, the quad layout's strips and the patch network). The
 * square is convex, so the map is bijective for positive weights; a triangle
 * with all three corners on one side of the square has no area there.
 * `mMethod` 5.
 */
DiscParam param_square(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       const std::vector<std::size_t>& rLoop,
                       const std::array<std::size_t, 4>& rCorners);

/**
 * @brief Mean-value coordinates of an annulus (two boundary loops, Euler
 * characteristic 0) between two concentric circles (`mesh_cad_surface`'s
 * band faces, roadmap 6.4.4): the longer loop on the circle of radius its
 * length over 2π, the other on the circle `exp(-2π w / L)` times as large,
 * the conformal modulus of a cylinder of circumference `L` (the loops' mean
 * length) and height `w` (the area over `L`), clamped to [1e-3, 0.9]. Each
 * loop is laid by arc length, the inner starting at its vertex nearest the
 * outer's first so the map does not twist, and running the other way round
 * as the triangles do. Floater's interior solve, as `param_mean_value`.
 * `mOk` false for a patch that is not an annulus; `mMethod` 6.
 */
DiscParam param_annulus(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris);

/// The triangles whose signed area in `(u, v)` has the minority sign or is zero.
std::size_t param_flips(const std::vector<PatchTri>& rTris, const std::vector<double>& rU,
                        const std::vector<double>& rV);

/**
 * @brief The parametrisation `fit_cad` uses: least-squares conformal maps,
 * then mean-value coordinates if those flip a triangle; `mOk` false for a
 * patch that is not a disc or that both flip.
 */
DiscParam param_disc(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris);

/// Which rungs `param_ladder` may take (`FitCadParameterization`'s values).
enum class ParamLadder : int {
    Auto = 0,        ///< projection, then conformal, then mean-value
    Projection = 1,  ///< the best-fit plane or the unrolled cylinder only
    Lscm = 2,        ///< least-squares conformal maps only
    Harmonic = 3,    ///< mean-value coordinates on a circle only
};

/**
 * @brief The parametrisation ladder of `fit_cad`'s freeform modes (v17.12.0):
 * the best-fit plane, then the unrolled cylinder (`patch_flatten`), each kept
 * only when it flips no triangle; then, for a disc, least-squares conformal
 * maps, then mean-value coordinates on a circle -- the rungs `Which` allows,
 * in that order. `mOk` false when every allowed rung flips a triangle or
 * needs a disc the patch is not; the caller then splits it.
 */
DiscParam param_ladder(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       ParamLadder Which);

}  // namespace detail
}  // namespace meshioplusplus
