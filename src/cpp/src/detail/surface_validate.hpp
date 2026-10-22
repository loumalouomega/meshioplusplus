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
 * @file detail/surface_validate.hpp
 * @brief Checking a fitted freeform B-spline surface against the patch it
 * stands for, in both directions, and the noise floor of a patch -- for
 * `fit_cad`'s freeform modes (v17.12.0).
 *
 *  - **Noise.** A scanned skin is no more accurate than its own noise, so a
 *    freeform face is held to `min(3 tol, max(tol, 4.5 sigma))`, the analogue
 *    of the chordal-sag allowance analytic faces already get. `sigma` is
 *    estimated from the skin alone, never from the fit: each interior
 *    vertex's residual from the height quadric fitted, by least squares in
 *    its one ring's plane, to its two-ring without it -- the quadric absorbs
 *    the surface's curvature -- and `sigma = median |residual| / 0.6745 /
 *    sqrt(1 + 6 / k)` (the normal MAD constant, and the prediction's own
 *    variance for `k` neighbours), robust to a few spikes. Noise reads the
 *    same against the one ring's plane; a patch that reads more than 1.5
 *    times larger there is curvature the tessellation under-resolves, and
 *    reads as clean (0).
 *  - **Forward.** The surface sampled at its Greville points and on a `K x K`
 *    grid per knot span, the samples kept where their parameters fall inside
 *    the patch's parameter-space boundary (the trim), and each sample's
 *    distance to the patch's triangles: the surface must not stray where no
 *    vertex pins it (between samples, past the boundary within the trim).
 *  - **Reverse.** Each vertex projected onto the surface (Newton, from its
 *    fitted parameters): the patch must not stray from the surface.
 *
 * A **core-private** header: nothing installed names it. Deterministic: the
 * sample order is fixed and the distance query is thread-count independent.
 */

// System includes
#include <cstddef>
#include <vector>

// Project includes (private, not installed)
#include "bspline_lsq.hpp"
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/// The noise deviation of a triangle patch, estimated from its vertex rings
/// (see the file comment); 0 for a patch too small to tell.
double patch_noise_sigma(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris);

/// What `surface_validate` measured.
struct SurfaceValidation {
    /// Samples of the surface inside the trim, and their distances to the patch.
    std::size_t mNumSamples = 0;
    double mMaxForward = 0.0;
    double mRmsForward = 0.0;
    /// The patch's vertices' distances to the surface.
    double mMaxReverse = 0.0;
    double mRmsReverse = 0.0;
};

/**
 * @brief Sample a fitted surface against its patch in both directions.
 * @param rS the surface; @param rU, rV every vertex's parameters on it (the
 *        fit's foot points; vertices the triangles do not use are ignored).
 * @param SamplesPerSpan the grid per knot span (`K`, at least 1).
 * @param pDistanceTris the triangles the samples are measured against, when
 *        not the patch's own (a network patch is measured against its
 *        neighbours' too: its smooth boundary runs across their zigzag one);
 *        the trim is the patch's either way.
 */
SurfaceValidation surface_validate(const BsplineSurface& rS, const std::vector<PatchVec>& rPts,
                                   const std::vector<PatchTri>& rTris,
                                   const std::vector<double>& rU, const std::vector<double>& rV,
                                   int SamplesPerSpan = 3,
                                   const std::vector<PatchTri>* pDistanceTris = nullptr);

}  // namespace detail
}  // namespace meshioplusplus
