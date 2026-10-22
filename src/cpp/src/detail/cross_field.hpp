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
 * @file detail/cross_field.hpp
 * @brief The smoothest cross field (4-RoSy direction field) on a triangle
 * patch, aligned with its boundary and feature edges, and its singularities
 * -- the guide `quad_layout` traces a patch network along (`fit_cad`'s
 * network mode, v17.12.0).
 *
 * Each triangle carries a frame (its first edge, its normal) and a cross,
 * stored as the complex number `z = exp(4 i theta)` so the four branches
 * `theta + k pi / 2` are one value. The field minimises the discrete
 * connection Dirichlet energy `sum_e |z_f - c_e z_g|^2` over interior edges,
 * `c_e = exp(4 i (alpha_f - alpha_g))` carrying the cross across the shared
 * edge (`alpha` the edge's angle in each frame): the smoothest field of
 * Knöppel, Crane, Pinkall and Schröder ("Globally optimal direction fields",
 * SIGGRAPH 2013), with the boundary and feature edges as hard constraints
 * (a constrained triangle's cross holds that edge's direction), plus, where
 * asked, a soft pull towards the principal curvature directions weighted by
 * their anisotropy. The Hermitian system is solved as its real symmetric
 * positive definite twin by the deterministic conjugate gradients of
 * `sparse_solve.hpp`; with nothing to constrain it (a closed patch with no
 * curvature term) the field is the smallest eigenvector, by inverse iteration
 * from a fixed start.
 *
 * A vertex's index is the cross's turning around it, `(sum of the wrapped
 * mismatches + the angle defect) / 2 pi`, a multiple of `1/4`: `+1/4` a
 * valence-3 singularity, `-1/4` a valence-5 one. Interior vertices only; on
 * the boundary the layout reads corners from the boundary's turning instead.
 *
 * A **core-private** header: nothing installed names it. Deterministic:
 * fixed loop orders, a fixed start and iteration count, thread-count
 * independent solves.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/// What `cross_field` aligns with.
struct CrossFieldOptions {
    /// Interior edges (vertex pairs) the field must follow, besides the
    /// boundary: the patch's feature edges.
    std::vector<std::pair<std::size_t, std::size_t>> mFeatureEdges;
    /// Whether the boundary edges constrain the field.
    bool mConstrainBoundary = true;
    /// The weight of the pull towards the principal curvature directions,
    /// times each triangle's anisotropy `(|k1| - |k2|) / (|k1| + |k2|)`; 0 for none.
    double mCurvatureWeight = 0.0;
    /// Inverse iterations for an unconstrained field.
    int mPowerIterations = 30;
};

/// A cross field on a patch.
struct CrossField {
    bool mOk = false;
    /// Per triangle: its frame (unit first edge, unit in-plane normal to it,
    /// unit normal) and the angle of the cross's first branch in it.
    std::vector<PatchVec> mE1, mE2, mNormal;
    std::vector<double> mTheta;
    /// Per vertex: its index in quarters (`round(4 index)`), 0 on the
    /// boundary and at regular vertices.
    std::vector<int> mIndex;
    /// The vertices with a nonzero index, ascending.
    std::vector<std::size_t> mSingular;
    /// Whether every triangle was constrained, solved or both (false for a
    /// patch the solve could not reach).
    std::size_t mNumConstrained = 0;
};

/// The smoothest cross field of a patch (see the file comment).
CrossField cross_field(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       const CrossFieldOptions& rOptions);

/// Branch `K` (taken mod 4) of the cross in triangle `F`, a unit 3-D vector.
PatchVec cross_branch(const CrossField& rField, std::size_t F, int K);

/// A singular vertex's valence in a quad layout: `4 - index`.
inline int cross_valence(const CrossField& rField, std::size_t V) {
    return 4 - rField.mIndex[V];
}

}  // namespace detail
}  // namespace meshioplusplus
