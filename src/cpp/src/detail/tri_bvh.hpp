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
 * @file detail/tri_bvh.hpp
 * @brief A bounding-volume hierarchy over a triangle soup, for the fast
 * generalised winding number and for ray casts, used by `octree_mesh`.
 *
 * **The winding number.** `w(q) = (1/4π) Σ Ω_t(q)`, the solid angles of the
 * triangles seen from `q` (Van Oosterom and Strackee, IEEE Trans. Biomed.
 * Eng. 30(2), 1983), is 1 inside a closed, outward-wound surface and 0
 * outside, and degrades gracefully on holes, self-intersections and
 * non-manifold edges (Jacobson, Kavan and Sorkine-Hornung, "Robust
 * inside-outside segmentation using generalized winding numbers", SIGGRAPH
 * 2013). Evaluated directly it is O(T) per query; the tree makes it
 * O(log T) for points away from the surface (Barill, Dickson, Schmidt, Levin
 * and Jacobson, "Fast winding numbers for soups and clouds", SIGGRAPH 2018):
 * a node whose centre `p` is farther than `β` times its radius from `q` is
 * replaced by the far-field expansion of its triangles' dipoles,
 *
 *     Ω ≈ N·F(r) + J(r):M + ½ H(r)⋮C,   r = p - q,   F(r) = r / |r|^3,
 *
 * the Taylor expansion of `F` to second order (`J` and `H` its first and
 * second derivatives), where `N = Σ a_t n_t` is the area-weighted normal,
 * `M = Σ_t a_t n_t (c_t - p)^T` the first moment and
 * `C = Σ_t n_t ⊗ ∫_t (x - p)(x - p)^T dA` the second, integrated exactly over
 * each triangle (`∫ (x - c)(x - c)^T = (a / 12) Σ_v (v - c)(v - c)^T`); nearer
 * nodes are opened, and leaves sum the exact solid angles.
 *
 * **Determinism.** The tree is built over the triangles sorted by the Morton
 * key of their centroids (ties by triangle id) and split at the median count,
 * so its shape is a function of the triangle count alone; node data is
 * accumulated bottom-up in a fixed (left, right) order, and every query
 * traverses left first and adds in traversal order. A query's value is
 * therefore the same bits on every backend and thread count. The ray cast
 * keeps the nearest hit with ties broken by triangle id, so which node is
 * visited first is unobservable (the reason `surface_distance.hpp` gives for
 * not using a BVH for its own nearest-point search does not apply).
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"

namespace meshioplusplus {
namespace detail {

/// One node in preorder: the left child is the next node, the right child is
/// `mRight` (-1 for a leaf, whose triangles are `mOrder[mBegin, mBegin + mCount)`).
struct TriBvhNode {
    Vec3 mLo{0.0, 0.0, 0.0};
    Vec3 mHi{0.0, 0.0, 0.0};
    std::int32_t mBegin = 0;
    std::int32_t mCount = 0;
    std::int32_t mRight = -1;
};

/// The tree, and the far-field data of every node.
struct TriBvh {
    const TriangleSoup* mpSoup = nullptr;
    /// Soup triangle ids in leaf order.
    std::vector<std::int32_t> mOrder;
    std::vector<TriBvhNode> mNodes;
    /// Per node, the area-weighted centroid of its triangles.
    std::vector<Vec3> mCentre;
    /// Per node, `Σ a_t n_t` = `Σ ½ (b - a) × (c - a)`.
    std::vector<Vec3> mNormalSum;
    /// Per node, `Σ a_t n_t (c_t - p)^T`, row-major.
    std::vector<std::array<double, 9>> mMoment;
    /// Per node, `Σ n_t ⊗ ∫ (x - p)(x - p)^T dA`, index `9 i + 3 j + k`.
    std::vector<std::array<double, 27>> mSecond;
    /// Per node, an upper bound of the distance from `mCentre` to its corners.
    std::vector<double> mRadius;
};

/// Build the tree over @p rSoup (kept by pointer: it must outlive the tree).
TriBvh build_tri_bvh(const TriangleSoup& rSoup, std::int32_t LeafSize = 8);

/// The signed solid angle of triangle `(A, B, C)` seen from @p rQ: positive
/// when `rQ` lies behind it, on the side opposite `cross(B - A, C - A)`.
double tri_solid_angle(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rQ);

/// The generalised winding number at @p rQ; @p Beta <= 0 sums every triangle
/// exactly.
double fast_winding(const TriBvh& rTree, const Vec3& rQ, double Beta);

/// `fast_winding` at every point, in parallel.
std::vector<double> fast_windings(const TriBvh& rTree, const std::vector<Vec3>& rPoints,
                                  double Beta);

/// The nearest hit of a ray.
struct TriRayHit {
    double mT = 0.0;              ///< Ray parameter of the hit.
    std::int64_t mTriangle = -1;  ///< Soup triangle id; -1 when nothing was hit.
    Vec3 mNormal{0.0, 0.0, 0.0};  ///< The triangle's unnormalised normal.
};

/// The nearest triangle hit by `rOrigin + t rDir` with `TMin < t <= TMax`
/// (Möller and Trumbore, JGT 2(1), 1997), skipping triangle @p Skip; ties
/// go to the lower triangle id. Barycentric coordinates within 1e-10 of a
/// triangle's edge count as a hit, so a ray through a shared edge or vertex
/// does not slip between its triangles.
TriRayHit tri_bvh_ray(const TriBvh& rTree, const Vec3& rOrigin, const Vec3& rDir, double TMin,
                      double TMax, std::int64_t Skip = -1);

/// Where a segment crosses the skin.
struct TriSegmentHit {
    double mT = 0.0;  ///< Segment parameter in [0, 1].
    /// +1 when the segment leaves through the triangle's front (its normal
    /// side), -1 when it enters: the winding number drops or rises by one.
    double mSide = 0.0;
};

/// Every crossing of the segment `rA` to `rB` with the skin, ascending in `t`;
/// crossings of several triangles at one point (a shared edge or vertex) count
/// once, with the side their triangles agree on.
std::vector<TriSegmentHit> tri_bvh_segment_hits(const TriBvh& rTree, const Vec3& rA,
                                                const Vec3& rB);

}  // namespace detail
}  // namespace meshioplusplus
