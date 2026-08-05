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
 * @file detail/tri_box.hpp
 * @brief Exact triangle / axis-aligned-box overlap, by the separating axis
 * theorem.
 *
 * `voxelize(VoxelFill::Surface)` needs to know which voxels a triangle passes
 * through, and "does the triangle's bounding box overlap the voxel" is not that
 * question -- a long thin diagonal triangle overlaps a great many boxes it never
 * actually enters. This is the exact test.
 *
 * ### The thirteen axes
 *
 * Two convex bodies are disjoint exactly when some axis separates their
 * projections, and for a triangle against an axis-aligned box it suffices to try
 * thirteen: the box's three face normals, the triangle's own normal, and the
 * nine cross products of the triangle's three edges with the box's three axes
 * (Akenine-Moller, *Fast 3D Triangle-Box Overlap Testing*). Dropping any of the
 * nine cross-product axes gives a test that is right for most configurations and
 * wrong for the edge-on ones, which is a bad failure to ship because it looks
 * like an off-by-one in the grid rather than a wrong predicate.
 *
 * ### Touching counts as overlapping
 *
 * The separation test is `min > rad || max < -rad`, so a triangle lying exactly
 * on a voxel's face is *not* separated and marks both neighbouring voxels. That
 * is deliberate: the alternative -- picking one side -- needs a tie-break rule
 * that has no geometric justification, and marking both keeps the occupied set a
 * closed cover of the surface, which is what a voxelization is for.
 *
 * Nothing here is transcendental: it is all `+ - *` and comparisons, so a numpy
 * transcription is bit-identical. Header-only and inline, because this runs once
 * per (triangle, candidate voxel) pair.
 */

// System includes
#include <cmath>
#include <cstddef>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

namespace tri_box_impl {

/// One separating-axis test: the three projections against the box's radius.
inline bool separated(double p0, double p1, double p2, double rad) {
    const double lo = p0 < p1 ? (p0 < p2 ? p0 : p2) : (p1 < p2 ? p1 : p2);
    const double hi = p0 > p1 ? (p0 > p2 ? p0 : p2) : (p1 > p2 ? p1 : p2);
    return lo > rad || hi < -rad;
}

/// Does the plane through @p rV with normal @p rNormal reach the box?
inline bool plane_box_overlap(const Vec3& rNormal, const Vec3& rV, const Vec3& rHalf) {
    Vec3 vmin{0.0, 0.0, 0.0};
    Vec3 vmax{0.0, 0.0, 0.0};
    for (std::size_t k = 0; k < 3; ++k) {
        if (rNormal[k] > 0.0) {
            vmin[k] = -rHalf[k] - rV[k];
            vmax[k] = rHalf[k] - rV[k];
        } else {
            vmin[k] = rHalf[k] - rV[k];
            vmax[k] = -rHalf[k] - rV[k];
        }
    }
    if (vec3_dot(rNormal, vmin) > 0.0)
        return false;
    return vec3_dot(rNormal, vmax) >= 0.0;
}

}  // namespace tri_box_impl

/**
 * @brief Does triangle (@p rA, @p rB, @p rC) overlap the axis-aligned box
 *        centred at @p rCentre with half-extents @p rHalf?
 *
 * Touching counts as overlapping -- see the file comment.
 */
inline bool tri_box_overlap(const Vec3& rCentre, const Vec3& rHalf, const Vec3& rA, const Vec3& rB,
                            const Vec3& rC) {
    using tri_box_impl::plane_box_overlap;
    using tri_box_impl::separated;

    // Work in the box's frame: the box becomes [-half, half] about the origin.
    const Vec3 v0 = vec3_sub(rA, rCentre);
    const Vec3 v1 = vec3_sub(rB, rCentre);
    const Vec3 v2 = vec3_sub(rC, rCentre);

    const Vec3 e0 = vec3_sub(v1, v0);
    const Vec3 e1 = vec3_sub(v2, v1);
    const Vec3 e2 = vec3_sub(v0, v2);

    // The nine edge x axis tests. Each cross product with a unit axis has a
    // closed form, written out rather than built through vec3_cross so the two
    // zero components are never computed or compared.
    const Vec3 edges[3] = {e0, e1, e2};
    const Vec3 verts[3] = {v0, v1, v2};
    for (std::size_t e = 0; e < 3; ++e) {
        const double ex = std::fabs(edges[e][0]);
        const double ey = std::fabs(edges[e][1]);
        const double ez = std::fabs(edges[e][2]);

        // axis = cross(x_hat, edge) = (0, -ez, ey)
        {
            double p[3];
            for (std::size_t i = 0; i < 3; ++i)
                p[i] = -edges[e][2] * verts[i][1] + edges[e][1] * verts[i][2];
            if (separated(p[0], p[1], p[2], ez * rHalf[1] + ey * rHalf[2]))
                return false;
        }
        // axis = cross(y_hat, edge) = (ez, 0, -ex)
        {
            double p[3];
            for (std::size_t i = 0; i < 3; ++i)
                p[i] = edges[e][2] * verts[i][0] - edges[e][0] * verts[i][2];
            if (separated(p[0], p[1], p[2], ez * rHalf[0] + ex * rHalf[2]))
                return false;
        }
        // axis = cross(z_hat, edge) = (-ey, ex, 0)
        {
            double p[3];
            for (std::size_t i = 0; i < 3; ++i)
                p[i] = -edges[e][1] * verts[i][0] + edges[e][0] * verts[i][1];
            if (separated(p[0], p[1], p[2], ey * rHalf[0] + ex * rHalf[1]))
                return false;
        }
    }

    // The three box face normals: the triangle's own bounding box against the box.
    for (std::size_t k = 0; k < 3; ++k)
        if (separated(v0[k], v1[k], v2[k], rHalf[k]))
            return false;

    // The triangle's plane.
    return plane_box_overlap(vec3_cross(e0, e1), v0, rHalf);
}

}  // namespace detail
}  // namespace meshioplusplus
