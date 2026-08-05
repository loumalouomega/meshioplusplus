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
 * @file detail/point_triangle.hpp
 * @brief The closest point on a triangle to a query point, and **which feature
 * of the triangle it lies on**.
 *
 * The feature is not diagnostic decoration: it is what selects the normal used
 * to sign a distance. A closest point in the triangle's interior takes the
 * triangle's own normal; one on an edge or at a vertex must take a blend of the
 * incident faces' normals, or the sign comes out wrong on the concave side of
 * every crease. See `operations/sdf.hpp`.
 *
 * ### The method, and the one that looks like it but is wrong
 *
 * This is Ericson's barycentric region classification (*Real-Time Collision
 * Detection*, §5.1.5). The tempting alternative -- project onto the triangle's
 * plane and clamp the barycentric coordinates into range -- is **wrong for
 * obtuse triangles**, where the clamped point can land on the far edge rather
 * than the near one, and is the commonest way this gets shipped broken.
 *
 * The seven regions are mutually exclusive and exhaustive with the comparisons
 * written exactly as they are. Every branch tests a quantity built by the same
 * sequence of multiplications and subtractions on both sides of its boundary, so
 * there is no configuration in which two regions both claim a point or neither
 * does.
 *
 * ### Why it is arranged for a bit-exact numpy twin
 *
 * There is not one transcendental function here: the whole computation is
 * `+ - * /` and comparisons, and `sqrt` is left to the caller so that a
 * min-reduction over many triangles can compare squared distances. IEEE requires
 * both `/` and `sqrt` to be correctly rounded, so a numpy transcription that
 * follows this expression order produces bit-identical results -- which is what
 * `tests/python/test_surface_distance.py::test_cpp_matches_python` asserts.
 *
 * A degenerate (zero-area) triangle is handled explicitly rather than left to
 * chance: it makes the denominator zero, and dividing by it would poison a
 * min-reduction with a NaN that then wins every comparison.
 *
 * Header-only and inline: this runs once per (query, candidate triangle) pair,
 * which is the hottest loop in the whole feature. `detail::` helpers are exempt
 * from the unique-prefix rule that `src/cpp/src/**.cpp` follows.
 */

// System includes
#include <cstddef>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

/// Which feature of a triangle a closest point lies on.
enum class TriangleFeature {
    VertexA = 0,
    VertexB = 1,
    VertexC = 2,
    EdgeAB = 3,
    EdgeBC = 4,
    EdgeCA = 5,
    Face = 6,
};

/// The closest point on a triangle, and where on it.
struct PointTriangleHit {
    Vec3 mPoint{0.0, 0.0, 0.0};   ///< The closest point itself.
    double mDistanceSq = 0.0;     ///< Its squared distance from the query.
    TriangleFeature mFeature = TriangleFeature::Face;  ///< Which feature it is on.
};

/**
 * @brief The closest point on triangle (A, B, C) to P.
 *
 * The expression order is load-bearing -- see the file comment. Do not
 * "simplify" it, and in particular do not reorder the region tests: they are
 * arranged so that each one can assume the preceding ones failed.
 */
inline PointTriangleHit closest_point_on_triangle(const Vec3& rP, const Vec3& rA, const Vec3& rB,
                                                  const Vec3& rC) {
    PointTriangleHit hit;
    const Vec3 ab = vec3_sub(rB, rA);
    const Vec3 ac = vec3_sub(rC, rA);
    const Vec3 ap = vec3_sub(rP, rA);
    const double d1 = vec3_dot(ab, ap);
    const double d2 = vec3_dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        hit.mPoint = rA;
        hit.mFeature = TriangleFeature::VertexA;
        hit.mDistanceSq = vec3_norm_sq(vec3_sub(rP, hit.mPoint));
        return hit;
    }

    const Vec3 bp = vec3_sub(rP, rB);
    const double d3 = vec3_dot(ab, bp);
    const double d4 = vec3_dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) {
        hit.mPoint = rB;
        hit.mFeature = TriangleFeature::VertexB;
        hit.mDistanceSq = vec3_norm_sq(vec3_sub(rP, hit.mPoint));
        return hit;
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        hit.mPoint = vec3_add(rA, vec3_scale(ab, v));
        hit.mFeature = TriangleFeature::EdgeAB;
        hit.mDistanceSq = vec3_norm_sq(vec3_sub(rP, hit.mPoint));
        return hit;
    }

    const Vec3 cp = vec3_sub(rP, rC);
    const double d5 = vec3_dot(ab, cp);
    const double d6 = vec3_dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) {
        hit.mPoint = rC;
        hit.mFeature = TriangleFeature::VertexC;
        hit.mDistanceSq = vec3_norm_sq(vec3_sub(rP, hit.mPoint));
        return hit;
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        hit.mPoint = vec3_add(rA, vec3_scale(ac, w));
        hit.mFeature = TriangleFeature::EdgeCA;
        hit.mDistanceSq = vec3_norm_sq(vec3_sub(rP, hit.mPoint));
        return hit;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        hit.mPoint = vec3_add(rB, vec3_scale(vec3_sub(rC, rB), w));
        hit.mFeature = TriangleFeature::EdgeBC;
        hit.mDistanceSq = vec3_norm_sq(vec3_sub(rP, hit.mPoint));
        return hit;
    }

    const double den = va + vb + vc;
    if (!(den > 0.0)) {
        // Degenerate triangle: no interior to project into, and the face branch
        // below would divide by zero and hand a NaN to a min-reduction, where it
        // would then win every comparison. Fall back to the three edges, taking
        // the nearest and breaking ties by edge ordinal so the answer does not
        // depend on how the compiler ordered the comparisons.
        const Vec3 ends[3][2] = {{rA, rB}, {rB, rC}, {rC, rA}};
        const TriangleFeature feats[3] = {TriangleFeature::EdgeAB, TriangleFeature::EdgeBC,
                                          TriangleFeature::EdgeCA};
        bool first = true;
        for (std::size_t e = 0; e < 3; ++e) {
            const Vec3 d = vec3_sub(ends[e][1], ends[e][0]);
            const double len2 = vec3_norm_sq(d);
            double t = 0.0;
            if (len2 > 0.0) {
                t = vec3_dot(vec3_sub(rP, ends[e][0]), d) / len2;
                t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
            }
            const Vec3 q = vec3_add(ends[e][0], vec3_scale(d, t));
            const double d2q = vec3_norm_sq(vec3_sub(rP, q));
            if (first || d2q < hit.mDistanceSq) {
                first = false;
                hit.mPoint = q;
                hit.mDistanceSq = d2q;
                hit.mFeature = feats[e];
            }
        }
        return hit;
    }

    const double inv = 1.0 / den;
    const double v = vb * inv;
    const double w = vc * inv;
    hit.mPoint = vec3_add(rA, vec3_add(vec3_scale(ab, v), vec3_scale(ac, w)));
    hit.mFeature = TriangleFeature::Face;
    hit.mDistanceSq = vec3_norm_sq(vec3_sub(rP, hit.mPoint));
    return hit;
}

}  // namespace detail
}  // namespace meshioplusplus
