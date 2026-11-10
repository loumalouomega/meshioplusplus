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

#pragma once

/**
 * @file detail/fem_kernels.hpp
 * @brief Shared per-element kernels: the formulas several operations used to
 * carry private copies of (doc/geometry_kernels.md).
 *
 * Everything here is a plain inline function or a `constexpr` table over
 * `detail::Vec3` corners already gathered through the uniform mesh API, so it
 * works under every mesh backend and costs nothing over an inline copy. The
 * header is private (not installed), so a change here is not an ABI event.
 * Kernels over a whole mesh live beside it (`cell_centroids.hpp`).
 *
 * Each kernel keeps the exact arithmetic of the copies it replaced: the
 * operations' outputs are pinned bit for bit by the determinism digests
 * (doc/benchmarks.md, "Determinism check"). A copy whose arithmetic differs
 * (another determinant expansion, another degenerate-corner value) is not
 * folded in here until its own output is proven unchanged.
 */

// System includes
#include <algorithm>
#include <cmath>
#include <limits>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus::detail {

/**
 * @brief Six times the signed volume of the tetrahedron (A, B, C, D): the
 * triple product of its three edges out of A. Positive when D lies on the side
 * of (A, B, C) its right-hand normal points to (meshio's order).
 */
inline double tet_volume6(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) {
    return triple_product(vec3_sub(rB, rA), vec3_sub(rC, rA), vec3_sub(rD, rA));
}

/** @brief The signed volume of the tetrahedron (A, B, C, D), `tet_volume6 / 6`. */
inline double tet_volume(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) {
    return tet_volume6(rA, rB, rC, rD) / 6.0;
}

/**
 * @brief The determinant of the three edge vectors at a corner, each scaled to
 * unit length first: the per-corner factor of the scaled Jacobian (Verdict's
 * definition, as `quality` reports it).
 * @return NaN when an edge is shorter than @p Eps.
 */
inline double det_unit3(const Vec3& rA, const Vec3& rB, const Vec3& rC, double Eps) {
    const double na = vec3_norm(rA);
    const double nb = vec3_norm(rB);
    const double nc = vec3_norm(rC);
    if (na < Eps || nb < Eps || nc < Eps)
        return std::numeric_limits<double>::quiet_NaN();
    return det3(vec3_scale(rA, 1.0 / na), vec3_scale(rB, 1.0 / nb), vec3_scale(rC, 1.0 / nc));
}

/**
 * @brief The corner's scaled Jacobian as a ratio, `det(a, b, c) / (|a||b||c|)`
 * (the determinant is not normalised edge by edge, unlike `det_unit3`).
 * @return @p Degenerate when an edge has zero (or NaN) length.
 */
inline double corner_jacobian_ratio(const Vec3& rA, const Vec3& rB, const Vec3& rC,
                                    double Degenerate) {
    const double l = vec3_norm(rA) * vec3_norm(rB) * vec3_norm(rC);
    if (!(l > 0.0))
        return Degenerate;
    return triple_product(rA, rB, rC) / l;
}

/**
 * @brief Each corner of a linear hexahedron (VTK/meshio order) and its three
 * edge neighbours, right-handed for a well-oriented cell: the corner frames of
 * the scaled Jacobian (Knupp's hexahedron corners).
 */
inline constexpr int kHexCornerEdges[8][3] = {{1, 3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2, 7},
                                              {7, 5, 0}, {4, 6, 1}, {5, 7, 2}, {6, 4, 3}};

/**
 * @brief Each corner of a linear wedge (meshio order) and its three edge
 * neighbours, in the order the wedge's scaled Jacobian reads them.
 */
inline constexpr int kWedgeCornerEdges[6][3] = {{1, 2, 3}, {2, 0, 4}, {0, 1, 5},
                                                {5, 4, 0}, {3, 5, 1}, {4, 3, 2}};

/// Degrees per radian, as `optimize_volume` has always spelled it.
inline constexpr double kTetRadToDeg = 57.295779513082320876798;

/**
 * @brief The least dihedral angle of the tetrahedron (A, B, C, D), in degrees,
 * whatever its orientation; 0 for a degenerate one.
 *
 * The angle through each edge is the one between the normals of the two faces
 * through it. A cosine more than 1e-12 below the largest one already evaluated
 * cannot give the least angle (acos is decreasing with slope at least 1 and
 * accurate to an ulp or so), so its acos is skipped: the result is bit for bit
 * what evaluating all six gives (roadmap 5.5.2). Moved verbatim from
 * `optimize_volume`, whose digests pin it.
 */
inline double tet_min_dihedral(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) {
    const Vec3 p[4] = {rA, rB, rC, rD};
    // The edge (i, j) has the opposite edge (k, l); its dihedral angle is the
    // angle between the normals of faces (i, j, k) and (i, j, l).
    constexpr int kEdge[6][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                 {1, 2, 0, 3}, {1, 3, 0, 2}, {2, 3, 0, 1}};
    constexpr double kAcosSkip = 1e-12;
    double least = 180.0;
    double c_most = -2.0;  // the largest cosine whose acos was evaluated
    for (const auto& e : kEdge) {
        const Vec3& a = p[e[0]];
        const Vec3 ab = vec3_sub(p[e[1]], a);
        const Vec3 n1 = vec3_cross(ab, vec3_sub(p[e[2]], a));
        const Vec3 n2 = vec3_cross(ab, vec3_sub(p[e[3]], a));
        const double den = std::sqrt(vec3_dot(n1, n1) * vec3_dot(n2, n2));
        if (!(den > 0.0))
            return 0.0;
        const double c = std::clamp(vec3_dot(n1, n2) / den, -1.0, 1.0);
        if (!(c > c_most - kAcosSkip))
            continue;
        least = std::min(least, std::acos(c) * kTetRadToDeg);
        c_most = std::max(c_most, c);
    }
    return least;
}

/**
 * @brief The least and the largest dihedral angle of (A, B, C, D), in degrees,
 * from all six; both 0 for a degenerate tetrahedron.
 */
inline void tet_dihedral_range(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD,
                               double& rMin, double& rMax) {
    const Vec3 p[4] = {rA, rB, rC, rD};
    constexpr int kEdge[6][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                 {1, 2, 0, 3}, {1, 3, 0, 2}, {2, 3, 0, 1}};
    rMin = 180.0;
    rMax = 0.0;
    for (const auto& e : kEdge) {
        const Vec3& a = p[e[0]];
        const Vec3 ab = vec3_sub(p[e[1]], a);
        const Vec3 n1 = vec3_cross(ab, vec3_sub(p[e[2]], a));
        const Vec3 n2 = vec3_cross(ab, vec3_sub(p[e[3]], a));
        const double den = std::sqrt(vec3_dot(n1, n1) * vec3_dot(n2, n2));
        if (!(den > 0.0)) {
            rMin = rMax = 0.0;
            return;
        }
        const double angle = std::acos(std::clamp(vec3_dot(n1, n2) / den, -1.0, 1.0)) * kTetRadToDeg;
        rMin = std::min(rMin, angle);
        rMax = std::max(rMax, angle);
    }
}

/**
 * @brief The mean ratio `12 (3 V)^(2/3) / sum(l^2)` of (A, B, C, D): 1 for the
 * regular tetrahedron, falling to 0 as it flattens; 0 when the signed volume is
 * not positive. Liu and Joe's measure, the reciprocal of the conformal AMIPS
 * energy (`tet_amips`).
 */
inline double tet_mean_ratio(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) {
    const double v = tet_volume6(rA, rB, rC, rD) / 6.0;
    if (!(v > 0.0))
        return 0.0;
    const double l2 = vec3_norm_sq(vec3_sub(rB, rA)) + vec3_norm_sq(vec3_sub(rC, rA)) +
                      vec3_norm_sq(vec3_sub(rD, rA)) + vec3_norm_sq(vec3_sub(rC, rB)) +
                      vec3_norm_sq(vec3_sub(rD, rB)) + vec3_norm_sq(vec3_sub(rD, rC));
    if (!(l2 > 0.0))
        return 0.0;
    return 12.0 * std::cbrt(9.0 * v * v) / l2;
}

/**
 * @brief The conformal AMIPS energy of (A, B, C, D) against the regular
 * tetrahedron, `tr(J^T J) / (3 det(J)^(2/3))`: 1 for the regular tetrahedron,
 * growing without bound as it flattens, +infinity when the signed volume is
 * not positive. Equal to `1 / tet_mean_ratio` (Hu et al., "Fast Tetrahedral
 * Meshing in the Wild", 2020, use it as their quality).
 */
inline double tet_amips(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) {
    const double q = tet_mean_ratio(rA, rB, rC, rD);
    return q > 0.0 ? 1.0 / q : std::numeric_limits<double>::infinity();
}

}  // namespace meshioplusplus::detail
