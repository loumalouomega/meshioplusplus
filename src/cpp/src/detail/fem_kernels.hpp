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

}  // namespace meshioplusplus::detail
