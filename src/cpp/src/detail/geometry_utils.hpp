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
 * @file detail/geometry_utils.hpp
 * @brief `GeometryUtils`: per-triangle measures that operations used to write
 * out by hand (doc/geometry_kernels.md, roadmap §7).
 *
 * A Kratos-style utility class (`Kratos::GeometryUtils`): static inline members
 * over corners already gathered as `Vec3` or as flat `double` triples. Each
 * keeps the arithmetic of the copies it replaced, operation for operation, so
 * a call site produces the same bits (see `math_utils.hpp`). Volumes and
 * scaled Jacobians stay in `fem_kernels.hpp`, polygon and polyhedron measures
 * in `polyhedron.hpp`.
 *
 * A **core-private** header (not installed): no API or ABI surface.
 */

// System includes
#include <array>
#include <cstddef>
#include <vector>

// Project includes
#include "math_utils.hpp"

namespace meshioplusplus::detail {

/// Static triangle measures.
struct GeometryUtils {
    ///@name Triangles in 3D
    ///@{

    /**
     * @brief The unnormalised normal `(b - a) x (c - a)`: right-handed about
     * the corner order, of length twice the triangle's area.
     */
    static Vec3 TriangleNormal(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
        return vec3_cross(vec3_sub(rB, rA), vec3_sub(rC, rA));
    }

    /// `TriangleNormal` of corners held as flat triples.
    static Vec3 TriangleNormal(const double* pA, const double* pB, const double* pC) {
        return vec3_cross(MathUtils::Sub(pB, pA), MathUtils::Sub(pC, pA));
    }

    /// `|(b - a) x (c - a)|`, twice the triangle's area.
    static double TriangleTwiceArea(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
        return vec3_norm(TriangleNormal(rA, rB, rC));
    }

    /// `TriangleTwiceArea` of corners held as flat triples.
    static double TriangleTwiceArea(const double* pA, const double* pB, const double* pC) {
        return vec3_norm(TriangleNormal(pA, pB, pC));
    }

    /**
     * @brief The triangle's three corner angles in degrees, corner `k` between
     * its edges to corners `k + 1` and `k + 2` (`MathUtils::AngleDeg`).
     */
    static std::array<double, 3> TriangleAnglesDeg(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
        const Vec3* p[3] = {&rA, &rB, &rC};
        std::array<double, 3> out{};
        for (std::size_t k = 0; k < 3; ++k) {
            const Vec3& o = *p[k];
            out[k] =
                MathUtils::AngleDeg(vec3_sub(*p[(k + 1) % 3], o), vec3_sub(*p[(k + 2) % 3], o));
        }
        return out;
    }

    /// The triangle's area, `0.5 * |(b - a) x (c - a)|`.
    static double TriangleArea(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
        return 0.5 * TriangleTwiceArea(rA, rB, rC);
    }

    /// `TriangleArea` of corners held as flat triples.
    static double TriangleArea(const double* pA, const double* pB, const double* pC) {
        return 0.5 * TriangleTwiceArea(pA, pB, pC);
    }

    ///@}
    ///@name Triangles in a parameter plane
    ///@{

    /**
     * @brief Twice the signed area of the triangle `(u0, v0), (u1, v1),
     * (u2, v2)`: `(u1 - u0)(v2 - v0) - (u2 - u0)(v1 - v0)`, positive when
     * counter-clockwise. Plain floating point; `predicates.hpp`'s `orient2d`
     * is the exact sign.
     */
    static double SignedTwiceArea2D(double U0, double V0, double U1, double V1, double U2,
                                    double V2) {
        return (U1 - U0) * (V2 - V0) - (U2 - U0) * (V1 - V0);
    }

    /**
     * @brief Whether `(u, v)` lies inside the closed polygon `(rPu[i], rPv[i])`
     * by the even-odd rule (a crossing count along `+u`). Several loops of one
     * region (an outer one and its holes) combine by XOR-ing their results.
     * A point on an edge may come out either way.
     */
    static bool PointInPolygon2D(const std::vector<double>& rPu, const std::vector<double>& rPv,
                                 double U, double V) {
        bool in = false;
        const std::size_t n = rPu.size();
        for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
            if ((rPv[i] > V) != (rPv[j] > V) &&
                U < (rPu[j] - rPu[i]) * (V - rPv[i]) / (rPv[j] - rPv[i]) + rPu[i])
                in = !in;
        }
        return in;
    }

    ///@}
};

}  // namespace meshioplusplus::detail
