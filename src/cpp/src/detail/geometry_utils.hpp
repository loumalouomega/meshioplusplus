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
//                   SPDX-License-Identifier: MIT
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
#include <algorithm>
#include <array>
#include <cmath>
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

    /**
     * @brief The six dihedral angles of a tetrahedron `(rA, rB, rC, rD)` in
     * degrees, in the fixed edge order `(0,1), (0,2), (0,3), (1,2), (1,3),
     * (2,3)`. When `GuardZero` is true, a degenerate face pair (`den <= 0`)
     * yields `0.0`.
     */
    static std::array<double, 6> TetDihedralAnglesDeg(const Vec3& rA, const Vec3& rB,
                                                      const Vec3& rC, const Vec3& rD,
                                                      bool GuardZero = true) {
        const Vec3 p[4] = {rA, rB, rC, rD};
        constexpr int kEdge[6][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                     {1, 2, 0, 3}, {1, 3, 0, 2}, {2, 3, 0, 1}};
        constexpr double kRadToDeg = 180.0 / MathUtils::kPi;
        std::array<double, 6> out{};
        for (std::size_t e = 0; e < 6; ++e) {
            const Vec3& a = p[kEdge[e][0]];
            const Vec3 ab = vec3_sub(p[kEdge[e][1]], a);
            const Vec3 n1 = vec3_cross(ab, vec3_sub(p[kEdge[e][2]], a));
            const Vec3 n2 = vec3_cross(ab, vec3_sub(p[kEdge[e][3]], a));
            const double den = std::sqrt(vec3_norm_sq(n1) * vec3_norm_sq(n2));
            if (GuardZero && !(den > 0.0)) {
                out[e] = 0.0;
            } else {
                out[e] = std::acos(std::clamp(vec3_dot(n1, n2) / den, -1.0, 1.0)) * kRadToDeg;
            }
        }
        return out;
    }

    /**
     * @brief Sorted, deduplicated 1-ring vertex neighbours for a triangle
     * patch `rTris` over `NumPts` vertices.
     */
    static std::vector<std::vector<std::size_t>> TriangleVertexRings(
        const std::vector<std::array<std::size_t, 3>>& rTris, std::size_t NumPts) {
        std::vector<std::vector<std::size_t>> nb(NumPts);
        for (const auto& t : rTris) {
            for (int k = 0; k < 3; ++k) {
                nb[t[static_cast<std::size_t>(k)]].push_back(
                    t[static_cast<std::size_t>((k + 1) % 3)]);
                nb[t[static_cast<std::size_t>((k + 1) % 3)]].push_back(
                    t[static_cast<std::size_t>(k)]);
            }
        }
        for (auto& v : nb) {
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        }
        return nb;
    }

    ///@}
};

}  // namespace meshioplusplus::detail
