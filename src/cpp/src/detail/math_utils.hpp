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
 * @file detail/math_utils.hpp
 * @brief `MathUtils`: the small vector and 3x3 algebra that operations used to
 * write out by hand (doc/geometry_kernels.md, roadmap §8).
 *
 * A Kratos-style utility class (`Kratos::MathUtils`): static inline members, no
 * state, called as `MathUtils::Cross(a, b)` in place of the copied block. The
 * `Vec3` members forward to `detail/geometry.hpp`'s `vec3_*` kernels; the
 * `const double*` overloads compute the same expressions in the same order, for
 * call sites that hold flat coordinate buffers.
 *
 * Every member keeps the exact floating-point arithmetic of the copies it
 * replaced, so a call site that switches to it produces the same bits (the core
 * builds with `-ffp-contract=off`, so an inlined expression rounds as the
 * hand-written one did). A copy whose arithmetic differs (a reciprocal where
 * this divides, another guard, another determinant expansion) is not folded in.
 *
 * A **core-private** header (not installed), so it adds nothing to the API or
 * the ABI.
 */

// System includes
#include <algorithm>
#include <cmath>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus::detail {

/// Static vector and small-matrix helpers over `Vec3` and flat `double` triples.
struct MathUtils {
    /// pi, as the copies spelled it (`3.14159265358979323846`).
    static constexpr double kPi = 3.14159265358979323846;

    ///@name Vectors
    ///@{

    /// `a - b`.
    static Vec3 Sub(const Vec3& rA, const Vec3& rB) { return vec3_sub(rA, rB); }

    /// `a - b` of two flat triples.
    static Vec3 Sub(const double* pA, const double* pB) {
        return {pA[0] - pB[0], pA[1] - pB[1], pA[2] - pB[2]};
    }

    /// `a + b`.
    static Vec3 Add(const Vec3& rA, const Vec3& rB) { return vec3_add(rA, rB); }

    /// `s * a`.
    static Vec3 Scale(const Vec3& rA, double S) { return vec3_scale(rA, S); }

    /// `a . b`, summed left to right.
    static double Dot(const Vec3& rA, const Vec3& rB) { return vec3_dot(rA, rB); }

    /// `a . b` of two flat triples, summed left to right.
    static double Dot(const double* pA, const double* pB) {
        return pA[0] * pB[0] + pA[1] * pB[1] + pA[2] * pB[2];
    }

    /// `a . b` of a flat triple and a `Vec3`, summed left to right.
    static double Dot(const double* pA, const Vec3& rB) {
        return pA[0] * rB[0] + pA[1] * rB[1] + pA[2] * rB[2];
    }

    /// `a x b`.
    static Vec3 Cross(const Vec3& rA, const Vec3& rB) { return vec3_cross(rA, rB); }

    /// `a x b` of two flat triples.
    static Vec3 Cross(const double* pA, const double* pB) {
        return {pA[1] * pB[2] - pA[2] * pB[1], pA[2] * pB[0] - pA[0] * pB[2],
                pA[0] * pB[1] - pA[1] * pB[0]};
    }

    /// `a x b` of two flat triples, written to `pOut` (which must not alias them).
    static void Cross(const double* pA, const double* pB, double* pOut) {
        pOut[0] = pA[1] * pB[2] - pA[2] * pB[1];
        pOut[1] = pA[2] * pB[0] - pA[0] * pB[2];
        pOut[2] = pA[0] * pB[1] - pA[1] * pB[0];
    }

    /// `|a|^2`.
    static double NormSq(const Vec3& rA) { return vec3_norm_sq(rA); }

    /// `|a|`.
    static double Norm(const Vec3& rA) { return vec3_norm(rA); }

    /// `|a - b|`.
    static double Distance(const Vec3& rA, const Vec3& rB) { return vec3_norm(vec3_sub(rA, rB)); }

    /// `|a - b|` of two flat triples.
    static double Distance(const double* pA, const double* pB) { return vec3_norm(Sub(pA, pB)); }

    /**
     * @brief `a / |a|`, each component divided by the length, or the zero
     * vector when `!(|a| > 0)` (a zero or NaN length).
     *
     * The division form, which `vec3_normalize`'s multiplication by the
     * reciprocal does not match bit for bit; it is the form `quad_layout`,
     * `feature_fit`, `cross_field`, `boundary_layer` and the external mesher's
     * STL normals used.
     */
    static Vec3 UnitByDivision(const Vec3& rA) {
        const double l = vec3_norm(rA);
        if (!(l > 0.0))
            return {0.0, 0.0, 0.0};
        return {rA[0] / l, rA[1] / l, rA[2] / l};
    }

    /**
     * @brief The angle between `a` and `b` in degrees, `atan2(|a x b|, a . b)`:
     * well conditioned near 0 and 180 degrees, where `acos` of the normalised
     * dot product is not.
     */
    static double AngleDeg(const Vec3& rA, const Vec3& rB) {
        return std::atan2(vec3_norm(vec3_cross(rA, rB)), vec3_dot(rA, rB)) * (180.0 / kPi);
    }

    ///@}
    ///@name Small matrices
    ///@{

    /**
     * @brief Solves the symmetric 3x3 system `M x = b` in closed cofactor
     * (adjugate) form.
     *
     * `M` is given by its upper triangle. The system counts as singular unless
     * `|det M| > 1e-12 * s^3`, `s` the largest magnitude of the six entries
     * (written `s * s * s`: `std::pow` rounds differently). On success writes
     * `x` to `pOut` and returns `true`; otherwise leaves `pOut` alone and
     * returns `false`. `decimate`'s optimal placement and `gradient`'s least
     * squares share it; the arithmetic is theirs, operation for operation.
     */
    static bool SolveSymmetric3(double M00, double M01, double M02, double M11, double M12,
                                double M22, const double* pB, double* pOut) {
        const double c00 = M11 * M22 - M12 * M12;
        const double c01 = M02 * M12 - M01 * M22;
        const double c02 = M01 * M12 - M02 * M11;
        const double det = M00 * c00 + M01 * c01 + M02 * c02;
        double scale = std::abs(M00);
        scale = std::max(scale, std::abs(M01));
        scale = std::max(scale, std::abs(M02));
        scale = std::max(scale, std::abs(M11));
        scale = std::max(scale, std::abs(M12));
        scale = std::max(scale, std::abs(M22));
        if (!(std::abs(det) > 1e-12 * (scale * scale * scale)))
            return false;
        const double c11 = M00 * M22 - M02 * M02;
        const double c12 = M01 * M02 - M00 * M12;
        const double c22 = M00 * M11 - M01 * M01;
        const double inv = 1.0 / det;
        pOut[0] = (c00 * pB[0] + c01 * pB[1] + c02 * pB[2]) * inv;
        pOut[1] = (c01 * pB[0] + c11 * pB[1] + c12 * pB[2]) * inv;
        pOut[2] = (c02 * pB[0] + c12 * pB[1] + c22 * pB[2]) * inv;
        return true;
    }

    ///@}
};

}  // namespace meshioplusplus::detail
