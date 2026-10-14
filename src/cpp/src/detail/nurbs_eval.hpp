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
 * @file detail/nurbs_eval.hpp
 * @brief Point evaluation of B-spline and NURBS curves and tensor-product
 * surfaces at arbitrary parameters (The NURBS Book, algorithms A2.1, A2.2,
 * A3.1/A4.1 and A3.5/A4.3).
 *
 * Knot vectors are in the full form, `n + p + 1` knots for `n` control
 * points. Control nets list direction 0 fastest (Kratos's and MFEM's order),
 * each point `Dim` doubles. A null weight pointer evaluates the polynomial
 * B-spline, with no division; otherwise the rational map. The arithmetic is
 * the order the Python twin (`meshioplusplus/_nurbs_eval.py`) follows, so
 * the two give the same bits.
 *
 * Core-private (not installed). The MFEM reader keeps its own kernel, which
 * works on MFEM's compressed knot spans.
 */

// System includes
#include <cstddef>
#include <utility>
#include <vector>

namespace meshioplusplus::detail {

/// The knot span index `i` (`p <= i <= n - 1`) with `U[i] <= u < U[i+1]`;
/// the last span for `u >= U[n]`, the first for `u <= U[p]`.
inline std::size_t nurbs_find_span(int P, const std::vector<double>& rKnots, std::size_t N,
                                   double U) {
    const auto p = static_cast<std::size_t>(P);
    const std::size_t n = N - 1;
    if (U >= rKnots[n + 1]) {
        std::size_t s = n;
        while (s > p && rKnots[s] == rKnots[s + 1])
            --s;
        return s;
    }
    if (U <= rKnots[p]) {
        std::size_t s = p;
        while (s < n && rKnots[s] == rKnots[s + 1])
            ++s;
        return s;
    }
    std::size_t low = p, high = n + 1, mid = (low + high) / 2;
    while (U < rKnots[mid] || U >= rKnots[mid + 1]) {
        if (U < rKnots[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

/// The `p + 1` nonzero basis functions of span `Span` at `U` into `pN`.
inline void nurbs_basis(std::size_t Span, double U, int P, const std::vector<double>& rKnots,
                        double* pN) {
    const auto p = static_cast<std::size_t>(P);
    std::vector<double> left(p + 1), right(p + 1);
    pN[0] = 1.0;
    for (std::size_t j = 1; j <= p; ++j) {
        left[j] = U - rKnots[Span + 1 - j];
        right[j] = rKnots[Span + j] - U;
        double saved = 0.0;
        for (std::size_t r = 0; r < j; ++r) {
            const double temp = pN[r] / (right[r + 1] + left[j - r]);
            pN[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        pN[j] = saved;
    }
}

/// The spans of nonzero length as `(knot index, knot index + 1)` pairs:
/// `i` from `p` to `n - 1` with `U[i] < U[i+1]`.
inline std::vector<std::size_t> nurbs_spans(int P, const std::vector<double>& rKnots,
                                            std::size_t N) {
    std::vector<std::size_t> spans;
    for (std::size_t i = static_cast<std::size_t>(P); i + 1 <= N; ++i)
        if (rKnots[i] < rKnots[i + 1])
            spans.push_back(i);
    return spans;
}

/// A curve point at `U` in span `Span` (from `nurbs_find_span` or a cell's
/// own span): `pOut[0 .. Dim)`.
inline void nurbs_curve_point_in_span(int P, const std::vector<double>& rKnots,
                                      const double* pCps, const double* pW, std::size_t Dim,
                                      std::size_t Span, double U, double* pOut) {
    const auto p = static_cast<std::size_t>(P);
    std::vector<double> nb(p + 1);
    nurbs_basis(Span, U, P, rKnots, nb.data());
    for (std::size_t d = 0; d < Dim; ++d)
        pOut[d] = 0.0;
    double w = 0.0;
    for (std::size_t k = 0; k <= p; ++k) {
        const std::size_t i = Span - p + k;
        const double c = pW ? nb[k] * pW[i] : nb[k];
        for (std::size_t d = 0; d < Dim; ++d)
            pOut[d] += c * pCps[i * Dim + d];
        w += c;
    }
    if (pW)
        for (std::size_t d = 0; d < Dim; ++d)
            pOut[d] /= w;
}

/// A surface point at `(U, V)` in spans `(SpanU, SpanV)`.
inline void nurbs_surface_point_in_span(int P, int Q, const std::vector<double>& rKnotsU,
                                        const std::vector<double>& rKnotsV, const double* pCps,
                                        const double* pW, std::size_t Nu, std::size_t Dim,
                                        std::size_t SpanU, std::size_t SpanV, double U, double V,
                                        double* pOut) {
    const auto p = static_cast<std::size_t>(P);
    const auto q = static_cast<std::size_t>(Q);
    std::vector<double> nu(p + 1), nv(q + 1), temp(Dim + 1);
    nurbs_basis(SpanU, U, P, rKnotsU, nu.data());
    nurbs_basis(SpanV, V, Q, rKnotsV, nv.data());
    for (std::size_t d = 0; d < Dim; ++d)
        pOut[d] = 0.0;
    double w = 0.0;
    for (std::size_t l = 0; l <= q; ++l) {
        for (std::size_t d = 0; d <= Dim; ++d)
            temp[d] = 0.0;
        const std::size_t j = SpanV - q + l;
        for (std::size_t k = 0; k <= p; ++k) {
            const std::size_t i = SpanU - p + k;
            const std::size_t idx = j * Nu + i;
            const double c = pW ? nu[k] * pW[idx] : nu[k];
            for (std::size_t d = 0; d < Dim; ++d)
                temp[d] += c * pCps[idx * Dim + d];
            temp[Dim] += c;
        }
        for (std::size_t d = 0; d < Dim; ++d)
            pOut[d] += nv[l] * temp[d];
        w += nv[l] * temp[Dim];
    }
    if (pW)
        for (std::size_t d = 0; d < Dim; ++d)
            pOut[d] /= w;
}

}  // namespace meshioplusplus::detail
