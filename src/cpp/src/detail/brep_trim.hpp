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
 * @file detail/brep_trim.hpp
 * @brief The parametric helpers a trimmed B-rep face is sampled and evaluated
 * with, shared by the CAD JSON sampler (`formats/cad_json.cpp`) and the IGA
 * layer (`operations/iga.cpp`): the cells of one knot-vector direction, their
 * lattice, a trimming loop as a polygon in `(u, v)`, and the even-odd
 * inside test over all loops of a face.
 *
 * The Python twins are `_axis`, `_lattice`, `_loop_polygon` and `_inside` in
 * `meshioplusplus/cad_json/_cad_json.py`, which repeat these bit for bit.
 *
 * Core-private (not installed).
 */

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"

// Project includes (private, not installed)
#include "nurbs_eval.hpp"

namespace meshioplusplus::detail {

/// The largest degree a B-rep curve or surface may have.
inline constexpr int kBrepMaxDegree = 32;

/// Check a curve's degree, control points, weights, knots and ids; the error
/// text is `pPrefix + rWhat + " ..."`.
inline void brep_check_curve(const BrepCurve& rC, const std::string& rWhat, const char* pPrefix) {
    const auto bad = [&](const char* pText) {
        throw std::invalid_argument(std::string(pPrefix) + rWhat + pText);
    };
    const std::size_t n = rC.mWeights.Size();
    if (rC.mDegree < 1 || rC.mDegree > kBrepMaxDegree)
        bad(" degree is out of range");
    if (rC.mControlPoints.Dtype() != DType::Float64 || rC.mControlPoints.Shape().size() != 2 ||
        rC.mControlPoints.Shape()[0] != n || rC.mControlPoints.Shape()[1] != 3 ||
        rC.mWeights.Dtype() != DType::Float64)
        bad(" control points are not Float64 (n, 3) beside n weights");
    if (n < static_cast<std::size_t>(rC.mDegree) + 1 ||
        rC.mKnots.size() != n + static_cast<std::size_t>(rC.mDegree) + 1)
        bad(" knot vector does not match its control points");
    for (std::size_t i = 1; i < rC.mKnots.size(); ++i)
        if (!(rC.mKnots[i] >= rC.mKnots[i - 1]))
            bad(" knot vector decreases");
    if (rC.mIds.Size() != 0 && (rC.mIds.Size() != n || rC.mIds.Dtype() != DType::Int64))
        bad(" ids do not match its control points");
}

/// The cell boundaries along one direction: every nonzero span split in `S`.
struct BrepAxis {
    std::vector<double> mLo, mHi;       // per cell
    std::vector<std::size_t> mSpan;     // knot index of each cell's span
    std::vector<std::int64_t> mSpanNo;  // the span's ordinal among nonzero spans
};

/// The cells of one direction (`N` control points), clipped to `pRange` when it
/// holds two values.
inline BrepAxis brep_axis(int P, const std::vector<double>& rKnots, std::size_t N, int S,
                          const std::vector<double>* pRange) {
    BrepAxis a;
    const std::vector<std::size_t> spans = nurbs_spans(P, rKnots, N);
    for (std::size_t s = 0; s < spans.size(); ++s) {
        double k0 = rKnots[spans[s]], k1 = rKnots[spans[s] + 1];
        if (pRange && pRange->size() == 2) {
            const double r0 = std::min((*pRange)[0], (*pRange)[1]);
            const double r1 = std::max((*pRange)[0], (*pRange)[1]);
            k0 = std::max(k0, r0);
            k1 = std::min(k1, r1);
            if (!(k0 < k1))
                continue;
        }
        for (int j = 0; j < S; ++j) {
            const double lo = k0 + (k1 - k0) * (static_cast<double>(j) / S);
            const double hi = j + 1 == S ? k1 : k0 + (k1 - k0) * (static_cast<double>(j + 1) / S);
            a.mLo.push_back(lo);
            a.mHi.push_back(hi);
            a.mSpan.push_back(spans[s]);
            a.mSpanNo.push_back(static_cast<std::int64_t>(s));
        }
    }
    return a;
}

/// The lattice coordinates of one axis: cell c's node i at index `c * K + i`.
inline std::vector<double> brep_lattice(const BrepAxis& rA, int K) {
    const std::size_t n = rA.mLo.size();
    std::vector<double> g(n * static_cast<std::size_t>(K) + 1, 0.0);
    for (std::size_t c = 0; c < n; ++c)
        for (int i = 0; i < K; ++i)
            g[c * static_cast<std::size_t>(K) + static_cast<std::size_t>(i)] =
                rA.mLo[c] + (rA.mHi[c] - rA.mLo[c]) * (static_cast<double>(i) / K);
    if (n)
        g[n * static_cast<std::size_t>(K)] = rA.mHi[n - 1];
    return g;
}

/// A loop as a closed polygon in `(u, v)`: every trim sampled at `Samples`
/// points per knot span of its active range, in loop order.
inline std::vector<std::array<double, 2>> brep_loop_polygon(const BrepLoop& rLoop, int Samples) {
    std::vector<std::array<double, 2>> poly;
    for (const BrepTrim& t : rLoop.mTrims) {
        const BrepCurve& c = t.mCurve;
        const std::size_t n = c.mWeights.Size();
        const BrepAxis axis = brep_axis(c.mDegree, c.mKnots, n, 1, &c.mActiveRange);
        std::vector<std::array<double, 2>> pts;
        const double* w = c.mRational ? c.mWeights.As<double>() : nullptr;
        for (std::size_t s = 0; s < axis.mLo.size(); ++s)
            for (int j = 0; j < Samples; ++j) {
                const double u =
                    axis.mLo[s] + (axis.mHi[s] - axis.mLo[s]) * (static_cast<double>(j) / Samples);
                double xyz[3];
                nurbs_curve_point_in_span(c.mDegree, c.mKnots, c.mControlPoints.As<double>(), w, 3,
                                          axis.mSpan[s], u, xyz);
                pts.push_back({xyz[0], xyz[1]});
            }
        if (!axis.mLo.empty()) {
            const std::size_t s = axis.mLo.size() - 1;
            double xyz[3];
            nurbs_curve_point_in_span(c.mDegree, c.mKnots, c.mControlPoints.As<double>(), w, 3,
                                      axis.mSpan[s], axis.mHi[s], xyz);
            pts.push_back({xyz[0], xyz[1]});
        }
        if (!t.mDirection)
            std::reverse(pts.begin(), pts.end());
        poly.insert(poly.end(), pts.begin(), pts.end());
    }
    return poly;
}

/// Even-odd point-in-polygon over every loop of a face.
inline bool brep_inside(const std::vector<std::vector<std::array<double, 2>>>& rLoops, double U,
                        double V) {
    bool inside = false;
    for (const auto& poly : rLoops) {
        const std::size_t n = poly.size();
        for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
            const auto& a = poly[i];
            const auto& b = poly[j];
            if ((a[1] > V) != (b[1] > V) && U < (b[0] - a[0]) * (V - a[1]) / (b[1] - a[1]) + a[0])
                inside = !inside;
        }
    }
    return inside;
}

}  // namespace meshioplusplus::detail
