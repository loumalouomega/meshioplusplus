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
 * @file detail/bspline_lsq.hpp
 * @brief The OpenCASCADE-free pieces of a tensor-product B-spline surface
 * least-squares fit: clamped knot vectors, basis functions with derivatives,
 * point evaluation with first and second derivatives, the collocation row of
 * a sample, a thin-plate smoothing term integrated by Gauss-Legendre
 * quadrature, the regularised least-squares solve of the poles, and the
 * Newton projection of a point onto the surface (parameter correction).
 *
 * The surface is polynomial (no weights), `mNu x mNv` poles listed `u`
 * fastest (pole `(i, j)` at index `i + mNu * j`), three coordinates each, with
 * full clamped knot vectors (`n + p + 1` knots for `n` poles). It reuses
 * `nurbs_find_span` from `nurbs_eval.hpp`; the derivatives are The NURBS
 * Book's A2.3 to any order (`iga_kernels.hpp` stops at the first).
 *
 * A **core-private** header (the `lsq_fit.hpp` precedent): nothing installed
 * names it. Serial and deterministic: fixed loop orders, fixed iteration caps.
 */

// System includes
#include <cstddef>
#include <vector>

// Project includes (private, not installed)
#include "dense_chol.hpp"

namespace meshioplusplus {
namespace detail {

/// The highest degree these routines take (`fit_cad`'s freeform limit).
inline constexpr int kBsplineMaxDegree = 8;

/// A polynomial tensor-product B-spline surface in 3D.
struct BsplineSurface {
    int mDegreeU = 3;
    int mDegreeV = 3;
    std::vector<double> mKnotsU;
    std::vector<double> mKnotsV;
    std::size_t mNu = 0;
    std::size_t mNv = 0;
    /// `mNu * mNv * 3` coordinates, pole `(i, j)` at `3 * (i + mNu * j)`.
    std::vector<double> mPoles;

    std::size_t NumPoles() const { return mNu * mNv; }
    /// The parameter domain `[U0, U1] x [V0, V1]`.
    double U0() const { return mKnotsU[static_cast<std::size_t>(mDegreeU)]; }
    double U1() const { return mKnotsU[mNu]; }
    double V0() const { return mKnotsV[static_cast<std::size_t>(mDegreeV)]; }
    double V1() const { return mKnotsV[mNv]; }
};

/**
 * @brief A clamped knot vector for `N` poles of degree `P` on `[Lo, Hi]`:
 * `P + 1` copies of each end and `N - P - 1` interior knots evenly spaced.
 * @throws std::invalid_argument unless `1 <= P <= kBsplineMaxDegree`,
 * `N >= P + 1` and `Lo < Hi`.
 */
std::vector<double> bspline_clamped_knots(int P, std::size_t N, double Lo, double Hi);

/**
 * @brief A surface shape with clamped uniform knots on the given domain and
 * zero poles.
 */
BsplineSurface bspline_shape(int P, int Q, std::size_t Nu, std::size_t Nv, double U0, double U1,
                             double V0, double V1);

/**
 * @brief The basis functions of span `Span` at `U` and their derivatives up to
 * order `NDer` (A2.3): `pDers` is `(NDer + 1)` rows of `P + 1` values, row `k`
 * the `k`-th derivatives (zero for `k > P`).
 */
void bspline_basis_ders(std::size_t Span, double U, int P, int NDer,
                        const std::vector<double>& rKnots, double* pDers);

/// A surface point and its partial derivatives.
struct BsplineEval {
    double mS[3] = {0.0, 0.0, 0.0};
    double mSu[3] = {0.0, 0.0, 0.0};
    double mSv[3] = {0.0, 0.0, 0.0};
    double mSuu[3] = {0.0, 0.0, 0.0};
    double mSuv[3] = {0.0, 0.0, 0.0};
    double mSvv[3] = {0.0, 0.0, 0.0};
};

/**
 * @brief Evaluate the surface at `(U, V)` (clamped to the domain) with
 * derivatives up to order `NDer` (0, 1 or 2; higher ones stay zero).
 */
void bspline_eval(const BsplineSurface& rS, double U, double V, int NDer, BsplineEval& rOut);

/**
 * @brief The collocation row of a sample at `(U, V)`: the `(p + 1)(q + 1)`
 * pole indices whose basis functions are non-zero there, and those values.
 */
void bspline_row(const BsplineSurface& rS, double U, double V, std::vector<std::size_t>& rCols,
                 std::vector<double>& rValues);

/**
 * @brief The envelope of the normal matrix: pole `(i, j)` couples with every
 * pole within `p` columns and `q` rows of it, so row `i + Nu j` starts at
 * `max(0, i - p) + Nu max(0, j - q)`.
 */
SkylineMatrix bspline_normal_envelope(const BsplineSurface& rS);

/**
 * @brief Add `Weight` times the thin-plate energy
 * `∫∫ |S_uu|² + 2 |S_uv|² + |S_vv|² du dv` (a quadratic form in each
 * coordinate of the poles, the same for all three) to `rA`, integrated
 * exactly by Gauss-Legendre quadrature over every knot span.
 */
void bspline_add_thin_plate(const BsplineSurface& rS, double Weight, SkylineMatrix& rA);

/**
 * @brief The same with a weight per knot span: `rSpanWeight` holds one value
 * per pair of non-empty spans, `u` fastest (`nurbs_spans` order in each
 * direction), so a fit can bend freely where its samples are and keep its
 * extension past them smooth.
 */
void bspline_add_thin_plate_spans(const BsplineSurface& rS, const std::vector<double>& rSpanWeight,
                                  SkylineMatrix& rA);

/**
 * @brief Fit the poles of `rS` (its degrees, knots and pole counts given) to
 * scattered samples by regularised least squares:
 * `min Σ w_k |S(u_k, v_k) - x_k|² + Smoothing · thin-plate energy`.
 * @param rU, rV the sample parameters (clamped to the domain).
 * @param rXyz `3 * rU.size()` coordinates.
 * @param pWeights per sample, or null for all ones.
 * @param pSpanSmoothing when given, per knot span its thin-plate weight
 *        (see `bspline_add_thin_plate_spans`), replacing `Smoothing`.
 * @return false (poles left as they were) when the normal matrix is singular:
 *         poles no sample touches and no smoothing to hold them.
 */
bool bspline_lsq_fit(BsplineSurface& rS, const std::vector<double>& rU,
                     const std::vector<double>& rV, const std::vector<double>& rXyz,
                     const std::vector<double>* pWeights, double Smoothing,
                     const std::vector<double>* pSpanSmoothing = nullptr);

/**
 * @brief Project `pP` onto the surface by Newton's method on the
 * closest-point conditions `(S - P)·S_u = (S - P)·S_v = 0`, falling back to
 * a Gauss-Newton step where the Hessian is not positive definite, the
 * parameters clamped to the domain.
 * @param rU, rV the start on entry, the foot point's parameters on exit.
 * @return the distance from `pP` to the surface point at the final `(U, V)`.
 */
double bspline_project(const BsplineSurface& rS, const double* pP, double& rU, double& rV,
                       int MaxIterations = 20);

/// The same with the foot point's parameters held inside `[U0, U1] x [V0, V1]`
/// (a part of the domain, such as where a fit's samples lie).
double bspline_project(const BsplineSurface& rS, const double* pP, double& rU, double& rV,
                       double U0, double U1, double V0, double V1, int MaxIterations = 20);

/**
 * @brief The `N` Gauss-Legendre nodes and weights on `[-1, 1]` (Newton on the
 * Legendre polynomial), exact for polynomials of degree `2N - 1`.
 */
void gauss_legendre(int N, std::vector<double>& rNodes, std::vector<double>& rWeights);

}  // namespace detail
}  // namespace meshioplusplus
