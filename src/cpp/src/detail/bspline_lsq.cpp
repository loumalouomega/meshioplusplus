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
// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

// Project includes (private, not installed)
#include "bspline_lsq.hpp"
#include "nurbs_eval.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kBslPi = 3.14159265358979323846;

double bsl_clamp(double X, double Lo, double Hi) {
    return std::min(std::max(X, Lo), Hi);
}

}  // namespace

std::vector<double> bspline_clamped_knots(int P, std::size_t N, double Lo, double Hi) {
    if (P < 1 || P > kBsplineMaxDegree || N < static_cast<std::size_t>(P) + 1 || !(Lo < Hi))
        throw std::invalid_argument(
            "meshio++: bspline_clamped_knots: needs a degree of 1 to 8, poles > degree and "
            "Lo < Hi");
    const auto p = static_cast<std::size_t>(P);
    std::vector<double> k(N + p + 1);
    const std::size_t spans = N - p;
    for (std::size_t i = 0; i <= p; ++i) {
        k[i] = Lo;
        k[N + i] = Hi;
    }
    for (std::size_t i = 1; i < spans; ++i)
        k[p + i] = Lo + (Hi - Lo) * static_cast<double>(i) / static_cast<double>(spans);
    return k;
}

BsplineSurface bspline_shape(int P, int Q, std::size_t Nu, std::size_t Nv, double U0, double U1,
                             double V0, double V1) {
    BsplineSurface s;
    s.mDegreeU = P;
    s.mDegreeV = Q;
    s.mNu = Nu;
    s.mNv = Nv;
    s.mKnotsU = bspline_clamped_knots(P, Nu, U0, U1);
    s.mKnotsV = bspline_clamped_knots(Q, Nv, V0, V1);
    s.mPoles.assign(Nu * Nv * 3, 0.0);
    return s;
}

void bspline_basis_ders(std::size_t Span, double U, int P, int NDer,
                        const std::vector<double>& rKnots, double* pDers) {
    // Fixed-size work arrays: no allocation on a path every projection takes.
    const auto p = static_cast<std::size_t>(P);
    const auto nd = static_cast<std::size_t>(std::max(NDer, 0));
    std::array<double, (kBsplineMaxDegree + 1) * (kBsplineMaxDegree + 1)> ndu{};
    std::array<double, kBsplineMaxDegree + 1> left{}, right{};
    const auto at = [p](std::size_t R, std::size_t C) { return R * (p + 1) + C; };
    ndu[at(0, 0)] = 1.0;
    for (std::size_t j = 1; j <= p; ++j) {
        left[j] = U - rKnots[Span + 1 - j];
        right[j] = rKnots[Span + j] - U;
        double saved = 0.0;
        for (std::size_t r = 0; r < j; ++r) {
            ndu[at(j, r)] = right[r + 1] + left[j - r];
            const double temp = ndu[at(r, j - 1)] / ndu[at(j, r)];
            ndu[at(r, j)] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[at(j, j)] = saved;
    }
    for (std::size_t k = 0; k <= nd; ++k)
        for (std::size_t j = 0; j <= p; ++j)
            pDers[k * (p + 1) + j] = 0.0;
    for (std::size_t j = 0; j <= p; ++j)
        pDers[j] = ndu[at(j, p)];
    const std::size_t n = std::min(nd, p);
    std::array<double, 2 * (kBsplineMaxDegree + 1)> a{};
    for (std::size_t r = 0; r <= p; ++r) {
        double* a0 = a.data();
        double* a1 = a.data() + (p + 1);
        a0[0] = 1.0;
        for (std::size_t k = 1; k <= n; ++k) {
            double d = 0.0;
            const auto rk = static_cast<long>(r) - static_cast<long>(k);
            const long pk = static_cast<long>(p) - static_cast<long>(k);
            if (rk >= 0) {
                a1[0] =
                    a0[0] / ndu[at(static_cast<std::size_t>(pk + 1), static_cast<std::size_t>(rk))];
                d = a1[0] * ndu[at(static_cast<std::size_t>(rk), static_cast<std::size_t>(pk))];
            }
            const long j1 = rk >= -1 ? 1 : -rk;
            const long j2 = (static_cast<long>(r) - 1 <= pk)
                                ? static_cast<long>(k) - 1
                                : static_cast<long>(p) - static_cast<long>(r);
            for (long j = j1; j <= j2; ++j) {
                const auto ju = static_cast<std::size_t>(j);
                a1[ju] =
                    (a0[ju] - a0[ju - 1]) /
                    ndu[at(static_cast<std::size_t>(pk + 1), static_cast<std::size_t>(rk + j))];
                d += a1[ju] *
                     ndu[at(static_cast<std::size_t>(rk + j), static_cast<std::size_t>(pk))];
            }
            if (static_cast<long>(r) <= pk) {
                a1[k] = -a0[k - 1] / ndu[at(static_cast<std::size_t>(pk + 1), r)];
                d += a1[k] * ndu[at(r, static_cast<std::size_t>(pk))];
            }
            pDers[k * (p + 1) + r] = d;
            std::swap(a0, a1);
        }
    }
    double factor = static_cast<double>(P);
    for (std::size_t k = 1; k <= n; ++k) {
        for (std::size_t j = 0; j <= p; ++j)
            pDers[k * (p + 1) + j] *= factor;
        factor *= static_cast<double>(P - static_cast<int>(k));
    }
}

void bspline_eval(const BsplineSurface& rS, double U, double V, int NDer, BsplineEval& rOut) {
    rOut = BsplineEval{};
    const auto p = static_cast<std::size_t>(rS.mDegreeU);
    const auto q = static_cast<std::size_t>(rS.mDegreeV);
    U = bsl_clamp(U, rS.U0(), rS.U1());
    V = bsl_clamp(V, rS.V0(), rS.V1());
    const std::size_t su = nurbs_find_span(rS.mDegreeU, rS.mKnotsU, rS.mNu, U);
    const std::size_t sv = nurbs_find_span(rS.mDegreeV, rS.mKnotsV, rS.mNv, V);
    const int nd = std::min(std::max(NDer, 0), 2);
    const auto ndz = static_cast<std::size_t>(nd);
    std::array<double, 3 * (kBsplineMaxDegree + 1)> bu{}, bv{};
    bspline_basis_ders(su, U, rS.mDegreeU, nd, rS.mKnotsU, bu.data());
    bspline_basis_ders(sv, V, rS.mDegreeV, nd, rS.mKnotsV, bv.data());
    double* outs[6] = {rOut.mS, rOut.mSu, rOut.mSv, rOut.mSuu, rOut.mSuv, rOut.mSvv};
    // (du order, dv order) of each output.
    static const std::size_t kOrd[6][2] = {{0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2}};
    for (std::size_t l = 0; l <= q; ++l) {
        const std::size_t j = sv - q + l;
        for (std::size_t k = 0; k <= p; ++k) {
            const std::size_t i = su - p + k;
            const double* pole = rS.mPoles.data() + 3 * (i + rS.mNu * j);
            for (std::size_t o = 0; o < 6; ++o) {
                if (kOrd[o][0] + kOrd[o][1] > ndz)
                    continue;
                const double w = bu[kOrd[o][0] * (p + 1) + k] * bv[kOrd[o][1] * (q + 1) + l];
                for (std::size_t c = 0; c < 3; ++c)
                    outs[o][c] += w * pole[c];
            }
        }
    }
}

void bspline_row(const BsplineSurface& rS, double U, double V, std::vector<std::size_t>& rCols,
                 std::vector<double>& rValues) {
    const auto p = static_cast<std::size_t>(rS.mDegreeU);
    const auto q = static_cast<std::size_t>(rS.mDegreeV);
    U = bsl_clamp(U, rS.U0(), rS.U1());
    V = bsl_clamp(V, rS.V0(), rS.V1());
    const std::size_t su = nurbs_find_span(rS.mDegreeU, rS.mKnotsU, rS.mNu, U);
    const std::size_t sv = nurbs_find_span(rS.mDegreeV, rS.mKnotsV, rS.mNv, V);
    std::array<double, kBsplineMaxDegree + 1> bu{}, bv{};
    bspline_basis_ders(su, U, rS.mDegreeU, 0, rS.mKnotsU, bu.data());
    bspline_basis_ders(sv, V, rS.mDegreeV, 0, rS.mKnotsV, bv.data());
    rCols.clear();
    rValues.clear();
    for (std::size_t l = 0; l <= q; ++l)
        for (std::size_t k = 0; k <= p; ++k) {
            rCols.push_back((su - p + k) + rS.mNu * (sv - q + l));
            rValues.push_back(bu[k] * bv[l]);
        }
}

SkylineMatrix bspline_normal_envelope(const BsplineSurface& rS) {
    const auto p = static_cast<std::size_t>(rS.mDegreeU);
    const auto q = static_cast<std::size_t>(rS.mDegreeV);
    std::vector<std::size_t> first(rS.NumPoles());
    for (std::size_t j = 0; j < rS.mNv; ++j)
        for (std::size_t i = 0; i < rS.mNu; ++i)
            first[i + rS.mNu * j] = (i > p ? i - p : 0) + rS.mNu * (j > q ? j - q : 0);
    return skyline_zero(first);
}

void gauss_legendre(int N, std::vector<double>& rNodes, std::vector<double>& rWeights) {
    const auto n = static_cast<std::size_t>(std::max(N, 1));
    rNodes.assign(n, 0.0);
    rWeights.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double x =
            std::cos(kBslPi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
        double dp = 1.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = x;
            for (std::size_t k = 2; k <= n; ++k) {
                const double kk = static_cast<double>(k);
                const double p2 = ((2.0 * kk - 1.0) * x * p1 - (kk - 1.0) * p0) / kk;
                p0 = p1;
                p1 = p2;
            }
            dp = static_cast<double>(n) * (x * p1 - p0) / (x * x - 1.0);
            const double dx = p1 / dp;
            x -= dx;
            if (std::abs(dx) < 1e-16)
                break;
        }
        rNodes[n - 1 - i] = x;
        rWeights[n - 1 - i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
}

void bspline_add_thin_plate_spans(const BsplineSurface& rS, const std::vector<double>& rSpanWeight,
                                  SkylineMatrix& rA) {
    const auto p = static_cast<std::size_t>(rS.mDegreeU);
    const auto q = static_cast<std::size_t>(rS.mDegreeV);
    // The integrand is a polynomial of degree 2p in u per span: p + 1 points are exact.
    std::vector<double> xu, wu, xv, wv;
    gauss_legendre(rS.mDegreeU + 1, xu, wu);
    gauss_legendre(rS.mDegreeV + 1, xv, wv);
    const std::vector<std::size_t> spans_u = nurbs_spans(rS.mDegreeU, rS.mKnotsU, rS.mNu);
    const std::vector<std::size_t> spans_v = nurbs_spans(rS.mDegreeV, rS.mKnotsV, rS.mNv);
    std::vector<double> bu(3 * (p + 1)), bv(3 * (q + 1));
    const std::size_t nloc = (p + 1) * (q + 1);
    std::vector<std::size_t> idx(nloc);
    std::vector<double> duu(nloc), duv(nloc), dvv(nloc);
    for (std::size_t jv = 0; jv < spans_v.size(); ++jv) {
        const std::size_t sv = spans_v[jv];
        const double v0 = rS.mKnotsV[sv], v1 = rS.mKnotsV[sv + 1];
        for (std::size_t ju = 0; ju < spans_u.size(); ++ju) {
            const std::size_t su = spans_u[ju];
            const double Weight = rSpanWeight[ju + spans_u.size() * jv];
            if (Weight == 0.0)
                continue;
            const double u0 = rS.mKnotsU[su], u1 = rS.mKnotsU[su + 1];
            const double jac = 0.25 * (u1 - u0) * (v1 - v0);
            for (std::size_t gv = 0; gv < xv.size(); ++gv) {
                const double v = 0.5 * (v0 + v1) + 0.5 * (v1 - v0) * xv[gv];
                bspline_basis_ders(sv, v, rS.mDegreeV, 2, rS.mKnotsV, bv.data());
                for (std::size_t gu = 0; gu < xu.size(); ++gu) {
                    const double u = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * xu[gu];
                    bspline_basis_ders(su, u, rS.mDegreeU, 2, rS.mKnotsU, bu.data());
                    const double w = Weight * jac * wu[gu] * wv[gv];
                    std::size_t m = 0;
                    for (std::size_t l = 0; l <= q; ++l)
                        for (std::size_t k = 0; k <= p; ++k, ++m) {
                            idx[m] = (su - p + k) + rS.mNu * (sv - q + l);
                            duu[m] = bu[2 * (p + 1) + k] * bv[l];
                            duv[m] = bu[(p + 1) + k] * bv[(q + 1) + l];
                            dvv[m] = bu[k] * bv[2 * (q + 1) + l];
                        }
                    for (std::size_t a = 0; a < nloc; ++a)
                        for (std::size_t b = 0; b <= a; ++b) {
                            const double e =
                                w * (duu[a] * duu[b] + 2.0 * duv[a] * duv[b] + dvv[a] * dvv[b]);
                            rA.At(idx[a], idx[b]) += e;
                        }
                }
            }
        }
    }
}

void bspline_add_thin_plate(const BsplineSurface& rS, double Weight, SkylineMatrix& rA) {
    if (Weight == 0.0)
        return;
    const std::size_t n = nurbs_spans(rS.mDegreeU, rS.mKnotsU, rS.mNu).size() *
                          nurbs_spans(rS.mDegreeV, rS.mKnotsV, rS.mNv).size();
    bspline_add_thin_plate_spans(rS, std::vector<double>(n, Weight), rA);
}

bool bspline_lsq_fit(BsplineSurface& rS, const std::vector<double>& rU,
                     const std::vector<double>& rV, const std::vector<double>& rXyz,
                     const std::vector<double>* pWeights, double Smoothing,
                     const std::vector<double>* pSpanSmoothing) {
    const std::size_t n = rS.NumPoles();
    SkylineMatrix a = bspline_normal_envelope(rS);
    std::vector<double> rhs(3 * n, 0.0);
    std::vector<std::size_t> cols;
    std::vector<double> vals;
    for (std::size_t s = 0; s < rU.size(); ++s) {
        const double w = pWeights ? (*pWeights)[s] : 1.0;
        if (w == 0.0)
            continue;
        bspline_row(rS, rU[s], rV[s], cols, vals);
        for (std::size_t a1 = 0; a1 < cols.size(); ++a1) {
            for (std::size_t b1 = 0; b1 < cols.size(); ++b1)
                if (cols[b1] <= cols[a1])
                    a.At(cols[a1], cols[b1]) += w * vals[a1] * vals[b1];
            for (std::size_t c = 0; c < 3; ++c)
                rhs[c * n + cols[a1]] += w * vals[a1] * rXyz[3 * s + c];
        }
    }
    if (pSpanSmoothing)
        bspline_add_thin_plate_spans(rS, *pSpanSmoothing, a);
    else
        bspline_add_thin_plate(rS, Smoothing, a);
    if (!skyline_ldlt(a))
        return false;
    for (std::size_t c = 0; c < 3; ++c)
        skyline_solve(a, rhs.data() + c * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t c = 0; c < 3; ++c) {
            const double x = rhs[c * n + i];
            if (!std::isfinite(x))
                return false;
        }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t c = 0; c < 3; ++c)
            rS.mPoles[3 * i + c] = rhs[c * n + i];
    return true;
}

double bspline_project(const BsplineSurface& rS, const double* pP, double& rU, double& rV,
                       int MaxIterations) {
    return bspline_project(rS, pP, rU, rV, rS.U0(), rS.U1(), rS.V0(), rS.V1(), MaxIterations);
}

double bspline_project(const BsplineSurface& rS, const double* pP, double& rU, double& rV,
                       double U0, double U1, double V0, double V1, int MaxIterations) {
    const double u0 = std::max(U0, rS.U0()), u1 = std::min(U1, rS.U1());
    const double v0 = std::max(V0, rS.V0()), v1 = std::min(V1, rS.V1());
    const double eps_u = 1e-14 * std::max(1.0, u1 - u0), eps_v = 1e-14 * std::max(1.0, v1 - v0);
    rU = bsl_clamp(rU, u0, u1);
    rV = bsl_clamp(rV, v0, v1);
    BsplineEval e;
    const auto dist = [&](double U, double V) {
        BsplineEval f;
        bspline_eval(rS, U, V, 0, f);
        double d2 = 0.0;
        for (int c = 0; c < 3; ++c)
            d2 += (f.mS[c] - pP[c]) * (f.mS[c] - pP[c]);
        return std::sqrt(d2);
    };
    double best = dist(rU, rV);
    for (int it = 0; it < MaxIterations; ++it) {
        bspline_eval(rS, rU, rV, 2, e);
        double r[3];
        for (int c = 0; c < 3; ++c)
            r[c] = e.mS[c] - pP[c];
        const auto dot = [](const double* pA, const double* pB) {
            return pA[0] * pB[0] + pA[1] * pB[1] + pA[2] * pB[2];
        };
        const double g0 = dot(r, e.mSu), g1 = dot(r, e.mSv);
        double h00 = dot(e.mSu, e.mSu) + dot(r, e.mSuu);
        double h01 = dot(e.mSu, e.mSv) + dot(r, e.mSuv);
        double h11 = dot(e.mSv, e.mSv) + dot(r, e.mSvv);
        double det = h00 * h11 - h01 * h01;
        if (!(h00 > 0.0) || !(det > 0.0)) {
            // Gauss-Newton: the first fundamental form, positive definite on a
            // regular surface.
            h00 = dot(e.mSu, e.mSu);
            h01 = dot(e.mSu, e.mSv);
            h11 = dot(e.mSv, e.mSv);
            det = h00 * h11 - h01 * h01;
            if (!(det > 0.0))
                break;
        }
        double du = -(h11 * g0 - h01 * g1) / det;
        double dv = -(h00 * g1 - h01 * g0) / det;
        // Halve the step until the distance does not grow (at most 30 times).
        double nu = bsl_clamp(rU + du, u0, u1), nv = bsl_clamp(rV + dv, v0, v1);
        double d = dist(nu, nv);
        for (int h = 0; h < 30 && d > best; ++h) {
            du *= 0.5;
            dv *= 0.5;
            nu = bsl_clamp(rU + du, u0, u1);
            nv = bsl_clamp(rV + dv, v0, v1);
            d = dist(nu, nv);
        }
        if (d > best)
            break;
        const double step_u = std::abs(nu - rU), step_v = std::abs(nv - rV);
        rU = nu;
        rV = nv;
        best = d;
        if (step_u <= eps_u && step_v <= eps_v)
            break;
    }
    return best;
}

}  // namespace detail
}  // namespace meshioplusplus
