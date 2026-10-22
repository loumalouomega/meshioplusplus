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
#include <cmath>

// Project includes (private, not installed)
#include "bspline_curve.hpp"
#include "bspline_lsq.hpp"
#include "lsq_fit.hpp"
#include "nurbs_eval.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

double bc_dist(const PatchVec& rA, const PatchVec& rB) {
    const double x = rA[0] - rB[0], y = rA[1] - rB[1], z = rA[2] - rB[2];
    return std::sqrt(x * x + y * y + z * z);
}

// The basis values (and first derivatives) of span `Span` at `T`.
void bc_basis(const BsplineCurve& rC, std::size_t Span, double T, double* pN, double* pD) {
    const int p = rC.mDegree;
    double ders[2 * (kBsplineMaxDegree + 1)];
    bspline_basis_ders(Span, T, p, 1, rC.mKnots, ders);
    for (int k = 0; k <= p; ++k) {
        pN[k] = ders[k];
        if (pD)
            pD[k] = ders[p + 1 + k];
    }
}

// `T`'s multiplicity in the knot vector.
int bc_multiplicity(const std::vector<double>& rKnots, double T) {
    int m = 0;
    for (double k : rKnots)
        m += k == T ? 1 : 0;
    return m;
}

// Split at an interior knot of full multiplicity (degree + 1): the left and
// right curves.
void bc_split_at(const BsplineCurve& rC, double T, BsplineCurve& rLeft, BsplineCurve& rRight) {
    const std::size_t p = static_cast<std::size_t>(rC.mDegree);
    std::size_t r0 = 0;
    while (rC.mKnots[r0] != T)
        ++r0;
    rLeft.mDegree = rRight.mDegree = rC.mDegree;
    rLeft.mKnots.assign(rC.mKnots.begin(),
                        rC.mKnots.begin() + static_cast<std::ptrdiff_t>(r0 + p + 1));
    rLeft.mPoles.assign(rC.mPoles.begin(), rC.mPoles.begin() + static_cast<std::ptrdiff_t>(r0));
    rRight.mKnots.assign(rC.mKnots.begin() + static_cast<std::ptrdiff_t>(r0), rC.mKnots.end());
    rRight.mPoles.assign(rC.mPoles.begin() + static_cast<std::ptrdiff_t>(r0), rC.mPoles.end());
}

}  // namespace

PatchVec bspline_curve_eval(const BsplineCurve& rC, double T, PatchVec* pD) {
    const std::size_t n = rC.mPoles.size();
    T = std::max(rC.T0(), std::min(rC.T1(), T));
    const std::size_t span = nurbs_find_span(rC.mDegree, rC.mKnots, n, T);
    double nb[kBsplineMaxDegree + 1], db[kBsplineMaxDegree + 1];
    bc_basis(rC, span, T, nb, db);
    PatchVec out{0.0, 0.0, 0.0}, d{0.0, 0.0, 0.0};
    const std::size_t p = static_cast<std::size_t>(rC.mDegree);
    for (std::size_t k = 0; k <= p; ++k) {
        const PatchVec& q = rC.mPoles[span - p + k];
        for (std::size_t j = 0; j < 3; ++j) {
            out[j] += nb[k] * q[j];
            d[j] += db[k] * q[j];
        }
    }
    if (pD)
        *pD = d;
    return out;
}

void bspline_curve_insert_knot(BsplineCurve& rC, double T, int Times) {
    const int p = rC.mDegree;
    for (int r = 0; r < Times; ++r) {
        const std::size_t n = rC.mPoles.size();
        const int s = bc_multiplicity(rC.mKnots, T);
        if (s > p)
            return;
        if (s == p) {
            // Full multiplicity: duplicate the pole the curve passes through.
            std::size_t r0 = 0;
            while (rC.mKnots[r0] != T)
                ++r0;
            rC.mKnots.insert(rC.mKnots.begin() + static_cast<std::ptrdiff_t>(r0), T);
            rC.mPoles.insert(rC.mPoles.begin() + static_cast<std::ptrdiff_t>(r0),
                             rC.mPoles[r0 - 1]);
            continue;
        }
        // Boehm: one insertion in span k.
        std::size_t k = 0;
        while (k + 1 < rC.mKnots.size() && rC.mKnots[k + 1] <= T)
            ++k;
        k = std::min(k, n - 1);
        std::vector<PatchVec> q(n + 1);
        const std::size_t pp = static_cast<std::size_t>(p);
        for (std::size_t i = 0; i <= k - pp; ++i)
            q[i] = rC.mPoles[i];
        for (std::size_t i = k + 1 - static_cast<std::size_t>(s); i <= n; ++i)
            q[i] = rC.mPoles[i - 1];
        for (std::size_t i = k - pp + 1; i <= k - static_cast<std::size_t>(s); ++i) {
            const double a = (T - rC.mKnots[i]) / (rC.mKnots[i + pp] - rC.mKnots[i]);
            for (std::size_t j = 0; j < 3; ++j)
                q[i][j] = a * rC.mPoles[i][j] + (1.0 - a) * rC.mPoles[i - 1][j];
        }
        rC.mKnots.insert(rC.mKnots.begin() + static_cast<std::ptrdiff_t>(k + 1), T);
        rC.mPoles = std::move(q);
    }
}

void bspline_curve_reparam(BsplineCurve& rC, double T0, double T1) {
    const double a = rC.T0(), b = rC.T1();
    for (double& k : rC.mKnots)
        k = T0 + (T1 - T0) * (k - a) / (b - a);
    // Keep the ends exact.
    const std::size_t p = static_cast<std::size_t>(rC.mDegree);
    for (std::size_t i = 0; i <= p; ++i) {
        rC.mKnots[i] = T0;
        rC.mKnots[rC.mKnots.size() - 1 - i] = T1;
    }
}

BsplineCurve bspline_curve_reverse(const BsplineCurve& rC) {
    BsplineCurve out;
    out.mDegree = rC.mDegree;
    const double a = rC.T0(), b = rC.T1();
    out.mKnots.resize(rC.mKnots.size());
    for (std::size_t i = 0; i < rC.mKnots.size(); ++i)
        out.mKnots[i] = a + b - rC.mKnots[rC.mKnots.size() - 1 - i];
    out.mPoles.assign(rC.mPoles.rbegin(), rC.mPoles.rend());
    const std::size_t p = static_cast<std::size_t>(rC.mDegree);
    for (std::size_t i = 0; i <= p; ++i) {
        out.mKnots[i] = a;
        out.mKnots[out.mKnots.size() - 1 - i] = b;
    }
    return out;
}

BsplineCurve bspline_curve_segment(const BsplineCurve& rC, double A, double B) {
    if (A > B) {
        BsplineCurve s = bspline_curve_segment(rC, B, A);
        return bspline_curve_reverse(s);
    }
    BsplineCurve c = rC;
    const int p = c.mDegree;
    A = std::max(c.T0(), A);
    B = std::min(c.T1(), B);
    if (A > c.T0()) {
        bspline_curve_insert_knot(c, A, p + 1 - bc_multiplicity(c.mKnots, A));
        BsplineCurve l, r;
        bc_split_at(c, A, l, r);
        c = std::move(r);
    }
    if (B < c.T1()) {
        bspline_curve_insert_knot(c, B, p + 1 - bc_multiplicity(c.mKnots, B));
        BsplineCurve l, r;
        bc_split_at(c, B, l, r);
        c = std::move(l);
    }
    bspline_curve_reparam(c, 0.0, 1.0);
    return c;
}

BsplineCurve bspline_curve_join(const std::vector<BsplineCurve>& rParts,
                                const std::vector<double>& rLengths) {
    if (rParts.size() == 1) {
        BsplineCurve c = rParts[0];
        bspline_curve_reparam(c, 0.0, 1.0);
        return c;
    }
    double total = 0.0;
    for (double l : rLengths)
        total += l;
    BsplineCurve out;
    const int p = rParts[0].mDegree;
    const std::size_t pp = static_cast<std::size_t>(p);
    out.mDegree = p;
    double at = 0.0;
    for (std::size_t i = 0; i < rParts.size(); ++i) {
        BsplineCurve c = rParts[i];
        const double next = i + 1 == rParts.size() ? 1.0 : at + rLengths[i] / total;
        bspline_curve_reparam(c, at, next);
        if (i == 0) {
            out.mKnots.assign(c.mKnots.begin(), c.mKnots.end() - 1);  // end at multiplicity p
            out.mPoles = c.mPoles;
        } else {
            // Average the shared end, then append the rest.
            PatchVec& last = out.mPoles.back();
            for (std::size_t j = 0; j < 3; ++j)
                last[j] = 0.5 * (last[j] + c.mPoles[0][j]);
            const bool final_part = i + 1 == rParts.size();
            out.mKnots.insert(out.mKnots.end(),
                              c.mKnots.begin() + static_cast<std::ptrdiff_t>(pp + 1),
                              final_part ? c.mKnots.end() : c.mKnots.end() - 1);
            out.mPoles.insert(out.mPoles.end(), c.mPoles.begin() + 1, c.mPoles.end());
        }
        at = next;
    }
    return out;
}

std::vector<double> bspline_knot_union(const std::vector<double>& rA,
                                       const std::vector<double>& rB) {
    std::vector<double> out;
    std::size_t i = 0, j = 0;
    while (i < rA.size() || j < rB.size()) {
        if (j >= rB.size() || (i < rA.size() && rA[i] < rB[j])) {
            out.push_back(rA[i++]);
        } else if (i >= rA.size() || rB[j] < rA[i]) {
            out.push_back(rB[j++]);
        } else {
            out.push_back(rA[i]);
            ++i;
            ++j;
        }
    }
    return out;
}

void bspline_curve_refine_to(BsplineCurve& rC, const std::vector<double>& rKnots) {
    std::size_t i = 0;
    while (i < rKnots.size()) {
        std::size_t j = i;
        while (j < rKnots.size() && rKnots[j] == rKnots[i])
            ++j;
        const int want = static_cast<int>(j - i);
        const int have = bc_multiplicity(rC.mKnots, rKnots[i]);
        if (want > have && rKnots[i] > rC.T0() && rKnots[i] < rC.T1())
            bspline_curve_insert_knot(rC, rKnots[i], want - have);
        i = j;
    }
}

CurveFitResult bspline_curve_fit(const std::vector<PatchVec>& rPts,
                                 const CurveFitOptions& rOptions) {
    CurveFitResult out;
    const std::size_t m = rPts.size();
    const int p = rOptions.mDegree;
    if (m < 2 || p < 1 || p > kBsplineMaxDegree)
        return out;
    const bool aniso = rOptions.mNormals.size() == m;
    // Chord-length parameters.
    std::vector<double> t(m, 0.0);
    for (std::size_t k = 1; k < m; ++k)
        t[k] = t[k - 1] + bc_dist(rPts[k], rPts[k - 1]);
    const double len = t.back();
    if (!(len > 0.0))
        return out;
    for (double& v : t)
        v /= len;
    t.back() = 1.0;
    // The points held exactly: the ends and those asked for.
    std::vector<std::size_t> held{0};
    for (std::size_t k : rOptions.mInterpolate)
        if (k > 0 && k + 1 < m && k > held.back())
            held.push_back(k);
    held.push_back(m - 1);
    const std::size_t nc = held.size();
    const std::size_t min_poles = std::max<std::size_t>(static_cast<std::size_t>(p) + 1, nc);
    const std::size_t cap = std::max(min_poles, std::min(rOptions.mMaxPoles, m + nc));
    std::size_t n = std::min(cap, std::max(min_poles, rOptions.mStartPoles));
    // A point's residual measure: along its normal, and lightly across, when
    // normals are given; plainly otherwise.
    const auto metric = [&](std::size_t K, double* pM) {
        for (int i = 0; i < 9; ++i)
            pM[i] = 0.0;
        if (!aniso) {
            pM[0] = pM[4] = pM[8] = 1.0;
            return;
        }
        const PatchVec& q = rOptions.mNormals[K];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j)
                pM[3 * i + j] = q[static_cast<std::size_t>(i)] * q[static_cast<std::size_t>(j)];
            pM[3 * i + i] += rOptions.mTangentWeight;
        }
    };
    const auto distance = [&](std::size_t K, const PatchVec& rAt) {
        const PatchVec r{rAt[0] - rPts[K][0], rAt[1] - rPts[K][1], rAt[2] - rPts[K][2]};
        if (!aniso)
            return std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        const PatchVec& q = rOptions.mNormals[K];
        return std::fabs(r[0] * q[0] + r[1] * q[1] + r[2] * q[2]);
    };
    for (;;) {
        BsplineCurve c;
        c.mDegree = p;
        c.mKnots = bspline_clamped_knots(p, n, 0.0, 1.0);
        c.mPoles.assign(n, PatchVec{0.0, 0.0, 0.0});
        std::vector<double> par = t;
        bool solved = false;
        for (int pass = 0; pass <= rOptions.mParameterCorrections; ++pass) {
            // Normal equations (3 x 3 blocks per pole pair) plus smoothing,
            // bordered by the constraints (3 rows per held point).
            const std::size_t dim = 3 * n + 3 * nc;
            std::vector<double> a(dim * dim, 0.0), rhs(dim, 0.0);
            double nb[kBsplineMaxDegree + 1];
            double mk[9];
            for (std::size_t k = 0; k < m; ++k) {
                const std::size_t span = nurbs_find_span(p, c.mKnots, n, par[k]);
                bc_basis(c, span, par[k], nb, nullptr);
                metric(k, mk);
                for (int x = 0; x <= p; ++x) {
                    const std::size_t i =
                        span - static_cast<std::size_t>(p) + static_cast<std::size_t>(x);
                    for (int y = 0; y <= p; ++y) {
                        const std::size_t j =
                            span - static_cast<std::size_t>(p) + static_cast<std::size_t>(y);
                        for (std::size_t r = 0; r < 3; ++r)
                            for (std::size_t q = 0; q < 3; ++q)
                                a[(3 * i + r) * dim + 3 * j + q] += nb[x] * nb[y] * mk[3 * r + q];
                    }
                    for (std::size_t r = 0; r < 3; ++r)
                        for (std::size_t q = 0; q < 3; ++q)
                            rhs[3 * i + r] += nb[x] * mk[3 * r + q] * rPts[k][q];
                }
            }
            const double lambda = rOptions.mSmoothing * static_cast<double>(m);
            for (std::size_t i = 1; i + 1 < n; ++i) {
                const std::size_t idx[3] = {i - 1, i, i + 1};
                const double w[3] = {1.0, -2.0, 1.0};
                for (int x = 0; x < 3; ++x)
                    for (int y = 0; y < 3; ++y)
                        for (std::size_t r = 0; r < 3; ++r)
                            a[(3 * idx[x] + r) * dim + 3 * idx[y] + r] += lambda * w[x] * w[y];
            }
            for (std::size_t q = 0; q < nc; ++q) {
                const std::size_t span = nurbs_find_span(p, c.mKnots, n, par[held[q]]);
                bc_basis(c, span, par[held[q]], nb, nullptr);
                for (int x = 0; x <= p; ++x) {
                    const std::size_t i =
                        span - static_cast<std::size_t>(p) + static_cast<std::size_t>(x);
                    for (std::size_t r = 0; r < 3; ++r) {
                        a[(3 * n + 3 * q + r) * dim + 3 * i + r] = nb[x];
                        a[(3 * i + r) * dim + 3 * n + 3 * q + r] = nb[x];
                    }
                }
                for (std::size_t r = 0; r < 3; ++r)
                    rhs[3 * n + 3 * q + r] = rPts[held[q]][r];
            }
            solved = lsq_solve(a, rhs, dim);
            if (!solved)
                break;
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t r = 0; r < 3; ++r)
                    c.mPoles[i][r] = rhs[3 * i + r];
            if (pass == rOptions.mParameterCorrections)
                break;
            // Parameter correction: each free point to its foot (Newton).
            std::vector<std::uint8_t> fixed(m, 0);
            for (std::size_t k : held)
                fixed[k] = 1;
            for (std::size_t k = 0; k < m; ++k) {
                if (fixed[k])
                    continue;
                double u = par[k];
                for (int it = 0; it < 8; ++it) {
                    PatchVec d1;
                    const PatchVec q = bspline_curve_eval(c, u, &d1);
                    const PatchVec r{q[0] - rPts[k][0], q[1] - rPts[k][1], q[2] - rPts[k][2]};
                    const double g = r[0] * d1[0] + r[1] * d1[1] + r[2] * d1[2];
                    const double h = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2];
                    if (!(h > 0.0))
                        break;
                    const double step = g / h;
                    u = std::max(0.0, std::min(1.0, u - step));
                    if (std::abs(step) < 1e-12)
                        break;
                }
                par[k] = u;
            }
            // Keep the parameters in order (a foot cannot pass a neighbour's).
            for (std::size_t k = 1; k < m; ++k)
                par[k] = std::max(par[k], par[k - 1]);
        }
        if (solved) {
            double worst = 0.0;
            for (std::size_t k = 0; k < m; ++k)
                worst = std::max(worst, distance(k, bspline_curve_eval(c, par[k])));
            out.mCurve = c;
            out.mParams = par;
            out.mMaxDist = worst;
            out.mOk = worst <= rOptions.mTolerance;
            if (out.mOk || n >= cap)
                return out;
        } else if (n >= cap) {
            return out;
        }
        n = std::min(cap, std::max(n + 1, 2 * n - static_cast<std::size_t>(p)));
    }
}

}  // namespace detail
}  // namespace meshioplusplus
