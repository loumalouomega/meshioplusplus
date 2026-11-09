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
// A small dense two-phase simplex. See detail/dense_lp.hpp.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "dense_lp.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kLpEps = 1e-12;

// The tableau: rows 0..m-1 the constraints, row m the reduced costs; the last
// column the right-hand side (the objective row holds minus the objective).
struct LpTableau {
    std::size_t mRows = 0;
    std::size_t mCols = 0;  // without the right-hand side
    std::vector<double> mT;
    std::vector<std::size_t> mBasis;

    double& At(std::size_t R, std::size_t C) { return mT[R * (mCols + 1) + C]; }
    double Rhs(std::size_t R) const { return mT[R * (mCols + 1) + mCols]; }

    void Pivot(std::size_t P, std::size_t Q) {
        const std::size_t w = mCols + 1;
        const double piv = mT[P * w + Q];
        for (std::size_t c = 0; c < w; ++c)
            mT[P * w + c] /= piv;
        for (std::size_t r = 0; r <= mRows; ++r) {
            if (r == P)
                continue;
            const double f = mT[r * w + Q];
            if (f == 0.0)
                continue;
            for (std::size_t c = 0; c < w; ++c)
                mT[r * w + c] -= f * mT[P * w + c];
            mT[r * w + Q] = 0.0;
        }
        mBasis[P] = Q;
    }

    // Bland's rule: the lowest allowed column with a negative reduced cost
    // enters; the leaving row by the minimum ratio, the lowest basic index on
    // a tie. Returns false at the optimum; sets rUnbounded when a column
    // enters that nothing bounds.
    bool Step(std::size_t Allowed, double CostEps, double PivEps, bool& rUnbounded) {
        std::size_t q = Allowed;
        for (std::size_t c = 0; c < Allowed; ++c)
            if (At(mRows, c) < -CostEps) {
                q = c;
                break;
            }
        if (q == Allowed)
            return false;
        std::size_t p = mRows;
        double best = 0.0;
        for (std::size_t r = 0; r < mRows; ++r) {
            const double a = At(r, q);
            if (!(a > PivEps))
                continue;
            const double ratio = Rhs(r) / a;
            if (p == mRows || ratio < best || (ratio == best && mBasis[r] < mBasis[p])) {
                p = r;
                best = ratio;
            }
        }
        if (p == mRows) {
            rUnbounded = true;
            return false;
        }
        Pivot(p, q);
        return true;
    }
};

}  // namespace

DenseLpResult dense_lp_equality(const std::vector<double>& rC, const std::vector<double>& rA,
                                const std::vector<double>& rB, std::int64_t MaxPivots) {
    DenseLpResult out;
    const std::size_t n = rC.size(), m = rB.size();
    if (rA.size() != n * m)
        return out;
    double scale_a = 0.0, scale_c = 1.0, scale_b = 1.0;
    for (const double v : rA)
        scale_a = std::max(scale_a, std::fabs(v));
    for (const double v : rC)
        scale_c = std::max(scale_c, std::fabs(v));
    for (const double v : rB)
        scale_b = std::max(scale_b, std::fabs(v));
    const double piv_eps = kLpEps * std::max(scale_a, 1e-300);

    LpTableau t;
    t.mRows = m;
    t.mCols = n + m;
    t.mT.assign((m + 1) * (n + m + 1), 0.0);
    t.mBasis.resize(m);
    std::vector<double> sign(m, 1.0);
    for (std::size_t r = 0; r < m; ++r) {
        sign[r] = rB[r] < 0.0 ? -1.0 : 1.0;
        for (std::size_t c = 0; c < n; ++c)
            t.At(r, c) = sign[r] * rA[r * n + c];
        t.At(r, n + r) = 1.0;
        t.At(r, n + m) = sign[r] * rB[r];
        t.mBasis[r] = n + r;
    }
    // Phase one: minimise the artificials' sum.
    for (std::size_t c = 0; c <= n + m; ++c) {
        if (c >= n && c < n + m)
            continue;
        double s = 0.0;
        for (std::size_t r = 0; r < m; ++r)
            s += t.At(r, c);
        t.At(m, c) = -s;
    }
    bool unbounded = false;
    while (t.Step(n, kLpEps, piv_eps, unbounded)) {
        if (++out.mPivots > MaxPivots) {
            out.mStatus = DenseLpStatus::IterationLimit;
            return out;
        }
    }
    if (-t.Rhs(m) > kLpEps * scale_b * static_cast<double>(std::max<std::size_t>(m, 1))) {
        out.mStatus = DenseLpStatus::Infeasible;
        return out;
    }
    // Drive the artificials still basic (at zero) out where a column allows.
    for (std::size_t r = 0; r < m; ++r) {
        if (t.mBasis[r] < n)
            continue;
        for (std::size_t c = 0; c < n; ++c)
            if (std::fabs(t.At(r, c)) > piv_eps) {
                t.Pivot(r, c);
                break;
            }
    }
    // Phase two: the objective, the artificials kept out.
    for (std::size_t c = 0; c <= n + m; ++c) {
        double z = c < n ? rC[c] : 0.0;
        for (std::size_t r = 0; r < m; ++r) {
            const std::size_t b = t.mBasis[r];
            if (b < n)
                z -= rC[b] * t.At(r, c);
        }
        t.At(m, c) = z;
    }
    unbounded = false;
    while (t.Step(n, kLpEps * scale_c, piv_eps, unbounded)) {
        if (++out.mPivots > MaxPivots) {
            out.mStatus = DenseLpStatus::IterationLimit;
            return out;
        }
    }
    if (unbounded) {
        out.mStatus = DenseLpStatus::Unbounded;
        return out;
    }
    out.mStatus = DenseLpStatus::Optimal;
    out.mX.assign(n, 0.0);
    for (std::size_t r = 0; r < m; ++r)
        if (t.mBasis[r] < n)
            out.mX[t.mBasis[r]] = t.Rhs(r);
    out.mObjective = 0.0;
    for (std::size_t c = 0; c < n; ++c)
        out.mObjective += rC[c] * out.mX[c];
    // The duals: an artificial column's reduced cost is minus its row's dual.
    out.mDuals.assign(m, 0.0);
    for (std::size_t r = 0; r < m; ++r)
        out.mDuals[r] = -sign[r] * t.At(m, n + r);
    return out;
}

DenseLpResult dense_lp_inequality(const std::vector<double>& rC, const std::vector<double>& rA,
                                  const std::vector<double>& rB, std::int64_t MaxPivots) {
    const std::size_t k = rC.size(), m = rB.size();
    DenseLpResult out;
    if (rA.size() != k * m)
        return out;
    // The dual: min b.y subject to A^T y = -c, y >= 0.
    std::vector<double> at(k * m);
    for (std::size_t r = 0; r < m; ++r)
        for (std::size_t c = 0; c < k; ++c)
            at[c * m + r] = rA[r * k + c];
    std::vector<double> rhs(k);
    for (std::size_t c = 0; c < k; ++c)
        rhs[c] = -rC[c];
    const DenseLpResult dual = dense_lp_equality(rB, at, rhs, MaxPivots);
    out.mPivots = dual.mPivots;
    if (dual.mStatus == DenseLpStatus::Infeasible) {
        out.mStatus = DenseLpStatus::Unbounded;
        return out;
    }
    if (dual.mStatus == DenseLpStatus::Unbounded) {
        out.mStatus = DenseLpStatus::Infeasible;
        return out;
    }
    if (dual.mStatus != DenseLpStatus::Optimal) {
        out.mStatus = dual.mStatus;
        return out;
    }
    out.mStatus = DenseLpStatus::Optimal;
    out.mX = dual.mDuals;
    out.mObjective = 0.0;
    for (std::size_t c = 0; c < k; ++c)
        out.mObjective += rC[c] * out.mX[c];
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
