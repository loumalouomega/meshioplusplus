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
 * @file detail/iga_kernels.hpp
 * @brief The isogeometric kernels behind `operations/iga.hpp`: exact knot
 * insertion and degree elevation of a tensor-product NURBS net, and its
 * evaluation with first derivatives (The NURBS Book, algorithms A2.3, A5.4
 * and A5.9).
 *
 * A net is held *homogeneous*: every control point is a row of `C` doubles
 * `(w x, w f..., w)`, so the geometry, the weights and every NURBS field
 * refine with one kernel and one rounding. A non-rational net drops the
 * weight column and the weights stay exactly 1. Control nets list direction 0
 * fastest. The arithmetic is the order the Python twin
 * (`meshioplusplus/_iga.py`) follows, so the two give the same bits; the
 * kernels are scalar per column, which is what lets the twin vectorise over
 * lines and columns without changing a bit (`-ffp-contract=off`).
 *
 * Core-private (not installed).
 */

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes (private, not installed)
#include "nurbs_eval.hpp"

namespace meshioplusplus::detail {

/// A tensor-product control net in homogeneous form.
struct IgaNet {
    std::size_t mDim = 0;
    std::array<int, 3> mDegrees{0, 0, 0};
    std::array<std::vector<double>, 3> mKnots;
    std::array<std::size_t, 3> mCounts{1, 1, 1};
    /// Columns per control point.
    std::size_t mColumns = 0;
    /// `prod(mCounts)` rows of `mColumns` values, direction 0 fastest.
    std::vector<double> mValues;

    std::size_t NumPoints() const { return mCounts[0] * mCounts[1] * mCounts[2]; }
};

inline std::int64_t iga_binomial(int N, int K) {
    if (K < 0 || K > N)
        return 0;
    std::int64_t r = 1;
    for (int i = 1; i <= K; ++i)
        r = r * (N - K + i) / i;
    return r;
}

/// The knots of one direction that split every nonzero span in `S` equal
/// parts (the interior ones only), ascending.
inline std::vector<double> iga_split_knots(int P, const std::vector<double>& rKnots, int S) {
    std::vector<double> x;
    if (S <= 1)
        return x;
    const std::size_t n = rKnots.size() - static_cast<std::size_t>(P) - 1;
    for (const std::size_t i : nurbs_spans(P, rKnots, n)) {
        const double k0 = rKnots[i], k1 = rKnots[i + 1];
        for (int j = 1; j < S; ++j)
            x.push_back(k0 + (k1 - k0) * (static_cast<double>(j) / S));
    }
    return x;
}

/// A5.4: refine the knot vector of one curve of `C`-column control rows with
/// the ascending knots `rX` (at least one).
inline void iga_refine_line(int P, const std::vector<double>& rU, const double* pPw,
                            std::size_t Count, std::size_t C, const std::vector<double>& rX,
                            std::vector<double>& rUbar, std::vector<double>& rQw) {
    using I = std::ptrdiff_t;
    const I p = P, n = static_cast<I>(Count) - 1, m = n + p + 1;
    const I r = static_cast<I>(rX.size()) - 1;
    const I a = static_cast<I>(nurbs_find_span(P, rU, Count, rX[0]));
    const I b = static_cast<I>(nurbs_find_span(P, rU, Count, rX[static_cast<std::size_t>(r)])) + 1;
    rUbar.assign(static_cast<std::size_t>(m + r + 2), 0.0);
    rQw.assign((Count + static_cast<std::size_t>(r) + 1) * C, 0.0);
    const auto row = [C](std::vector<double>& rV, I Index) {
        return rV.data() + static_cast<std::size_t>(Index) * C;
    };
    for (I j = 0; j <= a - p; ++j)
        for (std::size_t c = 0; c < C; ++c)
            row(rQw, j)[c] = pPw[static_cast<std::size_t>(j) * C + c];
    for (I j = b - 1; j <= n; ++j)
        for (std::size_t c = 0; c < C; ++c)
            row(rQw, j + r + 1)[c] = pPw[static_cast<std::size_t>(j) * C + c];
    for (I j = 0; j <= a; ++j)
        rUbar[static_cast<std::size_t>(j)] = rU[static_cast<std::size_t>(j)];
    for (I j = b + p; j <= m; ++j)
        rUbar[static_cast<std::size_t>(j + r + 1)] = rU[static_cast<std::size_t>(j)];
    I i = b + p - 1, k = b + p + r;
    for (I j = r; j >= 0; --j) {
        const double xj = rX[static_cast<std::size_t>(j)];
        while (xj <= rU[static_cast<std::size_t>(i)] && i > a) {
            for (std::size_t c = 0; c < C; ++c)
                row(rQw, k - p - 1)[c] = pPw[static_cast<std::size_t>(i - p - 1) * C + c];
            rUbar[static_cast<std::size_t>(k)] = rU[static_cast<std::size_t>(i)];
            --k;
            --i;
        }
        for (std::size_t c = 0; c < C; ++c)
            row(rQw, k - p - 1)[c] = row(rQw, k - p)[c];
        for (I l = 1; l <= p; ++l) {
            const I ind = k - p + l;
            double alfa = rUbar[static_cast<std::size_t>(k + l)] - xj;
            if (std::fabs(alfa) == 0.0) {
                for (std::size_t c = 0; c < C; ++c)
                    row(rQw, ind - 1)[c] = row(rQw, ind)[c];
            } else {
                alfa = alfa / (rUbar[static_cast<std::size_t>(k + l)] -
                               rU[static_cast<std::size_t>(i - p + l)]);
                for (std::size_t c = 0; c < C; ++c)
                    row(rQw, ind - 1)[c] =
                        alfa * row(rQw, ind - 1)[c] + (1.0 - alfa) * row(rQw, ind)[c];
            }
        }
        rUbar[static_cast<std::size_t>(k)] = xj;
        --k;
    }
}

/// A5.9: raise the degree of one curve of `C`-column control rows by `T`
/// (at least one), keeping every knot's continuity.
inline void iga_elevate_line(int P, const std::vector<double>& rU, const double* pPw,
                             std::size_t Count, std::size_t C, int T, std::vector<double>& rUh,
                             std::vector<double>& rQw) {
    using I = std::ptrdiff_t;
    const I p = P, t = T, ph = p + t, ph2 = ph / 2;
    const I n = static_cast<I>(Count) - 1, m = n + p + 1;
    // The distinct interior knots bound the segments: the output holds at most
    // (Count + t * (segments + 1)) rows and that many + ph + 1 knots.
    I distinct = 0;
    for (I i = p + 1; i <= n; ++i)
        if (rU[static_cast<std::size_t>(i)] != rU[static_cast<std::size_t>(i - 1)])
            ++distinct;
    const I max_rows = n + 1 + t * (distinct + 1);
    rQw.assign(static_cast<std::size_t>(max_rows) * C, 0.0);
    rUh.assign(static_cast<std::size_t>(max_rows + ph + 1), 0.0);
    std::vector<std::vector<double>> bezalfs(
        static_cast<std::size_t>(ph + 1),
        std::vector<double>(static_cast<std::size_t>(p + 1), 0.0));
    bezalfs[0][0] = 1.0;
    bezalfs[static_cast<std::size_t>(ph)][static_cast<std::size_t>(p)] = 1.0;
    for (I i = 1; i <= ph2; ++i) {
        const double inv =
            1.0 / static_cast<double>(iga_binomial(static_cast<int>(ph), static_cast<int>(i)));
        const I mpi = std::min(p, i);
        for (I j = std::max<I>(0, i - t); j <= mpi; ++j)
            bezalfs[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                inv * static_cast<double>(iga_binomial(static_cast<int>(p), static_cast<int>(j))) *
                static_cast<double>(iga_binomial(static_cast<int>(t), static_cast<int>(i - j)));
    }
    for (I i = ph2 + 1; i <= ph - 1; ++i) {
        const I mpi = std::min(p, i);
        for (I j = std::max<I>(0, i - t); j <= mpi; ++j)
            bezalfs[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                bezalfs[static_cast<std::size_t>(ph - i)][static_cast<std::size_t>(p - j)];
    }
    std::vector<double> bpts(static_cast<std::size_t>(p + 1) * C),
        ebpts(static_cast<std::size_t>(ph + 1) * C),
        nextbpts(static_cast<std::size_t>(std::max<I>(p - 1, 1)) * C),
        alfs(static_cast<std::size_t>(std::max<I>(p - 1, 1)));
    const auto brow = [C](std::vector<double>& rV, I Index) {
        return rV.data() + static_cast<std::size_t>(Index) * C;
    };
    const auto U = [&rU](I Index) { return rU[static_cast<std::size_t>(Index)]; };
    const auto Uh = [&rUh](I Index) -> double& { return rUh[static_cast<std::size_t>(Index)]; };
    I mh = ph, kind = ph + 1, r = -1, a = p, b = p + 1, cind = 1;
    double ua = U(0);
    for (std::size_t c = 0; c < C; ++c)
        brow(rQw, 0)[c] = pPw[c];
    for (I i = 0; i <= ph; ++i)
        Uh(i) = ua;
    for (I i = 0; i <= p; ++i)
        for (std::size_t c = 0; c < C; ++c)
            brow(bpts, i)[c] = pPw[static_cast<std::size_t>(i) * C + c];
    while (b < m) {
        I i = b;
        while (b < m && U(b) == U(b + 1))
            ++b;
        const I mul = b - i + 1;
        mh = mh + mul + t;
        const double ub = U(b);
        const I oldr = r;
        r = p - mul;
        const I lbz = oldr > 0 ? (oldr + 2) / 2 : 1;
        const I rbz = r > 0 ? ph - (r + 1) / 2 : ph;
        if (r > 0) {
            const double numer = ub - ua;
            for (I k = p; k > mul; --k)
                alfs[static_cast<std::size_t>(k - mul - 1)] = numer / (U(a + k) - ua);
            for (I j = 1; j <= r; ++j) {
                const I save = r - j, s = mul + j;
                for (I k = p; k >= s; --k) {
                    const double al = alfs[static_cast<std::size_t>(k - s)];
                    for (std::size_t c = 0; c < C; ++c)
                        brow(bpts, k)[c] =
                            al * brow(bpts, k)[c] + (1.0 - al) * brow(bpts, k - 1)[c];
                }
                for (std::size_t c = 0; c < C; ++c)
                    brow(nextbpts, save)[c] = brow(bpts, p)[c];
            }
        }
        for (I k = lbz; k <= ph; ++k) {
            for (std::size_t c = 0; c < C; ++c)
                brow(ebpts, k)[c] = 0.0;
            const I mpi = std::min(p, k);
            for (I j = std::max<I>(0, k - t); j <= mpi; ++j) {
                const double bz = bezalfs[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
                for (std::size_t c = 0; c < C; ++c)
                    brow(ebpts, k)[c] = brow(ebpts, k)[c] + bz * brow(bpts, j)[c];
            }
        }
        if (oldr > 1) {
            I first = kind - 2, last = kind;
            const double den = ub - ua;
            const double bet = (ub - Uh(kind - 1)) / den;
            for (I tr = 1; tr < oldr; ++tr) {
                I ii = first, j = last, kj = j - kind + 1;
                while (j - ii > tr) {
                    if (ii < cind) {
                        const double alf = (ub - Uh(ii)) / (ua - Uh(ii));
                        for (std::size_t c = 0; c < C; ++c)
                            brow(rQw, ii)[c] =
                                alf * brow(rQw, ii)[c] + (1.0 - alf) * brow(rQw, ii - 1)[c];
                    }
                    if (j >= lbz) {
                        if (j - tr <= kind - ph + oldr) {
                            const double gam = (ub - Uh(j - tr)) / den;
                            for (std::size_t c = 0; c < C; ++c)
                                brow(ebpts, kj)[c] =
                                    gam * brow(ebpts, kj)[c] + (1.0 - gam) * brow(ebpts, kj + 1)[c];
                        } else {
                            for (std::size_t c = 0; c < C; ++c)
                                brow(ebpts, kj)[c] =
                                    bet * brow(ebpts, kj)[c] + (1.0 - bet) * brow(ebpts, kj + 1)[c];
                        }
                    }
                    ++ii;
                    --j;
                    --kj;
                }
                --first;
                ++last;
            }
        }
        if (a != p)
            for (I k = 0; k < ph - oldr; ++k) {
                Uh(kind) = ua;
                ++kind;
            }
        for (I j = lbz; j <= rbz; ++j) {
            for (std::size_t c = 0; c < C; ++c)
                brow(rQw, cind)[c] = brow(ebpts, j)[c];
            ++cind;
        }
        if (b < m) {
            for (I j = 0; j < r; ++j)
                for (std::size_t c = 0; c < C; ++c)
                    brow(bpts, j)[c] = brow(nextbpts, j)[c];
            for (I j = r; j <= p; ++j)
                for (std::size_t c = 0; c < C; ++c)
                    brow(bpts, j)[c] = pPw[static_cast<std::size_t>(b - p + j) * C + c];
            a = b;
            ++b;
            ua = ub;
        } else {
            for (I k = 0; k <= ph; ++k)
                Uh(kind + k) = ub;
        }
    }
    const I nh = mh - ph - 1;
    rQw.resize(static_cast<std::size_t>(nh + 1) * C);
    rUh.resize(static_cast<std::size_t>(nh + ph + 2));
}

/// Apply a line operation `Op(line values, count, out knots, out values)` along
/// direction `Dir` of every line of the net.
template <class Op>
inline void iga_apply_direction(IgaNet& rNet, std::size_t Dir, std::size_t NewCount, Op&& rOp) {
    const std::size_t C = rNet.mColumns;
    std::size_t inner = 1, outer = 1;
    for (std::size_t e = 0; e < Dir; ++e)
        inner *= rNet.mCounts[e];
    for (std::size_t e = Dir + 1; e < 3; ++e)
        outer *= rNet.mCounts[e];
    const std::size_t n = rNet.mCounts[Dir];
    std::vector<double> out(inner * NewCount * outer * C);
    std::vector<double> line(n * C);
    for (std::size_t o = 0; o < outer; ++o)
        for (std::size_t in = 0; in < inner; ++in) {
            for (std::size_t k = 0; k < n; ++k)
                for (std::size_t c = 0; c < C; ++c)
                    line[k * C + c] = rNet.mValues[(in + inner * (k + n * o)) * C + c];
            const std::vector<double>& res = rOp(line);
            for (std::size_t k = 0; k < NewCount; ++k)
                for (std::size_t c = 0; c < C; ++c)
                    out[(in + inner * (k + NewCount * o)) * C + c] = res[k * C + c];
        }
    rNet.mValues = std::move(out);
    rNet.mCounts[Dir] = NewCount;
}

/// Insert `rX` (ascending) into direction `Dir`'s knot vector.
inline void iga_net_insert(IgaNet& rNet, std::size_t Dir, const std::vector<double>& rX) {
    if (rX.empty())
        return;
    const int p = rNet.mDegrees[Dir];
    const std::size_t n = rNet.mCounts[Dir];
    const double lo = rNet.mKnots[Dir][static_cast<std::size_t>(p)], hi = rNet.mKnots[Dir][n];
    for (std::size_t i = 0; i < rX.size(); ++i)
        if (!(rX[i] > lo && rX[i] < hi) || (i && rX[i] < rX[i - 1]))
            throw std::invalid_argument(
                "meshio++: iga: knots to insert must be ascending and inside the parameter range");
    std::vector<double> ubar, res;
    iga_apply_direction(rNet, Dir, n + rX.size(),
                        [&](const std::vector<double>& rLine) -> const std::vector<double>& {
                            iga_refine_line(p, rNet.mKnots[Dir], rLine.data(), n, rNet.mColumns, rX,
                                            ubar, res);
                            return res;
                        });
    rNet.mKnots[Dir] = ubar;
}

/// Raise direction `Dir`'s degree by `T`.
inline void iga_net_elevate(IgaNet& rNet, std::size_t Dir, int T) {
    if (T <= 0)
        return;
    const int p = rNet.mDegrees[Dir];
    const std::size_t n = rNet.mCounts[Dir];
    std::vector<double> uh, res;
    std::vector<double> uh_first;
    std::size_t new_count = 0;
    // The new count and knot vector depend on the knots only: take them from
    // a first, zero line.
    {
        std::vector<double> zero(n * rNet.mColumns, 0.0);
        iga_elevate_line(p, rNet.mKnots[Dir], zero.data(), n, rNet.mColumns, T, uh_first, res);
        new_count = res.size() / rNet.mColumns;
    }
    iga_apply_direction(
        rNet, Dir, new_count, [&](const std::vector<double>& rLine) -> const std::vector<double>& {
            iga_elevate_line(p, rNet.mKnots[Dir], rLine.data(), n, rNet.mColumns, T, uh, res);
            return res;
        });
    rNet.mKnots[Dir] = uh_first;
    rNet.mDegrees[Dir] = p + T;
}

/// Scratch doubles for the per-point evaluators: on the stack up to `N`,
/// past it on the heap. An allocation per call was most of their cost.
template <std::size_t N>
class IgaScratch {
public:
    IgaScratch() = default;
    explicit IgaScratch(std::size_t Size, double Value = 0.0) { Assign(Size, Value); }
    IgaScratch(const IgaScratch&) = delete;
    IgaScratch& operator=(const IgaScratch&) = delete;

    void Assign(std::size_t Size, double Value) {
        if (Size > N) {
            mHeap.assign(Size, Value);
            mData = mHeap.data();
        } else {
            std::fill(mFixed.begin(), mFixed.begin() + static_cast<std::ptrdiff_t>(Size), Value);
            mData = mFixed.data();
        }
    }
    double* data() { return mData; }
    double& operator[](std::size_t I) { return mData[I]; }
    double operator[](std::size_t I) const { return mData[I]; }

private:
    std::array<double, N> mFixed;
    std::vector<double> mHeap;
    double* mData = mFixed.data();
};

/// A2.3 for first derivatives: `pDers` is `(NDer + 1)` rows of `P + 1`
/// values, row 0 the basis functions, row 1 their derivatives.
inline void iga_basis_ders(std::size_t Span, double U, int P, int NDer,
                           const std::vector<double>& rKnots, double* pDers) {
    const auto p = static_cast<std::size_t>(P);
    IgaScratch<64> ndu((p + 1) * (p + 1)), left(p + 1), right(p + 1);
    const auto at = [p](std::size_t R, std::size_t Cc) { return R * (p + 1) + Cc; };
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
    for (std::size_t j = 0; j <= p; ++j)
        pDers[j] = ndu[at(j, p)];
    if (NDer < 1)
        return;
    // The first derivatives (A2.3 with k = 1: its inner loop is empty, and the
    // work table's entry (0, 0) is exactly 1).
    const double factor = static_cast<double>(P);
    for (std::size_t r = 0; r <= p; ++r) {
        double d = 0.0;
        if (r >= 1)
            d = (1.0 / ndu[at(p, r - 1)]) * ndu[at(r - 1, p - 1)];
        if (r + 1 <= p)
            d += (-1.0 / ndu[at(p, r)]) * ndu[at(r, p - 1)];
        pDers[(p + 1) + r] = d;
    }
    for (std::size_t j = 0; j <= p; ++j)
        pDers[(p + 1) + j] *= factor;
}

/// A tensor-product net to evaluate: `NV` columns per control point, weights
/// optional (null: polynomial).
struct IgaEvalNet {
    std::size_t mDim = 0;
    std::array<int, 3> mDegrees{0, 0, 0};
    std::array<const std::vector<double>*, 3> mKnots{nullptr, nullptr, nullptr};
    std::array<std::size_t, 3> mCounts{1, 1, 1};
    const double* mValues = nullptr;
    std::size_t mColumns = 0;
    const double* mWeights = nullptr;
};

/// The `NV` values at the parameters `pU[0 .. dim)` in knot spans `pSpan`, and
/// with `NDer == 1` their derivatives along each direction (`pDer`, row-major
/// `(dim, NV)`).
inline void iga_eval_point(const IgaEvalNet& rNet, const std::size_t* pSpan, const double* pU,
                           int NDer, double* pOut, double* pDer) {
    const std::size_t dim = rNet.mDim, NV = rNet.mColumns;
    std::array<IgaScratch<32>, 3> ders;
    for (std::size_t d = 0; d < dim; ++d) {
        ders[d].Assign(
            static_cast<std::size_t>(NDer + 1) * (static_cast<std::size_t>(rNet.mDegrees[d]) + 1),
            0.0);
        iga_basis_ders(pSpan[d], pU[d], rNet.mDegrees[d], NDer, *rNet.mKnots[d], ders[d].data());
    }
    IgaScratch<16> acc(NV);
    IgaScratch<48> dacc(dim * NV);
    double wsum = 0.0;
    std::array<double, 3> dwsum{0.0, 0.0, 0.0};
    std::array<std::size_t, 3> np1{1, 1, 1}, base{0, 0, 0};
    for (std::size_t d = 0; d < dim; ++d) {
        np1[d] = static_cast<std::size_t>(rNet.mDegrees[d]) + 1;
        base[d] = pSpan[d] - static_cast<std::size_t>(rNet.mDegrees[d]);
    }
    for (std::size_t j2 = 0; j2 < np1[2]; ++j2)
        for (std::size_t j1 = 0; j1 < np1[1]; ++j1)
            for (std::size_t j0 = 0; j0 < np1[0]; ++j0) {
                const std::array<std::size_t, 3> j{j0, j1, j2};
                const std::size_t idx =
                    (base[0] + j0) +
                    rNet.mCounts[0] * ((base[1] + j1) + rNet.mCounts[1] * (base[2] + j2));
                const double w = rNet.mWeights ? rNet.mWeights[idx] : 1.0;
                const double* v = rNet.mValues + idx * NV;
                double b = 1.0;
                for (std::size_t d = 0; d < dim; ++d)
                    b *= ders[d][j[d]];
                const double cc = w * b;
                for (std::size_t c = 0; c < NV; ++c)
                    acc[c] += cc * v[c];
                wsum += cc;
                if (NDer >= 1)
                    for (std::size_t k = 0; k < dim; ++k) {
                        double bk = 1.0;
                        for (std::size_t d = 0; d < dim; ++d)
                            bk *= d == k ? ders[d][np1[d] + j[d]] : ders[d][j[d]];
                        const double ck = w * bk;
                        for (std::size_t c = 0; c < NV; ++c)
                            dacc[k * NV + c] += ck * v[c];
                        dwsum[k] += ck;
                    }
            }
    if (rNet.mWeights) {
        for (std::size_t c = 0; c < NV; ++c)
            pOut[c] = acc[c] / wsum;
        if (NDer >= 1)
            for (std::size_t k = 0; k < dim; ++k)
                for (std::size_t c = 0; c < NV; ++c)
                    pDer[k * NV + c] = (dacc[k * NV + c] - dwsum[k] * pOut[c]) / wsum;
    } else {
        for (std::size_t c = 0; c < NV; ++c)
            pOut[c] = acc[c];
        if (NDer >= 1)
            for (std::size_t k = 0; k < dim; ++k)
                for (std::size_t c = 0; c < NV; ++c)
                    pDer[k * NV + c] = dacc[k * NV + c];
    }
}

/// A2.3 to the second derivatives: `pDers` is three rows of `P + 1` values,
/// the basis functions, their first and their second derivatives (zero above
/// the degree). Kept apart from `iga_basis_ders`, whose path the Python twin
/// matches bit for bit.
inline void iga_basis_ders2(std::size_t Span, double U, int P, const std::vector<double>& rKnots,
                            double* pDers) {
    const auto p = static_cast<std::size_t>(P);
    const std::size_t w = p + 1;
    IgaScratch<64> ndu(w * w), left(w), right(w), a(2 * w);
    const auto at = [w](std::size_t R, std::size_t Cc) { return R * w + Cc; };
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
    for (std::size_t j = 0; j <= p; ++j) {
        pDers[j] = ndu[at(j, p)];
        pDers[w + j] = 0.0;
        pDers[2 * w + j] = 0.0;
    }
    const int n = std::min(2, P);
    for (int r = 0; r <= P; ++r) {
        std::size_t s1 = 0, s2 = 1;
        a[0] = 1.0;
        for (int k = 1; k <= n; ++k) {
            double d = 0.0;
            const int rk = r - k, pk = P - k;
            const auto row = [&](std::size_t S, int J) -> double& {
                return a[S * w + static_cast<std::size_t>(J)];
            };
            if (r >= k) {
                row(s2, 0) =
                    row(s1, 0) /
                    ndu[at(static_cast<std::size_t>(pk + 1), static_cast<std::size_t>(rk))];
                d = row(s2, 0) *
                    ndu[at(static_cast<std::size_t>(rk), static_cast<std::size_t>(pk))];
            }
            const int j1 = rk >= -1 ? 1 : -rk;
            const int j2 = r - 1 <= pk ? k - 1 : P - r;
            for (int j = j1; j <= j2; ++j) {
                row(s2, j) =
                    (row(s1, j) - row(s1, j - 1)) /
                    ndu[at(static_cast<std::size_t>(pk + 1), static_cast<std::size_t>(rk + j))];
                d += row(s2, j) *
                     ndu[at(static_cast<std::size_t>(rk + j), static_cast<std::size_t>(pk))];
            }
            if (r <= pk) {
                row(s2, k) = -row(s1, k - 1) /
                             ndu[at(static_cast<std::size_t>(pk + 1), static_cast<std::size_t>(r))];
                d +=
                    row(s2, k) * ndu[at(static_cast<std::size_t>(r), static_cast<std::size_t>(pk))];
            }
            pDers[static_cast<std::size_t>(k) * w + static_cast<std::size_t>(r)] = d;
            std::swap(s1, s2);
        }
    }
    double factor = static_cast<double>(P);
    for (int k = 1; k <= n; ++k) {
        for (std::size_t j = 0; j <= p; ++j)
            pDers[static_cast<std::size_t>(k) * w + j] *= factor;
        factor *= static_cast<double>(P - k);
    }
}

/// A surface net (`mDim == 2`) at the parameters `pU` in spans `pSpan`: the
/// `NV` values, their first derivatives (`pD1`, rows u and v) and second
/// derivatives (`pD2`, rows uu, uv and vv), the rational ones by the quotient
/// rule of A4.4 (Piegl and Tiller 1997).
inline void iga_eval_surface_d2(const IgaEvalNet& rNet, const std::size_t* pSpan, const double* pU,
                                double* pOut, double* pD1, double* pD2) {
    const std::size_t NV = rNet.mColumns;
    std::array<IgaScratch<48>, 2> ders;
    std::array<std::size_t, 2> np1{}, base{};
    for (std::size_t d = 0; d < 2; ++d) {
        np1[d] = static_cast<std::size_t>(rNet.mDegrees[d]) + 1;
        base[d] = pSpan[d] - static_cast<std::size_t>(rNet.mDegrees[d]);
        ders[d].Assign(3 * np1[d], 0.0);
        iga_basis_ders2(pSpan[d], pU[d], rNet.mDegrees[d], *rNet.mKnots[d], ders[d].data());
    }
    // Homogeneous sums: index 0 value, 1 u, 2 v, 3 uu, 4 uv, 5 vv.
    constexpr std::array<std::array<int, 2>, 6> kOrder{
        {{0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2}}};
    IgaScratch<24> acc(6 * NV);
    std::array<double, 6> wsum{};
    for (std::size_t j1 = 0; j1 < np1[1]; ++j1)
        for (std::size_t j0 = 0; j0 < np1[0]; ++j0) {
            const std::size_t idx = (base[0] + j0) + rNet.mCounts[0] * (base[1] + j1);
            const double wt = rNet.mWeights ? rNet.mWeights[idx] : 1.0;
            const double* v = rNet.mValues + idx * NV;
            for (std::size_t q = 0; q < 6; ++q) {
                const double b = ders[0][static_cast<std::size_t>(kOrder[q][0]) * np1[0] + j0] *
                                 ders[1][static_cast<std::size_t>(kOrder[q][1]) * np1[1] + j1];
                const double cc = wt * b;
                for (std::size_t c = 0; c < NV; ++c)
                    acc[q * NV + c] += cc * v[c];
                wsum[q] += cc;
            }
        }
    const double w0 = rNet.mWeights ? wsum[0] : 1.0;
    for (std::size_t c = 0; c < NV; ++c) {
        if (!rNet.mWeights) {
            pOut[c] = acc[c];
            pD1[c] = acc[NV + c];
            pD1[NV + c] = acc[2 * NV + c];
            pD2[c] = acc[3 * NV + c];
            pD2[NV + c] = acc[4 * NV + c];
            pD2[2 * NV + c] = acc[5 * NV + c];
            continue;
        }
        const double x = acc[c] / w0;
        const double xu = (acc[NV + c] - wsum[1] * x) / w0;
        const double xv = (acc[2 * NV + c] - wsum[2] * x) / w0;
        pOut[c] = x;
        pD1[c] = xu;
        pD1[NV + c] = xv;
        pD2[c] = (acc[3 * NV + c] - 2.0 * wsum[1] * xu - wsum[3] * x) / w0;
        pD2[NV + c] = (acc[4 * NV + c] - wsum[1] * xv - wsum[2] * xu - wsum[4] * x) / w0;
        pD2[2 * NV + c] = (acc[5 * NV + c] - 2.0 * wsum[2] * xv - wsum[5] * x) / w0;
    }
}

}  // namespace meshioplusplus::detail
