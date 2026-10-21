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

// Truncated hierarchical B-splines. See detail/thb_spline.hpp.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/log.hpp"

// Project includes (private, not installed)
#include "dense_chol.hpp"
#include "thb_spline.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The fewest samples a cell needs before it is refined: a spline finer than
// its data is free between the samples and oscillates.
constexpr std::size_t kThbMinCellSamples = 8;

// The 1-D refinement matrix from level `l` (`S` spans) to level `l + 1`
// (`2 S`), clamped uniform knots on [0, 1]: column `i` holds the level-`l + 1`
// coefficients of the level-`l` B-spline `i` (Boehm's knot insertion, The
// NURBS Book A5.1, of every new midpoint), row-major `(2 S + p) x (S + p)`.
std::vector<double> thb_refine_1d(int P, std::size_t S) {
    const std::size_t p = static_cast<std::size_t>(P);
    const std::size_t n = S + p;
    std::vector<double> knots = bspline_clamped_knots(P, n, 0.0, 1.0);
    // The control "points": the identity, one column per coarse function.
    std::vector<std::vector<double>> pts(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i)
        pts[i][i] = 1.0;
    for (std::size_t k = 0; k < S; ++k) {
        const double t = (static_cast<double>(k) + 0.5) / static_cast<double>(S);
        // The span of t: knots[span] <= t < knots[span + 1].
        std::size_t span = p;
        while (span + 1 < knots.size() && knots[span + 1] <= t)
            ++span;
        std::vector<std::vector<double>> q(pts.size() + 1);
        for (std::size_t i = 0; i <= pts.size(); ++i) {
            if (i + p <= span) {
                q[i] = pts[i];
            } else if (i > span) {
                q[i] = pts[i - 1];
            } else {
                const double a = (t - knots[i]) / (knots[i + p] - knots[i]);
                q[i].resize(n);
                for (std::size_t c = 0; c < n; ++c)
                    q[i][c] = a * pts[i][c] + (1.0 - a) * pts[i - 1][c];
            }
        }
        knots.insert(knots.begin() + static_cast<std::ptrdiff_t>(span) + 1, t);
        pts = std::move(q);
    }
    std::vector<double> a(pts.size() * n);
    for (std::size_t r = 0; r < pts.size(); ++r)
        for (std::size_t c = 0; c < n; ++c)
            a[r * n + c] = std::fabs(pts[r][c]) < 1e-15 ? 0.0 : pts[r][c];
    return a;
}

// Whether every cell of the support of level-`L` B-spline `(I, J)` passes `In`.
template <class F>
bool thb_support_all(int P, std::size_t S, std::size_t I, std::size_t J, F&& rIn) {
    const std::size_t p = static_cast<std::size_t>(P);
    const std::size_t i0 = I >= p ? I - p : 0, i1 = std::min(I, S - 1);
    const std::size_t j0 = J >= p ? J - p : 0, j1 = std::min(J, S - 1);
    for (std::size_t j = j0; j <= j1; ++j)
        for (std::size_t i = i0; i <= i1; ++i)
            if (!rIn(i, j))
                return false;
    return true;
}

}  // namespace

ThbBasis thb_hierarchy(int Degree, std::size_t Spans0, int Levels) {
    ThbBasis b;
    b.mDegree = Degree;
    b.mSpans0 = std::max<std::size_t>(Spans0, 1);
    b.mRefined.resize(static_cast<std::size_t>(std::max(Levels, 1)));
    for (int l = 0; l < b.NumLevels(); ++l)
        b.mRefined[static_cast<std::size_t>(l)].assign(b.Spans(l) * b.Spans(l), 0);
    return b;
}

void thb_build(ThbBasis& rB) {
    rB.mFunctions.clear();
    const int p = rB.mDegree;
    const int levels = rB.NumLevels();
    // The finest level the hierarchy uses: past it nothing is refined.
    int top = 0;
    for (int l = 0; l + 1 < levels; ++l)
        if (std::any_of(rB.mRefined[static_cast<std::size_t>(l)].begin(),
                        rB.mRefined[static_cast<std::size_t>(l)].end(),
                        [](std::uint8_t X) { return X != 0; }))
            top = l + 1;
    // Cells of level l inside Omega^l: level 0 all, else a refined parent.
    const auto in_domain = [&](int L, std::size_t I, std::size_t J) {
        if (L == 0)
            return true;
        const std::size_t ps = rB.Spans(L - 1);
        return rB.mRefined[static_cast<std::size_t>(L - 1)][(I / 2) + ps * (J / 2)] != 0;
    };
    const auto refined = [&](int L, std::size_t I, std::size_t J) {
        return L < top && rB.mRefined[static_cast<std::size_t>(L)][I + rB.Spans(L) * J] != 0;
    };
    std::vector<std::vector<double>> refine(static_cast<std::size_t>(top));
    for (int l = 0; l < top; ++l)
        refine[static_cast<std::size_t>(l)] = thb_refine_1d(p, rB.Spans(l));
    // Normalise: the hierarchy the basis is built on ends at `top`.
    rB.mRefined.resize(static_cast<std::size_t>(top) + 1);
    std::fill(rB.mRefined.back().begin(), rB.mRefined.back().end(), 0);
    for (int l = 0; l <= top; ++l) {
        const std::size_t s = rB.Spans(l);
        const std::size_t n = s + static_cast<std::size_t>(p);
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                const bool inside = thb_support_all(
                    p, s, i, j, [&](std::size_t A, std::size_t B) { return in_domain(l, A, B); });
                const bool deeper = thb_support_all(
                    p, s, i, j, [&](std::size_t A, std::size_t B) { return refined(l, A, B); });
                if (!inside || deeper)
                    continue;
                // Truncate level by level, on the block the function covers.
                std::size_t r0 = i, r1 = i, c0 = j, c1 = j, nn = n;
                std::vector<double> block(1, 1.0);  // (r1 - r0 + 1) x (c1 - c0 + 1), row = u
                for (int k = l; k < top; ++k) {
                    const std::size_t sk = rB.Spans(k), nk = sk + static_cast<std::size_t>(p);
                    const std::size_t nf = 2 * sk + static_cast<std::size_t>(p);
                    const std::vector<double>& a = refine[static_cast<std::size_t>(k)];
                    // The fine rows the coarse block reaches.
                    std::size_t f0 = nf, f1 = 0;
                    for (std::size_t r = 0; r < nf; ++r)
                        for (std::size_t c = r0; c <= r1; ++c)
                            if (a[r * nk + c] != 0.0) {
                                f0 = std::min(f0, r);
                                f1 = std::max(f1, r);
                            }
                    std::size_t g0 = nf, g1 = 0;
                    for (std::size_t r = 0; r < nf; ++r)
                        for (std::size_t c = c0; c <= c1; ++c)
                            if (a[r * nk + c] != 0.0) {
                                g0 = std::min(g0, r);
                                g1 = std::max(g1, r);
                            }
                    const std::size_t bu = r1 - r0 + 1, bv = c1 - c0 + 1;
                    const std::size_t fu = f1 - f0 + 1, fv = g1 - g0 + 1;
                    // next = A_u block A_v^T.
                    std::vector<double> half(fu * bv, 0.0);
                    for (std::size_t x = 0; x < fu; ++x)
                        for (std::size_t y = 0; y < bv; ++y) {
                            double sum = 0.0;
                            for (std::size_t z = 0; z < bu; ++z)
                                sum += a[(f0 + x) * nk + (r0 + z)] * block[z + bu * y];
                            half[x + fu * y] = sum;
                        }
                    std::vector<double> next(fu * fv, 0.0);
                    for (std::size_t x = 0; x < fu; ++x)
                        for (std::size_t y = 0; y < fv; ++y) {
                            double sum = 0.0;
                            for (std::size_t z = 0; z < bv; ++z)
                                sum += half[x + fu * z] * a[(g0 + y) * nk + (c0 + z)];
                            next[x + fu * y] = sum;
                        }
                    // Truncation: the finer B-splines whose support lies in
                    // Omega^{k+1} drop out.
                    const std::size_t sf = 2 * sk;
                    for (std::size_t y = 0; y < fv; ++y)
                        for (std::size_t x = 0; x < fu; ++x)
                            if (thb_support_all(p, sf, f0 + x, g0 + y,
                                                [&](std::size_t A, std::size_t B) {
                                                    return in_domain(k + 1, A, B);
                                                }))
                                next[x + fu * y] = 0.0;
                    block = std::move(next);
                    r0 = f0;
                    r1 = f1;
                    c0 = g0;
                    c1 = g1;
                    nn = nf;
                }
                ThbBasis::Function fn;
                fn.mLevel = l;
                fn.mI = i;
                fn.mJ = j;
                const std::size_t bu = r1 - r0 + 1;
                for (std::size_t y = 0; y <= c1 - c0; ++y)
                    for (std::size_t x = 0; x < bu; ++x) {
                        const double c = block[x + bu * y];
                        if (std::fabs(c) > 1e-15)
                            fn.mFine.emplace_back((r0 + x) + nn * (c0 + y), c);
                    }
                std::sort(fn.mFine.begin(), fn.mFine.end());
                rB.mFunctions.push_back(std::move(fn));
            }
    }
}

BsplineSurface thb_flatten(const ThbBasis& rB, const std::vector<PatchVec>& rCoef) {
    const std::size_t n = rB.FineSize();
    BsplineSurface s = bspline_shape(rB.mDegree, rB.mDegree, n, n, 0.0, 1.0, 0.0, 1.0);
    for (std::size_t k = 0; k < rB.mFunctions.size() && k < rCoef.size(); ++k)
        for (const auto& f : rB.mFunctions[k].mFine)
            for (std::size_t c = 0; c < 3; ++c)
                s.mPoles[3 * f.first + c] += f.second * rCoef[k][c];
    return s;
}

BsplineFitResult thb_fit_patch(const std::vector<PatchVec>& rPts,
                               const std::vector<PatchTri>& rTris, const std::vector<double>& rU,
                               const std::vector<double>& rV, const BsplineFitOptions& rOptions,
                               ThbFitInfo* pInfo) {
    BsplineFitResult out;
    const int p = rOptions.mDegree;
    if (p < 1 || !(rOptions.mTolerance > 0.0) || rPts.empty() || rTris.empty())
        return out;
    const double tol = rOptions.mTolerance;
    const double ctol = rOptions.mCentroidTolerance > 0.0 ? rOptions.mCentroidTolerance : 2.0 * tol;

    // The samples on the unit square, as `bspline_fit_patch` takes them:
    // the vertices, then the facet centres.
    double lo[2] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    std::vector<std::uint8_t> used(rPts.size(), 0);
    for (const PatchTri& t : rTris)
        for (std::size_t i : t)
            used[i] = 1;
    for (std::size_t i = 0; i < rPts.size(); ++i)
        if (used[i]) {
            lo[0] = std::min(lo[0], rU[i]);
            hi[0] = std::max(hi[0], rU[i]);
            lo[1] = std::min(lo[1], rV[i]);
            hi[1] = std::max(hi[1], rV[i]);
        }
    double ext[2] = {hi[0] - lo[0], hi[1] - lo[1]};
    if (!(ext[0] > 0.0) || !(ext[1] > 0.0))
        return out;
    for (int d = 0; d < 2; ++d) {
        lo[d] -= rOptions.mMargin * ext[d];
        ext[d] *= 1.0 + 2.0 * rOptions.mMargin;
    }
    std::vector<double> xyz, su, sv, sw;
    std::vector<std::size_t> vert_of;
    for (std::size_t i = 0; i < rPts.size(); ++i)
        if (used[i]) {
            vert_of.push_back(i);
            xyz.insert(xyz.end(), rPts[i].begin(), rPts[i].end());
            su.push_back((rU[i] - lo[0]) / ext[0]);
            sv.push_back((rV[i] - lo[1]) / ext[1]);
            sw.push_back(1.0);
        }
    const std::size_t nvert = su.size();
    const double box[4] = {
        *std::min_element(su.begin(), su.end()), *std::max_element(su.begin(), su.end()),
        *std::min_element(sv.begin(), sv.end()), *std::max_element(sv.begin(), sv.end())};
    for (const PatchTri& t : rTris) {
        for (std::size_t c = 0; c < 3; ++c)
            xyz.push_back((rPts[t[0]][c] + rPts[t[1]][c] + rPts[t[2]][c]) / 3.0);
        su.push_back(((rU[t[0]] + rU[t[1]] + rU[t[2]]) / 3.0 - lo[0]) / ext[0]);
        sv.push_back(((rV[t[0]] + rV[t[1]] + rV[t[2]]) / 3.0 - lo[1]) / ext[1]);
        sw.push_back(rOptions.mCentroidWeight);
    }
    const std::size_t ns = su.size();
    const double lambda = rOptions.mSmoothing * static_cast<double>(ns);
    const double empty = rOptions.mEmptySmoothing * static_cast<double>(ns);

    // The hierarchy: level 0 from the starting poles, as many levels as keep
    // the finest within the poles per direction.
    const std::size_t pp = static_cast<std::size_t>(p);
    const std::size_t s0 = std::max<std::size_t>(1, std::max(rOptions.mStartPoles, pp + 1) - pp);
    // (The finest level may reach twice the poles per direction: it is
    // refined only where needed, the hierarchical freedom capped below.)
    int levels = 1;
    while ((s0 << levels) + pp <= 2 * std::max(rOptions.mMaxPoles, pp + 1) && levels < 8)
        ++levels;
    ThbBasis basis = thb_hierarchy(p, s0, levels);

    // A fit at the current parameters: the THB coefficients, flattened.
    const auto fit = [&](const ThbBasis& rB, std::vector<PatchVec>& rCoef, BsplineSurface& rS) {
        const std::size_t nd = rB.mFunctions.size();
        const std::size_t nfine = rB.FineSize();
        // The THB functions each finest B-spline is part of.
        std::vector<std::vector<std::pair<std::size_t, double>>> of_fine(nfine * nfine);
        for (std::size_t k = 0; k < nd; ++k)
            for (const auto& f : rB.mFunctions[k].mFine)
                of_fine[f.first].emplace_back(k, f.second);
        const BsplineSurface shape = bspline_shape(p, p, nfine, nfine, 0.0, 1.0, 0.0, 1.0);
        std::vector<double> m(nd * nd, 0.0), rhs(3 * nd, 0.0), row(nd, 0.0);
        std::vector<std::uint8_t> mark(nd, 0);
        std::vector<std::size_t> touched, cols;
        std::vector<double> vals;
        const auto add_row = [&](double W, const double* pX) {
            for (std::size_t a : touched) {
                const double ra = row[a];
                if (pX != nullptr)
                    for (std::size_t c = 0; c < 3; ++c)
                        rhs[3 * a + c] += W * ra * pX[c];
                for (std::size_t b : touched)
                    m[a * nd + b] += W * ra * row[b];
            }
            for (std::size_t a : touched) {
                row[a] = 0.0;
                mark[a] = 0;
            }
            touched.clear();
        };
        const auto gather = [&](std::size_t Fine, double V) {
            for (const auto& t : of_fine[Fine]) {
                if (!mark[t.first]) {
                    mark[t.first] = 1;
                    touched.push_back(t.first);
                }
                row[t.first] += V * t.second;
            }
        };
        // Which finest cells hold samples.
        const std::size_t sf = nfine - pp;
        std::vector<std::uint8_t> occupied(sf * sf, 0);
        for (std::size_t k = 0; k < ns; ++k) {
            const std::size_t i = std::min(
                sf - 1, static_cast<std::size_t>(std::max(0.0, su[k]) * static_cast<double>(sf)));
            const std::size_t j = std::min(
                sf - 1, static_cast<std::size_t>(std::max(0.0, sv[k]) * static_cast<double>(sf)));
            occupied[i + sf * j] = 1;
            if (sw[k] == 0.0)
                continue;
            bspline_row(shape, su[k], sv[k], cols, vals);
            for (std::size_t c = 0; c < cols.size(); ++c)
                gather(cols[c], vals[c]);
            add_row(sw[k], &xyz[3 * k]);
        }
        // Smoothing: the finest coefficients' second differences, scaled as
        // the thin-plate energy, strong over cells no sample reaches.
        const double scale = static_cast<double>(sf) * static_cast<double>(sf);
        const auto cell_weight = [&](std::size_t I, std::size_t J) {
            const std::size_t i = std::min(I > pp / 2 ? I - pp / 2 : 0, sf - 1);
            const std::size_t j = std::min(J > pp / 2 ? J - pp / 2 : 0, sf - 1);
            return (occupied[i + sf * j] ? lambda : empty) * scale;
        };
        for (std::size_t j = 0; j < nfine; ++j)
            for (std::size_t i = 0; i < nfine; ++i) {
                const double w = cell_weight(i, j);
                if (!(w > 0.0))
                    continue;
                if (i > 0 && i + 1 < nfine) {
                    gather((i - 1) + nfine * j, 1.0);
                    gather(i + nfine * j, -2.0);
                    gather((i + 1) + nfine * j, 1.0);
                    add_row(w, nullptr);
                }
                if (j > 0 && j + 1 < nfine) {
                    gather(i + nfine * (j - 1), 1.0);
                    gather(i + nfine * j, -2.0);
                    gather(i + nfine * (j + 1), 1.0);
                    add_row(w, nullptr);
                }
                if (i + 1 < nfine && j + 1 < nfine) {
                    gather(i + nfine * j, 1.0);
                    gather((i + 1) + nfine * j, -1.0);
                    gather(i + nfine * (j + 1), -1.0);
                    gather((i + 1) + nfine * (j + 1), 1.0);
                    add_row(2.0 * w, nullptr);
                }
            }
        // The normal equations, by Cholesky; a tiny ridge for any function
        // neither the data nor the smoothing reaches.
        double diag = 0.0;
        for (std::size_t a = 0; a < nd; ++a)
            diag = std::max(diag, m[a * nd + a]);
        for (std::size_t a = 0; a < nd; ++a)
            m[a * nd + a] += 1e-12 * diag + 1e-300;
        if (!chol_factor(m, nd))
            return false;
        rCoef.assign(nd, PatchVec{0.0, 0.0, 0.0});
        std::vector<double> b(nd);
        for (std::size_t c = 0; c < 3; ++c) {
            for (std::size_t a = 0; a < nd; ++a)
                b[a] = rhs[3 * a + c];
            chol_solve(m, nd, b.data());
            for (std::size_t a = 0; a < nd; ++a) {
                if (!std::isfinite(b[a]))
                    return false;
                rCoef[a][c] = b[a];
            }
        }
        rS = thb_flatten(rB, rCoef);
        return true;
    };
    // The samples projected onto the surface, their parameters updated.
    struct Errors {
        double mVertex = 0.0, mCentroid = 0.0, mRms = 0.0;
        std::vector<double> mDist;
        double Worst(double Ratio) const { return std::max(mVertex, Ratio * mCentroid); }
    };
    const auto project = [&](const BsplineSurface& rS, std::vector<double>& rSu,
                             std::vector<double>& rSv) {
        Errors e;
        e.mDist.resize(ns);
        double sum = 0.0;
        for (std::size_t k = 0; k < ns; ++k) {
            const double d =
                bspline_project(rS, &xyz[3 * k], rSu[k], rSv[k], box[0], box[1], box[2], box[3]);
            e.mDist[k] = d;
            sum += d * d;
            if (k < nvert)
                e.mVertex = std::max(e.mVertex, d);
            else
                e.mCentroid = std::max(e.mCentroid, d);
        }
        e.mRms = std::sqrt(sum / static_cast<double>(std::max<std::size_t>(ns, 1)));
        return e;
    };
    const double ratio = tol / ctol;
    const auto met = [&](const Errors& rE) { return rE.mVertex <= tol && rE.mCentroid <= ctol; };

    thb_build(basis);
    std::vector<PatchVec> coef;
    BsplineSurface s;
    if (!fit(basis, coef, s))
        return out;
    Errors err = project(s, su, sv);
    BsplineSurface best = s;
    std::vector<double> best_u = su, best_v = sv;
    Errors best_err = err;
    std::size_t best_dofs = basis.mFunctions.size();
    int best_levels = basis.NumLevels();
    for (int it = 1; it <= rOptions.mMaxIterations && !met(best_err); ++it) {
        out.mIterations = it;
        // Parameter correction at this hierarchy, while it helps.
        {
            std::vector<PatchVec> c2;
            BsplineSurface t;
            if (fit(basis, c2, t)) {
                std::vector<double> u2 = su, v2 = sv;
                const Errors e2 = project(t, u2, v2);
                if (e2.Worst(ratio) < 0.9 * err.Worst(ratio)) {
                    s = std::move(t);
                    coef = std::move(c2);
                    su = std::move(u2);
                    sv = std::move(v2);
                    err = e2;
                    if (err.Worst(ratio) < best_err.Worst(ratio)) {
                        best = s;
                        best_u = su;
                        best_v = sv;
                        best_err = err;
                        best_dofs = basis.mFunctions.size();
                        best_levels = basis.NumLevels();
                    }
                    continue;
                }
            }
        }
        // Refine the cells that hold samples still too far, at the deepest
        // level each lies in, and their neighbours there.
        ThbBasis next = basis;
        next.mRefined.resize(static_cast<std::size_t>(levels));
        for (int l = 0; l < levels; ++l)
            if (next.mRefined[static_cast<std::size_t>(l)].size() != next.Spans(l) * next.Spans(l))
                next.mRefined[static_cast<std::size_t>(l)].assign(next.Spans(l) * next.Spans(l), 0);
        bool any = false;
        std::vector<std::vector<std::size_t>> count(static_cast<std::size_t>(levels));
        for (int l = 0; l < levels; ++l)
            count[static_cast<std::size_t>(l)].assign(next.Spans(l) * next.Spans(l), 0);
        std::vector<std::pair<int, std::size_t>> cell_of(ns);
        for (std::size_t k = 0; k < ns; ++k) {
            int l = 0;
            std::size_t ci = 0;
            for (;; ++l) {
                const std::size_t sl = next.Spans(l);
                const std::size_t i = std::min(
                    sl - 1,
                    static_cast<std::size_t>(std::max(0.0, su[k]) * static_cast<double>(sl)));
                const std::size_t j = std::min(
                    sl - 1,
                    static_cast<std::size_t>(std::max(0.0, sv[k]) * static_cast<double>(sl)));
                ci = i + sl * j;
                if (l + 1 >= levels || !next.mRefined[static_cast<std::size_t>(l)][ci])
                    break;
            }
            cell_of[k] = {l, ci};
            ++count[static_cast<std::size_t>(l)][ci];
        }
        // The worst first: a cell is refined only for a sample at least half
        // as far out of its target as the worst, so the hierarchy grows
        // where the misfit is and the fit is redone before the rest follow.
        double worst = 0.0;
        for (std::size_t k = 0; k < ns; ++k)
            worst = std::max(worst, err.mDist[k] / (k < nvert ? tol : ctol));
        for (std::size_t k = 0; k < ns; ++k) {
            const double limit = k < nvert ? tol : ctol;
            if (!(err.mDist[k] > limit) || err.mDist[k] / limit < 0.5 * worst)
                continue;
            const int l = cell_of[k].first;
            const std::size_t ci = cell_of[k].second;
            if (l + 1 >= levels || count[static_cast<std::size_t>(l)][ci] < kThbMinCellSamples)
                continue;
            // With its neighbours: a level-`l + 1` cubic needs four of its
            // cells a side, two refined cells of level `l`.
            const std::size_t sl = next.Spans(l);
            const std::size_t i = ci % sl, j = ci / sl;
            for (std::size_t jj = j > 0 ? j - 1 : 0; jj <= std::min(sl - 1, j + 1); ++jj)
                for (std::size_t ii = i > 0 ? i - 1 : 0; ii <= std::min(sl - 1, i + 1); ++ii) {
                    const bool inside =
                        l == 0 || next.mRefined[static_cast<std::size_t>(l - 1)]
                                               [(ii / 2) + next.Spans(l - 1) * (jj / 2)] != 0;
                    std::uint8_t& r = next.mRefined[static_cast<std::size_t>(l)][ii + sl * jj];
                    if (inside && !r) {
                        r = 1;
                        any = true;
                    }
                }
        }
        if (!any) {
            log::debug("thb: no cell left to refine at {} function(s), {} level(s), worst {}",
                       basis.mFunctions.size(), basis.NumLevels(), err.Worst(ratio));
            break;
        }
        thb_build(next);
        if (rOptions.mMaxTotalPoles > 0 && next.mFunctions.size() > rOptions.mMaxTotalPoles) {
            log::debug("thb: {} function(s) would pass the budget of {}; worst {}",
                       next.mFunctions.size(), rOptions.mMaxTotalPoles, err.Worst(ratio));
            break;  // the budget is spent
        }
        std::vector<PatchVec> c2;
        BsplineSurface t;
        if (!fit(next, c2, t))
            break;
        basis = std::move(next);
        coef = std::move(c2);
        s = std::move(t);
        err = project(s, su, sv);
        ++out.mRefinements;
        if (err.Worst(ratio) < best_err.Worst(ratio)) {
            best = s;
            best_u = su;
            best_v = sv;
            best_err = err;
            best_dofs = basis.mFunctions.size();
            best_levels = basis.NumLevels();
        }
    }
    out.mOk = met(best_err);
    out.mSurface = std::move(best);
    out.mMaxVertexError = best_err.mVertex;
    out.mMaxCentroidError = best_err.mCentroid;
    out.mRmsError = best_err.mRms;
    out.mU.assign(rPts.size(), 0.0);
    out.mV.assign(rPts.size(), 0.0);
    for (std::size_t k = 0; k < vert_of.size(); ++k) {
        out.mU[vert_of[k]] = best_u[k];
        out.mV[vert_of[k]] = best_v[k];
    }
    if (pInfo != nullptr) {
        pInfo->mNumDofs = best_dofs;
        pInfo->mNumLevels = best_levels;
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
