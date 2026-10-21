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
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "bspline_fit.hpp"
#include "dense_chol.hpp"
#include "nurbs_eval.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// TDM steps couple the three coordinates, a system three times the poles':
// only taken while the surface has at most this many poles.
constexpr std::size_t kBsfTdmMaxPoles = 400;
// The fewest samples a span cell needs before it may be halved.
constexpr std::size_t kBsfMinCellSamples = 8;

// The samples: vertices then centroids, their points and parameters.
struct BsfSamples {
    std::vector<double> mXyz;
    std::vector<double> mU, mV;
    std::vector<double> mWeight;
    std::size_t mNumVertices = 0;
    std::size_t mNumData = 0;               // vertices and centroids; ghost samples follow
    double mBox[4] = {0.0, 1.0, 0.0, 1.0};  // the samples' (u, v) box

    std::size_t Size() const { return mU.size(); }
};

// The largest distances of the vertices and of the centroids, and the RMS.
struct BsfErrors {
    double mVertex = 0.0;
    double mCentroid = 0.0;
    double mRms = 0.0;
    std::vector<double> mDist;

    double mRatio = 0.5;  // the vertex tolerance over the centroid one

    double Worst() const { return std::max(mVertex, mRatio * mCentroid); }
};

// Project every sample onto the surface, updating its parameters, the foot
// points held inside the samples' own box: a boundary sample sliding into the
// margin would claim spans the extension needs to keep smooth.
BsfErrors bsf_project(const BsplineSurface& rS, BsfSamples& rSamples, double Ratio) {
    BsfErrors e;
    e.mRatio = Ratio;
    e.mDist.resize(rSamples.mNumData);
    double sum = 0.0;
    for (std::size_t k = 0; k < rSamples.mNumData; ++k) {
        const double d =
            bspline_project(rS, &rSamples.mXyz[3 * k], rSamples.mU[k], rSamples.mV[k],
                            rSamples.mBox[0], rSamples.mBox[1], rSamples.mBox[2], rSamples.mBox[3]);
        e.mDist[k] = d;
        sum += d * d;
        if (k < rSamples.mNumVertices)
            e.mVertex = std::max(e.mVertex, d);
        else
            e.mCentroid = std::max(e.mCentroid, d);
    }
    e.mRms = rSamples.mNumData ? std::sqrt(sum / static_cast<double>(rSamples.mNumData)) : 0.0;
    return e;
}

// A shape with the given knot vectors and zero poles.
BsplineSurface bsf_shape(int P, std::vector<double> rKu, std::vector<double> rKv) {
    BsplineSurface s;
    s.mDegreeU = s.mDegreeV = P;
    s.mNu = rKu.size() - static_cast<std::size_t>(P) - 1;
    s.mNv = rKv.size() - static_cast<std::size_t>(P) - 1;
    s.mKnotsU = std::move(rKu);
    s.mKnotsV = std::move(rKv);
    s.mPoles.assign(s.mNu * s.mNv * 3, 0.0);
    return s;
}

// Point distance minimisation: the poles fitted at the samples' parameters,
// the thin-plate weight raised tenfold up to three times if the system is
// singular.
// Per knot span its thin-plate weight: `Lambda` where samples are, `Empty`
// where none is.
std::vector<double> bsf_span_weights(const BsplineSurface& rS, const BsfSamples& rSamples,
                                     double Lambda, double Empty) {
    const std::vector<std::size_t> su = nurbs_spans(rS.mDegreeU, rS.mKnotsU, rS.mNu);
    const std::vector<std::size_t> sv = nurbs_spans(rS.mDegreeV, rS.mKnotsV, rS.mNv);
    std::vector<double> w(su.size() * sv.size(), Empty);
    for (std::size_t k = 0; k < rSamples.Size(); ++k) {
        const double u = std::min(std::max(rSamples.mU[k], rS.U0()), rS.U1());
        const double v = std::min(std::max(rSamples.mV[k], rS.V0()), rS.V1());
        const std::size_t iu = static_cast<std::size_t>(
            std::lower_bound(su.begin(), su.end(),
                             nurbs_find_span(rS.mDegreeU, rS.mKnotsU, rS.mNu, u)) -
            su.begin());
        const std::size_t iv = static_cast<std::size_t>(
            std::lower_bound(sv.begin(), sv.end(),
                             nurbs_find_span(rS.mDegreeV, rS.mKnotsV, rS.mNv, v)) -
            sv.begin());
        if (iu < su.size() && iv < sv.size())
            w[iu + su.size() * iv] = Lambda;
    }
    return w;
}

bool bsf_pdm(BsplineSurface& rS, const BsfSamples& rSamples, double Lambda, double Empty) {
    for (int attempt = 0; attempt < 4; ++attempt) {
        const std::vector<double> w = bsf_span_weights(rS, rSamples, Lambda, Empty);
        if (bspline_lsq_fit(rS, rSamples.mU, rSamples.mV, rSamples.mXyz, &rSamples.mWeight, Lambda,
                            &w))
            return true;
        Lambda = Lambda > 0.0 ? Lambda * 10.0 : 1e-9 * static_cast<double>(rSamples.Size());
    }
    return false;
}

// One tangent distance minimisation step with damping Mu (Wang, Pottmann and
// Liu 2006): min sum ((S(u_k) - x_k) . n_k)^2 + Mu |S(u_k) - x_k|^2
// + Lambda * thin plate, the three coordinates coupled.
bool bsf_tdm(BsplineSurface& rS, const BsfSamples& rSamples, double Lambda, double Empty,
             double Mu) {
    const std::size_t np = rS.NumPoles();
    const SkylineMatrix pole_env = bspline_normal_envelope(rS);
    std::vector<std::size_t> first(3 * np);
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t c = 0; c < 3; ++c)
            first[3 * i + c] = 3 * pole_env.mFirst[i];
    SkylineMatrix a = skyline_zero(first);
    std::vector<double> rhs(3 * np, 0.0);
    std::vector<std::size_t> cols;
    std::vector<double> vals;
    for (std::size_t k = 0; k < rSamples.Size(); ++k) {
        BsplineEval e;
        bspline_eval(rS, rSamples.mU[k], rSamples.mV[k], 1, e);
        double n[3] = {e.mSu[1] * e.mSv[2] - e.mSu[2] * e.mSv[1],
                       e.mSu[2] * e.mSv[0] - e.mSu[0] * e.mSv[2],
                       e.mSu[0] * e.mSv[1] - e.mSu[1] * e.mSv[0]};
        const double ln = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (ln > 0.0)
            for (double& x : n)
                x /= ln;
        const double* x = &rSamples.mXyz[3 * k];
        const double xn = x[0] * n[0] + x[1] * n[1] + x[2] * n[2];
        const double w = rSamples.mWeight[k];
        if (w == 0.0)
            continue;
        bspline_row(rS, rSamples.mU[k], rSamples.mV[k], cols, vals);
        for (std::size_t ia = 0; ia < cols.size(); ++ia) {
            for (std::size_t c = 0; c < 3; ++c)
                rhs[3 * cols[ia] + c] += w * vals[ia] * (n[c] * xn + Mu * x[c]);
            for (std::size_t ib = 0; ib < cols.size(); ++ib) {
                if (cols[ib] > cols[ia])
                    continue;
                const double bb = w * vals[ia] * vals[ib];
                for (std::size_t c = 0; c < 3; ++c)
                    for (std::size_t d = 0; d < 3; ++d) {
                        const std::size_t r = 3 * cols[ia] + c, s = 3 * cols[ib] + d;
                        if (s > r)
                            continue;
                        a.At(r, s) += bb * (n[c] * n[d] + (c == d ? Mu : 0.0));
                    }
            }
        }
    }
    {
        SkylineMatrix k = bspline_normal_envelope(rS);
        bspline_add_thin_plate_spans(rS, bsf_span_weights(rS, rSamples, Lambda, Empty), k);
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = k.mFirst[i]; j <= i; ++j)
                for (std::size_t c = 0; c < 3; ++c)
                    a.At(3 * i + c, 3 * j + c) += k.At(i, j);
    }
    if (!skyline_ldlt(a))
        return false;
    skyline_solve(a, rhs.data());
    for (double v : rhs)
        if (!std::isfinite(v))
            return false;
    rS.mPoles = std::move(rhs);
    return true;
}

// The knot vector with the midpoint of each marked span inserted.
std::vector<double> bsf_insert(const std::vector<double>& rK, int P,
                               const std::vector<std::size_t>& rSpans) {
    std::vector<double> mids;
    for (std::size_t s : rSpans)
        mids.push_back(0.5 * (rK[s] + rK[s + 1]));
    std::vector<double> out = rK;
    out.insert(out.end(), mids.begin(), mids.end());
    std::sort(out.begin(), out.end());
    (void)P;
    return out;
}

}  // namespace

BsplineFitResult bspline_fit_patch(const std::vector<PatchVec>& rPts,
                                   const std::vector<PatchTri>& rTris,
                                   const std::vector<double>& rU, const std::vector<double>& rV,
                                   const BsplineFitOptions& rOptions) {
    BsplineFitResult out;
    const int p = rOptions.mDegree;
    if (p < 1 || !(rOptions.mTolerance > 0.0) || rPts.empty() || rTris.empty())
        return out;
    const double tol = rOptions.mTolerance;
    const double ctol =
        rOptions.mCentroidTolerance > 0.0 ? rOptions.mCentroidTolerance : 2.0 * rOptions.mTolerance;

    // The parameters' box, widened, mapped onto the unit square.
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
    BsfSamples smp;
    std::vector<std::size_t> vert_of;
    for (std::size_t i = 0; i < rPts.size(); ++i)
        if (used[i]) {
            vert_of.push_back(i);
            smp.mXyz.insert(smp.mXyz.end(), rPts[i].begin(), rPts[i].end());
            smp.mU.push_back((rU[i] - lo[0]) / ext[0]);
            smp.mV.push_back((rV[i] - lo[1]) / ext[1]);
            smp.mWeight.push_back(1.0);
        }
    smp.mNumVertices = smp.mU.size();
    smp.mBox[0] = *std::min_element(smp.mU.begin(), smp.mU.end());
    smp.mBox[1] = *std::max_element(smp.mU.begin(), smp.mU.end());
    smp.mBox[2] = *std::min_element(smp.mV.begin(), smp.mV.end());
    smp.mBox[3] = *std::max_element(smp.mV.begin(), smp.mV.end());
    for (const PatchTri& t : rTris) {
        for (std::size_t c = 0; c < 3; ++c)
            smp.mXyz.push_back((rPts[t[0]][c] + rPts[t[1]][c] + rPts[t[2]][c]) / 3.0);
        smp.mU.push_back(((rU[t[0]] + rU[t[1]] + rU[t[2]]) / 3.0 - lo[0]) / ext[0]);
        smp.mV.push_back(((rV[t[0]] + rV[t[1]] + rV[t[2]]) / 3.0 - lo[1]) / ext[1]);
        smp.mWeight.push_back(rOptions.mCentroidWeight);
    }
    smp.mNumData = smp.Size();
    const double lambda = rOptions.mSmoothing * static_cast<double>(smp.Size());
    const double empty = rOptions.mEmptySmoothing * static_cast<double>(smp.Size());

    // Ghost samples where the parameter square holds no data (the margin, a
    // disc's corners): a coarse, smooth fit evaluated there, lightly weighted,
    // so the surface runs on past its patch as that smooth surface does
    // instead of wherever poles no sample pins fall.
    {
        const std::size_t nc = std::min<std::size_t>(
            6, std::max(rOptions.mStartPoles, std::size_t(static_cast<std::size_t>(p) + 1)));
        BsplineSurface coarse = bsf_shape(p, bspline_clamped_knots(p, nc, 0.0, 1.0),
                                          bspline_clamped_knots(p, nc, 0.0, 1.0));
        if (bspline_lsq_fit(coarse, smp.mU, smp.mV, smp.mXyz, &smp.mWeight,
                            1e-6 * static_cast<double>(smp.Size()))) {
            // Covered: every cell within the data's own spacing (the median
            // edge length in (u, v)) of a sample, so a gap between samples
            // inside the patch never takes a ghost.
            constexpr std::size_t kG = 128;
            std::vector<double> lens;
            for (const PatchTri& tr : rTris)
                for (int c = 0; c < 3; ++c) {
                    const std::size_t a = tr[static_cast<std::size_t>(c)];
                    const std::size_t b = tr[static_cast<std::size_t>((c + 1) % 3)];
                    lens.push_back(std::hypot((rU[a] - rU[b]) / ext[0], (rV[a] - rV[b]) / ext[1]));
                }
            std::nth_element(lens.begin(),
                             lens.begin() + static_cast<std::ptrdiff_t>(lens.size() / 2),
                             lens.end());
            const double h = lens[lens.size() / 2];
            const std::size_t reach = static_cast<std::size_t>(std::ceil(h * kG)) + 1;
            std::vector<std::uint8_t> occupied(kG * kG, 0);
            for (std::size_t k = 0; k < smp.mNumData; ++k) {
                const std::size_t i =
                    std::min(kG - 1, static_cast<std::size_t>(std::max(0.0, smp.mU[k]) * kG));
                const std::size_t j =
                    std::min(kG - 1, static_cast<std::size_t>(std::max(0.0, smp.mV[k]) * kG));
                for (std::size_t jj = j > reach ? j - reach : 0; jj <= std::min(kG - 1, j + reach);
                     ++jj)
                    for (std::size_t ii = i > reach ? i - reach : 0;
                         ii <= std::min(kG - 1, i + reach); ++ii)
                        occupied[ii + kG * jj] = 1;
            }
            for (std::size_t j = 0; j < kG; ++j)
                for (std::size_t i = 0; i < kG; ++i) {
                    if (occupied[i + kG * j])
                        continue;
                    const double u = (static_cast<double>(i) + 0.5) / kG;
                    const double v = (static_cast<double>(j) + 0.5) / kG;
                    BsplineEval e;
                    bspline_eval(coarse, u, v, 0, e);
                    smp.mXyz.insert(smp.mXyz.end(), e.mS, e.mS + 3);
                    smp.mU.push_back(u);
                    smp.mV.push_back(v);
                    smp.mWeight.push_back(rOptions.mGhostWeight);
                }
        }
    }

    // The start: uniform knots, as many poles as the samples and the budget allow.
    const std::size_t pmin = static_cast<std::size_t>(p) + 1;
    const std::size_t maxp = std::max(rOptions.mMaxPoles, pmin);
    std::size_t n0 = std::min(std::max(rOptions.mStartPoles, pmin), maxp);
    BsplineSurface s = bsf_shape(p, bspline_clamped_knots(p, n0, 0.0, 1.0),
                                 bspline_clamped_knots(p, n0, 0.0, 1.0));
    if (!bsf_pdm(s, smp, lambda, empty))
        return out;
    BsfErrors err = bsf_project(s, smp, tol / ctol);
    out.mHistory.push_back(err.Worst());
    const auto met = [&](const BsfErrors& rE) { return rE.mVertex <= tol && rE.mCentroid <= ctol; };

    // The working state moves on; the best one seen is the result.
    BsplineSurface best = s;
    BsfSamples best_smp = smp;
    BsfErrors best_err = err;
    const auto keep_best = [&]() {
        if (err.Worst() < best_err.Worst()) {
            best = s;
            best_smp = smp;
            best_err = err;
        }
    };
    double mu = 0.1;
    for (int it = 1; it <= rOptions.mMaxIterations && !met(best_err); ++it) {
        out.mIterations = it;
        const double before = err.Worst();
        bool progressed = false;
        // 1. Point distance minimisation at the corrected parameters.
        {
            BsplineSurface t = s;
            BsfSamples ts = smp;
            if (bsf_pdm(t, ts, lambda, empty)) {
                const BsfErrors te = bsf_project(t, ts, tol / ctol);
                if (te.Worst() < err.Worst()) {
                    progressed = te.Worst() < 0.9 * before;
                    s = std::move(t);
                    smp = std::move(ts);
                    err = te;
                }
            }
        }
        // 2. Tangent distance minimisation once that stalls.
        if (!progressed && !met(err) && s.NumPoles() <= kBsfTdmMaxPoles) {
            for (int attempt = 0; attempt < 1 && !progressed; ++attempt) {
                BsplineSurface t = s;
                BsfSamples ts = smp;
                if (!bsf_tdm(t, ts, lambda, empty, mu)) {
                    mu *= 10.0;
                    continue;
                }
                const BsfErrors te = bsf_project(t, ts, tol / ctol);
                if (te.Worst() < err.Worst()) {
                    progressed = te.Worst() < 0.9 * before;
                    s = std::move(t);
                    smp = std::move(ts);
                    err = te;
                    mu = std::max(mu * 0.5, 1e-2);
                } else {
                    mu *= 10.0;
                }
            }
        }
        // 3. Knot insertion where the samples are still too far.
        if (!progressed && !met(err)) {
            // A span is halved only while its cell holds enough samples: a
            // spline finer than its data is free between them and oscillates.
            std::map<std::pair<std::size_t, std::size_t>, std::size_t> in_cell;
            std::vector<std::pair<std::size_t, std::size_t>> cell_of(smp.mNumData);
            for (std::size_t k = 0; k < smp.mNumData; ++k) {
                const double u = std::min(std::max(smp.mU[k], 0.0), 1.0);
                const double v = std::min(std::max(smp.mV[k], 0.0), 1.0);
                cell_of[k] = {nurbs_find_span(p, s.mKnotsU, s.mNu, u),
                              nurbs_find_span(p, s.mKnotsV, s.mNv, v)};
                ++in_cell[cell_of[k]];
            }
            std::map<std::size_t, double> bad_u, bad_v;
            for (std::size_t k = 0; k < smp.mNumData; ++k) {
                const double limit = k < smp.mNumVertices ? tol : ctol;
                if (!(err.mDist[k] > limit) || in_cell[cell_of[k]] < kBsfMinCellSamples)
                    continue;
                const std::size_t su = cell_of[k].first, sv = cell_of[k].second;
                bad_u[su] = std::max(bad_u[su], err.mDist[k]);
                bad_v[sv] = std::max(bad_v[sv], err.mDist[k]);
            }
            // The worst spans first, within the pole budget.
            const auto pick = [&](const std::map<std::size_t, double>& rBad, std::size_t N) {
                std::vector<std::pair<double, std::size_t>> by;
                for (const auto& kv : rBad)
                    by.emplace_back(-kv.second, kv.first);
                std::sort(by.begin(), by.end());
                std::vector<std::size_t> spans;
                for (const auto& b : by)
                    if (N + spans.size() < maxp)
                        spans.push_back(b.second);
                return spans;
            };
            const std::vector<std::size_t> su = pick(bad_u, s.mNu);
            const std::vector<std::size_t> sv = pick(bad_v, s.mNv);
            if (su.empty() && sv.empty())
                break;  // the pole budget is spent
            BsplineSurface t =
                bsf_shape(p, bsf_insert(s.mKnotsU, p, su), bsf_insert(s.mKnotsV, p, sv));
            BsfSamples ts = smp;
            if (!bsf_pdm(t, ts, lambda, empty))
                break;
            err = bsf_project(t, ts, tol / ctol);
            s = std::move(t);
            smp = std::move(ts);
            ++out.mRefinements;
            progressed = true;
        }
        keep_best();
        out.mHistory.push_back(best_err.Worst());
        if (!progressed)
            break;
    }
    s = std::move(best);
    smp = std::move(best_smp);
    err = best_err;

    out.mOk = met(err);
    out.mSurface = std::move(s);
    out.mMaxVertexError = err.mVertex;
    out.mMaxCentroidError = err.mCentroid;
    out.mRmsError = err.mRms;
    out.mU.assign(rPts.size(), 0.0);
    out.mV.assign(rPts.size(), 0.0);
    for (std::size_t k = 0; k < vert_of.size(); ++k) {
        out.mU[vert_of[k]] = smp.mU[k];
        out.mV[vert_of[k]] = smp.mV[k];
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
