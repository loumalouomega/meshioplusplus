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
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/refine_patches.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/crease_edges.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/primitive_geom.hpp"
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kRpPrefix = "meshio++: refine_patches: ";
constexpr double kRpPi = 3.14159265358979323846;
// RANSAC hypotheses per freeform patch and iteration, and the cells a
// hypothesis is fitted to.
constexpr int kRpHypotheses = 24;
constexpr std::size_t kRpNeighbourhood = 64;
// fit_primitives' default share of a patch's points within the tolerance.
constexpr double kRpMinInlierFraction = 0.9;
// Relative slack on the energy comparisons (Lloyd step, merge).
constexpr double kRpEnergySlack = 1e-12;

// splitmix64 (Steele, Lea and Flood 2014): the same sequence on every
// platform, unlike a standard-library distribution.
struct RpRng {
    std::uint64_t mState;
    std::uint64_t Next() {
        std::uint64_t z = (mState += 0x9e3779b97f4a7c15ull);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        return z ^ (z >> 31);
    }
};

// The surface as faces: rings, geometry, smooth neighbours and paired edges.
struct RpSurface {
    detail::SurfaceRings mRings;
    std::vector<double> mXyz;
    std::vector<double> mNormal;    // unit, 3 per face
    std::vector<double> mCentroid;  // 3 per face
    std::vector<double> mArea;
    std::vector<std::int64_t> mLabel;    // input label per face, -1 outside every patch
    std::vector<std::size_t> mAdjStart;  // smooth neighbours, CSR
    std::vector<std::int64_t> mAdj;      // ascending per face
    std::vector<std::pair<std::int64_t, std::int64_t>> mPairs;  // every two-face edge (a < b)
    double mDiag = 0.0;
    double mTol = 0.0;
    double mWeight = 1.0;

    std::size_t NumFaces() const { return mRings.NumFaces(); }
};

// A patch's proxy: a primitive, a quadric height field or nothing (facets).
struct RpProxy {
    RefineProxy mKind = RefineProxy::Facets;
    detail::PrimitiveFit mFit;
    double mO[3] = {0, 0, 0};
    double mE1[3] = {1, 0, 0};
    double mE2[3] = {0, 1, 0};
    double mN[3] = {0, 0, 1};
    double mQ[6] = {0, 0, 0, 0, 0, 0};
    double mRms = 0.0;
    double mMax = 0.0;
};

double rp_dot(const double* pA, const double* pB) {
    return pA[0] * pB[0] + pA[1] * pB[1] + pA[2] * pB[2];
}

// The distance of a point from the proxy and, when pN is given, the proxy's
// unit normal there (false: no normal, the normal term is then zero).
double rp_distance(const RpProxy& rP, const double* pX, double* pN, bool* pHasN) {
    if (pHasN)
        *pHasN = false;
    if (rP.mKind == RefineProxy::Facets)
        return 0.0;
    if (rP.mKind != RefineProxy::Freeform) {
        const double d = detail::primitive_point_distance(rP.mFit.mType, rP.mFit.mParams, pX);
        if (pN && pHasN)
            *pHasN = detail::primitive_normal(rP.mFit.mType, rP.mFit.mParams, pX, pN);
        return d;
    }
    const double r[3] = {pX[0] - rP.mO[0], pX[1] - rP.mO[1], pX[2] - rP.mO[2]};
    const double x = rp_dot(r, rP.mE1), y = rp_dot(r, rP.mE2), h = rp_dot(r, rP.mN);
    const double* q = rP.mQ;
    const double z = q[0] + q[1] * x + q[2] * y + q[3] * x * x + q[4] * x * y + q[5] * y * y;
    const double gx = q[1] + 2.0 * q[3] * x + q[4] * y, gy = q[2] + q[4] * x + 2.0 * q[5] * y;
    const double g = std::sqrt(1.0 + gx * gx + gy * gy);
    if (pN && pHasN) {
        for (int k = 0; k < 3; ++k)
            pN[k] = (rP.mN[k] - gx * rP.mE1[k] - gy * rP.mE2[k]) / g;
        *pHasN = true;
    }
    // The vertical residual over the slope: the first-order orthogonal distance.
    return std::abs(h - z) / g;
}

// A face's cost against a proxy: area * (d^2 / t^2 + w (1 - |n . n_p|)).
double rp_cost(const RpSurface& rS, std::size_t F, const RpProxy& rP) {
    if (rP.mKind == RefineProxy::Facets)
        return 0.0;
    double n[3];
    bool has_n = false;
    const double d = rp_distance(rP, &rS.mCentroid[F * 3], n, &has_n);
    double c = d * d / (rS.mTol * rS.mTol);
    if (has_n)
        c += rS.mWeight * std::max(0.0, 1.0 - std::abs(rp_dot(n, &rS.mNormal[F * 3])));
    return rS.mArea[F] * c;
}

// The distinct corners, normals and centroids of a set of faces, the samples
// fit_primitives fits.
void rp_samples(const RpSurface& rS, const std::vector<std::int64_t>& rFaces,
                std::vector<double>& rPts, std::vector<double>& rNormals,
                std::vector<double>& rCentres) {
    const std::vector<std::int64_t>& start = rS.mRings.mStart;
    const std::vector<std::int64_t>& nodes = rS.mRings.mNodes;
    std::vector<std::int64_t> verts;
    for (std::int64_t f : rFaces)
        for (std::int64_t k = start[static_cast<std::size_t>(f)];
             k < start[static_cast<std::size_t>(f) + 1]; ++k)
            verts.push_back(nodes[static_cast<std::size_t>(k)]);
    std::sort(verts.begin(), verts.end());
    verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
    rPts.resize(verts.size() * 3);
    for (std::size_t i = 0; i < verts.size(); ++i)
        for (std::size_t d = 0; d < 3; ++d)
            rPts[i * 3 + d] = rS.mXyz[static_cast<std::size_t>(verts[i]) * 3 + d];
    rNormals.clear();
    rCentres.clear();
    for (std::int64_t f : rFaces)
        for (std::size_t d = 0; d < 3; ++d) {
            rNormals.push_back(rS.mNormal[static_cast<std::size_t>(f) * 3 + d]);
            rCentres.push_back(rS.mCentroid[static_cast<std::size_t>(f) * 3 + d]);
        }
}

// The largest and RMS distance of the points from the proxy.
void rp_errors(const RpProxy& rP, const std::vector<double>& rPts, double& rRms, double& rMax) {
    double s = 0.0, m = 0.0;
    const std::size_t n = rPts.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const double d = rp_distance(rP, &rPts[i * 3], nullptr, nullptr);
        s += d * d;
        m = std::max(m, d);
    }
    rRms = n ? std::sqrt(s / static_cast<double>(n)) : 0.0;
    rMax = m;
}

// A primitive the samples lie on as `fit_primitives` accepts one (RMS within
// the tolerance, nine points in ten within it), so the refinement's verdict
// and the fit table it writes agree, and not a curved surface so flat or so
// large it stands in for a plane (segment's sanity rule). The points outside
// the tolerance are cells the primitive does not take in the Lloyd step.
bool rp_fit_primitive(const RpSurface& rS, const std::vector<double>& rPts,
                      const std::vector<double>& rNormals, const std::vector<double>& rCentres,
                      std::uint32_t Types, detail::PrimitiveFit& rFit) {
    if (rPts.size() < 9 || Types == 0)
        return false;
    rFit = detail::fit_primitive_samples(rPts, rNormals, rCentres, rS.mTol, Types,
                                         kRpMinInlierFraction, 50);
    if (rFit.mType <= 0)
        return false;
    for (int k = 0; k < detail::kPrimitiveParams; ++k)
        if (!std::isfinite(rFit.mParams[k]))
            return false;
    if (rFit.mType == 4)
        return rFit.mParams[8] < 85.0;
    if (rFit.mType >= 2)
        return rFit.mParams[6] <= rS.mDiag && rFit.mParams[7] <= rS.mDiag;
    return true;
}

// A quadric height field over the samples' best-fit plane.
void rp_fit_quadric(const RpSurface& rS, const std::vector<std::int64_t>& rFaces,
                    const std::vector<double>& rPts, RpProxy& rP) {
    const std::size_t n = rPts.size() / 3;
    double c[3] = {0, 0, 0};
    for (std::size_t i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            c[k] += rPts[i * 3 + static_cast<std::size_t>(k)];
    for (double& v : c)
        v /= static_cast<double>(std::max<std::size_t>(n, 1));
    double cov[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (std::size_t i = 0; i < n; ++i)
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                cov[a * 3 + b] += (rPts[i * 3 + static_cast<std::size_t>(a)] - c[a]) *
                                  (rPts[i * 3 + static_cast<std::size_t>(b)] - c[b]);
    double w[3], v[9];
    detail::lsq_sym3_eigen(cov, w, v);
    double nrm[3] = {v[0], v[1], v[2]};
    // Orient the plane's normal by the faces' area-weighted normal.
    double mean[3] = {0, 0, 0};
    for (std::int64_t f : rFaces)
        for (int k = 0; k < 3; ++k)
            mean[k] += rS.mArea[static_cast<std::size_t>(f)] *
                       rS.mNormal[static_cast<std::size_t>(f) * 3 + static_cast<std::size_t>(k)];
    if (rp_dot(mean, nrm) < 0.0)
        for (double& x : nrm)
            x = -x;
    rP.mKind = RefineProxy::Freeform;
    for (int k = 0; k < 3; ++k) {
        rP.mO[k] = c[k];
        rP.mN[k] = nrm[k];
    }
    detail::lsq_frame(rP.mN, rP.mE1, rP.mE2);
    std::fill(rP.mQ, rP.mQ + 6, 0.0);
    if (n < 6)
        return;  // a plane through the points
    // Least squares h = q0 + q1 x + q2 y + q3 x^2 + q4 xy + q5 y^2, in
    // coordinates scaled by the patch's extent for conditioning.
    double ext = 0.0;
    std::vector<double> xs(n), ys(n), hs(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double r[3] = {rPts[i * 3] - c[0], rPts[i * 3 + 1] - c[1], rPts[i * 3 + 2] - c[2]};
        xs[i] = rp_dot(r, rP.mE1);
        ys[i] = rp_dot(r, rP.mE2);
        hs[i] = rp_dot(r, rP.mN);
        ext = std::max(ext, std::max(std::abs(xs[i]), std::abs(ys[i])));
    }
    if (!(ext > 0.0))
        return;
    std::vector<double> a(36, 0.0), b(6, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = xs[i] / ext, y = ys[i] / ext;
        const double row[6] = {1.0, x, y, x * x, x * y, y * y};
        for (int r = 0; r < 6; ++r) {
            b[static_cast<std::size_t>(r)] += row[r] * hs[i];
            for (int s = 0; s < 6; ++s)
                a[static_cast<std::size_t>(r * 6 + s)] += row[r] * row[s];
        }
    }
    for (int r = 0; r < 6; ++r)
        a[static_cast<std::size_t>(r * 7)] *= 1.0 + 1e-12;  // a whisker of ridge
    if (!detail::lsq_solve(a, b, 6))
        return;
    // Back to unscaled coordinates.
    rP.mQ[0] = b[0];
    rP.mQ[1] = b[1] / ext;
    rP.mQ[2] = b[2] / ext;
    rP.mQ[3] = b[3] / (ext * ext);
    rP.mQ[4] = b[4] / (ext * ext);
    rP.mQ[5] = b[5] / (ext * ext);
}

// The proxy of a set of faces: a primitive within the tolerance when one is,
// else a quadric height field, else (too few points) facets. When the fresh
// fit finds no primitive, `pKeep` -- the primitive the faces were grown or
// found by -- stays the proxy if it still explains every point.
RpProxy rp_fit(const RpSurface& rS, const std::vector<std::int64_t>& rFaces, std::uint32_t Types,
               const RpProxy* pKeep = nullptr) {
    RpProxy p;
    std::vector<double> pts, nrm, ctr;
    rp_samples(rS, rFaces, pts, nrm, ctr);
    if (pts.size() < 9) {
        p.mKind = RefineProxy::Facets;
        return p;
    }
    if (rp_fit_primitive(rS, pts, nrm, ctr, Types, p.mFit)) {
        p.mKind = static_cast<RefineProxy>(p.mFit.mType);
    } else if (pKeep && pKeep->mKind != RefineProxy::Facets &&
               pKeep->mKind != RefineProxy::Freeform) {
        p = *pKeep;
        rp_errors(p, pts, p.mRms, p.mMax);
        if (p.mMax <= rS.mTol)
            return p;
        p = RpProxy{};
        rp_fit_quadric(rS, rFaces, pts, p);
    } else {
        rp_fit_quadric(rS, rFaces, pts, p);
    }
    rp_errors(p, pts, p.mRms, p.mMax);
    return p;
}

bool rp_is_primitive(const RpProxy& rP) {
    return rP.mKind != RefineProxy::Facets && rP.mKind != RefineProxy::Freeform;
}

// Whether every corner of a face lies within the tolerance of a primitive
// proxy (always true of a freeform or facets proxy, which take any cell).
bool rp_takes(const RpSurface& rS, std::size_t F, const RpProxy& rP) {
    if (!rp_is_primitive(rP))
        return true;
    const std::vector<std::int64_t>& start = rS.mRings.mStart;
    for (std::int64_t k = start[F]; k < start[F + 1]; ++k) {
        const double* p =
            &rS.mXyz[static_cast<std::size_t>(rS.mRings.mNodes[static_cast<std::size_t>(k)]) * 3];
        if (!(detail::primitive_point_distance(rP.mFit.mType, rP.mFit.mParams, p) <= rS.mTol))
            return false;
    }
    return true;
}

double rp_energy(const RpSurface& rS, const std::vector<std::int64_t>& rFaces, const RpProxy& rP) {
    double e = 0.0;
    for (std::int64_t f : rFaces)
        e += rp_cost(rS, static_cast<std::size_t>(f), rP);
    return e;
}

// The working partition: per face its patch, per patch its faces and proxy.
struct RpState {
    std::vector<std::int64_t> mOf;                  // per face, -1 outside every patch
    std::vector<std::vector<std::int64_t>> mFaces;  // per patch, ascending
    std::vector<RpProxy> mProxy;
    std::vector<double> mEnergy;
    std::vector<std::int64_t> mIterations;  // the iterations that changed the patch
};

// Refit the patches flagged dirty, in parallel, each serially.
void rp_refit(const RpSurface& rS, RpState& rState, const std::vector<std::uint8_t>& rDirty,
              std::uint32_t Types) {
    std::vector<std::size_t> todo;
    for (std::size_t r = 0; r < rState.mFaces.size(); ++r)
        if (rDirty[r] && !rState.mFaces[r].empty())
            todo.push_back(r);
    parallel_for(
        todo.size(),
        [&](std::size_t i) {
            const std::size_t r = todo[i];
            const RpProxy keep = rState.mProxy[r];
            rState.mProxy[r] = rp_fit(rS, rState.mFaces[r], Types, &keep);
            rState.mEnergy[r] = rp_energy(rS, rState.mFaces[r], rState.mProxy[r]);
        },
        1);
}

double rp_total_energy(const RpState& rState) {
    double e = 0.0;
    for (std::size_t r = 0; r < rState.mFaces.size(); ++r)
        if (!rState.mFaces[r].empty())
            e += rState.mEnergy[r];
    return e;
}

double rp_explained(const RpSurface& rS, const RpState& rState, double TotalArea) {
    double a = 0.0;
    for (std::size_t r = 0; r < rState.mFaces.size(); ++r)
        if (rp_is_primitive(rState.mProxy[r]))
            for (std::int64_t f : rState.mFaces[r])
                a += rS.mArea[static_cast<std::size_t>(f)];
    return TotalArea > 0.0 ? a / TotalArea : 0.0;
}

// Rebuild the per-patch face lists from the per-face assignment.
void rp_rebuild(RpState& rState, std::size_t NumPatches) {
    rState.mFaces.assign(NumPatches, {});
    for (std::size_t f = 0; f < rState.mOf.size(); ++f)
        if (rState.mOf[f] >= 0)
            rState.mFaces[static_cast<std::size_t>(rState.mOf[f])].push_back(
                static_cast<std::int64_t>(f));
}

// The edge-connected pieces of a set of faces across smooth edges, each
// ascending, in the order of their lowest face.
std::vector<std::vector<std::int64_t>> rp_components(const RpSurface& rS,
                                                     const std::vector<std::int64_t>& rFaces) {
    std::map<std::int64_t, std::size_t> index;
    for (std::size_t i = 0; i < rFaces.size(); ++i)
        index[rFaces[i]] = i;
    std::vector<std::uint8_t> seen(rFaces.size(), 0);
    std::vector<std::vector<std::int64_t>> out;
    for (std::size_t i = 0; i < rFaces.size(); ++i) {
        if (seen[i])
            continue;
        std::vector<std::int64_t> comp;
        std::vector<std::int64_t> stack{rFaces[i]};
        seen[i] = 1;
        while (!stack.empty()) {
            const std::int64_t f = stack.back();
            stack.pop_back();
            comp.push_back(f);
            const std::size_t uf = static_cast<std::size_t>(f);
            for (std::size_t k = rS.mAdjStart[uf]; k < rS.mAdjStart[uf + 1]; ++k) {
                const auto it = index.find(rS.mAdj[k]);
                if (it != index.end() && !seen[it->second]) {
                    seen[it->second] = 1;
                    stack.push_back(rS.mAdj[k]);
                }
            }
        }
        std::sort(comp.begin(), comp.end());
        out.push_back(std::move(comp));
    }
    return out;
}

// One Lloyd step: seeds, then flooding from all of them through one priority
// queue. Returns false when no face moved.
bool rp_lloyd(const RpSurface& rS, RpState& rState) {
    const std::size_t np = rState.mFaces.size();
    using Entry = std::tuple<double, std::int64_t, std::int64_t>;  // (cost, patch, face)
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    std::vector<std::int64_t> next(rState.mOf.size(), -1);
    std::vector<std::uint8_t> fixed(np, 0);
    for (std::size_t r = 0; r < np; ++r) {
        const std::vector<std::int64_t>& faces = rState.mFaces[r];
        if (faces.empty())
            continue;
        if (rState.mProxy[r].mKind == RefineProxy::Facets) {
            // Nothing to grow from: the patch keeps its faces.
            fixed[r] = 1;
            for (std::int64_t f : faces)
                next[static_cast<std::size_t>(f)] = static_cast<std::int64_t>(r);
            continue;
        }
        std::int64_t seed = faces.front();
        double best = std::numeric_limits<double>::infinity();
        for (std::int64_t f : faces) {
            const double c =
                rp_cost(rS, static_cast<std::size_t>(f), rState.mProxy[r]) /
                std::max(rS.mArea[static_cast<std::size_t>(f)], std::numeric_limits<double>::min());
            if (c < best) {
                best = c;
                seed = f;
            }
        }
        queue.emplace(0.0, static_cast<std::int64_t>(r), seed);
    }
    while (!queue.empty()) {
        const auto [cost, r, f] = queue.top();
        queue.pop();
        (void)cost;
        const std::size_t uf = static_cast<std::size_t>(f);
        if (next[uf] >= 0)
            continue;
        next[uf] = r;
        const RpProxy& proxy = rState.mProxy[static_cast<std::size_t>(r)];
        for (std::size_t k = rS.mAdjStart[uf]; k < rS.mAdjStart[uf + 1]; ++k) {
            const std::size_t g = static_cast<std::size_t>(rS.mAdj[k]);
            if (next[g] >= 0 || rState.mOf[g] < 0 ||
                fixed[static_cast<std::size_t>(rState.mOf[g])] || !rp_takes(rS, g, proxy))
                continue;
            queue.emplace(rp_cost(rS, g, proxy), r, static_cast<std::int64_t>(g));
        }
    }
    // A face no seed reached (cannot happen on a partition into smooth-connected
    // patches, but keep it safe) keeps its patch.
    for (std::size_t f = 0; f < next.size(); ++f)
        if (rState.mOf[f] >= 0 && next[f] < 0)
            next[f] = rState.mOf[f];
    const bool moved = next != rState.mOf;
    rState.mOf = std::move(next);
    return moved;
}

}  // namespace

RefinePatchesResult refine_patches(const Mesh& rMesh, const RefinePatchesOptions& rOptions) {
    if (!(rOptions.mTolerance > 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kRpPrefix) + "tolerance must be finite and > 0");
    if (!(rOptions.mNormalWeight >= 0.0) || !std::isfinite(rOptions.mNormalWeight))
        throw std::invalid_argument(std::string(kRpPrefix) +
                                    "normal_weight must be finite and >= 0");
    if (!(rOptions.mEpsilon >= 0.0) || !std::isfinite(rOptions.mEpsilon))
        throw std::invalid_argument(std::string(kRpPrefix) + "epsilon must be finite and >= 0");
    if (!(rOptions.mFeatureAngleDeg >= 0.0 && rOptions.mFeatureAngleDeg <= 180.0))
        throw std::invalid_argument(std::string(kRpPrefix) +
                                    "feature_angle must be in [0, 180] degrees");
    if (rOptions.mMaxIterations < 0)
        throw std::invalid_argument(std::string(kRpPrefix) + "max_iterations must be >= 0");
    if (rOptions.mMinPatchCells < 1)
        throw std::invalid_argument(std::string(kRpPrefix) + "min_patch_cells must be >= 1");
    if (rOptions.mTypes & ~kFitAllTypes)
        throw std::invalid_argument(std::string(kRpPrefix) + "types holds an unknown primitive");
    if (rOptions.mPrefix.empty())
        throw std::invalid_argument(std::string(kRpPrefix) + "prefix must not be empty");
    if (!rMesh.HasCellData(kSegmentPatchName))
        throw std::invalid_argument(std::string(kRpPrefix) + "the mesh has no '" +
                                    kSegmentPatchName + "' cell data; run segment first");

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(bases));

    // Per global cell, its input label.
    std::vector<std::int64_t> cell_label(ncells, -1);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const NDArray& a = rMesh.CellData(kSegmentPatchName, b);
        const std::size_t nc = rMesh.Cells(b).NumCells();
        if (a.Size() != nc)
            throw std::invalid_argument(std::string(kRpPrefix) + "'" + kSegmentPatchName +
                                        "' must hold one label per cell");
        if (a.Dtype() == DType::Float32 || a.Dtype() == DType::Float64)
            throw std::invalid_argument(std::string(kRpPrefix) + "'" + kSegmentPatchName +
                                        "' must be an integer label");
        for (std::size_t c = 0; c < nc; ++c)
            cell_label[static_cast<std::size_t>(bases[b]) + c] = detail::read_int(a, c);
    }

    RpSurface s;
    s.mRings = detail::surface_rings(rMesh, std::vector<std::uint8_t>(ncells, 1), kRpPrefix);
    const std::size_t nf = s.NumFaces();
    const std::vector<std::int64_t>& start = s.mRings.mStart;
    const std::vector<std::int64_t>& nodes = s.mRings.mNodes;
    s.mXyz = detail::surface_xyz(rMesh);
    s.mNormal.assign(nf * 3, 0.0);
    s.mCentroid.assign(nf * 3, 0.0);
    s.mArea.assign(nf, 0.0);
    s.mLabel.assign(nf, -1);
    parallel_for(nf, [&](std::size_t f) {
        const std::int64_t* ring = nodes.data() + start[f];
        const std::size_t size = static_cast<std::size_t>(start[f + 1] - start[f]);
        detail::ring_unit_normal(s.mXyz.data(), ring, size, &s.mNormal[f * 3]);
        double c[3] = {0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < size; ++k)
            for (std::size_t d = 0; d < 3; ++d)
                c[d] += s.mXyz[static_cast<std::size_t>(ring[k]) * 3 + d];
        for (std::size_t d = 0; d < 3; ++d)
            s.mCentroid[f * 3 + d] = c[d] / static_cast<double>(size);
        // The area by the fan from the first corner.
        double ax = 0.0, ay = 0.0, az = 0.0;
        const double* p0 = &s.mXyz[static_cast<std::size_t>(ring[0]) * 3];
        for (std::size_t k = 1; k + 1 < size; ++k) {
            const double* p1 = &s.mXyz[static_cast<std::size_t>(ring[k]) * 3];
            const double* p2 = &s.mXyz[static_cast<std::size_t>(ring[k + 1]) * 3];
            const double u[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const double v[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            ax += u[1] * v[2] - u[2] * v[1];
            ay += u[2] * v[0] - u[0] * v[2];
            az += u[0] * v[1] - u[1] * v[0];
        }
        s.mArea[f] = 0.5 * std::sqrt(ax * ax + ay * ay + az * az);
        s.mLabel[f] = cell_label[static_cast<std::size_t>(s.mRings.mCell[f])];
    });

    // The tolerance and the bounding-box diagonal of the surface's points.
    {
        double lo[3] = {std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
        double hi[3] = {-lo[0], -lo[1], -lo[2]};
        for (std::int64_t v : nodes)
            for (std::size_t d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], s.mXyz[static_cast<std::size_t>(v) * 3 + d]);
                hi[d] = std::max(hi[d], s.mXyz[static_cast<std::size_t>(v) * 3 + d]);
            }
        double d2 = 0.0;
        for (std::size_t d = 0; d < 3; ++d)
            d2 += nodes.empty() ? 0.0 : (hi[d] - lo[d]) * (hi[d] - lo[d]);
        s.mDiag = std::sqrt(d2);
    }
    s.mTol = rOptions.mRelative ? rOptions.mTolerance * s.mDiag : rOptions.mTolerance;
    if (!(s.mTol > 0.0))
        s.mTol = rOptions.mTolerance;
    s.mWeight = rOptions.mNormalWeight;

    // Paired edges and the smooth neighbours (segment's crease rule).
    {
        struct Use {
            std::int64_t mLo, mHi, mFace;
            bool mForward;
        };
        std::vector<Use> uses;
        for (std::size_t f = 0; f < nf; ++f) {
            const std::int64_t b = start[f], n = start[f + 1] - b;
            for (std::int64_t k = 0; k < n; ++k) {
                const std::int64_t u = nodes[static_cast<std::size_t>(b + k)];
                const std::int64_t v = nodes[static_cast<std::size_t>(b + (k + 1) % n)];
                if (u == v)
                    continue;
                const std::int64_t face = static_cast<std::int64_t>(f);
                uses.push_back(u < v ? Use{u, v, face, true} : Use{v, u, face, false});
            }
        }
        parallel_sort(uses.begin(), uses.end(), [](const Use& rA, const Use& rB) {
            return std::tie(rA.mLo, rA.mHi, rA.mFace, rA.mForward) <
                   std::tie(rB.mLo, rB.mHi, rB.mFace, rB.mForward);
        });
        const double cos_feature = std::cos(rOptions.mFeatureAngleDeg * kRpPi / 180.0);
        std::vector<std::pair<std::int64_t, std::int64_t>> smooth;
        for (std::size_t b = 0; b < uses.size();) {
            std::size_t e = b + 1;
            while (e < uses.size() && uses[e].mLo == uses[b].mLo && uses[e].mHi == uses[b].mHi)
                ++e;
            if (e - b == 2 && uses[b].mFace != uses[b + 1].mFace) {
                const std::int64_t fa = uses[b].mFace, fb = uses[b + 1].mFace;
                s.mPairs.emplace_back(fa, fb);
                const double* na = &s.mNormal[static_cast<std::size_t>(fa) * 3];
                const double* nb = &s.mNormal[static_cast<std::size_t>(fb) * 3];
                const bool za = na[0] == 0.0 && na[1] == 0.0 && na[2] == 0.0;
                const bool zb = nb[0] == 0.0 && nb[1] == 0.0 && nb[2] == 0.0;
                bool ok = true;
                if (!za && !zb) {
                    double dot = rp_dot(na, nb);
                    if (uses[b].mForward == uses[b + 1].mForward)
                        dot = -dot;
                    ok = !(dot < cos_feature);
                }
                const std::int64_t la = s.mLabel[static_cast<std::size_t>(fa)];
                const std::int64_t lb = s.mLabel[static_cast<std::size_t>(fb)];
                if (ok && la >= 0 && lb >= 0)
                    smooth.emplace_back(fa, fb);
            }
            b = e;
        }
        s.mAdjStart.assign(nf + 1, 0);
        for (const auto& p : smooth) {
            ++s.mAdjStart[static_cast<std::size_t>(p.first) + 1];
            ++s.mAdjStart[static_cast<std::size_t>(p.second) + 1];
        }
        for (std::size_t f = 0; f < nf; ++f)
            s.mAdjStart[f + 1] += s.mAdjStart[f];
        s.mAdj.assign(s.mAdjStart[nf], 0);
        std::vector<std::size_t> fill(s.mAdjStart.begin(), s.mAdjStart.end() - 1);
        for (const auto& p : smooth) {
            s.mAdj[fill[static_cast<std::size_t>(p.first)]++] = p.second;
            s.mAdj[fill[static_cast<std::size_t>(p.second)]++] = p.first;
        }
        for (std::size_t f = 0; f < nf; ++f)
            std::sort(s.mAdj.begin() + static_cast<std::ptrdiff_t>(s.mAdjStart[f]),
                      s.mAdj.begin() + static_cast<std::ptrdiff_t>(s.mAdjStart[f + 1]));
    }

    // The initial partition: each input label split into its smooth-connected
    // pieces (a label is never joined across a crease).
    RpState st;
    st.mOf.assign(nf, -1);
    double total_area = 0.0;
    {
        std::map<std::int64_t, std::vector<std::int64_t>> by_label;
        for (std::size_t f = 0; f < nf; ++f)
            if (s.mLabel[f] >= 0) {
                by_label[s.mLabel[f]].push_back(static_cast<std::int64_t>(f));
                total_area += s.mArea[f];
            }
        std::vector<std::vector<std::int64_t>> pieces;
        for (const auto& kv : by_label)
            for (std::vector<std::int64_t>& c : rp_components(s, kv.second))
                pieces.push_back(std::move(c));
        // Patches numbered in the order of their lowest face.
        std::sort(pieces.begin(), pieces.end(),
                  [](const auto& rA, const auto& rB) { return rA.front() < rB.front(); });
        for (std::size_t r = 0; r < pieces.size(); ++r)
            for (std::int64_t f : pieces[r])
                st.mOf[static_cast<std::size_t>(f)] = static_cast<std::int64_t>(r);
        st.mFaces = std::move(pieces);
    }
    st.mProxy.assign(st.mFaces.size(), RpProxy{});
    st.mEnergy.assign(st.mFaces.size(), 0.0);
    st.mIterations.assign(st.mFaces.size(), 0);
    rp_refit(s, st, std::vector<std::uint8_t>(st.mFaces.size(), 1), rOptions.mTypes);

    RefinePatchesResult result;
    result.mTolerance = s.mTol;
    result.mInitialEnergy = rp_total_energy(st);
    result.mInitialExplainedAreaFraction = rp_explained(s, st, total_area);
    double explained = result.mInitialExplainedAreaFraction;
    result.mStop = rOptions.mMaxIterations == 0 ? RefineStop::IterationCap : RefineStop::Unchanged;
    const std::size_t min_cells = static_cast<std::size_t>(rOptions.mMinPatchCells);

    for (std::int64_t it = 1; it <= rOptions.mMaxIterations; ++it) {
        const std::vector<std::int64_t> before_of = st.mOf;
        const std::vector<std::vector<std::int64_t>> before_faces = st.mFaces;

        const std::int64_t moves_before =
            result.mNumSplits + result.mNumMerges + result.mNumAbsorbed;
        // Lloyd step, kept only when the energy does not rise.
        bool lloyd_undone = false;
        {
            RpState saved = st;
            const double e0 = rp_total_energy(st);
            if (rp_lloyd(s, st)) {
                rp_rebuild(st, st.mFaces.size());
                // A cell no patch could take kept its old patch, which may now
                // be in pieces: each further piece is a patch of its own.
                for (std::size_t r = 0, np = st.mFaces.size(); r < np; ++r) {
                    if (st.mFaces[r].size() < 2)
                        continue;
                    std::vector<std::vector<std::int64_t>> comps = rp_components(s, st.mFaces[r]);
                    for (std::size_t c = 1; c < comps.size(); ++c) {
                        const std::int64_t id = static_cast<std::int64_t>(st.mFaces.size());
                        for (std::int64_t f : comps[c])
                            st.mOf[static_cast<std::size_t>(f)] = id;
                        st.mFaces.push_back(std::move(comps[c]));
                        st.mProxy.emplace_back();
                        st.mEnergy.push_back(0.0);
                        st.mIterations.push_back(st.mIterations[r]);
                    }
                    if (comps.size() > 1)
                        st.mFaces[r] = std::move(comps[0]);
                }
                std::vector<std::uint8_t> dirty(st.mFaces.size(), 0);
                for (std::size_t r = 0; r < st.mFaces.size(); ++r)
                    dirty[r] = r >= saved.mFaces.size() || st.mFaces[r] != saved.mFaces[r];
                rp_refit(s, st, dirty, rOptions.mTypes);
                const double e1 = rp_total_energy(st);
                if (e1 > e0 + kRpEnergySlack * (e0 + total_area)) {
                    // Undone; the moves below may still change the partition.
                    st = std::move(saved);
                    lloyd_undone = true;
                } else {
                    ++result.mNumIterations;
                }
            }
        }

        // Split: a primitive found inside a freeform patch by seeded RANSAC.
        if (rOptions.mAllowSplit) {
            const std::size_t np = st.mFaces.size();
            for (std::size_t r = 0; r < np; ++r) {
                const std::vector<std::int64_t> faces = st.mFaces[r];
                if (st.mProxy[r].mKind != RefineProxy::Freeform ||
                    faces.size() < std::max<std::size_t>(2 * min_cells, 6))
                    continue;
                std::map<std::int64_t, std::size_t> in;
                for (std::size_t i = 0; i < faces.size(); ++i)
                    in[faces[i]] = i;
                RpRng rng{rOptions.mRansacSeed ^ (0x9e3779b97f4a7c15ull * (r + 1)) ^
                          (static_cast<std::uint64_t>(it) << 32)};
                std::vector<std::int64_t> best_inliers;
                RpProxy best_proxy;
                double best_area = 0.0;
                std::vector<double> pts, nrm, ctr;
                for (int h = 0; h < kRpHypotheses; ++h) {
                    const std::int64_t seed = faces[rng.Next() % faces.size()];
                    // The seed's neighbourhood breadth first, ring by ring: its
                    // 1-ring, its 2-ring and up to kRpNeighbourhood cells are
                    // three hypotheses (a plane needs few cells, a curved
                    // surface more, and a large one may span two surfaces).
                    std::vector<std::int64_t> hood{seed};
                    std::vector<std::size_t> ring_end{1};
                    std::vector<std::uint8_t> seen(faces.size(), 0);
                    seen[in[seed]] = 1;
                    for (std::size_t q = 0; q < hood.size() && hood.size() < kRpNeighbourhood;
                         ++q) {
                        if (q == ring_end.back())
                            ring_end.push_back(hood.size());
                        const std::size_t uf = static_cast<std::size_t>(hood[q]);
                        for (std::size_t k = s.mAdjStart[uf];
                             k < s.mAdjStart[uf + 1] && hood.size() < kRpNeighbourhood; ++k) {
                            const auto at = in.find(s.mAdj[k]);
                            if (at != in.end() && !seen[at->second]) {
                                seen[at->second] = 1;
                                hood.push_back(s.mAdj[k]);
                            }
                        }
                    }
                    ring_end.push_back(hood.size());
                    // ring_end[i] is the size of the hood through ring i.
                    std::vector<std::size_t> sizes;
                    for (std::size_t i : {std::size_t{1}, std::size_t{2}})
                        if (i < ring_end.size())
                            sizes.push_back(ring_end[i]);
                    sizes.push_back(std::min<std::size_t>(hood.size(), 24));
                    sizes.push_back(hood.size());
                    std::sort(sizes.begin(), sizes.end());
                    sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());
                    for (std::size_t size : sizes) {
                        std::vector<std::int64_t> sub(
                            hood.begin(), hood.begin() + static_cast<std::ptrdiff_t>(size));
                        std::sort(sub.begin(), sub.end());
                        rp_samples(s, sub, pts, nrm, ctr);
                        detail::PrimitiveFit fit;
                        if (!rp_fit_primitive(s, pts, nrm, ctr, rOptions.mTypes, fit))
                            continue;
                        // The connected cells it explains, grown from the seed.
                        const auto explains = [&](std::int64_t F) {
                            for (std::int64_t k = start[static_cast<std::size_t>(F)];
                                 k < start[static_cast<std::size_t>(F) + 1]; ++k) {
                                const double* p = &s.mXyz[static_cast<std::size_t>(
                                                              nodes[static_cast<std::size_t>(k)]) *
                                                          3];
                                if (!(detail::primitive_point_distance(fit.mType, fit.mParams, p) <=
                                      s.mTol))
                                    return false;
                            }
                            return true;
                        };
                        if (!explains(seed))
                            continue;
                        std::vector<std::int64_t> inl{seed};
                        std::vector<std::uint8_t> took(faces.size(), 0);
                        took[in[seed]] = 1;
                        double area = 0.0;
                        for (std::size_t q = 0; q < inl.size(); ++q) {
                            const std::size_t uf = static_cast<std::size_t>(inl[q]);
                            area += s.mArea[uf];
                            for (std::size_t k = s.mAdjStart[uf]; k < s.mAdjStart[uf + 1]; ++k) {
                                const auto at = in.find(s.mAdj[k]);
                                if (at != in.end() && !took[at->second]) {
                                    took[at->second] = 1;
                                    if (explains(s.mAdj[k]))
                                        inl.push_back(s.mAdj[k]);
                                }
                            }
                        }
                        if (area > best_area) {
                            best_area = area;
                            std::sort(inl.begin(), inl.end());
                            best_inliers = std::move(inl);
                            best_proxy = RpProxy{};
                            best_proxy.mKind = static_cast<RefineProxy>(fit.mType);
                            best_proxy.mFit = fit;
                        }
                    }
                }
                if (best_inliers.size() < std::max<std::size_t>(min_cells, 3) ||
                    best_inliers.size() >= faces.size())
                    continue;
                // The pieces: the inliers, and the rest's connected parts.
                std::vector<std::int64_t> rest;
                std::set_difference(faces.begin(), faces.end(), best_inliers.begin(),
                                    best_inliers.end(), std::back_inserter(rest));
                std::vector<std::vector<std::int64_t>> pieces{best_inliers};
                for (std::vector<std::int64_t>& c : rp_components(s, rest))
                    pieces.push_back(std::move(c));
                std::vector<RpProxy> proxies(pieces.size());
                std::vector<double> energies(pieces.size());
                for (std::size_t i = 0; i < pieces.size(); ++i) {
                    proxies[i] =
                        rp_fit(s, pieces[i], rOptions.mTypes, i == 0 ? &best_proxy : nullptr);
                    energies[i] = rp_energy(s, pieces[i], proxies[i]);
                }
                // The cells the primitive takes must cost less under it than
                // under the freeform proxy they leave.
                if (!rp_is_primitive(proxies[0]) ||
                    !(energies[0] <= rp_energy(s, pieces[0], st.mProxy[r])))
                    continue;
                // Keep the piece holding the patch's lowest face under its id.
                std::size_t keep = 0;
                for (std::size_t i = 1; i < pieces.size(); ++i)
                    if (pieces[i].front() < pieces[keep].front())
                        keep = i;
                for (std::size_t i = 0; i < pieces.size(); ++i) {
                    const std::size_t id = i == keep ? r : st.mFaces.size();
                    if (id == st.mFaces.size()) {
                        st.mFaces.emplace_back();
                        st.mProxy.emplace_back();
                        st.mEnergy.push_back(0.0);
                        st.mIterations.push_back(st.mIterations[r]);
                    }
                    for (std::int64_t f : pieces[i])
                        st.mOf[static_cast<std::size_t>(f)] = static_cast<std::int64_t>(id);
                    st.mFaces[id] = pieces[i];
                    st.mProxy[id] = proxies[i];
                    st.mEnergy[id] = energies[i];
                }
                ++result.mNumSplits;
            }
        }

        // Merge: neighbours one primitive explains together, cheapest first.
        if (rOptions.mAllowMerge) {
            std::map<std::pair<std::int64_t, std::int64_t>, int> pairs;
            for (std::size_t f = 0; f < nf; ++f) {
                const std::int64_t a = st.mOf[f];
                if (a < 0)
                    continue;
                for (std::size_t k = s.mAdjStart[f]; k < s.mAdjStart[f + 1]; ++k) {
                    const std::int64_t b = st.mOf[static_cast<std::size_t>(s.mAdj[k])];
                    if (b > a)
                        pairs[{a, b}] = 1;
                }
            }
            std::vector<std::pair<std::int64_t, std::int64_t>> cand;
            for (const auto& kv : pairs) {
                const RpProxy& pa = st.mProxy[static_cast<std::size_t>(kv.first.first)];
                const RpProxy& pb = st.mProxy[static_cast<std::size_t>(kv.first.second)];
                if (pa.mKind != RefineProxy::Facets && pb.mKind != RefineProxy::Facets)
                    cand.push_back(kv.first);
            }
            std::vector<RpProxy> unions(cand.size());
            std::vector<double> union_energy(cand.size(), 0.0);
            parallel_for(
                cand.size(),
                [&](std::size_t i) {
                    std::vector<std::int64_t> u =
                        st.mFaces[static_cast<std::size_t>(cand[i].first)];
                    const std::vector<std::int64_t>& v =
                        st.mFaces[static_cast<std::size_t>(cand[i].second)];
                    u.insert(u.end(), v.begin(), v.end());
                    std::sort(u.begin(), u.end());
                    unions[i] = rp_fit(s, u, rOptions.mTypes);
                    union_energy[i] = rp_energy(s, u, unions[i]);
                },
                1);
            std::vector<std::size_t> order;
            for (std::size_t i = 0; i < cand.size(); ++i)
                if (rp_is_primitive(unions[i]))
                    order.push_back(i);
            std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) {
                return unions[x].mRms < unions[y].mRms;
            });
            std::vector<std::uint8_t> used(st.mFaces.size(), 0);
            for (std::size_t i : order) {
                const std::size_t a = static_cast<std::size_t>(cand[i].first);
                const std::size_t b = static_cast<std::size_t>(cand[i].second);
                if (used[a] || used[b])
                    continue;
                const double sum = st.mEnergy[a] + st.mEnergy[b];
                if (!(union_energy[i] <= sum * (1.0 + 1e-9) + 1e-12 * s.mTol * s.mTol))
                    continue;
                used[a] = used[b] = 1;
                std::vector<std::int64_t>& fa = st.mFaces[a];
                fa.insert(fa.end(), st.mFaces[b].begin(), st.mFaces[b].end());
                std::sort(fa.begin(), fa.end());
                for (std::int64_t f : st.mFaces[b])
                    st.mOf[static_cast<std::size_t>(f)] = static_cast<std::int64_t>(a);
                st.mFaces[b].clear();
                st.mProxy[a] = unions[i];
                st.mEnergy[a] = union_energy[i];
                st.mIterations[a] = std::max(st.mIterations[a], st.mIterations[b]);
                ++result.mNumMerges;
            }
        }

        // Absorb: a small patch joins the smooth neighbour it shares most edges
        // with -- a fragment the flooding left, whatever proxy two triangles fit.
        if (min_cells > 1) {
            for (std::size_t r = 0; r < st.mFaces.size(); ++r) {
                if (st.mFaces[r].empty() || st.mFaces[r].size() >= min_cells)
                    continue;
                std::map<std::int64_t, std::int64_t> shared;
                for (std::int64_t f : st.mFaces[r]) {
                    const std::size_t uf = static_cast<std::size_t>(f);
                    for (std::size_t k = s.mAdjStart[uf]; k < s.mAdjStart[uf + 1]; ++k) {
                        const std::int64_t b = st.mOf[static_cast<std::size_t>(s.mAdj[k])];
                        if (b >= 0 && b != static_cast<std::int64_t>(r))
                            ++shared[b];
                    }
                }
                if (shared.empty())
                    continue;
                std::int64_t to = shared.begin()->first, most = shared.begin()->second;
                for (const auto& kv : shared)
                    if (kv.second > most) {
                        to = kv.first;
                        most = kv.second;
                    }
                std::vector<std::int64_t>& ft = st.mFaces[static_cast<std::size_t>(to)];
                ft.insert(ft.end(), st.mFaces[r].begin(), st.mFaces[r].end());
                std::sort(ft.begin(), ft.end());
                for (std::int64_t f : st.mFaces[r])
                    st.mOf[static_cast<std::size_t>(f)] = to;
                st.mFaces[r].clear();
                std::vector<std::uint8_t> dirty(st.mFaces.size(), 0);
                dirty[static_cast<std::size_t>(to)] = 1;
                rp_refit(s, st, dirty, rOptions.mTypes);
                ++result.mNumAbsorbed;
            }
        }

        // Count the iteration on every patch it changed.
        rp_rebuild(st, st.mFaces.size());
        for (std::size_t r = 0; r < st.mFaces.size(); ++r)
            if (!st.mFaces[r].empty() &&
                (r >= before_faces.size() || st.mFaces[r] != before_faces[r]))
                ++st.mIterations[r];
        const double now = rp_explained(s, st, total_area);
        const double gain = now - explained;
        explained = now;
        if (st.mOf == before_of) {
            result.mStop = lloyd_undone ? RefineStop::EnergyRose : RefineStop::Unchanged;
            break;
        }
        const bool moved =
            result.mNumSplits + result.mNumMerges + result.mNumAbsorbed > moves_before;
        if (gain < rOptions.mEpsilon && !moved) {
            result.mStop = RefineStop::SmallGain;
            break;
        }
        if (it == rOptions.mMaxIterations)
            result.mStop = RefineStop::IterationCap;
    }

    // Final labels, in the order of each patch's lowest face.
    std::vector<std::int64_t> final_id(st.mFaces.size(), -1);
    std::vector<std::size_t> order;  // final id -> patch
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t r = st.mOf[f];
        if (r >= 0 && final_id[static_cast<std::size_t>(r)] < 0) {
            final_id[static_cast<std::size_t>(r)] = static_cast<std::int64_t>(order.size());
            order.push_back(static_cast<std::size_t>(r));
        }
    }
    const std::size_t m = order.size();
    std::vector<std::int64_t> cell_patch(ncells, -1);
    for (std::size_t f = 0; f < nf; ++f)
        if (st.mOf[f] >= 0)
            cell_patch[static_cast<std::size_t>(s.mRings.mCell[f])] =
                final_id[static_cast<std::size_t>(st.mOf[f])];

    // Per patch its errors, over its own points; and the global ones.
    std::vector<double> rms(m), mx(m);
    double h_max = 0.0, h_sum = 0.0;
    std::size_t h_n = 0;
    for (std::size_t i = 0; i < m; ++i) {
        const std::size_t r = order[i];
        std::vector<double> pts, nrm, ctr;
        rp_samples(s, st.mFaces[r], pts, nrm, ctr);
        rp_errors(st.mProxy[r], pts, rms[i], mx[i]);
        for (std::size_t k = 0; k < pts.size() / 3; ++k) {
            const double d = rp_distance(st.mProxy[r], &pts[k * 3], nullptr, nullptr);
            h_max = std::max(h_max, d);
            h_sum += d * d;
            ++h_n;
        }
        const RefineProxy kind = st.mProxy[r].mKind;
        if (rp_is_primitive(st.mProxy[r]))
            ++result.mNumPrimitive;
        else if (kind == RefineProxy::Freeform)
            ++result.mNumFreeform;
        else
            ++result.mNumFacets;
    }
    result.mNumPatches = static_cast<std::int64_t>(m);
    result.mHausdorffMax = h_max;
    result.mHausdorffRms = h_n ? std::sqrt(h_sum / static_cast<double>(h_n)) : 0.0;
    result.mExplainedAreaFraction = rp_explained(s, st, total_area);
    result.mEnergy = rp_total_energy(st);

    // The mesh: the input with the refined labels, adjacency and regions.
    Mesh out = detail::clone_mesh(
        rMesh, [](DataLocation, const std::string&, std::string&) { return true; });
    std::vector<NDArray> blocks;
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const std::size_t base = static_cast<std::size_t>(bases[bi]);
        const std::size_t nc = rMesh.Cells(bi).NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {nc});
        std::copy(cell_patch.begin() + static_cast<std::ptrdiff_t>(base),
                  cell_patch.begin() + static_cast<std::ptrdiff_t>(base + nc),
                  a.As<std::int64_t>());
        blocks.push_back(std::move(a));
    }
    out.AddCellData(kSegmentPatchName, std::move(blocks));
    // A cell that is not a surface face is -1, as segment writes it.
    std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t> adjacency;
    for (const auto& p : s.mPairs) {
        std::int64_t a = st.mOf[static_cast<std::size_t>(p.first)];
        std::int64_t b = st.mOf[static_cast<std::size_t>(p.second)];
        if (a < 0 || b < 0 || a == b)
            continue;
        a = final_id[static_cast<std::size_t>(a)];
        b = final_id[static_cast<std::size_t>(b)];
        if (a > b)
            std::swap(a, b);
        ++adjacency[{a, b}];
    }
    NDArray table = NDArray::Uninit(DType::Int64, {adjacency.size(), std::size_t{3}});
    {
        std::int64_t* t = table.As<std::int64_t>();
        std::size_t r = 0;
        for (const auto& kv : adjacency) {
            t[r * 3 + 0] = kv.first.first;
            t[r * 3 + 1] = kv.first.second;
            t[r * 3 + 2] = kv.second;
            ++r;
        }
    }
    out.AddFieldData(kSegmentAdjacencyName, std::move(table));
    const std::string head = rOptions.mPrefix + ":";
    for (std::size_t i = out.NumRegions(); i-- > 0;) {
        const meshioplusplus::Region& r = out.Region(i);
        if (r.mKind != RegionKind::Cell || r.mName.size() <= head.size() ||
            r.mName.compare(0, head.size(), head) != 0)
            continue;
        const bool digits = std::all_of(r.mName.begin() + static_cast<std::ptrdiff_t>(head.size()),
                                        r.mName.end(), [](char c) { return c >= '0' && c <= '9'; });
        if (digits)
            out.RemoveRegion(i);
    }
    std::vector<std::vector<std::int64_t>> members(m);
    for (std::size_t c = 0; c < ncells; ++c)
        if (cell_patch[c] >= 0)
            members[static_cast<std::size_t>(cell_patch[c])].push_back(
                static_cast<std::int64_t>(c));
    for (std::size_t p = 0; p < m; ++p) {
        NDArray entries = NDArray::Uninit(DType::Int64, {members[p].size()});
        std::copy(members[p].begin(), members[p].end(), entries.As<std::int64_t>());
        out.AddRegion(meshioplusplus::Region(head + std::to_string(p), RegionKind::Cell, 2, -1,
                                             std::move(entries)));
    }

    // The fit table fit_cad reads, at the absolute tolerance.
    FitPrimitivesOptions fo;
    fo.mPatchData = kSegmentPatchName;
    fo.mTolerance = s.mTol;
    fo.mRelative = false;
    fo.mTypes = rOptions.mTypes;
    Mesh fitted = fit_primitives(out, fo).mMesh;

    NDArray pid = NDArray::Uninit(DType::Int64, {m});
    NDArray proxy = NDArray::Uninit(DType::Int64, {m});
    NDArray iters = NDArray::Uninit(DType::Int64, {m});
    NDArray prms = NDArray::Uninit(DType::Float64, {m});
    NDArray pmax = NDArray::Uninit(DType::Float64, {m});
    NDArray cells = NDArray::Uninit(DType::Int64, {m});
    NDArray demotion = NDArray::Uninit(DType::Int64, {m});
    for (std::size_t i = 0; i < m; ++i) {
        const std::size_t r = order[i];
        const RpProxy& p = st.mProxy[r];
        pid.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
        proxy.As<std::int64_t>()[i] = static_cast<std::int64_t>(p.mKind);
        iters.As<std::int64_t>()[i] = st.mIterations[r];
        prms.As<double>()[i] = rms[i];
        pmax.As<double>()[i] = mx[i];
        cells.As<std::int64_t>()[i] = static_cast<std::int64_t>(st.mFaces[r].size());
        RefineDemotion why = RefineDemotion::None;
        if (p.mKind == RefineProxy::Freeform)
            why = RefineDemotion::NoPrimitiveWithinTolerance;
        else if (p.mKind == RefineProxy::Facets)
            why = RefineDemotion::TooSmall;
        demotion.As<std::int64_t>()[i] = static_cast<std::int64_t>(why);
    }
    fitted.AddFieldData(kRefinePatchIdName, std::move(pid));
    fitted.AddFieldData(kRefinePatchProxyName, std::move(proxy));
    fitted.AddFieldData(kRefinePatchIterationsName, std::move(iters));
    fitted.AddFieldData(kRefinePatchRmsName, std::move(prms));
    fitted.AddFieldData(kRefinePatchMaxName, std::move(pmax));
    fitted.AddFieldData(kRefinePatchCellsName, std::move(cells));
    fitted.AddFieldData(kRefinePatchDemotionName, std::move(demotion));
    NDArray ea = NDArray::Uninit(DType::Float64, {1});
    ea.As<double>()[0] = result.mExplainedAreaFraction;
    fitted.AddFieldData(kRefineExplainedAreaName, std::move(ea));
    NDArray hd = NDArray::Uninit(DType::Float64, {2});
    hd.As<double>()[0] = result.mHausdorffMax;
    hd.As<double>()[1] = result.mHausdorffRms;
    fitted.AddFieldData(kRefineHausdorffName, std::move(hd));
    result.mMesh = std::move(fitted);
    return result;
}

}  // namespace meshioplusplus
