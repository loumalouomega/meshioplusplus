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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// Fixtures shared by the patch network tests (`test_patch_network.cpp`, the full ones, and
// `test_quick_patch_network.cpp`, their small versions that run in the pull-request coverage
// job): a skin as a patch of points and triangles, and the gap between the patches that meet
// along a shared edge.

#pragma once

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "../../src/cpp/src/detail/bspline_curve.hpp"
#include "../../src/cpp/src/detail/patch_network.hpp"
#include "../../src/cpp/src/detail/quad_layout.hpp"
#include "meshioplusplus/mesh.hpp"

namespace pnfix {

using namespace meshioplusplus;
using detail::PatchTri;
using detail::PatchVec;

struct PnTestPatch {
    std::vector<PatchVec> mPts;
    std::vector<PatchTri> mTris;
};

inline PnTestPatch pn_patch(const Mesh& rMesh, bool TopOnly) {
    PnTestPatch p;
    const double* x = rMesh.Points().As<double>();
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
        p.mPts.push_back({x[i * 3], x[i * 3 + 1], x[i * 3 + 2]});
    const std::int64_t* c = rMesh.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < rMesh.Cells(0).NumCells(); ++t) {
        const PatchTri tri{static_cast<std::size_t>(c[t * 3]),
                           static_cast<std::size_t>(c[t * 3 + 1]),
                           static_cast<std::size_t>(c[t * 3 + 2])};
        if (TopOnly) {
            const PatchVec& a = p.mPts[tri[0]];
            const PatchVec& b = p.mPts[tri[1]];
            const PatchVec& d = p.mPts[tri[2]];
            const double nx = (b[1] - a[1]) * (d[2] - a[2]) - (b[2] - a[2]) * (d[1] - a[1]);
            const double ny = (b[2] - a[2]) * (d[0] - a[0]) - (b[0] - a[0]) * (d[2] - a[2]);
            const double nz = (b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0]);
            if (!(nz > 0.3 * std::sqrt(nx * nx + ny * ny + nz * nz)))
                continue;
        }
        p.mTris.push_back(tri);
    }
    return p;
}

inline double pn_diag(const PnTestPatch& rP) {
    PatchVec lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const PatchTri& t : rP.mTris)
        for (std::size_t v : t)
            for (std::size_t j = 0; j < 3; ++j) {
                lo[j] = std::min(lo[j], rP.mPts[v][j]);
                hi[j] = std::max(hi[j], rP.mPts[v][j]);
            }
    return std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                     (hi[2] - lo[2]) * (hi[2] - lo[2]));
}

inline double pn_dist(const PatchVec& rA, const PatchVec& rB) {
    return std::sqrt((rA[0] - rB[0]) * (rA[0] - rB[0]) + (rA[1] - rB[1]) * (rA[1] - rB[1]) +
                     (rA[2] - rB[2]) * (rA[2] - rB[2]));
}

inline PatchVec pn_eval(const detail::BsplineSurface& rS, double U, double V) {
    detail::BsplineEval e;
    detail::bspline_eval(rS, U, V, 0, e);
    return {e.mS[0], e.mS[1], e.mS[2]};
}

inline void pn_side_uv(std::size_t Side, double S, double& rU, double& rV) {
    const double u[4] = {S, 1.0, 1.0 - S, 0.0}, v[4] = {0.0, S, 1.0, 1.0 - S};
    rU = u[Side];
    rV = v[Side];
}

// The largest gap between the two patches along every shared stretch.
inline double pn_g0_gap(const detail::QuadLayout& rL, const detail::PatchNetwork& rN) {
    double worst = 0.0;
    for (std::size_t c = 0; c < rL.mChains.size(); ++c) {
        struct Use {
            std::size_t mPatch, mSide;
            double mTb, mTe, mS0, mS1;
        };
        std::vector<Use> uses;
        for (std::size_t i = 0; i < rL.mPatches.size(); ++i)
            for (std::size_t j = 0; j < 4; ++j) {
                const auto& segs = rL.mSides[i][j].mSegments;
                if (segs.size() != 1)
                    continue;
                if (segs[0].mChain != c)
                    continue;
                uses.push_back({i, j, rN.mChainParams[c][segs[0].mBegin],
                                rN.mChainParams[c][segs[0].mEnd], 0.0, 1.0});
            }
        for (std::size_t a = 0; a < uses.size(); ++a)
            for (std::size_t b = a + 1; b < uses.size(); ++b) {
                if (rN.mPatches[uses[a].mPatch].mSurface.NumPoles() == 0 ||
                    rN.mPatches[uses[b].mPatch].mSurface.NumPoles() == 0)
                    continue;
                const double lo = std::max(std::min(uses[a].mTb, uses[a].mTe),
                                           std::min(uses[b].mTb, uses[b].mTe));
                const double hi = std::min(std::max(uses[a].mTb, uses[a].mTe),
                                           std::max(uses[b].mTb, uses[b].mTe));
                for (int k = 0; k <= 10 && hi > lo; ++k) {
                    const double t = lo + (hi - lo) * k / 10.0;
                    double ua, va, ub, vb;
                    pn_side_uv(uses[a].mSide, (t - uses[a].mTb) / (uses[a].mTe - uses[a].mTb), ua,
                               va);
                    pn_side_uv(uses[b].mSide, (t - uses[b].mTb) / (uses[b].mTe - uses[b].mTb), ub,
                               vb);
                    worst = std::max(
                        worst, pn_dist(pn_eval(rN.mPatches[uses[a].mPatch].mSurface, ua, va),
                                       pn_eval(rN.mPatches[uses[b].mPatch].mSurface, ub, vb)));
                }
            }
    }
    return worst;
}

}  // namespace pnfix
