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
// The discrete parameterisation of mesh_cad_surface: see
// detail/mesh_cad_reparam.hpp. Anonymous-namespace helpers are prefixed
// `mcr_` (the amalgamation concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

// Project includes
#include "cdt2d.hpp"
#include "mesh_cad_reparam.hpp"
#include "param_disc.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

using McrV2 = std::array<double, 2>;
using McrV3 = std::array<double, 3>;

constexpr double kMcrPi = 3.14159265358979323846;

McrV3 mcr_sub(const McrV3& rA, const McrV3& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

double mcr_dot(const McrV3& rA, const McrV3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

/// The angle at `rO` of the triangle `rO, rA, rB`, in radians; 0 when an
/// edge has no length.
double mcr_angle(const McrV3& rO, const McrV3& rA, const McrV3& rB) {
    const McrV3 a = mcr_sub(rA, rO);
    const McrV3 b = mcr_sub(rB, rO);
    const McrV3 c{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    const double s = std::sqrt(mcr_dot(c, c));
    const double d = mcr_dot(a, b);
    if (!(s > 0.0) && !(std::fabs(d) > 0.0))
        return 0.0;
    return std::atan2(s, d);
}

/// The eigenvalues of `(J L^-1)^T (J L^-1)`, `J = (Su, Sv)` and `L` the
/// upper Cholesky factor `(L11, L12; 0, L22)`: {smallest, largest}.
std::array<double, 2> mcr_eigen(const McrV3& rSu, const McrV3& rSv, double L11, double L12,
                                double L22) {
    McrV3 c1{};
    McrV3 c2{};
    for (std::size_t i = 0; i < 3; ++i) {
        c1[i] = rSu[i] / L11;
        c2[i] = rSv[i] / L22 - rSu[i] * L12 / (L11 * L22);
    }
    const double a = mcr_dot(c1, c1);
    const double b = mcr_dot(c1, c2);
    const double d = mcr_dot(c2, c2);
    const double root = std::sqrt((a - d) * (a - d) + 4.0 * b * b);
    return {0.5 * (a + d - root), 0.5 * (a + d + root)};
}

/// `J D`, the surface's derivative along the disc from its derivative along
/// the parameters and `D = d(u, v)/d(p, q)`.
std::array<McrV3, 2> mcr_chain(const McrV3& rSu, const McrV3& rSv,
                               const std::array<double, 4>& rD) {
    std::array<McrV3, 2> out{};
    for (std::size_t i = 0; i < 3; ++i) {
        out[0][i] = rSu[i] * rD[0] + rSv[i] * rD[2];
        out[1][i] = rSu[i] * rD[1] + rSv[i] * rD[3];
    }
    return out;
}

/// The derivative of the parameters along the disc on triangle `T` at
/// barycentrics `rL`: exact on an ordinary triangle, a central difference of
/// the collapsed coordinates on one at a pole.
std::array<double, 4> mcr_d(const McadDisc& rD, std::size_t T, const std::array<double, 3>& rL) {
    const McadCoarseMesh& r_c = *rD.pCoarse;
    const auto& r_t = r_c.mTris[T];
    const McrV2& p0 = rD.mP[r_t[0]];
    const McrV2& p1 = rD.mP[r_t[1]];
    const McrV2& p2 = rD.mP[r_t[2]];
    const double a11 = p1[0] - p0[0], a12 = p2[0] - p0[0];
    const double a21 = p1[1] - p0[1], a22 = p2[1] - p0[1];
    const double det = a11 * a22 - a12 * a21;
    if (!(std::fabs(det) > 0.0))
        return {0.0, 0.0, 0.0, 0.0};
    // d(lambda_1, lambda_2)/d(p, q): the inverse of A.
    const double i11 = a22 / det, i12 = -a12 / det;
    const double i21 = -a21 / det, i22 = a11 / det;
    if (r_c.mPoleCorner[T] < 0) {
        const auto& r_uv = r_c.mCornerUv[T];
        const double b11 = r_uv[1][0] - r_uv[0][0], b12 = r_uv[2][0] - r_uv[0][0];
        const double b21 = r_uv[1][1] - r_uv[0][1], b22 = r_uv[2][1] - r_uv[0][1];
        return {b11 * i11 + b12 * i21, b11 * i12 + b12 * i22, b21 * i11 + b22 * i21,
                b21 * i12 + b22 * i22};
    }
    // Along lambda_1 and lambda_2 (lambda_0 taking up the rest), a step of
    // a millionth, then through the inverse of A.
    constexpr double kStep = 1e-6;
    const auto at = [&](double D1, double D2) {
        return rD.UvAt(T, {rL[0] - D1 - D2, rL[1] + D1, rL[2] + D2});
    };
    const McrV2 u1p = at(kStep, 0.0), u1m = at(-kStep, 0.0);
    const McrV2 u2p = at(0.0, kStep), u2m = at(0.0, -kStep);
    const double b11 = (u1p[0] - u1m[0]) / (2.0 * kStep), b12 = (u2p[0] - u2m[0]) / (2.0 * kStep);
    const double b21 = (u1p[1] - u1m[1]) / (2.0 * kStep), b22 = (u2p[1] - u2m[1]) / (2.0 * kStep);
    return {b11 * i11 + b12 * i21, b11 * i12 + b12 * i22, b21 * i11 + b22 * i21,
            b21 * i12 + b22 * i22};
}

}  // namespace

McadDisc::V2 McadDisc::UvAt(std::size_t T, const std::array<double, 3>& rL) const {
    const McadCoarseMesh& r_c = *pCoarse;
    const auto& r_uv = r_c.mCornerUv[T];
    const int c = r_c.mPoleCorner[T];
    if (c < 0)
        return {rL[0] * r_uv[0][0] + rL[1] * r_uv[1][0] + rL[2] * r_uv[2][0],
                rL[0] * r_uv[0][1] + rL[1] * r_uv[1][1] + rL[2] * r_uv[2][1]};
    // Collapsed coordinates (Duffy 1982): the pole's corner at the point of
    // its degenerate segment nearest the opposite side's point.
    const std::size_t i = static_cast<std::size_t>(c);
    const std::size_t j = (i + 1) % 3;
    const std::size_t k = (i + 2) % 3;
    const auto& r_seg = r_c.mPoleSegment.at({r_c.mOwner[T], r_c.mTris[T][i]});
    const double s = rL[j] + rL[k];
    V2 b = r_seg[0];
    if (s > 0.0)
        b = {(rL[j] * r_uv[j][0] + rL[k] * r_uv[k][0]) / s,
             (rL[j] * r_uv[j][1] + rL[k] * r_uv[k][1]) / s};
    const double dx = r_seg[1][0] - r_seg[0][0];
    const double dy = r_seg[1][1] - r_seg[0][1];
    const double len2 = dx * dx + dy * dy;
    V2 pole = r_seg[0];
    if (len2 > 0.0) {
        const double t =
            std::clamp(((b[0] - r_seg[0][0]) * dx + (b[1] - r_seg[0][1]) * dy) / len2, 0.0, 1.0);
        pole = {r_seg[0][0] + t * dx, r_seg[0][1] + t * dy};
    }
    return {rL[i] * pole[0] + rL[j] * r_uv[j][0] + rL[k] * r_uv[k][0],
            rL[i] * pole[1] + rL[j] * r_uv[j][1] + rL[k] * r_uv[k][1]};
}

std::array<double, 3> McadDisc::Bary(std::size_t T, const V2& rP) const {
    const auto& r_t = pCoarse->mTris[T];
    const V2& a = mP[r_t[0]];
    const V2& b = mP[r_t[1]];
    const V2& c = mP[r_t[2]];
    const double det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
    if (!(std::fabs(det) > 0.0))
        return {-1.0, -1.0, -1.0};
    const double l1 = ((rP[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (rP[1] - a[1])) / det;
    const double l2 = ((b[0] - a[0]) * (rP[1] - a[1]) - (rP[0] - a[0]) * (b[1] - a[1])) / det;
    return {1.0 - l1 - l2, l1, l2};
}

void McadDisc::Index() {
    const std::size_t nt = pCoarse->mTris.size();
    V2 hi{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    mLo = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    for (const auto& r_t : pCoarse->mTris)
        for (const std::size_t v : r_t)
            for (std::size_t d = 0; d < 2; ++d) {
                mLo[d] = std::min(mLo[d], mP[v][d]);
                hi[d] = std::max(hi[d], mP[v][d]);
            }
    const double w = std::max(hi[0] - mLo[0], 0.0);
    const double hgt = std::max(hi[1] - mLo[1], 0.0);
    const double area = std::max(w * hgt, 1e-300);
    // About two triangles a bucket, at most a million buckets.
    mCell = std::sqrt(2.0 * area / double(std::max<std::size_t>(nt, 1)));
    if (!(mCell > 0.0))
        mCell = std::max({w, hgt, 1.0});
    mNx = std::min<std::size_t>(1000, static_cast<std::size_t>(w / mCell) + 1);
    mNy = std::min<std::size_t>(1000, static_cast<std::size_t>(hgt / mCell) + 1);
    mCell = std::max({mCell, w / double(mNx), hgt / double(mNy)});
    mBuckets.assign(mNx * mNy, {});
    const auto cell_of = [&](double X, std::size_t D, std::size_t N) {
        const double c = std::floor((X - mLo[D]) / mCell);
        return static_cast<std::size_t>(std::clamp(c, 0.0, double(N - 1)));
    };
    for (std::size_t t = 0; t < nt; ++t) {
        const auto& r_t = pCoarse->mTris[t];
        V2 lo{mP[r_t[0]]}, up{mP[r_t[0]]};
        for (const std::size_t v : r_t)
            for (std::size_t d = 0; d < 2; ++d) {
                lo[d] = std::min(lo[d], mP[v][d]);
                up[d] = std::max(up[d], mP[v][d]);
            }
        for (std::size_t j = cell_of(lo[1], 1, mNy); j <= cell_of(up[1], 1, mNy); ++j)
            for (std::size_t i = cell_of(lo[0], 0, mNx); i <= cell_of(up[0], 0, mNx); ++i)
                mBuckets[j * mNx + i].push_back(t);
    }
}

McadDisc::Hit McadDisc::Locate(const V2& rP) const {
    Hit hit;
    const auto cell_of = [&](double X, std::size_t D, std::size_t N) {
        const double c = std::floor((X - mLo[D]) / mCell);
        return static_cast<std::size_t>(std::clamp(c, 0.0, double(N - 1)));
    };
    const std::size_t ci = cell_of(rP[0], 0, mNx);
    const std::size_t cj = cell_of(rP[1], 1, mNy);
    std::size_t best = 0;
    double best_min = -std::numeric_limits<double>::infinity();
    std::array<double, 3> best_l{1.0, 0.0, 0.0};
    const auto consider = [&](const std::vector<std::size_t>& rCands) {
        for (const std::size_t t : rCands) {
            const std::array<double, 3> l = Bary(t, rP);
            const double m = std::min({l[0], l[1], l[2]});
            if (m >= -1e-12) {
                best = t;
                best_l = l;
                hit.mInside = true;
                return;
            }
            if (m > best_min) {
                best_min = m;
                best = t;
                best_l = l;
            }
        }
    };
    consider(mBuckets[cj * mNx + ci]);
    // An empty bucket lies outside every triangle (each is filed in every
    // bucket its box meets), so only the clamped fallback is wanted: the
    // nearest ring of buckets that holds any, searched outward in a fixed
    // order, rather than every triangle, which made each lattice point
    // outside a curved disc cost the whole coarse mesh (found by the fuzz
    // target).
    if (mBuckets[cj * mNx + ci].empty()) {
        const std::size_t rings = std::max(mNx, mNy);
        bool found = false;
        for (std::size_t r = 1; r <= rings && !found; ++r) {
            const auto ring = [&](std::size_t I, std::size_t J) {
                const auto& r_b = mBuckets[J * mNx + I];
                if (!r_b.empty()) {
                    found = true;
                    consider(r_b);
                }
            };
            const std::size_t i0 = ci >= r ? ci - r : 0;
            const std::size_t i1 = std::min(ci + r, mNx - 1);
            const std::size_t j0 = cj >= r ? cj - r : 0;
            const std::size_t j1 = std::min(cj + r, mNy - 1);
            for (std::size_t j = j0; j <= j1; ++j)
                for (std::size_t i = i0; i <= i1; ++i)
                    if (std::max(i > ci ? i - ci : ci - i, j > cj ? j - cj : cj - j) == r)
                        ring(i, j);
        }
        hit.mInside = false;
    }
    double sum = 0.0;
    for (double& r_l : best_l) {
        r_l = std::max(r_l, 0.0);
        sum += r_l;
    }
    if (sum > 0.0)
        for (double& r_l : best_l)
            r_l /= sum;
    else
        best_l = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    hit.mTri = best;
    hit.mOwner = pCoarse->mOwner[best];
    hit.mUv = UvAt(best, best_l);
    hit.mD = mcr_d(*this, best, best_l);
    return hit;
}

double mcad_coarse_condition(const McadCoarseMesh& rCoarse,
                             const std::vector<const BrepSurfaceEval*>& rSurfaces, double L11,
                             double L12, double L22) {
    McadDisc probe;  // only its parameters at a centroid
    probe.pCoarse = &rCoarse;
    double smax = 0.0;
    double smin = std::numeric_limits<double>::infinity();
    for (std::size_t t = 0; t < rCoarse.mTris.size(); ++t) {
        const McrV2 uv = probe.UvAt(t, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
        McrV3 x{}, su{}, sv{};
        rSurfaces[static_cast<std::size_t>(rCoarse.mOwner[t])]->Eval(uv[0], uv[1], x.data(),
                                                                     su.data(), sv.data());
        const auto e = mcr_eigen(su, sv, L11, L12, L22);
        smax = std::max(smax, std::sqrt(std::max(e[1], 0.0)));
        smin = std::min(smin, std::sqrt(std::max(e[0], 0.0)));
    }
    if (!(smax > 0.0))
        return 1.0;
    return smin > 0.0 ? smax / smin : std::numeric_limits<double>::infinity();
}

double mcad_disc_condition(const McadDisc& rDisc) {
    double smax = 0.0;
    double smin = std::numeric_limits<double>::infinity();
    for (const auto& r_j : rDisc.mJac) {
        const auto e = mcr_eigen(r_j[0], r_j[1], rDisc.mL11, rDisc.mL12, rDisc.mL22);
        smax = std::max(smax, std::sqrt(std::max(e[1], 0.0)));
        smin = std::min(smin, std::sqrt(std::max(e[0], 0.0)));
    }
    if (!(smax > 0.0))
        return 1.0;
    return smin > 0.0 ? smax / smin : std::numeric_limits<double>::infinity();
}

std::string mcad_build_disc(const McadCoarseMesh& rCoarse,
                            const std::vector<const BrepSurfaceEval*>& rSurfaces, McadDisc& rDisc) {
    const McadCoarseMesh& r_c = rCoarse;
    const std::size_t nv = r_c.mX.size();
    const std::size_t nt = r_c.mTris.size();
    if (nt == 0)
        return "its coarse mesh has no triangles";
    rDisc = McadDisc{};
    rDisc.pCoarse = &rCoarse;
    rDisc.mSurfaces = rSurfaces;
    // The orientation each triangle must keep in the disc: its owner's.
    const auto want = [&](std::size_t T) {
        return r_c.mReversed[static_cast<std::size_t>(r_c.mOwner[T])] ? -1.0 : 1.0;
    };
    const auto signed_area = [&](const std::vector<double>& rU, const std::vector<double>& rV,
                                 const std::array<std::size_t, 3>& rT) {
        return (rU[rT[1]] - rU[rT[0]]) * (rV[rT[2]] - rV[rT[0]]) -
               (rU[rT[2]] - rU[rT[0]]) * (rV[rT[1]] - rV[rT[0]]);
    };
    // Every real triangle the orientation it must keep, after reflecting v
    // when all are the other way round.
    const auto oriented = [&](std::vector<double>& rU, std::vector<double>& rV) {
        std::size_t pos = 0, neg = 0;
        for (std::size_t t = 0; t < nt; ++t) {
            const double s = signed_area(rU, rV, r_c.mTris[t]) * want(t);
            if (s > 0.0)
                ++pos;
            else if (s < 0.0)
                ++neg;
            else
                return false;
        }
        if (pos > 0 && neg > 0)
            return false;
        if (neg > 0)
            for (double& r_v : rV)
                r_v = -r_v;
        return true;
    };
    bool loops_ok = true;
    const std::vector<std::vector<std::size_t>> loops = param_boundary_loops(r_c.mTris, loops_ok);
    if (!loops_ok || loops.empty())
        return "its coarse mesh's boundary is pinched or not closed";
    // The loops neither cross nor touch in the disc (exact).
    const auto simple = [&](const std::vector<double>& rU, const std::vector<double>& rV) {
        std::vector<Cdt2d::Vec2> pts(nv);
        for (std::size_t i = 0; i < nv; ++i)
            pts[i] = {rU[i], rV[i]};
        std::vector<std::array<std::int64_t, 2>> segs;
        for (const auto& r_l : loops)
            for (std::size_t k = 0; k < r_l.size(); ++k)
                segs.push_back({static_cast<std::int64_t>(r_l[k]),
                                static_cast<std::int64_t>(r_l[(k + 1) % r_l.size()])});
        for (const auto& r_p : pts)
            if (!std::isfinite(r_p[0]) || !std::isfinite(r_p[1]))
                return false;
        return cdt2d_find_crossings(pts, segs).empty();
    };
    std::vector<double> u, v;
    {
        DiscParam lscm = param_lscm(r_c.mX, r_c.mTris);
        if (lscm.mOk && lscm.mNumFlips == 0 && oriented(lscm.mU, lscm.mV) &&
            simple(lscm.mU, lscm.mV)) {
            u = std::move(lscm.mU);
            v = std::move(lscm.mV);
            rDisc.mMethod = 1;
        }
    }
    if (rDisc.mMethod == 0) {
        // Mean-value coordinates on a circle, the holes filled by fans: only
        // where no boundary corner is under 60 degrees on the surface, which
        // the circle would flatten out of refinement's sight.
        std::map<std::size_t, double> corner;
        std::set<std::size_t> on_loop;
        for (const auto& r_l : loops)
            on_loop.insert(r_l.begin(), r_l.end());
        for (const auto& r_t : r_c.mTris)
            for (std::size_t k = 0; k < 3; ++k)
                if (on_loop.count(r_t[k]))
                    corner[r_t[k]] += mcr_angle(r_c.mX[r_t[k]], r_c.mX[r_t[(k + 1) % 3]],
                                                r_c.mX[r_t[(k + 2) % 3]]);
        for (const auto& [vtx, a] : corner)
            if (a < kMcrPi / 3.0)
                return "the conformal map folds it, and a corner under 60 degrees keeps it off "
                       "a circle";
        // The outer loop: the longest on the surface.
        std::size_t outer = 0;
        double longest = -1.0;
        for (std::size_t l = 0; l < loops.size(); ++l) {
            double len = 0.0;
            for (std::size_t k = 0; k < loops[l].size(); ++k) {
                const McrV3 d =
                    mcr_sub(r_c.mX[loops[l][(k + 1) % loops[l].size()]], r_c.mX[loops[l][k]]);
                len += std::sqrt(mcr_dot(d, d));
            }
            if (len > longest) {
                longest = len;
                outer = l;
            }
        }
        std::vector<PatchVec> pts = r_c.mX;
        std::vector<PatchTri> tris = r_c.mTris;
        for (std::size_t l = 0; l < loops.size(); ++l) {
            if (l == outer)
                continue;
            McrV3 centre{0.0, 0.0, 0.0};
            for (const std::size_t w : loops[l])
                for (std::size_t d = 0; d < 3; ++d)
                    centre[d] += r_c.mX[w][d] / double(loops[l].size());
            const std::size_t c = pts.size();
            pts.push_back(centre);
            for (std::size_t k = 0; k < loops[l].size(); ++k)
                tris.push_back({loops[l][(k + 1) % loops[l].size()], loops[l][k], c});
        }
        DiscParam mvc = param_mean_value(pts, tris);
        if (!mvc.mOk)
            return "neither the conformal map nor mean-value coordinates flatten it";
        mvc.mU.resize(nv);
        mvc.mV.resize(nv);
        if (param_flips(r_c.mTris, mvc.mU, mvc.mV) != 0 || !oriented(mvc.mU, mvc.mV) ||
            !simple(mvc.mU, mvc.mV))
            return "neither the conformal map nor mean-value coordinates flatten it without "
                   "folding";
        u = std::move(mvc.mU);
        v = std::move(mvc.mV);
        rDisc.mMethod = 2;
    }
    rDisc.mP.resize(nv);
    for (std::size_t i = 0; i < nv; ++i)
        rDisc.mP[i] = {u[i], v[i]};
    rDisc.Index();
    // The surface's derivative along the disc at each centroid, and the
    // disc's plane: the area-weighted mean of J'^T J'.
    double me = 0.0, mf = 0.0, mg = 0.0, total = 0.0;
    rDisc.mJac.reserve(nt);
    rDisc.mJacAt.reserve(nt);
    for (std::size_t t = 0; t < nt; ++t) {
        const std::array<double, 3> third{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        const McrV2 uv = rDisc.UvAt(t, third);
        McrV3 x{}, su{}, sv{};
        rSurfaces[static_cast<std::size_t>(r_c.mOwner[t])]->Eval(uv[0], uv[1], x.data(), su.data(),
                                                                 sv.data());
        const auto j = mcr_chain(su, sv, mcr_d(rDisc, t, third));
        rDisc.mJac.push_back(j);
        rDisc.mJacAt.push_back(x);
        const auto& r_t = r_c.mTris[t];
        const double area = 0.5 * std::fabs((u[r_t[1]] - u[r_t[0]]) * (v[r_t[2]] - v[r_t[0]]) -
                                            (u[r_t[2]] - u[r_t[0]]) * (v[r_t[1]] - v[r_t[0]]));
        me += area * mcr_dot(j[0], j[0]);
        mf += area * mcr_dot(j[0], j[1]);
        mg += area * mcr_dot(j[1], j[1]);
        total += area;
    }
    if (total > 0.0) {
        me /= total;
        mf /= total;
        mg /= total;
    }
    if (me > 0.0 && mg > 0.0 && me - mf * mf / mg > 1e-12 * me) {
        rDisc.mL11 = std::sqrt(me);
        rDisc.mL12 = mf / rDisc.mL11;
        rDisc.mL22 = std::sqrt(mg - rDisc.mL12 * rDisc.mL12);
    } else {
        return "its flattening has no metric to speak of";
    }
    return {};
}

}  // namespace detail
}  // namespace meshioplusplus
