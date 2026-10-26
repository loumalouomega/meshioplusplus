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
#include <limits>
#include <vector>

// Project includes
#include "meshioplusplus/nurbs.hpp"
#include "surface_project.hpp"

namespace meshioplusplus::detail {

namespace {

double spj_dot(const double* pA, const double* pB) {
    return pA[0] * pB[0] + pA[1] * pB[1] + pA[2] * pB[2];
}

double spj_dist(const EdgeVec3& rA, const EdgeVec3& rB) {
    return std::sqrt((rA[0] - rB[0]) * (rA[0] - rB[0]) + (rA[1] - rB[1]) * (rA[1] - rB[1]) +
                     (rA[2] - rB[2]) * (rA[2] - rB[2]));
}

// Seeds per direction: (degree + 1) per span of nonzero length, 8 to 32.
int spj_count(const NurbsPatch& rPatch, std::size_t D) {
    if (D >= rPatch.mKnots.size() || D >= rPatch.mDegrees.size())
        return 8;
    const auto& r_k = rPatch.mKnots[D];
    int spans = 0;
    for (std::size_t i = 0; i + 1 < r_k.size(); ++i)
        spans += r_k[i + 1] > r_k[i] ? 1 : 0;
    return std::clamp((rPatch.mDegrees[D] + 1) * std::max(spans, 1), 8, 32);
}

}  // namespace

SurfaceProjector::SurfaceProjector(const BrepSurfaceEval& rSurface, const NurbsPatch& rPatch)
    : pSurface(&rSurface) {
    const auto& r_b = rSurface.Box();
    const int nu = spj_count(rPatch, 0), nv = spj_count(rPatch, 1);
    for (int j = 0; j <= nv; ++j)
        for (int i = 0; i <= nu; ++i) {
            const std::array<double, 2> uv{r_b[0] + (r_b[1] - r_b[0]) * double(i) / nu,
                                           r_b[2] + (r_b[3] - r_b[2]) * double(j) / nv};
            EdgeVec3 x{};
            rSurface.Eval(uv[0], uv[1], x.data(), nullptr, nullptr);
            mSeedUv.push_back(uv);
            mSeedX.push_back(x);
        }
}

double SurfaceProjector::Polish(const EdgeVec3& rP, std::array<double, 2>& rUv) const {
    const auto& r_b = pSurface->Box();
    const double eps_u = 1e-14 * std::max(1.0, r_b[1] - r_b[0]);
    const double eps_v = 1e-14 * std::max(1.0, r_b[3] - r_b[2]);
    const auto dist = [&](double U, double V) {
        EdgeVec3 x{};
        pSurface->Eval(U, V, x.data(), nullptr, nullptr);
        return spj_dist(x, rP);
    };
    double best = dist(rUv[0], rUv[1]);
    for (int it = 0; it < 50; ++it) {
        double s[3], su[3], sv[3], suu[3], suv[3], svv[3];
        pSurface->Eval2(rUv[0], rUv[1], s, su, sv, suu, suv, svv);
        const double r[3] = {s[0] - rP[0], s[1] - rP[1], s[2] - rP[2]};
        const double g0 = spj_dot(r, su), g1 = spj_dot(r, sv);
        double h00 = spj_dot(su, su) + spj_dot(r, suu);
        double h01 = spj_dot(su, sv) + spj_dot(r, suv);
        double h11 = spj_dot(sv, sv) + spj_dot(r, svv);
        double det = h00 * h11 - h01 * h01;
        if (!(h00 > 0.0) || !(det > 0.0)) {
            h00 = spj_dot(su, su);
            h01 = spj_dot(su, sv);
            h11 = spj_dot(sv, sv);
            det = h00 * h11 - h01 * h01;
            if (!(det > 0.0))
                break;
        }
        double du = -(h11 * g0 - h01 * g1) / det;
        double dv = -(h00 * g1 - h01 * g0) / det;
        double nu = std::clamp(rUv[0] + du, r_b[0], r_b[1]);
        double nv = std::clamp(rUv[1] + dv, r_b[2], r_b[3]);
        double d = dist(nu, nv);
        for (int h = 0; h < 30 && d > best; ++h) {
            du *= 0.5;
            dv *= 0.5;
            nu = std::clamp(rUv[0] + du, r_b[0], r_b[1]);
            nv = std::clamp(rUv[1] + dv, r_b[2], r_b[3]);
            d = dist(nu, nv);
        }
        if (d > best)
            break;
        const double step_u = std::fabs(nu - rUv[0]), step_v = std::fabs(nv - rUv[1]);
        rUv = {nu, nv};
        best = d;
        if (step_u <= eps_u && step_v <= eps_v)
            break;
    }
    return best;
}

std::array<double, 2> SurfaceProjector::Project(const EdgeVec3& rP, EdgeVec3& rX) const {
    // The three nearest seeds, ties to the lower index.
    constexpr std::size_t kStarts = 3;
    std::array<std::size_t, kStarts> near{};
    std::array<double, kStarts> near_d{};
    near_d.fill(std::numeric_limits<double>::infinity());
    for (std::size_t s = 0; s < mSeedX.size(); ++s) {
        const double d = spj_dist(mSeedX[s], rP);
        for (std::size_t k = 0; k < kStarts; ++k)
            if (d < near_d[k]) {
                for (std::size_t m = kStarts - 1; m > k; --m) {
                    near_d[m] = near_d[m - 1];
                    near[m] = near[m - 1];
                }
                near_d[k] = d;
                near[k] = s;
                break;
            }
    }
    std::array<double, 2> best_uv =
        mSeedUv.empty() ? std::array<double, 2>{0.0, 0.0} : mSeedUv[near[0]];
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < kStarts && k < mSeedUv.size(); ++k) {
        std::array<double, 2> uv = mSeedUv[near[k]];
        const double d = Polish(rP, uv);
        if (d < best) {
            best = d;
            best_uv = uv;
        }
    }
    pSurface->Eval(best_uv[0], best_uv[1], rX.data(), nullptr, nullptr);
    return best_uv;
}

}  // namespace meshioplusplus::detail
