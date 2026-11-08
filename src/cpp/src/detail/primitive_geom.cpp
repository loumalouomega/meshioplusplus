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

// Point queries against a fitted primitive. See detail/primitive_geom.hpp.

// System includes
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

// Project includes (private, not installed)
#include "primitive_geom.hpp"
#include "lsq_fit.hpp"
#include "math_utils.hpp"

namespace meshioplusplus {
namespace detail {
namespace {

constexpr double kPgPi = 3.14159265358979323846;

// The radial unit vector of pP about the axis (pX, pD), its height along the
// axis and its distance from it; false on the axis.
bool pg_radial(const double* pX, const double* pD, const double* pP, double* pR, double& rH,
               double& rRho) {
    const double v[3] = {pP[0] - pX[0], pP[1] - pX[1], pP[2] - pX[2]};
    rH = detail::MathUtils::Dot(v, pD);
    for (int k = 0; k < 3; ++k)
        pR[k] = v[k] - rH * pD[k];
    rRho = std::sqrt(detail::MathUtils::Dot(pR, pR));
    if (!(rRho > 0.0))
        return false;
    for (int k = 0; k < 3; ++k)
        pR[k] /= rRho;
    return true;
}

}  // namespace

double primitive_point_distance(std::int64_t Type, const double* pParams, const double* pP) {
    const double* x = pParams;
    const double* d = pParams + 3;
    const double v[3] = {pP[0] - x[0], pP[1] - x[1], pP[2] - x[2]};
    switch (Type) {
        case 1:
            return std::fabs(detail::MathUtils::Dot(v, d));
        case 2:
            return std::fabs(std::sqrt(detail::MathUtils::Dot(v, v)) - pParams[6]);
        case 3:
        case 4:
        case 5: {
            const double h = detail::MathUtils::Dot(v, d);
            const double r[3] = {v[0] - h * d[0], v[1] - h * d[1], v[2] - h * d[2]};
            const double rho = std::sqrt(detail::MathUtils::Dot(r, r));
            if (Type == 3)
                return std::fabs(rho - pParams[6]);
            if (Type == 4) {
                const double a = pParams[8] * kPgPi / 180.0;
                if (h * std::cos(a) + rho * std::sin(a) < 0.0)
                    return std::sqrt(detail::MathUtils::Dot(v, v));  // behind the apex
                return std::fabs(rho * std::cos(a) - h * std::sin(a));
            }
            const double radial = rho - pParams[6];
            return std::fabs(std::sqrt(radial * radial + h * h) - pParams[7]);
        }
        default:
            return std::numeric_limits<double>::infinity();
    }
}

bool primitive_project(std::int64_t Type, const double* pParams, const double* pP, double* pOut) {
    for (int k = 0; k < 3; ++k)
        pOut[k] = pP[k];
    const double* x = pParams;
    const double* d = pParams + 3;
    switch (Type) {
        case 1: {
            const double v[3] = {pP[0] - x[0], pP[1] - x[1], pP[2] - x[2]};
            const double h = detail::MathUtils::Dot(v, d);
            for (int k = 0; k < 3; ++k)
                pOut[k] = pP[k] - h * d[k];
            return true;
        }
        case 2: {
            const double v[3] = {pP[0] - x[0], pP[1] - x[1], pP[2] - x[2]};
            const double len = std::sqrt(detail::MathUtils::Dot(v, v));
            if (!(len > 0.0))
                return false;
            for (int k = 0; k < 3; ++k)
                pOut[k] = x[k] + pParams[6] * v[k] / len;
            return true;
        }
        case 3: {
            double r[3], h, rho;
            if (!pg_radial(x, d, pP, r, h, rho))
                return false;
            for (int k = 0; k < 3; ++k)
                pOut[k] = x[k] + h * d[k] + pParams[6] * r[k];
            return true;
        }
        case 4: {
            double r[3], h, rho;
            if (!pg_radial(x, d, pP, r, h, rho))
                return false;
            const double a = pParams[8] * kPgPi / 180.0;
            const double t = h * std::cos(a) + rho * std::sin(a);
            if (t <= 0.0) {
                for (int k = 0; k < 3; ++k)
                    pOut[k] = x[k];
                return true;
            }
            for (int k = 0; k < 3; ++k)
                pOut[k] = x[k] + t * (std::cos(a) * d[k] + std::sin(a) * r[k]);
            return true;
        }
        case 5: {
            double r[3], h, rho;
            if (!pg_radial(x, d, pP, r, h, rho))
                return false;
            double q[3], w[3];
            for (int k = 0; k < 3; ++k) {
                q[k] = x[k] + pParams[6] * r[k];
                w[k] = pP[k] - q[k];
            }
            const double len = std::sqrt(detail::MathUtils::Dot(w, w));
            if (!(len > 0.0))
                return false;
            for (int k = 0; k < 3; ++k)
                pOut[k] = q[k] + pParams[7] * w[k] / len;
            return true;
        }
        default:
            return false;
    }
}

bool primitive_normal(std::int64_t Type, const double* pParams, const double* pP, double* pOut) {
    const double* x = pParams;
    const double* d = pParams + 3;
    switch (Type) {
        case 1:
            for (int k = 0; k < 3; ++k)
                pOut[k] = d[k];
            return true;
        case 2: {
            const double v[3] = {pP[0] - x[0], pP[1] - x[1], pP[2] - x[2]};
            const double len = std::sqrt(detail::MathUtils::Dot(v, v));
            if (!(len > 0.0))
                return false;
            for (int k = 0; k < 3; ++k)
                pOut[k] = v[k] / len;
            return true;
        }
        case 3: {
            double h, rho;
            return pg_radial(x, d, pP, pOut, h, rho);
        }
        case 4: {
            double r[3], h, rho;
            if (!pg_radial(x, d, pP, r, h, rho))
                return false;
            const double a = pParams[8] * kPgPi / 180.0;
            for (int k = 0; k < 3; ++k)
                pOut[k] = std::cos(a) * r[k] - std::sin(a) * d[k];
            return true;
        }
        case 5: {
            double r[3], h, rho;
            if (!pg_radial(x, d, pP, r, h, rho))
                return false;
            double w[3];
            for (int k = 0; k < 3; ++k)
                w[k] = pP[k] - (x[k] + pParams[6] * r[k]);
            const double len = std::sqrt(detail::MathUtils::Dot(w, w));
            if (!(len > 0.0))
                return false;
            for (int k = 0; k < 3; ++k)
                pOut[k] = w[k] / len;
            return true;
        }
        default:
            return false;
    }
}

double primitive_signed_distance(std::int64_t Type, const double* pParams, const double* pP) {
    const double* x = pParams;
    const double* d = pParams + 3;
    const double v[3] = {pP[0] - x[0], pP[1] - x[1], pP[2] - x[2]};
    switch (Type) {
        case 1:
            return detail::MathUtils::Dot(v, d);
        case 2:
            return std::sqrt(detail::MathUtils::Dot(v, v)) - pParams[6];
        case 3:
        case 4:
        case 5: {
            const double h = detail::MathUtils::Dot(v, d);
            const double r[3] = {v[0] - h * d[0], v[1] - h * d[1], v[2] - h * d[2]};
            const double rho = std::sqrt(detail::MathUtils::Dot(r, r));
            if (Type == 3)
                return rho - pParams[6];
            if (Type == 4) {
                const double a = pParams[8] * kPgPi / 180.0;
                return rho * std::cos(a) - h * std::sin(a);
            }
            const double radial = rho - pParams[6];
            return std::sqrt(radial * radial + h * h) - pParams[7];
        }
        default:
            return std::numeric_limits<double>::quiet_NaN();
    }
}

PrimitiveFormMetrics primitive_form_metrics(std::int64_t Type, const double* pParams,
                                            const double* pPts, std::size_t N) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    PrimitiveFormMetrics m{nan, nan, nan, nan, nan, {nan, nan}};
    if (Type < 1 || Type > 5 || N == 0)
        return m;
    double lo = std::numeric_limits<double>::infinity();
    double hi = -lo;
    for (std::size_t i = 0; i < N; ++i) {
        const double s = primitive_signed_distance(Type, pParams, pPts + 3 * i);
        lo = s < lo ? s : lo;
        hi = s > hi ? s : hi;
    }
    m.mFormError = hi - lo;
    const double* x = pParams;
    const double* d = pParams + 3;
    if (Type == 1) {
        // The in-plane principal directions of the points' scatter.
        double mean[3] = {0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < N; ++i)
            for (int k = 0; k < 3; ++k)
                mean[k] += pPts[3 * i + static_cast<std::size_t>(k)];
        for (int k = 0; k < 3; ++k)
            mean[k] /= static_cast<double>(N);
        double cov[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        for (std::size_t i = 0; i < N; ++i) {
            const double q[3] = {pPts[3 * i] - mean[0], pPts[3 * i + 1] - mean[1],
                                 pPts[3 * i + 2] - mean[2]};
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    cov[a * 3 + b] += q[a] * q[b];
        }
        double w[3], vec[9];
        detail::lsq_sym3_eigen(cov, w, vec);
        for (int e = 0; e < 2; ++e) {
            const double* u = vec + 3 * (2 - e);  // largest eigenvalue first
            double a = std::numeric_limits<double>::infinity(), b = -a;
            for (std::size_t i = 0; i < N; ++i) {
                const double t = detail::MathUtils::Dot(pPts + 3 * i, u);
                a = t < a ? t : a;
                b = t > b ? t : b;
            }
            m.mExtent[e] = b - a;
        }
        return m;
    }
    if (Type == 2) {
        double rlo = std::numeric_limits<double>::infinity(), rhi = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            const double* p = pPts + 3 * i;
            const double v[3] = {p[0] - x[0], p[1] - x[1], p[2] - x[2]};
            const double r = std::sqrt(detail::MathUtils::Dot(v, v));
            rlo = r < rlo ? r : rlo;
            rhi = r > rhi ? r : rhi;
        }
        m.mDiameterMin = 2.0 * rlo;
        m.mDiameterMax = 2.0 * rhi;
        return m;
    }
    if (Type == 5)
        return m;
    // Cylinder and cone: heights, radial distances and angles about the axis.
    double e1[3], e2[3];
    detail::lsq_frame(d, e1, e2);
    double hlo = std::numeric_limits<double>::infinity(), hhi = -hlo;
    double rlo = hlo, rhi = 0.0;
    std::vector<double> angle;
    angle.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const double* p = pPts + 3 * i;
        const double v[3] = {p[0] - x[0], p[1] - x[1], p[2] - x[2]};
        const double h = detail::MathUtils::Dot(v, d);
        const double r[3] = {v[0] - h * d[0], v[1] - h * d[1], v[2] - h * d[2]};
        const double rho = std::sqrt(detail::MathUtils::Dot(r, r));
        hlo = h < hlo ? h : hlo;
        hhi = h > hhi ? h : hhi;
        rlo = rho < rlo ? rho : rlo;
        rhi = rho > rhi ? rho : rhi;
        if (rho > 0.0) {
            double a = std::atan2(detail::MathUtils::Dot(r, e2), detail::MathUtils::Dot(r, e1)) *
                       180.0 / kPgPi;
            angle.push_back(a < 0.0 ? a + 360.0 : a);
        }
    }
    m.mLength = hhi - hlo;
    if (Type == 3) {
        m.mDiameterMin = 2.0 * rlo;
        m.mDiameterMax = 2.0 * rhi;
    }
    // Points on several rims share their angles: keep one of each.
    std::sort(angle.begin(), angle.end());
    angle.erase(
        std::unique(angle.begin(), angle.end(), [](double A, double B) { return B - A <= 1e-9; }),
        angle.end());
    if (angle.size() >= 2) {
        std::vector<double> gap(angle.size());
        for (std::size_t i = 0; i + 1 < angle.size(); ++i)
            gap[i] = angle[i + 1] - angle[i];
        gap.back() = angle.front() + 360.0 - angle.back();
        const double largest = *std::max_element(gap.begin(), gap.end());
        std::vector<double> sorted = gap;
        std::nth_element(sorted.begin(),
                         sorted.begin() + static_cast<std::ptrdiff_t>(sorted.size() / 2),
                         sorted.end());
        const double median = sorted[sorted.size() / 2];
        m.mArcDeg = largest <= 2.0 * median ? 360.0 : 360.0 - largest;
    }
    return m;
}

}  // namespace detail
}  // namespace meshioplusplus
