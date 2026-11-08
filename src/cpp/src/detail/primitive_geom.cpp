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
#include <cstdint>
#include <limits>

// Project includes (private, not installed)
#include "primitive_geom.hpp"
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

}  // namespace detail
}  // namespace meshioplusplus
