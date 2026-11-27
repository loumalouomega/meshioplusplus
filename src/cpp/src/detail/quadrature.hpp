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
#pragma once

/**
 * @file detail/quadrature.hpp
 * @brief Shared Gauss quadrature rules on reference elements (roadmap 3.4.3 / 15.1).
 *
 * All element rules live in the `[0, 1]` reference coordinates used by
 * `detail/reference_element.hpp` (`lagrange::vtk_lattice`):
 *   - Line: `[0, 1]`, measure 1.
 *   - Triangle: `u, v >= 0, u + v <= 1`, measure 1/2.
 *   - Quad: `[0, 1]^2`, measure 1.
 *   - Tetra: `u, v, w >= 0, u + v + w <= 1`, measure 1/6.
 *   - Hexahedron: `[0, 1]^3`, measure 1.
 *   - Wedge: `Triangle x [0, 1]`, measure 1/2.
 *   - Pyramid: `0 <= w <= 1, 0 <= u, v <= 1 - w`, measure 1/3.
 *
 * Weights sum to the reference element's measure so that
 * `sum_q w_q * |det J(xi_q)| * f(xi_q)` approximates `\int_K f dK`.
 */

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus::detail {

/// A quadrature rule in `[0, 1]` reference coordinates.
struct QuadratureRule {
    std::vector<std::array<double, 3>> mPoints;
    std::vector<double> mWeights;
};

/**
 * @brief @p N-point Gauss-Legendre nodes and weights on `[-1, 1]`, sorted
 * ascending by node. Exact for polynomials of degree `2N - 1`.
 */
inline void gauss_legendre_1d(int N, std::vector<double>& rNodes, std::vector<double>& rWeights) {
    rNodes.resize(static_cast<std::size_t>(N));
    rWeights.resize(static_cast<std::size_t>(N));
    if (N <= 0)
        return;
    constexpr double kPi = 3.14159265358979323846;
    const int m = (N + 1) / 2;
    for (int i = 0; i < m; ++i) {
        double z = std::cos(kPi * (static_cast<double>(i) + 0.75) / (static_cast<double>(N) + 0.5));
        double pp = 0.0;
        for (int iter = 0; iter < 50; ++iter) {
            double p1 = 1.0;
            double p2 = 0.0;
            for (int j = 1; j <= N; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / static_cast<double>(j);
            }
            pp = static_cast<double>(N) * (z * p1 - p2) / (z * z - 1.0);
            const double dz = p1 / pp;
            z -= dz;
            if (std::abs(dz) <= 1e-16)
                break;
        }
        const double w = 2.0 / ((1.0 - z * z) * pp * pp);
        rNodes[static_cast<std::size_t>(i)] = -z;
        rNodes[static_cast<std::size_t>(N - 1 - i)] = z;
        rWeights[static_cast<std::size_t>(i)] = w;
        rWeights[static_cast<std::size_t>(N - 1 - i)] = w;
    }
}

/// @p N-point Gauss-Legendre nodes and weights mapped to `[0, 1]`.
inline void gauss_legendre_unit(int N, std::vector<double>& rNodes, std::vector<double>& rWeights) {
    gauss_legendre_1d(N, rNodes, rWeights);
    for (std::size_t i = 0; i < rNodes.size(); ++i) {
        rNodes[i] = 0.5 * (rNodes[i] + 1.0);
        rWeights[i] *= 0.5;
    }
}

/// Symmetric triangle quadrature on the unit triangle (`u, v >= 0, u+v <= 1`, area 0.5).
inline QuadratureRule quadrature_triangle(int ExactOrder) {
    QuadratureRule rule;
    const auto add_sym3 = [&](double a, double b, double w_area_frac) {
        const double w = 0.5 * w_area_frac;
        const double c = 1.0 - a - b;
        rule.mPoints.push_back({a, b, 0.0});
        rule.mWeights.push_back(w);
        rule.mPoints.push_back({b, c, 0.0});
        rule.mWeights.push_back(w);
        rule.mPoints.push_back({c, a, 0.0});
        rule.mWeights.push_back(w);
    };
    const auto add_sym6 = [&](double a, double b, double w_area_frac) {
        const double w = 0.5 * w_area_frac;
        const double c = 1.0 - a - b;
        const double p[6][2] = {{a, b}, {b, a}, {b, c}, {c, b}, {c, a}, {a, c}};
        for (const auto& q : p) {
            rule.mPoints.push_back({q[0], q[1], 0.0});
            rule.mWeights.push_back(w);
        }
    };

    if (ExactOrder <= 1) {
        rule.mPoints.push_back({1.0 / 3.0, 1.0 / 3.0, 0.0});
        rule.mWeights.push_back(0.5);
        return rule;
    }
    if (ExactOrder <= 2) {
        add_sym3(1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0);
        return rule;
    }
    if (ExactOrder <= 4) {
        // Dunavant 6-point rule, exact to degree 4, all positive weights.
        add_sym3(0.44594849091596488632, 0.44594849091596488632, 0.22338158967801146570);
        add_sym3(0.09157621350977074346, 0.09157621350977074346, 0.10995174365532186764);
        return rule;
    }
    if (ExactOrder <= 5) {
        // Dunavant 7-point rule, exact to degree 5, all positive weights.
        rule.mPoints.push_back({1.0 / 3.0, 1.0 / 3.0, 0.0});
        rule.mWeights.push_back(0.5 * 0.225);
        add_sym3(0.05971587178976982046, 0.47014206410511508977, 0.13239415278850618074);
        add_sym3(0.79742698535308732240, 0.10128650732345633880, 0.12593918054482715260);
        return rule;
    }
    if (ExactOrder <= 6) {
        // Dunavant 12-point rule, exact to degree 6, all positive weights.
        add_sym3(0.50142650965817915742, 0.24928674517091042129, 0.11678627572637936603);
        add_sym3(0.87382197101699554332, 0.06308901449150222834, 0.05084490637020681692);
        add_sym6(0.05314504984481694735, 0.31035245103378440542, 0.08285107561837357519);
        return rule;
    }
    // General higher-order fallback: Duffy transform from [0, 1]^2 with u = (1-a)b, v = a,
    // Jacobian = (1 - a). Polynomial of degree p in (u,v) has degree <= p+1 in a, <= p in b.
    const int n = (ExactOrder + 2) / 2 + 1;
    std::vector<double> gx, gw;
    gauss_legendre_unit(n, gx, gw);
    for (int i = 0; i < n; ++i) {
        const double a = gx[static_cast<std::size_t>(i)];
        const double wa = gw[static_cast<std::size_t>(i)] * (1.0 - a);
        for (int j = 0; j < n; ++j) {
            const double b = gx[static_cast<std::size_t>(j)];
            rule.mPoints.push_back({(1.0 - a) * b, a, 0.0});
            rule.mWeights.push_back(wa * gw[static_cast<std::size_t>(j)]);
        }
    }
    return rule;
}

/// Symmetric tetrahedron quadrature on the unit tet (`u,v,w >= 0, u+v+w <= 1`, volume 1/6).
inline QuadratureRule quadrature_tetra(int ExactOrder) {
    QuadratureRule rule;
    constexpr double kVol = 1.0 / 6.0;
    const auto add_sym4 = [&](double a, double w_vol_frac) {
        const double b = 1.0 - 3.0 * a;
        const double w = kVol * w_vol_frac;
        rule.mPoints.push_back({a, a, a});
        rule.mWeights.push_back(w);
        rule.mPoints.push_back({b, a, a});
        rule.mWeights.push_back(w);
        rule.mPoints.push_back({a, b, a});
        rule.mWeights.push_back(w);
        rule.mPoints.push_back({a, a, b});
        rule.mWeights.push_back(w);
    };
    const auto add_sym6 = [&](double a, double w_vol_frac) {
        const double b = 0.5 - a;
        const double w = kVol * w_vol_frac;
        const double p[6][3] = {{a, a, b}, {a, b, a}, {b, a, a}, {a, b, b}, {b, a, b}, {b, b, a}};
        for (const auto& q : p) {
            rule.mPoints.push_back({q[0], q[1], q[2]});
            rule.mWeights.push_back(w);
        }
    };

    if (ExactOrder <= 1) {
        rule.mPoints.push_back({0.25, 0.25, 0.25});
        rule.mWeights.push_back(kVol);
        return rule;
    }
    if (ExactOrder <= 2) {
        // 4-point Keast rule, exact to degree 2, positive weights.
        const double a = (5.0 - std::sqrt(5.0)) / 20.0;
        add_sym4(a, 0.25);
        return rule;
    }
    if (ExactOrder <= 4) {
        // 11-point Keast / Shunn-Ham rule, exact to degree 4, all positive weights.
        rule.mPoints.push_back({0.25, 0.25, 0.25});
        rule.mWeights.push_back(kVol * (-0.07893333333333333));
        // Note: use 14-point positive-weight rule for degree 4 (Shunn & Ham 2012)
        // so all weights are strictly positive!
        rule.mPoints.clear();
        rule.mWeights.clear();
        add_sym4(0.09273525031089122640, 0.07349304311636194954);
        add_sym4(0.31088591926330060980, 0.11268792571801585080);
        add_sym6(0.04550370412564964949, 0.04254602077708146644);
        return rule;
    }
    // General higher-order fallback: Duffy transform from [0, 1]^3 with
    // w = a, v = (1 - a) b, u = (1 - a)(1 - b) c, Jacobian = (1 - a)^2 (1 - b).
    const int n = (ExactOrder + 3) / 2 + 1;
    std::vector<double> gx, gw;
    gauss_legendre_unit(n, gx, gw);
    for (int i = 0; i < n; ++i) {
        const double a = gx[static_cast<std::size_t>(i)];
        const double oma = 1.0 - a;
        const double wa = gw[static_cast<std::size_t>(i)] * oma * oma;
        for (int j = 0; j < n; ++j) {
            const double b = gx[static_cast<std::size_t>(j)];
            const double omb = 1.0 - b;
            const double wab = wa * gw[static_cast<std::size_t>(j)] * omb;
            for (int k = 0; k < n; ++k) {
                const double c = gx[static_cast<std::size_t>(k)];
                rule.mPoints.push_back({oma * omb * c, oma * b, a});
                rule.mWeights.push_back(wab * gw[static_cast<std::size_t>(k)]);
            }
        }
    }
    return rule;
}

/// Quadrature rule for a base cell family (`CellType`) exact to polynomial degree @p ExactOrder.
inline QuadratureRule quadrature_rule(CellType Type, int ExactOrder) {
    const int order = std::max(1, ExactOrder);
    const int n1d = (order + 2) / 2;  // 2*n1d - 1 >= order
    QuadratureRule rule;
    switch (Type) {
        case CellType::Line:
        case CellType::Line3:
        case CellType::VtkLagrangeCurve:
        case CellType::VtkBezierCurve: {
            std::vector<double> gx, gw;
            gauss_legendre_unit(n1d, gx, gw);
            for (int i = 0; i < n1d; ++i) {
                rule.mPoints.push_back({gx[static_cast<std::size_t>(i)], 0.0, 0.0});
                rule.mWeights.push_back(gw[static_cast<std::size_t>(i)]);
            }
            return rule;
        }
        case CellType::Triangle:
        case CellType::Triangle6:
        case CellType::Triangle7:
        case CellType::VtkLagrangeTriangle:
        case CellType::VtkBezierTriangle:
            return quadrature_triangle(order);

        case CellType::Quad:
        case CellType::Quad8:
        case CellType::Quad9:
        case CellType::VtkLagrangeQuadrilateral:
        case CellType::VtkBezierQuadrilateral: {
            std::vector<double> gx, gw;
            gauss_legendre_unit(n1d, gx, gw);
            for (int j = 0; j < n1d; ++j) {
                for (int i = 0; i < n1d; ++i) {
                    rule.mPoints.push_back(
                        {gx[static_cast<std::size_t>(i)], gx[static_cast<std::size_t>(j)], 0.0});
                    rule.mWeights.push_back(gw[static_cast<std::size_t>(i)] *
                                            gw[static_cast<std::size_t>(j)]);
                }
            }
            return rule;
        }
        case CellType::Tetra:
        case CellType::Tetra10:
        case CellType::VtkLagrangeTetrahedron:
        case CellType::VtkBezierTetrahedron:
            return quadrature_tetra(order);

        case CellType::Hexahedron:
        case CellType::Hexahedron20:
        case CellType::Hexahedron27:
        case CellType::VtkLagrangeHexahedron:
        case CellType::VtkBezierHexahedron: {
            std::vector<double> gx, gw;
            gauss_legendre_unit(n1d, gx, gw);
            for (int k = 0; k < n1d; ++k) {
                for (int j = 0; j < n1d; ++j) {
                    for (int i = 0; i < n1d; ++i) {
                        rule.mPoints.push_back({gx[static_cast<std::size_t>(i)],
                                                gx[static_cast<std::size_t>(j)],
                                                gx[static_cast<std::size_t>(k)]});
                        rule.mWeights.push_back(gw[static_cast<std::size_t>(i)] *
                                                gw[static_cast<std::size_t>(j)] *
                                                gw[static_cast<std::size_t>(k)]);
                    }
                }
            }
            return rule;
        }
        case CellType::Wedge:
        case CellType::Wedge15:
        case CellType::Wedge18:
        case CellType::VtkLagrangeWedge:
        case CellType::VtkBezierWedge: {
            const QuadratureRule tri = quadrature_triangle(order);
            std::vector<double> gz, wz;
            gauss_legendre_unit(n1d, gz, wz);
            for (int k = 0; k < n1d; ++k) {
                for (std::size_t t = 0; t < tri.mPoints.size(); ++t) {
                    rule.mPoints.push_back(
                        {tri.mPoints[t][0], tri.mPoints[t][1], gz[static_cast<std::size_t>(k)]});
                    rule.mWeights.push_back(tri.mWeights[t] * wz[static_cast<std::size_t>(k)]);
                }
            }
            return rule;
        }
        case CellType::Pyramid:
        case CellType::Pyramid13:
        case CellType::Pyramid14:
        case CellType::VtkLagrangePyramid:
        case CellType::VtkBezierPyramid: {
            // Duffy mapping from [0, 1]^3: w = c, u = (1-c)*a, v = (1-c)*b, J = (1-c)^2.
            const int nz = (order + 4) / 2;
            std::vector<double> gx, gw, gz, wz;
            gauss_legendre_unit(n1d, gx, gw);
            gauss_legendre_unit(nz, gz, wz);
            for (int k = 0; k < nz; ++k) {
                const double c = gz[static_cast<std::size_t>(k)];
                const double omc = 1.0 - c;
                const double wk = wz[static_cast<std::size_t>(k)] * omc * omc;
                for (int j = 0; j < n1d; ++j) {
                    for (int i = 0; i < n1d; ++i) {
                        rule.mPoints.push_back({omc * gx[static_cast<std::size_t>(i)],
                                                omc * gx[static_cast<std::size_t>(j)], c});
                        rule.mWeights.push_back(wk * gw[static_cast<std::size_t>(i)] *
                                                gw[static_cast<std::size_t>(j)]);
                    }
                }
            }
            return rule;
        }
        default:
            return rule;
    }
}

}  // namespace meshioplusplus::detail
