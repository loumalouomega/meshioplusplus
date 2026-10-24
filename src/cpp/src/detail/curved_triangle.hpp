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
 * @file detail/curved_triangle.hpp
 * @brief Lagrange triangles of order 2 and 3 in 3-D, and whether one is valid
 * everywhere inside, not only at its nodes (roadmap 6.6.1).
 *
 *  - **The element.** Nodes in the layout of `triangle6` and `triangle10`
 *    (Gmsh's and meshio's): the corners 0, 1, 2 at reference (0, 0), (1, 0),
 *    (0, 1); then each edge's nodes from its first corner, edges (0, 1),
 *    (1, 2), (2, 0); for order 3 then the centre. The map is
 *    `x(xi, eta) = sum_i N_i(xi, eta) X_i` with the Lagrange shape functions.
 *  - **The Jacobian.** A surface element has no square Jacobian matrix; its
 *    signed area element against a reference direction `N` (the surface
 *    normal of the face it lies on) is `J = (x_xi x x_eta) . N`, a
 *    polynomial of degree `2 (p - 1)` in the reference coordinates.
 *  - **Validity** (Johnen, Remacle and Geuzaine 2013). `J` positive at the
 *    nodes is necessary, not sufficient: a curved element can fold between
 *    them. `J` is expanded in the Bernstein basis of its degree on the
 *    triangle; its values lie in the convex hull of the coefficients, so all
 *    coefficients positive proves `J > 0` everywhere, and any value not
 *    positive proves the element invalid. Otherwise the triangle is split
 *    into four and each part judged, to `MaxDepth`; undecided at that depth
 *    counts as invalid. Exact up to the depth and floating-point rounding.
 */

// System includes
#include <array>

namespace meshioplusplus::detail {

using CurvedVec3 = std::array<double, 3>;

/// Nodes of a Lagrange triangle of order 1, 2 or 3 (3, 6 or 10).
int curved_triangle_num_nodes(int Order);

/// The element's point at reference `(Xi, Eta)`.
CurvedVec3 curved_triangle_point(int Order, const CurvedVec3* pNodes, double Xi, double Eta);

/// `J = (x_xi x x_eta) . N` at reference `(Xi, Eta)`; `pNodes` holds
/// `curved_triangle_num_nodes(Order)` points.
double curved_triangle_jacobian(int Order, const CurvedVec3* pNodes, const CurvedVec3& rNormal,
                                double Xi, double Eta);

/// The verdict on one element.
struct CurvedTriangleCheck {
    bool mValid = false;   ///< J > 0 everywhere, proved
    bool mDecided = true;  ///< false: undecided at the depth (counted invalid)
    double mMinNodal = 0;  ///< the smallest J at the element's own nodes
};

/// Whether the element is valid everywhere inside (see the file comment).
CurvedTriangleCheck curved_triangle_check(int Order, const CurvedVec3* pNodes,
                                          const CurvedVec3& rNormal, int MaxDepth = 8);

/// The Bernstein coefficients of the polynomial of degree `Degree` (2 or 4)
/// on the reference triangle whose values at the degree's lattice points are
/// `pValues`, in the order `(i, j)` with `i` (the `xi` index) outer, `j`
/// ascending, `i + j <= Degree`. Exposed for tests.
void curved_triangle_bernstein(int Degree, const double* pValues, double* pCoefficients);

}  // namespace meshioplusplus::detail
