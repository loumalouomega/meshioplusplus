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
 * @file detail/curved_quad.hpp
 * @brief Biquadratic quadrilaterals (`quad9`) in 3-D, and whether one is
 * valid everywhere inside (roadmap 6.8.3).
 *
 *  - **The element.** Nodes in the layout of `quad9` (Gmsh's and VTK's):
 *    corners 0..3 at reference (0, 0), (1, 0), (1, 1), (0, 1); the edges'
 *    midpoints 4..7 on edges (0, 1), (1, 2), (2, 3), (3, 0); the centre 8.
 *    The map is the tensor product of the quadratic Lagrange basis at 0,
 *    1/2 and 1 in each direction.
 *  - **The Jacobian.** As for the curved triangles (`curved_triangle.hpp`),
 *    `J = (x_xi x x_eta) . N` against the face's normal `N`. `x_xi` is of
 *    degree (1, 2) and `x_eta` of degree (2, 1), so `J` is a tensor-product
 *    polynomial of degree (3, 3).
 *  - **Validity** (Johnen, Remacle and Geuzaine 2013, whose bounds hold for
 *    the tensor-product Bézier basis as for the simplex one). `J`'s 16
 *    Bernstein coefficients bound it: all positive proves `J > 0`
 *    everywhere, a value not positive proves the element invalid, and
 *    otherwise the square is split into four and each part judged, to
 *    `MaxDepth`; undecided at that depth counts as invalid.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>

// Project includes (private, not installed)
#include "curved_triangle.hpp"

namespace meshioplusplus::detail {

/// The element's point at reference `(Xi, Eta)`; `pNodes` holds 9 points.
CurvedVec3 curved_quad_point(const CurvedVec3* pNodes, double Xi, double Eta);

/// `J = (x_xi x x_eta) . N` at reference `(Xi, Eta)`.
double curved_quad_jacobian(const CurvedVec3* pNodes, const CurvedVec3& rNormal, double Xi,
                            double Eta);

/// Whether the element is valid everywhere inside (see the file comment).
CurvedTriangleCheck curved_quad_check(const CurvedVec3* pNodes, const CurvedVec3& rNormal,
                                      int MaxDepth = 8);

}  // namespace meshioplusplus::detail
