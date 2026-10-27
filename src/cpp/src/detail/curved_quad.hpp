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
 * @brief Biquadratic and bicubic quadrilaterals (`quad9`, `quad16`) in 3-D,
 * and whether one is valid everywhere inside (roadmap 6.8.3).
 *
 *  - **The element.** Nodes in the layout of `quad9` (Gmsh's and VTK's):
 *    corners 0..3 at reference (0, 0), (1, 0), (1, 1), (0, 1); the edges'
 *    midpoints 4..7 on edges (0, 1), (1, 2), (2, 3), (3, 0); the centre 8.
 *    `quad16` in Gmsh's layout: the corners, each edge's two nodes from its
 *    first corner (4..11), then the inside (12..15) at (1/3, 1/3),
 *    (2/3, 1/3), (2/3, 2/3), (1/3, 2/3). The map is the tensor product of
 *    the Lagrange basis of degree `Order` at equispaced points.
 *  - **The Jacobian.** As for the curved triangles (`curved_triangle.hpp`),
 *    `J = (x_xi x x_eta) . N` against the face's normal `N`. `x_xi` is of
 *    degree (p - 1, p) and `x_eta` of degree (p, p - 1), so `J` is a
 *    tensor-product polynomial of degree (2p - 1, 2p - 1): (3, 3) for
 *    `quad9`, (5, 5) for `quad16`.
 *  - **Validity** (Johnen, Remacle and Geuzaine 2013, whose bounds hold for
 *    the tensor-product Bézier basis as for the simplex one). `J`'s 16 (36)
 *    Bernstein coefficients bound it: all positive proves `J > 0`
 *    everywhere, a value not positive proves the element invalid, and
 *    otherwise the square is split into four and each part judged, to
 *    `MaxDepth`; undecided at that depth counts as invalid.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>

// Project includes (private, not installed)
#include "curved_triangle.hpp"

namespace meshioplusplus::detail {

/// Nodes of a quad of order `Order` (2: 9, 3: 16).
std::size_t curved_quad_num_nodes(int Order);

/// Node `Node`'s reference position.
std::array<double, 2> curved_quad_node_position(int Order, std::size_t Node);

/// The element's point at reference `(Xi, Eta)`; `pNodes` holds
/// `curved_quad_num_nodes(Order)` points.
CurvedVec3 curved_quad_point(int Order, const CurvedVec3* pNodes, double Xi, double Eta);

double curved_quad_jacobian(int Order, const CurvedVec3* pNodes, const CurvedVec3& rNormal,
                            double Xi, double Eta);

CurvedTriangleCheck curved_quad_check(int Order, const CurvedVec3* pNodes,
                                      const CurvedVec3& rNormal, int MaxDepth = 8);

/// A lower bound of the scaled Jacobian, as `curved_triangle_scaled_bound`:
/// the least Bernstein coefficient of `J` over the largest in magnitude.
double curved_quad_scaled_bound(int Order, const CurvedVec3* pNodes, const CurvedVec3& rNormal);

/// `quad9` (order 2), as before.
CurvedVec3 curved_quad_point(const CurvedVec3* pNodes, double Xi, double Eta);

double curved_quad_jacobian(const CurvedVec3* pNodes, const CurvedVec3& rNormal, double Xi,
                            double Eta);

CurvedTriangleCheck curved_quad_check(const CurvedVec3* pNodes, const CurvedVec3& rNormal,
                                      int MaxDepth = 8);

}  // namespace meshioplusplus::detail
