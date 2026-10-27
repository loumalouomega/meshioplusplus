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
 * @file detail/curved_untangle.hpp
 * @brief Untangling curved surface triangles and quadrilaterals by moving
 * their nodes on the surface (roadmap 6.6.3).
 *
 * A curved element whose Jacobian is not positive everywhere (the Bézier
 * check of `curved_triangle.hpp`) loses the order its nodes bought if it is
 * made straight-sided. Toulorge, Geuzaine, Remacle and Lambrechts ("Robust
 * untangling of curvilinear meshes", J. Comput. Phys. 254, 2013) instead move
 * the high-order nodes that are not on the boundary to raise the worst
 * scaled Jacobian of the elements around them, the boundary nodes staying on
 * the CAD. Here the free nodes are a face's own (interior edges' nodes,
 * cubic centres and, for split quads, the corners inside the face) and
 * each moves in the face's parameters, lifted by the surface, so every node
 * stays on the geometry by construction; the measure raised is
 * `curved_triangle_scaled_bound` (`curved_quad_scaled_bound` for a quad), a
 * lower bound of `min J / max |J|`. Any node of an element may be free,
 * corners included. The optimiser is a deterministic Gauss-Seidel pattern
 * search (the precedent of `hex_untangle`), not the paper's barrier method:
 * each free node of an invalid element, in ascending id order, takes the
 * best of eight steps around it that raises the least bound of its elements,
 * and halves its step when none does; sweeps repeat until every element is
 * proved valid, no step is left, or the sweep cap. Only elements around a
 * moved node can change, so the result does not depend on threads.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

// Project includes (private, not installed)
#include "curved_triangle.hpp"

namespace meshioplusplus::detail {

/// One element: its node ids in the `triangle6` / `triangle10` layout, or
/// with `mQuad` in the `quad9` / `quad16` one (`curved_quad.hpp`), and the
/// face's normal the sign of its Jacobian is taken against.
struct CurvedUntangleElement {
    std::array<std::int64_t, 16> mNodes{};
    CurvedVec3 mNormal{0.0, 0.0, 0.0};
    bool mQuad = false;
};

/// A curved triangle as an element to untangle.
inline CurvedUntangleElement curved_untangle_triangle(const std::array<std::int64_t, 10>& rNodes,
                                                      const CurvedVec3& rNormal) {
    CurvedUntangleElement e;
    for (std::size_t i = 0; i < rNodes.size(); ++i)
        e.mNodes[i] = rNodes[i];
    e.mNormal = rNormal;
    return e;
}

struct CurvedUntangleReport {
    std::int64_t mInvalidBefore = 0;  ///< elements not proved valid on entry
    std::int64_t mInvalidAfter = 0;   ///< and on return
    std::int64_t mMoves = 0;          ///< node moves kept
    std::int64_t mSweeps = 0;
};

/// A free node's lift: the point at parameters `rUv` of free node `Free`.
using CurvedUntangleLift =
    std::function<CurvedVec3(std::size_t Free, const std::array<double, 2>& rUv)>;

/**
 * @brief Move the free nodes `rFree` (ids into `pX`, three coordinates each)
 * whose parameters are `rUv` (updated in place) to untangle `rElems`;
 * `rLift` gives the surface point at parameters. Every other node of the
 * elements is fixed. Elements with no free node are judged, not changed.
 */
CurvedUntangleReport curved_untangle(
    int Order, const std::vector<CurvedUntangleElement>& rElems, double* pX,
    const std::vector<std::int64_t>& rFree, std::vector<std::array<double, 2>>& rUv,
    const std::function<CurvedVec3(const std::array<double, 2>&)>& rLift, int MaxSweeps = 60);

/**
 * @brief The same with a lift of its own per free node and, when `pBounds`
 * is given, a box `{lo0, hi0, lo1, hi1}` per free node its parameters stay
 * in. A node whose lift ignores its second parameter (a node on a curve,
 * `{t, 0}`) has no speed along it, so it moves along the curve only.
 */
CurvedUntangleReport curved_untangle(int Order, const std::vector<CurvedUntangleElement>& rElems,
                                     double* pX, const std::vector<std::int64_t>& rFree,
                                     std::vector<std::array<double, 2>>& rUv,
                                     const CurvedUntangleLift& rLift,
                                     const std::vector<std::array<double, 4>>* pBounds,
                                     int MaxSweeps = 60);

}  // namespace meshioplusplus::detail
