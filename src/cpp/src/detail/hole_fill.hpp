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
 * @file detail/hole_fill.hpp
 * @brief A faired patch over one boundary loop (roadmap 7.8.6), after Liepa
 * (2003): the minimum-area triangulation of the loop (dynamic programming
 * over its vertices, O(n^3)), refined to the rim's edge density (a triangle
 * is split at its centroid while `sqrt(2)` times the centroid's distance to
 * each corner exceeds both the centroid's and that corner's spacing; after
 * each pass, the patch's interior edges are flipped to the Delaunay
 * criterion, opposite angles summing to at most pi), and faired with the rim
 * fixed: a membrane (`Membrane`, each new point the mean of its neighbours,
 * the uniform Laplacian) or a thin plate (`ThinPlate`, the umbrella vectors
 * of the new points and the rim, whose rings run into the surrounding
 * surface, minimised in least squares, so the patch continues the
 * surface's curvature). Each linear system is solved by `csr_pcg`.
 *
 * A **core-private** header: no installed header names it. Serial and
 * deterministic: ties go to the lowest index, the flips run in triangle
 * order.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus::detail {

enum class HoleFairing { None, Membrane, ThinPlate };

struct HoleFillInput {
    /// Every point of the surface; the loop's and the rings' ids index it.
    const std::vector<Vec3>* pPoints = nullptr;
    /// The loop, in the order whose fan `(v_i, v_{i+1}, c)` winds as the
    /// surface does.
    std::vector<std::int64_t> mLoop;
    /// Per loop vertex, its neighbours on the surface (outside the hole).
    std::vector<std::vector<std::int64_t>> mRimRing;
    /// Per loop vertex, the mean length of its surface edges.
    std::vector<double> mRimSpacing;
    HoleFairing mFairing = HoleFairing::Membrane;
    bool mRefine = true;
};

struct HoleFillPatch {
    /// The new points; the k-th has the id `Base + k`.
    std::vector<Vec3> mPoints;
    /// The patch's triangles, wound as the surface, in point ids.
    std::vector<std::array<std::int64_t, 3>> mTriangles;
};

/// The faired patch over @p rIn's loop (see the file comment); the new points
/// are numbered from @p Base.
HoleFillPatch hole_fill(const HoleFillInput& rIn, std::int64_t Base);

}  // namespace meshioplusplus::detail
