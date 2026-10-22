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
//  Provenance:      built from the published description of Labelle and
//                   Shewchuk, "Isosurface Stuffing: Fast Tetrahedral Meshes
//                   with Good Dihedral Angles" (SIGGRAPH 2007) only; the
//                   paper's reference implementation (Stellar) is not read.
//
#pragma once

/**
 * @file detail/isosurface_stuffing.hpp
 * @brief Isosurface stuffing of one interface through a tetrahedral mesh, for
 * `octree_mesh`.
 *
 * Every vertex carries a sign: -1 on the near side of the interface (inside),
 * +1 on the far side, 0 on it. Stuffing (Labelle and Shewchuk 2007) then:
 *
 *  1. places a *cut point* on every edge whose endpoints have opposite signs,
 *     where the interface crosses it (the caller's `StuffCut`); the point is
 *     keyed by the edge, so every tetrahedron around the edge sees the same one;
 *  2. *warps*: a cut point closer to an endpoint than `α` times the edge
 *     (`α_short` on an edge between a lattice centre and another vertex,
 *     `α_long` otherwise) violates that endpoint, and each violated vertex, in
 *     ascending order, moves onto its nearest violating cut point and takes the
 *     sign 0, which removes the cut points on its edges; a warp that would
 *     invert a tetrahedron around the vertex is refused and counted;
 *  3. cuts every tetrahedron of the material being split by a stencil chosen
 *     by the counts of its negative, zero and positive vertices: a tetrahedron
 *     with no positive vertex is kept whole, one with no negative vertex goes
 *     to the far side, and the others are cut into at most three tetrahedra on
 *     each side. A quadrilateral piece of a mesh face is split by the diagonal
 *     through the lower-numbered of its two mesh vertices, which both
 *     tetrahedra sharing the face choose alike; a quadrilateral on the
 *     interface by its shorter diagonal (ties by cut-point number), which both
 *     sides choose alike.
 *
 * The near side's pieces take material `Material + 1`, the far side's keep
 * `Material`; tetrahedra of other materials pass through unchanged (their
 * vertices can still be warped, under the same inversion guard, unless they
 * are pinned). Running it once per interface of a nested family of
 * interfaces gives a conforming multi-material mesh. Everything is serial or
 * writes disjoint slots in a fixed order: the result does not depend on the
 * thread count.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

/// A tetrahedral mesh with a material per tetrahedron.
struct StuffMesh {
    std::vector<Vec3> mPoints;
    std::vector<std::array<std::int64_t, 4>> mTets;
    /// Per tetrahedron, its material.
    std::vector<std::int32_t> mMaterial;
    /// Per tetrahedron, the lattice tetrahedron it was cut from.
    std::vector<std::int64_t> mParent;
    /// Per tetrahedron, 1 when stuffing cut it (it is no longer a lattice shape).
    std::vector<std::uint8_t> mCut;
    /// Per point, 1 for a lattice centre (the short-edge rule).
    std::vector<std::uint8_t> mCentre;
    /// Per point, 1 when it may not be warped (it lies on an earlier interface).
    std::vector<std::uint8_t> mPinned;
};

/// Where the interface crosses the edge from point @p A (sign -1) to @p B
/// (sign +1).
using StuffCut = std::function<Vec3(std::int64_t A, std::int64_t B)>;

/// Options for `isosurface_stuff`.
struct StuffOptions {
    double mAlphaLong = 0.24799;
    double mAlphaShort = 0.40173;
    /// The material being split.
    std::int32_t mMaterial = 0;
    /// Refuse warps that would leave a tetrahedron with a volume below this.
    double mMinVolume = 0.0;
};

/// What `isosurface_stuff` did.
struct StuffReport {
    std::int64_t mNumCutPoints = 0;    ///< cut points that survived warping
    std::int64_t mNumWarped = 0;       ///< vertices moved onto the interface
    std::int64_t mNumWarpRefused = 0;  ///< warps refused by the inversion guard
    std::int64_t mNumDegenerate = 0;   ///< output tetrahedra with a volume below the floor
};

/**
 * @brief Stuff one interface through @p rMesh, in place: warped points move,
 * cut points are appended, and the tetrahedra of `rOptions.mMaterial` are
 * replaced by their pieces. @p rSign is per point (resized for the cut points,
 * which are 0); points that end on the interface are pinned.
 */
StuffReport isosurface_stuff(StuffMesh& rMesh, std::vector<std::int8_t>& rSign,
                             const StuffCut& rCut, const StuffOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus
