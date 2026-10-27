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

#pragma once

/**
 * @file detail/hex_fit.hpp
 * @brief Core-private pieces of an all-hexahedral mesh fitted to a skin:
 * manifold repair of a set of grid cells, the boundary of a hexahedral mesh,
 * pillowing, the corner scaled Jacobian and an untangler.
 *
 * The grid cells kept by a voxel test need not have a manifold boundary: two
 * cells may touch along an edge or at a corner only. `hex_grid_make_manifold`
 * adds cells until, around every grid vertex, the kept and the empty cells of
 * its 2x2x2 block are each face-connected, which makes the boundary a closed
 * 2-manifold (the condition the "manifold voxelisation" of Sculpt-style
 * hexahedral meshers enforces). A pillow (Mitchell and Tautges 1995) then
 * puts one layer of hexahedra along the boundary, so no hexahedron has two
 * faces on the skin. The untangler is a Gauss-Seidel pattern search on each
 * node's worst corner scaled Jacobian, in ascending node order, so it is
 * deterministic.
 */

#include <array>
#include <cstdint>
#include <functional>
#include <vector>

#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Add cells to @p rKept (an `Nx x Ny x Nz` grid, cell `(k * Ny + j) *
 * Nx + i`) until every vertex's 2x2x2 block has face-connected kept and empty
 * cells. Each failing block takes the empty cell with the largest @p rScore
 * first (ties: the lowest index); cells outside the grid count as empty and
 * are never added.
 * @return the cells added.
 */
std::int64_t hex_grid_make_manifold(const std::array<std::int64_t, 3>& rN,
                                    std::vector<std::uint8_t>& rKept,
                                    const std::vector<double>& rScore);

/// A hexahedral mesh, VTK node order.
struct HexMesh {
    std::vector<Vec3> mPoints;
    std::vector<std::array<std::int64_t, 8>> mHexes;
};

/**
 * @brief `hex_grid_make_manifold` on an unstructured hexahedral mesh: add
 * hexahedra to @p rKept until, around every node, the kept hexahedra and the
 * rest (with the outside of the mesh) are each connected through the faces at
 * that node. A failing node takes the empty hexahedron around it with the
 * largest @p rScore first (ties: the lowest index); sweeps the nodes in
 * ascending order until none fails.
 * @return the hexahedra added.
 */
std::int64_t hex_make_manifold(const HexMesh& rMesh, std::vector<std::uint8_t>& rKept,
                               const std::vector<double>& rScore);

/// The faces of a hexahedral mesh used by one hexahedron, ordered so their
/// normal (right-hand rule) points out of it, with the (hexahedron, facet)
/// each comes from (facets in the `cell_faces` order of `hexahedron`).
struct HexBoundary {
    std::vector<std::array<std::int64_t, 4>> mQuads;
    std::vector<std::array<std::int64_t, 2>> mSides;
    /// Edges used by one boundary face, and by more than two.
    std::int64_t mNumOpenEdges = 0;
    std::int64_t mNumNonManifoldEdges = 0;
};

/// The boundary of @p rMesh.
HexBoundary hex_boundary(const HexMesh& rMesh);

/**
 * @brief Pillow @p rMesh: every boundary vertex gets a copy (appended to the
 * points, at the same place), and every boundary face `(a, b, c, d)` a new
 * hexahedron `(a, b, c, d, a', b', c', d')` (appended to the hexahedra), so
 * the copies are the new boundary.
 * @return per point of the pillowed mesh, its copy, or -1 (copies and
 *         interior points).
 */
std::vector<std::int64_t> hex_pillow(HexMesh& rMesh, const HexBoundary& rBoundary);

/// The smallest scaled Jacobian of a hexahedron's eight corners (-1 for a
/// corner with a zero-length edge).
double hex_min_scaled_jacobian(const std::vector<Vec3>& rPoints,
                               const std::array<std::int64_t, 8>& rHex);

/// How a node may move in `hex_untangle`.
enum class HexNodeMode : std::uint8_t {
    Free = 0,   ///< anywhere
    Slide = 1,  ///< along the skin (every candidate is projected)
    Fixed = 2,  ///< never
};

/// `hex_untangle`'s options.
struct HexUntangleOptions {
    std::int64_t mIterations = 50;  ///< sweeps at most
    double mTarget = 0.6;           ///< nodes whose worst corner reaches this are left
};

/// What `hex_untangle` did.
struct HexUntangleReport {
    std::int64_t mMoves = 0;
    std::int64_t mSweeps = 0;
};

/**
 * @brief Move the nodes of @p rHexes to raise the worst corner scaled
 * Jacobian around each: per sweep, every movable node below the target tries
 * the mean of its edge neighbours and six axis steps (a fifth of its mean
 * edge, halving every ten sweeps), projected by @p rProject when it slides,
 * and takes the best candidate if it improves.
 */
HexUntangleReport hex_untangle(std::vector<Vec3>& rPoints,
                               const std::vector<std::array<std::int64_t, 8>>& rHexes,
                               const std::vector<HexNodeMode>& rMode,
                               const std::function<Vec3(const Vec3&)>& rProject,
                               const HexUntangleOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus
