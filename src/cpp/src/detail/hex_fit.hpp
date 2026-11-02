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
#include <cstddef>
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
 * largest @p rScore first (ties: the lowest index), among @p pAllowed when
 * given; sweeps the nodes in ascending order until none fails.
 * @return the hexahedra added.
 */
std::int64_t hex_make_manifold(const HexMesh& rMesh, std::vector<std::uint8_t>& rKept,
                               const std::vector<double>& rScore,
                               const std::vector<std::uint8_t>* pAllowed = nullptr);

/**
 * @brief The converse of `hex_make_manifold`: remove kept hexahedra (the
 * smallest @p rScore first, ties: the lowest index) until every node passes.
 * Removing every kept hexahedron around a node always makes it pass, so this
 * terminates. @return the hexahedra removed.
 */
std::int64_t hex_shrink_manifold(const HexMesh& rMesh, std::vector<std::uint8_t>& rKept,
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

/// What `hex_pillow_interface` did.
struct HexInterface {
    /// Per point of the mesh after it, the shell node and the inner copy made
    /// for each interface node (-1 elsewhere).
    std::vector<std::int64_t> mShell, mInner;
    /// The interface faces on their original (outer) nodes, each pointing out
    /// of its inner hexahedron, and that hexahedron.
    std::vector<std::array<std::int64_t, 4>> mQuads;
    std::vector<std::int64_t> mInnerHex;
    /// The first hexahedron added: two per face, the inner layer then the outer.
    std::size_t mFirstHex = 0;
};

/**
 * @brief Pillow both sides of the interface between the hexahedra in
 * @p rInner and the others: every interface node `v` gets a shell node and an
 * inner copy (appended), the inner hexahedra move onto the inner copies, and
 * every interface face gets two hexahedra, one from the inner copies to the
 * shell nodes and one from the shell nodes to the original nodes, so no
 * hexahedron has two faces on the interface. The inner set's boundary must be
 * a manifold of interior faces (`hex_make_manifold`), away from the mesh's
 * boundary.
 */
HexInterface hex_pillow_interface(HexMesh& rMesh, const std::vector<std::uint8_t>& rInner);

/// The smallest scaled Jacobian of a hexahedron's eight corners (-1 for a
/// corner with a zero-length edge).
double hex_min_scaled_jacobian(const std::vector<Vec3>& rPoints,
                               const std::array<std::int64_t, 8>& rHex);

/// The scaled Jacobian of each of a hexahedron's eight corners (-1 for a
/// corner with a zero-length edge).
std::array<double, 8> hex_corner_scaled_jacobians(const std::vector<Vec3>& rPoints,
                                                  const std::array<std::int64_t, 8>& rHex);

/**
 * @brief The unit outward normal of the boundary @p rQuads at each point:
 * the sum of the area vectors of the quadrilaterals around it, then @p Smooth
 * times replaced by the sum over itself and its boundary neighbours (in
 * ascending order, so the result is deterministic). Zero for a point on no
 * quadrilateral, or where the sum vanishes.
 */
std::vector<Vec3> hex_boundary_normals(const std::vector<Vec3>& rPoints,
                                       const std::vector<std::array<std::int64_t, 4>>& rQuads,
                                       int Smooth);

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
    /// Also try each node's trilinear reflection, 26 directions and a line
    /// search for a node below 0.2, and restart the step once when a sweep
    /// stalls with a corner below 0.2. False is the plain pattern search.
    bool mRobust = false;
};

/// What `hex_untangle` did.
struct HexUntangleReport {
    std::int64_t mMoves = 0;
    std::int64_t mSweeps = 0;
};

/// Projects a candidate position of node `Node` (a sliding one) onto the skin.
using HexProjector = std::function<Vec3(std::size_t Node, const Vec3& rP)>;

/**
 * @brief Move the nodes of @p rHexes to raise the worst corner scaled
 * Jacobian around each: per sweep, every movable node below the target tries
 * the mean of its edge neighbours and six axis steps (a fifth of its mean
 * edge, halving every ten sweeps), projected by @p rProject when it slides,
 * and takes the best candidate if it improves. With `mRobust`, a node also
 * tries the mean of its trilinear reflections (the position that would make
 * each corner of it a parallelepiped's), a node below 0.2 steps along the 26
 * directions of a cube and doubles or halves its best step, and a sweep that
 * moves nothing while a corner is below 0.2 restarts the step size once.
 */
HexUntangleReport hex_untangle(std::vector<Vec3>& rPoints,
                               const std::vector<std::array<std::int64_t, 8>>& rHexes,
                               const std::vector<HexNodeMode>& rMode, const HexProjector& rProject,
                               const HexUntangleOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus
