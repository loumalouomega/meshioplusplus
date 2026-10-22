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
//  Provenance:      our own construction, after the body-centred cubic
//                   lattice of Labelle and Shewchuk, "Isosurface Stuffing:
//                   Fast Tetrahedral Meshes with Good Dihedral Angles"
//                   (SIGGRAPH 2007), extended to a graded octree; the angle
//                   table is derived here and pinned by the tests.
//
#pragma once

/**
 * @file detail/bcc_octree_tets.hpp
 * @brief A conforming, graded body-centred cubic tetrahedralisation of a
 * 2:1-balanced (faces, edges and corners) linear octree, for `octree_mesh`.
 *
 * Coordinates are integers on the **doubled** finest lattice, so that every
 * vertex is exact: a leaf with lo corner `X` and span `s` has its corners at
 * `2X + {0, 2s}` and its centre at `2X + s`.
 *
 * **Construction.** Every leaf is a cone, from its centre, over fans of its
 * face portions:
 *
 *  - a face on the domain boundary, or one the leaf shares with a coarser
 *    neighbour, is one portion; a face shared with four finer neighbours is
 *    four portions, their faces;
 *  - each portion is fanned from its own centre (a *face-centre* vertex), its
 *    boundary edges cut at their midpoints where a leaf corner lies there
 *    (2:1 balance across edges means at most one cut per edge);
 *  - a face shared with a neighbour of the same size is not fanned: the two
 *    cones over it form a bipyramid around the segment between the two
 *    centres, which is tetrahedralised around that axis instead.
 *
 * A portion's fan depends only on the portion and on which edge midpoints are
 * vertices, which both sides of a face see alike, so the tetrahedra conform
 * across every face; a leaf is convex and its centre sees its whole boundary,
 * so each cone fills its leaf. Where neighbouring leaves are the same size the
 * tetrahedra are exactly the BCC lattice's. There are six shapes:
 *
 * | class | where | dihedral angles (degrees) |
 * |---|---|---|
 * | T1 | equal face, whole edge | 60, 60, 60, 60, 90, 90 |
 * | T2 | equal face, half edge | 45, 60, 60, 90, 90, 90 |
 * | T3 | boundary or finer side, whole edge | 45, 60, 60, 90, 90, 90 |
 * | T4 | boundary or finer side, half edge | 45, 45, 60, 90, 90, 90 |
 * | T5 | coarser side, half of the leaf's edge | 35.26, 45, 60, 90, 90, 125.26 |
 * | T6 | coarser side, from an edge midpoint to the face centre | 45, 54.74, 54.74, 90, 90, 90 |
 *
 * so every tetrahedron's dihedral angles lie in [35.26°, 125.26°], and a
 * region of equal leaves is T1 only. Vertices are numbered by the Morton key
 * of their doubled coordinates and tetrahedra come out leaf by leaf in key
 * order, so the result is a function of the tree alone. Every tetrahedron is
 * positively oriented (`(b - a) · ((c - a) × (d - a)) > 0`).
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "linear_octree.hpp"
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

/// What a lattice vertex is.
enum class BccVertexType : std::uint8_t {
    Corner = 0,      ///< a leaf corner
    Centre = 1,      ///< a leaf centre
    FaceCentre = 2,  ///< the centre of a fanned face portion
};

/// The vertices and tetrahedra of the lattice.
struct BccLattice {
    /// Per vertex, the Morton key of its doubled coordinates, ascending.
    std::vector<std::uint64_t> mKey;
    /// Per vertex, its doubled coordinates.
    std::vector<std::array<std::int64_t, 3>> mIjk2;
    /// Per vertex, what it is.
    std::vector<BccVertexType> mType;
    /// The tetrahedra, as vertex indices, positively oriented.
    std::vector<std::array<std::int64_t, 4>> mTets;
    /// Per tetrahedron, its shape class (1-6, as in the table above).
    std::vector<std::uint8_t> mShape;
    /// Per tetrahedron, the leaf whose cone (or, across an equal face, whose
    /// lower side) it belongs to.
    std::vector<std::int64_t> mLeaf;
};

/**
 * @brief Tetrahedralise @p rTree, which must be 2:1-balanced across faces,
 * edges and corners and have `Extent(k) <= 2^20` on every axis.
 * @throws std::invalid_argument when it would have more than @p MaxTets
 *         tetrahedra (0 disables the check), and std::logic_error when the
 *         tree is not balanced.
 */
BccLattice bcc_octree_lattice(const LinearOctree& rTree, std::int64_t MaxTets);

/// The six dihedral angles of a tetrahedron, in degrees, ascending.
std::array<double, 6> tet_dihedral_angles(const Vec3& rA, const Vec3& rB, const Vec3& rC,
                                          const Vec3& rD);

}  // namespace detail
}  // namespace meshioplusplus
