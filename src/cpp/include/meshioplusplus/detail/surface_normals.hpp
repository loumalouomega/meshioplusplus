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
 * @file detail/surface_normals.hpp
 * @brief Vertex normals of a triangle soup, with the option to split a vertex
 * into one normal group per smooth patch.
 *
 * Two consumers share this kernel and must not disagree: the signed-distance
 * pseudonormal table in `surface_distance.cpp`, which needs one angle-weighted
 * normal per welded vertex, and `compute_normals` / the glTF writer, which need
 * per-*corner* normals because a crease has two normals at one position.
 *
 * ### Why the split is a graph problem, not a per-vertex loop
 *
 * A vertex is split by grouping the triangle corners around it into smooth
 * fans. Two corners of the same vertex belong to one fan when the triangles
 * that own them are joined along an edge through that vertex. An edge joins
 * its two triangles only when the pair is a *proper* manifold pair -- exactly
 * two users, walked in opposite directions -- and the dihedral angle between
 * their face normals is within the split angle. Boundary edges, non-manifold
 * edges and wound-the-same-way pairs therefore always cut. Triangles fanned
 * from the same polygon are joined unconditionally, so a non-planar polygon
 * never splits along its own diagonals.
 *
 * The joins are resolved with a union-find whose root is always the smallest
 * corner index in the set, and the groups are numbered by ascending root, so
 * the partition and its numbering depend only on the soup and the angle -- not
 * on the order the edges happened to be visited in. The numpy twin
 * (`_normals.py`) reproduces the same numbering.
 *
 * Normals are accumulated by a **serial** scatter in ascending
 * (triangle, corner) order: the order in which unit normals are summed changes
 * the last bits, and the SDF sign test can flip on a last-bit change.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/operations/sdf.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief The angle triangle (a, b, c) subtends at corner @p rA, in radians.
 *
 * 0 when either edge is degenerate. Used as a positive weight on a unit normal.
 */
MESHIOPLUSPLUS_API double corner_angle(const Vec3& rA, const Vec3& rB, const Vec3& rC);

/// Per triangle, the unnormalized normal `cross(b - a, c - a)`. Parallel; each
/// entry is independent.
MESHIOPLUSPLUS_API std::vector<Vec3> soup_face_normals(const TriangleSoup& rSoup);

/**
 * @brief The weighted sum of incident unit face normals at every soup point.
 *
 * @param rSoup the soup; the result is indexed by the soup's own point ids.
 * @param rFaceNormal `soup_face_normals(rSoup)`.
 * @param Weight `Angle` weights each incident face by the corner angle, `Area`
 *        by `|cross|`. Degenerate triangles contribute nothing.
 * @return one unnormalized sum per entry of `rSoup.mPoints`.
 */
MESHIOPLUSPLUS_API std::vector<Vec3> accumulate_vertex_normals(const TriangleSoup& rSoup,
                                                               const std::vector<Vec3>& rFaceNormal,
                                                               SdfPseudonormalWeight Weight);

/// The corner groups of a soup and one unit normal per group.
struct VertexNormalGroups {
    /// Per triangle corner (`3 * triangle + corner`), the group it belongs to.
    std::vector<std::int64_t> mCornerGroup;
    /// Per group, the soup point id every corner of the group sits on.
    std::vector<std::int64_t> mGroupPoint;
    /// Per group, its smallest corner index. Groups are numbered by this.
    std::vector<std::int64_t> mGroupRoot;
    /// Per group, the unit normal, or {0, 0, 0} when the incident faces sum to
    /// nothing (only degenerate triangles, or an exactly cancelling pair).
    std::vector<Vec3> mGroupNormal;
    /// Triangles whose normal has no direction.
    std::int64_t mNumDegenerate = 0;

    std::size_t NumGroups() const { return mGroupPoint.size(); }
};

/**
 * @brief Group the corners of @p rSoup into smooth fans and compute a unit
 * normal for each.
 *
 * @param SplitAngleDeg an edge joins its two triangles only when the angle
 *        between their face normals does not exceed this. A negative value
 *        disables splitting: every corner of a vertex is one group, so the
 *        result is one normal per touched point. Values of 180 and above keep
 *        every proper manifold pair joined, so only boundary, non-manifold and
 *        inconsistently wound edges still cut.
 */
MESHIOPLUSPLUS_API VertexNormalGroups vertex_normal_groups(const TriangleSoup& rSoup,
                                                           SdfPseudonormalWeight Weight,
                                                           double SplitAngleDeg);

}  // namespace detail
}  // namespace meshioplusplus
