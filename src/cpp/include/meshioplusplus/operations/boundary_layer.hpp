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
 * @file boundary_layer.hpp
 * @brief Prismatic boundary layers grown from a surface (roadmap 6.9; the
 * design is `doc/boundary_layer_plan.md`).
 *
 * From a triangle and quad surface, `mNumLayers` layers of cells grow along
 * the surface's normals: wedges from triangles, hexahedra from quads, the
 * first `mFirstHeight` thick and each next `mGrowth` times the one before.
 *
 *  - **Normals.** Each point grows along its *most normal* normal, the unit
 *    vector `n` maximising the least `n . n_f` over its faces (A. Aubry and
 *    R. Löhner, "On the 'most normal' normal", Commun. Numer. Meth. Engng.
 *    25, 2009), found from the angle-weighted mean by ascent on the face
 *    that sees it least; its offsets are divided by that least cosine, so
 *    every face's layers keep their thickness at a ridge or a corner. A point
 *    no direction sees all its faces from (a saddle folded past 90 degrees)
 *    does not grow.
 *  - **Collisions.** A ray from each point along its direction against the
 *    surface: the point's layers take at most `mCollisionFraction` of the
 *    distance to what it hits, so layers growing toward each other from both
 *    sides of a gap stop before they meet. The fraction a point keeps is
 *    smoothed `mSmoothIterations` times over its neighbours, never raised,
 *    so the layers taper rather than step (the advancing-layer practice of
 *    S. Pirzadeh, AIAA J. 34, 1996, and R. V. Garimella and M. S. Shephard,
 *    Int. J. Numer. Meth. Engng. 49, 2000).
 *  - **Validity.** Every cell is judged by its least corner Jacobian; the
 *    points of a cell that folds have their layers halved, at most
 *    `mMaxShrink` times, and then collapsed: a point whose layers collapse
 *    stays on the surface, and its cells become pyramids and tetrahedra
 *    (`detail/layer_sweep.hpp`). With `mCollapse` = `Layers` a point instead
 *    drops its top layer, so it grows fewer layers at their full thickness
 *    and its cells above end in pyramids and tetrahedra (Garimella and
 *    Shephard 2000). No cell is written folded.
 *  - **Fronts.** Once every cell is valid, each face of the last layer is
 *    tested against every other one and the surface (`mFrontCheck`), the
 *    advancing-layer front test of Pirzadeh (1996), with the exact
 *    triangle–triangle kernel of `detail/tri_tri.hpp` over a bounding-volume
 *    hierarchy; faces whose columns share a surface point are neighbours and
 *    are not tested. The points of a face that crosses or touches another
 *    give way as a folding cell's do, and the layers are judged again. A ray
 *    misses what lies between its point's neighbours; the front test does
 *    not.
 *  - **Smoothing the normals.** `mNormalSmoothIterations` Laplacian passes
 *    move each direction toward its neighbours' mean, a move kept only while
 *    every face around the point still sees it at more than half the least
 *    cosine its most normal normal had; the stretch follows the new least
 *    cosine. On a coarse surface it turns the directions less sharply from
 *    point to point.
 *
 * The result holds the cells, the surface (`boundary_layer:wall`) and the
 * last layer (`boundary_layer:cap`) as faces, named Cell regions and the
 * cell data `boundary_layer:layer` (1..L for the cells, 0 the wall, L + 1
 * the cap); every cell data array of the surface's faces (`cad:face` among
 * them) is carried to the cells and faces grown from them. Point data is
 * carried to every copy of a point. The volume the cap closes off is the
 * caller's (`octree_mesh` fills it, fitted rather than conforming; roadmap
 * 6.7 has the conforming route). Deterministic: points and cells are
 * visited in index order.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

inline constexpr const char* kBoundaryLayerName = "boundary_layer:layer";
inline constexpr const char* kBoundaryLayerWall = "boundary_layer:wall";
inline constexpr const char* kBoundaryLayerCap = "boundary_layer:cap";
inline constexpr const char* kBoundaryLayerCells = "boundary_layer:cells";

/// How a point whose cells fold gives way (roadmap 6.9.1).
enum class BoundaryLayerCollapse : std::int32_t {
    Points = 0,  ///< its layers halved, then the point collapsed
    Layers = 1,  ///< its top layer dropped: fewer layers at the point
};

struct BoundaryLayerOptions {
    double mFirstHeight = 0.0;         ///< the first layer's thickness (required, > 0)
    double mGrowth = 1.2;              ///< each layer over the one before (>= 1)
    std::int32_t mNumLayers = 3;       ///< 1 to 1000
    bool mReverse = false;             ///< grow against the faces' normals
    double mCollisionFraction = 0.45;  ///< of the distance a point's ray hits, (0, 1]
    std::int32_t mSmoothIterations = 5;
    std::int32_t mMaxShrink = 6;  ///< halvings before a point collapses
    BoundaryLayerCollapse mCollapse = BoundaryLayerCollapse::Points;
    /// Laplacian passes over the directions (0 to 1000).
    std::int32_t mNormalSmoothIterations = 0;
    /// Test the last layer's faces against each other and the surface.
    bool mFrontCheck = true;
};

struct BoundaryLayerResult {
    Mesh mMesh;
    std::int64_t mNumCells = 0;  ///< volume cells
    std::int64_t mNumWedges = 0;
    std::int64_t mNumHexahedra = 0;
    std::int64_t mNumPyramids = 0;
    std::int64_t mNumTetrahedra = 0;
    std::int64_t mNumShrunk = 0;      ///< points whose layers are thinner than asked
    std::int64_t mNumCollapsed = 0;   ///< of those, points that do not grow at all
    std::int64_t mNumNotVisible = 0;  ///< of those, points no direction sees all faces from
    /// With `Layers`, the points that grow fewer layers than asked (but some).
    std::int64_t mNumFewerLayers = 0;
    /// Pairs of last-layer faces found crossing or touching another face,
    /// over every round of the front test.
    std::int64_t mNumFrontCrossings = 0;
    double mMinScaledJacobian = 0.0;  ///< least corner Jacobian over the cells
    double mMinThickness = 0.0;       ///< the thinnest and thickest total, over points
    double mMaxThickness = 0.0;
    /// Whether every cap edge is shared by two cap faces (a closed surface
    /// gives a closed cap).
    bool mCapClosed = false;
    std::vector<std::string> mWarnings;
};

/**
 * @brief Grow boundary layers from the triangle and quad surface `rSurface`
 * (see the file comment).
 * @throws std::invalid_argument for an option out of range (`mCollapse`
 * among them), a surface with
 * no triangle or quad, a cell block that is neither (lines and vertices are
 * skipped), an edge used by more than two faces, or two faces that disagree
 * on its direction.
 */
MESHIOPLUSPLUS_API BoundaryLayerResult boundary_layer(const Mesh& rSurface,
                                                      const BoundaryLayerOptions& rOptions);

}  // namespace meshioplusplus
