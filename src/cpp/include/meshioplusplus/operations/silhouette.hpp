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
 * @file operations/silhouette.hpp
 * @brief The outline of a part seen along a direction (roadmap 7.6.3): its
 * surface projected onto a plane, unioned, as closed chains with its holes,
 * optionally offset by a clearance or made a prism (a nest, a jig, a foam
 * insert).
 *
 * The surface's triangles (quads and polygons fanned, a volume mesh's skin)
 * are projected onto the plane through the origin normal to `mNormal`, in
 * the frame `lsq_frame` completes it to, and their union is found exactly
 * (`detail/planar_union.hpp`): the candidate edges are the surface's
 * boundary edges and its contour generators, split where they cross, and
 * the region's boundary is walked into loops, outer boundaries
 * counter-clockwise about the normal and holes clockwise.
 *
 * `mClearance` > 0 grows the outline by that distance (a loose fit) and < 0
 * shrinks it (a press fit): the Minkowski sum or difference with a disk,
 * through the same union, of the region and a strip along its boundary. A
 * corner the offset opens is rounded with an arc (polygonal,
 * `mArcSegments` a full turn, its sides tangent to the true arc so the
 * offset is never less than the clearance) when it turns by less than
 * `mRoundAngleDeg`, and kept sharp (a mitre, bevelled beyond four times the
 * clearance) otherwise.
 *
 * The outline is a chained line mesh in the plane at the part's lowest
 * point along the normal: the Int64 cell data `section:chain` and
 * `section:closed` and the Float64 `section:offset` that `fit_curves`
 * reads, and `silhouette:hole` (1 for a hole). With `mExtrude`, the result is
 * instead the closed triangle surface of the prism over the outline from the
 * part's lowest to its highest point along the normal, outward wound. Both
 * carry the field data `silhouette:normal` and `silhouette:extent` (the
 * lowest and highest heights).
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The Int64 cell data marking a silhouette's holes (1) and outer loops (0).
inline constexpr const char* kSilhouetteHoleName = "silhouette:hole";
/// The Float64 field data: the unit viewing direction.
inline constexpr const char* kSilhouetteNormalName = "silhouette:normal";
/// The Float64 field data: the lowest and highest heights along the normal.
inline constexpr const char* kSilhouetteExtentName = "silhouette:extent";

struct SilhouetteOptions {
    std::array<double, 3> mNormal{0.0, 0.0, 1.0};  ///< the viewing direction
    double mClearance = 0.0;                       ///< > 0 grows the outline, < 0 shrinks it
    double mRoundAngleDeg = 180.0;                 ///< corners turning less than this are rounded
    std::int64_t mArcSegments = 64;                ///< polygon sides of a full turn of arc
    bool mExtrude = false;                         ///< a prism instead of the outline
    std::string mRegion;                           ///< restrict to this cell region
    /// Points closer than this times the projected box's diagonal are one.
    double mMergeTolerance = 1e-9;
};

struct SilhouetteResult {
    Mesh mMesh;
    std::int64_t mNumLoops = 0;
    std::int64_t mNumOuter = 0;
    std::int64_t mNumHoles = 0;
    std::int64_t mNumTriangles = 0;  ///< projected triangles of non-zero area
    double mArea = 0.0;              ///< of the outline's region
    double mPerimeter = 0.0;         ///< the loops' total length
    double mMinHeight = 0.0;         ///< along the normal
    double mMaxHeight = 0.0;
    std::vector<std::string> mWarnings;
};

/**
 * @brief The outline of @p rMesh along `rOptions.mNormal` (see the file comment).
 * @throws std::invalid_argument on a zero normal, a non-finite clearance or
 *         angle, fewer than 4 arc segments, a line or point mesh, an unknown
 *         region, or coordinates the exact predicates cannot hold.
 */
MESHIOPLUSPLUS_API SilhouetteResult silhouette(const Mesh& rMesh,
                                               const SilhouetteOptions& rOptions = {});

}  // namespace meshioplusplus
