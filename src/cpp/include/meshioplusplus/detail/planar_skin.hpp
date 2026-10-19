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
 * @file detail/planar_skin.hpp
 * @brief The 2-D skin shared by `voxelize` (pixels) and `octree` (quadtrees): a
 * planar polyline, taken from `line` cells or from the boundary edges of
 * `triangle`/`quad` cells, with the segment/box overlap and signed-distance
 * kernels that go with it.
 *
 * The skin is read in the xy plane; z is carried but never used by the kernels.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// Which coordinate plane a planar skin lies in.
enum class PlanarPlane {
    /// The axis along which the points have zero extent (z, then y, then x), else z.
    Auto = 0,
    XY = 1,
    XZ = 2,
    YZ = 3,
};

/// Parse a plane name (`auto`, `xy`, `xz`, `yz`).
MESHIOPLUSPLUS_API PlanarPlane planar_plane_from_name(const std::string& rName);

/// An axis permutation that maps a planar skin's plane to xy: frame axis `k`
/// is world axis `mAxis[k]`, so the in-plane axes come first and the normal
/// last. `Permute`/`Unpermute` convert between the two.
struct PlanarFrame {
    std::array<int, 3> mAxis{{0, 1, 2}};
    Vec3 Permute(const Vec3& rWorld) const {
        return Vec3{{rWorld[static_cast<std::size_t>(mAxis[0])],
                     rWorld[static_cast<std::size_t>(mAxis[1])],
                     rWorld[static_cast<std::size_t>(mAxis[2])]}};
    }
    Vec3 Unpermute(const Vec3& rFrame) const {
        Vec3 out{};
        for (std::size_t k = 0; k < 3; ++k)
            out[static_cast<std::size_t>(mAxis[k])] = rFrame[k];
        return out;
    }
};

/// The frame for @p Plane, deciding `Auto` from @p rMesh's point extents.
MESHIOPLUSPLUS_API PlanarFrame planar_frame(const Mesh& rMesh, PlanarPlane Plane);

/// A planar skin: segments, and the mesh point ids of their endpoints.
struct PlanarSkin {
    std::vector<std::array<Vec3, 2>> mSegments;
    std::vector<std::array<std::int64_t, 2>> mIds;
    /// True when at least one segment came from a `line` cell (such a polyline
    /// can be open; the boundary of a surface cannot).
    bool mHasLines = false;
};

/**
 * @brief Collect the segments of @p rMesh: every `line` cell, and the edges used
 * exactly once by `triangle`/`quad` cells (in first-occurrence order, after the
 * lines). Other blocks are ignored.
 * @param rPrefix error-message prefix, e.g. `"meshio++: octree: "`.
 * @param rFrame permutation applied to every segment point (identity: xy).
 * @throws std::invalid_argument when a cell references a missing point.
 */
MESHIOPLUSPLUS_API PlanarSkin collect_planar_skin(const Mesh& rMesh, const char* pPrefix,
                                                  const PlanarFrame& rFrame = {});

/// True when the mesh has `line` cells and nothing but `line`/`vertex` blocks.
MESHIOPLUSPLUS_API bool mesh_is_polyline_only(const Mesh& rMesh);

/// Throw, naming up to eight open ends, when an endpoint is used an odd number
/// of times (the polyline does not close, so it has no inside).
MESHIOPLUSPLUS_API void require_closed_polyline(const PlanarSkin& rSkin, const char* pPrefix);

/// Segment / box overlap in the xy plane by the separating axis test. The box is
/// centred at @p rCentre with half-extents @p rHalf; touching counts.
MESHIOPLUSPLUS_API bool segment_box_overlap_2d(const Vec3& rCentre, const Vec3& rHalf,
                                               const Vec3& rA, const Vec3& rB);

/// The point of the polyline nearest to @p rP, in the skin's xy plane (z is
/// copied from @p rP). An empty skin returns @p rP.
MESHIOPLUSPLUS_API Vec3 polyline_closest_point_2d(const PlanarSkin& rSkin, const Vec3& rP);

/// Signed distance from @p rP to the polyline, negative inside (winding number
/// nonzero, or odd crossing count with @p EvenOdd).
MESHIOPLUSPLUS_API double polyline_signed_distance_2d(const PlanarSkin& rSkin, const Vec3& rP,
                                                      bool EvenOdd);

}  // namespace detail
}  // namespace meshioplusplus
