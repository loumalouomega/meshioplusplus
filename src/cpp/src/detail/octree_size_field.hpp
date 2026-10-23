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
 * @file detail/octree_size_field.hpp
 * @brief The skin of `octree_mesh` (welded, with normals and features), its
 * size field, and the adaptive linear octree that size field refines.
 *
 * The size field is evaluated per skin primitive, as the smallest size its
 * vertices ask for:
 *
 *     h = clamp(min(h_max, h_curv, h_feat, h_gap, h_user), h_min, h_max)
 *
 *  - `h_curv = 2 sin(θ/2) / κ_max`, the chord that subtends the angle θ on a
 *    circle of the largest principal curvature `κ_max = |H| + sqrt(H² - K)`;
 *  - `h_feat`, the size along sharp and non-manifold edges (segment primitives);
 *  - `h_gap = t / gap_cells`, `t` the distance an inward ray from the vertex
 *    travels before it leaves the skin through its far wall;
 *  - `h_user`, a `point_data` array.
 *
 * The root lattice has `L = ceil(log2(h_max / h_min))` levels and a root edge
 * `R = h_min 2^L`, so the finest leaf is `h_min` exactly and a leaf at level
 * `lev(h) = ceil(log2(R / h))` is no larger than `h`. The octree is built top
 * down as `octree` builds its own: each leaf keeps the primitives whose exact
 * triangle-box overlap (Akenine-Möller, JGT 6(1), 2001) hits it, and splits
 * while it is coarser than the finest level its primitives ask for. Leaves
 * come out sorted by Morton key, so the result does not depend on the order of
 * the work.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "linear_octree.hpp"
#include "tri_bvh.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// The skin as welded triangles.
struct OctreeSkin {
    /// The triangle soup (as `build_triangle_soup` reduces the input).
    TriangleSoup mSoup;
    /// Welded points: exactly coincident input points are one.
    std::vector<Vec3> mPoints;
    /// Per soup triangle, its welded vertices.
    std::vector<std::array<std::int64_t, 3>> mTris;
    /// Per welded point, the input point it was first seen as.
    std::vector<std::int64_t> mSourcePoint;
    /// Per welded point, the area-weighted unit normal of its triangles.
    std::vector<Vec3> mNormal;
    /// The bounding box of the triangles' corners.
    Vec3 mLo{0.0, 0.0, 0.0};
    Vec3 mHi{0.0, 0.0, 0.0};
    /// The boundary of a volume input (empty for a surface input).
    Mesh mExtracted;
    bool mFromVolume = false;
    /// The surface the soup was taken from: @p rInput itself, or its boundary.
    const Mesh& Surface(const Mesh& rInput) const { return mFromVolume ? mExtracted : rInput; }
};

/// Reduce @p rMesh (restricted to @p rRegion) to an `OctreeSkin`.
/// @throws std::invalid_argument when no surface triangle remains.
OctreeSkin octree_skin(const Mesh& rMesh, const std::string& rRegion, const char* pPrefix);

/// A sharp or non-manifold edge of the skin, between welded points.
struct OctreeFeatureEdge {
    std::int64_t mA = 0;
    std::int64_t mB = 0;
};

/// The sharp (dihedral above @p FeatureAngleDeg) and non-manifold edges of
/// the skin, ascending; empty when @p FeatureAngleDeg is not positive.
std::vector<OctreeFeatureEdge> octree_skin_features(const OctreeSkin& rSkin,
                                                    double FeatureAngleDeg);

/// The root lattice of the octree.
struct OctreeFrame {
    std::array<double, 3> mOrigin{{0.0, 0.0, 0.0}};
    double mRootSize = 0.0;
    std::array<std::int64_t, 3> mRoots{{1, 1, 1}};
    std::int32_t mMaxLevel = 0;
    double mMinSize = 0.0;  ///< the finest leaf edge, `mRootSize / 2^mMaxLevel`
    double mMaxSize = 0.0;  ///< the largest leaf edge allowed
    /// Edge of a finest cell.
    double Fine() const { return mRootSize / static_cast<double>(std::int64_t{1} << mMaxLevel); }
    /// The level whose leaves are no larger than @p Size.
    std::int32_t LevelFor(double Size) const;
    /// World coordinates of a finest-lattice point.
    Vec3 Position(std::int64_t X, std::int64_t Y, std::int64_t Z) const;
};

/**
 * @brief The root lattice over the box `[rLo, rHi]` (already padded): root
 * edge `MinSize * 2^L` with `L = ceil(log2(MaxSize / MinSize))` (capped at
 * @p MaxLevelCap when it is not negative), roots per axis covering the box,
 * centred on it.
 * @throws std::invalid_argument when the lattice needs `2^20` or more finest
 *         cells on an axis.
 */
OctreeFrame octree_frame(const Vec3& rLo, const Vec3& rHi, double MinSize, double MaxSize,
                         std::int64_t MaxLevelCap, const char* pPrefix);

/// A skin primitive with its target size: a triangle, or a feature segment
/// `(A, B, B)`.
struct SizedPrim {
    Vec3 mA{0.0, 0.0, 0.0};
    Vec3 mB{0.0, 0.0, 0.0};
    Vec3 mC{0.0, 0.0, 0.0};
    double mSize = 0.0;
};

/// What drives the size field.
struct SizeFieldSpec {
    double mMinSize = 0.0;
    double mMaxSize = 0.0;
    double mCurvatureAngleDeg = 0.0;  ///< 0: off
    double mFeatureSize = 0.0;        ///< 0: off (no feature primitives)
    double mGapCells = 0.0;           ///< 0: off
    bool mGapOutside = false;
    const double* pUserSize = nullptr;  ///< per input point, or null
};

/// The size each welded point asks for (before clamping), and the primitives.
struct SizeField {
    std::vector<double> mPointSize;
    std::vector<SizedPrim> mPrims;
    /// The smallest size any primitive asks for (after clamping).
    double mSkinMinSize = 0.0;
};

/// Evaluate the size field on the skin.
SizeField octree_size_field(const OctreeSkin& rSkin, const TriBvh& rTree,
                            const std::vector<OctreeFeatureEdge>& rFeatures,
                            const SizeFieldSpec& rSpec);

/// A refined octree with, per leaf, the primitives that reach it.
struct AdaptiveOctree {
    LinearOctree mTree;
    /// Per leaf, 1 when a primitive overlaps it.
    std::vector<std::uint8_t> mCut;
    /// Per leaf, the indices of the primitives that overlap it.
    std::vector<std::vector<std::int32_t>> mCand;
};

/**
 * @brief Refine the frame's roots top down until every leaf is at the level
 * its primitives ask for: `max(lev(size))` over the primitives that overlap
 * it (or @p SurfaceLevel for every cut leaf when it is not negative), and at
 * least `lev(MaxSize)` everywhere.
 */
AdaptiveOctree octree_adaptive_build(const OctreeFrame& rFrame,
                                     const std::vector<SizedPrim>& rPrims,
                                     std::int32_t SurfaceLevel, std::int64_t MaxCells);

/// After the tree was refined further (`rOrigin[i]` the old leaf new leaf `i`
/// came from), recompute the candidates and the cut flags of the new leaves.
void octree_adaptive_refresh(AdaptiveOctree& rOctree, const OctreeFrame& rFrame,
                             const std::vector<SizedPrim>& rPrims,
                             const std::vector<std::int64_t>& rOrigin,
                             std::vector<std::vector<std::int32_t>> OldCand);

/// Whether primitive @p rPrim overlaps the closed box of a leaf.
bool octree_prim_hits(const OctreeFrame& rFrame, const SizedPrim& rPrim,
                      const std::array<std::int64_t, 3>& rLo, std::int64_t Span);

}  // namespace detail
}  // namespace meshioplusplus
