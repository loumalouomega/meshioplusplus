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
 * @file operations/octree.hpp
 * @brief Octree (3-D) and quadtree (2-D) meshes from a skin: leaves are fine
 * only where the surface is.
 *
 * The output is an ordinary `Mesh` of `hexahedron` (3-D) or `quad` (2-D) leaves,
 * with hanging nodes where a fine leaf meets a coarse one, flagged in the
 * `octree:hanging` point data. See doc/octree.md for the contract.
 *
 * Every point is an integer coordinate on the finest level's lattice, so shared
 * corners are deduplicated arithmetically -- no hash, no tolerance -- and
 * coordinates are `origin + spacing * index`, never accumulated. Leaves come out
 * in Morton order, so the result is byte-identical across backends and thread
 * counts, and across the C++ core and the numpy twin.
 *
 * @see doc/octree.md, operations/voxelize.hpp
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

namespace meshioplusplus {

/// Per-leaf Int64 level (0 = a root cell) cell data.
inline constexpr const char* kOctreeLevelName = "octree:level";
/// Per-leaf Int64 0/1 cell data: 1 when the skin passes through the leaf.
inline constexpr const char* kOctreeCutName = "octree:cut";
/// Per-point Int64 0/1 point data: 1 for a node on the interior of a face or
/// edge of a leaf rather than at its corners.
inline constexpr const char* kOctreeHangingName = "octree:hanging";

/// How strictly neighbouring leaves must agree in size.
enum class OctreeBalance {
    /// No balancing: neighbours may differ by any number of levels.
    None = 0,
    /// 2:1 across faces (edges in 2-D).
    Face = 1,
    /// 2:1 across faces, edges and corners.
    Full = 2,
};

/// Parse an `OctreeBalance` name (`none`, `face`, `full`).
MESHIOPLUSPLUS_API OctreeBalance octree_balance_from_name(const std::string& rName);

/// Options for `octree`.
struct OctreeOptions {
    /// Root cells per axis (a forest of cubic roots). The third entry is ignored
    /// in 2-D.
    std::array<std::int64_t, 3> mRootCells{{1, 1, 1}};
    /// Every leaf is split down to at least this level.
    std::int64_t mMinLevel = 0;
    /// No leaf is split beyond this level (at most 20).
    std::int64_t mMaxLevel = 4;
    /// Explicit bounds `{lo[3], hi[3]}`; unset uses the skin's bounding box.
    std::optional<std::array<double, 6>> mBounds;
    /// Padding added to every side, in world units.
    double mPadding = 0.0;
    /// Padding added to every side, as a fraction of the bounding-box diagonal.
    double mPaddingRelative = 0.0;
    /// 2 for a planar skin of `line` cells (or the boundary of a 2-D surface
    /// mesh) in the plane of its first point's z; 3 for a surface; 0 decides
    /// from the cell types (only `line`/`vertex` blocks mean 2).
    std::int64_t mDimension = 0;
    /// Refine further where the skin is within this distance of a leaf (0: only
    /// where the skin cuts it).
    double mSampleDistance = 0.0;
    /// Neighbour size agreement.
    OctreeBalance mBalance = OctreeBalance::Face;
    /// Which leaves to keep, as for `voxelize`.
    VoxelFill mFill = VoxelFill::All;
    /// Attach `sdf:distance` (signed distance at leaf centres).
    bool mAttachDistance = false;
    /// 2-D only: use the even-odd rule instead of the winding number.
    bool mEvenOdd = false;
    /// Refuse to generate more leaves than this.
    std::int64_t mMaxCells = 20000000;
    /// How 3-D `VoxelFill::Inside` decides what is inside.
    SurfaceDistanceOptions mDistance;
};

/// The result of `octree`.
struct OctreeResult {
    Mesh mMesh;                                       ///< The leaves.
    std::int64_t mDimension = 3;                      ///< 2 or 3.
    std::array<double, 3> mOrigin{{0.0, 0.0, 0.0}};   ///< Lo corner of the root lattice.
    double mRootSize = 0.0;                           ///< Edge length of a root cell.
    std::int64_t mNumLeaves = 0;                      ///< Leaves kept.
    std::int64_t mNumHanging = 0;                     ///< Hanging nodes in the output.
    std::int64_t mFinestLevel = 0;                    ///< Deepest level present.
};

/**
 * @brief Build an octree (quadtree) around the skin @p rMesh.
 * @throws std::invalid_argument on a non-positive root count, a level range
 *         that is empty or above 20, a lattice too fine for 63-bit keys,
 *         inverted bounds, more leaves than `mMaxCells`, and (for an `Inside`
 *         fill in 2-D) a polyline that does not close, naming its open ends.
 */
MESHIOPLUSPLUS_API OctreeResult octree(const Mesh& rMesh, const OctreeOptions& rOptions = {});

}  // namespace meshioplusplus
