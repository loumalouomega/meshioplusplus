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
 * @file operations/octree_mesh_internal.hpp
 * @brief What the stages of `octree_mesh` share: the prepared skin, its trees,
 * the frame and the balanced octree, and one entry point per output.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <chrono>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "../detail/linear_octree.hpp"
#include "../detail/octree_size_field.hpp"
#include "../detail/tri_bvh.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// The prefix of every `octree_mesh` message.
inline constexpr const char* kOctreeMeshPrefix = "meshio++: octree_mesh: ";

/// Everything the outputs are built from. Not movable once built: the trees
/// point into the skin's soup.
struct OctreeMeshSetup {
    OctreeSkin mSkin;
    TriBvh mBvh;
    DistanceQuery mQuery;
    OctreeFrame mFrame;
    std::vector<OctreeFeatureEdge> mFeatures;
    SizeField mField;
    /// The balanced octree, with per-leaf candidates and cut flags.
    AdaptiveOctree mOctree;
    /// The skin's bounding-box diagonal.
    double mDiagonal = 0.0;
    /// -1 when the skin is wound inside out (its winding numbers are negated).
    double mWindingSign = 1.0;
    /// Points at most this far from the skin are on it.
    double mOnSkin = 0.0;
};

/// A stopwatch that records each stage's seconds in the result.
struct OctreeMeshClock {
    std::chrono::steady_clock::time_point mStart = std::chrono::steady_clock::now();
    void Lap(OctreeMeshResult& rResult, const char* pStage);
};

/// Log a warning and keep it in the result.
void octree_mesh_warn(OctreeMeshResult& rResult, const std::string& rMessage);

/// The winding number at every point (the skin's orientation applied), with
/// the values near the threshold re-evaluated more accurately.
std::vector<double> octree_mesh_winding(const OctreeMeshSetup& rSetup,
                                        const std::vector<Vec3>& rPoints,
                                        const OctreeMeshOptions& rOptions,
                                        OctreeMeshResult& rResult);

/// The centre of every leaf of a tree.
std::vector<Vec3> octree_mesh_leaf_centres(const OctreeFrame& rFrame, const LinearOctree& rTree);

/// The signed distance (negative inside, 0 on the skin) and the closest skin
/// point of every point.
void octree_mesh_signed_distance(const OctreeMeshSetup& rSetup, const std::vector<Vec3>& rPoints,
                                 const OctreeMeshOptions& rOptions, OctreeMeshResult& rResult,
                                 std::vector<double>& rDistance, std::vector<Vec3>& rClosest);

/// Fill the report's octree counters (leaves per level, finest level).
void octree_mesh_report_tree(const LinearOctree& rTree, OctreeMeshResult& rResult);

/// The castellated output.
void octree_mesh_castellated(OctreeMeshSetup& rSetup, const OctreeMeshOptions& rOptions,
                             OctreeMeshResult& rResult);

/// A tetrahedral mesh between the stages of the tetrahedral output.
struct OctreeMeshTets {
    std::vector<Vec3> mPoints;
    std::vector<std::array<std::int64_t, 4>> mTets;
    /// Per tetrahedron, its material (0 outside, 1 inside, higher for nested shells).
    std::vector<std::int32_t> mMaterial;
    /// Per tetrahedron, 1 when stuffing or optimisation changed it.
    std::vector<std::uint8_t> mCut;
    /// Per tetrahedron, the level of the leaf it lies in.
    std::vector<std::int64_t> mLevel;
    /// Per tetrahedron, its lattice shape class (0 when cut).
    std::vector<std::int64_t> mShape;
    /// Per point, 1 when it lies on the skin or an interface.
    std::vector<std::uint8_t> mOnSkin;
    /// Per point, 1 on a feature corner, 2 on a feature chain.
    std::vector<std::uint8_t> mFeature;
    std::int32_t mNumMaterials = 2;
};

/// The lattice of the balanced octree, stuffed (or cleaved) to the skin.
void octree_mesh_stuffed_lattice(OctreeMeshSetup& rSetup, const OctreeMeshOptions& rOptions,
                                 OctreeMeshResult& rResult, OctreeMeshTets& rOut);

/// Snap the tetrahedra's boundary onto the skin's feature corners and chains.
void octree_mesh_snap_features(OctreeMeshSetup& rSetup, const OctreeMeshOptions& rOptions,
                               OctreeMeshResult& rResult, OctreeMeshTets& rTets);

/// Fit the skin's nodes along it and optimise the tetrahedra (ODT relocation,
/// flips) under the dihedral-floor guard.
void octree_mesh_optimize(OctreeMeshSetup& rSetup, const OctreeMeshOptions& rOptions,
                          OctreeMeshResult& rResult, OctreeMeshTets& rTets);

/// The two-sided distance between the triangles @p rFaces over @p rPoints
/// (an output's boundary) and the skin.
double octree_mesh_hausdorff_to_skin(const OctreeMeshSetup& rSetup,
                                     const std::vector<Vec3>& rPoints,
                                     const std::vector<std::array<std::int64_t, 3>>& rFaces);

/// The hexahedral output.
void octree_mesh_hexahedra(OctreeMeshSetup& rSetup, const OctreeMeshOptions& rOptions,
                           OctreeMeshResult& rResult);

/// The edge of the leaf holding @p rP (clamped into the forest).
double octree_mesh_local_size(const OctreeMeshSetup& rSetup, const Vec3& rP);

/// Write the tetrahedra (the inside materials, and the outside on request) to
/// the result's mesh, with their regions and data.
void octree_mesh_emit_tets(OctreeMeshSetup& rSetup, const OctreeMeshOptions& rOptions,
                           OctreeMeshResult& rResult, const OctreeMeshTets& rIn);

}  // namespace detail
}  // namespace meshioplusplus
