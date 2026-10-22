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

}  // namespace detail
}  // namespace meshioplusplus
