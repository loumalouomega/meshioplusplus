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
 * @file operations/filter_points.hpp
 * @brief Conditioning of a point cloud: outlier removal and thinning (roadmap 9.9).
 *
 * Scans read from PLY, PCD, XYZ or LAS carry noise and far more points than a
 * reconstruction or a registration needs. `filter_points(mesh, options)`
 * applies up to four steps, in this order, to the mesh's points, each step
 * seeing the points the one before kept; a step is off at its default of 0:
 *
 * 1. **Statistical outlier removal** (`mStatisticalK` > 0; R. B. Rusu et al.,
 *    "Towards 3D point cloud based object maps for household environments",
 *    Robotics and Autonomous Systems 56(11), 2008): each point's mean distance
 *    to its `mStatisticalK` nearest other points is compared with the mean
 *    `mu` and the (population) standard deviation `sigma` of those means over
 *    all points; a point whose mean exceeds `mu + mStdRatio * sigma` goes.
 * 2. **Radius outlier removal** (`mRadius` > 0): a point with fewer than
 *    `mMinNeighbors` other points within `mRadius` goes.
 * 3. **Voxel-grid thinning** (`mVoxelSize` > 0): the points are binned into
 *    cubes of that side anchored at the bounding box's minimum corner, and one
 *    **real point** per occupied cube survives: the one nearest the mean of its
 *    cube's points (ties to the lower id), so the output is a subset of the
 *    input with every array intact. Cubes are visited in ascending (x, y, z)
 *    key order, which is the output order of this step.
 * 4. **Poisson-disk thinning** (`mMinDistance` > 0; R. Bridson, "Fast Poisson
 *    disk sampling in arbitrary dimensions", SIGGRAPH Sketches 2007, here as
 *    greedy dart throwing in index order, with no random numbers): a point is
 *    kept when no already-kept point lies closer than `mMinDistance`. The kept
 *    points are at least that far apart, and every dropped point is within it
 *    of a kept one.
 *
 * The result is a point cloud: the surviving points (a step keeps the order of
 * the one before it; only the voxel step reorders), every per-point array
 * gathered with them, one `vertex` block, and `mPointMap` giving each output
 * row's input point. Point regions follow the surviving points (a group that
 * loses every point is kept, empty); cells and Cell and Side regions are
 * dropped, with a warning. Nothing else is changed: the input is never modified.
 *
 * Every search is exact and deterministic (`detail/point_knn.hpp`: ties to the
 * lower id), the per-point work runs in parallel and every reduction is serial
 * in input order, so the result is the same on every backend and thread count.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// How voxel-grid thinning represents each occupied cube (roadmap 9.9.1).
enum class VoxelMode : std::int32_t {
    Nearest = 0,  ///< Keep the real input point nearest the cube's mean (default).
    Average = 1,  ///< Emit the mean of each cube's coordinates and float arrays (`mPointMap` = -1).
};

/// Which steps `filter_points` runs and how strictly.
struct PointFilterOptions {
    /// Statistical outlier removal: neighbours per mean distance; 0 is off.
    std::int64_t mStatisticalK = 0;
    /// Statistical outlier removal: standard deviations above the mean.
    double mStdRatio = 1.0;
    /// Radius outlier removal: the search radius; 0 is off.
    double mRadius = 0.0;
    /// Radius outlier removal: other points needed within `mRadius`.
    std::int64_t mMinNeighbors = 4;
    /// Voxel thinning: the cube side; 0 is off.
    double mVoxelSize = 0.0;
    /// Poisson-disk thinning: the minimum distance between kept points; 0 is off.
    double mMinDistance = 0.0;
    /// Voxel thinning representative: nearest real point or cube mean.
    VoxelMode mVoxelMode = VoxelMode::Nearest;
};

/// What `filter_points` kept.
struct PointFilterResult {
    /// The point cloud.
    Mesh mMesh;
    /// Per output point, the input point it is.
    std::vector<std::int64_t> mPointMap;
    std::int64_t mNumInput = 0;
    std::int64_t mNumRemovedStatistical = 0;
    std::int64_t mNumRemovedRadius = 0;
    std::int64_t mNumRemovedVoxel = 0;
    std::int64_t mNumRemovedPoisson = 0;
    std::vector<std::string> mWarnings;
};

/**
 * @brief Remove outliers and thin a point cloud (see the file comment).
 * @throws std::invalid_argument on a negative or non-finite parameter, a
 *         `mMinNeighbors` below 0, a non-finite coordinate, no step enabled, or
 *         a voxel size too small for the cloud's extent.
 */
MESHIOPLUSPLUS_API PointFilterResult filter_points(const Mesh& rMesh,
                                                   const PointFilterOptions& rOptions = {});

}  // namespace meshioplusplus
