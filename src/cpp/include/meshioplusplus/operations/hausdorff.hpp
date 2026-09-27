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
 * @file operations/hausdorff.hpp
 * @brief The Hausdorff distance between two surfaces: how far apart they are
 * at their worst.
 *
 * `hausdorff_distance(a, b)` samples each surface, finds every sample's
 * unsigned distance to the *other* surface with the `distance_to_surface`
 * kernel (the bucket-grid nearest-triangle search, so the two cannot
 * disagree), and reduces: the one-sided distances are the largest sample
 * distance each way, and the Hausdorff distance is the larger of the two. It is
 * the scalar a remeshing, decimation or format round trip wants to assert on --
 * "no point moved further than this".
 *
 * ### Sampling, and why the answer is a lower bound
 *
 * The distance from a surface to another is attained somewhere on the first
 * surface, not necessarily at a vertex. With `mFaceSamples == 0` only the
 * vertices of each surface are sampled, which is exact when the farthest point
 * is a vertex (a surface compared with a refinement of itself, a rigid motion)
 * and a lower bound otherwise. `mFaceSamples = s > 0` also samples the centroid
 * of each of the `s * s` sub-triangles a triangle splits into, which tightens
 * the bound as `s` grows.
 *
 * ### Inputs
 *
 * Each mesh is a surface (triangles, quads, polygons), or a volume mesh, whose
 * boundary (`extract_surface`) is compared. A region names the `Cell` region of
 * surface cells to compare, and is refused on a volume mesh.
 *
 * Deterministic: samples are generated and reduced in a fixed order, and the
 * worst sample is the first one attaining the maximum.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// How `hausdorff_distance` samples the two surfaces.
struct HausdorffOptions {
    /// 0 samples the vertices only; `s > 0` also samples the centroids of the
    /// `s * s` sub-triangles of every triangle.
    std::int64_t mFaceSamples = 0;
    /// Restrict mesh A to this named `Cell` region of surface cells; empty
    /// takes every surface cell.
    std::string mRegionA;
    /// The same, for mesh B.
    std::string mRegionB;
    /// The bucket size of the nearest-triangle search, 0 for the automatic one.
    double mGridCellSize = 0.0;
};

/// The distances between two surfaces.
struct HausdorffResult {
    /// The Hausdorff distance, `max(mAtoB, mBtoA)`.
    double mDistance = 0.0;
    /// The largest distance from a sample of A to B.
    double mAtoB = 0.0;
    /// The largest distance from a sample of B to A.
    double mBtoA = 0.0;
    /// The mean and root-mean-square distance of A's samples to B.
    double mMeanAtoB = 0.0;
    double mRmsAtoB = 0.0;
    /// The mean and root-mean-square distance of B's samples to A.
    double mMeanBtoA = 0.0;
    double mRmsBtoA = 0.0;
    /// How many points of each surface were sampled.
    std::int64_t mNumSamplesA = 0;
    std::int64_t mNumSamplesB = 0;
    /// The sample of A farthest from B, and of B farthest from A.
    std::array<double, 3> mWorstPointA{{0.0, 0.0, 0.0}};
    std::array<double, 3> mWorstPointB{{0.0, 0.0, 0.0}};
};

/**
 * @brief The (sampled) Hausdorff distance between the surfaces of @p rA and @p rB.
 *
 * @throws std::invalid_argument when either mesh has no surface triangles, on a
 *         higher-order surface block (pointing at `linearize`), an unknown
 *         region name, a region on a volume mesh, or a negative `mFaceSamples`.
 */
MESHIOPLUSPLUS_API HausdorffResult hausdorff_distance(const Mesh& rA, const Mesh& rB,
                                                      const HausdorffOptions& rOptions = {});

}  // namespace meshioplusplus
