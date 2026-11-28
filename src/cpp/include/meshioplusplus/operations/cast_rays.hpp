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
 * @file operations/cast_rays.hpp
 * @brief The nearest surface hit of a batch of rays (roadmap 9.3).
 *
 * `cast_rays(mesh, origins, directions, options)` shoots every ray
 * `origin + t * direction` (`t` in length units: the direction is
 * normalised) at the surface cells of `mesh` and reports the nearest hit
 * with `0 < t <= max_distance`: whether it hit, the distance, the point, the
 * triangle's unit normal (its winding side, `(b - a) x (c - a)`) and the
 * global index of the **input cell** it belongs to (a quad or polygon is
 * fanned into triangles that all report their own cell).
 *
 * - The surface is the mesh's triangle, quad and rectangular polygon cells,
 *   optionally restricted to a named `Cell` region (`mRegion`). A volume
 *   mesh is refused by name, pointing at `extract_surface`; a higher-order
 *   block at `linearize`.
 * - A ray whose origin is on the surface does not hit it at `t = 0`; a hit
 *   needs `t > 0`.
 * - A ray through a shared edge or vertex hits once, at the lowest-numbered
 *   triangle (the traversal's barycentric tolerance of 1e-10 closes the seams).
 * - `directions` holds one direction (three numbers) for every ray, or one per
 *   ray; a zero or non-finite direction throws by name.
 * - A miss reports `hit = 0`, `NaN` distance, point and normal and cell `-1`.
 *
 * Nearest hit only: all crossings along a ray and a reusable query handle
 * are the remainder of roadmap 9.3. The rays are traced in parallel over a
 * bounding-volume hierarchy of the triangles; the result is the same on every
 * backend and thread count. The mesh is never modified.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// What `cast_rays` shoots at, and how far.
struct CastRaysOptions {
    /// Ignore hits farther than this; 0 (or any value not above 0) means no limit.
    double mMaxDistance = 0.0;
    /// Restrict the surface to this named `Cell` region; empty takes every surface cell.
    std::string mRegion;
    /// When `true`, return every crossing along each ray in CSR layout (`mRayOffset`).
    bool mAllHits = false;
};

/// One entry per ray (or per hit when `mAllHits` is `true`, indexed by `mRayOffset`).
struct CastRaysResult {
    std::vector<std::uint8_t> mHit;
    std::vector<double> mDistance;    ///< `NaN` on a miss (nearest mode); per hit when `mAllHits`.
    std::vector<double> mPoint;       ///< Three per ray (or per hit when `mAllHits`).
    std::vector<double> mNormal;      ///< Three per ray (or per hit), unit.
    std::vector<std::int64_t> mCell;  ///< Global input cell; -1 on a miss (nearest mode).
    std::int64_t mNumHits = 0;
    /// CSR offsets of length `num_rays + 1` when `mAllHits` is `true`; empty otherwise.
    std::vector<std::int64_t> mRayOffset;
};

/**
 * @brief Cast rays at the surface cells of @p rMesh (see the file comment).
 * @param rOrigins three numbers per ray.
 * @param rDirections three numbers, or three per ray.
 * @throws std::invalid_argument on a size that is not a multiple of three, a
 *         direction count that is neither one nor the ray count, a zero or
 *         non-finite direction or origin, a volume or higher-order block, an
 *         unknown region, or a mesh with no surface cell.
 */
MESHIOPLUSPLUS_API CastRaysResult cast_rays(const Mesh& rMesh, const std::vector<double>& rOrigins,
                                            const std::vector<double>& rDirections,
                                            const CastRaysOptions& rOptions = {});

}  // namespace meshioplusplus
