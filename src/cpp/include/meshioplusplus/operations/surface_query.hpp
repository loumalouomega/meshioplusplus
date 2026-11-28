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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file operations/surface_query.hpp
 * @brief A reusable surface query handle (roadmap 9.3.2): builds a surface's
 * triangle soup and its bounding-volume hierarchies once, then answers many
 * batches of ray casts and closest-point queries against them.
 *
 * `cast_rays(mesh, ...)` is `SurfaceQuery(mesh, region).CastRays(...)`, so a
 * batch answered through the handle is bit-identical to the one-shot call.
 * The handle owns copies of everything it needs: the mesh may change or go
 * away after construction. It is safe to query from several threads at once.
 */

// System includes
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/cast_rays.hpp"

namespace meshioplusplus {

/// The nearest point on the surface to each query point.
struct SurfaceClosestResult {
    std::vector<double> mPoint;       ///< Three per query point.
    std::vector<double> mDistance;    ///< Unsigned distance, one per query point.
    std::vector<std::int64_t> mCell;  ///< Global input cell of the nearest point.
};

/// A surface's triangle soup and hierarchies, built once (roadmap 9.3.2).
class MESHIOPLUSPLUS_API SurfaceQuery {
public:
    /// The triangle, quad and polygon cells of @p rMesh, or only those of the
    /// cell region @p rRegion. Throws `std::invalid_argument` on an unknown
    /// region or a mesh with no surface cell.
    explicit SurfaceQuery(const Mesh& rMesh, const std::string& rRegion = "");
    ~SurfaceQuery();
    SurfaceQuery(SurfaceQuery&& rOther) noexcept;
    SurfaceQuery& operator=(SurfaceQuery&& rOther) noexcept;
    SurfaceQuery(const SurfaceQuery&) = delete;
    SurfaceQuery& operator=(const SurfaceQuery&) = delete;

    /// The surface's triangles (quads and polygons fanned).
    std::size_t NumTriangles() const noexcept;

    /// `cast_rays` against this surface; `rOptions.mRegion` must be empty (the
    /// region was chosen at construction).
    CastRaysResult CastRays(const std::vector<double>& rOrigins,
                            const std::vector<double>& rDirections,
                            const CastRaysOptions& rOptions = {}) const;

    /// The nearest surface point to each point of @p rPoints (three numbers each).
    SurfaceClosestResult ClosestPoints(const std::vector<double>& rPoints) const;

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

}  // namespace meshioplusplus
