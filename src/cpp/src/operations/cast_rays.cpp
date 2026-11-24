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
// Nearest-hit ray casting: see operations/cast_rays.hpp for the contract. One
// bounding-volume hierarchy over the surface's triangle soup, one traversal
// per ray (`detail::tri_bvh_ray`), in parallel.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/cast_rays.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/tri_bvh.hpp"

namespace meshioplusplus {

CastRaysResult cast_rays(const Mesh& rMesh, const std::vector<double>& rOrigins,
                         const std::vector<double>& rDirections, const CastRaysOptions& rOptions) {
    static const char* const kPrefix = "meshio++: cast_rays: ";
    if (rOrigins.size() % 3 != 0)
        throw std::invalid_argument(std::string(kPrefix) +
                                    "origins must hold three numbers per ray");
    if (rDirections.size() % 3 != 0 || rDirections.empty())
        throw std::invalid_argument(std::string(kPrefix) +
                                    "directions must hold three numbers, or three per ray");
    const std::size_t n = rOrigins.size() / 3;
    const bool common = rDirections.size() == 3;
    if (!common && rDirections.size() != rOrigins.size())
        throw std::invalid_argument(std::string(kPrefix) +
                                    "directions must be one direction or one per ray (" +
                                    std::to_string(n) + " rays, " +
                                    std::to_string(rDirections.size() / 3) + " directions)");

    // Unit directions, validated before any work.
    std::vector<detail::Vec3> dirs(common ? 1 : n);
    for (std::size_t i = 0; i < dirs.size(); ++i) {
        const detail::Vec3 d{rDirections[3 * i], rDirections[3 * i + 1], rDirections[3 * i + 2]};
        const double len = detail::vec3_norm(d);
        if (!std::isfinite(len) || len == 0.0)
            throw std::invalid_argument(std::string(kPrefix) + "direction " + std::to_string(i) +
                                        " is zero or not finite");
        dirs[i] = detail::vec3_scale(d, 1.0 / len);
    }
    for (std::size_t k = 0; k < rOrigins.size(); ++k)
        if (!std::isfinite(rOrigins[k]))
            throw std::invalid_argument(std::string(kPrefix) + "origin " +
                                        std::to_string(k / 3) + " is not finite");

    const detail::TriangleSoup soup = detail::build_triangle_soup(rMesh, rOptions.mRegion);
    if (soup.NumTriangles() == 0)
        throw std::invalid_argument(std::string(kPrefix) + "the mesh has no surface cell to hit");
    const detail::TriBvh tree = detail::build_tri_bvh(soup);
    const double t_max = rOptions.mMaxDistance > 0.0 ? rOptions.mMaxDistance
                                                     : std::numeric_limits<double>::infinity();
    const double nan = std::nan("");

    CastRaysResult out;
    out.mHit.assign(n, 0);
    out.mDistance.assign(n, nan);
    out.mPoint.assign(3 * n, nan);
    out.mNormal.assign(3 * n, nan);
    out.mCell.assign(n, -1);
    parallel_for(n, [&](std::size_t i) {
        const detail::Vec3 origin{rOrigins[3 * i], rOrigins[3 * i + 1], rOrigins[3 * i + 2]};
        const detail::Vec3& dir = dirs[common ? 0 : i];
        const detail::TriRayHit hit = detail::tri_bvh_ray(tree, origin, dir, 0.0, t_max);
        if (hit.mTriangle < 0)
            return;
        out.mHit[i] = 1;
        out.mDistance[i] = hit.mT;
        for (std::size_t k = 0; k < 3; ++k)
            out.mPoint[3 * i + k] = origin[k] + hit.mT * dir[k];
        const double len = detail::vec3_norm(hit.mNormal);
        if (len > 0.0)
            for (std::size_t k = 0; k < 3; ++k)
                out.mNormal[3 * i + k] = hit.mNormal[k] / len;
        out.mCell[i] = soup.mSourceCell[static_cast<std::size_t>(hit.mTriangle)];
    });
    for (std::size_t i = 0; i < n; ++i)
        out.mNumHits += out.mHit[i];
    return out;
}

}  // namespace meshioplusplus
