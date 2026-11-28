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
#include "meshioplusplus/operations/surface_query.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/tri_bvh.hpp"

namespace meshioplusplus {

namespace {

const char* const kCastRaysPrefix = "meshio++: cast_rays: ";

// The unit directions of a batch, validated (with the origins) before any work.
std::vector<detail::Vec3> cast_rays_unit_directions(const std::vector<double>& rOrigins,
                                                    const std::vector<double>& rDirections) {
    const char* const kPrefix = kCastRaysPrefix;
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
            throw std::invalid_argument(std::string(kPrefix) + "origin " + std::to_string(k / 3) +
                                        " is not finite");
    return dirs;
}

}  // namespace

struct SurfaceQuery::Impl {
    detail::TriangleSoup mSoup;
    detail::TriBvh mTree;  // points into mSoup, which never moves (Impl is heap-held)
    detail::DistanceQuery mDistance;
};

SurfaceQuery::SurfaceQuery(const Mesh& rMesh, const std::string& rRegion)
    : mpImpl(std::make_unique<Impl>()) {
    mpImpl->mSoup = detail::build_triangle_soup(rMesh, rRegion);
    if (mpImpl->mSoup.NumTriangles() == 0)
        throw std::invalid_argument(std::string(kCastRaysPrefix) +
                                    "the mesh has no surface cell to hit");
    mpImpl->mTree = detail::build_tri_bvh(mpImpl->mSoup);
    mpImpl->mDistance = detail::build_distance_query(mpImpl->mSoup, SurfaceDistanceOptions{});
}

SurfaceQuery::~SurfaceQuery() = default;
SurfaceQuery::SurfaceQuery(SurfaceQuery&& rOther) noexcept = default;
SurfaceQuery& SurfaceQuery::operator=(SurfaceQuery&& rOther) noexcept = default;

std::size_t SurfaceQuery::NumTriangles() const noexcept {
    return mpImpl ? mpImpl->mSoup.NumTriangles() : 0;
}

SurfaceClosestResult SurfaceQuery::ClosestPoints(const std::vector<double>& rPoints) const {
    if (rPoints.size() % 3 != 0)
        throw std::invalid_argument(
            "meshio++: SurfaceQuery: points must hold three numbers per point");
    const std::size_t n = rPoints.size() / 3;
    std::vector<detail::Vec3> pts(n);
    for (std::size_t i = 0; i < n; ++i) {
        pts[i] = {rPoints[3 * i], rPoints[3 * i + 1], rPoints[3 * i + 2]};
        for (const double v : pts[i])
            if (!std::isfinite(v))
                throw std::invalid_argument("meshio++: SurfaceQuery: point " + std::to_string(i) +
                                            " is not finite");
    }
    const std::vector<detail::ClosestPointHit> hits =
        detail::query_closest_points(mpImpl->mDistance, pts);
    SurfaceClosestResult out;
    out.mPoint.resize(3 * n);
    out.mDistance.resize(n);
    out.mCell.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < 3; ++k)
            out.mPoint[3 * i + k] = hits[i].mPoint[k];
        out.mDistance[i] = hits[i].mDistance;
        out.mCell[i] = hits[i].mSourceCell;
    }
    return out;
}

CastRaysResult SurfaceQuery::CastRays(const std::vector<double>& rOrigins,
                                      const std::vector<double>& rDirections,
                                      const CastRaysOptions& rOptions) const {
    if (!rOptions.mRegion.empty())
        throw std::invalid_argument(std::string(kCastRaysPrefix) +
                                    "a SurfaceQuery's region is chosen at construction");
    const std::vector<detail::Vec3> dirs = cast_rays_unit_directions(rOrigins, rDirections);
    const std::size_t n = rOrigins.size() / 3;
    const bool common = rDirections.size() == 3;
    const detail::TriangleSoup& soup = mpImpl->mSoup;
    const detail::TriBvh& tree = mpImpl->mTree;
    const double t_max = rOptions.mMaxDistance > 0.0 ? rOptions.mMaxDistance
                                                     : std::numeric_limits<double>::infinity();
    const double nan = std::nan("");

    CastRaysResult out;
    if (rOptions.mAllHits) {
        std::vector<std::vector<detail::TriRayHit>> per_ray(n);
        parallel_for(n, [&](std::size_t i) {
            const detail::Vec3 origin{rOrigins[3 * i], rOrigins[3 * i + 1], rOrigins[3 * i + 2]};
            const detail::Vec3& dir = dirs[common ? 0 : i];
            per_ray[i] = detail::tri_bvh_ray_all_hits(tree, origin, dir, 0.0, t_max);
        });
        out.mHit.assign(n, 0);
        out.mRayOffset.assign(n + 1, 0);
        for (std::size_t i = 0; i < n; ++i) {
            out.mHit[i] = per_ray[i].empty() ? 0 : 1;
            out.mRayOffset[i + 1] =
                out.mRayOffset[i] + static_cast<std::int64_t>(per_ray[i].size());
        }
        const auto total = static_cast<std::size_t>(out.mRayOffset[n]);
        out.mNumHits = static_cast<std::int64_t>(total);
        out.mDistance.assign(total, nan);
        out.mPoint.assign(3 * total, nan);
        out.mNormal.assign(3 * total, nan);
        out.mCell.assign(total, -1);
        parallel_for(n, [&](std::size_t i) {
            const detail::Vec3 origin{rOrigins[3 * i], rOrigins[3 * i + 1], rOrigins[3 * i + 2]};
            const detail::Vec3& dir = dirs[common ? 0 : i];
            const auto base = static_cast<std::size_t>(out.mRayOffset[i]);
            const auto& hits = per_ray[i];
            for (std::size_t h = 0; h < hits.size(); ++h) {
                const std::size_t dst = base + h;
                out.mDistance[dst] = hits[h].mT;
                for (std::size_t k = 0; k < 3; ++k)
                    out.mPoint[3 * dst + k] = origin[k] + hits[h].mT * dir[k];
                const double len = detail::vec3_norm(hits[h].mNormal);
                if (len > 0.0)
                    for (std::size_t k = 0; k < 3; ++k)
                        out.mNormal[3 * dst + k] = hits[h].mNormal[k] / len;
                out.mCell[dst] = soup.mSourceCell[static_cast<std::size_t>(hits[h].mTriangle)];
            }
        });
        return out;
    }

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

CastRaysResult cast_rays(const Mesh& rMesh, const std::vector<double>& rOrigins,
                         const std::vector<double>& rDirections, const CastRaysOptions& rOptions) {
    // The rays are validated before the surface is built, as they always were.
    cast_rays_unit_directions(rOrigins, rDirections);
    const SurfaceQuery query(rMesh, rOptions.mRegion);
    CastRaysOptions options = rOptions;
    options.mRegion.clear();
    return query.CastRays(rOrigins, rDirections, options);
}

}  // namespace meshioplusplus
