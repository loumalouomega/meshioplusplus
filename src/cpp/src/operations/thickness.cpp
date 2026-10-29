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
// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/thickness.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
// Project includes (private, not installed)
#include "../detail/tri_bvh.hpp"

namespace meshioplusplus {

namespace {

// One measured cell: where its ray starts and which way its normal points.
struct ThkFacet {
    std::int64_t mCell = -1;
    detail::Vec3 mOrigin{0.0, 0.0, 0.0};
    detail::Vec3 mNormal{0.0, 0.0, 0.0};  // unit; zero for a degenerate cell
};

// The cells of a soup, from its runs of triangles: the area-weighted centroid
// and the summed normal of each cell's fan.
std::vector<ThkFacet> thk_facets(const detail::TriangleSoup& rSoup) {
    std::vector<std::size_t> starts;
    for (std::size_t t = 0; t < rSoup.NumTriangles(); ++t)
        if (t == 0 || rSoup.mSourceCell[t] != rSoup.mSourceCell[t - 1])
            starts.push_back(t);
    starts.push_back(rSoup.NumTriangles());
    std::vector<ThkFacet> facets(starts.size() - 1);
    parallel_for(facets.size(), [&](std::size_t f) {
        ThkFacet& r_facet = facets[f];
        r_facet.mCell = rSoup.mSourceCell[starts[f]];
        detail::Vec3 normal{0.0, 0.0, 0.0}, moment{0.0, 0.0, 0.0}, mean{0.0, 0.0, 0.0};
        double area = 0.0;
        for (std::size_t t = starts[f]; t < starts[f + 1]; ++t) {
            const detail::Vec3* c = &rSoup.mCorners[3 * t];
            const detail::Vec3 n =
                detail::vec3_cross(detail::vec3_sub(c[1], c[0]), detail::vec3_sub(c[2], c[0]));
            const double a = 0.5 * detail::vec3_norm(n);
            const detail::Vec3 centre =
                detail::vec3_scale(detail::vec3_add(detail::vec3_add(c[0], c[1]), c[2]), 1.0 / 3.0);
            normal = detail::vec3_add(normal, n);
            moment = detail::vec3_add(moment, detail::vec3_scale(centre, a));
            mean = detail::vec3_add(mean, centre);
            area += a;
        }
        const double len = detail::vec3_norm(normal);
        r_facet.mOrigin =
            area > 0.0
                ? detail::vec3_scale(moment, 1.0 / area)
                : detail::vec3_scale(mean, 1.0 / static_cast<double>(starts[f + 1] - starts[f]));
        if (len > 0.0)
            r_facet.mNormal = detail::vec3_scale(normal, 1.0 / len);
    });
    return facets;
}

ThicknessResult thk_run(const Mesh& rMesh, const Mesh* pTarget, const ThicknessOptions& rOptions) {
    if (!std::isfinite(rOptions.mMaxDistance) || rOptions.mMaxDistance < 0.0)
        throw std::invalid_argument(
            "meshio++: thickness: max_distance must be finite and non-negative");
    if (!std::isfinite(rOptions.mThreshold) || rOptions.mThreshold < 0.0)
        throw std::invalid_argument(
            "meshio++: thickness: threshold must be finite and non-negative");
    const bool self = pTarget == nullptr;
    ThicknessResult result;
    result.mDirection = rOptions.mDirection != ThicknessDirection::Auto
                            ? rOptions.mDirection
                            : (self ? ThicknessDirection::Inward : ThicknessDirection::Outward);
    const bool inward = result.mDirection == ThicknessDirection::Inward;

    const detail::TriangleSoup soup = detail::build_triangle_soup(rMesh, rOptions.mRegion);
    const detail::TriangleSoup target =
        detail::build_triangle_soup(self ? rMesh : *pTarget, rOptions.mRegionTarget);
    if (soup.NumTriangles() == 0)
        result.mWarnings.emplace_back("the mesh has no 2-D cells to measure");
    if (target.NumTriangles() == 0)
        result.mWarnings.emplace_back("the target surface has no 2-D cells to reach");

    // The ray length: given, or the diagonal of both surfaces.
    double ray = rOptions.mMaxDistance;
    if (ray == 0.0) {
        const double inf = std::numeric_limits<double>::infinity();
        detail::Vec3 lo{inf, inf, inf}, hi{-inf, -inf, -inf};
        for (const detail::TriangleSoup* p_soup : {&soup, &target})
            for (const detail::Vec3& r_p : p_soup->mCorners)
                for (std::size_t d = 0; d < 3; ++d) {
                    lo[d] = std::min(lo[d], r_p[d]);
                    hi[d] = std::max(hi[d], r_p[d]);
                }
        ray = soup.mCorners.empty() ? 0.0 : detail::vec3_norm(detail::vec3_sub(hi, lo));
    }
    result.mMaxDistance = ray;

    const std::vector<ThkFacet> facets = thk_facets(soup);
    const detail::TriBvh tree = detail::build_tri_bvh(target);
    std::vector<double> dist(facets.size(), std::numeric_limits<double>::quiet_NaN());
    std::vector<std::int64_t> hit_cell(facets.size(), -1);
    parallel_for(
        facets.size(),
        [&](std::size_t f) {
            const ThkFacet& r_facet = facets[f];
            if (detail::vec3_norm_sq(r_facet.mNormal) == 0.0 || !(ray > 0.0))
                return;
            const detail::Vec3 dir =
                inward ? detail::vec3_scale(r_facet.mNormal, -1.0) : r_facet.mNormal;
            const std::int64_t own = self ? r_facet.mCell : -1;
            const detail::TriRayHit hit = detail::tri_bvh_ray_filtered(
                tree, r_facet.mOrigin, dir, 0.0, ray,
                [&](std::int64_t T, const detail::Vec3& rNormal) {
                    if (target.mSourceCell[static_cast<std::size_t>(T)] == own)
                        return false;
                    const double facing = detail::vec3_dot(rNormal, dir);
                    return inward ? facing > 0.0 : facing < 0.0;
                });
            if (hit.mTriangle < 0)
                return;
            dist[f] = hit.mT;
            hit_cell[f] = target.mSourceCell[static_cast<std::size_t>(hit.mTriangle)];
        },
        64);

    // The summary, in cell order.
    std::int64_t degenerate = 0;
    double sum = 0.0;
    for (std::size_t f = 0; f < facets.size(); ++f) {
        ++result.mNumCells;
        degenerate += detail::vec3_norm_sq(facets[f].mNormal) == 0.0;
        if (hit_cell[f] < 0) {
            ++result.mNumMissed;
            continue;
        }
        const double d = dist[f];
        if (result.mNumHit == 0 || d < result.mMin) {
            result.mMin = d;
            result.mMinCell = facets[f].mCell;
            result.mMinPoint = facets[f].mOrigin;
        }
        result.mMax = result.mNumHit == 0 ? d : std::max(result.mMax, d);
        sum += d;
        ++result.mNumHit;
        if (rOptions.mThreshold > 0.0 && d < rOptions.mThreshold)
            ++result.mNumBelow;
    }
    if (result.mNumHit > 0)
        result.mMean = sum / static_cast<double>(result.mNumHit);
    if (degenerate > 0)
        result.mWarnings.emplace_back(std::to_string(degenerate) +
                                      " cells have no area and no normal; they read as misses");

    result.mMesh = detail::clone_mesh(
        rMesh, [&](DataLocation Location, const std::string& rName, std::string&) {
            return !(rOptions.mAttach && Location == DataLocation::Cell &&
                     (rName == kThicknessDistanceName || rName == kThicknessHitCellName));
        });
    if (rOptions.mAttach) {
        const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
        const auto total = static_cast<std::size_t>(detail::total_cells(bases));
        std::vector<double> all_dist(total, std::numeric_limits<double>::quiet_NaN());
        std::vector<std::int64_t> all_hit(total, -1);
        for (std::size_t f = 0; f < facets.size(); ++f) {
            const auto c = static_cast<std::size_t>(facets[f].mCell);
            all_dist[c] = dist[f];
            all_hit[c] = hit_cell[f];
        }
        std::vector<NDArray> dist_blocks, hit_blocks;
        for (std::size_t b = 0; b + 1 < bases.size(); ++b) {
            const auto n = static_cast<std::size_t>(bases[b + 1] - bases[b]);
            NDArray d = NDArray::Uninit(DType::Float64, {n});
            std::copy(all_dist.begin() + bases[b], all_dist.begin() + bases[b + 1], d.As<double>());
            NDArray h = NDArray::Uninit(DType::Int64, {n});
            std::copy(all_hit.begin() + bases[b], all_hit.begin() + bases[b + 1],
                      h.As<std::int64_t>());
            dist_blocks.push_back(std::move(d));
            hit_blocks.push_back(std::move(h));
        }
        result.mMesh.AddCellData(kThicknessDistanceName, std::move(dist_blocks));
        result.mMesh.AddCellData(kThicknessHitCellName, std::move(hit_blocks));
    }
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        result.mMesh.AddRegion(rMesh.Region(i));
    return result;
}

}  // namespace

ThicknessResult thickness(const Mesh& rMesh, const ThicknessOptions& rOptions) {
    return thk_run(rMesh, nullptr, rOptions);
}

ThicknessResult thickness(const Mesh& rMesh, const Mesh& rTarget,
                          const ThicknessOptions& rOptions) {
    return thk_run(rMesh, &rTarget, rOptions);
}

}  // namespace meshioplusplus
