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
#include "../detail/lsq_fit.hpp"
#include "../detail/render_field.hpp"
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

constexpr double kThkPi = 3.14159265358979323846;

// The fan of a cone: the central direction, then a Fibonacci spiral of
// `Count - 1` directions within `HalfAngle` radians of it, in the frame
// `lsq_frame` completes around it (deterministic).
std::vector<detail::Vec3> thk_fan(const detail::Vec3& rAxis, double HalfAngle, std::int64_t Count) {
    std::vector<detail::Vec3> out{rAxis};
    if (Count <= 1)
        return out;
    detail::Vec3 e1{0.0, 0.0, 0.0}, e2{0.0, 0.0, 0.0};
    detail::lsq_frame(rAxis.data(), e1.data(), e2.data());
    const double golden = kThkPi * (3.0 - std::sqrt(5.0));
    const auto n = static_cast<double>(Count - 1);
    for (std::int64_t k = 0; k < Count - 1; ++k) {
        const double polar = HalfAngle * std::sqrt((static_cast<double>(k) + 0.5) / n);
        const double az = golden * static_cast<double>(k);
        const double sp = std::sin(polar);
        out.push_back(
            detail::vec3_add(detail::vec3_scale(rAxis, std::cos(polar)),
                             detail::vec3_add(detail::vec3_scale(e1, sp * std::cos(az)),
                                              detail::vec3_scale(e2, sp * std::sin(az)))));
    }
    return out;
}

// The vertex ids of each facet's triangles, sorted and unique.
std::vector<std::vector<std::int64_t>> thk_facet_vertices(const detail::TriangleSoup& rSoup,
                                                          std::size_t NumFacets) {
    std::vector<std::vector<std::int64_t>> out(NumFacets);
    std::size_t f = 0;
    for (std::size_t t = 0; t < rSoup.NumTriangles(); ++t) {
        if (t > 0 && rSoup.mSourceCell[t] != rSoup.mSourceCell[t - 1])
            ++f;
        for (const std::int64_t v : rSoup.mVertices[t])
            out[f].push_back(v);
    }
    for (std::vector<std::int64_t>& r_v : out) {
        std::sort(r_v.begin(), r_v.end());
        r_v.erase(std::unique(r_v.begin(), r_v.end()), r_v.end());
    }
    return out;
}

ThicknessResult thk_run(const Mesh& rMesh, const Mesh* pTarget, const ThicknessOptions& rOptions) {
    if (!std::isfinite(rOptions.mMaxDistance) || rOptions.mMaxDistance < 0.0)
        throw std::invalid_argument(
            "meshio++: thickness: max_distance must be finite and non-negative");
    if (!std::isfinite(rOptions.mThreshold) || rOptions.mThreshold < 0.0)
        throw std::invalid_argument(
            "meshio++: thickness: threshold must be finite and non-negative");
    if (!(rOptions.mMaxNormalDeviationDeg > 0.0 && rOptions.mMaxNormalDeviationDeg <= 90.0))
        throw std::invalid_argument(
            "meshio++: thickness: max_normal_deviation must lie in (0, 90] degrees");
    if (rOptions.mMethod == ThicknessMethod::Cone &&
        (!(rOptions.mConeAngleDeg > 0.0 && rOptions.mConeAngleDeg < 90.0) ||
         rOptions.mConeRays < 1))
        throw std::invalid_argument(
            "meshio++: thickness: the cone needs an angle in (0, 90) degrees and at least one "
            "ray");
    if (rOptions.mMethod != ThicknessMethod::Ray && rOptions.mMethod != ThicknessMethod::Sphere &&
        rOptions.mMethod != ThicknessMethod::Cone)
        throw std::invalid_argument("meshio++: thickness: unknown method");
    const bool self = pTarget == nullptr;
    ThicknessResult result;
    result.mMethod = rOptions.mMethod;
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
        if (rOptions.mTwoSided)
            ray *= 0.2;  // a scan's walls: a fifth of the diagonal
    }
    result.mMaxDistance = ray;

    const std::vector<ThkFacet> facets = thk_facets(soup);
    const detail::TriBvh tree = detail::build_tri_bvh(target);
    std::vector<double> dist(facets.size(), std::numeric_limits<double>::quiet_NaN());
    std::vector<std::int64_t> hit_cell(facets.size(), -1);
    const bool legacy = rOptions.mMethod == ThicknessMethod::Ray && !rOptions.mTwoSided &&
                        rOptions.mMaxNormalDeviationDeg >= 90.0;
    if (legacy) {
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
    } else {
        // Exit (inward), entry (outward) or either face (two-sided), within the
        // normal deviation; 90 degrees keeps the strict sign test.
        const bool strict = rOptions.mMaxNormalDeviationDeg >= 90.0;
        const double cos_dev = std::cos(rOptions.mMaxNormalDeviationDeg * kThkPi / 180.0);
        const auto facing_ok = [&](double Facing) {
            return strict ? Facing > 0.0 : Facing >= cos_dev;
        };
        std::vector<detail::Vec3> hit_dir(facets.size(), detail::Vec3{0.0, 0.0, 0.0});
        parallel_for(
            facets.size(),
            [&](std::size_t f) {
                const ThkFacet& r_facet = facets[f];
                if (detail::vec3_norm_sq(r_facet.mNormal) == 0.0 || !(ray > 0.0))
                    return;
                const detail::Vec3 axis =
                    inward ? detail::vec3_scale(r_facet.mNormal, -1.0) : r_facet.mNormal;
                const std::int64_t own = self ? r_facet.mCell : -1;
                std::vector<detail::Vec3> dirs =
                    rOptions.mMethod == ThicknessMethod::Cone
                        ? thk_fan(axis, rOptions.mConeAngleDeg * kThkPi / 180.0, rOptions.mConeRays)
                        : std::vector<detail::Vec3>{axis};
                if (rOptions.mTwoSided) {
                    const std::size_t n = dirs.size();
                    for (std::size_t k = 0; k < n; ++k)
                        dirs.push_back(detail::vec3_scale(dirs[k], -1.0));
                }
                for (const detail::Vec3& r_dir : dirs) {
                    const detail::TriRayHit hit = detail::tri_bvh_ray_filtered(
                        tree, r_facet.mOrigin, r_dir, 0.0, ray,
                        [&](std::int64_t T, const detail::Vec3& rNormal) {
                            if (target.mSourceCell[static_cast<std::size_t>(T)] == own)
                                return false;
                            const double len = detail::vec3_norm(rNormal);
                            if (!(len > 0.0))
                                return false;
                            const double facing = detail::vec3_dot(rNormal, r_dir) / len;
                            if (rOptions.mTwoSided)
                                return facing_ok(std::abs(facing));
                            return facing_ok(inward ? facing : -facing);
                        });
                    if (hit.mTriangle < 0 || !(hit.mT < dist[f] || std::isnan(dist[f])))
                        continue;
                    dist[f] = hit.mT;
                    hit_cell[f] = target.mSourceCell[static_cast<std::size_t>(hit.mTriangle)];
                    hit_dir[f] = r_dir;
                }
            },
            16);
        if (rOptions.mMethod == ThicknessMethod::Sphere) {
            // Twice the distance from each ray's midpoint to the nearest
            // triangle outside the facet's one-ring, never more than the ray.
            const detail::DistanceQuery query =
                detail::build_distance_query(target, SurfaceDistanceOptions{});
            const std::vector<std::vector<std::int64_t>> ring =
                thk_facet_vertices(soup, facets.size());
            std::vector<std::size_t> which;
            std::vector<detail::Vec3> mids;
            for (std::size_t f = 0; f < facets.size(); ++f)
                if (hit_cell[f] >= 0) {
                    which.push_back(f);
                    mids.push_back(detail::vec3_add(facets[f].mOrigin,
                                                    detail::vec3_scale(hit_dir[f], 0.5 * dist[f])));
                }
            const detail::TriangleFilter accept = [&](std::size_t Q, std::int64_t T) {
                if (!self)
                    return true;
                const std::vector<std::int64_t>& r_ring = ring[which[Q]];
                for (const std::int64_t v : target.mVertices[static_cast<std::size_t>(T)])
                    if (std::binary_search(r_ring.begin(), r_ring.end(), v))
                        return false;
                return true;
            };
            const std::vector<detail::SurfaceProjection> proj =
                detail::query_surface_projections_filtered(query, mids, accept, ray);
            for (std::size_t q = 0; q < which.size(); ++q) {
                const std::size_t f = which[q];
                if (!proj[q].mFound)
                    continue;
                const double d = 2.0 * proj[q].mDistance;
                if (d < dist[f]) {
                    dist[f] = d;
                    hit_cell[f] = proj[q].mSourceCell;
                }
            }
        }
    }

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
    if (result.mNumHit > 0) {
        result.mMean = sum / static_cast<double>(result.mNumHit);
        std::vector<double> hits;
        hits.reserve(static_cast<std::size_t>(result.mNumHit));
        for (std::size_t f = 0; f < facets.size(); ++f)
            if (hit_cell[f] >= 0)
                hits.push_back(dist[f]);
        result.mP2 = detail::render_percentile(hits, 2.0);
        result.mP95 = detail::render_percentile(hits, 95.0);
    }
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

ThicknessMethod thickness_method_from_name(const std::string& rName) {
    if (rName == "ray")
        return ThicknessMethod::Ray;
    if (rName == "sphere")
        return ThicknessMethod::Sphere;
    if (rName == "cone")
        return ThicknessMethod::Cone;
    throw std::invalid_argument("meshio++: thickness: unknown method '" + rName +
                                "' (expected 'ray', 'sphere' or 'cone')");
}

ThicknessResult thickness(const Mesh& rMesh, const ThicknessOptions& rOptions) {
    return thk_run(rMesh, nullptr, rOptions);
}

ThicknessResult thickness(const Mesh& rMesh, const Mesh& rTarget,
                          const ThicknessOptions& rOptions) {
    return thk_run(rMesh, &rTarget, rOptions);
}

}  // namespace meshioplusplus
