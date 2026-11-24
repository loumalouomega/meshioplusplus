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
// Point-cloud conditioning: see operations/filter_points.hpp for the contract.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/filter_points.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/bounding_box_utils.hpp"
#include "../detail/bucket_build.hpp"
#include "../detail/point_knn.hpp"
#include "../detail/point_locator.hpp"
#include "../detail/subset_gather.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;
using FpIds = std::vector<std::int64_t>;

// The coordinates of the points `rAlive` names, in that order.
std::vector<Vec3> fp_cloud(const std::vector<Vec3>& rXyz, const FpIds& rAlive) {
    std::vector<Vec3> out(rAlive.size());
    for (std::size_t i = 0; i < rAlive.size(); ++i)
        out[i] = rXyz[static_cast<std::size_t>(rAlive[i])];
    return out;
}

// The grid cell size of `rCloud` by the shared rule.
double fp_cell_size(const std::vector<Vec3>& rCloud) {
    double lo[3], hi[3];
    detail::BoundingBoxUtils::FromCoords(
        rCloud.size(), [&](std::size_t i) { return rCloud[i]; }, lo, hi);
    return detail::point_grid_cell_size(lo, hi, rCloud.size());
}

// Step 1: statistical outlier removal.
FpIds fp_statistical(const std::vector<Vec3>& rXyz, const FpIds& rAlive, std::int64_t K,
                     double StdRatio, std::vector<std::string>& rWarnings) {
    const std::size_t m = rAlive.size();
    if (m < 2) {
        rWarnings.push_back("statistical outlier removal skipped: fewer than two points");
        return rAlive;
    }
    const std::size_t k = std::min<std::size_t>(static_cast<std::size_t>(K), m - 1);
    const std::vector<Vec3> cloud = fp_cloud(rXyz, rAlive);
    const detail::PointKnn knn(cloud, fp_cell_size(cloud));
    std::vector<double> mean_dist(m, 0.0);
    parallel_for(m, [&](std::size_t i) {
        std::vector<std::pair<double, std::int64_t>> buf(k + 1);
        const std::size_t found = knn.KNearest(cloud[i], k + 1, buf.data());
        double sum = 0.0;
        std::size_t used = 0;
        for (std::size_t j = 0; j < found && used < k; ++j) {
            if (buf[j].second == static_cast<std::int64_t>(i))
                continue;  // the point itself
            sum += std::sqrt(buf[j].first);
            ++used;
        }
        mean_dist[i] = used > 0 ? sum / static_cast<double>(used) : 0.0;
    });
    double mu = 0.0;
    for (double d : mean_dist)
        mu += d;
    mu /= static_cast<double>(m);
    double var = 0.0;
    for (double d : mean_dist)
        var += (d - mu) * (d - mu);
    const double sigma = std::sqrt(var / static_cast<double>(m));
    const double threshold = mu + StdRatio * sigma;
    // A relative guard so equal means (sigma 0) are not split by rounding.
    const double limit = threshold + 1e-12 * std::abs(threshold);
    FpIds out;
    out.reserve(m);
    for (std::size_t i = 0; i < m; ++i)
        if (mean_dist[i] <= limit)
            out.push_back(rAlive[i]);
    return out;
}

// Step 2: radius outlier removal.
FpIds fp_radius(const std::vector<Vec3>& rXyz, const FpIds& rAlive, double Radius,
                std::int64_t MinNeighbors) {
    const std::size_t m = rAlive.size();
    if (m == 0)
        return rAlive;
    const std::vector<Vec3> cloud = fp_cloud(rXyz, rAlive);
    const detail::PointKnn knn(cloud, Radius);
    std::vector<std::uint8_t> keep(m, 0);
    const std::size_t need = static_cast<std::size_t>(MinNeighbors);
    parallel_for(m, [&](std::size_t i) {
        // The point counts itself (distance 0), so it needs one more.
        keep[i] = knn.CountWithin(cloud[i], Radius, need + 1) >= need + 1 ? 1 : 0;
    });
    FpIds out;
    out.reserve(m);
    for (std::size_t i = 0; i < m; ++i)
        if (keep[i])
            out.push_back(rAlive[i]);
    return out;
}

struct FpVoxelRow {
    std::int64_t mKx, mKy, mKz, mId;  // mId indexes the alive list
    bool operator<(const FpVoxelRow& r) const {
        if (mKx != r.mKx)
            return mKx < r.mKx;
        if (mKy != r.mKy)
            return mKy < r.mKy;
        if (mKz != r.mKz)
            return mKz < r.mKz;
        return mId < r.mId;
    }
    bool SameVoxel(const FpVoxelRow& r) const {
        return mKx == r.mKx && mKy == r.mKy && mKz == r.mKz;
    }
};

// Step 3: one real point per occupied voxel, the one nearest the voxel's mean.
FpIds fp_voxel(const std::vector<Vec3>& rXyz, const FpIds& rAlive, double Size) {
    const std::size_t m = rAlive.size();
    if (m == 0)
        return rAlive;
    const std::vector<Vec3> cloud = fp_cloud(rXyz, rAlive);
    double lo[3], hi[3];
    detail::BoundingBoxUtils::FromCoords(
        m, [&](std::size_t i) { return cloud[i]; }, lo, hi);
    for (int d = 0; d < 3; ++d)
        if (!((hi[d] - lo[d]) / Size < 1e15))
            throw std::invalid_argument(
                "meshio++: filter_points: the voxel size is too small for the extent of the cloud");
    std::vector<FpVoxelRow> rows(m);
    parallel_for(m, [&](std::size_t i) {
        rows[i] = FpVoxelRow{static_cast<std::int64_t>(std::floor((cloud[i][0] - lo[0]) / Size)),
                             static_cast<std::int64_t>(std::floor((cloud[i][1] - lo[1]) / Size)),
                             static_cast<std::int64_t>(std::floor((cloud[i][2] - lo[2]) / Size)),
                             static_cast<std::int64_t>(i)};
    });
    std::sort(rows.begin(), rows.end());
    FpIds out;
    for (std::size_t a = 0; a < m;) {
        std::size_t b = a;
        Vec3 mean{0.0, 0.0, 0.0};
        while (b < m && rows[b].SameVoxel(rows[a])) {
            mean = detail::vec3_add(mean, cloud[static_cast<std::size_t>(rows[b].mId)]);
            ++b;
        }
        mean = detail::vec3_scale(mean, 1.0 / static_cast<double>(b - a));
        std::int64_t best = rows[a].mId;
        double best_d2 = std::numeric_limits<double>::infinity();
        for (std::size_t j = a; j < b; ++j) {
            const double d2 =
                detail::vec3_norm_sq(detail::vec3_sub(cloud[static_cast<std::size_t>(rows[j].mId)], mean));
            if (d2 < best_d2) {  // ascending id: a tie keeps the lower id
                best_d2 = d2;
                best = rows[j].mId;
            }
        }
        out.push_back(rAlive[static_cast<std::size_t>(best)]);
        a = b;
    }
    return out;
}

// Step 4: greedy Poisson-disk thinning in index order.
FpIds fp_poisson(const std::vector<Vec3>& rXyz, const FpIds& rAlive, double MinDistance) {
    const std::size_t m = rAlive.size();
    if (m == 0)
        return rAlive;
    const std::vector<Vec3> cloud = fp_cloud(rXyz, rAlive);
    std::vector<detail::GridKey> keys(m);
    for (std::size_t i = 0; i < m; ++i)
        keys[i] = detail::bucket_key(cloud[i].data(), MinDistance);
    const detail::BucketTable grid = detail::build_bucket_table(MinDistance, keys, keys);
    const double d2_min = MinDistance * MinDistance;
    std::vector<std::uint8_t> kept(m, 0);
    FpIds out;
    for (std::size_t i = 0; i < m; ++i) {
        const detail::GridKey k = keys[i];
        const detail::GridKey lo{k.x - 1, k.y - 1, k.z - 1};
        const detail::GridKey hi{k.x + 1, k.y + 1, k.z + 1};
        bool blocked = false;
        grid.ForEachInBox(lo, hi, [&](const detail::BucketView& rIds) {
            if (blocked)
                return;
            for (const std::int64_t j : rIds) {
                if (!kept[static_cast<std::size_t>(j)])
                    continue;
                if (detail::vec3_norm_sq(detail::vec3_sub(cloud[i], cloud[static_cast<std::size_t>(j)])) <
                    d2_min) {
                    blocked = true;
                    return;
                }
            }
        });
        if (!blocked) {
            kept[i] = 1;
            out.push_back(rAlive[i]);
        }
    }
    return out;
}

}  // namespace

PointFilterResult filter_points(const Mesh& rMesh, const PointFilterOptions& rOptions) {
    static const char* const kPrefix = "meshio++: filter_points: ";
    const auto bad = [&](const char* pWhat) {
        throw std::invalid_argument(std::string(kPrefix) + pWhat);
    };
    if (rOptions.mStatisticalK < 0)
        bad("statistical_k must be zero or above");
    if (!std::isfinite(rOptions.mStdRatio) || rOptions.mStdRatio < 0.0)
        bad("std_ratio must be a finite number, zero or above");
    if (!std::isfinite(rOptions.mRadius) || rOptions.mRadius < 0.0)
        bad("radius must be a finite number, zero or above");
    if (rOptions.mMinNeighbors < 0)
        bad("min_neighbors must be zero or above");
    if (!std::isfinite(rOptions.mVoxelSize) || rOptions.mVoxelSize < 0.0)
        bad("voxel_size must be a finite number, zero or above");
    if (!std::isfinite(rOptions.mMinDistance) || rOptions.mMinDistance < 0.0)
        bad("min_distance must be a finite number, zero or above");
    const bool do_stat = rOptions.mStatisticalK > 0;
    const bool do_radius = rOptions.mRadius > 0.0;
    const bool do_voxel = rOptions.mVoxelSize > 0.0;
    const bool do_poisson = rOptions.mMinDistance > 0.0;
    if (!do_stat && !do_radius && !do_voxel && !do_poisson)
        bad("no step is enabled (give statistical_k, radius, voxel_size or min_distance)");

    const std::size_t n = rMesh.NumPoints();
    std::vector<Vec3> xyz(n);
    for (std::size_t i = 0; i < n; ++i) {
        xyz[i] = detail::read_point(rMesh.Points(), rMesh.PointDim(), static_cast<std::int64_t>(i));
        if (!std::isfinite(xyz[i][0] + xyz[i][1] + xyz[i][2]))
            throw std::invalid_argument(std::string(kPrefix) + "point " + std::to_string(i) +
                                        " is not finite");
    }

    PointFilterResult result;
    result.mNumInput = static_cast<std::int64_t>(n);
    FpIds alive(n);
    std::iota(alive.begin(), alive.end(), std::int64_t{0});
    if (do_stat) {
        const std::size_t before = alive.size();
        alive = fp_statistical(xyz, alive, rOptions.mStatisticalK, rOptions.mStdRatio,
                               result.mWarnings);
        result.mNumRemovedStatistical = static_cast<std::int64_t>(before - alive.size());
    }
    if (do_radius) {
        const std::size_t before = alive.size();
        alive = fp_radius(xyz, alive, rOptions.mRadius, rOptions.mMinNeighbors);
        result.mNumRemovedRadius = static_cast<std::int64_t>(before - alive.size());
    }
    if (do_voxel) {
        const std::size_t before = alive.size();
        alive = fp_voxel(xyz, alive, rOptions.mVoxelSize);
        result.mNumRemovedVoxel = static_cast<std::int64_t>(before - alive.size());
    }
    if (do_poisson) {
        const std::size_t before = alive.size();
        alive = fp_poisson(xyz, alive, rOptions.mMinDistance);
        result.mNumRemovedPoisson = static_cast<std::int64_t>(before - alive.size());
    }

    // The point cloud: points and per-point arrays gathered, one vertex block.
    Mesh out;
    detail::subset_gather_points(rMesh, alive, out);
    const std::size_t kept = alive.size();
    if (kept > 0) {
        NDArray conn = NDArray::Uninit(DType::Int64, {kept, std::size_t{1}});
        for (std::size_t i = 0; i < kept; ++i)
            conn.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
        out.AddCellBlock("vertex", std::move(conn));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    // Point regions follow the surviving points; the rest cannot.
    std::vector<std::int64_t> inverse(n, -1);
    for (std::size_t j = 0; j < kept; ++j)
        inverse[static_cast<std::size_t>(alive[j])] = static_cast<std::int64_t>(j);
    std::size_t dropped_regions = 0;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& region = rMesh.Region(r);
        if (region.mKind != RegionKind::Point) {
            ++dropped_regions;
            continue;
        }
        std::vector<std::int64_t> mapped;
        const std::int64_t* entries = region.Entries();
        for (std::size_t e = 0; e < region.NumEntries(); ++e) {
            const std::int64_t p = entries[e];
            if (p >= 0 && static_cast<std::size_t>(p) < n && inverse[static_cast<std::size_t>(p)] >= 0)
                mapped.push_back(inverse[static_cast<std::size_t>(p)]);
        }
        NDArray arr = NDArray::Uninit(DType::Int64, {mapped.size()});
        std::copy(mapped.begin(), mapped.end(), arr.As<std::int64_t>());
        out.AddRegion(Region(region.mName, RegionKind::Point, region.mDim, region.mTag, std::move(arr)));
    }
    std::size_t input_cells = 0;
    for (const auto cb : rMesh.CellRange())
        input_cells += cb.NumCells();
    if (input_cells > 0)
        result.mWarnings.push_back(std::to_string(input_cells) +
                                   " cell(s) dropped: the result is a point cloud");
    if (dropped_regions > 0)
        result.mWarnings.push_back(std::to_string(dropped_regions) +
                                   " cell or side region(s) dropped: the result is a point cloud");
    result.mMesh = std::move(out);
    result.mPointMap = std::move(alive);
    return result;
}

}  // namespace meshioplusplus
