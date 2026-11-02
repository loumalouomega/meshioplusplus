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
 * @file detail/point_knn.hpp
 * @brief Exact nearest and k-nearest queries of arbitrary points against a
 * fixed point cloud, over `detail/spatial_hash.hpp`'s bucket grid. Hoisted
 * from `register_rigid`'s nearest-point search (v18.x), which it serves
 * unchanged with `k = 1`, for `reconstruct_surface`'s tangent planes.
 *
 * The answer is exact and deterministic: candidates are ordered by squared
 * distance (`dx*dx + dy*dy + dz*dz`, in that order), ties to the lower id,
 * and the shells expand until the next shell's lower bound exceeds the k-th
 * best. The cell size is a speed choice, never the answer. Queries are const
 * and may run concurrently.
 */

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/spatial_hash.hpp"

namespace meshioplusplus {
namespace detail {

class PointKnn {
public:
    using Vec = std::array<double, 3>;

    PointKnn() = default;

    /// Index (a copy of) @p Points with buckets of side @p CellSize (positive).
    PointKnn(std::vector<Vec> Points, double CellSize)
        : mPoints(std::move(Points)), mGrid(CellSize) {
        for (std::size_t i = 0; i < mPoints.size(); ++i)
            mGrid.Insert(mGrid.KeyOf(mPoints[i].data()), static_cast<std::int64_t>(i));
    }

    /// The indexed points.
    const std::vector<Vec>& Points() const { return mPoints; }

    /// The nearest point's id (ties to the lower id), or -1 for an empty cloud.
    std::int64_t Nearest(const Vec& rP) const {
        std::pair<double, std::int64_t> best{0.0, -1};
        KNearest(rP, 1, &best);
        return best.second;
    }

    /**
     * @brief The @p K nearest points to @p rP, ascending by (squared
     * distance, id), into @p pOut (room for @p K). Returns how many were
     * found (fewer than @p K only when the cloud is smaller).
     */
    std::size_t KNearest(const Vec& rP, std::size_t K,
                         std::pair<double, std::int64_t>* pOut) const {
        if (K == 0 || mPoints.empty())
            return 0;
        const GridKey k = mGrid.KeyOf(rP.data());
        const GridKey& lo = mGrid.OccupiedLo();
        const GridKey& hi = mGrid.OccupiedHi();
        const auto cap = [](std::int64_t v, std::int64_t a, std::int64_t b) {
            return std::max(v > a ? v - a : a - v, v > b ? v - b : b - v);
        };
        const auto start = [](std::int64_t v, std::int64_t a, std::int64_t b) {
            return v < a ? a - v : (v > b ? v - b : std::int64_t{0});
        };
        const std::int64_t rcap =
            std::max({cap(k.x, lo.x, hi.x), cap(k.y, lo.y, hi.y), cap(k.z, lo.z, hi.z)});
        const std::int64_t r0 =
            std::max({start(k.x, lo.x, hi.x), start(k.y, lo.y, hi.y), start(k.z, lo.z, hi.z)});
        const double cell = mGrid.CellSize();
        const std::vector<Vec>& pts = mPoints;
        std::size_t found = 0;
        // pOut[0, found) stays sorted; an insertion sort, K is small.
        const auto offer = [&](double D2, std::int64_t Id) {
            const std::pair<double, std::int64_t> c{D2, Id};
            if (found == K && !(c < pOut[K - 1]))
                return;
            std::size_t at = found < K ? found++ : K - 1;
            while (at > 0 && c < pOut[at - 1]) {
                pOut[at] = pOut[at - 1];
                --at;
            }
            pOut[at] = c;
        };
        for (std::int64_t r = r0; r <= rcap; ++r) {
            if (found == K && r >= 1) {
                const double ring = static_cast<double>(r - 1) * cell;
                if (ring * ring > pOut[K - 1].first)
                    break;
            }
            mGrid.ForEachInShell(k, r, [&](const std::vector<std::int64_t>& rIds) {
                for (const std::int64_t id : rIds) {
                    const Vec& q = pts[static_cast<std::size_t>(id)];
                    const double dx = rP[0] - q[0], dy = rP[1] - q[1], dz = rP[2] - q[2];
                    offer(dx * dx + dy * dy + dz * dz, id);
                }
            });
        }
        return found;
    }

private:
    std::vector<Vec> mPoints;
    SpatialGrid mGrid{1.0};
};

}  // namespace detail
}  // namespace meshioplusplus
