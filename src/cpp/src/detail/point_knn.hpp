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
 * @file detail/point_knn.hpp
 * @brief Exact nearest and k-nearest queries of arbitrary points against a
 * fixed point cloud, over a frozen bucket grid (`detail/bucket_build.hpp`). Hoisted
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
#include "meshioplusplus/detail/bucket_table.hpp"

#include "bucket_build.hpp"

namespace meshioplusplus {
namespace detail {

class PointKnn {
public:
    using Vec = std::array<double, 3>;

    PointKnn() = default;

    /// Index (a copy of) @p Points with buckets of side @p CellSize (positive).
    PointKnn(std::vector<Vec> Points, double CellSize) : mPoints(std::move(Points)) {
        std::vector<GridKey> keys(mPoints.size());
        for (std::size_t i = 0; i < mPoints.size(); ++i)
            keys[i] = bucket_key(mPoints[i].data(), CellSize);
        mGrid = build_bucket_table(CellSize, keys, keys);
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
            mGrid.ForEachInShell(k, r, [&](const BucketView& rIds) {
                for (const std::int64_t id : rIds) {
                    const Vec& q = pts[static_cast<std::size_t>(id)];
                    const double dx = rP[0] - q[0], dy = rP[1] - q[1], dz = rP[2] - q[2];
                    offer(dx * dx + dy * dy + dz * dz, id);
                }
            });
        }
        return found;
    }

    /**
     * @brief How many indexed points lie within @p Radius of @p rP
     * (`dx*dx + dy*dy + dz*dz <= Radius * Radius`, the query point itself
     * included when it is indexed), counting no further than @p Limit.
     *
     * The grid cells of the box `rP +- Radius` are visited, so the cost grows
     * with `(Radius / CellSize)^3`: build the index with a cell size of about
     * the radius being asked about.
     */
    std::size_t CountWithin(const Vec& rP, double Radius,
                            std::size_t Limit = std::numeric_limits<std::size_t>::max()) const {
        if (mPoints.empty() || Limit == 0 || !(Radius >= 0.0))
            return 0;
        const Vec lo{rP[0] - Radius, rP[1] - Radius, rP[2] - Radius};
        const Vec hi{rP[0] + Radius, rP[1] + Radius, rP[2] + Radius};
        const double r2 = Radius * Radius;
        std::size_t count = 0;
        mGrid.ForEachInBox(
            mGrid.KeyOf(lo.data()), mGrid.KeyOf(hi.data()), [&](const BucketView& rIds) {
                for (const std::int64_t id : rIds) {
                    if (count >= Limit)
                        return;
                    const Vec& q = mPoints[static_cast<std::size_t>(id)];
                    const double dx = rP[0] - q[0], dy = rP[1] - q[1], dz = rP[2] - q[2];
                    if (dx * dx + dy * dy + dz * dz <= r2)
                        ++count;
                }
            });
        return count;
    }

    /**
     * @brief All indexed points within @p Radius of @p rP
     * (`dx*dx + dy*dy + dz*dz <= Radius * Radius`), appended to @p rOut as
     * `(squared_distance, id)` pairs sorted deterministically by ascending id.
     */
    void WithinRadius(const Vec& rP, double Radius,
                      std::vector<std::pair<double, std::int64_t>>& rOut) const {
        rOut.clear();
        if (mPoints.empty() || !(Radius >= 0.0))
            return;
        const Vec lo{rP[0] - Radius, rP[1] - Radius, rP[2] - Radius};
        const Vec hi{rP[0] + Radius, rP[1] + Radius, rP[2] + Radius};
        const double r2 = Radius * Radius;
        mGrid.ForEachInBox(
            mGrid.KeyOf(lo.data()), mGrid.KeyOf(hi.data()), [&](const BucketView& rIds) {
                for (const std::int64_t id : rIds) {
                    const Vec& q = mPoints[static_cast<std::size_t>(id)];
                    const double dx = rP[0] - q[0], dy = rP[1] - q[1], dz = rP[2] - q[2];
                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 <= r2)
                        rOut.emplace_back(d2, id);
                }
            });
        std::sort(rOut.begin(), rOut.end(), [](const auto& rA, const auto& rB) {
            return rA.second < rB.second;
        });
    }

private:
    std::vector<Vec> mPoints;
    BucketTable mGrid{1.0};
};

}  // namespace detail
}  // namespace meshioplusplus
