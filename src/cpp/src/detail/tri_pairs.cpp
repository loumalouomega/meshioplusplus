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
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

// Project includes
#include "tri_pairs.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus::detail {

namespace {

using TpBox = std::array<double, 6>;

std::vector<TpBox> tp_boxes(const TriangleSoup& rSoup) {
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<TpBox> boxes(rSoup.NumTriangles());
    parallel_for(boxes.size(), [&](std::size_t t) {
        TpBox& r_box = boxes[t];
        r_box = {{inf, inf, inf, -inf, -inf, -inf}};
        for (std::size_t k = 0; k < 3; ++k) {
            const Vec3& p = rSoup.mCorners[3 * t + k];
            for (std::size_t d = 0; d < 3; ++d) {
                r_box[d] = std::min(r_box[d], p[d]);
                r_box[3 + d] = std::max(r_box[3 + d], p[d]);
            }
        }
    });
    return boxes;
}

bool tp_apart(const TpBox& rA, const TpBox& rB) {
    return rA[0] > rB[3] || rB[0] > rA[3] || rA[1] > rB[4] || rB[1] > rA[4] || rA[2] > rB[5] ||
           rB[2] > rA[5];
}

}  // namespace

void soup_box_pairs(const TriangleSoup& rA, const TriangleSoup* pB, const TriPairVisitor& rVisit) {
    const std::vector<TpBox> box_a = tp_boxes(rA);
    const std::vector<TpBox> box_b = pB != nullptr ? tp_boxes(*pB) : std::vector<TpBox>{};
    const std::vector<TpBox>& r_grid_boxes = pB != nullptr ? box_b : box_a;
    if (box_a.empty() || r_grid_boxes.empty())
        return;

    // The grid's cell from the gridded triangles' mean extent and the diagonal
    // of both soups.
    const double inf = std::numeric_limits<double>::infinity();
    std::array<double, 3> lo{{inf, inf, inf}}, hi{{-inf, -inf, -inf}};
    double extent_sum = 0.0;
    for (const std::vector<TpBox>* p_boxes : {&box_a, &box_b})
        for (const TpBox& r_box : *p_boxes)
            for (std::size_t d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], r_box[d]);
                hi[d] = std::max(hi[d], r_box[3 + d]);
            }
    for (const TpBox& r_box : r_grid_boxes)
        extent_sum += std::max({r_box[3] - r_box[0], r_box[4] - r_box[1], r_box[5] - r_box[2]});
    const double diag =
        std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                  (hi[2] - lo[2]) * (hi[2] - lo[2]));
    double cell = std::max(extent_sum / static_cast<double>(r_grid_boxes.size()), diag / 256.0);
    if (!(cell > 0.0) || !std::isfinite(cell))
        cell = 1.0;
    SpatialGrid grid(cell);
    for (std::size_t t = 0; t < r_grid_boxes.size(); ++t)
        grid.InsertBox(grid.KeyOf(r_grid_boxes[t].data()), grid.KeyOf(r_grid_boxes[t].data() + 3),
                       static_cast<std::int64_t>(t));

    const bool self = pB == nullptr;
    parallel_for(
        box_a.size(),
        [&](std::size_t i) {
            const TpBox& bi = box_a[i];
            std::vector<std::int64_t> near;
            grid.ForEachInBox(grid.KeyOf(bi.data()), grid.KeyOf(bi.data() + 3),
                              [&](const std::vector<std::int64_t>& rIds) {
                                  for (std::int64_t id : rIds)
                                      if (!self || static_cast<std::size_t>(id) > i)
                                          near.push_back(id);
                              });
            std::sort(near.begin(), near.end());
            near.erase(std::unique(near.begin(), near.end()), near.end());
            near.erase(std::remove_if(near.begin(), near.end(),
                                      [&](std::int64_t j) {
                                          return tp_apart(
                                              bi, r_grid_boxes[static_cast<std::size_t>(j)]);
                                      }),
                       near.end());
            if (!near.empty())
                rVisit(i, near);
        },
        64);
}

}  // namespace meshioplusplus::detail
