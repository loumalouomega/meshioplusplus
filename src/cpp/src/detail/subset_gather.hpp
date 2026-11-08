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
 * @file detail/subset_gather.hpp
 * @brief The point gather `build_cell_subset` and `clean` share: the points and
 * the point_data of a mesh, gathered by a list of source rows.
 *
 * A **core-private** header. The body is forced inline (`MIO_ALWAYS_INLINE`)
 * into both callers so each compiles as it did with its own copy (roadmap
 * 3.4.11): as a plain inline function GCC kept it out of line, and the extra
 * function in `subset.cpp` cost `build_cell_subset` its inlining of
 * `read_int` (`split` +10%); an out-of-line entry beside a forced-inline body
 * in `subset.cpp` did the same.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/parallel.hpp"

#include "always_inline.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Assigns @p rOut the points of @p rMesh gathered by @p rPointSrc (output
 * row -> source point), as byte rows keeping the dtype, and adds each point_data
 * array: gathered by the same rows when it has one row per point, an owned
 * verbatim copy otherwise.
 */
MIO_ALWAYS_INLINE void subset_gather_points(const Mesh& rMesh,
                                            const std::vector<std::int64_t>& rPointSrc,
                                            Mesh& rOut) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();
    {
        NDArray newpts = NDArray::Uninit(points.Dtype(), {rPointSrc.size(), dim});
        if (!rPointSrc.empty() && dim > 0) {
            const std::size_t rb = dim * dtype_size(points.Dtype());
            const std::byte* s = points.Data();
            std::byte* d = newpts.Data();
            parallel_for_bw(rPointSrc.size(), [&](std::size_t j) {
                std::memcpy(d + j * rb, s + static_cast<std::size_t>(rPointSrc[j]) * rb, rb);
            });
        }
        rOut.AssignPoints(std::move(newpts));
    }
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        const std::size_t rows = a.Shape().empty() ? 0 : a.Shape()[0];
        if (rows == n)
            rOut.AddPointData(name, subset_gather_rows(a, rPointSrc));
        else {
            NDArray c = a;
            c.MakeOwned();
            rOut.AddPointData(name, std::move(c));
        }
    }
}

}  // namespace detail
}  // namespace meshioplusplus
