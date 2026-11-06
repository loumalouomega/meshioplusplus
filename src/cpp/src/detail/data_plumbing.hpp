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
 * @file detail/data_plumbing.hpp
 * @brief Row plumbing shared by the operations that rebuild a mesh from a
 * subset or a refinement of another: gathering the point table, compacting it
 * by a `used` mask, replicating parent cell rows to children, and copying the
 * geometry under new coordinates.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): it lives beside
 * the `.cpp` files, no installed header names it, and it adds nothing to the
 * API or the ABI. Built on the uniform mesh API only, so it works under every
 * backend. Each body is the byte-copy loop its callers used to repeat, moved
 * verbatim, so outputs are unchanged.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/ragged_csr.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Assigns @p rOut the points of @p rMesh whose @p rUsed flag is set, in
 * ascending order (`NumUsed` of them), as a byte-row subset keeping the dtype.
 */
inline void compact_points_by_mask(const Mesh& rMesh, const std::vector<char>& rUsed,
                                   std::size_t NumUsed, Mesh& rOut) {
    const std::size_t num_points = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t dim = detail::cols(points);
    const std::size_t row_bytes = dim * dtype_size(points.Dtype());
    NDArray new_points = NDArray::Uninit(points.Dtype(), {NumUsed, dim});
    const std::byte* src = points.Data();
    std::byte* dst = new_points.Data();
    std::size_t w = 0;
    for (std::size_t i = 0; i < num_points; ++i)
        if (rUsed[i])
            std::memcpy(dst + (w++) * row_bytes, src + i * row_bytes, row_bytes);
    rOut.AssignPoints(std::move(new_points));
}

/**
 * @brief Adds to @p rOut each `point_data` array of @p rMesh subset to the rows
 * whose @p rUsed flag is set (`NumUsed` of them), keeping the dtype.
 *
 * An array whose row count is not the point count is not per-point data: with
 * @p KeepMismatched it is copied verbatim (owned), otherwise it is dropped.
 */
inline void compact_point_data_by_mask(const Mesh& rMesh, const std::vector<char>& rUsed,
                                       std::size_t NumUsed, bool KeepMismatched, Mesh& rOut) {
    const std::size_t num_points = rMesh.NumPoints();
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (detail::rows(a) != num_points) {
            if (KeepMismatched)
                rOut.AddPointData(name, detail::data_owned_copy(a));
            continue;
        }
        std::vector<std::size_t> shape = a.Shape();
        shape[0] = NumUsed;
        const std::size_t row_bytes = num_points == 0 ? 0 : a.Nbytes() / num_points;
        NDArray b = NDArray::Uninit(a.Dtype(), std::move(shape));
        const std::byte* src = a.Data();
        std::byte* dst = b.Data();
        std::size_t w = 0;
        for (std::size_t i = 0; i < num_points; ++i)
            if (rUsed[i])
                std::memcpy(dst + (w++) * row_bytes, src + i * row_bytes, row_bytes);
        rOut.AddPointData(name, std::move(b));
    }
}

/**
 * @brief The blocks of cell_data array @p rName of @p rMesh with each parent's
 * row replicated to its children (`rParentOfChild[b][i]` is child `i`'s parent
 * in block `b`), by raw byte copy (dtype-agnostic, no precision loss).
 *
 * A block is copied verbatim (owned) instead when the array's block count is
 * not `NumBlocks` (it is not per-cell data), when it has no rows, or when
 * `rCopyVerbatim(b, in_rows)` says so -- the caller's own skip test.
 */
template <class TCopyVerbatim>
std::vector<NDArray> replicate_parent_rows(
    const Mesh& rMesh, const std::string& rName, std::size_t NumBlocks,
    const std::vector<std::vector<std::int64_t>>& rParentOfChild, TCopyVerbatim&& rCopyVerbatim) {
    const std::size_t ndata = rMesh.CellDataNumBlocks(rName);
    std::vector<NDArray> blocks;
    blocks.reserve(ndata);
    for (std::size_t b = 0; b < ndata; ++b) {
        const NDArray& a = rMesh.CellData(rName, b);
        const std::size_t in_rows = detail::rows(a);
        if (ndata != NumBlocks || in_rows == 0 || rCopyVerbatim(b, in_rows)) {
            blocks.push_back(detail::data_owned_copy(a));
            continue;
        }
        const std::vector<std::int64_t>& parents = rParentOfChild[b];
        std::vector<std::size_t> shape = a.Shape();
        shape[0] = parents.size();
        const std::size_t row_bytes = a.Nbytes() / in_rows;
        NDArray replicated = NDArray::Uninit(a.Dtype(), std::move(shape));
        const std::byte* src = a.Data();
        std::byte* dst = replicated.Data();
        parallel_for_bw(parents.size(), [&](std::size_t i) {
            std::memcpy(dst + i * row_bytes, src + static_cast<std::size_t>(parents[i]) * row_bytes,
                        row_bytes);
        });
        blocks.push_back(std::move(replicated));
    }
    return blocks;
}

/**
 * @brief Fills the empty mesh @p rOut with @p Points and owned copies of every
 * cell block (ragged included), named region and property set of @p rMesh --
 * `clone_geometry` with new coordinates. No data arrays are carried.
 */
inline void clone_geometry_into(const Mesh& rMesh, NDArray Points, Mesh& rOut) {
    rOut.AssignPoints(std::move(Points));
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsRagged()) {
            detail::append_ragged_copy(cb, rOut);
        } else {
            rOut.AddCellBlock(std::string(cb.Type()), data_owned_copy(cb.Conn()));
        }
    }
    // Named regions ride through verbatim. Every caller leaves the point and
    // cell numbering untouched (the data operations, `transform`, `smooth`,
    // `interpolate`'s target clone, `attach_quality`), so there is nothing to
    // remap — see doc/regions.md's table of which operations remap and which
    // pass through.
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        rOut.AddRegion(rMesh.Region(i));
    // Property sets ride along too: these operations preserve the mesh's shape,
    // so dropping a deck's material data here would be new lossiness. They are
    // keyed by id, not by entity index, so there is nothing to remap.
    for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i)
        rOut.AddPropertySet(rMesh.GetPropertySet(i));
}

}  // namespace detail
}  // namespace meshioplusplus
