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
// Cross-mesh field transfer: sample a source mesh's data arrays onto a target
// mesh (nearest source point, or barycentric in a simplexified source), built
// entirely through the uniform mesh API so it compiles under every mesh
// backend. See operations/interpolate.hpp for the contract, in particular the
// determinism rules the pure-numpy twin (_interpolate.py) replicates
// expression for expression.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/bucket_table.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/bounding_box_utils.hpp"
#include "../detail/bucket_build.hpp"
#include "../detail/cell_centroids.hpp"
#include "../detail/point_locator.hpp"
#include "../detail/transfer_names.hpp"
#include "../detail/complex_pairs.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

// Total cell count across all blocks (the global cell numbering size).
std::size_t interp_total_cells(const Mesh& rMesh) {
    std::size_t total = 0;
    for (const auto cb : rMesh.CellRange())
        total += cb.NumCells();
    return total;
}

// The grid cell size for `n` items spread over the bbox `[pLo, pHi]`: the
// largest axis extent divided by the smallest integer R with R*R*R >= n.
// Deliberately no cbrt/pow — the integer loop is exact in both C++ and Python,
// and the single division is IEEE-identical, which is what keeps the numpy
// fallback's grid byte-compatible.
double interp_cell_size(const double* pLo, const double* pHi, std::size_t n) {
    return detail::point_grid_cell_size(pLo, pHi, n);
}

// The point grid of `n` coordinates, id `i` for coordinate `i`: keys filled
// in parallel, then grouped by `build_bucket_table` -- every bucket ascending,
// as serial inserts in id order left the map this replaced (roadmap 5.4.1).
template <class TCoord>
detail::BucketTable interp_point_grid(double Cell, std::size_t n, TCoord&& coord) {
    std::vector<detail::GridKey> keys(n);
    parallel_for(n, [&](std::size_t i) {
        const detail::Vec3 c = coord(i);
        keys[i] =
            detail::GridKey{detail::grid_quantize(c[0], Cell), detail::grid_quantize(c[1], Cell),
                            detail::grid_quantize(c[2], Cell)};
    });
    return detail::build_bucket_table(Cell, keys, keys);
}

// Nearest id in `rGrid` to `rP` by expanding-shell search: any candidate in a
// cell at Chebyshev key radius r >= 1 lies at least (r-1)*cell away, so once
// ((r-1)*cell)^2 exceeds the best squared distance no farther shell can win;
// the occupied-key bbox bounds the expansion (the grid is never empty here).
// Tie-break: smallest squared distance, then lowest id — enforced by the
// strict comparison plus the ascending id order within every bucket.
template <class TCoord>
std::int64_t interp_nearest(const detail::BucketTable& rGrid, TCoord&& coord,
                            const detail::Vec3& rP) {
    const detail::GridKey k = rGrid.KeyOf(rP.data());
    const detail::GridKey& klo = rGrid.OccupiedLo();
    const detail::GridKey& khi = rGrid.OccupiedHi();
    auto axis_cap = [](std::int64_t v, std::int64_t lo, std::int64_t hi) {
        const std::int64_t a = v > lo ? v - lo : lo - v;
        const std::int64_t b = v > hi ? v - hi : hi - v;
        return a > b ? a : b;
    };
    std::int64_t rcap = axis_cap(k.x, klo.x, khi.x);
    rcap = std::max(rcap, axis_cap(k.y, klo.y, khi.y));
    rcap = std::max(rcap, axis_cap(k.z, klo.z, khi.z));
    // Shells below the Chebyshev distance to the occupied bbox are provably
    // empty; starting there skips them without affecting any candidate.
    auto axis_start = [](std::int64_t v, std::int64_t lo, std::int64_t hi) {
        if (v < lo)
            return lo - v;
        if (v > hi)
            return v - hi;
        return std::int64_t{0};
    };
    std::int64_t r0 = axis_start(k.x, klo.x, khi.x);
    r0 = std::max(r0, axis_start(k.y, klo.y, khi.y));
    r0 = std::max(r0, axis_start(k.z, klo.z, khi.z));

    const double cell = rGrid.CellSize();
    double best_d2 = std::numeric_limits<double>::infinity();
    std::int64_t best = -1;
    for (std::int64_t r = r0;; ++r) {
        if (best >= 0 && r >= 1) {
            const double ring = static_cast<double>(r - 1) * cell;
            if (ring * ring > best_d2)
                break;
        }
        if (r > rcap)
            break;
        rGrid.ForEachInShell(k, r, [&](const detail::BucketView& rIds) {
            for (std::int64_t id : rIds) {
                const detail::Vec3 q = coord(id);
                double d2 = 0.0;
                for (int d = 0; d < 3; ++d) {
                    const double delta = rP[d] - q[d];
                    d2 += delta * delta;
                }
                if (d2 < best_d2 || (d2 == best_d2 && id < best)) {
                    best_d2 = d2;
                    best = id;
                }
            }
        });
    }
    return best;
}

// The output shape for `nrows` sampled rows mirroring `rSrc`'s trailing
// dimensions (a 1-D scalar array stays 1-D, the merge convention).
std::vector<std::size_t> interp_out_shape(const NDArray& rSrc, std::size_t nrows) {
    std::vector<std::size_t> shape = rSrc.Shape();
    if (shape.empty())
        shape = {nrows};
    else
        shape[0] = nrows;
    return shape;
}

// Gather rows of `rSrc` by index, dtype-preserving (a raw byte gather — never
// through double, so int64 payloads beyond 2^53 survive bit-for-bit).
NDArray interp_gather_rows(const NDArray& rSrc, const std::vector<std::int64_t>& rIdx) {
    NDArray dst = NDArray::Uninit(rSrc.Dtype(), interp_out_shape(rSrc, rIdx.size()));
    const std::size_t src_rows = detail::rows(rSrc);
    if (rIdx.empty() || src_rows == 0)
        return dst;
    const std::size_t rb = rSrc.Nbytes() / src_rows;
    const std::byte* src = rSrc.Data();
    std::byte* out = dst.Data();
    parallel_for_bw(rIdx.size(), [&](std::size_t n) {
        std::memcpy(out + n * rb, src + static_cast<std::size_t>(rIdx[n]) * rb, rb);
    });
    return dst;
}

}  // namespace

InterpolateMethod interpolate_method_from_name(const std::string& rName) {
    if (rName == "nearest")
        return InterpolateMethod::Nearest;
    if (rName == "barycentric")
        return InterpolateMethod::Barycentric;
    throw std::invalid_argument("meshio++: interpolate: unknown method '" + rName +
                                "' (expected 'nearest' or 'barycentric')");
}

InterpolateConflict interpolate_conflict_from_name(const std::string& rName) {
    return detail::transfer_conflict_from_name<InterpolateConflict>("interpolate", rName);
}

Mesh interpolate(const Mesh& rSource, const Mesh& rTarget, const InterpolateOptions& rOptions) {
    if (rSource.NumPoints() == 0)
        throw std::invalid_argument("meshio++: interpolate: the source mesh has no points");

    // --- resolve which arrays transfer (validation before any work) ---------
    std::vector<std::string> point_names;
    std::vector<std::string> cell_names;
    if (rOptions.mArrays.empty()) {
        point_names = rSource.PointDataNames();  // sorted on every backend
    } else {
        for (const std::string& name : rOptions.mArrays) {
            bool any = false;
            if (rSource.HasPointData(name)) {
                point_names.push_back(name);
                any = true;
            }
            if (rSource.HasCellData(name)) {
                cell_names.push_back(name);
                any = true;
            }
            if (!any) {
                std::string re_name;
                std::string im_name;
                if (detail::find_complex_pair_for_base(rSource, DataLocation::Point, name, re_name,
                                                       im_name)) {
                    point_names.push_back(re_name);
                    point_names.push_back(im_name);
                    any = true;
                }
                if (detail::find_complex_pair_for_base(rSource, DataLocation::Cell, name, re_name,
                                                       im_name)) {
                    cell_names.push_back(re_name);
                    cell_names.push_back(im_name);
                    any = true;
                }
            }
            if (!any)
                throw std::invalid_argument(
                    detail::transfer_unknown_array_message("interpolate", rSource, name));
        }
    }

    for (const std::string& name : cell_names) {
        if (rSource.CellDataNumBlocks(name) != rSource.NumCellBlocks())
            throw std::invalid_argument("meshio++: interpolate: source cell_data '" + name +
                                        "' does not have one array per cell block");
        const std::size_t nblocks = rSource.CellDataNumBlocks(name);
        for (std::size_t b = 1; b < nblocks; ++b)
            if (data_num_components(rSource.CellData(name, b)) !=
                data_num_components(rSource.CellData(name, 0)))
                throw std::invalid_argument("meshio++: interpolate: source cell_data '" + name +
                                            "' has inconsistent components across blocks");
    }

    std::vector<std::string> point_out(point_names.size());
    for (std::size_t i = 0; i < point_names.size(); ++i)
        point_out[i] = detail::transfer_resolve_name("interpolate", rTarget, true, point_names[i],
                                                     rOptions.mOnConflict);
    std::vector<std::string> cell_out(cell_names.size());
    for (std::size_t i = 0; i < cell_names.size(); ++i)
        cell_out[i] = detail::transfer_resolve_name("interpolate", rTarget, false, cell_names[i],
                                                    rOptions.mOnConflict);

    // KRATOS's Mesh is not copy-constructible — clone through the uniform API.
    Mesh out = detail::clone_mesh(rTarget);

    const NDArray& tpts = rTarget.Points();
    const std::size_t tdim = rTarget.PointDim();
    const std::size_t nt = rTarget.NumPoints();

    // --- point_data: sampled at the target's points -------------------------
    if (!point_names.empty()) {
        if (rOptions.mMethod == InterpolateMethod::Nearest) {
            const NDArray& spts = rSource.Points();
            const std::size_t sdim = rSource.PointDim();
            const std::size_t ns = rSource.NumPoints();
            auto scoord = [&](std::int64_t i) { return detail::read_point(spts, sdim, i); };
            double lo[3], hi[3];
            detail::BoundingBoxUtils::FromCoords(
                ns, [&](std::size_t i) { return scoord(static_cast<std::int64_t>(i)); }, lo, hi);
            const detail::BucketTable grid = interp_point_grid(
                interp_cell_size(lo, hi, ns), ns,
                [&](std::size_t i) { return scoord(static_cast<std::int64_t>(i)); });

            std::vector<std::int64_t> nearest(nt);
            parallel_for(nt, [&](std::size_t t) {
                const detail::Vec3 p = detail::read_point(tpts, tdim, static_cast<std::int64_t>(t));
                nearest[t] = interp_nearest(grid, scoord, p);
            });
            for (std::size_t i = 0; i < point_names.size(); ++i)
                out.AddPointData(point_out[i],
                                 interp_gather_rows(rSource.PointData(point_names[i]), nearest));
        } else {
            // Simplexify so containment/weights are exact and uniform; points
            // are untouched for a linear source, and point_data rides through
            // (row-subset under the higher-order linearization), so sampling
            // the simplexified mesh's own arrays is always index-consistent.
            ConvertCellsOptions cco;
            cco.mMode = ConvertCellsMode::Simplexify;
            const ConvertCellsResult simp = convert_cells(rSource, cco);
            const Mesh& smesh = simp.mMesh;

            // The search structure (flat simplex connectivity in ascending
            // block-then-cell order, the quantised-box grid, the first
            // containing simplex in ascending index): detail/point_locator.hpp,
            // shared with locate_points.
            const detail::PointLocator locator(smesh, detail::kPointLocatorTolerance,
                                               "interpolate");
            const std::size_t nv = locator.NodesPerSimplex();
            const std::vector<std::int64_t>& sconn = locator.Connectivity();

            const NDArray& spts = smesh.Points();
            const std::size_t sdim = smesh.PointDim();
            const std::size_t ns = smesh.NumPoints();
            auto scoord = [&](std::int64_t i) { return detail::read_point(spts, sdim, i); };
            const double* lo = locator.Lo();
            const double* hi = locator.Hi();

            std::vector<std::int64_t> simplex_of(nt, -1);
            std::vector<std::array<double, 4>> wts(nt);
            parallel_for(nt, [&](std::size_t t) {
                const detail::Vec3 p = detail::read_point(tpts, tdim, static_cast<std::int64_t>(t));
                double w[4] = {0.0, 0.0, 0.0, 0.0};
                const std::int64_t s = locator.Locate(p, w);
                if (s < 0)
                    return;
                simplex_of[t] = s;
                for (std::size_t k = 0; k < nv; ++k)
                    wts[t][k] = w[k];
            });

            std::size_t n_outside = 0;
            for (std::size_t t = 0; t < nt; ++t)
                if (simplex_of[t] < 0)
                    ++n_outside;

            std::vector<std::int64_t> near_of;
            if (n_outside > 0 && rOptions.mExtrapolate) {
                const detail::BucketTable pgrid = interp_point_grid(
                    interp_cell_size(lo, hi, ns), ns,
                    [&](std::size_t i) { return scoord(static_cast<std::int64_t>(i)); });
                near_of.assign(nt, -1);
                parallel_for(nt, [&](std::size_t t) {
                    if (simplex_of[t] >= 0)
                        return;
                    const detail::Vec3 p =
                        detail::read_point(tpts, tdim, static_cast<std::int64_t>(t));
                    near_of[t] = interp_nearest(pgrid, scoord, p);
                });
            } else if (n_outside > 0) {
                log::warn(
                    "meshio++: interpolate: {} of {} target points lie outside the source "
                    "domain; filled with default_value = {}",
                    n_outside, nt, rOptions.mDefaultValue);
            }

            for (std::size_t i = 0; i < point_names.size(); ++i) {
                const NDArray& src = smesh.PointData(point_names[i]);
                const std::size_t comps = data_num_components(src);
                NDArray arr = NDArray::Uninit(DType::Float64, interp_out_shape(src, nt));
                double* o = arr.As<double>();
                parallel_for(nt, [&](std::size_t t) {
                    if (simplex_of[t] >= 0) {
                        const std::int64_t* nodes =
                            &sconn[static_cast<std::size_t>(simplex_of[t]) * nv];
                        for (std::size_t c = 0; c < comps; ++c) {
                            // Corner-by-corner in connectivity order (fixed —
                            // the numpy twin sums in the same order).
                            double acc = wts[t][0] *
                                         detail::read_double(
                                             src, static_cast<std::size_t>(nodes[0]) * comps + c);
                            for (std::size_t k = 1; k < nv; ++k)
                                acc += wts[t][k] *
                                       detail::read_double(
                                           src, static_cast<std::size_t>(nodes[k]) * comps + c);
                            o[t * comps + c] = acc;
                        }
                    } else if (rOptions.mExtrapolate) {
                        for (std::size_t c = 0; c < comps; ++c)
                            o[t * comps + c] = detail::read_double(
                                src, static_cast<std::size_t>(near_of[t]) * comps + c);
                    } else {
                        for (std::size_t c = 0; c < comps; ++c)
                            o[t * comps + c] = rOptions.mDefaultValue;
                    }
                });
                out.AddPointData(point_out[i], std::move(arr));
            }
        }
    }

    // --- cell_data: always nearest source-cell centroid ---------------------
    if (!cell_names.empty()) {
        const std::size_t src_total = interp_total_cells(rSource);
        if (src_total == 0)
            throw std::invalid_argument(
                "meshio++: interpolate: cell_data was requested but the source mesh has no cells");
        if (rTarget.NumCellBlocks() == 0) {
            log::warn(
                "meshio++: interpolate: cell_data requested but the target has no cell blocks; "
                "skipping {} array(s)",
                cell_names.size());
        } else {
            const std::size_t tgt_total = interp_total_cells(rTarget);
            const std::vector<double> scent = detail::cell_centroids(rSource, src_total);
            const std::vector<double> tcent = detail::cell_centroids(rTarget, tgt_total);
            auto ccoord = [&](std::int64_t i) {
                const std::size_t j = static_cast<std::size_t>(i);
                return detail::Vec3{scent[3 * j], scent[3 * j + 1], scent[3 * j + 2]};
            };
            double lo[3], hi[3];
            detail::BoundingBoxUtils::FromCoords(
                src_total, [&](std::size_t i) { return ccoord(static_cast<std::int64_t>(i)); }, lo,
                hi);
            const detail::BucketTable grid = interp_point_grid(
                interp_cell_size(lo, hi, src_total), src_total,
                [&](std::size_t i) { return ccoord(static_cast<std::int64_t>(i)); });

            std::vector<std::int64_t> nearest_cell(tgt_total);
            parallel_for(tgt_total, [&](std::size_t t) {
                const detail::Vec3 p{tcent[3 * t], tcent[3 * t + 1], tcent[3 * t + 2]};
                nearest_cell[t] = interp_nearest(grid, ccoord, p);
            });

            // Source block offsets in the global cell numbering.
            std::vector<std::size_t> soff;
            soff.reserve(rSource.NumCellBlocks() + 1);
            soff.push_back(0);
            for (const auto cb : rSource.CellRange())
                soff.push_back(soff.back() + cb.NumCells());

            for (std::size_t i = 0; i < cell_names.size(); ++i) {
                const std::string& name = cell_names[i];
                const std::size_t nblocks = rSource.CellDataNumBlocks(name);
                const std::size_t comps = data_num_components(rSource.CellData(name, 0));
                DType dt = rSource.CellData(name, 0).Dtype();
                for (std::size_t b = 1; b < nblocks; ++b)
                    if (rSource.CellData(name, b).Dtype() != dt)
                        dt = DType::Float64;  // mixed dtypes widen

                std::vector<NDArray> blocks_out;
                blocks_out.reserve(rTarget.NumCellBlocks());
                std::size_t tbase = 0;
                for (const auto cb : rTarget.CellRange()) {
                    const std::size_t ncells = cb.NumCells();
                    NDArray arr =
                        NDArray::Uninit(dt, interp_out_shape(rSource.CellData(name, 0), ncells));
                    parallel_for_bw(ncells, [&, tbase](std::size_t c) {
                        const std::size_t g = static_cast<std::size_t>(nearest_cell[tbase + c]);
                        std::size_t sb = 0;
                        while (g >= soff[sb + 1])
                            ++sb;
                        const std::size_t sl = g - soff[sb];
                        const NDArray& src = rSource.CellData(name, sb);
                        if (src.Dtype() == dt) {
                            const std::size_t rb = comps * dtype_size(dt);
                            std::memcpy(arr.Data() + c * rb, src.Data() + sl * rb, rb);
                        } else {
                            for (std::size_t k = 0; k < comps; ++k)
                                detail::store_double_truncating(
                                    arr, c * comps + k, detail::read_double(src, sl * comps + k));
                        }
                    });
                    blocks_out.push_back(std::move(arr));
                    tbase += ncells;
                }
                out.AddCellData(cell_out[i], std::move(blocks_out));
            }
        }
    }

    return out;
}

}  // namespace meshioplusplus
