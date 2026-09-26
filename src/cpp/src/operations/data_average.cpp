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
// Point <-> cell data averaging. The point->cell direction is a pure gather
// (one output row per cell, reading only that cell's nodes) and is therefore
// parallel over cells. The cell->point direction is a scatter, so its
// accumulation pass is deliberately SERIAL: floating-point addition is not
// associative, and a parallel scatter would make the result depend on the
// thread count and on the scheduling, which would flake the cross-backend
// tests. Only the final divide is parallelised. This is the same phase-split
// philosophy as operations/surface.cpp.
//
// Geometry is never modified. See operations/data_average.hpp for the contract.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/slot_runs.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

/// Output shape for `rows` entries of `ncomp` components: 1-D when scalar, so
/// scalar data stays scalar (the same rule merge/clean follow).
std::vector<std::size_t> davg_shape(std::size_t rows, std::size_t ncomp) {
    if (ncomp <= 1)
        return {rows};
    return {rows, ncomp};
}

/// Applies the NaN policy to one produced value.
double davg_apply_nan_policy(double v, const DataAverageOptions& rOpts, const char* pWhat,
                             const std::string& rName, std::size_t index) {
    if (std::isfinite(v))
        return v;
    switch (rOpts.nan_policy) {
        case NanPolicy::Ignore:
            return v;
        case NanPolicy::Replace:
            return rOpts.nan_replacement;
        case NanPolicy::Fail:
            throw std::invalid_argument(
                std::string("meshio++: data_average: non-finite value for ") + pWhat + " " +
                std::to_string(index) + " of '" + rName + "'");
    }
    return v;
}

/// The target name for one input array, per prefix/suffix.
std::string davg_target_name(const std::string& rName, const DataAverageOptions& rOpts) {
    return rOpts.prefix + rName + rOpts.suffix;
}

/// Resolves the requested names, validating each against the source location.
std::vector<std::string> davg_resolve_names(const Mesh& rMesh, DataLocation source,
                                            const DataAverageOptions& rOpts) {
    if (rOpts.names.empty())
        return data_names(rMesh, source);
    for (const std::string& n : rOpts.names)
        if (!data_has(rMesh, source, n))
            throw std::invalid_argument(data_unknown_key_message(rMesh, source, n));
    return rOpts.names;
}

/// Rejects an output name that already exists when overwriting is disabled.
void davg_check_target(const Mesh& rMesh, DataLocation target, const std::string& rTarget,
                       const DataAverageOptions& rOpts) {
    if (!rOpts.overwrite && data_has(rMesh, target, rTarget))
        throw std::invalid_argument(std::string("meshio++: data_average: ") +
                                    data_location_name(target) + " '" + rTarget +
                                    "' already exists (pass overwrite=true to replace it)");
}

/// The node ids of one cell, covering rectangular, polygon and polyhedron
/// blocks. For a polyhedron the *distinct* nodes across all faces are used.
void davg_cell_nodes(const Mesh::CellView& rCell, std::size_t index,
                     std::vector<std::int64_t>& rNodes) {
    rNodes.clear();
    if (rCell.IsPolyhedron()) {
        std::unordered_set<std::int64_t> seen;
        for (std::size_t f = 0; f < rCell.NumFaces(index); ++f) {
            const auto face = rCell.Face(index, f);
            for (std::size_t i = 0; i < face.second; ++i)
                if (seen.insert(face.first[i]).second)
                    rNodes.push_back(face.first[i]);
        }
        return;
    }
    if (rCell.IsRagged()) {
        const std::int64_t* row = rCell.Row(index);
        rNodes.assign(row, row + rCell.RowSize(index));
        return;
    }
    // Connectivity dtype is not guaranteed to be Int64 (a MESHIO-backed mesh
    // often carries Int32 straight from numpy), so read it dtype-agnostically.
    const std::size_t npc = rCell.NodesPerCell();
    const NDArray& conn = rCell.Conn();
    rNodes.reserve(npc);
    for (std::size_t k = 0; k < npc; ++k)
        rNodes.push_back(detail::read_int(conn, index * npc + k));
}

/// Every point's incident cells (global, block-major ids), one entry per
/// occurrence of the point in `davg_cell_nodes`, in (block, cell, node)
/// order: the order in which the serial scatter added their contributions.
/// Entries naming no point of the mesh are dropped, as the scatter skipped
/// them.
struct DavgIncidence {
    std::vector<std::uint64_t> mStart;   ///< `npoints + 1` offsets into `mCellOf`
    std::vector<std::uint64_t> mCellOf;  ///< per entry, its global cell
};

DavgIncidence davg_incidence(const Mesh& rMesh, const std::vector<std::size_t>& rBlockBase,
                             std::size_t NumPoints) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::size_t ncells = rBlockBase[nblocks];
    std::vector<std::uint64_t> count(ncells);
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        parallel_for(cb.NumCells(), [&](std::size_t c) {
            thread_local std::vector<std::int64_t> nodes;
            davg_cell_nodes(cb, c, nodes);
            count[rBlockBase[b] + c] = nodes.size();
        });
    }
    std::vector<std::uint64_t> offset(ncells);
    const std::uint64_t nentries =
        parallel_exclusive_scan(count.data(), ncells, offset.data(), std::uint64_t{0});
    std::vector<std::uint64_t> point_of(nentries);
    std::vector<std::uint64_t> cell_of(nentries);
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        parallel_for(cb.NumCells(), [&](std::size_t c) {
            thread_local std::vector<std::int64_t> nodes;
            davg_cell_nodes(cb, c, nodes);
            const std::size_t gc = rBlockBase[b] + c;
            std::uint64_t e = offset[gc];
            for (std::int64_t node : nodes) {
                point_of[e] = node < 0 || static_cast<std::size_t>(node) >= NumPoints
                                  ? NumPoints
                                  : static_cast<std::uint64_t>(node);
                cell_of[e] = gc;
                ++e;
            }
        });
    }
    // A stable counting sort by point keeps each point's entries in entry
    // order; bucket `NumPoints` collects the dropped ones.
    std::vector<std::uint64_t> order;
    DavgIncidence out;
    detail::slot_runs_impl::counting_sort(
        nentries, NumPoints + 1, [&](std::uint64_t e) { return point_of[e]; }, order, out.mStart);
    out.mStart.pop_back();  // drop the invalid bucket's end
    const std::uint64_t kept = out.mStart[NumPoints];
    out.mCellOf.resize(kept);
    parallel_for_bw(kept, [&](std::size_t j) { out.mCellOf[j] = cell_of[order[j]]; });
    return out;
}

}  // namespace

CellPointWeight cell_point_weight_from_name(const std::string& rName) {
    if (rName == "uniform" || rName == "simple")
        return CellPointWeight::Uniform;
    if (rName == "measure" || rName == "area" || rName == "volume")
        return CellPointWeight::Measure;
    throw std::invalid_argument("meshio++: unknown averaging weight '" + rName +
                                "' (expected 'uniform' or 'measure')");
}

Mesh point_data_to_cell_data(const Mesh& rMesh, const DataAverageOptions& rOpts) {
    const std::vector<std::string> names = davg_resolve_names(rMesh, DataLocation::Point, rOpts);
    for (const std::string& n : names)
        davg_check_target(rMesh, DataLocation::Cell, davg_target_name(n, rOpts), rOpts);

    Mesh out = detail::clone_mesh(rMesh);
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (nblocks == 0)
        return out;

    for (const std::string& name : names) {
        const NDArray& src = rMesh.PointData(name);
        const std::size_t ncomp = data_num_components(src);
        const std::size_t npoints = rMesh.NumPoints();

        // One output array per cell block — the uniform API requires exactly
        // NumCellBlocks() of them.
        std::vector<NDArray> blocks;
        blocks.reserve(nblocks);
        for (std::size_t b = 0; b < nblocks; ++b) {
            const auto cb = rMesh.Cells(b);
            const std::size_t nc = cb.NumCells();
            NDArray dst(DType::Float64, davg_shape(nc, ncomp));
            double* pdst = dst.As<double>();

            // Pure gather: each cell reads only its own nodes, so this is safe
            // to run in parallel and is genuine per-element compute.
            parallel_for(nc, [&](std::size_t c) {
                std::vector<std::int64_t> nodes;
                davg_cell_nodes(cb, c, nodes);
                for (std::size_t k = 0; k < ncomp; ++k) {
                    double sum = 0.0;
                    std::int64_t count = 0;
                    for (std::int64_t node : nodes) {
                        if (node < 0 || static_cast<std::size_t>(node) >= npoints)
                            continue;
                        const double v =
                            detail::read_double(src, static_cast<std::size_t>(node) * ncomp + k);
                        if (std::isfinite(v)) {
                            sum += v;
                            ++count;
                        }
                    }
                    pdst[c * ncomp + k] =
                        count > 0 ? sum / static_cast<double>(count) : std::nan("");
                }
            });

            // The NaN policy is applied serially afterwards so that a `Fail`
            // throw is deterministic rather than racing between threads.
            for (std::size_t i = 0; i < nc * ncomp; ++i)
                pdst[i] =
                    davg_apply_nan_policy(pdst[i], rOpts, "cell", name, i / (ncomp ? ncomp : 1));
            blocks.push_back(std::move(dst));
        }
        out.AddCellData(davg_target_name(name, rOpts), std::move(blocks));
    }
    return out;
}

Mesh cell_data_to_point_data(const Mesh& rMesh, const DataAverageOptions& rOpts) {
    const std::vector<std::string> names = davg_resolve_names(rMesh, DataLocation::Cell, rOpts);
    for (const std::string& n : names)
        davg_check_target(rMesh, DataLocation::Point, davg_target_name(n, rOpts), rOpts);

    Mesh out = detail::clone_mesh(rMesh);
    const std::size_t npoints = rMesh.NumPoints();
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (npoints == 0 || nblocks == 0)
        return out;

    // Per-cell weights, computed once and reused by every array.
    std::vector<std::vector<double>> weights(nblocks);
    bool warned_measure = false;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        weights[b].assign(cb.NumCells(), 1.0);
        if (rOpts.weight != CellPointWeight::Measure)
            continue;
        const NDArray& points = rMesh.Points();
        const std::size_t pdim = rMesh.PointDim();
        std::vector<double>& w = weights[b];
        parallel_for(cb.NumCells(),
                     [&](std::size_t c) { w[c] = detail::cell_measure(points, pdim, cb, c); });
        for (double& v : w) {
            if (!std::isfinite(v) || v <= 0.0) {
                v = 1.0;
                if (!warned_measure) {
                    warned_measure = true;
                    log::warn(
                        "data_average: cell measure unavailable for some cells (ragged or "
                        "unsupported type); falling back to a unit weight there");
                }
            }
        }
    }

    // Built on first use and shared by every array.
    std::unique_ptr<DavgIncidence> incidence;
    std::vector<double> cell_weight;
    std::vector<std::size_t> cell_block_base;
    for (const std::string& name : names) {
        if (rMesh.CellDataNumBlocks(name) != nblocks)
            throw std::invalid_argument("meshio++: data_average: cell_data '" + name + "' has " +
                                        std::to_string(rMesh.CellDataNumBlocks(name)) +
                                        " block(s) but the mesh has " + std::to_string(nblocks) +
                                        " cell block(s)");

        const std::size_t ncomp = data_num_components(rMesh.CellData(name, 0));
        for (std::size_t b = 0; b < nblocks; ++b)
            if (data_num_components(rMesh.CellData(name, b)) != ncomp)
                throw std::invalid_argument("meshio++: data_average: cell_data '" + name +
                                            "' has inconsistent component counts across blocks");
        NDArray dst(DType::Float64, davg_shape(npoints, ncomp));
        double* pdst = dst.As<double>();
        if (detail::slot_runs_impl::workers() <= 1) {
            // One worker (a SEQ build, or one thread): the plain scatter, which
            // adds the same terms in the same order as the gather below and is
            // cheaper than building the incidence table first.
            std::vector<double> sum(npoints * ncomp, 0.0);
            std::vector<double> wsum(npoints * ncomp, 0.0);
            std::vector<std::int64_t> nodes;
            for (std::size_t b = 0; b < nblocks; ++b) {
                const auto cb = rMesh.Cells(b);
                const detail::DoubleView src(rMesh.CellData(name, b));
                for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                    davg_cell_nodes(cb, c, nodes);
                    const double w = weights[b][c];
                    for (std::size_t k = 0; k < ncomp; ++k) {
                        const double v = src[c * ncomp + k];
                        if (!std::isfinite(v))
                            continue;
                        for (std::int64_t node : nodes) {
                            if (node < 0 || static_cast<std::size_t>(node) >= npoints)
                                continue;
                            const std::size_t idx = static_cast<std::size_t>(node) * ncomp + k;
                            sum[idx] += w * v;
                            wsum[idx] += w;
                        }
                    }
                }
            }
            for (std::size_t i = 0; i < npoints * ncomp; ++i)
                pdst[i] = wsum[i] > 0.0 ? sum[i] / wsum[i] : std::nan("");
            for (std::size_t i = 0; i < npoints * ncomp; ++i)
                pdst[i] =
                    davg_apply_nan_policy(pdst[i], rOpts, "point", name, i / (ncomp ? ncomp : 1));
            out.AddPointData(davg_target_name(name, rOpts), std::move(dst));
            continue;
        }
        if (!incidence) {
            std::vector<std::size_t> block_base(nblocks + 1, 0);
            for (std::size_t b = 0; b < nblocks; ++b)
                block_base[b + 1] = block_base[b] + rMesh.Cells(b).NumCells();
            incidence = std::make_unique<DavgIncidence>(davg_incidence(rMesh, block_base, npoints));
            cell_weight.resize(block_base[nblocks]);
            for (std::size_t b = 0; b < nblocks; ++b)
                std::copy(weights[b].begin(), weights[b].end(),
                          cell_weight.begin() + static_cast<std::ptrdiff_t>(block_base[b]));
            cell_block_base = std::move(block_base);
        }
        // The cell values as doubles, block-major (read_double's conversion).
        std::vector<double> values(cell_weight.size() * ncomp);
        for (std::size_t b = 0; b < nblocks; ++b) {
            const detail::DoubleView src(rMesh.CellData(name, b));
            const std::size_t n = (cell_block_base[b + 1] - cell_block_base[b]) * ncomp;
            double* to = values.data() + cell_block_base[b] * ncomp;
            parallel_for_bw(n, [&](std::size_t i) { to[i] = src[i]; });
        }

        // --- accumulation, gathered per point: each point adds its incident
        // --- cells' weighted values in (block, cell, node) order -- the order
        // --- the serial scatter added them in (FP add is not associative) --
        // --- so the sums are bit-identical on every backend and thread count.
        const DavgIncidence& inc = *incidence;
        parallel_for(npoints, [&](std::size_t pnt) {
            for (std::size_t k = 0; k < ncomp; ++k) {
                double sum = 0.0;
                double wsum = 0.0;
                for (std::uint64_t j = inc.mStart[pnt]; j < inc.mStart[pnt + 1]; ++j) {
                    const std::uint64_t gc = inc.mCellOf[j];
                    const double v = values[gc * ncomp + k];
                    if (!std::isfinite(v))
                        continue;
                    const double w = cell_weight[gc];
                    sum += w * v;
                    wsum += w;
                }
                pdst[pnt * ncomp + k] = wsum > 0.0 ? sum / wsum : std::nan("");
            }
        });
        for (std::size_t i = 0; i < npoints * ncomp; ++i)
            pdst[i] = davg_apply_nan_policy(pdst[i], rOpts, "point", name, i / (ncomp ? ncomp : 1));

        out.AddPointData(davg_target_name(name, rOpts), std::move(dst));
    }
    return out;
}

}  // namespace meshioplusplus
