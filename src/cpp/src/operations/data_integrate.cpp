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
// Cell-measure-weighted field integration: see operations/data_integrate.hpp
// for the contract. Built entirely through the uniform mesh API, so it
// compiles under every mesh backend.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/data_integrate.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// The cells belonging to one named Cell region, resolved to (block, row)
// pairs once and reused across every requested array.
struct DintRegionCells {
    std::string mName;
    std::vector<std::pair<std::size_t, std::int64_t>> mCells;
};

// Turns per-component accumulators into the public report shape.
FieldIntegralRegion dint_build_region(std::string name,
                                      const std::vector<detail::WeightedSum>& rStats,
                                      std::int64_t NumCells, std::int64_t NumSkipped,
                                      std::size_t Comps) {
    FieldIntegralRegion out;
    out.mName = std::move(name);
    out.mNumCells = NumCells;
    out.mNumSkipped = NumSkipped;
    out.mDomainMeasurePerComponent.resize(Comps);
    out.mTotalPerComponent.resize(Comps);
    out.mMeanPerComponent.resize(Comps);
    out.mNumNanPerComponent.resize(Comps);
    for (std::size_t k = 0; k < Comps; ++k) {
        out.mDomainMeasurePerComponent[k] = rStats[k].mSumW;
        out.mTotalPerComponent[k] = rStats[k].mSumVW;
        out.mMeanPerComponent[k] = rStats[k].Mean();
        out.mNumNanPerComponent[k] = rStats[k].mNumNan;
    }
    return out;
}

std::string dint_unknown_or_point_message(const Mesh& rMesh, const std::string& rName) {
    if (rMesh.HasPointData(rName))
        return "meshio++: data_integrate: '" + rName +
               "' is a point_data array; convert it first with point_data_to_cell_data (CLI: "
               "`data to-cell`)";
    return data_unknown_key_message(rMesh, DataLocation::Cell, rName);
}

}  // namespace

DataIntegrateReport data_integrate(const Mesh& rMesh, const DataIntegrateOptions& rOptions) {
    // --- resolve which arrays to integrate (validation before any work) -----
    std::vector<std::string> names = rOptions.mArrayNames;
    if (names.empty()) {
        names = rMesh.CellDataNames();  // sorted on every backend
    } else {
        for (const std::string& name : names)
            if (!rMesh.HasCellData(name))
                throw std::invalid_argument(dint_unknown_or_point_message(rMesh, name));
    }
    for (const std::string& name : names)
        if (rMesh.CellDataNumBlocks(name) != rMesh.NumCellBlocks())
            throw std::invalid_argument("meshio++: data_integrate: cell_data '" + name +
                                        "' does not have one array per cell block");

    DataIntegrateReport report;
    if (names.empty())
        return report;

    const std::size_t nblocks = rMesh.NumCellBlocks();

    // --- per-cell measures, computed once and shared by every array/region --
    // A negative sentinel marks a cell whose measure could not be computed
    // (ragged/unsupported/degenerate), so every array's accumulation excludes
    // it identically without recomputing anything geometric.
    std::vector<std::vector<double>> measures(nblocks);
    std::int64_t domain_num_cells = 0;
    std::int64_t domain_num_skipped = 0;
    if (nblocks > 0) {
        const NDArray& points = rMesh.Points();
        const std::size_t pdim = rMesh.PointDim();
        for (std::size_t b = 0; b < nblocks; ++b) {
            const auto cb = rMesh.Cells(b);
            std::vector<double>& w = measures[b];
            w.assign(cb.NumCells(), 0.0);
            parallel_for(cb.NumCells(), [&](std::size_t c) {
                w[c] = std::abs(detail::cell_measure(points, pdim, cb, c));
            });
            for (double& v : w) {
                if (std::isfinite(v) && v > 0.0)
                    ++domain_num_cells;
                else {
                    v = -1.0;  // sentinel: unmeasurable
                    ++domain_num_skipped;
                }
            }
        }
    }

    // --- named Cell regions, resolved once -----------------------------------
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    std::vector<DintRegionCells> region_defs;
    for (const std::string& rname : rMesh.RegionNames()) {
        const std::size_t idx = rMesh.FindRegion(rname, RegionKind::Cell);
        if (idx == Mesh::npos)
            continue;
        const Region& region = rMesh.Region(idx);
        DintRegionCells rc;
        rc.mName = rname;
        const std::int64_t* entries = region.Entries();
        const std::size_t n = region.NumEntries();
        rc.mCells.reserve(n);
        for (std::size_t e = 0; e < n; ++e) {
            const auto [b, row] = detail::global_to_block_row(bases, entries[e]);
            if (b != static_cast<std::size_t>(-1))
                rc.mCells.push_back({b, row});
        }
        region_defs.push_back(std::move(rc));
    }
    std::vector<std::int64_t> region_num_cells(region_defs.size(), 0);
    std::vector<std::int64_t> region_num_skipped(region_defs.size(), 0);
    for (std::size_t r = 0; r < region_defs.size(); ++r)
        for (const auto& [b, row] : region_defs[r].mCells) {
            if (measures[b][static_cast<std::size_t>(row)] > 0.0)
                ++region_num_cells[r];
            else
                ++region_num_skipped[r];
        }

    // --- per array: whole-mesh, then every region ----------------------------
    report.mArrays.reserve(names.size());
    for (const std::string& name : names) {
        FieldIntegralArray arr;
        arr.mName = name;
        const std::size_t comps = data_num_components(rMesh.CellData(name, 0));
        arr.mNumComponents = static_cast<std::int64_t>(comps);

        std::vector<detail::WeightedSum> domain_stats(comps);
        for (std::size_t b = 0; b < nblocks; ++b)
            detail::accumulate_weighted(rMesh.CellData(name, b), comps, measures[b], domain_stats);
        arr.mDomain =
            dint_build_region("", domain_stats, domain_num_cells, domain_num_skipped, comps);

        arr.mRegions.reserve(region_defs.size());
        for (std::size_t r = 0; r < region_defs.size(); ++r) {
            // Regions are typically few and small relative to the whole mesh,
            // so a plain serial fold (in the region's own canonical, sorted
            // entry order) is simpler than chunked parallelism here and is
            // trivially deterministic.
            std::vector<detail::WeightedSum> region_stats(comps);
            for (const auto& [b, row] : region_defs[r].mCells) {
                const double w = measures[b][static_cast<std::size_t>(row)];
                if (!(w > 0.0))
                    continue;
                const NDArray& block_arr = rMesh.CellData(name, b);
                for (std::size_t k = 0; k < comps; ++k)
                    region_stats[k].Add(
                        detail::read_double(block_arr, static_cast<std::size_t>(row) * comps + k),
                        w);
            }
            arr.mRegions.push_back(dint_build_region(region_defs[r].mName, region_stats,
                                                     region_num_cells[r], region_num_skipped[r],
                                                     comps));
        }

        report.mArrays.push_back(std::move(arr));
    }

    return report;
}

}  // namespace meshioplusplus
