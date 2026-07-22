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
// Read-only per-array summary of a mesh's data maps -- the data-side
// complement to the topological `info` verb and the geometric compute_stats.
// Reductions go through the shared detail::FiniteStats accumulator (parallel
// chunks, serial combine, non-finite values excluded and counted), so the
// numbers match what data_condition would compute from the same array.
//
// This operation never throws on non-finite values: counting them is the point.
// See operations/data_info.hpp for the contract.

// System includes
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/data_info.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {

namespace {

/// Fills the reduction-derived fields of `rInfo` from per-component stats.
void dinfo_finish(DataArrayInfo& rInfo, const std::vector<detail::FiniteStats>& rStats) {
    rInfo.mMinPerComponent.reserve(rStats.size());
    rInfo.mMaxPerComponent.reserve(rStats.size());
    rInfo.mMeanPerComponent.reserve(rStats.size());
    for (const detail::FiniteStats& s : rStats) {
        rInfo.mMinPerComponent.push_back(s.Min());
        rInfo.mMaxPerComponent.push_back(s.Max());
        rInfo.mMeanPerComponent.push_back(s.Mean());
    }
    const detail::FiniteStats all = detail::combine_components(rStats);
    rInfo.mMin = all.Min();
    rInfo.mMax = all.Max();
    rInfo.mMean = all.Mean();
    rInfo.mNumNan = all.mNumNan;
    rInfo.mNumInf = all.mNumInf;
    rInfo.mNumFinite = all.mNumFinite;
}

/// Summarizes a single (non-cell) array.
DataArrayInfo dinfo_simple(DataLocation location, const std::string& rName, const NDArray& rArray) {
    DataArrayInfo info;
    info.mLocation = location;
    info.mName = rName;
    info.mDtype = rArray.Dtype();
    info.mShape = rArray.Shape();
    info.mNumBlocks = 1;
    const std::size_t ncomp = data_num_components(rArray);
    info.mNumComponents = static_cast<std::int64_t>(ncomp);
    info.mNumValues = static_cast<std::int64_t>(rArray.Size());
    info.mNumEntries = ncomp > 0 ? static_cast<std::int64_t>(rArray.Size() / ncomp) : 0;

    std::vector<detail::FiniteStats> stats;
    detail::accumulate_stats(rArray, ncomp, stats);
    if (stats.empty())
        stats.resize(ncomp);
    dinfo_finish(info, stats);
    return info;
}

/// Summarizes one cell_data array across all of its blocks.
DataArrayInfo dinfo_cell(const Mesh& rMesh, const std::string& rName) {
    DataArrayInfo info;
    info.mLocation = DataLocation::Cell;
    info.mName = rName;
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
    info.mNumBlocks = static_cast<std::int64_t>(nblocks);
    if (nblocks == 0) {
        info.mNumComponents = 1;
        dinfo_finish(info, std::vector<detail::FiniteStats>(1));
        return info;
    }

    const NDArray& first = rMesh.CellData(rName, 0);
    info.mDtype = first.Dtype();
    info.mShape = first.Shape();
    const std::size_t ncomp = data_num_components(first);
    info.mNumComponents = static_cast<std::int64_t>(ncomp);

    std::vector<detail::FiniteStats> stats;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const NDArray& block = rMesh.CellData(rName, b);
        if (data_num_components(block) != ncomp) {
            // Diagnostic, not an error: report it and skip the odd block so the
            // rest of the summary is still usable.
            info.mInconsistentBlocks = true;
            continue;
        }
        info.mNumValues += static_cast<std::int64_t>(block.Size());
        info.mNumEntries += static_cast<std::int64_t>(block.Size() / (ncomp > 0 ? ncomp : 1));
        detail::accumulate_stats(block, ncomp, stats);
    }
    if (stats.empty())
        stats.resize(ncomp);
    dinfo_finish(info, stats);
    return info;
}

}  // namespace

DataInfoReport data_info(const Mesh& rMesh) {
    DataInfoReport report;
    for (const std::string& name : rMesh.PointDataNames()) {
        report.mArrays.push_back(dinfo_simple(DataLocation::Point, name, rMesh.PointData(name)));
        ++report.mNumPointData;
    }
    for (const std::string& name : rMesh.CellDataNames()) {
        report.mArrays.push_back(dinfo_cell(rMesh, name));
        ++report.mNumCellData;
    }
    for (const std::string& name : rMesh.FieldDataNames()) {
        report.mArrays.push_back(dinfo_simple(DataLocation::Field, name, rMesh.FieldData(name)));
        ++report.mNumFieldData;
    }
    return report;
}

DataArrayInfo data_array_info(const Mesh& rMesh, DataLocation Location, const std::string& rName) {
    if (!data_has(rMesh, Location, rName))
        throw std::invalid_argument(data_unknown_key_message(rMesh, Location, rName));
    switch (Location) {
        case DataLocation::Point:
            return dinfo_simple(Location, rName, rMesh.PointData(rName));
        case DataLocation::Cell:
            return dinfo_cell(rMesh, rName);
        case DataLocation::Field:
            return dinfo_simple(Location, rName, rMesh.FieldData(rName));
    }
    return DataArrayInfo{};  // unreachable
}

}  // namespace meshioplusplus
