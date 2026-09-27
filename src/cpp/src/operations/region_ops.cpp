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

// edit_regions. See operations/region_ops.hpp for the contract.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/region_ops.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kRoPrefix = "meshio++: edit_regions: ";

std::string ro_describe(const meshioplusplus::Region& rR) {
    return "'" + rR.mName + "' (" + region_kind_name(rR.mKind) + ", dim " +
           std::to_string(rR.mDim) + ", tag " + std::to_string(rR.mTag) + ")";
}

std::string ro_describe(const RegionSelector& rS) {
    std::string s = "'" + rS.mName + "'";
    if (rS.mKind >= 0 && rS.mKind <= 2)
        s += std::string(" kind ") + region_kind_name(static_cast<RegionKind>(rS.mKind));
    if (rS.mDim != kRegionAny)
        s += " dim " + std::to_string(rS.mDim);
    if (rS.mTag != kRegionAny)
        s += " tag " + std::to_string(rS.mTag);
    return s;
}

bool ro_matches(const meshioplusplus::Region& rR, const RegionSelector& rS) {
    return rR.mName == rS.mName &&
           (rS.mKind < 0 || static_cast<std::int32_t>(rR.mKind) == rS.mKind) &&
           (rS.mDim == kRegionAny || rR.mDim == rS.mDim) &&
           (rS.mTag == kRegionAny || rR.mTag == rS.mTag);
}

// A region's entries as sortable rows: (entry, 0) for point/cell, (cell, facet)
// for side. Stored entries are canonical, so the rows come out sorted.
std::vector<std::pair<std::int64_t, std::int64_t>> ro_rows(const meshioplusplus::Region& rR) {
    const std::size_t n = rR.NumEntries();
    const std::size_t stride = rR.Stride();
    std::vector<std::pair<std::int64_t, std::int64_t>> out(n);
    for (std::size_t i = 0; i < n; ++i)
        out[i] = {detail::read_int(rR.mEntries, i * stride),
                  stride == 2 ? detail::read_int(rR.mEntries, i * stride + 1) : 0};
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

NDArray ro_entries(const std::vector<std::pair<std::int64_t, std::int64_t>>& rRows,
                   RegionKind Kind) {
    const bool side = Kind == RegionKind::Side;
    NDArray a = side ? NDArray::Uninit(DType::Int64, {rRows.size(), std::size_t{2}})
                     : NDArray::Uninit(DType::Int64, {rRows.size()});
    std::int64_t* d = a.As<std::int64_t>();
    for (std::size_t i = 0; i < rRows.size(); ++i) {
        if (side) {
            d[i * 2] = rRows[i].first;
            d[i * 2 + 1] = rRows[i].second;
        } else {
            d[i] = rRows[i].first;
        }
    }
    return a;
}

// The key an edit's output will have must not belong to an unrelated region.
void ro_check_free(const Mesh& rMesh, const meshioplusplus::Region& rOut,
                   const std::vector<std::size_t>& rInputs, const char* pOp) {
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        if (rMesh.Region(i).Key() != rOut.Key())
            continue;
        if (std::find(rInputs.begin(), rInputs.end(), i) != rInputs.end())
            return;  // replacing one of the edit's own inputs
        throw std::invalid_argument(std::string(kRoPrefix) + pOp + ": the result " +
                                    ro_describe(rOut) + " would replace an existing region");
    }
}

// Remove regions @p rIdx (indices into the current list), highest first so the
// lower indices stay valid.
void ro_remove(Mesh& rMesh, std::vector<std::size_t> rIdx) {
    std::sort(rIdx.begin(), rIdx.end());
    rIdx.erase(std::unique(rIdx.begin(), rIdx.end()), rIdx.end());
    for (auto it = rIdx.rbegin(); it != rIdx.rend(); ++it)
        rMesh.RemoveRegion(*it);
}

void ro_apply(Mesh& rMesh, const RegionEdit& rEdit) {
    const char* op = region_op_name(rEdit.mOp);
    std::vector<std::size_t> in;
    in.reserve(rEdit.mInputs.size());
    for (const RegionSelector& s : rEdit.mInputs)
        in.push_back(find_region(rMesh, s));
    const std::size_t nin = in.size();

    switch (rEdit.mOp) {
        case RegionOp::Delete: {
            if (nin == 0)
                throw std::invalid_argument(std::string(kRoPrefix) + "delete: no region given");
            ro_remove(rMesh, in);
            return;
        }
        case RegionOp::Rename:
        case RegionOp::Retag: {
            if (nin != 1)
                throw std::invalid_argument(std::string(kRoPrefix) + op +
                                            ": takes exactly one region, got " +
                                            std::to_string(nin));
            const bool rename = rEdit.mOp == RegionOp::Rename;
            if (rename && rEdit.mOutputName.empty())
                throw std::invalid_argument(std::string(kRoPrefix) + "rename: no new name given");
            if (!rename && rEdit.mOutputTag == kRegionAny && rEdit.mOutputDim == kRegionAny)
                throw std::invalid_argument(std::string(kRoPrefix) +
                                            "retag: give a new tag and/or dimension");
            meshioplusplus::Region out = rMesh.Region(in[0]);
            if (!rEdit.mOutputName.empty())
                out.mName = rEdit.mOutputName;
            if (rEdit.mOutputDim != kRegionAny)
                out.mDim = static_cast<int>(rEdit.mOutputDim);
            if (rEdit.mOutputTag != kRegionAny)
                out.mTag = rEdit.mOutputTag;
            ro_check_free(rMesh, out, in, op);
            if (!rename && out.mKind == RegionKind::Cell && rMesh.HasCellData("gmsh:physical"))
                log::warn(
                    "edit_regions: retagged cell region '{}', but the mesh carries "
                    "'gmsh:physical' cell data, which a gmsh write takes the physical tag from; "
                    "that array still holds the old tag",
                    out.mName);
            ro_remove(rMesh, in);
            rMesh.AddRegion(std::move(out));
            return;
        }
        case RegionOp::Union:
        case RegionOp::Intersection:
        case RegionOp::Difference: {
            if (nin < 2)
                throw std::invalid_argument(std::string(kRoPrefix) + op +
                                            ": takes two or more regions, got " +
                                            std::to_string(nin));
            if (rEdit.mOutputName.empty())
                throw std::invalid_argument(std::string(kRoPrefix) + op +
                                            ": no name given for the result");
            const meshioplusplus::Region& first = rMesh.Region(in[0]);
            int dim = first.mDim;
            for (std::size_t k = 1; k < nin; ++k) {
                const meshioplusplus::Region& r = rMesh.Region(in[k]);
                if (r.mKind != first.mKind)
                    throw std::invalid_argument(std::string(kRoPrefix) + op + ": " +
                                                ro_describe(first) + " and " + ro_describe(r) +
                                                " are of different kinds");
                if (r.mDim != dim)
                    dim = -1;
            }
            auto acc = ro_rows(first);
            for (std::size_t k = 1; k < nin; ++k) {
                const auto rows = ro_rows(rMesh.Region(in[k]));
                std::vector<std::pair<std::int64_t, std::int64_t>> next;
                if (rEdit.mOp == RegionOp::Union)
                    std::set_union(acc.begin(), acc.end(), rows.begin(), rows.end(),
                                   std::back_inserter(next));
                else if (rEdit.mOp == RegionOp::Intersection)
                    std::set_intersection(acc.begin(), acc.end(), rows.begin(), rows.end(),
                                          std::back_inserter(next));
                else
                    std::set_difference(acc.begin(), acc.end(), rows.begin(), rows.end(),
                                        std::back_inserter(next));
                acc = std::move(next);
            }
            meshioplusplus::Region out;
            out.mName = rEdit.mOutputName;
            out.mKind = first.mKind;
            out.mDim = rEdit.mOutputDim != kRegionAny ? static_cast<int>(rEdit.mOutputDim) : dim;
            out.mTag = rEdit.mOutputTag != kRegionAny ? rEdit.mOutputTag : -1;
            out.mEntries = ro_entries(acc, out.mKind);
            ro_check_free(rMesh, out, in, op);
            if (!rEdit.mKeepInputs)
                ro_remove(rMesh, in);
            rMesh.AddRegion(std::move(out));
            return;
        }
    }
    throw std::invalid_argument(std::string(kRoPrefix) + "unknown operation");
}

}  // namespace

const char* region_op_name(RegionOp Op) {
    switch (Op) {
        case RegionOp::Union:
            return "union";
        case RegionOp::Intersection:
            return "intersection";
        case RegionOp::Difference:
            return "difference";
        case RegionOp::Rename:
            return "rename";
        case RegionOp::Retag:
            return "retag";
        case RegionOp::Delete:
            return "delete";
    }
    return "unknown";
}

RegionOp region_op_from_name(const std::string& rName) {
    if (rName == "union")
        return RegionOp::Union;
    if (rName == "intersection" || rName == "intersect")
        return RegionOp::Intersection;
    if (rName == "difference")
        return RegionOp::Difference;
    if (rName == "rename")
        return RegionOp::Rename;
    if (rName == "retag")
        return RegionOp::Retag;
    if (rName == "delete")
        return RegionOp::Delete;
    throw std::invalid_argument(std::string(kRoPrefix) + "unknown operation '" + rName +
                                "' (expected union, intersection, difference, rename, retag "
                                "or delete)");
}

std::size_t find_region(const Mesh& rMesh, const RegionSelector& rSelector) {
    std::vector<std::size_t> hits;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        if (ro_matches(rMesh.Region(i), rSelector))
            hits.push_back(i);
    if (hits.size() == 1)
        return hits[0];
    std::string list;
    if (hits.empty()) {
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
            list += (list.empty() ? "" : ", ") + ro_describe(rMesh.Region(i));
        throw std::invalid_argument(std::string(kRoPrefix) + "no region matches " +
                                    ro_describe(rSelector) +
                                    " (regions: " + (list.empty() ? "none" : list) + ")");
    }
    for (std::size_t i : hits)
        list += (list.empty() ? "" : ", ") + ro_describe(rMesh.Region(i));
    throw std::invalid_argument(std::string(kRoPrefix) + ro_describe(rSelector) + " matches " +
                                std::to_string(hits.size()) + " regions (" + list +
                                "); pin the kind, dim or tag");
}

Mesh edit_regions(const Mesh& rMesh, const std::vector<RegionEdit>& rEdits) {
    Mesh out = detail::clone_mesh(rMesh);
    for (const RegionEdit& e : rEdits)
        ro_apply(out, e);
    return out;
}

}  // namespace meshioplusplus
