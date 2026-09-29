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

// The VTK XML region convention. See region_field_data.hpp.

// System includes
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "region_field_data.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr const char* kRfdPrefix = "region:";
constexpr const char* kRfdMetaPrefix = "region-meta:";

bool rfd_starts_with(const std::string& rS, const char* pPrefix) {
    return rS.rfind(pPrefix, 0) == 0;
}

const char* rfd_kind_name(RegionKind Kind) {
    switch (Kind) {
        case RegionKind::Point:
            return "point";
        case RegionKind::Cell:
            return "cell";
        default:
            return "side";
    }
}

/// `<kind>:<name>` after a prefix -> kind and name; false if malformed.
bool rfd_parse(const std::string& rRest, RegionKind& rKind, std::string& rName) {
    const std::size_t colon = rRest.find(':');
    if (colon == std::string::npos)
        return false;
    const std::string kind = rRest.substr(0, colon);
    if (kind == "point")
        rKind = RegionKind::Point;
    else if (kind == "cell")
        rKind = RegionKind::Cell;
    else if (kind == "side")
        rKind = RegionKind::Side;
    else
        return false;
    rName = rRest.substr(colon + 1);
    return true;
}

bool rfd_is_int(DType Dtype) {
    return Dtype != DType::Float32 && Dtype != DType::Float64;
}

}  // namespace

bool is_region_field_name(const std::string& rName) {
    return rfd_starts_with(rName, kRfdPrefix) || rfd_starts_with(rName, kRfdMetaPrefix);
}

std::vector<std::pair<std::string, NDArray>> regions_to_field_arrays(
    const Mesh& rMesh, const std::vector<std::int64_t>* pGlobalToFile) {
    std::vector<std::pair<std::string, NDArray>> out;
    auto file_cell = [&](std::int64_t Global) -> std::int64_t {
        if (!pGlobalToFile)
            return Global;
        return Global >= 0 && static_cast<std::size_t>(Global) < pGlobalToFile->size()
                   ? (*pGlobalToFile)[static_cast<std::size_t>(Global)]
                   : -1;
    };
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r_region = rMesh.Region(i);
        const std::size_t n = r_region.NumEntries();
        const std::size_t stride = r_region.Stride();
        const std::int64_t* e = r_region.Entries();
        NDArray a = stride == 1 ? NDArray(DType::Int64, {n}) : NDArray(DType::Int64, {n, 2});
        std::int64_t* p = a.As<std::int64_t>();
        for (std::size_t k = 0; k < n; ++k) {
            if (r_region.mKind == RegionKind::Point) {
                p[k] = e[k];
            } else if (r_region.mKind == RegionKind::Cell) {
                p[k] = file_cell(e[k]);
            } else {
                p[2 * k] = file_cell(e[2 * k]);
                p[2 * k + 1] = e[2 * k + 1];
            }
        }
        const std::string key = std::string(rfd_kind_name(r_region.mKind)) + ":" + r_region.mName;
        out.emplace_back(kRfdPrefix + key, std::move(a));
        if (r_region.mDim != -1 || r_region.mTag != -1) {
            NDArray meta(DType::Int64, {2});
            meta.As<std::int64_t>()[0] = r_region.mDim;
            meta.As<std::int64_t>()[1] = r_region.mTag;
            out.emplace_back(kRfdMetaPrefix + key, std::move(meta));
        }
    }
    return out;
}

void regions_from_field_arrays(Mesh& rMesh, std::vector<std::pair<std::string, NDArray>>& rArrays,
                               const std::vector<std::int64_t>* pFileToGlobal,
                               const char* pFormat) {
    // [dim, tag] by "<kind>:<name>", read first so a region array can use it.
    std::map<std::string, std::pair<int, std::int64_t>> metas;
    for (const auto& [name, arr] : rArrays) {
        if (!rfd_starts_with(name, kRfdMetaPrefix))
            continue;
        if (rfd_is_int(arr.Dtype()) && arr.Size() == 2)
            metas[name.substr(std::string(kRfdMetaPrefix).size())] = {
                static_cast<int>(read_int(arr, 0)), read_int(arr, 1)};
    }
    const auto n_points = static_cast<std::int64_t>(rMesh.NumPoints());
    std::int64_t n_cells = 0;
    for (const auto cb : rMesh.CellRange())
        n_cells += static_cast<std::int64_t>(cb.NumCells());
    auto global_cell = [&](std::int64_t File) -> std::int64_t {
        if (!pFileToGlobal)
            return File >= 0 && File < n_cells ? File : -1;
        return File >= 0 && static_cast<std::size_t>(File) < pFileToGlobal->size()
                   ? (*pFileToGlobal)[static_cast<std::size_t>(File)]
                   : -1;
    };

    for (auto& [name, arr] : rArrays) {
        if (rfd_starts_with(name, kRfdMetaPrefix)) {
            if (metas.count(name.substr(std::string(kRfdMetaPrefix).size())))
                continue;  // consumed
            log::warn("{}: malformed region meta array '{}' kept as field data", pFormat, name);
            rMesh.AddFieldData(name, std::move(arr));
            continue;
        }
        RegionKind kind{};
        std::string region_name;
        if (!rfd_starts_with(name, kRfdPrefix) ||
            !rfd_parse(name.substr(std::string(kRfdPrefix).size()), kind, region_name)) {
            rMesh.AddFieldData(name, std::move(arr));
            continue;
        }
        const std::size_t stride = kind == RegionKind::Side ? 2 : 1;
        const std::vector<std::size_t>& shape = arr.Shape();
        const std::size_t comps = shape.size() >= 2 ? shape[1] : 1;
        bool ok = rfd_is_int(arr.Dtype()) && comps == stride && arr.Size() % stride == 0;
        const std::size_t n = ok ? arr.Size() / stride : 0;
        std::vector<std::int64_t> flat(arr.Size());
        for (std::size_t k = 0; ok && k < n; ++k) {
            if (kind == RegionKind::Point) {
                flat[k] = read_int(arr, k);
                ok = flat[k] >= 0 && flat[k] < n_points;
            } else {
                const std::int64_t g = global_cell(read_int(arr, k * stride));
                ok = g >= 0;
                flat[k * stride] = g;
                if (stride == 2) {
                    flat[k * 2 + 1] = read_int(arr, k * 2 + 1);
                    ok = ok && flat[k * 2 + 1] >= 0;
                }
            }
        }
        if (!ok) {
            log::warn("{}: '{}' is not a well-formed region array; kept as field data", pFormat,
                      name);
            rMesh.AddFieldData(name, std::move(arr));
            continue;
        }
        NDArray entries = stride == 1 ? NDArray(DType::Int64, {n}) : NDArray(DType::Int64, {n, 2});
        std::copy(flat.begin(), flat.end(), entries.As<std::int64_t>());
        int dim = -1;
        std::int64_t tag = -1;
        const auto it = metas.find(name.substr(std::string(kRfdPrefix).size()));
        if (it != metas.end()) {
            dim = it->second.first;
            tag = it->second.second;
        }
        rMesh.AddRegion(Region(region_name, kind, dim, tag, std::move(entries)));
    }
    rArrays.clear();
}

std::vector<RegionSummary> region_summaries_from_field_names(
    const std::vector<std::pair<std::string, std::size_t>>& rArrays) {
    std::vector<RegionSummary> out;
    for (const auto& [name, tuples] : rArrays) {
        RegionSummary rs;
        if (!rfd_starts_with(name, kRfdPrefix) ||
            !rfd_parse(name.substr(std::string(kRfdPrefix).size()), rs.mKind, rs.mName))
            continue;
        rs.mNumEntries = tuples;
        out.push_back(std::move(rs));
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
