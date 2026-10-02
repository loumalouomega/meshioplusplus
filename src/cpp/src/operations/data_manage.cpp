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
// Array management for a mesh's data arrays: keep-only, drop, rename, applied
// in that order. Validation runs against the input mesh up front so a bad key
// costs no work; the rewrite itself is a single detail::clone_mesh pass whose
// filter consults the precomputed per-location decision tables. Geometry is
// copied verbatim. See operations/data_manage.hpp for the contract.

// System includes
#include <array>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <map>
#include <cstring>
#include <type_traits>

// Project includes
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"

namespace meshioplusplus {

namespace {

/// The three locations, in the order they are reported.
constexpr std::array<DataLocation, 3> DMANAGE_LOCATIONS = {DataLocation::Point, DataLocation::Cell,
                                                           DataLocation::Field};

/// Index of a location within `DMANAGE_LOCATIONS`.
std::size_t dmanage_index(DataLocation loc) {
    switch (loc) {
        case DataLocation::Point:
            return 0;
        case DataLocation::Cell:
            return 1;
        case DataLocation::Field:
            return 2;
    }
    return 0;  // unreachable
}

/// "point_data:T", the form used in the dropped/renamed report lists.
std::string dmanage_qualified(DataLocation loc, const std::string& rName) {
    return std::string(data_location_name(loc)) + ":" + rName;
}

/// Throws unless `rName` exists at `loc` (or the caller asked to ignore that).
/// @return whether the key exists.
bool dmanage_require(const Mesh& rMesh, DataLocation loc, const std::string& rName,
                     bool ignore_missing) {
    if (data_has(rMesh, loc, rName))
        return true;
    if (ignore_missing)
        return false;
    throw std::invalid_argument(data_unknown_key_message(rMesh, loc, rName));
}

/// Per-location decision tables built once from the options.
struct DmanagePlan {
    /// Whether a whitelist applies at this location at all.
    std::array<bool, 3> mHasKeep = {false, false, false};
    /// The whitelist, when `mHasKeep` is set.
    std::array<std::unordered_set<std::string>, 3> mKeep;
    /// Names to remove.
    std::array<std::unordered_set<std::string>, 3> mDrop;
    /// Renames, old name -> new name.
    std::array<std::unordered_map<std::string, std::string>, 3> mRename;
};

}  // namespace

DataManageResult data_manage(const Mesh& rMesh, const DataManageOptions& rOpts) {
    DmanagePlan plan;

    // --- validate + build the keep whitelists -------------------------------
    for (const DataKey& k : rOpts.keep) {
        const std::size_t li = dmanage_index(k.location);
        // A location named at all switches on whitelisting there, even if the
        // key itself turns out to be missing and ignored — "keep nothing from
        // point_data" must still drop everything in point_data.
        plan.mHasKeep[li] = true;
        if (dmanage_require(rMesh, k.location, k.name, rOpts.ignore_missing))
            plan.mKeep[li].insert(k.name);
    }

    // --- validate + build the drop sets --------------------------------------
    for (const DataKey& k : rOpts.drop) {
        if (dmanage_require(rMesh, k.location, k.name, rOpts.ignore_missing))
            plan.mDrop[dmanage_index(k.location)].insert(k.name);
    }

    // --- validate + build the rename maps ------------------------------------
    // Sources and targets are both checked for collisions. A rename whose
    // source survives keep/drop but whose target collides with another array
    // that is *not* itself being renamed away would silently clobber it.
    std::array<std::unordered_set<std::string>, 3> rename_targets;
    for (const DataRename& r : rOpts.rename) {
        const std::size_t li = dmanage_index(r.location);
        const char* loc = data_location_name(r.location);
        if (!dmanage_require(rMesh, r.location, r.from, rOpts.ignore_missing))
            continue;
        if (r.to.empty())
            throw std::invalid_argument(std::string("meshio++: data_manage: cannot rename ") + loc +
                                        " '" + r.from + "' to an empty name");
        if (plan.mRename[li].count(r.from) != 0)
            throw std::invalid_argument(std::string("meshio++: data_manage: ") + loc + " '" +
                                        r.from + "' is renamed twice");
        if (!rename_targets[li].insert(r.to).second)
            throw std::invalid_argument(std::string("meshio++: data_manage: two renames both "
                                                    "target ") +
                                        loc + " '" + r.to + "'");
        plan.mRename[li].emplace(r.from, r.to);
    }
    // Target-collision check, once every rename source is known.
    for (const DataRename& r : rOpts.rename) {
        const std::size_t li = dmanage_index(r.location);
        if (plan.mRename[li].count(r.from) == 0)
            continue;  // skipped as missing
        if (r.to == r.from)
            continue;
        const bool target_exists = data_has(rMesh, r.location, r.to);
        const bool target_renamed_away = plan.mRename[li].count(r.to) != 0;
        const bool target_dropped = plan.mDrop[li].count(r.to) != 0 ||
                                    (plan.mHasKeep[li] && plan.mKeep[li].count(r.to) == 0);
        if (target_exists && !target_renamed_away && !target_dropped)
            throw std::invalid_argument(std::string("meshio++: data_manage: cannot rename ") +
                                        data_location_name(r.location) + " '" + r.from + "' to '" +
                                        r.to + "' (that name already exists)");
    }

    // --- rewrite -------------------------------------------------------------
    DataManageResult result;
    std::vector<std::string> dropped;
    std::vector<std::pair<std::string, std::string>> renamed;
    result.mMesh = detail::clone_mesh(rMesh, [&](DataLocation loc, const std::string& rName,
                                                 std::string& rTarget) {
        const std::size_t li = dmanage_index(loc);
        if (plan.mHasKeep[li] && plan.mKeep[li].count(rName) == 0) {
            dropped.push_back(dmanage_qualified(loc, rName));
            return false;
        }
        if (plan.mDrop[li].count(rName) != 0) {
            dropped.push_back(dmanage_qualified(loc, rName));
            return false;
        }
        const auto it = plan.mRename[li].find(rName);
        if (it != plan.mRename[li].end() && it->second != rName) {
            rTarget = it->second;
            renamed.emplace_back(dmanage_qualified(loc, rName), dmanage_qualified(loc, it->second));
        }
        return true;
    });

    std::sort(dropped.begin(), dropped.end());
    result.mDropped = std::move(dropped);
    result.mRenamed = std::move(renamed);
    return result;
}

Mesh data_drop(const Mesh& rMesh, DataLocation Location, const std::vector<std::string>& rNames,
               bool IgnoreMissing) {
    DataManageOptions opts;
    opts.ignore_missing = IgnoreMissing;
    for (const std::string& n : rNames)
        opts.drop.push_back(DataKey{Location, n});
    return data_manage(rMesh, opts).mMesh;
}

Mesh data_keep(const Mesh& rMesh, DataLocation Location, const std::vector<std::string>& rNames,
               bool IgnoreMissing) {
    // Validate the requested names against the input first, so an unknown key
    // is reported the same way `drop` reports one.
    std::unordered_set<std::string> wanted;
    for (const std::string& n : rNames)
        if (dmanage_require(rMesh, Location, n, IgnoreMissing))
            wanted.insert(n);

    // Expressed as the complement — drop everything at this location that was
    // not asked for. This keeps an empty `rNames` meaningful (drop the lot)
    // without needing a "whitelist is active" sentinel key.
    DataManageOptions opts;
    for (const std::string& n : data_names(rMesh, Location))
        if (wanted.count(n) == 0)
            opts.drop.push_back(DataKey{Location, n});
    return data_manage(rMesh, opts).mMesh;
}

Mesh data_rename(const Mesh& rMesh, DataLocation Location, const std::string& rFrom,
                 const std::string& rTo) {
    DataManageOptions opts;
    opts.rename.push_back(DataRename{Location, rFrom, rTo});
    return data_manage(rMesh, opts).mMesh;
}

namespace {

RegionKind dmanage_set_kind(DataLocation location) {
    if (location == DataLocation::Point)
        return RegionKind::Point;
    if (location == DataLocation::Cell)
        return RegionKind::Cell;
    throw std::invalid_argument("meshio++: sets/data conversions require point or cell location");
}

NDArray dmanage_indices(const std::vector<std::int64_t>& rEntries) {
    NDArray data(DType::Int64, {rEntries.size()});
    if (!rEntries.empty())
        std::memcpy(data.Data(), rEntries.data(), rEntries.size() * sizeof(std::int64_t));
    return data;
}

// Exact integer ordering, including UInt64 values above INT64_MAX, without
// converting labels through double or overflowing on INT64_MIN.
struct DmanageTag {
    bool mNegative;
    std::uint64_t mMagnitude;
    bool operator<(const DmanageTag& rOther) const {
        if (mNegative != rOther.mNegative)
            return mNegative;
        return mNegative ? mMagnitude > rOther.mMagnitude : mMagnitude < rOther.mMagnitude;
    }
    std::string Name() const { return (mNegative ? "-" : "") + std::to_string(mMagnitude); }
};

template <class T>
DmanageTag dmanage_tag(T value) {
    if constexpr (std::is_signed_v<T>) {
        if (value < 0)
            return {true, static_cast<std::uint64_t>(-(value + 1)) + 1};
    }
    return {false, static_cast<std::uint64_t>(value)};
}

}  // namespace

Mesh sets_to_data(const Mesh& rMesh, DataLocation Location, const std::optional<std::string>& rName,
                  const std::string& rJoin, const std::vector<std::string>& rOrder) {
    const auto kind = dmanage_set_kind(Location);
    std::vector<std::string> names;
    std::unordered_map<std::string, std::size_t> regions;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const auto& region = rMesh.Region(i);
        if (region.mKind != kind)
            continue;
        if (regions.emplace(region.mName, i).second)
            names.push_back(region.mName);
        else
            regions[region.mName] = i;  // Compatibility view uses the last same-name region.
    }
    if (!rOrder.empty()) {
        std::unordered_set<std::string> seen;
        if (rOrder.size() != names.size())
            throw std::invalid_argument("meshio++: set order must name every set exactly once");
        for (const auto& name : rOrder)
            if (!regions.count(name) || !seen.insert(name).second)
                throw std::invalid_argument("meshio++: invalid/duplicate set name in order: " +
                                            name);
        names = rOrder;
    }
    if (names.empty())
        return detail::clone_mesh(rMesh);
    std::string key;
    if (rName)
        key = *rName;
    else
        for (const auto& name : names) {
            if (&name != &names.front())
                key += rJoin;
            key += name;
        }
    const auto bases = detail::block_bases(rMesh);
    const std::size_t n = Location == DataLocation::Point ? rMesh.NumPoints()
                                                          : static_cast<std::size_t>(bases.back());
    std::vector<std::int64_t> labels(n, -1);
    for (std::size_t label = 0; label < names.size(); ++label) {
        const auto& region = rMesh.Region(regions.at(names[label]));
        for (std::size_t i = 0; i < region.NumEntries(); ++i) {
            auto index = region.Entries()[i];
            if (Location == DataLocation::Point && index < 0 &&
                index >= -static_cast<std::int64_t>(n))
                index += static_cast<std::int64_t>(n);
            if (index < 0 || static_cast<std::uint64_t>(index) >= n) {
                if (Location == DataLocation::Cell)
                    continue;  // Python's global-to-block compatibility view drops these.
                throw std::invalid_argument("meshio++: point set index is out of range");
            }
            labels[static_cast<std::size_t>(index)] = static_cast<std::int64_t>(label);
        }
    }
    if (std::find(labels.begin(), labels.end(), -1) != labels.end())
        log::warn("sets_to_data: not all entities belong to a set; using default value -1");
    Mesh out = detail::clone_mesh(rMesh);
    if (Location == DataLocation::Point)
        out.AddPointData(key, dmanage_indices(labels));
    else {
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            NDArray data(DType::Int64, {rMesh.Cells(b).NumCells()});
            if (data.Size())
                std::memcpy(data.Data(), labels.data() + bases[b],
                            data.Size() * sizeof(std::int64_t));
            blocks.push_back(std::move(data));
        }
        out.AddCellData(key, std::move(blocks));
    }
    for (std::size_t i = out.NumRegions(); i > 0; --i)
        if (out.Region(i - 1).mKind == kind)
            out.RemoveRegion(i - 1);
    return out;
}

Mesh data_to_sets(const Mesh& rMesh, DataLocation Location, const std::string& rKey) {
    const auto kind = dmanage_set_kind(Location);
    dmanage_require(rMesh, Location, rKey, false);
    const auto bases = detail::block_bases(rMesh);
    const std::size_t num_blocks = Location == DataLocation::Point ? 1 : rMesh.NumCellBlocks();
    if (Location == DataLocation::Cell && rMesh.CellDataNumBlocks(rKey) != num_blocks)
        throw std::invalid_argument("meshio++: cell data must have one array per block");
    std::map<DmanageTag, std::vector<std::int64_t>> tags;
    for (std::size_t b = 0; b < num_blocks; ++b) {
        const auto& array =
            Location == DataLocation::Point ? rMesh.PointData(rKey) : rMesh.CellData(rKey, b);
        const auto rows =
            Location == DataLocation::Point ? rMesh.NumPoints() : rMesh.Cells(b).NumCells();
        if (array.Ndim() != 1 || array.Size() != rows)
            throw Unsupported("meshio++: data_to_sets requires scalar one-dimensional fields");
        detail::dispatch_dtype(array.Dtype(), [&]<class T>() {
            if constexpr (!std::is_integral_v<T>)
                throw std::invalid_argument("meshio++: data_to_sets array '" + rKey +
                                            "' is not int data");
            else {
                const T* values = array.As<T>();
                const auto base = Location == DataLocation::Point ? 0 : bases[b];
                for (std::size_t i = 0; i < rows; ++i)
                    tags[dmanage_tag(values[i])].push_back(base + static_cast<std::int64_t>(i));
            }
        });
    }
    std::vector<std::string> names;
    std::unordered_set<std::string> seen;
    for (std::size_t start = 0;;) {
        const auto end = rKey.find('-', start);
        const auto name = rKey.substr(start, end == std::string::npos ? end : end - start);
        if (seen.insert(name).second)
            names.push_back(name);
        if (end == std::string::npos)
            break;
        start = end + 1;
    }
    if (names.size() != tags.size()) {
        names.clear();
        for (const auto& [tag, entries] : tags)
            names.push_back("set-" + (Location == DataLocation::Point ? std::string("key") : rKey) +
                            "-" + tag.Name());
    }
    Mesh out = detail::clone_mesh(
        rMesh, [&](DataLocation location, const std::string& rName, std::string&) {
            return location != Location || rName != rKey;
        });
    std::size_t i = 0;
    for (const auto& [tag, entries] : tags) {
        const auto& name = names[i++];
        Region region(name, kind, dmanage_indices(entries));
        for (std::size_t j = 0; j < out.NumRegions(); ++j)
            if (out.Region(j).mKind == kind && out.Region(j).mName == name) {
                region.mDim = out.Region(j).mDim;
                region.mTag = out.Region(j).mTag;
                out.RemoveRegion(j);
                break;
            }
        out.AddRegion(std::move(region));
    }
    return out;
}

}  // namespace meshioplusplus
