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
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/operations/data_common.hpp"

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

}  // namespace meshioplusplus
