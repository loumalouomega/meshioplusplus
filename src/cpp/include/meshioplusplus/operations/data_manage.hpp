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
 * @file operations/data_manage.hpp
 * @brief Array management for a mesh's data: rename, drop, and keep-only.
 *
 * `data_manage` rewrites *which* data arrays a mesh carries and under what
 * names. Values are copied verbatim — no array's contents, dtype or shape is
 * ever reinterpreted — and the geometry (points and every cell block) comes
 * through bit-identical.
 *
 * The three phases run in a fixed, documented order: **keep**, then **drop**,
 * then **rename**. `keep` is a whitelist that only applies to the locations it
 * actually mentions, so `keep({{Point, "T"}})` prunes `point_data` down to `T`
 * while leaving `cell_data` and `field_data` completely untouched.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format — it is not
 * in the format registry.
 */

// System includes
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {

/// One (location, name) key identifying a single data array.
struct DataKey {
    DataLocation location = DataLocation::Point;  ///< Which data map.
    std::string name;                             ///< The array name.
};

/// One rename request: `from` -> `to`, within a single location.
struct DataRename {
    DataLocation location = DataLocation::Point;  ///< Which data map.
    std::string from;                             ///< Existing array name.
    std::string to;                               ///< New array name.
};

/// Options for `data_manage` (see the file comment for the phase order).
struct DataManageOptions {
    /// Whitelist, applied first. Only affects the locations that appear in it:
    /// a location mentioned by no entry passes through untouched.
    std::vector<DataKey> keep;
    /// Arrays to remove, applied after `keep`.
    std::vector<DataKey> drop;
    /// Renames, applied last.
    std::vector<DataRename> rename;
    /// Skip entries whose key does not exist instead of throwing.
    bool ignore_missing = false;
};

/// Result of `data_manage`.
struct DataManageResult {
    Mesh mMesh;  ///< The rewritten mesh (geometry bit-identical to the input).
    /// Arrays that were removed, as `"point_data:T"`, sorted.
    std::vector<std::string> mDropped;
    /// Renames that were applied, as `("point_data:old", "point_data:new")`.
    std::vector<std::pair<std::string, std::string>> mRenamed;
};

/**
 * @brief Rewrites which data arrays @p rMesh carries.
 *
 * All validation runs against the *input* mesh before anything is cloned, so a
 * bad key costs nothing.
 * @param rMesh the mesh to rewrite (unmodified).
 * @param rOpts the keep / drop / rename requests.
 * @return the rewritten mesh plus what was dropped and renamed.
 * @throws std::invalid_argument on an unknown key (unless `ignore_missing`), on
 *         a rename whose target already exists and is not itself being renamed
 *         away, on two renames targeting the same name, or on two renames of
 *         the same source.
 */
MESHIOPLUSPLUS_API DataManageResult data_manage(const Mesh& rMesh, const DataManageOptions& rOpts);

/**
 * @brief Drops the named arrays at one location.
 * @param rMesh the mesh to rewrite (unmodified).
 * @param Location which data map to drop from.
 * @param rNames the array names to remove.
 * @param IgnoreMissing skip names that do not exist instead of throwing.
 * @return the rewritten mesh.
 * @throws std::invalid_argument on an unknown name (unless @p IgnoreMissing).
 */
MESHIOPLUSPLUS_API Mesh data_drop(const Mesh& rMesh, DataLocation Location, const std::vector<std::string>& rNames,
               bool IgnoreMissing = false);

/**
 * @brief Keeps only the named arrays at one location, dropping the rest there.
 *
 * The other two locations are left untouched.
 * @param rMesh the mesh to rewrite (unmodified).
 * @param Location which data map to prune.
 * @param rNames the array names to retain.
 * @param IgnoreMissing skip names that do not exist instead of throwing.
 * @return the rewritten mesh.
 * @throws std::invalid_argument on an unknown name (unless @p IgnoreMissing).
 */
MESHIOPLUSPLUS_API Mesh data_keep(const Mesh& rMesh, DataLocation Location, const std::vector<std::string>& rNames,
               bool IgnoreMissing = false);

/**
 * @brief Renames one array, preserving its values, dtype and shape.
 * @param rMesh the mesh to rewrite (unmodified).
 * @param Location which data map the array lives in.
 * @param rFrom the existing name.
 * @param rTo the new name.
 * @return the rewritten mesh.
 * @throws std::invalid_argument if @p rFrom does not exist or @p rTo already does.
 */
MESHIOPLUSPLUS_API Mesh data_rename(const Mesh& rMesh, DataLocation Location, const std::string& rFrom,
                 const std::string& rTo);

}  // namespace meshioplusplus
