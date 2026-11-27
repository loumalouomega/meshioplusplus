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
#pragma once

/**
 * @file detail/region_tree.hpp
 * @brief Region hierarchy helpers (`/` path names + `field_data["meshio:tree"]`, roadmap §8.5).
 *
 * Regions remain flat on `Mesh`, with `/` separating hierarchy levels (e.g.
 * `Vehicle/Wing/Skin`). Internal assembly nodes that own no direct entities
 * are recorded as a UTF-8 string list under `field_data["meshio:tree"]` (using
 * `detail/string_field.hpp`).
 *
 * Core-private: no installed header includes this.
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "string_field.hpp"

namespace meshioplusplus::detail {

inline constexpr std::string_view kRegionTreeKey = "meshio:tree";

/// Split a `/`-separated hierarchy path into non-empty segments.
inline std::vector<std::string> split_region_path(std::string_view rPath) {
    std::vector<std::string> parts;
    std::size_t pos = 0;
    while (pos < rPath.size()) {
        const std::size_t slash = rPath.find('/', pos);
        const std::string_view seg =
            slash == std::string_view::npos ? rPath.substr(pos) : rPath.substr(pos, slash - pos);
        if (!seg.empty())
            parts.emplace_back(seg);
        if (slash == std::string_view::npos)
            break;
        pos = slash + 1;
    }
    return parts;
}

/// Join path segments with `/`.
inline std::string join_region_path(const std::vector<std::string>& rParts, std::size_t Count) {
    std::string out;
    for (std::size_t i = 0; i < Count && i < rParts.size(); ++i) {
        if (i != 0)
            out += '/';
        out += rParts[i];
    }
    return out;
}

/// Parent path of a `/`-separated name (`"a/b/c"` -> `"a/b"`, `"a"` -> `""`).
inline std::string parent_region_path(std::string_view rPath) {
    const auto parts = split_region_path(rPath);
    return parts.size() <= 1 ? std::string() : join_region_path(parts, parts.size() - 1);
}

/// Leaf segment of a `/`-separated name (`"a/b/c"` -> `"c"`).
inline std::string leaf_region_name(std::string_view rPath) {
    const auto parts = split_region_path(rPath);
    return parts.empty() ? std::string() : parts.back();
}

/// Store sorted, deduplicated empty internal assembly paths in `field_data["meshio:tree"]`.
inline void set_region_tree_nodes(Mesh& rMesh, const std::vector<std::string>& rPaths) {
    std::set<std::string> uniq;
    for (const auto& p : rPaths) {
        const auto segs = split_region_path(p);
        if (!segs.empty())
            uniq.insert(join_region_path(segs, segs.size()));
    }
    if (uniq.empty())
        return;
    std::vector<std::string> sorted(uniq.begin(), uniq.end());
    set_field_strings(rMesh, std::string(kRegionTreeKey), sorted);
}

/// Read internal assembly paths from `field_data["meshio:tree"]` (empty if absent).
inline std::vector<std::string> get_region_tree_nodes(const Mesh& rMesh) {
    auto opt = get_field_strings(rMesh, std::string(kRegionTreeKey));
    return opt ? std::move(*opt) : std::vector<std::string>{};
}

struct RegionTreeNode {
    std::string mPath;
    std::string mName;
    std::vector<std::string> mChildren;       ///< Direct child full paths, insertion/sorted order.
    std::vector<std::size_t> mRegionIndices;  ///< Indices into `rMesh.Regions()` at this node.
};

struct RegionTree {
    std::vector<std::string> mRoots;
    std::unordered_map<std::string, RegionTreeNode> mNodes;

    bool HasHierarchy() const {
        for (const auto& [path, node] : mNodes)
            if (!node.mChildren.empty() || path.find('/') != std::string::npos)
                return true;
        return false;
    }
};

/// Build the region hierarchy tree from `rMesh.Regions()` and `field_data["meshio:tree"]`.
inline RegionTree build_region_tree(const Mesh& rMesh) {
    RegionTree tree;
    auto ensure_path = [&](const std::string& rPath) -> std::string {
        const auto parts = split_region_path(rPath);
        if (parts.empty())
            return {};
        std::string cur;
        for (std::size_t i = 0; i < parts.size(); ++i) {
            const std::string parent = cur;
            cur = i == 0 ? parts[i] : cur + "/" + parts[i];
            auto it = tree.mNodes.find(cur);
            if (it == tree.mNodes.end()) {
                RegionTreeNode node;
                node.mPath = cur;
                node.mName = parts[i];
                tree.mNodes.emplace(cur, std::move(node));
                if (parent.empty())
                    tree.mRoots.push_back(cur);
                else
                    tree.mNodes.at(parent).mChildren.push_back(cur);
            }
        }
        return cur;
    };

    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const std::string canonical = ensure_path(rMesh.Region(i).mName);
        if (!canonical.empty())
            tree.mNodes.at(canonical).mRegionIndices.push_back(i);
    }
    for (const std::string& p : get_region_tree_nodes(rMesh))
        (void)ensure_path(p);

    return tree;
}

}  // namespace meshioplusplus::detail
