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
// The shared mesh-carving helper. See detail/mesh_carve.hpp for the contract; this
// used to be VTKHDF's own private vtkhdf_carve(), generalized so EnSight Gold's
// multi-part writer can reuse it rather than reimplementing the same region-vs-block
// fallback logic (roadmap §1.1).

// System includes
#include <algorithm>
#include <cstddef>
#include <set>
#include <string>
#include <tuple>
#include <vector>

// Project includes
#include "meshioplusplus/detail/mesh_carve.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

std::vector<MeshPart> carve_by_region(const Mesh& rMesh, const std::string& rCaller,
                                      bool StrictNames, bool RejectSlash) {
    const std::vector<std::int64_t> bases = block_bases(rMesh);
    const std::size_t total = static_cast<std::size_t>(bases.back());
    if (total == 0)
        throw WriteError("meshio++: " + rCaller +
                         ": a multi-part write needs at least one cell to carve into parts");

    std::vector<const meshioplusplus::Region*> regions;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        if (rMesh.Region(i).mKind == RegionKind::Cell)
            regions.push_back(&rMesh.Region(i));
    // Parts are written in (tag, name) order -- the mesh's own region storage order.
    std::stable_sort(regions.begin(), regions.end(),
                     [](const meshioplusplus::Region* a, const meshioplusplus::Region* b) {
                         return std::tie(a->mTag, a->mName) < std::tie(b->mTag, b->mName);
                     });
    if (!regions.empty()) {
        std::vector<char> hit(total, 0);
        std::size_t covered = 0;
        bool ok = true;
        for (const auto* r : regions)
            for (std::size_t i = 0; i < r->NumEntries(); ++i) {
                const std::int64_t g = r->Entries()[i];
                if (g < 0 || static_cast<std::size_t>(g) >= total ||
                    hit[static_cast<std::size_t>(g)]) {
                    ok = false;
                    continue;
                }
                hit[static_cast<std::size_t>(g)] = 1;
                ++covered;
            }
        bool names_ok = true;
        if (ok && covered == total) {
            std::vector<MeshPart> parts;
            std::set<std::string> names;
            for (std::size_t i = 0; i < regions.size(); ++i) {
                MeshPart p;
                p.mName =
                    regions[i]->mName.empty() ? "block_" + std::to_string(i) : regions[i]->mName;
                const bool has_slash = RejectSlash && p.mName.find('/') != std::string::npos;
                if (has_slash || !names.insert(p.mName).second) {
                    if (StrictNames)
                        throw WriteError("meshio++: " + rCaller +
                                         ": cell region names must be unique" +
                                         (RejectSlash ? " and free of '/'" : "") +
                                         " to name parts");
                    names_ok = false;
                    break;
                }
                for (std::size_t e = 0; e < regions[i]->NumEntries(); ++e)
                    p.mCells.push_back(static_cast<std::size_t>(regions[i]->Entries()[e]));
                std::sort(p.mCells.begin(), p.mCells.end());
                parts.push_back(std::move(p));
            }
            if (names_ok)
                return parts;
            log::warn(
                "meshio++: {}: cell regions do not have unique names; writing one part per "
                "cell block instead.",
                rCaller);
        } else {
            log::warn(
                "meshio++: {}: cell regions overlap or do not cover every cell; writing one "
                "part per cell block instead.",
                rCaller);
        }
    }

    std::vector<MeshPart> parts;
    for (std::size_t bi = 0; bi + 1 < bases.size(); ++bi) {
        if (bases[bi + 1] == bases[bi])
            continue;
        MeshPart p;
        p.mName = "block_" + std::to_string(bi);
        for (std::int64_t g = bases[bi]; g < bases[bi + 1]; ++g)
            p.mCells.push_back(static_cast<std::size_t>(g));
        parts.push_back(std::move(p));
    }
    return parts;
}

}  // namespace detail
}  // namespace meshioplusplus
