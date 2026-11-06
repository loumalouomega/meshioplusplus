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
 * @file detail/patch_output.hpp
 * @brief The patch labelling `segment` and `refine_patches` write onto their
 * output mesh: the per-cell patch id, the patch adjacency table and one
 * `<prefix>:<patch>` cell region per patch (replacing an earlier run's).
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

/// Write the patch labels onto `rOut` (a clone of `rMesh`): `rCellPatch` the
/// global block-major patch of each cell (-1 off the surface), `rBases` the
/// block bases, `rAdjacency` the shared-edge count of each patch pair and
/// `NumPatches` the patch count.
inline void patch_output_write(
    const Mesh& rMesh, Mesh& rOut, const std::vector<std::int64_t>& rBases,
    const std::vector<std::int64_t>& rCellPatch, std::size_t NumPatches,
    const std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t>& rAdjacency,
    const std::string& rPrefix) {
    std::vector<NDArray> blocks;
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const std::size_t base = static_cast<std::size_t>(rBases[bi]);
        const std::size_t nc = rMesh.Cells(bi).NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {nc});
        std::copy(rCellPatch.begin() + static_cast<std::ptrdiff_t>(base),
                  rCellPatch.begin() + static_cast<std::ptrdiff_t>(base + nc),
                  a.As<std::int64_t>());
        blocks.push_back(std::move(a));
    }
    rOut.AddCellData(kSegmentPatchName, std::move(blocks));

    NDArray table = NDArray::Uninit(DType::Int64, {rAdjacency.size(), std::size_t{3}});
    {
        std::int64_t* t = table.As<std::int64_t>();
        std::size_t r = 0;
        for (const auto& kv : rAdjacency) {
            t[r * 3 + 0] = kv.first.first;
            t[r * 3 + 1] = kv.first.second;
            t[r * 3 + 2] = kv.second;
            ++r;
        }
    }
    rOut.AddFieldData(kSegmentAdjacencyName, std::move(table));

    // Replace an earlier run's `<prefix>:<integer>` cell regions.
    const std::string head = rPrefix + ":";
    for (std::size_t i = rOut.NumRegions(); i-- > 0;) {
        const meshioplusplus::Region& r = rOut.Region(i);
        if (r.mKind != RegionKind::Cell || r.mName.size() <= head.size() ||
            r.mName.compare(0, head.size(), head) != 0)
            continue;
        const bool digits = std::all_of(r.mName.begin() + static_cast<std::ptrdiff_t>(head.size()),
                                        r.mName.end(), [](char c) { return c >= '0' && c <= '9'; });
        if (digits)
            rOut.RemoveRegion(i);
    }
    const std::size_t ncells = rCellPatch.size();
    std::vector<std::vector<std::int64_t>> members(NumPatches);
    for (std::size_t c = 0; c < ncells; ++c)
        if (rCellPatch[c] >= 0)
            members[static_cast<std::size_t>(rCellPatch[c])].push_back(
                static_cast<std::int64_t>(c));
    for (std::size_t p = 0; p < NumPatches; ++p) {
        NDArray entries = NDArray::Uninit(DType::Int64, {members[p].size()});
        std::copy(members[p].begin(), members[p].end(), entries.As<std::int64_t>());
        rOut.AddRegion(meshioplusplus::Region(head + std::to_string(p), RegionKind::Cell, 2, -1,
                                              std::move(entries)));
    }
}

}  // namespace detail
}  // namespace meshioplusplus
