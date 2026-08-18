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

// System includes
#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <utility>

// Project includes
#include "meshioplusplus/detail/refine_hierarchy.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/refine.hpp"

namespace meshioplusplus {
namespace detail {

RefineHierarchyState refine_read_hierarchy(const Mesh& rMesh,
                                           const std::vector<std::int64_t>& rBases,
                                           std::vector<std::int64_t>& rIds, std::int64_t& rIdBase) {
    if (!rMesh.HasCellData(kRefineCellIdName))
        return RefineHierarchyState::Absent;

    const std::size_t nblocks = rMesh.NumCellBlocks();
    const bool has_parent = rMesh.HasCellData(kRefineParentIdName);
    if (rMesh.CellDataNumBlocks(kRefineCellIdName) != nblocks ||
        (has_parent && rMesh.CellDataNumBlocks(kRefineParentIdName) != nblocks)) {
        log::warn("refine: ignoring '{}'/'{}': they do not cover every cell block.",
                  kRefineCellIdName, kRefineParentIdName);
        return RefineHierarchyState::Invalid;
    }

    const std::int64_t total = total_cells(rBases);
    std::vector<std::int64_t> ids(static_cast<std::size_t>(total));
    std::vector<std::int64_t> parent_ids;
    if (has_parent)
        parent_ids.resize(static_cast<std::size_t>(total));

    std::size_t bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t base = static_cast<std::size_t>(rBases[bi]);
        const std::size_t ncells = cb.NumCells();
        const NDArray& a = rMesh.CellData(kRefineCellIdName, bi);
        if (rows(a) != ncells || (ncells != 0 && a.Size() / ncells != 1)) {
            log::warn("refine: ignoring '{}': block {} is not one scalar value per cell.",
                      kRefineCellIdName, bi);
            return RefineHierarchyState::Invalid;
        }
        for (std::size_t c = 0; c < ncells; ++c)
            ids[base + c] = read_int(a, c);
        if (has_parent) {
            const NDArray& p = rMesh.CellData(kRefineParentIdName, bi);
            if (rows(p) != ncells || (ncells != 0 && p.Size() / ncells != 1)) {
                log::warn("refine: ignoring '{}': block {} is not one scalar value per cell.",
                          kRefineParentIdName, bi);
                return RefineHierarchyState::Invalid;
            }
            for (std::size_t c = 0; c < ncells; ++c)
                parent_ids[base + c] = read_int(p, c);
        }
        ++bi;
    }

    std::int64_t max_id = -1;
    std::unordered_map<std::int64_t, char> seen;
    seen.reserve(ids.size() * 2);
    for (std::int64_t id : ids) {
        if (id < 0) {
            log::warn("refine: ignoring '{}'/'{}': ids must be non-negative.", kRefineCellIdName,
                      kRefineParentIdName);
            return RefineHierarchyState::Invalid;
        }
        if (!seen.emplace(id, 0).second) {
            log::warn(
                "refine: ignoring '{}'/'{}': the id {} is not unique, so the mesh was merged or "
                "the array was replicated by another operation.",
                kRefineCellIdName, kRefineParentIdName, id);
            return RefineHierarchyState::Invalid;
        }
        max_id = std::max(max_id, id);
    }
    for (std::int64_t id : parent_ids)
        max_id = std::max(max_id, id);

    rIds = std::move(ids);
    rIdBase = max_id + 1;
    return RefineHierarchyState::Valid;
}

}  // namespace detail
}  // namespace meshioplusplus
