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
// Green-element undo: lookup and substitution against a caller-supplied coarse
// mesh -- see operations/undo_green.hpp for the contract and the reasoning
// for why no subdivision-table inversion is needed here.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/undo_green.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/refine_hierarchy.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kUgPrefix = "meshio++: undo_green: ";

// A required scalar cell_data array (one value per cell, covering every
// block), read in place: one view per block, zero-copy for the Int64 arrays
// refine() writes (a converted copy of that block otherwise), indexed by global
// cell. The flat copies this replaces were 2 GB at 83M fine cells (roadmap
// 6.1.1). Throws by name on any mismatch -- undo_green's fine-mesh
// preconditions are hard requirements, unlike detail::refine_read_hierarchy's
// warn-and-fall-back contract.
class UgColumn {
public:
    UgColumn(const Mesh& rMesh, const std::vector<std::int64_t>& rBases, const char* pName)
        : mrBases(rBases) {
        if (!rMesh.HasCellData(pName))
            throw std::invalid_argument(
                std::string(kUgPrefix) + "the fine mesh has no '" + pName +
                "' cell_data; run refine(..., record_hierarchy=True, record_levels=True) first");
        const std::size_t nblocks = rMesh.NumCellBlocks();
        if (rMesh.CellDataNumBlocks(pName) != nblocks)
            throw std::invalid_argument(std::string(kUgPrefix) + "'" + pName +
                                        "' does not cover every cell block");
        mViews = std::vector<std::optional<detail::Int64View>>(nblocks);
        mBlocks.assign(nblocks, nullptr);
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const std::size_t ncells = cb.NumCells();
            const NDArray& a = rMesh.CellData(pName, bi);
            if (detail::rows(a) != ncells || (ncells != 0 && a.Size() / ncells != 1))
                throw std::invalid_argument(std::string(kUgPrefix) + "'" + pName + "' block " +
                                            std::to_string(bi) +
                                            " is not one scalar value per cell");
            mViews[bi].emplace(a);
            mBlocks[bi] = mViews[bi]->Data();
            ++bi;
        }
    }
    UgColumn(const UgColumn&) = delete;
    UgColumn& operator=(const UgColumn&) = delete;

    /// Block @p bi's values, row-indexed.
    const std::int64_t* Block(std::size_t bi) const { return mBlocks[bi]; }

    /// The value of global cell @p g.
    std::int64_t At(std::size_t g) const {
        if (mBlocks.size() == 1)
            return mBlocks[0][g];
        const auto it =
            std::upper_bound(mrBases.begin() + 1, mrBases.end(), static_cast<std::int64_t>(g));
        const auto bi = static_cast<std::size_t>(it - mrBases.begin() - 1);
        return mBlocks[bi][g - static_cast<std::size_t>(mrBases[bi])];
    }

private:
    const std::vector<std::int64_t>& mrBases;
    std::vector<std::optional<detail::Int64View>> mViews;
    std::vector<const std::int64_t*> mBlocks;
};

// Same shape check, but returns empty rather than throwing when `pName` is
// absent or malformed -- used for the coarse mesh's OPTIONAL refine:level,
// where absence means "never refined" (implicit level 0 for every cell).
std::vector<std::int64_t> ug_read_optional(const Mesh& rMesh,
                                           const std::vector<std::int64_t>& rBases,
                                           const char* pName) {
    if (!rMesh.HasCellData(pName))
        return {};
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (rMesh.CellDataNumBlocks(pName) != nblocks)
        return {};
    std::vector<std::int64_t> out(static_cast<std::size_t>(detail::total_cells(rBases)));
    std::size_t bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t base = static_cast<std::size_t>(rBases[bi]);
        const std::size_t ncells = cb.NumCells();
        const NDArray& a = rMesh.CellData(pName, bi);
        if (detail::rows(a) != ncells || (ncells != 0 && a.Size() / ncells != 1))
            return {};
        for (std::size_t c = 0; c < ncells; ++c)
            out[base + c] = detail::read_int(a, c);
        ++bi;
    }
    return out;
}

// Whether a name is one of the six reserved refine:* arrays -- always
// excluded from undo_green's output (stale hierarchy bookkeeping).
bool ug_is_reserved(const std::string& rName) {
    return rName == kRefineParentCellName || rName == kRefineLevelName ||
           rName == kRefineHangingName || rName == kRefineEntityName ||
           rName == kRefineCellIdName || rName == kRefineParentIdName;
}

// A cell's role in the output, decided once per sibling group.
enum class UgRole : std::uint8_t { Keep, AnchorSubstitute, Suppressed };

}  // namespace

UndoGreenResult undo_green(const Mesh& rCoarse, const Mesh& rFine) {
    if (rCoarse.NumPoints() > rFine.NumPoints())
        throw std::invalid_argument(
            std::string(kUgPrefix) +
            "the coarse mesh has more points than the fine mesh, so they cannot be the "
            "coarse/fine pair of one refine() call");

    const std::vector<std::int64_t> fine_bases = detail::block_bases(rFine);
    const std::vector<std::int64_t> coarse_bases = detail::block_bases(rCoarse);
    const auto total_fine = static_cast<std::size_t>(detail::total_cells(fine_bases));

    // What the classification leaves for the output: every global fine cell's
    // role and green group, and per group its coarse row. Everything else it
    // reads (the id lookup, the sibling order) lives in the scope below and is
    // released before the output is built; the three refine:* columns are read
    // in place. A group id is read only for a substituted or suppressed cell,
    // so 32 bits each is enough: a green group has at least two cells.
    std::vector<UgRole> role(total_fine, UgRole::Keep);
    std::vector<std::uint32_t> group_id_of(total_fine, 0);
    std::vector<std::int64_t> group_coarse_row;  // indexed by green group id
    std::int64_t num_removed = 0;
    {
        const UgColumn fine_id(rFine, fine_bases, kRefineCellIdName);
        const UgColumn fine_parent(rFine, fine_bases, kRefineParentIdName);
        const UgColumn fine_level(rFine, fine_bases, kRefineLevelName);

        // --- resolve the coarse mesh's id space: its own recorded ids if valid,
        // else the implicit global-block-major ids -- refine_attach_hierarchy's
        // own "start fresh" fallback, where an id is its row and needs no table.
        const auto total_coarse = detail::total_cells(coarse_bases);
        std::vector<std::int64_t> coarse_ids;
        std::int64_t coarse_id_base = 0;
        const bool coarse_ids_recorded =
            detail::refine_read_hierarchy(rCoarse, coarse_bases, coarse_ids, coarse_id_base) ==
            detail::RefineHierarchyState::Valid;
        // (id, row) pairs sorted -- a total order -- searched by binary search.
        std::vector<std::pair<std::int64_t, std::int64_t>> id_to_coarse_row;
        if (coarse_ids_recorded) {
            id_to_coarse_row.resize(coarse_ids.size());
            parallel_for(coarse_ids.size(), [&](std::size_t i) {
                id_to_coarse_row[i] = {coarse_ids[i], static_cast<std::int64_t>(i)};
            });
            parallel_sort(id_to_coarse_row.begin(), id_to_coarse_row.end());
        }
        coarse_ids = std::vector<std::int64_t>();
        // The coarse row of an id, or -1. With recorded ids, the last row with
        // this id, as the former map's overwrite kept it.
        const auto coarse_row_of = [&](std::int64_t id) -> std::int64_t {
            if (!coarse_ids_recorded)
                return (id >= 0 && id < total_coarse) ? id : -1;
            auto it = std::upper_bound(
                id_to_coarse_row.begin(), id_to_coarse_row.end(), id,
                [](std::int64_t v, const std::pair<std::int64_t, std::int64_t>& rE) {
                    return v < rE.first;
                });
            if (it != id_to_coarse_row.begin())
                --it;
            if (it == id_to_coarse_row.end() || it->first != id)
                return -1;
            return it->second;
        };

        const std::vector<std::int64_t> coarse_level =
            ug_read_optional(rCoarse, coarse_bases, kRefineLevelName);
        auto coarse_level_at = [&](std::int64_t row) -> std::int64_t {
            return coarse_level.empty() ? 0 : coarse_level[static_cast<std::size_t>(row)];
        };

        // --- group fine cells by parent_id --------------------------------------
        // Sorted (parent_id, cell) pairs: each run is one sibling group, its cells
        // ascending, and the groups come in ascending parent_id -- the same on
        // every platform, where a hash map's iteration order was not (it decided
        // which of several malformed groups was reported).
        // refine() writes its cells parent-major, so parent_id usually arrives
        // non-decreasing: the sorted pairs are then (parent_id[g], g) for g in
        // order, and the sort (and its N pairs) is skipped.
        bool parents_sorted = true;
        {
            const std::int64_t* p_prev = nullptr;
            for (std::size_t bi = 0; parents_sorted && bi + 1 < fine_bases.size(); ++bi) {
                const auto n = static_cast<std::size_t>(fine_bases[bi + 1] - fine_bases[bi]);
                if (n == 0)
                    continue;
                const std::int64_t* p = fine_parent.Block(bi);
                parents_sorted = (p_prev == nullptr || *p_prev <= p[0]) && std::is_sorted(p, p + n);
                p_prev = p + n - 1;
            }
        }
        std::vector<std::pair<std::int64_t, std::int64_t>> by_parent;
        if (!parents_sorted) {
            by_parent.resize(total_fine);
            for (std::size_t bi = 0; bi + 1 < fine_bases.size(); ++bi) {
                const auto base = static_cast<std::size_t>(fine_bases[bi]);
                const auto n = static_cast<std::size_t>(fine_bases[bi + 1]) - base;
                const std::int64_t* p = fine_parent.Block(bi);
                parallel_for(n, [&](std::size_t r) {
                    by_parent[base + r] = {p[r], static_cast<std::int64_t>(base + r)};
                });
            }
            parallel_sort(by_parent.begin(), by_parent.end());
        }
        const auto parent_at = [&](std::size_t k) {
            return parents_sorted ? fine_parent.At(k) : by_parent[k].first;
        };
        const auto cell_at = [&](std::size_t k) {
            return parents_sorted ? static_cast<std::int64_t>(k) : by_parent[k].second;
        };

        // --- classify every global fine cell's role; for green groups, resolve
        // the substitution source once ------------------------------------------
        std::vector<std::int64_t> members;
        for (std::size_t r0 = 0, r1 = 0; r0 < total_fine; r0 = r1) {
            const std::int64_t parent_id = parent_at(r0);
            members.clear();
            for (r1 = r0; r1 < total_fine && parent_at(r1) == parent_id; ++r1)
                members.push_back(cell_at(r1));

            if (members.size() == 1) {
                const std::int64_t g = members.front();
                if (fine_id.At(static_cast<std::size_t>(g)) != parent_id)
                    throw std::invalid_argument(
                        std::string(kUgPrefix) +
                        "malformed hierarchy: a singleton sibling group's refine:cell_id does not "
                        "equal its refine:parent_id (cell " +
                        std::to_string(g) + ")");
                continue;  // untouched: role stays Keep
            }

            const std::int64_t coarse_row = coarse_row_of(parent_id);
            if (coarse_row < 0)
                throw std::invalid_argument(
                    std::string(kUgPrefix) + "refine:parent_id " + std::to_string(parent_id) +
                    " does not resolve in the coarse mesh's id space -- these two meshes are not "
                    "the coarse/fine pair of one refine() call");
            const std::int64_t clevel = coarse_level_at(coarse_row);

            const std::int64_t flevel = fine_level.At(static_cast<std::size_t>(members.front()));
            for (std::int64_t g : members)
                if (fine_level.At(static_cast<std::size_t>(g)) != flevel)
                    throw std::invalid_argument(
                        std::string(kUgPrefix) +
                        "malformed hierarchy: the sibling group under refine:parent_id " +
                        std::to_string(parent_id) + " does not agree on refine:level");

            if (flevel == clevel + 1)
                continue;  // red: a genuine refinement, kept unchanged (role stays Keep)

            if (flevel > clevel + 1)
                throw std::invalid_argument(
                    std::string(kUgPrefix) + "the sibling group under refine:parent_id " +
                    std::to_string(parent_id) + " has refine:level " + std::to_string(flevel) +
                    ", more than one deeper than its coarse parent's own level (" +
                    std::to_string(clevel) +
                    "); undo_green only supports a single-pass (levels=1) hierarchy");

            if (flevel != clevel)
                throw std::invalid_argument(
                    std::string(kUgPrefix) + "the sibling group under refine:parent_id " +
                    std::to_string(parent_id) + " has refine:level " + std::to_string(flevel) +
                    ", which is neither its coarse parent's own level (" + std::to_string(clevel) +
                    ", green) nor one more (" + std::to_string(clevel + 1) + ", red)");

            // green: substitute the whole group with one row from the coarse mesh
            if (group_coarse_row.size() >= std::numeric_limits<std::uint32_t>::max())
                throw Unsupported(std::string(kUgPrefix) + "more than " +
                                  std::to_string(std::numeric_limits<std::uint32_t>::max()) +
                                  " green groups");
            const auto gid = static_cast<std::uint32_t>(group_coarse_row.size());
            group_coarse_row.push_back(coarse_row);
            num_removed += static_cast<std::int64_t>(members.size()) - 1;
            role[static_cast<std::size_t>(members.front())] = UgRole::AnchorSubstitute;
            for (std::size_t k = 1; k < members.size(); ++k)
                role[static_cast<std::size_t>(members[k])] = UgRole::Suppressed;
            for (std::int64_t g : members)
                group_id_of[static_cast<std::size_t>(g)] = gid;
        }
    }  // the classification scratch is released here

    // --- build the output mesh: fine's own block structure, unchanged types
    // and order, rows compacted/substituted -----------------------------------
    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rFine.Points()));

    const std::size_t nblocks = rFine.NumCellBlocks();
    std::vector<std::size_t> out_ncells(nblocks, 0);
    {
        std::size_t bi = 0;
        for (const auto cb : rFine.CellRange()) {
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::size_t g = static_cast<std::size_t>(fine_bases[bi]) + r;
                if (role[g] != UgRole::Suppressed)
                    ++out_ncells[bi];
            }
            ++bi;
        }
    }

    std::vector<NDArray> cell_maps(nblocks);
    std::vector<NDArray> out_conn(nblocks);
    const std::size_t num_groups = group_coarse_row.size();
    std::vector<std::int64_t> group_output_global(num_groups, -1);
    std::vector<std::uint32_t> group_fine_block(num_groups, 0);

    {
        std::size_t bi = 0;
        std::int64_t out_block_base = 0;
        for (const auto cb : rFine.CellRange()) {
            const std::size_t npc = cb.NodesPerCell();
            NDArray conn = NDArray::Uninit(DType::Int64, {out_ncells[bi], npc});
            std::int64_t* dst = conn.As<std::int64_t>();
            NDArray cmap = NDArray::Uninit(DType::Int64, {cb.NumCells()});
            std::int64_t* cm = cmap.As<std::int64_t>();

            std::size_t out_row = 0;
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::size_t g = static_cast<std::size_t>(fine_bases[bi]) + r;
                if (role[g] == UgRole::Suppressed) {
                    cm[r] = group_output_global[group_id_of[g]];
                    continue;
                }
                if (role[g] == UgRole::Keep) {
                    const NDArray& src = cb.Conn();
                    for (std::size_t k = 0; k < npc; ++k)
                        dst[out_row * npc + k] = detail::read_int(src, r * npc + k);
                } else {  // AnchorSubstitute
                    const std::uint32_t gid = group_id_of[g];
                    const std::int64_t coarse_row = group_coarse_row[gid];
                    const auto blk_row = detail::global_to_block_row(coarse_bases, coarse_row);
                    const auto ccb = rCoarse.Cells(blk_row.first);
                    const NDArray& csrc = ccb.Conn();
                    for (std::size_t k = 0; k < npc; ++k)
                        dst[out_row * npc + k] = detail::read_int(
                            csrc, static_cast<std::size_t>(blk_row.second) * npc + k);
                    group_output_global[gid] = out_block_base + static_cast<std::int64_t>(out_row);
                    group_fine_block[gid] = static_cast<std::uint32_t>(bi);
                }
                cm[r] = out_block_base + static_cast<std::int64_t>(out_row);
                ++out_row;
            }

            out_conn[bi] = std::move(conn);
            cell_maps[bi] = std::move(cmap);
            out_block_base += static_cast<std::int64_t>(out_ncells[bi]);
            ++bi;
        }
    }

    {
        std::size_t bi = 0;
        for (const auto cb : rFine.CellRange()) {
            out.AddCellBlock(std::string(cb.Type()), std::move(out_conn[bi]));
            ++bi;
        }
    }

    // --- point_data / field_data: unchanged, no new or pruned points --------
    // (refine:entity / refine:hanging are POINT data, not cell_data, but are
    // just as much stale hierarchy bookkeeping as the four cell_data ones --
    // ug_is_reserved covers all six, so this loop must consult it too)
    for (const std::string& name : rFine.PointDataNames()) {
        if (ug_is_reserved(name))
            continue;
        out.AddPointData(name, detail::data_owned_copy(rFine.PointData(name)));
    }
    for (const std::string& name : rFine.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rFine.FieldData(name)));

    // --- cell_data: reserved refine:* arrays dropped; everything else kept
    // (fine's own row for Keep, coarse's for AnchorSubstitute) ---------------
    for (const std::string& name : rFine.CellDataNames()) {
        if (ug_is_reserved(name))
            continue;
        if (rFine.CellDataNumBlocks(name) != nblocks) {
            log::warn(
                "{}cell_data '{}' does not have one array per fine block; dropped rather than "
                "guessed at",
                kUgPrefix, name);
            continue;
        }

        // Pre-check every green group's substitution source before touching
        // anything, so a bad array is dropped whole rather than half-built.
        bool ok = true;
        if (num_groups > 0) {
            if (!rCoarse.HasCellData(name) ||
                rCoarse.CellDataNumBlocks(name) != rCoarse.NumCellBlocks()) {
                ok = false;
            } else {
                for (std::size_t ug = 0; ok && ug < num_groups; ++ug) {
                    const auto blk_row =
                        detail::global_to_block_row(coarse_bases, group_coarse_row[ug]);
                    const NDArray& fsrc = rFine.CellData(name, group_fine_block[ug]);
                    const NDArray& csrc = rCoarse.CellData(name, blk_row.first);
                    const std::size_t f_row_bytes =
                        fsrc.Nbytes() / std::max<std::size_t>(1, detail::rows(fsrc));
                    const std::size_t c_row_bytes =
                        csrc.Nbytes() / std::max<std::size_t>(1, detail::rows(csrc));
                    if (fsrc.Dtype() != csrc.Dtype() || f_row_bytes != c_row_bytes)
                        ok = false;
                }
            }
        }
        if (!ok) {
            log::warn(
                "{}cell_data '{}' cannot be honestly restored for a substituted cell (missing, "
                "incomplete or a different shape/dtype on the coarse mesh); dropped",
                kUgPrefix, name);
            continue;
        }

        std::vector<NDArray> out_blocks(nblocks);
        std::size_t bi = 0;
        for (const auto cb : rFine.CellRange()) {
            const NDArray& fsrc = rFine.CellData(name, bi);
            std::vector<std::size_t> shape = fsrc.Shape();
            shape[0] = out_ncells[bi];
            NDArray dst = NDArray::Uninit(fsrc.Dtype(), shape);
            const std::size_t row_bytes =
                fsrc.Nbytes() / std::max<std::size_t>(1, detail::rows(fsrc));
            std::byte* p_dst = dst.Data();
            const std::byte* p_fsrc = fsrc.Data();

            std::size_t out_row = 0;
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::size_t g = static_cast<std::size_t>(fine_bases[bi]) + r;
                if (role[g] == UgRole::Suppressed)
                    continue;
                if (role[g] == UgRole::Keep) {
                    std::memcpy(p_dst + out_row * row_bytes, p_fsrc + r * row_bytes, row_bytes);
                } else {  // AnchorSubstitute
                    const std::uint32_t ug = group_id_of[g];
                    const auto blk_row =
                        detail::global_to_block_row(coarse_bases, group_coarse_row[ug]);
                    const NDArray& csrc = rCoarse.CellData(name, blk_row.first);
                    std::memcpy(p_dst + out_row * row_bytes,
                                csrc.Data() + static_cast<std::size_t>(blk_row.second) * row_bytes,
                                row_bytes);
                }
                ++out_row;
            }
            out_blocks[bi] = std::move(dst);
            ++bi;
        }
        out.AddCellData(name, std::move(out_blocks));
    }

    UndoGreenResult res;
    res.mMesh = std::move(out);
    res.mCellMaps = std::move(cell_maps);
    res.mNumGroupsUndone = static_cast<std::int64_t>(num_groups);
    res.mNumCellsRemoved = num_removed;

    detail::RegionRemap rmap;
    rmap.mCellMapKind = detail::CellMapKind::Direct;
    rmap.pCellMaps = &res.mCellMaps;
    rmap.mOpName = "undo_green";
    detail::remap_regions(rFine, res.mMesh, rmap);

    return res;
}

namespace {

std::size_t ug_red_child_count(const std::string& rType) {
    if (rType == "line")
        return 2;
    if (rType == "triangle" || rType == "quad")
        return 4;
    if (rType == "tetra" || rType == "wedge" || rType == "hexahedron")
        return 8;
    return 0;
}

std::array<double, 3> ug_read_pt(const NDArray& rPts, std::size_t idx, std::size_t dim) {
    std::array<double, 3> p{0.0, 0.0, 0.0};
    for (std::size_t d = 0; d < std::min<std::size_t>(dim, 3); ++d)
        p[d] = detail::read_double(rPts, idx * dim + d);
    return p;
}

double ug_tet_signed_vol(const NDArray& rPts, std::size_t dim,
                         const std::vector<std::int64_t>& rN) {
    const auto a = ug_read_pt(rPts, static_cast<std::size_t>(rN[0]), dim);
    const auto b = ug_read_pt(rPts, static_cast<std::size_t>(rN[1]), dim);
    const auto c = ug_read_pt(rPts, static_cast<std::size_t>(rN[2]), dim);
    const auto d = ug_read_pt(rPts, static_cast<std::size_t>(rN[3]), dim);
    const double ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
    const double vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2];
    const double wx = d[0] - a[0], wy = d[1] - a[1], wz = d[2] - a[2];
    return ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx);
}

std::array<double, 3> ug_tri_normal(const NDArray& rPts, std::size_t dim, std::int64_t i0,
                                    std::int64_t i1, std::int64_t i2) {
    const auto a = ug_read_pt(rPts, static_cast<std::size_t>(i0), dim);
    const auto b = ug_read_pt(rPts, static_cast<std::size_t>(i1), dim);
    const auto c = ug_read_pt(rPts, static_cast<std::size_t>(i2), dim);
    const double ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
    const double vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2];
    return {uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx};
}

}  // namespace

UndoGreenResult undo_green(const Mesh& rFine) {
    if (!rFine.HasPointData(kRefineEntityName))
        throw std::invalid_argument(
            std::string(kUgPrefix) + "single-mesh undo_green requires '" + kRefineEntityName +
            "' point_data; run refine(..., record_hierarchy=True, record_levels=True) first");
    const NDArray& entity_arr = rFine.PointData(kRefineEntityName);
    const std::size_t npts = rFine.NumPoints();
    if (detail::rows(entity_arr) != npts || (npts != 0 && entity_arr.Size() / npts != 4))
        throw std::invalid_argument(std::string(kUgPrefix) + "'" + kRefineEntityName +
                                    "' must have shape (num_points, 4)");
    std::vector<std::array<std::int64_t, 4>> entity_keys(npts);
    for (std::size_t p = 0; p < npts; ++p)
        for (std::size_t k = 0; k < 4; ++k)
            entity_keys[p][k] = detail::read_int(entity_arr, p * 4 + k);

    const std::vector<std::int64_t> fine_bases = detail::block_bases(rFine);
    const std::size_t nblocks = rFine.NumCellBlocks();
    const std::size_t total_fine = static_cast<std::size_t>(fine_bases.back());
    UgColumn fine_cid(rFine, fine_bases, kRefineCellIdName);
    UgColumn fine_pid(rFine, fine_bases, kRefineParentIdName);
    UgColumn fine_lvl(rFine, fine_bases, kRefineLevelName);

    const NDArray& pts = rFine.Points();
    const std::size_t pdim = rFine.PointDim();

    enum class UgSingleRole : std::uint8_t { Keep = 0, Anchor = 1, Suppressed = 2 };
    std::vector<UgSingleRole> role(total_fine, UgSingleRole::Keep);
    std::vector<std::vector<std::int64_t>> anchor_corners(total_fine);
    std::vector<std::int64_t> anchor_pid(total_fine, -1);
    std::vector<std::size_t> group_anchor(total_fine, 0);
    std::int64_t num_groups = 0;
    std::int64_t num_removed = 0;

    std::size_t bi = 0;
    for (const auto cb : rFine.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        const std::size_t red_count = ug_red_child_count(cb.Type());
        const std::size_t nc = cb.IsRagged() ? 0 : cb.NodesPerCell();
        if (red_count <= 2 || nc == 0 || ncells == 0) {
            ++bi;
            continue;
        }
        // Group cells in this block by parent_id (where parent_id != cell_id).
        std::unordered_map<std::int64_t, std::vector<std::size_t>> groups;
        for (std::size_t r = 0; r < ncells; ++r) {
            const std::size_t g = static_cast<std::size_t>(fine_bases[bi]) + r;
            const std::int64_t pid = fine_pid.At(g);
            const std::int64_t cid = fine_cid.At(g);
            if (pid != cid)
                groups[pid].push_back(r);
        }
        const NDArray& sconn = cb.Conn();
        for (const auto& [pid, rows] : groups) {
            if (rows.size() <= 1 || rows.size() >= red_count)
                continue;
            // Collect all distinct vertices used by this green sibling group.
            std::vector<std::int64_t> verts;
            for (std::size_t r : rows) {
                for (std::size_t k = 0; k < nc; ++k)
                    verts.push_back(detail::read_int(sconn, r * nc + k));
            }
            std::sort(verts.begin(), verts.end());
            verts.erase(std::unique(verts.begin(), verts.end()), verts.end());

            // A vertex is a mid-entity node of THIS split iff its refine:entity has >= 2 valid
            // parent corners AND every valid corner in refine:entity is present in `verts`.
            auto in_verts = [&](std::int64_t v) {
                return std::binary_search(verts.begin(), verts.end(), v);
            };
            std::vector<std::int64_t> corners;
            for (std::int64_t v : verts) {
                bool is_mid = false;
                if (v >= 0 && static_cast<std::size_t>(v) < npts) {
                    const auto& ek = entity_keys[static_cast<std::size_t>(v)];
                    int valid = 0;
                    bool all_in = true;
                    for (int k = 0; k < 4; ++k) {
                        if (ek[k] >= 0) {
                            ++valid;
                            if (!in_verts(ek[k]))
                                all_in = false;
                        }
                    }
                    if (valid >= 2 && all_in)
                        is_mid = true;
                }
                if (!is_mid)
                    corners.push_back(v);
            }
            if (corners.size() != nc)
                continue;

            // Orient corners to match child 0's orientation.
            const std::size_t r0 = rows[0];
            if (cb.Type() == "triangle") {
                const auto n0 = ug_tri_normal(pts, pdim, detail::read_int(sconn, r0 * nc + 0),
                                              detail::read_int(sconn, r0 * nc + 1),
                                              detail::read_int(sconn, r0 * nc + 2));
                const auto np = ug_tri_normal(pts, pdim, corners[0], corners[1], corners[2]);
                if (n0[0] * np[0] + n0[1] * np[1] + n0[2] * np[2] < 0.0)
                    std::swap(corners[1], corners[2]);
            } else if (cb.Type() == "tetra") {
                std::vector<std::int64_t> c0 = {detail::read_int(sconn, r0 * nc + 0),
                                                detail::read_int(sconn, r0 * nc + 1),
                                                detail::read_int(sconn, r0 * nc + 2),
                                                detail::read_int(sconn, r0 * nc + 3)};
                const double v0 = ug_tet_signed_vol(pts, pdim, c0);
                const double vp = ug_tet_signed_vol(pts, pdim, corners);
                if (v0 * vp < 0.0)
                    std::swap(corners[2], corners[3]);
            } else if (cb.Type() == "quad") {
                // Order 4 corners cyclically around their centroid using child 0's normal.
                const auto n0 = ug_tri_normal(pts, pdim, detail::read_int(sconn, r0 * nc + 0),
                                              detail::read_int(sconn, r0 * nc + 1),
                                              detail::read_int(sconn, r0 * nc + 2));
                std::array<double, 3> ctr{0.0, 0.0, 0.0};
                for (std::int64_t v : corners) {
                    const auto p = ug_read_pt(pts, static_cast<std::size_t>(v), pdim);
                    ctr[0] += p[0] * 0.25;
                    ctr[1] += p[1] * 0.25;
                    ctr[2] += p[2] * 0.25;
                }
                const auto p0 = ug_read_pt(pts, static_cast<std::size_t>(corners[0]), pdim);
                std::array<double, 3> e1{p0[0] - ctr[0], p0[1] - ctr[1], p0[2] - ctr[2]};
                std::array<double, 3> e2{n0[1] * e1[2] - n0[2] * e1[1],
                                         n0[2] * e1[0] - n0[0] * e1[2],
                                         n0[0] * e1[1] - n0[1] * e1[0]};
                std::sort(corners.begin(), corners.end(), [&](std::int64_t a, std::int64_t b) {
                    const auto pa = ug_read_pt(pts, static_cast<std::size_t>(a), pdim);
                    const auto pb = ug_read_pt(pts, static_cast<std::size_t>(b), pdim);
                    const double ax = (pa[0] - ctr[0]) * e1[0] + (pa[1] - ctr[1]) * e1[1] +
                                      (pa[2] - ctr[2]) * e1[2];
                    const double ay = (pa[0] - ctr[0]) * e2[0] + (pa[1] - ctr[1]) * e2[1] +
                                      (pa[2] - ctr[2]) * e2[2];
                    const double bx = (pb[0] - ctr[0]) * e1[0] + (pb[1] - ctr[1]) * e1[1] +
                                      (pb[2] - ctr[2]) * e1[2];
                    const double by = (pb[0] - ctr[0]) * e2[0] + (pb[1] - ctr[1]) * e2[1] +
                                      (pb[2] - ctr[2]) * e2[2];
                    return std::atan2(ay, ax) < std::atan2(by, bx);
                });
            }

            const std::size_t anchor_g = static_cast<std::size_t>(fine_bases[bi]) + rows[0];
            role[anchor_g] = UgSingleRole::Anchor;
            anchor_corners[anchor_g] = std::move(corners);
            anchor_pid[anchor_g] = pid;
            group_anchor[anchor_g] = anchor_g;
            for (std::size_t k = 1; k < rows.size(); ++k) {
                const std::size_t sg = static_cast<std::size_t>(fine_bases[bi]) + rows[k];
                role[sg] = UgSingleRole::Suppressed;
                group_anchor[sg] = anchor_g;
            }
            ++num_groups;
            num_removed += static_cast<std::int64_t>(rows.size());
        }
        ++bi;
    }

    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rFine.Points()));
    for (const std::string& name : rFine.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rFine.PointData(name)));
    for (const std::string& name : rFine.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rFine.FieldData(name)));

    std::vector<NDArray> cell_maps(nblocks);
    std::vector<std::size_t> out_ncells(nblocks, 0);
    std::vector<std::int64_t> anchor_out_row(total_fine, -1);

    bi = 0;
    for (const auto cb : rFine.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        std::size_t kept = 0;
        for (std::size_t r = 0; r < ncells; ++r) {
            const std::size_t g = static_cast<std::size_t>(fine_bases[bi]) + r;
            if (role[g] != UgSingleRole::Suppressed)
                ++kept;
        }
        out_ncells[bi] = kept;
        NDArray cmap = NDArray::Uninit(DType::Int64, {ncells});
        auto* cm = cmap.As<std::int64_t>();
        if (cb.IsRagged()) {
            out.AddCellBlock(cb.Type(), detail::data_owned_copy(cb.Conn()));
            for (std::size_t r = 0; r < ncells; ++r)
                cm[r] = static_cast<std::int64_t>(r);
        } else {
            const std::size_t nc = cb.NodesPerCell();
            const NDArray& sconn = cb.Conn();
            NDArray conn = NDArray::Uninit(DType::Int64, {kept, nc});
            auto* dst = conn.As<std::int64_t>();
            std::size_t out_row = 0;
            for (std::size_t r = 0; r < ncells; ++r) {
                const std::size_t g = static_cast<std::size_t>(fine_bases[bi]) + r;
                if (role[g] == UgSingleRole::Suppressed)
                    continue;
                if (role[g] == UgSingleRole::Keep) {
                    for (std::size_t k = 0; k < nc; ++k)
                        dst[out_row * nc + k] = detail::read_int(sconn, r * nc + k);
                } else {
                    for (std::size_t k = 0; k < nc; ++k)
                        dst[out_row * nc + k] = anchor_corners[g][k];
                    anchor_out_row[g] = static_cast<std::int64_t>(out_row);
                }
                cm[r] = static_cast<std::int64_t>(out_row);
                ++out_row;
            }
            for (std::size_t r = 0; r < ncells; ++r) {
                const std::size_t g = static_cast<std::size_t>(fine_bases[bi]) + r;
                if (role[g] == UgSingleRole::Suppressed)
                    cm[r] = anchor_out_row[group_anchor[g]];
            }
            out.AddCellBlock(cb.Type(), std::move(conn));
        }
        cell_maps[bi] = std::move(cmap);
        ++bi;
    }

    // Carry cell_data (restoring cell_id = parent_id on undone anchors).
    for (const std::string& name : rFine.CellDataNames()) {
        if (rFine.CellDataNumBlocks(name) != nblocks)
            continue;
        std::vector<NDArray> out_blocks(nblocks);
        bi = 0;
        for (const auto cb : rFine.CellRange()) {
            const NDArray& fsrc = rFine.CellData(name, bi);
            std::vector<std::size_t> shape = fsrc.Shape();
            shape[0] = out_ncells[bi];
            NDArray dst = NDArray::Uninit(fsrc.Dtype(), shape);
            const std::size_t row_bytes =
                fsrc.Nbytes() / std::max<std::size_t>(1, detail::rows(fsrc));
            std::byte* p_dst = dst.Data();
            const std::byte* p_fsrc = fsrc.Data();
            std::size_t out_row = 0;
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::size_t g = static_cast<std::size_t>(fine_bases[bi]) + r;
                if (role[g] == UgSingleRole::Suppressed)
                    continue;
                if (role[g] == UgSingleRole::Anchor &&
                    (name == kRefineCellIdName || name == kRefineParentIdName) &&
                    fsrc.Dtype() == DType::Int64) {
                    dst.As<std::int64_t>()[out_row] = anchor_pid[g];
                } else if (row_bytes > 0) {
                    std::memcpy(p_dst + out_row * row_bytes, p_fsrc + r * row_bytes, row_bytes);
                }
                ++out_row;
            }
            out_blocks[bi] = std::move(dst);
            ++bi;
        }
        out.AddCellData(name, std::move(out_blocks));
    }

    UndoGreenResult res;
    res.mMesh = std::move(out);
    res.mCellMaps = std::move(cell_maps);
    res.mNumGroupsUndone = num_groups;
    res.mNumCellsRemoved = num_removed;

    detail::RegionRemap rmap;
    rmap.mCellMapKind = detail::CellMapKind::Direct;
    rmap.pCellMaps = &res.mCellMaps;
    rmap.mOpName = "undo_green";
    detail::remap_regions(rFine, res.mMesh, rmap);

    return res;
}

}  // namespace meshioplusplus
