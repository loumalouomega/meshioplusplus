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
// Green-element undo: lookup and substitution against a caller-supplied coarse
// mesh -- see operations/undo_green.hpp for the contract and the reasoning
// for why no subdivision-table inversion is needed here.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/operations/undo_green.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/refine_hierarchy.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/refine.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kUgPrefix = "meshio++: undo_green: ";

// Read a required Int64 scalar cell_data array (one value per cell, covering
// every block) into a flat global-cell-order vector. Throws by name on any
// mismatch -- undo_green's fine-mesh preconditions are hard requirements,
// unlike detail::refine_read_hierarchy's warn-and-fall-back contract.
std::vector<std::int64_t> ug_read_required(const Mesh& rMesh,
                                           const std::vector<std::int64_t>& rBases,
                                           const char* pName) {
    if (!rMesh.HasCellData(pName))
        throw std::invalid_argument(
            std::string(kUgPrefix) + "the fine mesh has no '" + pName +
            "' cell_data; run refine(..., record_hierarchy=True, record_levels=True) first");
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (rMesh.CellDataNumBlocks(pName) != nblocks)
        throw std::invalid_argument(std::string(kUgPrefix) + "'" + pName +
                                    "' does not cover every cell block");
    std::vector<std::int64_t> out(static_cast<std::size_t>(detail::total_cells(rBases)));
    std::size_t bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t base = static_cast<std::size_t>(rBases[bi]);
        const std::size_t ncells = cb.NumCells();
        const NDArray& a = rMesh.CellData(pName, bi);
        if (detail::rows(a) != ncells || (ncells != 0 && a.Size() / ncells != 1))
            throw std::invalid_argument(std::string(kUgPrefix) + "'" + pName + "' block " +
                                        std::to_string(bi) + " is not one scalar value per cell");
        for (std::size_t c = 0; c < ncells; ++c)
            out[base + c] = detail::read_int(a, c);
        ++bi;
    }
    return out;
}

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
enum class UgRole { Keep, AnchorSubstitute, Suppressed };

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

    const std::vector<std::int64_t> fine_id =
        ug_read_required(rFine, fine_bases, kRefineCellIdName);
    const std::vector<std::int64_t> fine_parent =
        ug_read_required(rFine, fine_bases, kRefineParentIdName);
    const std::vector<std::int64_t> fine_level =
        ug_read_required(rFine, fine_bases, kRefineLevelName);

    // --- resolve the coarse mesh's id space: its own recorded ids if valid,
    // else the implicit global-block-major ids -- refine_attach_hierarchy's
    // own "start fresh" fallback. --------------------------------------------
    std::vector<std::int64_t> coarse_ids;
    std::int64_t coarse_id_base = 0;
    const detail::RefineHierarchyState state =
        detail::refine_read_hierarchy(rCoarse, coarse_bases, coarse_ids, coarse_id_base);
    if (state != detail::RefineHierarchyState::Valid) {
        const auto total_coarse = static_cast<std::size_t>(detail::total_cells(coarse_bases));
        coarse_ids.resize(total_coarse);
        for (std::size_t i = 0; i < total_coarse; ++i)
            coarse_ids[i] = static_cast<std::int64_t>(i);
    }
    std::unordered_map<std::int64_t, std::int64_t> id_to_coarse_row;
    id_to_coarse_row.reserve(coarse_ids.size());
    for (std::size_t i = 0; i < coarse_ids.size(); ++i)
        id_to_coarse_row[coarse_ids[i]] = static_cast<std::int64_t>(i);

    const std::vector<std::int64_t> coarse_level =
        ug_read_optional(rCoarse, coarse_bases, kRefineLevelName);
    auto coarse_level_at = [&](std::int64_t row) -> std::int64_t {
        return coarse_level.empty() ? 0 : coarse_level[static_cast<std::size_t>(row)];
    };

    // --- group fine cells by parent_id --------------------------------------
    std::unordered_map<std::int64_t, std::vector<std::int64_t>> groups_by_parent;
    for (std::size_t g = 0; g < total_fine; ++g)
        groups_by_parent[fine_parent[g]].push_back(static_cast<std::int64_t>(g));

    // --- classify every global fine cell's role; for green groups, resolve
    // the substitution source once ------------------------------------------
    std::vector<UgRole> role(total_fine, UgRole::Keep);
    std::vector<std::int64_t> group_id_of(total_fine, -1);
    std::vector<std::int64_t> group_coarse_row;  // indexed by green group id
    std::vector<std::int64_t> group_size;        // indexed by green group id
    std::int64_t next_group_id = 0;

    for (auto& entry : groups_by_parent) {
        const std::int64_t parent_id = entry.first;
        std::vector<std::int64_t>& members = entry.second;
        std::sort(members.begin(), members.end());

        if (members.size() == 1) {
            const std::int64_t g = members.front();
            if (fine_id[static_cast<std::size_t>(g)] != parent_id)
                throw std::invalid_argument(
                    std::string(kUgPrefix) +
                    "malformed hierarchy: a singleton sibling group's refine:cell_id does not "
                    "equal its refine:parent_id (cell " +
                    std::to_string(g) + ")");
            continue;  // untouched: role stays Keep
        }

        const auto it = id_to_coarse_row.find(parent_id);
        if (it == id_to_coarse_row.end())
            throw std::invalid_argument(
                std::string(kUgPrefix) + "refine:parent_id " + std::to_string(parent_id) +
                " does not resolve in the coarse mesh's id space -- these two meshes are not "
                "the coarse/fine pair of one refine() call");
        const std::int64_t coarse_row = it->second;
        const std::int64_t clevel = coarse_level_at(coarse_row);

        const std::int64_t flevel = fine_level[static_cast<std::size_t>(members.front())];
        for (std::int64_t g : members)
            if (fine_level[static_cast<std::size_t>(g)] != flevel)
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
        const std::int64_t gid = next_group_id++;
        group_coarse_row.push_back(coarse_row);
        group_size.push_back(static_cast<std::int64_t>(members.size()));
        role[static_cast<std::size_t>(members.front())] = UgRole::AnchorSubstitute;
        group_id_of[static_cast<std::size_t>(members.front())] = gid;
        for (std::size_t k = 1; k < members.size(); ++k) {
            role[static_cast<std::size_t>(members[static_cast<std::ptrdiff_t>(k)])] =
                UgRole::Suppressed;
            group_id_of[static_cast<std::size_t>(members[static_cast<std::ptrdiff_t>(k)])] = gid;
        }
    }

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
    std::vector<std::int64_t> group_output_global(static_cast<std::size_t>(next_group_id), -1);
    std::vector<std::size_t> group_fine_block(static_cast<std::size_t>(next_group_id), 0);

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
                    const std::int64_t gid = group_id_of[g];
                    cm[r] = group_output_global[static_cast<std::size_t>(gid)];
                    continue;
                }
                if (role[g] == UgRole::Keep) {
                    const NDArray& src = cb.Conn();
                    for (std::size_t k = 0; k < npc; ++k)
                        dst[out_row * npc + k] = detail::read_int(src, r * npc + k);
                } else {  // AnchorSubstitute
                    const std::int64_t gid = group_id_of[g];
                    const std::int64_t coarse_row = group_coarse_row[static_cast<std::size_t>(gid)];
                    const auto blk_row = detail::global_to_block_row(coarse_bases, coarse_row);
                    const auto ccb = rCoarse.Cells(blk_row.first);
                    const NDArray& csrc = ccb.Conn();
                    for (std::size_t k = 0; k < npc; ++k)
                        dst[out_row * npc + k] = detail::read_int(
                            csrc, static_cast<std::size_t>(blk_row.second) * npc + k);
                    group_output_global[static_cast<std::size_t>(gid)] =
                        out_block_base + static_cast<std::int64_t>(out_row);
                    group_fine_block[static_cast<std::size_t>(gid)] = bi;
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
        if (next_group_id > 0) {
            if (!rCoarse.HasCellData(name) ||
                rCoarse.CellDataNumBlocks(name) != rCoarse.NumCellBlocks()) {
                ok = false;
            } else {
                for (std::int64_t gid = 0; ok && gid < next_group_id; ++gid) {
                    const auto ug = static_cast<std::size_t>(gid);
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
                    const std::int64_t gid = group_id_of[g];
                    const auto ug = static_cast<std::size_t>(gid);
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
    res.mNumGroupsUndone = next_group_id;
    std::int64_t removed = 0;
    for (std::int64_t sz : group_size)
        removed += sz - 1;
    res.mNumCellsRemoved = removed;

    detail::RegionRemap rmap;
    rmap.mCellMapKind = detail::CellMapKind::Direct;
    rmap.pCellMaps = &res.mCellMaps;
    rmap.mDropSideRegions = true;
    rmap.mOpName = "undo_green";
    detail::remap_regions(rFine, res.mMesh, rmap);

    return res;
}

}  // namespace meshioplusplus
