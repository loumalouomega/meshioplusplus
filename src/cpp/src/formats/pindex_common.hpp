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
 * @file formats/pindex_common.hpp
 * @brief Machinery shared by the ParaView index formats `.pvtu`, `.pvtp` and
 * `.pvd`: resolving a piece path, reading and merging pieces, dropping ghost
 * cells, and summarizing pieces' metadata.
 *
 * A **format-private** header (the `formats/gid_common.hpp` precedent): it sits
 * beside the `.cpp` files because nothing outside these formats needs it and no
 * installed header may name it. Its functions are in a *named* namespace
 * (`pidx`), not an anonymous one, because the amalgamation concatenates every
 * `.cpp` under `src/cpp/src` into one translation unit, and the amalgamator
 * inlines a quoted include exactly once.
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/pvtp.hpp"
#include "meshioplusplus/formats/pvtu.hpp"
#include "meshioplusplus/formats/vtm.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/operations/crop.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace pidx {

namespace fs = std::filesystem;

/// VTK's ghost vocabulary (`vtkDataSetAttributes`) and the name `partition`
/// tags its halo layers with.
inline constexpr const char* kGhostName = "vtkGhostType";
inline constexpr const char* kPartitionGhost = "partition:ghost";

/// The value of an XML attribute, escaped for a double-quoted context.
inline std::string escape_attr(const std::string& rValue) {
    std::string out;
    out.reserve(rValue.size());
    for (const char c : rValue) {
        switch (c) {
            case '&':
                out += "&amp;";
                break;
            case '<':
                out += "&lt;";
                break;
            case '>':
                out += "&gt;";
                break;
            case '"':
                out += "&quot;";
                break;
            default:
                out += c;
        }
    }
    return out;
}

/**
 * Where a piece named by @p rFile lives: resolved against the index's own
 * directory; an absolute path is used as written. There is deliberately no
 * search for the file elsewhere -- an absolute path from another machine that
 * does not exist here is an error naming the attribute and the index, because a
 * fallback could read the wrong file.
 */
inline fs::path resolve_path(const std::string& rIndexPath, const std::string& rFile,
                             const std::string& rFormat, const char* pAttr) {
    if (rFile.empty())
        throw ReadError("meshio++: " + rFormat + ": an entry is missing its '" + pAttr +
                        "' attribute: " + rIndexPath);
    const fs::path file(rFile);
    const fs::path base = fs::path(rIndexPath).parent_path();
    const fs::path path = file.is_absolute() || base.empty() ? file : base / file;
    std::error_code ec;
    if (!fs::exists(path, ec))
        throw ReadError("meshio++: " + rFormat + ": piece file '" + rFile + "' (" + pAttr +
                        "=) named by " + rIndexPath + " does not exist (looked for " +
                        path.string() +
                        "); a relative path is resolved against the index's own directory and an "
                        "absolute path is not searched for elsewhere");
    return path;
}

/// `<stem>_0000<ext>` ... zero-padded to `max(4, digits(count - 1))`, the width
/// `{step}` patterns use.
inline std::vector<std::string> piece_names(const std::string& rStem, std::size_t Count,
                                            const std::string& rExt) {
    std::size_t width = 4;
    if (Count > 1)
        width = std::max<std::size_t>(width, std::to_string(Count - 1).size());
    std::vector<std::string> names;
    names.reserve(Count);
    for (std::size_t i = 0; i < Count; ++i) {
        std::string idx = std::to_string(i);
        names.push_back(rStem + "_" + std::string(width - idx.size(), '0') + idx + rExt);
    }
    return names;
}

/// The options handed to a piece's own reader: the narrowing options apply, but
/// the index's step and piece selectors were consumed by the index itself.
inline ReadOptions child_options(const ReadOptions& rOpts) {
    ReadOptions child = rOpts;
    child.mTimeStep = 0;
    child.mPiece = 0;
    child.mPieceSet = false;
    return child;
}

/// A mesh with no points and no cells: what an index naming no piece reads as.
inline Mesh empty_mesh() {
    Mesh empty;
    empty.AssignPoints(NDArray::Uninit(DType::Float64, {0, 3}));
    return empty;
}

/**
 * Read one piece with the reader for its extension. @p Wide widens the set to
 * what a `.pvd` may name (`.vtm`, `.pvtu`, `.pvtp` as well); a parallel index
 * names serial files only, so the only nesting is a `.pvd` over them and no
 * cycle is expressible.
 */
inline Mesh read_child(const fs::path& rPath, const ReadOptions& rOpts,
                       const PvtuReadOptions& rGhost, bool Wide) {
    const std::string ext = rPath.extension().string();
    const ReadOptions child = child_options(rOpts);
    if (ext == ".vtu")
        return read_vtu(rPath.string(), child);
    if (ext == ".vtp")
        return read_vtp(rPath.string(), child);
    if (Wide) {
        if (ext == ".vtm")
            return read_vtm(rPath.string(), child);
        if (ext == ".pvtu")
            return read_pvtu(rPath.string(), child, rGhost);
        if (ext == ".pvtp")
            return read_pvtp(rPath.string(), child, rGhost);
    }
    throw ReadError("meshio++: unsupported piece '" + rPath.string() + "': only " +
                    (Wide ? ".vtu/.vtp/.vtm/.pvtu/.pvtp" : ".vtu/.vtp") + " pieces are read");
}

inline MeshMetadata read_child_metadata(const fs::path& rPath, const ReadOptions& rOpts,
                                        bool Wide) {
    const std::string ext = rPath.extension().string();
    const ReadOptions child = child_options(rOpts);
    if (ext == ".vtp")
        return read_vtp_metadata(rPath.string(), child);
    if (ext == ".vtu")
        return read_vtu_metadata(rPath.string(), child);
    if (Wide) {
        if (ext == ".vtm")
            return read_vtm_metadata(rPath.string(), child);
        if (ext == ".pvtu")
            return read_pvtu_metadata(rPath.string(), child);
        if (ext == ".pvtp")
            return read_pvtp_metadata(rPath.string(), child);
    }
    throw ReadError("meshio++: unsupported piece '" + rPath.string() + "': only " +
                    (Wide ? ".vtu/.vtp/.vtm/.pvtu/.pvtp" : ".vtu/.vtp") + " pieces are read");
}

/// Whether every entry of @p rArray is zero (any integer dtype).
inline bool all_zero(const NDArray& rArray) {
    for (const std::int64_t v : detail::vtu_to_int64(rArray))
        if (v != 0)
            return false;
    return true;
}

/**
 * Drop ghost cells (any `vtkGhostType` bit set) and the points only they used,
 * then the now meaningless ghost arrays. A piece with no `vtkGhostType` array is
 * returned unchanged; one with the array but no flagged cell keeps its points
 * and loses only the array. `partition:ghost` goes too, when it is all zero.
 */
inline Mesh drop_ghosts(Mesh Piece) {
    if (!Piece.HasCellData(kGhostName))
        return Piece;
    bool any_flagged = false;
    for (std::size_t b = 0; b < Piece.CellDataNumBlocks(kGhostName) && !any_flagged; ++b)
        any_flagged = !all_zero(Piece.CellData(kGhostName, b));

    Mesh kept;
    if (any_flagged)
        kept = crop_predicate(Piece, kGhostName, RefineCompare::Equal, 0.0).mMesh;
    else
        kept = std::move(Piece);

    bool drop_layers = kept.HasCellData(kPartitionGhost);
    for (std::size_t b = 0; drop_layers && b < kept.CellDataNumBlocks(kPartitionGhost); ++b)
        drop_layers = all_zero(kept.CellData(kPartitionGhost, b));
    return detail::clone_mesh(
        kept, [drop_layers](DataLocation Where, const std::string& rName, std::string&) {
            if (Where == DataLocation::Field)
                return true;
            if (rName == kGhostName)
                return false;
            return !(drop_layers && Where == DataLocation::Cell && rName == kPartitionGhost);
        });
}

/**
 * Merge @p rPieces without welding and add one `Cell` region per piece, named
 * `rNames[i]`. A single piece is returned as read: wrapping it in one all-cells
 * region would hide the regions it already has (a `.pvd` step that is one
 * `.pvtu`).
 */
inline Mesh merge_pieces(std::vector<Mesh> Pieces, const std::vector<std::string>& rNames) {
    if (Pieces.size() == 1)
        return std::move(Pieces.front());

    std::vector<const Mesh*> ptrs;
    ptrs.reserve(Pieces.size());
    for (const Mesh& m : Pieces)
        ptrs.push_back(&m);

    MergeOptions mopts;
    mopts.weld = false;
    mopts.source_tag = false;
    mopts.data_policy = MergeDataPolicy::Fill;
    MergeResult result = merge(ptrs, mopts);

    for (std::size_t i = 0; i < Pieces.size(); ++i)
        result.mMesh.AddRegion(Region(rNames[i], RegionKind::Cell, std::move(result.mCellMaps[i])));
    return std::move(result.mMesh);
}

/// The sum of the pieces' own metadata, in the same first-seen, same-type
/// consolidation order `merge_pieces` would produce.
inline MeshMetadata aggregate_metadata(const std::vector<MeshMetadata>& rParts,
                                       const std::string& rFormat) {
    MeshMetadata meta;
    meta.mFormat = rFormat;
    std::vector<CellBlockInfo> blocks;
    std::vector<std::string> point_names, cell_names, field_names;
    auto merge_names = [](std::vector<std::string>& rInto, const std::vector<std::string>& rFrom) {
        for (const std::string& n : rFrom)
            if (std::find(rInto.begin(), rInto.end(), n) == rInto.end())
                rInto.push_back(n);
    };
    for (const MeshMetadata& pm : rParts) {
        meta.mNumPoints += pm.mNumPoints;
        meta.mPointDim = std::max(meta.mPointDim, pm.mPointDim);
        meta.mFellBackToFullRead = meta.mFellBackToFullRead || pm.mFellBackToFullRead;
        merge_names(point_names, pm.mPointDataNames);
        merge_names(cell_names, pm.mCellDataNames);
        merge_names(field_names, pm.mFieldDataNames);
        for (const CellBlockInfo& cb : pm.mCellBlocks) {
            auto it = std::find_if(blocks.begin(), blocks.end(),
                                   [&](const CellBlockInfo& rB) { return rB.mType == cb.mType; });
            if (it == blocks.end()) {
                blocks.push_back(cb);
            } else {
                it->mNumCells += cb.mNumCells;
                it->mRagged = it->mRagged || cb.mRagged;
                if (it->mNodesPerCell != cb.mNodesPerCell)
                    it->mNodesPerCell = 0;
            }
        }
    }
    meta.mCellBlocks = std::move(blocks);
    std::sort(point_names.begin(), point_names.end());
    std::sort(cell_names.begin(), cell_names.end());
    std::sort(field_names.begin(), field_names.end());
    meta.mPointDataNames = std::move(point_names);
    meta.mCellDataNames = std::move(cell_names);
    meta.mFieldDataNames = std::move(field_names);
    return meta;
}

}  // namespace pidx
}  // namespace meshioplusplus
