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
#include <cstdint>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "slot_runs.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

/// Per-block bookkeeping for the two parallel passes: where this block's cells
/// land in the compact numbering, and where its faces land in the slot space.
struct FmBlockDesc {
    std::size_t mBlock = 0;      ///< index into the mesh's cell blocks
    std::size_t mFirstCell = 0;  ///< first compact cell id
    std::size_t mNumCells = 0;
};

/// Does this block bound a volume? Mirrors `cell_rings`' own acceptance rule --
/// deliberately, so a block can never be counted here and then declined there.
bool fm_block_has_volume(const Mesh::CellView& rBlock) {
    if (rBlock.IsPolyhedron())
        return true;
    if (rBlock.IsRagged())
        return false;  // 1-level (polygon): a surface cell, no enclosed volume
    return skin_supported(cell_type_from_name(rBlock.Type()));
}

}  // namespace

GlobalFaces build_global_faces(const Mesh& rMesh) {
    return build_global_faces(rMesh, {});
}

GlobalFaces build_global_faces(const Mesh& rMesh, const std::vector<std::size_t>& rBlocks) {
    GlobalFaces out;

    // Empty means "every block"; that is what the one-argument form passes, so
    // there is exactly one implementation rather than two that could drift.
    std::vector<bool> wanted;
    if (!rBlocks.empty()) {
        wanted.assign(rMesh.NumCellBlocks(), false);
        for (std::size_t b : rBlocks)
            if (b < wanted.size())
                wanted[b] = true;
    }

    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::vector<std::int64_t> bases = block_bases(rMesh);

    // ---- Phase 0: the compact cell space --------------------------------
    // Serial and cheap: one pass over the blocks, no per-cell work.
    std::vector<FmBlockDesc> descs;
    std::size_t n_cells = 0;
    for (std::size_t k = 0; k < rMesh.NumCellBlocks(); ++k) {
        const auto cb = rMesh.Cells(k);
        if ((!wanted.empty() && !wanted[k]) || !fm_block_has_volume(cb)) {
            out.mNonCellBlocks.push_back(k);
            continue;
        }
        FmBlockDesc d;
        d.mBlock = k;
        d.mFirstCell = n_cells;
        d.mNumCells = cb.NumCells();
        n_cells += d.mNumCells;
        descs.push_back(d);
    }

    out.mCellToGlobal.resize(n_cells);
    for (const FmBlockDesc& d : descs) {
        const std::int64_t base = bases[d.mBlock];
        for (std::size_t i = 0; i < d.mNumCells; ++i)
            out.mCellToGlobal[d.mFirstCell + i] = base + static_cast<std::int64_t>(i);
    }

    // ---- Phase 1: face counts and ring sizes ----------------------------
    // `cell_rings` is a gather, so running it twice is cheap; `orient_rings`
    // (the BFS) runs only in phase 2. Counting first is what lets phase 2 write
    // straight into pre-sized buffers instead of holding a per-cell vector pair
    // for the whole mesh -- on a ten-million-cell case that is gigabytes.
    std::vector<std::int64_t> n_cell_faces(n_cells, 0);
    for (const FmBlockDesc& d : descs) {
        const auto cb = rMesh.Cells(d.mBlock);
        parallel_for(d.mNumCells, [&](std::size_t i) {
            static thread_local CellRings rings;
            static thread_local std::vector<Vec3> coords;
            if (cell_rings(cb, i, points, dim, rings, coords))
                n_cell_faces[d.mFirstCell + i] = static_cast<std::int64_t>(rings.NumFaces());
        });
    }

    out.mCellFaceStart.resize(n_cells + 1);
    out.mCellFaceStart[n_cells] = parallel_exclusive_scan(
        n_cell_faces.data(), n_cells, out.mCellFaceStart.data(), std::int64_t{0});
    const std::size_t n_slots = static_cast<std::size_t>(out.mCellFaceStart[n_cells]);

    // Ring sizes per slot, then their prefix, so phase 2's writes are disjoint.
    std::vector<std::int64_t> slot_size(n_slots, 0);
    for (const FmBlockDesc& d : descs) {
        const auto cb = rMesh.Cells(d.mBlock);
        parallel_for(d.mNumCells, [&](std::size_t i) {
            static thread_local CellRings rings;
            static thread_local std::vector<Vec3> coords;
            if (!cell_rings(cb, i, points, dim, rings, coords))
                return;
            const std::size_t base = static_cast<std::size_t>(out.mCellFaceStart[d.mFirstCell + i]);
            for (std::size_t f = 0; f < rings.NumFaces(); ++f)
                slot_size[base + f] = static_cast<std::int64_t>(rings.FaceSize(f));
        });
    }

    std::vector<std::int64_t> slot_start(n_slots + 1);
    slot_start[n_slots] =
        parallel_exclusive_scan(slot_size.data(), n_slots, slot_start.data(), std::int64_t{0});

    // ---- Phase 2: build every cell's outward-wound rings ----------------
    std::vector<std::int64_t> slot_nodes(static_cast<std::size_t>(slot_start[n_slots]));
    std::vector<std::uint8_t> cell_flipped(n_cells, 0);
    std::vector<std::uint8_t> cell_unorientable(n_cells, 0);

    for (const FmBlockDesc& d : descs) {
        const auto cb = rMesh.Cells(d.mBlock);
        parallel_for(d.mNumCells, [&](std::size_t i) {
            static thread_local CellRings rings;
            static thread_local std::vector<Vec3> coords;
            const std::size_t cell = d.mFirstCell + i;
            if (!cell_rings(cb, i, points, dim, rings, coords))
                return;
            // Repair unconditionally: `cell_faces.hpp`'s rows are outward on the
            // REFERENCE element, so an inverted hexahedron yields six inward
            // normals and every one of its faces would be written misoriented.
            const RingOrientation ro = orient_rings(rings, coords.data());
            if (ro == RingOrientation::Repaired)
                cell_flipped[cell] = 1;
            else if (ro == RingOrientation::Unorientable)
                cell_unorientable[cell] = 1;

            const std::size_t base = static_cast<std::size_t>(out.mCellFaceStart[cell]);
            for (std::size_t f = 0; f < rings.NumFaces(); ++f) {
                const std::uint32_t* ring = rings.Face(f);
                const std::size_t n = rings.FaceSize(f);
                std::int64_t* dst = slot_nodes.data() + slot_start[base + f];
                for (std::size_t k = 0; k < n; ++k)
                    dst[k] = rings.mNodes[ring[k]];
            }
        });
    }

    for (std::size_t c = 0; c < n_cells; ++c) {
        out.mNumFlipped += cell_flipped[c];
        out.mNumUnorientable += cell_unorientable[c];
    }

    // ---- Phase 3: first-seen dedup ---------------------------------------
    // Slot order is ascending (compact cell, local face), and face ids are
    // handed out in the order a serial sweep over the slots first meets each
    // corner set -- so they, and therefore `owner < neighbour`, which the
    // OpenFOAM writer validates rather than assumes, do not depend on thread
    // count or hash order. The sort-based table (slot_runs.hpp) recovers that
    // sweep from a parallel sort: a face's owner is the cell of its first slot,
    // its neighbour the cell of its second, and every further slot is one more
    // non-manifold use.
    std::vector<std::int64_t> cell_of_slot(n_slots);
    parallel_for(n_cells, [&](std::size_t c) {
        for (auto s = out.mCellFaceStart[c]; s < out.mCellFaceStart[c + 1]; ++s)
            cell_of_slot[static_cast<std::size_t>(s)] = static_cast<std::int64_t>(c);
    });
    std::vector<FacetKey> keys(n_slots);
    parallel_for(n_slots, [&](std::size_t s) {
        keys[s] =
            FacetKey(slot_nodes.data() + slot_start[s], static_cast<std::size_t>(slot_size[s]));
    });
    const SlotRuns runs = group_facet_slots(
        keys, [](const FacetKey& rK) -> const FacetKey& { return rK; }, rMesh.NumPoints());
    std::vector<FacetKey>().swap(keys);
    const FirstSeen seen = number_first_seen(runs, n_slots);
    const std::size_t n_faces = seen.NumIds();

    out.mOwner.resize(n_faces);
    out.mNeighbour.resize(n_faces);
    out.mFaceStart.resize(n_faces + 1);
    std::vector<std::int64_t> face_size(n_faces);
    parallel_for(n_faces, [&](std::size_t fid) {
        const std::size_t r = static_cast<std::size_t>(seen.mRunOfId[fid]);
        const std::uint64_t head = runs.Head(r);
        out.mOwner[fid] = cell_of_slot[head];
        out.mNeighbour[fid] = runs.Size(r) > 1 ? cell_of_slot[runs.Begin(r)[1]] : -1;
        face_size[fid] = slot_size[head];
    });
    out.mFaceStart[n_faces] =
        parallel_exclusive_scan(face_size.data(), n_faces, out.mFaceStart.data(), std::int64_t{0});
    out.mFaceNodes.resize(static_cast<std::size_t>(out.mFaceStart[n_faces]));
    parallel_for(n_faces, [&](std::size_t fid) {
        const std::uint64_t head = runs.Head(static_cast<std::size_t>(seen.mRunOfId[fid]));
        std::copy_n(slot_nodes.data() + slot_start[head], face_size[fid],
                    out.mFaceNodes.data() + out.mFaceStart[fid]);
    });

    out.mCellFaces.resize(n_slots);
    parallel_for(n_slots, [&](std::size_t s) {
        const std::int64_t fid = seen.mIdOfSlot[s];
        // Positive: stored as this cell wound it (its first use); negative:
        // reversed from stored.
        const bool first = runs.Head(static_cast<std::size_t>(seen.mRunOfId[fid])) == s;
        out.mCellFaces[s] = first ? fid + 1 : -(fid + 1);
    });
    out.mNumNonManifold = parallel_reduce(
        runs.NumRuns(), 4096, std::int64_t{0},
        [&](std::size_t b, std::size_t e) {
            std::int64_t extra = 0;
            for (std::size_t r = b; r < e; ++r)
                extra += runs.Size(r) > 2 ? static_cast<std::int64_t>(runs.Size(r) - 2) : 0;
            return extra;
        },
        [](std::int64_t acc, std::int64_t part) { return acc + part; });

    return out;
}

FaceLookup::FaceLookup(const GlobalFaces& rFaces) {
    mSorted.resize(rFaces.NumFaces());
    parallel_for(rFaces.NumFaces(), [&](std::size_t f) {
        mSorted[f] = {FacetKey(rFaces.Face(f), rFaces.FaceSize(f)), static_cast<std::int64_t>(f)};
    });
    // Face corner sets are distinct, and the id breaks any tie a malformed
    // list could hold: a total order, so the table is the same on every build.
    const FacetKeyLess less;
    parallel_sort(mSorted.begin(), mSorted.end(), [&](const auto& rA, const auto& rB) {
        if (less(rA.first, rB.first))
            return true;
        if (less(rB.first, rA.first))
            return false;
        return rA.second < rB.second;
    });
}

std::int64_t FaceLookup::Find(const std::int64_t* pIds, std::size_t N) const {
    const FacetKey key(pIds, N);
    const FacetKeyLess less;
    const auto it = std::lower_bound(mSorted.begin(), mSorted.end(), key,
                                     [&](const std::pair<FacetKey, std::int64_t>& rE,
                                         const FacetKey& rK) { return less(rE.first, rK); });
    return it == mSorted.end() || !(it->first == key) ? -1 : it->second;
}

}  // namespace detail
}  // namespace meshioplusplus
