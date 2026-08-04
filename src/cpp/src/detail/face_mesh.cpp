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
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

/// Per-block bookkeeping for the two parallel passes: where this block's cells
/// land in the compact numbering, and where its faces land in the slot space.
struct FmBlockDesc {
    std::size_t mBlock = 0;        ///< index into the mesh's cell blocks
    std::size_t mFirstCell = 0;    ///< first compact cell id
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
    out.mCellFaceStart[0] = 0;
    for (std::size_t c = 0; c < n_cells; ++c)
        out.mCellFaceStart[c + 1] = out.mCellFaceStart[c] + n_cell_faces[c];
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
            const std::size_t base =
                static_cast<std::size_t>(out.mCellFaceStart[d.mFirstCell + i]);
            for (std::size_t f = 0; f < rings.NumFaces(); ++f)
                slot_size[base + f] = static_cast<std::int64_t>(rings.FaceSize(f));
        });
    }

    std::vector<std::int64_t> slot_start(n_slots + 1);
    slot_start[0] = 0;
    for (std::size_t s = 0; s < n_slots; ++s)
        slot_start[s + 1] = slot_start[s] + slot_size[s];

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

    // ---- Phase 3: SERIAL first-seen dedup -------------------------------
    // Slot order is ascending (compact cell, local face), so face ids -- and
    // therefore `owner < neighbour`, which the OpenFOAM writer validates rather
    // than assumes -- do not depend on thread count or hash order.
    out.mCellFaces.resize(n_slots);
    out.mFaceStart.push_back(0);

    std::unordered_map<FacetKey, std::int64_t, FacetKeyHash> first_use;
    first_use.reserve(n_slots);

    for (std::size_t c = 0; c < n_cells; ++c) {
        const std::size_t lo = static_cast<std::size_t>(out.mCellFaceStart[c]);
        const std::size_t hi = static_cast<std::size_t>(out.mCellFaceStart[c + 1]);
        for (std::size_t s = lo; s < hi; ++s) {
            const std::int64_t* ring = slot_nodes.data() + slot_start[s];
            const std::size_t n = static_cast<std::size_t>(slot_size[s]);
            const FacetKey key(ring, n);
            auto it = first_use.find(key);
            if (it == first_use.end()) {
                const std::int64_t fid = static_cast<std::int64_t>(out.NumFaces());
                out.mFaceNodes.insert(out.mFaceNodes.end(), ring, ring + n);
                out.mFaceStart.push_back(static_cast<std::int64_t>(out.mFaceNodes.size()));
                out.mOwner.push_back(static_cast<std::int64_t>(c));
                out.mNeighbour.push_back(-1);
                first_use.emplace(key, fid);
                out.mCellFaces[s] = fid + 1;  // positive: stored as this cell wound it
            } else {
                const std::int64_t fid = it->second;
                if (out.mNeighbour[static_cast<std::size_t>(fid)] < 0)
                    out.mNeighbour[static_cast<std::size_t>(fid)] = static_cast<std::int64_t>(c);
                else
                    ++out.mNumNonManifold;
                out.mCellFaces[s] = -(fid + 1);  // negative: reversed from stored
            }
        }
    }

    return out;
}

FaceLookup::FaceLookup(const GlobalFaces& rFaces) {
    mMap.reserve(rFaces.NumFaces());
    for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
        mMap.emplace(FacetKey(rFaces.Face(f), rFaces.FaceSize(f)), static_cast<std::int64_t>(f));
}

std::int64_t FaceLookup::Find(const std::int64_t* pIds, std::size_t N) const {
    const auto it = mMap.find(FacetKey(pIds, N));
    return it == mMap.end() ? -1 : it->second;
}

}  // namespace detail
}  // namespace meshioplusplus
