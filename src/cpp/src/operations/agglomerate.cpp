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
// Polyhedral coarsening: greedy seed-and-grow over the shared-face dual
// (detail::build_global_faces), one merged polyhedron per group -- see
// operations/agglomerate.hpp for the contract.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kAggPrefix = "meshio++: agglomerate: ";

/// A staged output block: either an unchanged copy of a pass-through input
/// block (any of its three storage shapes), or the one merged polyhedron
/// block. Mirrors subdivide.cpp's own staging shape; not shared with it,
/// since that one is file-private there too.
struct AggOutBlock {
    std::string mType;
    std::size_t mNodesPerCell = 0;  // rectangular pass-through only
    std::vector<std::int64_t> mConn;
    std::vector<std::vector<std::int64_t>> mPolygonRows;
    std::vector<std::vector<std::vector<std::int64_t>>> mPolyhedronCells;
    bool mIsRagged = false;
    bool mIsPolyhedron = false;
};

AggOutBlock agg_stage_passthrough(const Mesh::CellView& rBlock) {
    AggOutBlock out;
    out.mType = std::string(rBlock.Type());
    if (rBlock.IsPolyhedron()) {
        out.mIsRagged = true;
        out.mIsPolyhedron = true;
        out.mPolyhedronCells.resize(rBlock.NumCells());
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c) {
            for (std::size_t f = 0; f < rBlock.NumFaces(c); ++f) {
                const auto face = rBlock.Face(c, f);
                out.mPolyhedronCells[c].emplace_back(face.first, face.first + face.second);
            }
        }
    } else if (rBlock.IsRagged()) {
        out.mIsRagged = true;
        out.mPolygonRows.resize(rBlock.NumCells());
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c) {
            const std::int64_t* row = rBlock.Row(c);
            out.mPolygonRows[c].assign(row, row + rBlock.RowSize(c));
        }
    } else {
        const NDArray& conn = rBlock.Conn();
        const std::size_t npc = rBlock.NodesPerCell();
        out.mNodesPerCell = npc;
        out.mConn.resize(rBlock.NumCells() * npc);
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c)
            for (std::size_t k = 0; k < npc; ++k)
                out.mConn[c * npc + k] = detail::read_int(conn, c * npc + k);
    }
    return out;
}

void agg_emit_block(Mesh& rOut, AggOutBlock& rBlock) {
    if (rBlock.mIsPolyhedron) {
        rOut.AddPolyhedronBlock(rBlock.mType, std::move(rBlock.mPolyhedronCells));
    } else if (rBlock.mIsRagged) {
        rOut.AddPolygonBlock(rBlock.mType, std::move(rBlock.mPolygonRows));
    } else {
        const std::size_t npc = rBlock.mNodesPerCell;
        const std::size_t ncells = npc == 0 ? 0 : rBlock.mConn.size() / npc;
        NDArray conn = NDArray::Uninit(DType::Int64, {ncells, npc});
        std::int64_t* dst = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < rBlock.mConn.size(); ++i)
            dst[i] = rBlock.mConn[i];
        rOut.AddCellBlock(rBlock.mType, std::move(conn));
    }
}

/// Area of global face `f`, gathering its corner coordinates into a local
/// buffer first -- GlobalFaces::Face() returns GLOBAL node ids, while
/// polygon_area's ring-taking overload expects local indices, so the plain
/// (no-ring) overload over a freshly gathered buffer is the fit.
double agg_face_area(const detail::GlobalFaces& rFaces, std::size_t f, const NDArray& rPoints,
                     std::size_t PointDim) {
    const std::size_t n = rFaces.FaceSize(f);
    const std::int64_t* ring = rFaces.Face(f);
    std::vector<detail::Vec3> coords(n);
    for (std::size_t k = 0; k < n; ++k) {
        const auto pid = static_cast<std::size_t>(ring[k]);
        for (std::size_t d = 0; d < 3; ++d)
            coords[k][d] = d < PointDim ? detail::read_double(rPoints, pid * PointDim + d) : 0.0;
    }
    return detail::polygon_area(coords.data(), n);
}

/// Frontier entry ordered by DESCENDING accumulated shared-face area, ties
/// broken by ASCENDING compact cell id -- storing the negated area keeps a
/// plain ascending std::set a max-by-area, min-by-id priority structure.
struct FrontierKey {
    double mNegArea;
    std::int64_t mId;
    bool operator<(const FrontierKey& rOther) const {
        if (mNegArea != rOther.mNegArea)
            return mNegArea < rOther.mNegArea;
        return mId < rOther.mId;
    }
};

}  // namespace

AgglomerateResult agglomerate(const Mesh& rMesh, const AgglomerateOptions& rOptions) {
    if (rOptions.mTargetGroupSize == 0)
        throw std::invalid_argument(std::string(kAggPrefix) + "mTargetGroupSize must be >= 1");

    const detail::GlobalFaces gf = detail::build_global_faces(rMesh);
    if (gf.mNumNonManifold > 0)
        throw std::invalid_argument(
            std::string(kAggPrefix) + "mesh contains " + std::to_string(gf.mNumNonManifold) +
            " face(s) shared by three or more cells (non-manifold); refusing rather than "
            "guessing a boundary classification");

    const std::size_t n_compact = gf.NumCells();
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();

    // --- greedy seed-and-grow over the face dual --------------------------
    std::vector<std::int64_t> group_of(n_compact, -1);
    std::vector<std::vector<std::int64_t>> groups;

    for (std::size_t seed = 0; seed < n_compact; ++seed) {
        if (group_of[seed] != -1)
            continue;
        const auto gid = static_cast<std::int64_t>(groups.size());
        std::vector<std::int64_t> members{static_cast<std::int64_t>(seed)};
        group_of[seed] = gid;

        std::unordered_map<std::int64_t, double> pending;
        std::set<FrontierKey> frontier;

        auto push_neighbours = [&](std::int64_t c) {
            const std::size_t nf = gf.NumCellFaces(static_cast<std::size_t>(c));
            const std::int64_t* row = gf.CellFaces(static_cast<std::size_t>(c));
            for (std::size_t k = 0; k < nf; ++k) {
                const std::int64_t sid = row[k];
                const auto f = static_cast<std::size_t>((sid > 0 ? sid : -sid) - 1);
                const std::int64_t owner = gf.mOwner[f];
                const std::int64_t neigh = gf.mNeighbour[f];
                const std::int64_t other = (owner == c) ? neigh : owner;
                if (other < 0)
                    continue;  // mesh boundary
                if (group_of[static_cast<std::size_t>(other)] != -1)
                    continue;  // already claimed (by this group or would be a bug otherwise)
                const double a = agg_face_area(gf, f, points, pdim);
                auto it = pending.find(other);
                if (it != pending.end()) {
                    frontier.erase(FrontierKey{-it->second, other});
                    it->second += a;
                } else {
                    it = pending.emplace(other, a).first;
                }
                frontier.insert(FrontierKey{-it->second, other});
            }
        };

        push_neighbours(static_cast<std::int64_t>(seed));

        while (members.size() < rOptions.mTargetGroupSize && !frontier.empty()) {
            const auto fit = frontier.begin();
            const std::int64_t c = fit->mId;
            frontier.erase(fit);
            pending.erase(c);
            if (group_of[static_cast<std::size_t>(c)] != -1)
                continue;  // defensive; unreachable given the push-time check above
            group_of[static_cast<std::size_t>(c)] = gid;
            members.push_back(c);
            push_neighbours(c);
        }

        std::sort(members.begin(), members.end());
        groups.push_back(std::move(members));
    }

    // --- emit: one polyhedron cell per group, external faces only ---------
    std::vector<std::vector<std::vector<std::int64_t>>> merged_cells(groups.size());
    for (std::size_t g = 0; g < groups.size(); ++g) {
        auto& faces_out = merged_cells[g];
        for (std::int64_t c : groups[g]) {
            const std::size_t nf = gf.NumCellFaces(static_cast<std::size_t>(c));
            const std::int64_t* row = gf.CellFaces(static_cast<std::size_t>(c));
            for (std::size_t k = 0; k < nf; ++k) {
                const std::int64_t sid = row[k];
                const auto f = static_cast<std::size_t>((sid > 0 ? sid : -sid) - 1);
                const std::int64_t owner = gf.mOwner[f];
                const std::int64_t neigh = gf.mNeighbour[f];
                const std::int64_t other = (owner == c) ? neigh : owner;
                if (other >= 0 && group_of[static_cast<std::size_t>(other)] ==
                                      group_of[static_cast<std::size_t>(c)])
                    continue;  // internal to the group: dropped from both sides
                const std::size_t n = gf.FaceSize(f);
                const std::int64_t* ring = gf.Face(f);
                std::vector<std::int64_t> face_nodes(ring, ring + n);
                if (sid < 0)
                    std::reverse(face_nodes.begin(), face_nodes.end());
                faces_out.push_back(std::move(face_nodes));
            }
        }
    }

    // --- which original blocks carry volume, and the compact<->global bridge
    std::vector<bool> is_volume(rMesh.NumCellBlocks(), true);
    for (std::size_t b : gf.mNonCellBlocks)
        if (b < is_volume.size())
            is_volume[b] = false;

    const std::vector<std::int64_t> in_bases = detail::block_bases(rMesh);
    std::unordered_map<std::int64_t, std::int64_t> global_to_compact;
    global_to_compact.reserve(n_compact);
    for (std::size_t c = 0; c < n_compact; ++c)
        global_to_compact[gf.mCellToGlobal[c]] = static_cast<std::int64_t>(c);

    // --- build the output mesh's blocks, tracking where each original block
    // landed (pass-through) or whether it fed the merged block -----------
    Mesh out;
    out.AssignPoints(detail::data_owned_copy(points));

    const std::size_t nblocks_in = rMesh.NumCellBlocks();
    std::vector<std::int64_t> pass_out_block(nblocks_in, -1);
    std::int64_t merged_out_block = -1;
    bool merged_emitted = false;
    std::size_t bi = 0;
    std::int64_t out_block_idx = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (is_volume[bi]) {
            if (!merged_emitted) {
                AggOutBlock mb;
                mb.mType = "polyhedron";
                mb.mIsRagged = true;
                mb.mIsPolyhedron = true;
                mb.mPolyhedronCells = std::move(merged_cells);
                if (!mb.mPolyhedronCells.empty()) {
                    agg_emit_block(out, mb);
                    merged_out_block = out_block_idx++;
                }
                merged_emitted = true;
            }
        } else {
            AggOutBlock pb = agg_stage_passthrough(cb);
            agg_emit_block(out, pb);
            pass_out_block[bi] = out_block_idx++;
        }
        ++bi;
    }

    const std::vector<std::int64_t> out_bases = detail::block_bases(out);

    // --- the flat cell map --------------------------------------------------
    NDArray cell_map =
        NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(detail::total_cells(in_bases))});
    std::int64_t* cm = cell_map.As<std::int64_t>();
    bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t base_in = in_bases[bi];
        if (is_volume[bi]) {
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::int64_t global_in = base_in + static_cast<std::int64_t>(r);
                const std::int64_t compact = global_to_compact.at(global_in);
                cm[global_in] = merged_out_block < 0
                                    ? -1
                                    : out_bases[static_cast<std::size_t>(merged_out_block)] +
                                          group_of[static_cast<std::size_t>(compact)];
            }
        } else {
            const std::int64_t ob = pass_out_block[bi];
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::int64_t global_in = base_in + static_cast<std::int64_t>(r);
                cm[global_in] =
                    out_bases[static_cast<std::size_t>(ob)] + static_cast<std::int64_t>(r);
            }
        }
        ++bi;
    }

    // --- point_data / field_data: unchanged, no new or pruned points -------
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    // --- cell_data: pass-through blocks copied verbatim; the merged block's
    // row for group g is its first (ascending compact id) member's row -----
    for (const std::string& name : rMesh.CellDataNames()) {
        const std::size_t ndata = rMesh.CellDataNumBlocks(name);
        if (ndata != nblocks_in) {
            log::warn(
                "{}cell_data '{}' does not have one array per input block; dropped rather "
                "than guessed at",
                kAggPrefix, name);
            continue;
        }
        std::vector<NDArray> out_blocks(static_cast<std::size_t>(out_block_idx));
        bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            (void)cb;
            if (!is_volume[bi] && pass_out_block[bi] >= 0)
                out_blocks[static_cast<std::size_t>(pass_out_block[bi])] =
                    detail::data_owned_copy(rMesh.CellData(name, bi));
            ++bi;
        }
        if (merged_out_block >= 0) {
            // Resolve each group's first member's (block, row) and gather.
            const NDArray& first_src = rMesh.CellData(name, 0);
            std::vector<std::size_t> shape = first_src.Shape();
            shape[0] = groups.size();
            NDArray dst = NDArray::Uninit(first_src.Dtype(), std::move(shape));
            const std::size_t row_bytes =
                first_src.Nbytes() / std::max<std::size_t>(1, detail::rows(first_src));
            std::byte* p_dst = dst.Data();
            for (std::size_t g = 0; g < groups.size(); ++g) {
                const std::int64_t compact_first = groups[g].front();
                const std::int64_t global =
                    gf.mCellToGlobal[static_cast<std::size_t>(compact_first)];
                const auto [blk, row] = detail::global_to_block_row(in_bases, global);
                const NDArray& src = rMesh.CellData(name, blk);
                const std::byte* p_src = src.Data();
                std::memcpy(p_dst + g * row_bytes,
                            p_src + static_cast<std::size_t>(row) * row_bytes, row_bytes);
            }
            out_blocks[static_cast<std::size_t>(merged_out_block)] = std::move(dst);
        }
        out.AddCellData(name, std::move(out_blocks));
    }

    AgglomerateResult res;
    res.mMesh = std::move(out);
    res.mCellMap = std::move(cell_map);

    detail::RegionRemap rmap;
    rmap.mCellMapKind = detail::CellMapKind::Global;
    rmap.pGlobalCellMap = &res.mCellMap;
    rmap.mOpName = "agglomerate";
    detail::remap_regions(rMesh, res.mMesh, rmap);

    return res;
}

}  // namespace meshioplusplus
