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
// Surface/boundary extraction, following the facet-hashing algorithm of Kratos
// Multiphysics' SkinDetectionProcess: a facet (a face of a 3D cell, or an edge
// of a 2D cell) whose sorted corner-node key occurs exactly once across all
// cells of the chosen dimension is a boundary facet; facets seen twice are
// interior and cancel out. This single translation unit implements both the
// general `extract_surface` and the volume-only, linearize-able `extract_skin`
// (declared in skin.hpp) so the two never drift.

// System includes
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/skin.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// --- shared facet abstraction -----------------------------------------------

// Which kind of boundary facet we extract.
enum class SurfaceMode { Face, Edge };

// One facet of a cell, corners-first: unifies CellFaceDef (3D) and CellEdgeDef
// (2D) into a single shape so the two-phase extractor is dimension-agnostic.
struct SurfaceFacetDef {
    CellType mType;                      // output facet type
    std::uint8_t mNumCorners;            // leading entries of mNodes
    std::uint8_t mNumNodes;              // full node count of the facet
    std::array<std::uint8_t, 9> mNodes;  // local node indices, corners first
};

// The facets of a cell type for the given mode (empty if unsupported).
std::vector<SurfaceFacetDef> surface_facets_for(CellType type, SurfaceMode mode) {
    std::vector<SurfaceFacetDef> out;
    if (mode == SurfaceMode::Face) {
        for (const detail::CellFaceDef& fd : detail::cell_faces(type))
            out.push_back({fd.mFaceType, fd.mNumCorners, fd.mNumNodes, fd.mNodes});
    } else {
        for (const detail::CellEdgeDef& ed : detail::cell_edges(type)) {
            std::array<std::uint8_t, 9> nodes = {};
            nodes[0] = ed.mNodes[0];
            nodes[1] = ed.mNodes[1];
            nodes[2] = ed.mNodes[2];
            out.push_back({ed.mEdgeType, ed.mNumCorners, ed.mNumNodes, nodes});
        }
    }
    return out;
}

bool surface_mode_supported(CellType type, SurfaceMode mode) {
    return mode == SurfaceMode::Face ? detail::skin_supported(type)
                                     : detail::surface_edge_supported(type);
}

// Effective topological dimension for mode selection (polyhedron counts as 3).
int surface_effective_dim(const Mesh::CellView& rBlock) {
    if (rBlock.IsPolyhedron())
        return 3;
    return cell_type_dimension(cell_type_from_name(rBlock.Type()));
}

// A supported 3D block is a non-ragged volume block with a known face table.
bool surface_skin_block_supported(CellType type, bool is_ragged) {
    return cell_type_dimension(type) == 3 && !is_ragged && detail::skin_supported(type);
}

// --- facet key + hash -------------------------------------------------------
//
// detail::FacetKey rather than the fixed array<int64_t,4> this used before
// v9.16.0: a polyhedron's face can have any number of corners. It keeps a
// 4-entry inline fast path, so a triangle or quad still costs no allocation and
// the rectangular path is unchanged. Sharing ONE key type across both is what
// lets a hexahedron and a polyhedron that meet on a face cancel each other out
// -- with two key types they would each report that face as boundary.

using SurfaceFacetKey = detail::FacetKey;
using SurfaceFacetKeyHash = detail::FacetKeyHash;

SurfaceFacetKey surface_facet_key(const NDArray& rConn, std::size_t rowOffset,
                                  const SurfaceFacetDef& rFacet) {
    std::array<std::int64_t, 4> ids{};
    for (std::uint8_t k = 0; k < rFacet.mNumCorners; ++k)
        ids[k] = detail::read_int(rConn, rowOffset + rFacet.mNodes[k]);
    return SurfaceFacetKey(ids.data(), rFacet.mNumCorners);
}

// --- canonical output buckets -----------------------------------------------

// Face mode: triangle, triangle6, quad, quad8, quad9 (matches legacy skin).
// Edge mode: line, line3.
std::size_t surface_num_out_types(SurfaceMode mode) {
    return mode == SurfaceMode::Face ? 5 : 2;
}

std::size_t surface_out_type_index(SurfaceMode mode, CellType type) {
    if (mode == SurfaceMode::Face) {
        switch (type) {
            case CellType::Triangle:
                return 0;
            case CellType::Triangle6:
                return 1;
            case CellType::Quad:
                return 2;
            case CellType::Quad8:
                return 3;
            default:  // Quad9
                return 4;
        }
    }
    return type == CellType::Line ? 0 : 1;  // Line / Line3
}

const char* surface_out_type_name(SurfaceMode mode, std::size_t index) {
    static const char* face_names[5] = {"triangle", "triangle6", "quad", "quad8", "quad9"};
    static const char* edge_names[2] = {"line", "line3"};
    return mode == SurfaceMode::Face ? face_names[index] : edge_names[index];
}

std::size_t surface_out_type_nodes(SurfaceMode mode, std::size_t index) {
    static const std::size_t face_counts[5] = {3, 6, 4, 8, 9};
    static const std::size_t edge_counts[2] = {2, 3};
    return mode == SurfaceMode::Face ? face_counts[index] : edge_counts[index];
}

// --- per-block enumeration descriptor ---------------------------------------

struct SurfaceBlockDesc {
    const NDArray* mpConn = nullptr;  // rectangular blocks only
    std::size_t mNpc = 0;
    std::size_t mNumCells = 0;
    std::vector<SurfaceFacetDef> mFacets;  // rectangular blocks only
    std::size_t mFacetsPerCell = 0;        // rectangular blocks only
    std::size_t mFirstFacet = 0;           // offset into the flat record buffer
    std::int64_t mGlobalCellBase = 0;      // input-cell index of this block's cell 0
    // Polyhedron blocks: faces per cell vary, so the flat record offsets are
    // tabulated instead of computed as cell * mFacetsPerCell.
    bool mPolyhedron = false;
    std::size_t mBlock = 0;               // index into the mesh's cell blocks
    std::vector<std::size_t> mFaceStart;  // mNumCells + 1 record offsets
};

// One facet occurrence recorded in phase 1 (keys built in parallel).
struct SurfaceFacetRecord {
    SurfaceFacetKey mKey;
    std::uint32_t mDesc;
    std::uint32_t mCell;
    std::uint32_t mSlot;
    std::int64_t mParent;
};

}  // namespace

namespace detail {

// The shared core: extract boundary facets of the chosen mode.
//  - forceFaceMode true  => volume-only (extract_skin), honoring `linearize`.
//  - forceFaceMode false => auto (extract_surface): faces if any volume cell
//    exists, else edges.
Mesh surface_extract(const Mesh& rMesh, bool forceFaceMode, bool linearize, bool recordParentIds,
                     const char* pOpName) {
    // --- select the mode ---
    SurfaceMode mode = SurfaceMode::Face;
    if (!forceFaceMode) {
        int max_dim = -1;
        for (const auto cb : rMesh.CellRange())
            max_dim = std::max(max_dim, surface_effective_dim(cb));
        if (max_dim == 3)
            mode = SurfaceMode::Face;
        else if (max_dim == 2)
            mode = SurfaceMode::Edge;
        else
            throw std::invalid_argument(std::string(pOpName) +
                                        ": mesh contains no supported 2D or 3D cell block");
    }

    // --- pre-scan: warn per unsupported same-dimension block, require one ---
    // Face mode mirrors the legacy skin pre-scan (polyhedron counts as volume);
    // the wording ("volume"/"3D volume") is kept byte-identical to the old
    // extract_skin messages so its callers/tests are unaffected.
    const char* noun = mode == SurfaceMode::Face ? "volume" : "surface";
    const int scan_dim = mode == SurfaceMode::Face ? 3 : 2;
    bool any_supported = false;
    for (const auto cb : rMesh.CellRange()) {
        const CellType ct = cell_type_from_name(cb.Type());
        const bool same_dim = mode == SurfaceMode::Face
                                  ? (cell_type_dimension(ct) == 3 || cb.IsPolyhedron())
                                  : (cell_type_dimension(ct) == 2 && !cb.IsPolyhedron());
        if (!same_dim)
            continue;
        // A polyhedron block carries its own faces, so it is supported by
        // construction -- there is no table for it to be missing from.
        if (cb.IsPolyhedron() || (!cb.IsRagged() && surface_mode_supported(ct, mode))) {
            any_supported = true;
        } else {
            log::warn("{}: {} cell block '{}' is not supported; skipping it.", pOpName, noun,
                      cb.Type());
        }
    }
    if (!any_supported)
        throw std::invalid_argument(std::string(pOpName) + ": mesh contains no supported " +
                                    std::to_string(scan_dim) + "D " + noun + " cell block");

    // --- build per-block descriptors + a flat record buffer ---
    std::vector<SurfaceBlockDesc> descs;
    std::size_t total_facets = 0;
    std::int64_t global_cell_base = 0;
    std::size_t block_idx = 0;
    for (const auto cb : rMesh.CellRange()) {
        const CellType ct = cell_type_from_name(cb.Type());
        const std::int64_t base = global_cell_base;
        global_cell_base += static_cast<std::int64_t>(cb.NumCells());
        ++block_idx;
        const bool same_dim = mode == SurfaceMode::Face
                                  ? (cell_type_dimension(ct) == 3 || cb.IsPolyhedron())
                                  : (cell_type_dimension(ct) == 2 && !cb.IsPolyhedron());
        if (!same_dim)
            continue;
        if (cb.IsPolyhedron()) {
            // Variable faces per cell: tabulate the record offsets instead of
            // computing them from a constant facets-per-cell.
            SurfaceBlockDesc d;
            d.mPolyhedron = true;
            d.mBlock = block_idx - 1;
            d.mNumCells = cb.NumCells();
            d.mFirstFacet = total_facets;
            d.mGlobalCellBase = base;
            d.mFaceStart.reserve(d.mNumCells + 1);
            std::size_t at = 0;
            d.mFaceStart.push_back(0);
            for (std::size_t c = 0; c < d.mNumCells; ++c) {
                at += cb.NumFaces(c);
                d.mFaceStart.push_back(at);
            }
            total_facets += at;
            descs.push_back(std::move(d));
            continue;
        }
        if (cb.IsRagged() || !surface_mode_supported(ct, mode))
            continue;
        SurfaceBlockDesc d;
        d.mpConn = &cb.Conn();
        d.mNpc = cb.NodesPerCell();
        d.mNumCells = cb.NumCells();
        d.mFacets = surface_facets_for(ct, mode);
        d.mFacetsPerCell = d.mFacets.size();
        d.mFirstFacet = total_facets;
        d.mGlobalCellBase = base;
        total_facets += d.mNumCells * d.mFacetsPerCell;
        descs.push_back(std::move(d));
    }

    // --- phase 1: build facet keys into disjoint slots (parallel-safe) ---
    std::vector<SurfaceFacetRecord> recs(total_facets);
    for (std::uint32_t di = 0; di < descs.size(); ++di) {
        const SurfaceBlockDesc& d = descs[di];
        if (d.mPolyhedron) {
            const auto cb = rMesh.Cells(d.mBlock);
            parallel_for(d.mNumCells, [&, di](std::size_t c) {
                const std::size_t nf = cb.NumFaces(c);
                for (std::size_t f = 0; f < nf; ++f) {
                    const auto face = cb.Face(c, f);
                    SurfaceFacetRecord& r = recs[d.mFirstFacet + d.mFaceStart[c] + f];
                    r.mKey = SurfaceFacetKey(face.first, face.second);
                    r.mDesc = di;
                    r.mCell = static_cast<std::uint32_t>(c);
                    r.mSlot = static_cast<std::uint32_t>(f);
                    r.mParent = d.mGlobalCellBase + static_cast<std::int64_t>(c);
                }
            });
            continue;
        }
        const NDArray& conn = *d.mpConn;
        const std::size_t npc = d.mNpc;
        const std::size_t fpc = d.mFacetsPerCell;
        const std::size_t n = d.mNumCells * fpc;
        parallel_for(n, [&, di](std::size_t j) {
            const std::size_t cell = j / fpc;
            const std::size_t slot = j % fpc;
            SurfaceFacetRecord& r = recs[d.mFirstFacet + j];
            r.mKey = surface_facet_key(conn, cell * npc, d.mFacets[slot]);
            r.mDesc = di;
            r.mCell = static_cast<std::uint32_t>(cell);
            r.mSlot = static_cast<std::uint32_t>(slot);
            r.mParent = d.mGlobalCellBase + static_cast<std::int64_t>(cell);
        });
    }

    // --- phase 2, pass A: count key occurrences (serial → deterministic) ---
    std::unordered_map<SurfaceFacetKey, std::uint32_t, SurfaceFacetKeyHash> face_count;
    face_count.reserve(total_facets * 2);
    for (const SurfaceFacetRecord& r : recs)
        ++face_count[r.mKey];

    // --- phase 2, pass B: emit boundary facets (count == 1) in stored order ---
    const std::size_t num_out = surface_num_out_types(mode);
    std::vector<std::vector<std::int64_t>> out_conn(num_out);
    std::vector<std::vector<std::int64_t>> out_parent(num_out);
    // A polyhedral face may be an n-gon, which no fixed-width bucket can hold;
    // those go to a ragged `polygon` block emitted alongside the fixed ones.
    std::vector<std::vector<std::int64_t>> poly_rows;
    std::vector<std::int64_t> poly_parent;
    for (const SurfaceFacetRecord& r : recs) {
        if (face_count[r.mKey] != 1)
            continue;
        const SurfaceBlockDesc& d = descs[r.mDesc];
        if (d.mPolyhedron) {
            const auto cb = rMesh.Cells(d.mBlock);
            const auto face =
                cb.Face(static_cast<std::size_t>(r.mCell), static_cast<std::size_t>(r.mSlot));
            if (face.second == 3 || face.second == 4) {
                const std::size_t bucket = surface_out_type_index(
                    mode, face.second == 3 ? CellType::Triangle : CellType::Quad);
                std::vector<std::int64_t>& dst = out_conn[bucket];
                for (std::size_t k = 0; k < face.second; ++k)
                    dst.push_back(face.first[k]);
                if (recordParentIds)
                    out_parent[bucket].push_back(r.mParent);
            } else {
                poly_rows.emplace_back(face.first, face.first + face.second);
                if (recordParentIds)
                    poly_parent.push_back(r.mParent);
            }
            continue;
        }
        const NDArray& conn = *d.mpConn;
        const std::size_t row_offset = static_cast<std::size_t>(r.mCell) * d.mNpc;
        const SurfaceFacetDef& facet = d.mFacets[r.mSlot];
        const std::uint8_t n_out = linearize ? facet.mNumCorners : facet.mNumNodes;
        const CellType out_type =
            linearize ? (facet.mNumCorners == 3 ? CellType::Triangle : CellType::Quad)
                      : facet.mType;
        const std::size_t bucket = surface_out_type_index(mode, out_type);
        std::vector<std::int64_t>& dst = out_conn[bucket];
        for (std::uint8_t k = 0; k < n_out; ++k)
            dst.push_back(detail::read_int(conn, row_offset + facet.mNodes[k]));
        if (recordParentIds)
            out_parent[bucket].push_back(r.mParent);
    }

    // --- compaction: keep only referenced points, ascending original index ---
    const std::size_t num_points = rMesh.NumPoints();
    std::vector<char> used(num_points, 0);
    for (const std::vector<std::int64_t>& conn : out_conn)
        for (std::int64_t id : conn)
            used[static_cast<std::size_t>(id)] = 1;
    for (const std::vector<std::int64_t>& row : poly_rows)
        for (std::int64_t id : row)
            used[static_cast<std::size_t>(id)] = 1;
    std::vector<std::int64_t> remap(num_points, -1);
    std::size_t num_used = 0;
    for (std::size_t i = 0; i < num_points; ++i)
        if (used[i])
            remap[i] = static_cast<std::int64_t>(num_used++);

    Mesh surface;

    // Points: byte-row subset preserving the source dtype.
    {
        const NDArray& points = rMesh.Points();
        const std::size_t dim = detail::cols(points);
        const std::size_t row_bytes = dim * dtype_size(points.Dtype());
        NDArray new_points = NDArray::Uninit(points.Dtype(), {num_used, dim});
        const std::byte* src = points.Data();
        std::byte* dst = new_points.Data();
        std::size_t w = 0;
        for (std::size_t i = 0; i < num_points; ++i)
            if (used[i])
                std::memcpy(dst + (w++) * row_bytes, src + i * row_bytes, row_bytes);
        surface.AssignPoints(std::move(new_points));
    }

    // Cell blocks in canonical order, remapped connectivity; optional parent ids.
    std::vector<NDArray> parent_blocks;
    for (std::size_t t = 0; t < num_out; ++t) {
        const std::vector<std::int64_t>& conn = out_conn[t];
        if (conn.empty())
            continue;
        const std::size_t npc = surface_out_type_nodes(mode, t);
        NDArray block = NDArray::Uninit(DType::Int64, {conn.size() / npc, npc});
        std::int64_t* dst = block.As<std::int64_t>();
        for (std::size_t i = 0; i < conn.size(); ++i)
            dst[i] = remap[static_cast<std::size_t>(conn[i])];
        surface.AddCellBlock(surface_out_type_name(mode, t), std::move(block));
        if (recordParentIds) {
            const std::vector<std::int64_t>& par = out_parent[t];
            NDArray pblock = NDArray::Uninit(DType::Int64, {par.size(), 1});
            std::memcpy(pblock.Data(), par.data(), par.size() * sizeof(std::int64_t));
            parent_blocks.push_back(std::move(pblock));
        }
    }
    if (!poly_rows.empty()) {
        for (std::vector<std::int64_t>& row : poly_rows)
            for (std::int64_t& id : row)
                id = remap[static_cast<std::size_t>(id)];
        surface.AddPolygonBlock("polygon", std::move(poly_rows));
        if (recordParentIds) {
            NDArray pblock = NDArray::Uninit(DType::Int64, {poly_parent.size(), 1});
            std::memcpy(pblock.Data(), poly_parent.data(),
                        poly_parent.size() * sizeof(std::int64_t));
            parent_blocks.push_back(std::move(pblock));
        }
    }
    if (recordParentIds && !parent_blocks.empty())
        surface.AddCellData("surface:parent_cell", std::move(parent_blocks));

    // point_data: subset rows through the same remap; other data maps dropped.
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (detail::rows(a) != num_points)
            continue;  // shape does not match the point table; drop
        std::vector<std::size_t> new_shape = a.Shape();
        new_shape[0] = num_used;
        const std::size_t row_bytes = num_points == 0 ? 0 : a.Nbytes() / num_points;
        NDArray b = NDArray::Uninit(a.Dtype(), std::move(new_shape));
        const std::byte* src = a.Data();
        std::byte* dst = b.Data();
        std::size_t w = 0;
        for (std::size_t i = 0; i < num_points; ++i)
            if (used[i])
                std::memcpy(dst + (w++) * row_bytes, src + i * row_bytes, row_bytes);
        surface.AddPointData(name, std::move(b));
    }

    // The extracted facets are newly created cells one dimension below the
    // input's, so no named region can be carried across. Say so rather than
    // dropping them silently; `record_parent_ids` is the escape hatch for a
    // caller that wants to rebuild a group itself.
    detail::warn_regions_dropped(rMesh, pOpName ? pOpName : "extract_surface");

    return surface;
}

}  // namespace detail

// --- public API -------------------------------------------------------------

Mesh extract_surface(const Mesh& rMesh, bool recordParentIds) {
    return detail::surface_extract(rMesh, /*forceFaceMode=*/false, /*linearize=*/false,
                                   recordParentIds, "extract_surface");
}

bool has_surface_extractable_cells(const Mesh& rMesh) {
    int max_dim = -1;
    for (const auto cb : rMesh.CellRange())
        max_dim = std::max(max_dim, surface_effective_dim(cb));
    const SurfaceMode mode =
        max_dim == 3 ? SurfaceMode::Face : (max_dim == 2 ? SurfaceMode::Edge : SurfaceMode::Face);
    if (max_dim != 2 && max_dim != 3)
        return false;
    for (const auto cb : rMesh.CellRange()) {
        const CellType ct = cell_type_from_name(cb.Type());
        if (!cb.IsRagged() && !cb.IsPolyhedron() && surface_mode_supported(ct, mode))
            return true;
    }
    return false;
}

Mesh extract_skin(const Mesh& rMesh, bool linearize) {
    return detail::surface_extract(rMesh, /*forceFaceMode=*/true, linearize,
                                   /*recordParentIds=*/false, "extract_skin");
}

bool has_skinnable_cells(const Mesh& rMesh) {
    // A polyhedron block carries its own faces, so it is supported by
    // construction -- there is no table for it to be missing from, the same
    // exception `surface_extract`'s own pre-scan makes. This was missed here
    // until v10.3.0: `surface_skin_block_supported` excludes every ragged
    // block including polyhedra, so a mesh with ONLY a polyhedron block (e.g.
    // `subdivide`'s output) reported false and every `has_skinnable_cells`
    // consumer (convertSurfaceOps/convert_surface, and the `skin=true`
    // default on the STL/PLY/SVG/TikZ writers) silently fell back to its
    // unskinned path instead of actually skinning it.
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron() ||
            surface_skin_block_supported(cell_type_from_name(cb.Type()), cb.IsRagged()))
            return true;
    }
    return false;
}

// Declared in operations/surface.hpp; see there for the rationale. Arrays
// whose block layout does not match the mesh are skipped rather than guessed
// at, matching `_viewer._flatten_cell_data` on the Python side.
void gather_cell_data_onto_surface(const Mesh& rSource, Mesh& rSurface) {
    constexpr const char* kParent = "surface:parent_cell";
    if (!rSurface.HasCellData(kParent))
        return;

    const std::size_t n_src_blocks = rSource.NumCellBlocks();
    const std::size_t n_out_blocks = rSurface.NumCellBlocks();

    // Offset of each source block in the global cell numbering, which is what
    // `surface:parent_cell` indexes.
    std::vector<std::size_t> base(n_src_blocks + 1, 0);
    for (std::size_t b = 0; b < n_src_blocks; ++b)
        base[b + 1] = base[b] + rSource.Cells(b).NumCells();
    const std::size_t n_src_cells = base[n_src_blocks];

    for (const auto& name : rSource.CellDataNames()) {
        if (rSource.CellDataNumBlocks(name) != n_src_blocks)
            continue;

        // Flatten the per-block arrays into one buffer indexed by global cell
        // id, requiring a consistent dtype and component count.
        const NDArray& first = rSource.CellData(name, 0);
        const DType dt = first.Dtype();
        const std::size_t item = dtype_size(dt);
        const std::size_t comps = first.Size() / std::max<std::size_t>(1, first.Shape()[0]);
        const std::size_t row = comps * item;

        std::vector<std::byte> flat(n_src_cells * row);
        bool usable = true;
        for (std::size_t b = 0; b < n_src_blocks && usable; ++b) {
            const NDArray& a = rSource.CellData(name, b);
            const std::size_t n = rSource.Cells(b).NumCells();
            usable =
                a.Dtype() == dt && !a.Shape().empty() && a.Shape()[0] == n && a.Size() == n * comps;
            if (usable && n)
                std::memcpy(flat.data() + base[b] * row, a.Data(), n * row);
        }
        if (!usable)
            continue;

        // Gather one output array per boundary block, by owning cell.
        std::vector<NDArray> blocks;
        blocks.reserve(n_out_blocks);
        for (std::size_t b = 0; b < n_out_blocks; ++b) {
            const NDArray& parents = rSurface.CellData(kParent, b);
            const std::size_t n = parents.Size();
            std::vector<std::size_t> shape{n};
            if (comps > 1)
                shape.push_back(comps);
            NDArray out = NDArray::Uninit(dt, shape);
            const auto* pIds = reinterpret_cast<const std::int64_t*>(parents.Data());
            for (std::size_t i = 0; i < n; ++i) {
                const auto id = static_cast<std::size_t>(pIds[i]);
                if (id < n_src_cells)
                    std::memcpy(out.Data() + i * row, flat.data() + id * row, row);
                else
                    std::memset(out.Data() + i * row, 0, row);
            }
            blocks.push_back(std::move(out));
        }
        rSurface.AddCellData(name, std::move(blocks));
    }
}

}  // namespace meshioplusplus
