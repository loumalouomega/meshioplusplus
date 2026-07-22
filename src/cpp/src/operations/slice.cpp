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
// Planar cross-section of a mesh via marching tetrahedra on a simplexified
// input, built entirely through the uniform mesh API so it compiles under every
// mesh backend. See operations/slice.hpp for the contract, in particular the
// degeneracy/winding/determinism rules the pure-numpy twin (_slice.py)
// replicates expression for expression.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

// A canonical (low, high) simplexified-node pair identifying one cut edge; two
// simplices sharing an edge produce the same key, so the crossing point is
// deduped to one output node (the section is watertight).
using SliceEdgeKey = std::pair<std::int64_t, std::int64_t>;

struct SliceEdgeKeyHash {
    std::size_t operator()(const SliceEdgeKey& rKey) const {
        std::size_t h = std::hash<std::int64_t>{}(rKey.first);
        h ^= std::hash<std::int64_t>{}(rKey.second) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

// Local edges of a tetrahedron: e0=(0,1) e1=(0,2) e2=(0,3) e3=(1,2) e4=(1,3) e5=(2,3).
constexpr std::uint8_t kTetEdge[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

// Local edges of a triangle: e0=(0,1) e1=(1,2) e2=(0,2).
constexpr std::uint8_t kTriEdge[3][2] = {{0, 1}, {1, 2}, {0, 2}};

// The crossing-edge ring of a cut simplex, keyed by its below-plane sign mask
// (bit i set iff d_i < 0). The ring lists local edge indices in cyclic order;
// its winding is fixed only up to a global flip, which slice() resolves at
// runtime so every section face's Newell normal points toward the +normal side.
struct SliceRing {
    std::uint8_t mCount;
    std::array<std::uint8_t, 4> mEdges;
};

// Marching-tetrahedra table (16 masks). 0/full masks emit nothing; a single
// isolated vertex crosses 3 edges (triangle); a 2-2 split crosses 4 (quad). The
// quad rings are the cyclic orders derived from the classic two-triangle cases.
constexpr SliceRing kTetRing[16] = {
    /*0x0*/ {0, {{0, 0, 0, 0}}}, /*0x1*/ {3, {{0, 1, 2, 0}}},
    /*0x2*/ {3, {{0, 3, 4, 0}}}, /*0x3*/ {4, {{2, 1, 3, 4}}},
    /*0x4*/ {3, {{1, 3, 5, 0}}}, /*0x5*/ {4, {{0, 2, 5, 3}}},
    /*0x6*/ {4, {{0, 4, 5, 1}}}, /*0x7*/ {3, {{2, 4, 5, 0}}},
    /*0x8*/ {3, {{2, 4, 5, 0}}}, /*0x9*/ {4, {{0, 4, 5, 1}}},
    /*0xA*/ {4, {{0, 2, 5, 3}}}, /*0xB*/ {3, {{1, 3, 5, 0}}},
    /*0xC*/ {4, {{2, 1, 3, 4}}}, /*0xD*/ {3, {{0, 3, 4, 0}}},
    /*0xE*/ {3, {{0, 1, 2, 0}}}, /*0xF*/ {0, {{0, 0, 0, 0}}},
};

// Marching-triangle table (8 masks): a single isolated vertex crosses 2 edges
// (a segment); 0/full masks emit nothing.
constexpr SliceRing kTriRing[8] = {
    /*0x0*/ {0, {{0, 0, 0, 0}}}, /*0x1*/ {2, {{0, 2, 0, 0}}},
    /*0x2*/ {2, {{0, 1, 0, 0}}}, /*0x3*/ {2, {{1, 2, 0, 0}}},
    /*0x4*/ {2, {{1, 2, 0, 0}}}, /*0x5*/ {2, {{0, 1, 0, 0}}},
    /*0x6*/ {2, {{0, 2, 0, 0}}}, /*0x7*/ {0, {{0, 0, 0, 0}}},
};

// One cut simplex's section primitive: the ordered output node ids plus the
// provenance needed to replicate cell_data and record the input cell.
struct SliceFace {
    std::uint8_t mNumVerts;              // 2 (line), 3 (triangle) or 4 (quad)
    std::array<std::int64_t, 4> mNodes;  // section node ids (0-based)
    std::size_t mParentBlock;            // simplexified block index
    std::size_t mParentLocal;            // cell index within that block
    std::int64_t mParentGlobalCell;      // global (block-major) input cell index
};

// Newell normal of a polygon corner ring (twin of test_skin.cpp's helper).
Vec3 slice_newell(const std::vector<Vec3>& rRing) {
    Vec3 n = {0.0, 0.0, 0.0};
    const std::size_t k = rRing.size();
    for (std::size_t i = 0; i < k; ++i) {
        const Vec3& a = rRing[i];
        const Vec3& b = rRing[(i + 1) % k];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    return n;
}

// A staged output cell block: its type, connectivity, and per-cell provenance.
struct SliceOutBlock {
    std::string mType;
    std::size_t mNodesPerCell;
    std::vector<std::int64_t> mConn;
    std::vector<std::size_t> mParentBlock;
    std::vector<std::size_t> mParentLocal;
    std::vector<std::int64_t> mParentGlobalCell;
};

}  // namespace

Mesh slice(const Mesh& rMesh, const SliceOptions& rOptions) {
    const Vec3 origin = {rOptions.mOrigin[0], rOptions.mOrigin[1], rOptions.mOrigin[2]};
    const Vec3 normal = {rOptions.mNormal[0], rOptions.mNormal[1], rOptions.mNormal[2]};
    if (detail::vec3_norm_sq(normal) <= 0.0)
        throw std::invalid_argument("meshio++: slice: normal must be non-zero");

    // Simplexify: every 3D cell becomes a tetra, every 2D cell a triangle, with
    // convert:parent_cell recording the input cell (within its block).
    ConvertCellsResult cc =
        convert_cells(rMesh, {ConvertCellsMode::Simplexify, /*RecordParentIds=*/true});
    const Mesh& simp = cc.mMesh;

    const std::size_t nblocks = simp.NumCellBlocks();
    const std::size_t dim = simp.PointDim();

    // Original (input) global cell base per block: blocks map 1:1 through
    // simplexify, so simp block bi corresponds to input block bi.
    std::vector<std::int64_t> orig_block_base(nblocks, 0);
    {
        std::size_t bi = 0;
        std::int64_t base = 0;
        for (const auto cb : rMesh.CellRange()) {
            if (bi < nblocks)
                orig_block_base[bi] = base;
            base += static_cast<std::int64_t>(cb.NumCells());
            ++bi;
        }
    }

    // The maximum topological dimension decides the cut: a volume mesh (any
    // tetra) sections into surface faces, a 2D surface mesh into line segments.
    int maxdim = 0;
    for (const auto cb : simp.CellRange())
        maxdim = std::max(maxdim, cell_type_dimension(cell_type_from_name(std::string(cb.Type()))));
    const CellType want = maxdim == 3 ? CellType::Tetra : CellType::Triangle;

    // Signed distance of every simplexified node to the plane (the hot loop).
    const std::size_t nnodes_in = simp.NumPoints();
    std::vector<double> dist(nnodes_in, 0.0);
    {
        const NDArray& points = simp.Points();
        parallel_for(nnodes_in, [&](std::size_t g) {
            dist[g] = detail::vec3_dot(
                detail::vec3_sub(detail::read_point(points, dim, static_cast<std::int64_t>(g)),
                                 origin),
                normal);
        });
    }

    const bool has_parent = simp.HasCellData("convert:parent_cell") &&
                            simp.CellDataNumBlocks("convert:parent_cell") == nblocks;

    // --- serial cut + edge-key dedup (the determinism pin) -------------------
    // Iterated in a fixed (block, cell, ring-edge) order; node ids are handed
    // out by a single sweep, so output is identical across backends and threads.
    std::unordered_map<SliceEdgeKey, std::int64_t, SliceEdgeKeyHash> node_id;
    std::vector<SliceEdgeKey> node_edges;
    std::vector<SliceFace> faces;

    if (maxdim == 2 || maxdim == 3) {
        std::size_t bi = 0;
        for (const auto cb : simp.CellRange()) {
            const std::size_t block = bi++;
            const CellType type = cell_type_from_name(std::string(cb.Type()));
            if (type != want || cb.IsRagged())
                continue;
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            const std::size_t ncorners = maxdim == 3 ? 4u : 3u;
            const auto* edge_tbl = maxdim == 3 ? &kTetEdge[0] : &kTriEdge[0];
            const SliceRing* ring_tbl = maxdim == 3 ? kTetRing : kTriRing;
            const NDArray* parent =
                has_parent ? &simp.CellData("convert:parent_cell", block) : nullptr;
            const std::size_t nc = cb.NumCells();
            for (std::size_t c = 0; c < nc; ++c) {
                const std::size_t row = c * npc;
                std::array<std::int64_t, 4> nid{};
                unsigned mask = 0;
                for (std::size_t k = 0; k < ncorners; ++k) {
                    nid[k] = detail::read_int(conn, row + k);
                    if (dist[static_cast<std::size_t>(nid[k])] < 0.0)
                        mask |= (1u << k);
                }
                const SliceRing& ring = ring_tbl[mask];
                if (ring.mCount == 0)
                    continue;
                SliceFace f;
                f.mNumVerts = ring.mCount;
                f.mParentBlock = block;
                f.mParentLocal = c;
                const std::int64_t within =
                    parent ? detail::read_int(*parent, c) : static_cast<std::int64_t>(c);
                f.mParentGlobalCell = orig_block_base[block] + within;
                for (std::size_t e = 0; e < ring.mCount; ++e) {
                    const std::uint8_t ei = ring.mEdges[e];
                    const std::int64_t ga = nid[edge_tbl[ei][0]];
                    const std::int64_t gb = nid[edge_tbl[ei][1]];
                    const SliceEdgeKey key = ga < gb ? SliceEdgeKey{ga, gb} : SliceEdgeKey{gb, ga};
                    auto it = node_id.find(key);
                    std::int64_t id;
                    if (it == node_id.end()) {
                        id = static_cast<std::int64_t>(node_edges.size());
                        node_id.emplace(key, id);
                        node_edges.push_back(key);
                    } else {
                        id = it->second;
                    }
                    f.mNodes[e] = id;
                }
                faces.push_back(f);
            }
        }
    }

    const std::size_t nnodes = node_edges.size();

    Mesh out;

    // --- crossing-point coordinates (parallel) -------------------------------
    NDArray out_points = NDArray::Uninit(DType::Float64, {nnodes, dim});
    {
        const NDArray& points = simp.Points();
        double* dst = out_points.As<double>();
        parallel_for_bw(nnodes, [&](std::size_t i) {
            const std::int64_t lo = node_edges[i].first;
            const std::int64_t hi = node_edges[i].second;
            const double dl = dist[static_cast<std::size_t>(lo)];
            const double dh = dist[static_cast<std::size_t>(hi)];
            const double t = dl / (dl - dh);
            for (std::size_t k = 0; k < dim; ++k) {
                const double pl =
                    detail::read_double(points, static_cast<std::size_t>(lo) * dim + k);
                const double ph =
                    detail::read_double(points, static_cast<std::size_t>(hi) * dim + k);
                dst[i * dim + k] = pl + t * (ph - pl);
            }
        });
    }
    out.AssignPoints(std::move(out_points));

    // --- point_data: interpolate every source array at the same t (Float64) --
    for (const std::string& name : simp.PointDataNames()) {
        const NDArray& a = simp.PointData(name);
        const std::size_t src_rows = detail::rows(a);
        if (src_rows == 0)
            continue;
        const std::size_t comp = a.Size() / src_rows;
        std::vector<std::size_t> shape = a.Shape();
        if (shape.empty())
            shape = {nnodes};
        else
            shape[0] = nnodes;
        NDArray dst_arr = NDArray::Uninit(DType::Float64, std::move(shape));
        double* dst = dst_arr.As<double>();
        parallel_for_bw(nnodes, [&](std::size_t i) {
            const std::int64_t lo = node_edges[i].first;
            const std::int64_t hi = node_edges[i].second;
            const double dl = dist[static_cast<std::size_t>(lo)];
            const double dh = dist[static_cast<std::size_t>(hi)];
            const double t = dl / (dl - dh);
            for (std::size_t k = 0; k < comp; ++k) {
                const double vl = detail::read_double(a, static_cast<std::size_t>(lo) * comp + k);
                const double vh = detail::read_double(a, static_cast<std::size_t>(hi) * comp + k);
                dst[i * comp + k] = vl + t * (vh - vl);
            }
        });
        out.AddPointData(name, std::move(dst_arr));
    }

    if (faces.empty())
        return out;

    // --- bbox-relative degeneracy tolerance ----------------------------------
    double bbdiag = 0.0;
    {
        const double* p = out.Points().As<double>();
        std::vector<double> lo(dim, 0.0), hi(dim, 0.0);
        for (std::size_t k = 0; k < dim; ++k) {
            lo[k] = p[k];
            hi[k] = p[k];
        }
        for (std::size_t i = 0; i < nnodes; ++i)
            for (std::size_t k = 0; k < dim; ++k) {
                lo[k] = std::min(lo[k], p[i * dim + k]);
                hi[k] = std::max(hi[k], p[i * dim + k]);
            }
        double s = 0.0;
        for (std::size_t k = 0; k < dim; ++k)
            s += (hi[k] - lo[k]) * (hi[k] - lo[k]);
        bbdiag = std::sqrt(s);
    }
    const double rel = 1e-10;
    const double area_tol = rel * bbdiag * bbdiag;  // |Newell| ~ 2 * area
    const double len_tol = rel * bbdiag;

    // --- stage output blocks: orient by Newell flip, drop collapses ----------
    // Fixed output order: triangle, quad (volume) / line (surface).
    SliceOutBlock tri_blk{cell_type_name(CellType::Triangle), 3, {}, {}, {}, {}};
    SliceOutBlock quad_blk{cell_type_name(CellType::Quad), 4, {}, {}, {}, {}};
    SliceOutBlock line_blk{cell_type_name(CellType::Line), 2, {}, {}, {}, {}};
    const double* pts = out.Points().As<double>();
    for (const SliceFace& f : faces) {
        if (f.mNumVerts >= 3) {
            std::vector<Vec3> ring(f.mNumVerts);
            for (std::size_t v = 0; v < f.mNumVerts; ++v) {
                const std::size_t nd = static_cast<std::size_t>(f.mNodes[v]);
                ring[v] = detail::read_point(out.Points(), dim, static_cast<std::int64_t>(nd));
            }
            const Vec3 nrm = slice_newell(ring);
            if (detail::vec3_norm(nrm) < area_tol)
                continue;
            const bool flip = detail::vec3_dot(nrm, normal) < 0.0;
            SliceOutBlock& blk = f.mNumVerts == 3 ? tri_blk : quad_blk;
            for (std::size_t v = 0; v < f.mNumVerts; ++v) {
                const std::size_t src = flip ? (f.mNumVerts - 1 - v) : v;
                blk.mConn.push_back(f.mNodes[src]);
            }
            blk.mParentBlock.push_back(f.mParentBlock);
            blk.mParentLocal.push_back(f.mParentLocal);
            blk.mParentGlobalCell.push_back(f.mParentGlobalCell);
        } else {  // line segment
            const std::size_t a = static_cast<std::size_t>(f.mNodes[0]);
            const std::size_t b = static_cast<std::size_t>(f.mNodes[1]);
            double s = 0.0;
            for (std::size_t k = 0; k < dim; ++k) {
                const double dd = pts[a * dim + k] - pts[b * dim + k];
                s += dd * dd;
            }
            if (std::sqrt(s) < len_tol)
                continue;
            line_blk.mConn.push_back(f.mNodes[0]);
            line_blk.mConn.push_back(f.mNodes[1]);
            line_blk.mParentBlock.push_back(f.mParentBlock);
            line_blk.mParentLocal.push_back(f.mParentLocal);
            line_blk.mParentGlobalCell.push_back(f.mParentGlobalCell);
        }
    }

    std::vector<SliceOutBlock*> out_blocks;
    if (!tri_blk.mConn.empty())
        out_blocks.push_back(&tri_blk);
    if (!quad_blk.mConn.empty())
        out_blocks.push_back(&quad_blk);
    if (!line_blk.mConn.empty())
        out_blocks.push_back(&line_blk);

    for (SliceOutBlock* blk : out_blocks) {
        const std::size_t ncells = blk->mConn.size() / blk->mNodesPerCell;
        NDArray conn = NDArray::Uninit(DType::Int64, {ncells, blk->mNodesPerCell});
        std::memcpy(conn.Data(), blk->mConn.data(), blk->mConn.size() * sizeof(std::int64_t));
        out.AddCellBlock(blk->mType, std::move(conn));
    }

    // --- cell_data: replicate each parent's row to its section cells ----------
    for (const std::string& name : simp.CellDataNames()) {
        if (name == "convert:parent_cell")
            continue;
        if (simp.CellDataNumBlocks(name) != nblocks)
            continue;  // not per-cell data; dropped (new topology)
        const NDArray& a0 = simp.CellData(name, 0);
        const DType dt = a0.Dtype();
        std::size_t comp = 1;
        for (std::size_t d = 1; d < a0.Shape().size(); ++d)
            comp *= a0.Shape()[d];
        const std::size_t row_bytes = comp * dtype_size(dt);
        std::vector<NDArray> blocks;
        blocks.reserve(out_blocks.size());
        for (SliceOutBlock* blk : out_blocks) {
            const std::size_t ncells = blk->mParentGlobalCell.size();
            std::vector<std::size_t> shape = a0.Shape();
            if (shape.empty())
                shape = {ncells};
            else
                shape[0] = ncells;
            NDArray dst_arr = NDArray::Uninit(dt, std::move(shape));
            std::byte* dst = dst_arr.Data();
            for (std::size_t c = 0; c < ncells; ++c) {
                const NDArray& src = simp.CellData(name, blk->mParentBlock[c]);
                std::memcpy(dst + c * row_bytes, src.Data() + blk->mParentLocal[c] * row_bytes,
                            row_bytes);
            }
            blocks.push_back(std::move(dst_arr));
        }
        out.AddCellData(name, std::move(blocks));
    }

    // --- slice:parent_cell provenance ----------------------------------------
    if (rOptions.mRecordParentIds) {
        std::vector<NDArray> blocks;
        blocks.reserve(out_blocks.size());
        for (SliceOutBlock* blk : out_blocks) {
            NDArray a = NDArray::Uninit(DType::Int64, {blk->mParentGlobalCell.size()});
            if (!blk->mParentGlobalCell.empty())
                std::memcpy(a.Data(), blk->mParentGlobalCell.data(),
                            blk->mParentGlobalCell.size() * sizeof(std::int64_t));
            blocks.push_back(std::move(a));
        }
        out.AddCellData("slice:parent_cell", std::move(blocks));
    }

    return out;
}

}  // namespace meshioplusplus
