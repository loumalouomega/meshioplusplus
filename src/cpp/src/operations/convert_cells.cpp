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
// Element-representation conversion: linearize (drop higher-order nodes),
// simplexify (decompose into same-dimension simplices via fixed templates), and
// elevate (promote linear cells to serendipity quadratic, adding one node per
// unique edge). Built entirely through the uniform mesh API so it compiles under
// every mesh backend. See operations/convert_cells.hpp for the contract.
//
// Determinism: the decomposition templates are fixed, and the elevate mid-edge
// numbering comes from a serial pass over a parallel-filled record buffer --
// surface.cpp's phase-split idiom -- never from a concurrent hash insert. Output
// is therefore byte-identical across backends and thread counts.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// --- tables ------------------------------------------------------------------

// The linear base of a cell type (identity for types that are already linear or
// have no linear counterpart). Corner nodes always come first in meshio/VTK
// orderings, so linearizing is "keep the leading cell_type_num_nodes(base)
// columns" -- no permutation is ever needed.
CellType ccells_linear_base(CellType type) {
    switch (type) {
        case CellType::Line3:
        case CellType::Line4:
            return CellType::Line;
        case CellType::Triangle6:
        case CellType::Triangle10:
            return CellType::Triangle;
        case CellType::Quad8:
        case CellType::Quad9:
        case CellType::Quad16:
            return CellType::Quad;
        case CellType::Tetra10:
        case CellType::Tetra20:
            return CellType::Tetra;
        case CellType::Hexahedron20:
        case CellType::Hexahedron24:
        case CellType::Hexahedron27:
        case CellType::Hexahedron64:
            return CellType::Hexahedron;
        case CellType::Wedge15:
        case CellType::Wedge18:
            return CellType::Wedge;
        case CellType::Pyramid13:
        case CellType::Pyramid14:
            return CellType::Pyramid;
        default:
            return type;
    }
}

// One decomposition template: `mNumChildren` simplices of `mChildType`, each
// listing `mNodesPerChild` local node indices into the parent's connectivity.
struct CcellsSimplexTemplate {
    CellType mChildType;
    std::uint8_t mNodesPerChild;
    std::uint8_t mNumChildren;
    std::vector<std::uint8_t> mNodes;  // mNumChildren * mNodesPerChild entries
};

// The decomposition templates, all with consistently positive orientation for a
// well-oriented parent (pinned by tests/cpp/test_convert_cells.cpp).
//
// Conventions come from detail/cell_faces.hpp: tetra base (0,1,2) normal points
// toward apex 3; hexahedron base (0,1,2,3) normal points toward the top
// (4,5,6,7); wedge base (0,1,2) normal points toward the top (3,4,5); pyramid
// base (0,1,2,3) normal points toward apex 4.
//
// The hexahedron template is the canonical Freudenthal fan around the main
// diagonal 0-6: every one of the six tetra contains the edge (0,6), and the
// remaining two nodes walk the six faces not touching that diagonal. It is
// conforming across a shared face only if neighbouring hexahedra agree on the
// diagonal, which is why the diagonal is fixed at 0-6 rather than chosen per
// cell -- a deterministic, documented choice.
const std::unordered_map<CellType, CcellsSimplexTemplate>& ccells_simplex_table() {
    static const std::unordered_map<CellType, CcellsSimplexTemplate> table = [] {
        std::unordered_map<CellType, CcellsSimplexTemplate> t;
        t[CellType::Quad] = {CellType::Triangle, 3, 2, {0, 1, 2, 0, 2, 3}};
        t[CellType::Pyramid] = {CellType::Tetra, 4, 2, {0, 1, 2, 4, 0, 2, 3, 4}};
        t[CellType::Wedge] = {CellType::Tetra, 4, 3, {0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5}};
        t[CellType::Hexahedron] = {CellType::Tetra, 4, 6, {0, 1, 2, 6, 0, 2, 3, 6, 0, 3, 7, 6,
                                                           0, 7, 4, 6, 0, 4, 5, 6, 0, 5, 1, 6}};
        return t;
    }();
    return table;
}

const CcellsSimplexTemplate* ccells_simplex_template(CellType type) {
    const auto& table = ccells_simplex_table();
    auto it = table.find(type);
    return it == table.end() ? nullptr : &it->second;
}

// A cell type that is already a simplex of its own dimension needs no work.
bool ccells_is_simplex(CellType type) {
    return type == CellType::Line || type == CellType::Triangle || type == CellType::Tetra ||
           type == CellType::Vertex;
}

// The serendipity-quadratic promotion of a linear type: the target type, its
// corner count, and its edges **in the order the target's mid-node slots
// expect** (target node `mNumCorners + k` is the midpoint of edge `k`).
struct CcellsElevateSpec {
    CellType mTarget;
    std::uint8_t mNumCorners;
    std::vector<std::array<std::uint8_t, 2>> mEdges;
};

// Edge orderings come from detail/cell_subdivision.hpp, which is the single
// owner of "the edges of type T, in the order T's quadratic mid-node slots
// expect" -- shared with operations/refine.cpp, which needs the identical order
// for the same reason. (It in turn delegates the 2D rows to
// detail/cell_edges.hpp.) Keeping one table is what stops the two operations
// from drifting: a mismatched order still yields a valid, positively-oriented
// mesh, just one whose quadratic mid-nodes are permuted.
const std::unordered_map<CellType, CcellsElevateSpec>& ccells_elevate_table() {
    static const std::unordered_map<CellType, CcellsElevateSpec> table = [] {
        auto edges_of = [](CellType source) { return detail::cell_refine_edges(source); };
        std::unordered_map<CellType, CcellsElevateSpec> t;
        t[CellType::Line] = {CellType::Line3, 2, edges_of(CellType::Line)};
        t[CellType::Triangle] = {CellType::Triangle6, 3, edges_of(CellType::Triangle)};
        t[CellType::Quad] = {CellType::Quad8, 4, edges_of(CellType::Quad)};
        t[CellType::Tetra] = {CellType::Tetra10, 4, edges_of(CellType::Tetra)};
        t[CellType::Hexahedron] = {CellType::Hexahedron20, 8, edges_of(CellType::Hexahedron)};
        t[CellType::Wedge] = {CellType::Wedge15, 6, edges_of(CellType::Wedge)};
        t[CellType::Pyramid] = {CellType::Pyramid13, 5, edges_of(CellType::Pyramid)};
        return t;
    }();
    return table;
}

const CcellsElevateSpec* ccells_elevate_spec(CellType type) {
    const auto& table = ccells_elevate_table();
    auto it = table.find(type);
    return it == table.end() ? nullptr : &it->second;
}

// The full-Lagrange targets this version deliberately does not produce: they
// need face-centre and body-centre nodes, not just mid-edge ones.
bool ccells_is_full_lagrange(CellType type) {
    return type == CellType::Quad9 || type == CellType::Hexahedron27 || type == CellType::Wedge18 ||
           type == CellType::Pyramid14;
}

// --- output block staging ----------------------------------------------------

// One staged output cell block, in whichever of the three storage shapes the
// uniform API accepts. Emitted by ccells_emit_block, optionally through a point
// remap.
struct CcellsOutBlock {
    std::string mType;
    std::size_t mNodesPerCell = 0;  // rectangular blocks only
    std::vector<std::int64_t> mConn;
    std::vector<std::vector<std::int64_t>> mPolygonRows;
    std::vector<std::vector<std::vector<std::int64_t>>> mPolyhedronCells;
    bool mIsRagged = false;
    bool mIsPolyhedron = false;
};

// Stage an input block unchanged (used for every pass-through case).
CcellsOutBlock ccells_stage_passthrough(const Mesh::CellView& rBlock) {
    CcellsOutBlock out;
    out.mType = std::string(rBlock.Type());
    if (rBlock.IsPolyhedron()) {
        out.mIsRagged = true;
        out.mIsPolyhedron = true;
        out.mPolyhedronCells.resize(rBlock.NumCells());
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c) {
            for (std::size_t f = 0; f < rBlock.NumFaces(c); ++f) {
                auto face = rBlock.Face(c, f);
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
        parallel_for_bw(rBlock.NumCells(), [&](std::size_t c) {
            for (std::size_t k = 0; k < npc; ++k)
                out.mConn[c * npc + k] = detail::read_int(conn, c * npc + k);
        });
    }
    return out;
}

// Add a staged block to `rOut`, translating node ids through `pRemap` when given.
void ccells_emit_block(Mesh& rOut, CcellsOutBlock& rBlock,
                       const std::vector<std::int64_t>* pRemap) {
    auto map = [pRemap](std::int64_t id) {
        return pRemap ? (*pRemap)[static_cast<std::size_t>(id)] : id;
    };
    if (rBlock.mIsPolyhedron) {
        if (pRemap)
            for (auto& cell : rBlock.mPolyhedronCells)
                for (auto& face : cell)
                    for (std::int64_t& id : face)
                        id = map(id);
        rOut.AddPolyhedronBlock(rBlock.mType, std::move(rBlock.mPolyhedronCells));
    } else if (rBlock.mIsRagged) {
        if (pRemap)
            for (auto& row : rBlock.mPolygonRows)
                for (std::int64_t& id : row)
                    id = map(id);
        rOut.AddPolygonBlock(rBlock.mType, std::move(rBlock.mPolygonRows));
    } else {
        const std::size_t npc = rBlock.mNodesPerCell;
        const std::size_t ncells = npc == 0 ? 0 : rBlock.mConn.size() / npc;
        NDArray conn = NDArray::Uninit(DType::Int64, {ncells, npc});
        std::int64_t* dst = conn.As<std::int64_t>();
        parallel_for_bw(rBlock.mConn.size(), [&](std::size_t i) { dst[i] = map(rBlock.mConn[i]); });
        rOut.AddCellBlock(rBlock.mType, std::move(conn));
    }
}

// Visit every node id of a staged block (for the "which points survive" scan).
template <class TFunc>
void ccells_visit_block_nodes(const CcellsOutBlock& rBlock, TFunc&& rFunc) {
    if (rBlock.mIsPolyhedron) {
        for (const auto& cell : rBlock.mPolyhedronCells)
            for (const auto& face : cell)
                for (std::int64_t id : face)
                    rFunc(id);
    } else if (rBlock.mIsRagged) {
        for (const auto& row : rBlock.mPolygonRows)
            for (std::int64_t id : row)
                rFunc(id);
    } else {
        for (std::int64_t id : rBlock.mConn)
            rFunc(id);
    }
}

// --- shared data plumbing ----------------------------------------------------

NDArray ccells_int64_vector(const std::vector<std::int64_t>& rValues) {
    NDArray out = NDArray::Uninit(DType::Int64, {rValues.size()});
    if (!rValues.empty())
        std::memcpy(out.Data(), rValues.data(), rValues.size() * sizeof(std::int64_t));
    return out;
}

// Per-block identity cell maps (Linearize/Elevate: one child per parent).
std::vector<NDArray> ccells_identity_cell_maps(const Mesh& rMesh) {
    std::vector<NDArray> maps;
    maps.reserve(rMesh.NumCellBlocks());
    for (const auto cb : rMesh.CellRange()) {
        NDArray m = NDArray::Uninit(DType::Int64, {cb.NumCells()});
        std::int64_t* dst = m.As<std::int64_t>();
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            dst[c] = static_cast<std::int64_t>(c);
        maps.push_back(std::move(m));
    }
    return maps;
}

NDArray ccells_identity_point_map(std::size_t NumPoints) {
    NDArray m = NDArray::Uninit(DType::Int64, {NumPoints});
    std::int64_t* dst = m.As<std::int64_t>();
    for (std::size_t i = 0; i < NumPoints; ++i)
        dst[i] = static_cast<std::int64_t>(i);
    return m;
}

// Copy cell_data verbatim (Linearize/Elevate keep the cell count and the 1:1
// block correspondence, so the arrays transfer untouched).
void ccells_copy_cell_data(const Mesh& rIn, Mesh& rOut) {
    for (const std::string& name : rIn.CellDataNames()) {
        std::vector<NDArray> blocks;
        blocks.reserve(rIn.CellDataNumBlocks(name));
        for (std::size_t b = 0; b < rIn.CellDataNumBlocks(name); ++b)
            blocks.push_back(detail::data_owned_copy(rIn.CellData(name, b)));
        rOut.AddCellData(name, std::move(blocks));
    }
}

void ccells_copy_point_data(const Mesh& rIn, Mesh& rOut) {
    for (const std::string& name : rIn.PointDataNames())
        rOut.AddPointData(name, detail::data_owned_copy(rIn.PointData(name)));
}

void ccells_copy_field_data(const Mesh& rIn, Mesh& rOut) {
    for (const std::string& name : rIn.FieldDataNames())
        rOut.AddFieldData(name, detail::data_owned_copy(rIn.FieldData(name)));
}

// --- mode: linearize ---------------------------------------------------------

ConvertCellsResult ccells_linearize(const Mesh& rMesh, bool RecordParentIds) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<CcellsOutBlock> staged;
    staged.reserve(nblocks);

    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsRagged()) {
            // Polygon / polyhedron blocks are already linear.
            staged.push_back(ccells_stage_passthrough(cb));
            continue;
        }
        const CellType type = cell_type_from_name(std::string(cb.Type()));
        const CellType base = ccells_linear_base(type);
        if (base == type) {
            staged.push_back(ccells_stage_passthrough(cb));
            continue;
        }
        const std::size_t npc_in = cb.NodesPerCell();
        const std::size_t npc_out = static_cast<std::size_t>(cell_type_num_nodes(base));
        const NDArray& conn = cb.Conn();
        const std::size_t ncells = cb.NumCells();

        CcellsOutBlock out;
        out.mType = cell_type_name(base);
        out.mNodesPerCell = npc_out;
        out.mConn.resize(ncells * npc_out);
        std::vector<std::int64_t>& dst = out.mConn;
        parallel_for_bw(ncells, [&](std::size_t c) {
            for (std::size_t k = 0; k < npc_out; ++k)
                dst[c * npc_out + k] = detail::read_int(conn, c * npc_in + k);
        });
        staged.push_back(std::move(out));
    }

    // Which points survive, in ascending original order.
    const std::size_t num_points = rMesh.NumPoints();
    std::vector<char> used(num_points, 0);
    for (const CcellsOutBlock& block : staged)
        ccells_visit_block_nodes(block,
                                 [&](std::int64_t id) { used[static_cast<std::size_t>(id)] = 1; });
    std::vector<std::int64_t> remap(num_points, -1);
    std::size_t num_used = 0;
    for (std::size_t i = 0; i < num_points; ++i)
        if (used[i])
            remap[i] = static_cast<std::int64_t>(num_used++);

    Mesh out;

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
        out.AssignPoints(std::move(new_points));
    }

    for (CcellsOutBlock& block : staged)
        ccells_emit_block(out, block, &remap);

    // point_data: the same row subset; arrays that do not match the point table
    // are copied verbatim (they are not per-point data).
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (detail::rows(a) != num_points) {
            out.AddPointData(name, detail::data_owned_copy(a));
            continue;
        }
        std::vector<std::size_t> shape = a.Shape();
        shape[0] = num_used;
        const std::size_t row_bytes = num_points == 0 ? 0 : a.Nbytes() / num_points;
        NDArray b = NDArray::Uninit(a.Dtype(), std::move(shape));
        const std::byte* src = a.Data();
        std::byte* dst = b.Data();
        std::size_t w = 0;
        for (std::size_t i = 0; i < num_points; ++i)
            if (used[i])
                std::memcpy(dst + (w++) * row_bytes, src + i * row_bytes, row_bytes);
        out.AddPointData(name, std::move(b));
    }

    ccells_copy_cell_data(rMesh, out);
    ccells_copy_field_data(rMesh, out);

    if (RecordParentIds) {
        std::vector<NDArray> blocks;
        blocks.reserve(nblocks);
        for (const auto cb : rMesh.CellRange()) {
            NDArray a = NDArray::Uninit(DType::Int64, {cb.NumCells(), 1});
            std::int64_t* dst = a.As<std::int64_t>();
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                dst[c] = static_cast<std::int64_t>(c);
            blocks.push_back(std::move(a));
        }
        if (!blocks.empty())
            out.AddCellData("convert:parent_cell", std::move(blocks));
    }

    ConvertCellsResult res;
    res.mMesh = std::move(out);
    res.mPointMap = ccells_int64_vector(remap);
    res.mCellMaps = ccells_identity_cell_maps(rMesh);
    return res;
}

// --- mode: simplexify --------------------------------------------------------

ConvertCellsResult ccells_simplexify(const Mesh& rMesh, bool RecordParentIds) {
    // Higher-order input is linearized first (which also prunes the mid nodes
    // the decomposition would otherwise orphan); the decomposition then runs on
    // a purely linear mesh.
    bool needs_linearize = false;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            throw std::invalid_argument("convert_cells: cannot simplexify polyhedron cell block '" +
                                        std::string(cb.Type()) +
                                        "' (2-level ragged blocks have no simplex template)");
        if (cb.IsRagged())
            continue;
        const CellType type = cell_type_from_name(std::string(cb.Type()));
        if (ccells_linear_base(type) != type)
            needs_linearize = true;
    }

    ConvertCellsResult prepared;
    if (needs_linearize)
        prepared = ccells_linearize(rMesh, /*RecordParentIds=*/false);
    const Mesh& mesh = needs_linearize ? prepared.mMesh : rMesh;

    const std::size_t nblocks = mesh.NumCellBlocks();
    std::vector<CcellsOutBlock> staged;
    staged.reserve(nblocks);
    std::vector<std::vector<std::int64_t>> first_child(nblocks);
    std::vector<std::vector<std::int64_t>> parent_of_child(nblocks);

    std::size_t bi = 0;
    for (const auto cb : mesh.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        std::vector<std::int64_t>& firsts = first_child[bi];
        std::vector<std::int64_t>& parents = parent_of_child[bi];
        firsts.resize(ncells);

        const CellType type = cell_type_from_name(std::string(cb.Type()));

        // Polygons fan into (n-2) triangles around node 0. This covers both
        // storage shapes: a jagged block (rows of differing length) and a
        // rectangular block of uniform n-gons, which is how a polygon block
        // arrives whenever every cell happens to have the same node count.
        if (cb.IsRagged() || type == CellType::Polygon) {
            CcellsOutBlock out;
            out.mType = cell_type_name(CellType::Triangle);
            out.mNodesPerCell = 3;
            const bool ragged = cb.IsRagged();
            const NDArray* conn = ragged ? nullptr : &cb.Conn();
            const std::size_t npc = ragged ? 0 : cb.NodesPerCell();
            for (std::size_t c = 0; c < ncells; ++c) {
                firsts[c] = static_cast<std::int64_t>(out.mConn.size() / 3);
                const std::size_t n = ragged ? cb.RowSize(c) : npc;
                auto node = [&](std::size_t k) {
                    return ragged ? cb.Row(c)[k] : detail::read_int(*conn, c * npc + k);
                };
                for (std::size_t k = 1; k + 1 < n; ++k) {
                    out.mConn.push_back(node(0));
                    out.mConn.push_back(node(k));
                    out.mConn.push_back(node(k + 1));
                    parents.push_back(static_cast<std::int64_t>(c));
                }
            }
            staged.push_back(std::move(out));
            ++bi;
            continue;
        }

        const CcellsSimplexTemplate* tpl = ccells_simplex_template(type);
        if (tpl == nullptr) {
            if (!ccells_is_simplex(type))
                log::warn(
                    "convert_cells: cell block '{}' has no simplex template; passing it through "
                    "unchanged.",
                    cb.Type());
            staged.push_back(ccells_stage_passthrough(cb));
            for (std::size_t c = 0; c < ncells; ++c) {
                firsts[c] = static_cast<std::int64_t>(c);
                parents.push_back(static_cast<std::int64_t>(c));
            }
            ++bi;
            continue;
        }

        const std::size_t npc_in = cb.NodesPerCell();
        const std::size_t npc_out = tpl->mNodesPerChild;
        const std::size_t per_cell = tpl->mNumChildren;
        const NDArray& conn = cb.Conn();

        CcellsOutBlock out;
        out.mType = cell_type_name(tpl->mChildType);
        out.mNodesPerCell = npc_out;
        out.mConn.resize(ncells * per_cell * npc_out);
        parents.resize(ncells * per_cell);
        std::vector<std::int64_t>& dst = out.mConn;
        const std::vector<std::uint8_t>& nodes = tpl->mNodes;
        parallel_for_bw(ncells, [&](std::size_t c) {
            const std::size_t base_in = c * npc_in;
            const std::size_t base_out = c * per_cell * npc_out;
            for (std::size_t child = 0; child < per_cell; ++child) {
                for (std::size_t k = 0; k < npc_out; ++k)
                    dst[base_out + child * npc_out + k] =
                        detail::read_int(conn, base_in + nodes[child * npc_out + k]);
                parents[c * per_cell + child] = static_cast<std::int64_t>(c);
            }
        });
        for (std::size_t c = 0; c < ncells; ++c)
            firsts[c] = static_cast<std::int64_t>(c * per_cell);
        staged.push_back(std::move(out));
        ++bi;
    }

    Mesh out;
    out.AssignPoints(detail::data_owned_copy(mesh.Points()));
    for (CcellsOutBlock& block : staged)
        ccells_emit_block(out, block, nullptr);
    ccells_copy_point_data(mesh, out);
    ccells_copy_field_data(mesh, out);

    // cell_data: replicate each parent's row to its children. An array whose
    // block count does not match the mesh is not per-cell data, so it is copied
    // verbatim rather than mangled.
    for (const std::string& name : mesh.CellDataNames()) {
        const std::size_t ndata = mesh.CellDataNumBlocks(name);
        std::vector<NDArray> blocks;
        blocks.reserve(ndata);
        for (std::size_t b = 0; b < ndata; ++b) {
            const NDArray& a = mesh.CellData(name, b);
            const std::size_t in_rows = detail::rows(a);
            if (ndata != nblocks || in_rows == 0 || parent_of_child[b].empty()) {
                blocks.push_back(detail::data_owned_copy(a));
                continue;
            }
            const std::vector<std::int64_t>& parents = parent_of_child[b];
            std::vector<std::size_t> shape = a.Shape();
            shape[0] = parents.size();
            const std::size_t row_bytes = a.Nbytes() / in_rows;
            NDArray replicated = NDArray::Uninit(a.Dtype(), std::move(shape));
            const std::byte* src = a.Data();
            std::byte* dst = replicated.Data();
            parallel_for_bw(parents.size(), [&](std::size_t i) {
                std::memcpy(dst + i * row_bytes,
                            src + static_cast<std::size_t>(parents[i]) * row_bytes, row_bytes);
            });
            blocks.push_back(std::move(replicated));
        }
        out.AddCellData(name, std::move(blocks));
    }

    if (RecordParentIds) {
        std::vector<NDArray> blocks;
        blocks.reserve(nblocks);
        for (std::size_t b = 0; b < nblocks; ++b)
            blocks.push_back(ccells_int64_vector(parent_of_child[b]));
        if (!blocks.empty())
            out.AddCellData("convert:parent_cell", std::move(blocks));
    }

    ConvertCellsResult res;
    res.mMesh = std::move(out);
    res.mPointMap = needs_linearize ? std::move(prepared.mPointMap)
                                    : ccells_identity_point_map(rMesh.NumPoints());
    res.mCellMaps.reserve(nblocks);
    for (std::size_t b = 0; b < nblocks; ++b)
        res.mCellMaps.push_back(ccells_int64_vector(first_child[b]));
    return res;
}

// --- mode: elevate -----------------------------------------------------------

// A canonical (low, high) node pair identifying one edge.
using CcellsEdgeKey = std::pair<std::int64_t, std::int64_t>;

struct CcellsEdgeKeyHash {
    std::size_t operator()(const CcellsEdgeKey& rKey) const {
        std::size_t h = std::hash<std::int64_t>{}(rKey.first);
        h ^= std::hash<std::int64_t>{}(rKey.second) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

// Per-block bookkeeping for the elevate pass.
struct CcellsElevateBlock {
    const CcellsElevateSpec* mpSpec = nullptr;  // null => pass the block through
    std::size_t mFirstSlot = 0;                 // offset into the flat edge buffer
    std::size_t mNumEdges = 0;
};

ConvertCellsResult ccells_elevate(const Mesh& rMesh, bool RecordParentIds) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<CcellsElevateBlock> descs(nblocks);

    std::size_t total_slots = 0;
    std::size_t bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        CcellsElevateBlock& d = descs[bi];
        if (!cb.IsRagged()) {
            const CellType type = cell_type_from_name(std::string(cb.Type()));
            if (ccells_is_full_lagrange(type))
                throw std::invalid_argument(
                    "convert_cells: cannot elevate to the full-Lagrange type '" +
                    cell_type_name(type) +
                    "' (face/body-centre nodes are not supported in this version)");
            d.mpSpec = ccells_elevate_spec(type);
            if (d.mpSpec != nullptr) {
                d.mNumEdges = d.mpSpec->mEdges.size();
                d.mFirstSlot = total_slots;
                total_slots += cb.NumCells() * d.mNumEdges;
            }
        }
        ++bi;
    }

    // --- phase 1: fill the edge keys into disjoint slots (parallel-safe) ---
    std::vector<CcellsEdgeKey> keys(total_slots);
    bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const CcellsElevateBlock& d = descs[bi++];
        if (d.mpSpec == nullptr)
            continue;
        const NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t nedges = d.mNumEdges;
        const std::vector<std::array<std::uint8_t, 2>>& edges = d.mpSpec->mEdges;
        parallel_for(cb.NumCells() * nedges, [&](std::size_t j) {
            const std::size_t cell = j / nedges;
            const std::size_t slot = j % nedges;
            const std::int64_t a = detail::read_int(conn, cell * npc + edges[slot][0]);
            const std::int64_t b = detail::read_int(conn, cell * npc + edges[slot][1]);
            keys[d.mFirstSlot + j] = a < b ? CcellsEdgeKey{a, b} : CcellsEdgeKey{b, a};
        });
    }

    // --- phase 2: serial dedup in stored order -> deterministic numbering ---
    const std::size_t num_points = rMesh.NumPoints();
    std::unordered_map<CcellsEdgeKey, std::int64_t, CcellsEdgeKeyHash> edge_id;
    edge_id.reserve(total_slots * 2);
    std::vector<CcellsEdgeKey> new_edges;
    std::vector<std::int64_t> slot_id(total_slots);
    for (std::size_t i = 0; i < total_slots; ++i) {
        auto it = edge_id.find(keys[i]);
        if (it == edge_id.end()) {
            const std::int64_t id = static_cast<std::int64_t>(num_points + new_edges.size());
            edge_id.emplace(keys[i], id);
            new_edges.push_back(keys[i]);
            slot_id[i] = id;
        } else {
            slot_id[i] = it->second;
        }
    }

    Mesh out;

    // --- phase 3: points = originals + one midpoint per unique edge ---
    const NDArray& points = rMesh.Points();
    const std::size_t dim = detail::cols(points);
    {
        NDArray new_points = NDArray::Uninit(points.Dtype(), {num_points + new_edges.size(), dim});
        std::memcpy(new_points.Data(), points.Data(), points.Nbytes());
        parallel_for_bw(new_edges.size(), [&](std::size_t i) {
            const std::size_t a = static_cast<std::size_t>(new_edges[i].first);
            const std::size_t b = static_cast<std::size_t>(new_edges[i].second);
            for (std::size_t k = 0; k < dim; ++k) {
                const double v = 0.5 * (detail::read_double(points, a * dim + k) +
                                        detail::read_double(points, b * dim + k));
                detail::write_double(new_points, (num_points + i) * dim + k, v);
            }
        });
        out.AssignPoints(std::move(new_points));
    }

    // --- phase 4: connectivity = corners + the mid node of each edge ---
    bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const CcellsElevateBlock& d = descs[bi++];
        if (d.mpSpec == nullptr) {
            CcellsOutBlock passthrough = ccells_stage_passthrough(cb);
            ccells_emit_block(out, passthrough, nullptr);
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t npc_in = cb.NodesPerCell();
        const std::size_t ncorners = d.mpSpec->mNumCorners;
        const std::size_t nedges = d.mNumEdges;
        const std::size_t npc_out = ncorners + nedges;
        const std::size_t ncells = cb.NumCells();

        NDArray block = NDArray::Uninit(DType::Int64, {ncells, npc_out});
        std::int64_t* dst = block.As<std::int64_t>();
        parallel_for_bw(ncells, [&](std::size_t c) {
            for (std::size_t k = 0; k < ncorners; ++k)
                dst[c * npc_out + k] = detail::read_int(conn, c * npc_in + k);
            for (std::size_t e = 0; e < nedges; ++e)
                dst[c * npc_out + ncorners + e] = slot_id[d.mFirstSlot + c * nedges + e];
        });
        out.AddCellBlock(cell_type_name(d.mpSpec->mTarget), std::move(block));
    }

    // --- phase 5: point_data = originals + the endpoint mean per new node ---
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (detail::rows(a) != num_points) {
            out.AddPointData(name, detail::data_owned_copy(a));
            continue;
        }
        const std::size_t ncomp = num_points == 0 ? 0 : a.Size() / num_points;
        std::vector<std::size_t> shape = a.Shape();
        shape[0] = num_points + new_edges.size();
        NDArray b = NDArray::Uninit(a.Dtype(), std::move(shape));
        std::memcpy(b.Data(), a.Data(), a.Nbytes());
        parallel_for_bw(new_edges.size(), [&](std::size_t i) {
            const std::size_t p = static_cast<std::size_t>(new_edges[i].first);
            const std::size_t q = static_cast<std::size_t>(new_edges[i].second);
            for (std::size_t k = 0; k < ncomp; ++k) {
                const double v = 0.5 * (detail::read_double(a, p * ncomp + k) +
                                        detail::read_double(a, q * ncomp + k));
                detail::write_double(b, (num_points + i) * ncomp + k, v);
            }
        });
        out.AddPointData(name, std::move(b));
    }

    ccells_copy_cell_data(rMesh, out);
    ccells_copy_field_data(rMesh, out);

    if (RecordParentIds) {
        std::vector<NDArray> blocks;
        blocks.reserve(nblocks);
        for (const auto cb : rMesh.CellRange()) {
            NDArray a = NDArray::Uninit(DType::Int64, {cb.NumCells(), 1});
            std::int64_t* dst = a.As<std::int64_t>();
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                dst[c] = static_cast<std::int64_t>(c);
            blocks.push_back(std::move(a));
        }
        if (!blocks.empty())
            out.AddCellData("convert:parent_cell", std::move(blocks));
    }

    ConvertCellsResult res;
    res.mMesh = std::move(out);
    res.mPointMap = ccells_identity_point_map(num_points);
    res.mCellMaps = ccells_identity_cell_maps(rMesh);
    return res;
}

}  // namespace

// --- public API --------------------------------------------------------------

ConvertCellsMode convert_cells_mode_from_name(const std::string& rName) {
    if (rName == "linearize")
        return ConvertCellsMode::Linearize;
    if (rName == "simplexify")
        return ConvertCellsMode::Simplexify;
    if (rName == "elevate")
        return ConvertCellsMode::Elevate;
    throw std::invalid_argument("convert_cells: unknown mode '" + rName +
                                "' (expected 'linearize', 'simplexify', or 'elevate')");
}

namespace {

// Carry the input's named regions onto the output. The cell map is FirstChild:
// a parent's children occupy a contiguous run, which is also correct for the
// 1:1 modes (their maps are the monotone identity). Side regions are dropped —
// a child cell is a new cell of a subdivided or different topology, so its
// facets have no correspondence with the parent's. See detail/region_remap.hpp.
void ccells_carry_regions(const Mesh& rIn, ConvertCellsResult& rRes) {
    detail::RegionRemap rmap;
    rmap.pPointMap = &rRes.mPointMap;
    rmap.mCellMapKind = detail::CellMapKind::FirstChild;
    rmap.pCellMaps = &rRes.mCellMaps;
    rmap.mOpName = "convert_cells";
    detail::remap_regions(rIn, rRes.mMesh, rmap);
}

}  // namespace

ConvertCellsResult convert_cells(const Mesh& rMesh, const ConvertCellsOptions& rOptions) {
    ConvertCellsResult res;
    switch (rOptions.mMode) {
        case ConvertCellsMode::Linearize:
            res = ccells_linearize(rMesh, rOptions.mRecordParentIds);
            break;
        case ConvertCellsMode::Simplexify:
            res = ccells_simplexify(rMesh, rOptions.mRecordParentIds);
            break;
        case ConvertCellsMode::Elevate:
            res = ccells_elevate(rMesh, rOptions.mRecordParentIds);
            break;
        default:
            throw std::invalid_argument("convert_cells: unhandled mode");
    }
    ccells_carry_regions(rMesh, res);
    return res;
}

}  // namespace meshioplusplus
