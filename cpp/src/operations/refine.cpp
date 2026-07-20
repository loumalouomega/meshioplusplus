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

// Uniform mesh refinement: subdivide every cell into congruent same-type
// children, inserting nodes at edge midpoints, quad-face centres and (for the
// hexahedron) the body centre. Built entirely through the uniform mesh API so it
// compiles under every mesh backend. See operations/refine.hpp for the contract.
//
// Determinism: the subdivision templates are fixed, and the new-node numbering
// comes from a serial dedup pass over a parallel-filled disjoint-slot buffer --
// surface.cpp's phase-split idiom -- never from a concurrent hash insert. Output
// is therefore byte-identical across backends and thread counts.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// --- templates ---------------------------------------------------------------

// One cell type's subdivision template.
//
// Local node ids in `mChildren` address a single flat space per parent cell:
//   [0, mNumCorners)                          the parent's own corner nodes
//   [mNumCorners, +num_edges)                 one node per edge, in
//                                             cell_refine_edges() order
//   [.., +num_quad_faces)                     one node per quad face, in
//                                             cell_refine_quad_faces() order
//   [.., +1) when mHasBody                    the body centre
//
// That layout is not arbitrary: it coincides exactly with each type's own
// meshio/VTK full-Lagrange node numbering (line3, triangle6, quad9, tetra10,
// wedge18, hexahedron27), which is what makes the tables below readable against
// the reference elements.
struct RefineTemplate {
    CellType mType;
    std::uint8_t mNumCorners;
    bool mHasBody;
    std::vector<std::vector<std::uint8_t>> mChildren;
};

const std::unordered_map<CellType, RefineTemplate>& refine_table() {
    static const std::unordered_map<CellType, RefineTemplate> table = [] {
        std::unordered_map<CellType, RefineTemplate> t;
        // line3 layout: 2 = mid(0,1).
        t[CellType::Line] = {CellType::Line, 2, false, {{0, 2}, {2, 1}}};
        // triangle6 layout: 3 = m(0,1), 4 = m(1,2), 5 = m(2,0). The central
        // child keeps the parent's winding.
        t[CellType::Triangle] = {
            CellType::Triangle, 3, false, {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}}};
        // quad9 layout: 4..7 = edge mids, 8 = face centre. Each parent corner
        // stays in its own slot, which is what preserves the winding.
        t[CellType::Quad] = {
            CellType::Quad, 4, false, {{0, 4, 8, 7}, {4, 1, 5, 8}, {8, 5, 2, 6}, {7, 8, 6, 3}}};
        // tetra10 layout: 4=m(0,1) 5=m(1,2) 6=m(0,2) 7=m(0,3) 8=m(1,3) 9=m(2,3).
        // Four corner tetrahedra (each a half-scale homothety of the parent
        // about its own vertex), then the residual octahedron split along the
        // fixed interior diagonal 4-9 with the remaining ring 6->7->8->5. All
        // eight children have exactly one eighth of the parent's volume.
        t[CellType::Tetra] = {CellType::Tetra,
                              4,
                              false,
                              {
                                  {0, 4, 6, 7},
                                  {4, 1, 5, 8},
                                  {6, 5, 2, 9},
                                  {7, 8, 9, 3},
                                  {4, 9, 6, 7},
                                  {4, 9, 7, 8},
                                  {4, 9, 8, 5},
                                  {4, 9, 5, 6},
                              }};
        // wedge18 layout: 6..8 bottom-triangle mids, 9..11 top-triangle mids,
        // 12..14 vertical mids, 15..17 quad-face centres. A wedge refines as
        // "triangle 1-to-4 split x 2 vertical levels", and the mid-level
        // triangle's three edge midpoints ARE the three quad-face centres --
        // which is why a wedge needs no body node.
        t[CellType::Wedge] = {CellType::Wedge,
                              6,
                              false,
                              {
                                  {0, 6, 8, 12, 15, 17},
                                  {6, 1, 7, 15, 13, 16},
                                  {8, 7, 2, 17, 16, 14},
                                  {6, 7, 8, 15, 16, 17},
                                  {12, 15, 17, 3, 9, 11},
                                  {15, 13, 16, 9, 4, 10},
                                  {17, 16, 14, 11, 10, 5},
                                  {15, 16, 17, 9, 10, 11},
                              }};
        // hexahedron27 layout: 8..19 edge mids, 20..25 face centres, 26 body.
        // Rows follow the parent's own parametric (i,j,k) ordering over the
        // 3x3x3 lattice, so orientation is preserved by construction.
        t[CellType::Hexahedron] = {CellType::Hexahedron,
                                   8,
                                   true,
                                   {
                                       {0, 8, 24, 11, 16, 20, 26, 23},
                                       {8, 1, 9, 24, 20, 17, 21, 26},
                                       {11, 24, 10, 3, 23, 26, 22, 19},
                                       {24, 9, 2, 10, 26, 21, 18, 22},
                                       {16, 20, 26, 23, 4, 12, 25, 15},
                                       {20, 17, 21, 26, 12, 5, 13, 25},
                                       {23, 26, 22, 19, 15, 25, 14, 7},
                                       {26, 21, 18, 22, 25, 13, 6, 14},
                                   }};
        return t;
    }();
    return table;
}

const RefineTemplate* refine_template(CellType Type) {
    const auto& table = refine_table();
    auto it = table.find(Type);
    return it == table.end() ? nullptr : &it->second;
}

// --- new-node keys -----------------------------------------------------------

// A sorted, -1-padded node key identifying the entity a new node sits on: an
// edge is {-1, -1, a, b} and a quad face is {p, q, r, s}. Node ids are
// non-negative and the pad is -1, so an edge key can never equal a face key --
// the two key spaces share one map with no discriminator field. This is the
// same padding trick surface.cpp uses to let triangle and quad facets share a
// map. The arity is recoverable from the key itself (key[0] < 0 => edge), so
// the coordinate pass needs no side table.
using RefineNodeKey = std::array<std::int64_t, 4>;

struct RefineNodeKeyHash {
    std::size_t operator()(const RefineNodeKey& rKey) const {
        std::size_t h = 0;
        for (std::int64_t v : rKey)
            h ^= std::hash<std::int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

RefineNodeKey refine_edge_key(std::int64_t a, std::int64_t b) {
    return a < b ? RefineNodeKey{-1, -1, a, b} : RefineNodeKey{-1, -1, b, a};
}

RefineNodeKey refine_face_key(std::int64_t p, std::int64_t q, std::int64_t r, std::int64_t s) {
    RefineNodeKey key{p, q, r, s};
    std::sort(key.begin(), key.end());
    return key;
}

// --- per-block plumbing ------------------------------------------------------

// Everything one input block contributes to the shared slot buffer.
struct RefineBlockDesc {
    const RefineTemplate* mpTpl = nullptr;
    const std::vector<detail::CellEdgePair>* mpEdges = nullptr;
    const std::vector<detail::CellQuadFace>* mpFaces = nullptr;
    std::size_t mNumCells = 0;
    std::size_t mSlotsPerCell = 0;
    std::size_t mFirstSlot = 0;  // into the flat edge+face key buffer
    std::size_t mFirstBody = 0;  // into the body-centre range
};

// Reject, by name, every construct that cannot yield same-type children.
void refine_check_block(const Mesh::CellView& rBlock) {
    if (rBlock.IsRagged())
        throw std::invalid_argument(
            "refine: cannot refine ragged cell block '" + std::string(rBlock.Type()) +
            "' (polygon/polyhedron blocks have no same-type subdivision template)");
    const CellType type = cell_type_from_name(std::string(rBlock.Type()));
    if (refine_template(type) != nullptr)
        return;
    if (type == CellType::Pyramid)
        throw std::invalid_argument(
            "refine: cannot refine cell type 'pyramid' into same-type children (its uniform "
            "refinement is 6 pyramids + 4 tetrahedra; simplexify the mesh first)");
    throw std::invalid_argument("refine: cannot refine cell type '" + cell_type_name(type) +
                                "' into same-type children (higher-order cells have no same-type "
                                "subdivision template; linearize the mesh first)");
}

// The identity Int64 map [0, n).
NDArray refine_identity_map(std::size_t n) {
    NDArray a = NDArray::Uninit(DType::Int64, {n});
    std::int64_t* dst = a.As<std::int64_t>();
    for (std::size_t i = 0; i < n; ++i)
        dst[i] = static_cast<std::int64_t>(i);
    return a;
}

// A vector of Int64 as a 1-D NDArray.
NDArray refine_int64_vector(const std::vector<std::int64_t>& rValues) {
    NDArray a = NDArray::Uninit(DType::Int64, {rValues.size()});
    if (!rValues.empty())
        std::memcpy(a.Data(), rValues.data(), rValues.size() * sizeof(std::int64_t));
    return a;
}

// --- one refinement level ----------------------------------------------------

RefineResult refine_once(const Mesh& rMesh) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::size_t num_points = rMesh.NumPoints();

    // --- pre-scan: validate, size the shared slot buffer and the body range ---
    std::vector<RefineBlockDesc> descs(nblocks);
    std::size_t total_slots = 0;
    std::size_t total_bodies = 0;
    {
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            refine_check_block(cb);
            const CellType type = cell_type_from_name(std::string(cb.Type()));
            RefineBlockDesc& d = descs[bi];
            d.mpTpl = refine_template(type);
            d.mpEdges = &detail::cell_refine_edges(type);
            d.mpFaces = &detail::cell_refine_quad_faces(type);
            d.mNumCells = cb.NumCells();
            d.mSlotsPerCell = d.mpEdges->size() + d.mpFaces->size();
            d.mFirstSlot = total_slots;
            total_slots += d.mNumCells * d.mSlotsPerCell;
            if (d.mpTpl->mHasBody) {
                d.mFirstBody = total_bodies;
                total_bodies += d.mNumCells;
            }
            ++bi;
        }
    }

    // --- phase 1: fill entity keys into disjoint slots (parallel) ------------
    std::vector<RefineNodeKey> keys(total_slots);
    {
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const RefineBlockDesc& d = descs[bi++];
            if (d.mSlotsPerCell == 0 || d.mNumCells == 0)
                continue;
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            const std::vector<detail::CellEdgePair>& edges = *d.mpEdges;
            const std::vector<detail::CellQuadFace>& faces = *d.mpFaces;
            const std::size_t nedges = edges.size();
            const std::size_t nslots = d.mSlotsPerCell;
            parallel_for(d.mNumCells * nslots, [&](std::size_t j) {
                const std::size_t cell = j / nslots;
                const std::size_t slot = j % nslots;
                const std::size_t row = cell * npc;
                RefineNodeKey& key = keys[d.mFirstSlot + j];
                if (slot < nedges) {
                    const detail::CellEdgePair& e = edges[slot];
                    key = refine_edge_key(detail::read_int(conn, row + e[0]),
                                          detail::read_int(conn, row + e[1]));
                } else {
                    const detail::CellQuadFace& f = faces[slot - nedges];
                    key = refine_face_key(
                        detail::read_int(conn, row + f[0]), detail::read_int(conn, row + f[1]),
                        detail::read_int(conn, row + f[2]), detail::read_int(conn, row + f[3]));
                }
            });
        }
    }

    // --- phase 2: dedup in stored order (SERIAL -> deterministic) ------------
    // This pass is the determinism pin: ids are handed out by a single sweep
    // over the slot buffer, whose order is a pure function of (block, cell,
    // slot). It must never become a concurrent insert.
    std::unordered_map<RefineNodeKey, std::int64_t, RefineNodeKeyHash> node_id;
    node_id.reserve(total_slots * 2);
    std::vector<RefineNodeKey> new_nodes;
    std::vector<std::int64_t> slot_id(total_slots);
    for (std::size_t i = 0; i < total_slots; ++i) {
        auto it = node_id.find(keys[i]);
        if (it == node_id.end()) {
            const std::int64_t id = static_cast<std::int64_t>(num_points + new_nodes.size());
            node_id.emplace(keys[i], id);
            new_nodes.push_back(keys[i]);
            slot_id[i] = id;
        } else {
            slot_id[i] = it->second;
        }
    }
    keys.clear();
    keys.shrink_to_fit();  // phase 4 needs only slot_id; halves peak memory
    node_id.clear();

    // Body centres are unique per cell and so need no dedup, but their ids only
    // become known once the deduped range is closed.
    const std::size_t body_base = num_points + new_nodes.size();
    const std::size_t num_points_out = body_base + total_bodies;

    Mesh out;

    // --- phase 3: points -----------------------------------------------------
    // A new node's coordinate is the mean of its entity's corners. The mean is
    // order-independent, so two neighbours sharing an entity compute
    // bit-identical coordinates from the same sorted key -- no tie-break rule.
    {
        const NDArray& points = rMesh.Points();
        const std::size_t dim = detail::cols(points);
        NDArray new_points = NDArray::Uninit(points.Dtype(), {num_points_out, dim});
        std::memcpy(new_points.Data(), points.Data(), points.Nbytes());
        parallel_for_bw(new_nodes.size(), [&](std::size_t i) {
            const RefineNodeKey& key = new_nodes[i];
            const std::size_t first = key[0] < 0 ? 2 : 0;
            const double inv = 1.0 / static_cast<double>(4 - first);
            for (std::size_t k = 0; k < dim; ++k) {
                double sum = 0.0;
                for (std::size_t c = first; c < 4; ++c)
                    sum += detail::read_double(points, static_cast<std::size_t>(key[c]) * dim + k);
                detail::write_double(new_points, (num_points + i) * dim + k, sum * inv);
            }
        });
        // Body centres: the mean of the parent's eight corners.
        if (total_bodies > 0) {
            std::size_t bi = 0;
            for (const auto cb : rMesh.CellRange()) {
                const RefineBlockDesc& d = descs[bi++];
                if (!d.mpTpl->mHasBody || d.mNumCells == 0)
                    continue;
                const NDArray& conn = cb.Conn();
                const std::size_t npc = cb.NodesPerCell();
                const std::size_t ncorners = d.mpTpl->mNumCorners;
                const double inv = 1.0 / static_cast<double>(ncorners);
                parallel_for_bw(d.mNumCells, [&](std::size_t c) {
                    const std::size_t dst = body_base + d.mFirstBody + c;
                    for (std::size_t k = 0; k < dim; ++k) {
                        double sum = 0.0;
                        for (std::size_t n = 0; n < ncorners; ++n) {
                            const std::size_t p =
                                static_cast<std::size_t>(detail::read_int(conn, c * npc + n));
                            sum += detail::read_double(points, p * dim + k);
                        }
                        detail::write_double(new_points, dst * dim + k, sum * inv);
                    }
                });
            }
        }
        out.AssignPoints(std::move(new_points));
    }

    // --- phase 4: child connectivity ----------------------------------------
    std::vector<NDArray> cell_maps;
    cell_maps.reserve(nblocks);
    {
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const RefineBlockDesc& d = descs[bi++];
            const RefineTemplate& tpl = *d.mpTpl;
            const std::size_t nchildren = tpl.mChildren.size();
            const std::size_t out_npc = tpl.mChildren[0].size();
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            const std::size_t nedges = d.mpEdges->size();
            const std::size_t nfaces = d.mpFaces->size();
            const std::size_t ncorners = tpl.mNumCorners;

            NDArray out_conn = NDArray::Uninit(DType::Int64, {d.mNumCells * nchildren, out_npc});
            std::int64_t* dst = out_conn.As<std::int64_t>();
            parallel_for_bw(d.mNumCells, [&](std::size_t c) {
                // Resolve the parent's local node space once per cell.
                std::int64_t local[27];
                for (std::size_t n = 0; n < ncorners; ++n)
                    local[n] = detail::read_int(conn, c * npc + n);
                const std::size_t slot0 = d.mFirstSlot + c * d.mSlotsPerCell;
                for (std::size_t s = 0; s < nedges + nfaces; ++s)
                    local[ncorners + s] = slot_id[slot0 + s];
                if (tpl.mHasBody)
                    local[ncorners + nedges + nfaces] =
                        static_cast<std::int64_t>(body_base + d.mFirstBody + c);
                for (std::size_t k = 0; k < nchildren; ++k) {
                    const std::vector<std::uint8_t>& child = tpl.mChildren[k];
                    std::int64_t* row = dst + (c * nchildren + k) * out_npc;
                    for (std::size_t n = 0; n < out_npc; ++n)
                        row[n] = local[child[n]];
                }
            });
            out.AddCellBlock(cell_type_name(tpl.mType), std::move(out_conn));

            std::vector<std::int64_t> first_child(d.mNumCells);
            for (std::size_t c = 0; c < d.mNumCells; ++c)
                first_child[c] = static_cast<std::int64_t>(c * nchildren);
            cell_maps.push_back(refine_int64_vector(first_child));
        }
    }

    // --- phase 5: point_data -------------------------------------------------
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (detail::rows(a) != num_points) {
            out.AddPointData(name, detail::data_owned_copy(a));
            continue;
        }
        const std::size_t ncomp = num_points == 0 ? 0 : a.Size() / num_points;
        std::vector<std::size_t> shape = a.Shape();
        shape[0] = num_points_out;
        NDArray b = NDArray::Uninit(a.Dtype(), std::move(shape));
        std::memcpy(b.Data(), a.Data(), a.Nbytes());
        parallel_for_bw(new_nodes.size(), [&](std::size_t i) {
            const RefineNodeKey& key = new_nodes[i];
            const std::size_t first = key[0] < 0 ? 2 : 0;
            const double inv = 1.0 / static_cast<double>(4 - first);
            for (std::size_t k = 0; k < ncomp; ++k) {
                double sum = 0.0;
                for (std::size_t c = first; c < 4; ++c)
                    sum += detail::read_double(a, static_cast<std::size_t>(key[c]) * ncomp + k);
                detail::write_double(b, (num_points + i) * ncomp + k, sum * inv);
            }
        });
        if (total_bodies > 0) {
            std::size_t bi = 0;
            for (const auto cb : rMesh.CellRange()) {
                const RefineBlockDesc& d = descs[bi++];
                if (!d.mpTpl->mHasBody || d.mNumCells == 0)
                    continue;
                const NDArray& conn = cb.Conn();
                const std::size_t npc = cb.NodesPerCell();
                const std::size_t ncorners = d.mpTpl->mNumCorners;
                const double inv = 1.0 / static_cast<double>(ncorners);
                parallel_for_bw(d.mNumCells, [&](std::size_t c) {
                    const std::size_t dst = body_base + d.mFirstBody + c;
                    for (std::size_t k = 0; k < ncomp; ++k) {
                        double sum = 0.0;
                        for (std::size_t n = 0; n < ncorners; ++n) {
                            const std::size_t p =
                                static_cast<std::size_t>(detail::read_int(conn, c * npc + n));
                            sum += detail::read_double(a, p * ncomp + k);
                        }
                        detail::write_double(b, dst * ncomp + k, sum * inv);
                    }
                });
            }
        }
        out.AddPointData(name, std::move(b));
    }

    // --- phase 6: cell_data replicated parent -> children --------------------
    for (const std::string& name : rMesh.CellDataNames()) {
        const std::size_t ndata = rMesh.CellDataNumBlocks(name);
        std::vector<NDArray> blocks;
        blocks.reserve(ndata);
        for (std::size_t b = 0; b < ndata; ++b) {
            const NDArray& a = rMesh.CellData(name, b);
            const std::size_t in_rows = detail::rows(a);
            if (ndata != nblocks || in_rows == 0 || in_rows != descs[b].mNumCells) {
                blocks.push_back(detail::data_owned_copy(a));
                continue;
            }
            const std::size_t nchildren = descs[b].mpTpl->mChildren.size();
            std::vector<std::size_t> shape = a.Shape();
            shape[0] = in_rows * nchildren;
            const std::size_t row_bytes = a.Nbytes() / in_rows;
            NDArray replicated = NDArray::Uninit(a.Dtype(), std::move(shape));
            const std::byte* src = a.Data();
            std::byte* dst = replicated.Data();
            parallel_for_bw(in_rows * nchildren, [&](std::size_t i) {
                std::memcpy(dst + i * row_bytes, src + (i / nchildren) * row_bytes, row_bytes);
            });
            blocks.push_back(std::move(replicated));
        }
        out.AddCellData(name, std::move(blocks));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    RefineResult res;
    res.mMesh = std::move(out);
    res.mPointMap = refine_identity_map(num_points);
    res.mCellMaps = std::move(cell_maps);
    return res;
}

// Fold `rNext`'s maps into `rAcc`'s, so the accumulated maps always point from
// the ORIGINAL mesh into the current one. Must run before `rNext.mMesh` is
// moved from.
void refine_compose_maps(RefineResult& rAcc, const RefineResult& rNext) {
    std::int64_t* pm = rAcc.mPointMap.As<std::int64_t>();
    const std::int64_t* npm = rNext.mPointMap.As<std::int64_t>();
    for (std::size_t i = 0; i < rAcc.mPointMap.Size(); ++i)
        pm[i] = pm[i] < 0 ? -1 : npm[pm[i]];
    // Blocks correspond 1:1 at every level, and a cell's descendants stay
    // contiguous under the parent-major/child-minor emission order, so composing
    // "first child" maps is a plain lookup.
    for (std::size_t b = 0; b < rAcc.mCellMaps.size() && b < rNext.mCellMaps.size(); ++b) {
        std::int64_t* cm = rAcc.mCellMaps[b].As<std::int64_t>();
        const std::int64_t* ncm = rNext.mCellMaps[b].As<std::int64_t>();
        for (std::size_t c = 0; c < rAcc.mCellMaps[b].Size(); ++c)
            cm[c] = cm[c] < 0 ? -1 : ncm[cm[c]];
    }
}

// Attach refine:parent_cell from the composed maps, so it names the ORIGINAL
// ancestor rather than the immediate parent even at several levels.
void refine_attach_parent_ids(RefineResult& rResult) {
    const Mesh& mesh = rResult.mMesh;
    const std::size_t nblocks = mesh.NumCellBlocks();
    if (nblocks != rResult.mCellMaps.size() || nblocks == 0)
        return;
    std::vector<NDArray> blocks;
    blocks.reserve(nblocks);
    std::size_t bi = 0;
    for (const auto cb : mesh.CellRange()) {
        const NDArray& map = rResult.mCellMaps[bi++];
        const std::int64_t* first = map.As<std::int64_t>();
        const std::size_t nparents = map.Size();
        const std::size_t ncells = cb.NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {ncells, 1});
        std::int64_t* dst = a.As<std::int64_t>();
        std::fill(dst, dst + ncells, static_cast<std::int64_t>(-1));
        for (std::size_t p = 0; p < nparents; ++p) {
            if (first[p] < 0)
                continue;
            const std::size_t lo = static_cast<std::size_t>(first[p]);
            const std::size_t hi = p + 1 < nparents && first[p + 1] >= 0
                                       ? static_cast<std::size_t>(first[p + 1])
                                       : ncells;
            for (std::size_t c = lo; c < hi && c < ncells; ++c)
                dst[c] = static_cast<std::int64_t>(p);
        }
        blocks.push_back(std::move(a));
    }
    rResult.mMesh.AddCellData("refine:parent_cell", std::move(blocks));
}

// Cells the next level would produce, for the pre-flight size warning.
std::size_t refine_projected_cells(const Mesh& rMesh) {
    std::size_t total = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsRagged())
            continue;
        const RefineTemplate* tpl = refine_template(cell_type_from_name(std::string(cb.Type())));
        total += cb.NumCells() * (tpl == nullptr ? 1 : tpl->mChildren.size());
    }
    return total;
}

}  // namespace

// --- public API --------------------------------------------------------------

RefineResult refine(const Mesh& rMesh, const RefineOptions& rOptions) {
    if (rOptions.mLevels <= 0) {
        RefineResult res;
        res.mMesh = detail::clone_mesh(rMesh);
        res.mPointMap = refine_identity_map(rMesh.NumPoints());
        res.mCellMaps.reserve(rMesh.NumCellBlocks());
        for (const auto cb : rMesh.CellRange())
            res.mCellMaps.push_back(refine_identity_map(cb.NumCells()));
        if (rOptions.mRecordParentIds)
            refine_attach_parent_ids(res);
        return res;
    }

    // Refinement is exponential in `levels`; the failure mode at depth is a
    // bad_alloc rather than a wrong answer, so say so before it happens.
    static constexpr std::size_t refine_warn_cells = 20'000'000;
    if (refine_projected_cells(rMesh) > refine_warn_cells)
        log::warn("refine: one level of this mesh yields ~{} cells; {} level(s) requested.",
                  refine_projected_cells(rMesh), rOptions.mLevels);

    RefineResult acc = refine_once(rMesh);
    for (int level = 1; level < rOptions.mLevels; ++level) {
        RefineResult next = refine_once(acc.mMesh);  // borrows acc.mMesh
        refine_compose_maps(acc, next);              // reads next's maps first
        acc.mMesh = std::move(next.mMesh);           // sequenced after
    }
    if (rOptions.mRecordParentIds)
        refine_attach_parent_ids(acc);
    return acc;
}

}  // namespace meshioplusplus
