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

// Mesh refinement: subdivide cells into congruent same-type children, either
// every cell (uniform) or a selected subset closed up to a conforming mesh
// (selective). Built entirely through the uniform mesh API so it compiles under
// every mesh backend. See operations/refine.hpp for the contract and
// detail/refine_templates.hpp for the per-(type, split-edge mask) tables.
//
// The whole operation is driven by one set: which EDGES carry a new mid-edge
// node. Everything else is derived from it -- a quad face carries a centre iff
// all four of its edges are split, a hexahedron carries a body node iff all
// twelve are, and a cell's subdivision is the template for its own split-edge
// mask. Two cells sharing an entity read the same edges, so they cannot
// disagree, which is what makes conformity structural rather than sampled. With
// every edge split those rules collapse into the uniform templates, which is
// what makes uniform output byte-identical to the pre-selective implementation.
//
// Determinism: the split set grows under a monotone idempotent closure operator
// (detail/refine_templates.hpp), so its fixed point does not depend on the order
// cells are visited in; and the new-node numbering comes from a serial dedup
// pass over a parallel-filled disjoint-slot buffer -- surface.cpp's phase-split
// idiom -- never from a concurrent hash insert. Output is therefore
// byte-identical across backends and thread counts.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
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
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/refine_hierarchy.hpp"
#include "meshioplusplus/detail/refine_templates.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

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
    CellType mType = CellType::Custom;
    const std::vector<detail::CellEdgePair>* mpEdges = nullptr;
    const std::vector<detail::CellQuadFace>* mpFaces = nullptr;
    // Per quad face, the bits of the four edges bounding it: that face carries a
    // centre node exactly when all of them are split. Derived here once per
    // block from cell_refine_edges/cell_refine_quad_faces rather than tabulated,
    // so the two tables stay the only owners of the orderings.
    std::vector<std::uint16_t> mFaceEdgeMasks;
    std::uint8_t mNumCorners = 0;
    std::uint16_t mFullMask = 0;
    std::size_t mNumCells = 0;
    std::size_t mSlotsPerCell = 0;
    std::size_t mFirstSlot = 0;  // into the flat edge+face key buffer
    std::size_t mFirstCell = 0;  // this block's base in the global cell numbering
};

// Reject, by name, every construct that cannot yield same-type children.
void refine_check_block(const Mesh::CellView& rBlock) {
    if (rBlock.IsRagged())
        throw std::invalid_argument(
            "refine: cannot refine ragged cell block (run convert_cells(simplexify) first) '" +
            std::string(rBlock.Type()) +
            "' (polygon/polyhedron blocks have no same-type subdivision template)");
    const CellType type = cell_type_from_name(std::string(rBlock.Type()));
    if (detail::refine_type_supported(type))
        return;
    if (type == CellType::Pyramid)
        throw std::invalid_argument(
            "refine: cannot refine cell type 'pyramid' into same-type children (its uniform "
            "refinement is 6 pyramids + 4 tetrahedra; simplexify the mesh first)");
    throw std::invalid_argument("refine: cannot refine cell type '" + cell_type_name(type) +
                                "' into same-type children (higher-order cells have no same-type "
                                "subdivision template; linearize the mesh first)");
}

// The edge-bit mask of each quadrilateral face of a cell type.
std::vector<std::uint16_t> refine_face_edge_masks(CellType Type) {
    const std::vector<detail::CellEdgePair>& edges = detail::cell_refine_edges(Type);
    const std::vector<detail::CellQuadFace>& faces = detail::cell_refine_quad_faces(Type);
    std::vector<std::uint16_t> out(faces.size(), 0);
    for (std::size_t f = 0; f < faces.size(); ++f) {
        for (std::size_t i = 0; i < 4; ++i) {
            const std::uint8_t a = faces[f][i];
            const std::uint8_t b = faces[f][(i + 1) % 4];
            for (std::size_t k = 0; k < edges.size(); ++k) {
                if ((edges[k][0] == a && edges[k][1] == b) ||
                    (edges[k][0] == b && edges[k][1] == a)) {
                    out[f] = static_cast<std::uint16_t>(out[f] | (1u << k));
                    break;
                }
            }
        }
    }
    return out;
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

// --- refine:level ------------------------------------------------------------

// The per-cell refinement depth already recorded on a mesh, or an empty vector
// when it carries none. A malformed array is ignored with a warning rather than
// rejected: it is this operation's own bookkeeping, not user input, and
// refusing to refine because of it would be a worse answer than recomputing it.
std::vector<std::int64_t> refine_read_levels(const Mesh& rMesh,
                                             const std::vector<RefineBlockDesc>& rDescs) {
    if (!rMesh.HasCellData(kRefineLevelName))
        return {};
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (rMesh.CellDataNumBlocks(kRefineLevelName) != nblocks) {
        log::warn("refine: ignoring '{}': it does not cover every cell block.", kRefineLevelName);
        return {};
    }
    std::vector<std::int64_t> levels;
    levels.reserve(rDescs.empty() ? 0 : rDescs.back().mFirstCell + rDescs.back().mNumCells);
    for (std::size_t b = 0; b < nblocks; ++b) {
        const NDArray& a = rMesh.CellData(kRefineLevelName, b);
        if (detail::rows(a) != rDescs[b].mNumCells ||
            (rDescs[b].mNumCells != 0 && a.Size() / rDescs[b].mNumCells != 1)) {
            log::warn("refine: ignoring '{}': block {} is not one scalar value per cell.",
                      kRefineLevelName, b);
            return {};
        }
        for (std::size_t c = 0; c < rDescs[b].mNumCells; ++c)
            levels.push_back(detail::read_int(a, c));
    }
    return levels;
}

// --- refine:cell_id / refine:parent_id ----------------------------------------
//
// RefineHierarchyState / refine_read_hierarchy live in
// detail/refine_hierarchy.hpp now -- hoisted verbatim so undo_green.cpp can
// share the identical Absent/Valid/Invalid read against the *coarse* mesh
// rather than transcribing it a second time.
using detail::refine_read_hierarchy;
using detail::RefineHierarchyState;

// --- refine:entity -----------------------------------------------------------

// The per-point entity keys already recorded on a mesh, or an empty vector when
// it carries none (or carries a malformed or stale one).
//
// Staleness is the real hazard, not malformation: the entries name input point
// indices, which `refine` never renumbers but `reorder`/`clean`/`crop`/`merge`
// do, and which `transform`/`smooth` move. A stale key would silently graft one
// entity's node onto another's, so it is not enough to range-check it -- every
// key must still reproduce its own point's coordinates, by the same arithmetic
// phase 5 used to write them. The check runs through a scratch of the points'
// own dtype so a Float32 mesh compares what was actually stored rather than a
// widened double.
std::vector<RefineNodeKey> refine_read_entity_keys(const Mesh& rMesh) {
    if (!rMesh.HasPointData(kRefineEntityName))
        return {};
    const std::size_t num_points = rMesh.NumPoints();
    const NDArray& a = rMesh.PointData(kRefineEntityName);
    if (detail::rows(a) != num_points || (num_points != 0 && a.Size() / num_points != 4)) {
        log::warn("refine: ignoring '{}': it is not four values per point.", kRefineEntityName);
        return {};
    }

    std::vector<RefineNodeKey> keys(num_points);
    for (std::size_t p = 0; p < num_points; ++p)
        for (std::size_t k = 0; k < 4; ++k)
            keys[p][k] = detail::read_int(a, p * 4 + k);

    const NDArray& points = rMesh.Points();
    const std::size_t dim = detail::cols(points);
    NDArray scratch = NDArray::Uninit(points.Dtype(), {std::size_t{1}});
    for (std::size_t p = 0; p < num_points; ++p) {
        const RefineNodeKey& key = keys[p];
        if (key[3] < 0)
            continue;  // the {-1,-1,-1,-1} sentinel: a point refine did not create
        const std::size_t first = key[0] < 0 ? 2 : 0;
        bool ok = true;
        for (std::size_t c = first; c < 4 && ok; ++c)
            ok = key[c] >= 0 && static_cast<std::size_t>(key[c]) < num_points;
        // The phase 5 arithmetic, verbatim -- same order, same divisor, same
        // store-then-read round trip.
        const double inv = 1.0 / static_cast<double>(4 - first);
        for (std::size_t k = 0; k < dim && ok; ++k) {
            double sum = 0.0;
            for (std::size_t c = first; c < 4; ++c)
                sum += detail::read_double(points, static_cast<std::size_t>(key[c]) * dim + k);
            detail::write_double(scratch, 0, sum * inv);
            ok = detail::read_double(scratch, 0) == detail::read_double(points, p * dim + k);
        }
        if (!ok) {
            log::warn(
                "refine: ignoring '{}': point {} no longer sits on the entity it records, so the "
                "mesh was renumbered or moved since it was written.",
                kRefineEntityName, p);
            return {};
        }
    }
    return keys;
}

// --- one refinement level ----------------------------------------------------

// Refine `rMesh` once. `pRedSeed` is a per-global-cell flag choosing the cells
// to split fully; `nullptr` means every cell (the uniform path, which skips the
// closure entirely). `pOutRedChildren`, when given, receives the same flag for
// the OUTPUT cells -- set on the children of a fully-split parent -- which is
// what the next level uses as its own seed. `ForceEntity` attaches
// `refine:entity` unconditionally (the caller asked for the persistent
// hierarchy, whose node-level interpolation stencil this array already IS)
// even when this particular pass leaves no hanging node.
RefineResult refine_once(const Mesh& rMesh, const std::vector<char>* pRedSeed,
                         RefineClosure Closure, bool RecordLevels, bool ForceEntity,
                         std::vector<char>* pOutRedChildren) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::size_t num_points = rMesh.NumPoints();
    const bool propagate = Closure == RefineClosure::Propagate;

    // --- phase 0: validate, size the shared slot buffer ----------------------
    std::vector<RefineBlockDesc> descs(nblocks);
    std::size_t total_slots = 0;
    std::size_t total_cells = 0;
    {
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            refine_check_block(cb);
            const CellType type = cell_type_from_name(std::string(cb.Type()));
            RefineBlockDesc& d = descs[bi];
            d.mType = type;
            d.mpEdges = &detail::cell_refine_edges(type);
            d.mpFaces = &detail::cell_refine_quad_faces(type);
            d.mFaceEdgeMasks = refine_face_edge_masks(type);
            d.mNumCorners = static_cast<std::uint8_t>(cell_type_num_nodes(type));
            d.mFullMask = detail::refine_full_mask(type);
            d.mNumCells = cb.NumCells();
            d.mSlotsPerCell = d.mpEdges->size() + d.mpFaces->size();
            d.mFirstSlot = total_slots;
            d.mFirstCell = total_cells;
            total_slots += d.mNumCells * d.mSlotsPerCell;
            total_cells += d.mNumCells;
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
    // This pass is the determinism pin: entities are numbered by a single sweep
    // over the slot buffer, whose order is a pure function of (block, cell,
    // slot). It must never become a concurrent insert. Note it numbers
    // ENTITIES, not points -- which of them earn a node is decided in phase 4,
    // and the point ids are handed out in this same first-seen order there.
    std::vector<RefineNodeKey> entities;
    std::vector<std::int64_t> entity_of_slot(total_slots);
    {
        std::unordered_map<RefineNodeKey, std::int64_t, RefineNodeKeyHash> entity_id;
        entity_id.reserve(total_slots * 2);
        for (std::size_t i = 0; i < total_slots; ++i) {
            auto it = entity_id.find(keys[i]);
            if (it == entity_id.end()) {
                const std::int64_t id = static_cast<std::int64_t>(entities.size());
                entity_id.emplace(keys[i], id);
                entities.push_back(keys[i]);
                entity_of_slot[i] = id;
            } else {
                entity_of_slot[i] = it->second;
            }
        }
    }
    keys.clear();
    keys.shrink_to_fit();  // nothing downstream needs the per-slot keys

    // --- phase 3: the split-edge set and its closure -------------------------
    // `masks[c]` is the promoted split-edge mask of cell c. Uniform refinement
    // short-circuits the whole fixed point: every cell is red, every edge is
    // split, and the promotion of a full mask is itself.
    std::vector<std::uint16_t> masks(total_cells, 0);
    std::vector<char> split(entities.size(), 0);
    // Balanced needs the per-cell level up front rather than at emit time: the
    // whole rule is a comparison between neighbours' levels.
    const std::vector<std::int64_t> base_levels = refine_read_levels(rMesh, descs);
    if (pRedSeed != nullptr && Closure == RefineClosure::Balanced) {
        // 2:1 balance. A cell is split fully or not at all -- no transitional
        // templates -- and is drawn in only when a neighbour would otherwise end
        // up more than one level finer. Marking a cell only ever raises a post
        // level, so the rule is monotone and its least fixed point is unique: as
        // with the mask closure, determinism here is a property of the
        // formulation rather than of the traversal.
        //
        // Adjacency is by shared NODE, not by shared edge entity. That is not a
        // stylistic choice: across a hanging interface the coarse cell spans a
        // whole edge while the fine cell has only half of it, so the two are
        // *different* entities and an edge-keyed rule is blind to exactly the
        // coarse/fine adjacency it exists to police. The corner nodes are what
        // the two still share. (It is also the stronger, standard "corner
        // balance", so a diagonal neighbour counts too.)
        std::vector<char> red(*pRedSeed);
        const auto level_of = [&](std::size_t c) {
            return base_levels.empty() ? std::int64_t{0} : base_levels[c];
        };
        for (;;) {
            // Per input point, the largest post-refinement level among the cells
            // that touch it.
            std::vector<std::int64_t> node_max(num_points,
                                               std::numeric_limits<std::int64_t>::min());
            {
                std::size_t bi = 0;
                for (const auto cb : rMesh.CellRange()) {
                    const RefineBlockDesc& d = descs[bi++];
                    const NDArray& conn = cb.Conn();
                    const std::size_t npc = cb.NodesPerCell();
                    for (std::size_t c = 0; c < d.mNumCells; ++c) {
                        const std::size_t g = d.mFirstCell + c;
                        const std::int64_t post = level_of(g) + (red[g] ? 1 : 0);
                        for (std::size_t n = 0; n < npc; ++n) {
                            const std::int64_t p = detail::read_int(conn, c * npc + n);
                            if (p >= 0 && p < static_cast<std::int64_t>(num_points)) {
                                std::int64_t& m = node_max[static_cast<std::size_t>(p)];
                                m = std::max(m, post);
                            }
                        }
                    }
                }
            }
            bool changed = false;
            std::size_t bi = 0;
            for (const auto cb : rMesh.CellRange()) {
                const RefineBlockDesc& d = descs[bi++];
                const NDArray& conn = cb.Conn();
                const std::size_t npc = cb.NodesPerCell();
                for (std::size_t c = 0; c < d.mNumCells; ++c) {
                    const std::size_t g = d.mFirstCell + c;
                    if (red[g])
                        continue;
                    for (std::size_t n = 0; n < npc; ++n) {
                        const std::int64_t p = detail::read_int(conn, c * npc + n);
                        if (p >= 0 && p < static_cast<std::int64_t>(num_points) &&
                            node_max[static_cast<std::size_t>(p)] > level_of(g) + 1) {
                            red[g] = 1;
                            changed = true;
                            break;
                        }
                    }
                }
            }
            if (!changed)
                break;
        }
        // Only the fully split cells create nodes; an unrefined neighbour keeps
        // its whole shape and simply acquires hanging nodes on its edges.
        for (std::size_t b = 0; b < nblocks; ++b) {
            const RefineBlockDesc& d = descs[b];
            const std::size_t nedges = d.mpEdges->size();
            for (std::size_t c = 0; c < d.mNumCells; ++c) {
                const std::size_t g = d.mFirstCell + c;
                if (!red[g])
                    continue;
                masks[g] = d.mFullMask;
                const std::size_t slot0 = d.mFirstSlot + c * d.mSlotsPerCell;
                for (std::size_t k = 0; k < nedges; ++k)
                    split[static_cast<std::size_t>(entity_of_slot[slot0 + k])] = 1;
            }
        }
    } else if (pRedSeed == nullptr) {
        for (std::size_t b = 0; b < nblocks; ++b)
            std::fill(masks.begin() + static_cast<std::ptrdiff_t>(descs[b].mFirstCell),
                      masks.begin() +
                          static_cast<std::ptrdiff_t>(descs[b].mFirstCell + descs[b].mNumCells),
                      descs[b].mFullMask);
        std::fill(split.begin(), split.end(), static_cast<char>(1));
    } else {
        // Seed: every edge of every selected cell.
        for (std::size_t b = 0; b < nblocks; ++b) {
            const RefineBlockDesc& d = descs[b];
            const std::size_t nedges = d.mpEdges->size();
            for (std::size_t c = 0; c < d.mNumCells; ++c) {
                if (!(*pRedSeed)[d.mFirstCell + c])
                    continue;
                const std::size_t slot0 = d.mFirstSlot + c * d.mSlotsPerCell;
                for (std::size_t k = 0; k < nedges; ++k)
                    split[static_cast<std::size_t>(entity_of_slot[slot0 + k])] = 1;
            }
        }
        // Closure. `refine_promote_mask` is monotone and idempotent, so this
        // converges to a unique fixed point no matter what order cells are
        // visited in; the serial union below is what makes the *number* of
        // sweeps reproducible too. Bounded by the entity count, since bits are
        // only ever set.
        std::vector<std::uint16_t> promoted(total_cells, 0);
        for (;;) {
            for (std::size_t b = 0; b < nblocks; ++b) {
                const RefineBlockDesc& d = descs[b];
                if (d.mNumCells == 0)
                    continue;
                const std::size_t nedges = d.mpEdges->size();
                const CellType type = d.mType;
                parallel_for(d.mNumCells, [&](std::size_t c) {
                    const std::size_t slot0 = d.mFirstSlot + c * d.mSlotsPerCell;
                    std::uint16_t mask = 0;
                    for (std::size_t k = 0; k < nedges; ++k) {
                        if (split[static_cast<std::size_t>(entity_of_slot[slot0 + k])])
                            mask = static_cast<std::uint16_t>(mask | (1u << k));
                    }
                    promoted[d.mFirstCell + c] = detail::refine_promote_mask(type, mask, propagate);
                });
            }
            bool changed = false;
            for (std::size_t b = 0; b < nblocks; ++b) {
                const RefineBlockDesc& d = descs[b];
                const std::size_t nedges = d.mpEdges->size();
                for (std::size_t c = 0; c < d.mNumCells; ++c) {
                    const std::uint16_t mask = promoted[d.mFirstCell + c];
                    if (mask == 0)
                        continue;
                    const std::size_t slot0 = d.mFirstSlot + c * d.mSlotsPerCell;
                    for (std::size_t k = 0; k < nedges; ++k) {
                        if (!(mask & (1u << k)))
                            continue;
                        char& s = split[static_cast<std::size_t>(entity_of_slot[slot0 + k])];
                        if (!s) {
                            s = 1;
                            changed = true;
                        }
                    }
                }
            }
            if (!changed)
                break;
        }
        masks = std::move(promoted);
    }

    // The nodes a previous pass already placed on entities of this mesh, keyed by
    // the entity. The {-1,-1,-1,-1} sentinel is skipped: it is what an original
    // point and a body centre both carry, and admitting it would make every such
    // point collide on one key. First-seen wins, which only arises on an input
    // that was already torn.
    const std::vector<RefineNodeKey> base_entity_keys = refine_read_entity_keys(rMesh);
    std::unordered_map<RefineNodeKey, std::int64_t, RefineNodeKeyHash> persisted_node;
    for (std::size_t p = 0; p < base_entity_keys.size(); ++p) {
        if (base_entity_keys[p][3] < 0)
            continue;
        persisted_node.emplace(base_entity_keys[p], static_cast<std::int64_t>(p));
    }

    // --- phase 4: which entities earn a node, and their point ids ------------
    // An edge entity earns one iff it is split; a quad-face entity iff all four
    // of the bounding edges are. Both cells sharing an entity evaluate the same
    // predicate over the same edges, so they cannot disagree -- that is the
    // conformity argument, and it is why nothing here is tabulated per template.
    //
    // An entity whose key the input already recorded (see kRefineEntityName)
    // keeps the node it already has, **whether or not this pass splits it**.
    // Both halves of that matter. Reusing a split entity's node is what stops the
    // 2:1 balance rule from tearing the mesh, and resolving an *unsplit* one is
    // what lets the hanging rule below stay a single statement over entities: a
    // node still sitting on an edge its cell spans whole is reported by exactly
    // the same test that reports a freshly created one, so nothing has to be
    // carried forward or unmarked by hand.
    std::vector<std::int64_t> node_of_entity(entities.size(), -1);
    std::vector<std::int64_t> new_entities;  // entities needing a NEW point, in id order
    {
        std::vector<char> needs(entities.size(), 0);
        for (std::size_t e = 0; e < entities.size(); ++e)
            needs[e] = entities[e][0] < 0 ? split[e] : 0;
        for (std::size_t b = 0; b < nblocks; ++b) {
            const RefineBlockDesc& d = descs[b];
            const std::size_t nedges = d.mpEdges->size();
            const std::size_t nfaces = d.mpFaces->size();
            if (nfaces == 0 || d.mNumCells == 0)
                continue;
            for (std::size_t c = 0; c < d.mNumCells; ++c) {
                const std::uint16_t mask = masks[d.mFirstCell + c];
                const std::size_t slot0 = d.mFirstSlot + c * d.mSlotsPerCell;
                for (std::size_t f = 0; f < nfaces; ++f) {
                    const std::uint16_t need = d.mFaceEdgeMasks[f];
                    if ((mask & need) == need)
                        needs[static_cast<std::size_t>(entity_of_slot[slot0 + nedges + f])] = 1;
                }
            }
        }
        if (!persisted_node.empty()) {
            for (std::size_t e = 0; e < entities.size(); ++e) {
                auto it = persisted_node.find(entities[e]);
                if (it != persisted_node.end())
                    node_of_entity[e] = it->second;
            }
        }
        new_entities.reserve(entities.size());
        for (std::size_t e = 0; e < entities.size(); ++e) {
            if (!needs[e] || node_of_entity[e] >= 0)
                continue;
            node_of_entity[e] = static_cast<std::int64_t>(num_points + new_entities.size());
            new_entities.push_back(static_cast<std::int64_t>(e));
        }
    }

    // Body centres are unique per cell and so need no dedup, but their ids only
    // become known once the deduped range is closed. Only a fully split
    // hexahedron has one.
    const std::size_t body_base = num_points + new_entities.size();
    std::vector<std::int64_t> body_of_cell(total_cells, -1);
    std::size_t total_bodies = 0;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const RefineBlockDesc& d = descs[b];
        if (d.mType != CellType::Hexahedron)
            continue;
        for (std::size_t c = 0; c < d.mNumCells; ++c) {
            if (masks[d.mFirstCell + c] != d.mFullMask)
                continue;
            body_of_cell[d.mFirstCell + c] = static_cast<std::int64_t>(body_base + total_bodies);
            ++total_bodies;
        }
    }
    const std::size_t num_points_out = body_base + total_bodies;

    Mesh out;

    // --- phase 5: points -----------------------------------------------------
    // A new node's coordinate is the mean of its entity's corners. The mean is
    // order-independent, so two neighbours sharing an entity compute
    // bit-identical coordinates from the same sorted key -- no tie-break rule.
    {
        const NDArray& points = rMesh.Points();
        const std::size_t dim = detail::cols(points);
        NDArray new_points = NDArray::Uninit(points.Dtype(), {num_points_out, dim});
        std::memcpy(new_points.Data(), points.Data(), points.Nbytes());
        parallel_for_bw(new_entities.size(), [&](std::size_t i) {
            const RefineNodeKey& key = entities[static_cast<std::size_t>(new_entities[i])];
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
                if (d.mType != CellType::Hexahedron || d.mNumCells == 0)
                    continue;
                const NDArray& conn = cb.Conn();
                const std::size_t npc = cb.NodesPerCell();
                const std::size_t ncorners = d.mNumCorners;
                const double inv = 1.0 / static_cast<double>(ncorners);
                parallel_for_bw(d.mNumCells, [&](std::size_t c) {
                    const std::int64_t body = body_of_cell[d.mFirstCell + c];
                    if (body < 0)
                        return;
                    for (std::size_t k = 0; k < dim; ++k) {
                        double sum = 0.0;
                        for (std::size_t n = 0; n < ncorners; ++n) {
                            const std::size_t p =
                                static_cast<std::size_t>(detail::read_int(conn, c * npc + n));
                            sum += detail::read_double(points, p * dim + k);
                        }
                        detail::write_double(new_points, static_cast<std::size_t>(body) * dim + k,
                                             sum * inv);
                    }
                });
            }
        }
        out.AssignPoints(std::move(new_points));
    }

    // --- phase 6: child counts, then connectivity ----------------------------
    // A cell's child count depends on its mask, so the first-child map is a
    // serial ascending prefix sum rather than a constant stride. It stays
    // monotone with every entry non-negative and every run contiguous, which is
    // what keeps detail/region_remap.hpp's FirstChild contract valid.
    std::vector<NDArray> cell_maps;
    cell_maps.reserve(nblocks);
    std::vector<std::vector<std::int64_t>> parent_of_child(nblocks);
    std::vector<std::vector<char>> red_child(nblocks);
    {
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const std::size_t b = bi++;
            const RefineBlockDesc& d = descs[b];
            const std::size_t out_npc = static_cast<std::size_t>(cell_type_num_nodes(d.mType));
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            const std::size_t nedges = d.mpEdges->size();
            const std::size_t nfaces = d.mpFaces->size();
            const std::size_t ncorners = d.mNumCorners;

            std::vector<std::int64_t> first_child(d.mNumCells);
            std::size_t total_children = 0;
            for (std::size_t c = 0; c < d.mNumCells; ++c) {
                first_child[c] = static_cast<std::int64_t>(total_children);
                total_children +=
                    detail::refine_mask_template(d.mType, masks[d.mFirstCell + c]).NumChildren();
            }
            parent_of_child[b].resize(total_children);
            red_child[b].resize(total_children);

            NDArray out_conn = NDArray::Uninit(DType::Int64, {total_children, out_npc});
            std::int64_t* dst = out_conn.As<std::int64_t>();
            parallel_for_bw(d.mNumCells, [&](std::size_t c) {
                // Resolve the parent's local node space once per cell. A slot
                // whose entity earned no node stays -1; a template never
                // references one, by construction of the node rules in phase 4.
                std::int64_t local[27];
                for (std::size_t n = 0; n < ncorners; ++n)
                    local[n] = detail::read_int(conn, c * npc + n);
                const std::size_t slot0 = d.mFirstSlot + c * d.mSlotsPerCell;
                for (std::size_t s = 0; s < nedges + nfaces; ++s)
                    local[ncorners + s] =
                        node_of_entity[static_cast<std::size_t>(entity_of_slot[slot0 + s])];
                local[ncorners + nedges + nfaces] = body_of_cell[d.mFirstCell + c];

                const std::uint16_t mask = masks[d.mFirstCell + c];
                const detail::RefineMaskTemplate& tpl = detail::refine_mask_template(d.mType, mask);
                // The one choice made from GLOBAL node ids rather than the
                // template's local numbering, so that the cell across the
                // affected face resolves it the same way.
                const bool alt = tpl.HasVariant() && local[tpl.mTieB] < local[tpl.mTieA];
                const std::vector<std::vector<std::uint8_t>>& children =
                    alt ? tpl.mChildrenAlt : tpl.mChildren;
                const bool is_red = mask == d.mFullMask && mask != 0;
                for (std::size_t k = 0; k < children.size(); ++k) {
                    const std::vector<std::uint8_t>& child = children[k];
                    const std::size_t row_index = static_cast<std::size_t>(first_child[c]) + k;
                    std::int64_t* row = dst + row_index * out_npc;
                    for (std::size_t n = 0; n < out_npc; ++n)
                        row[n] = local[child[n]];
                    parent_of_child[b][row_index] = static_cast<std::int64_t>(c);
                    red_child[b][row_index] = is_red ? 1 : 0;
                }
            });
            out.AddCellBlock(cell_type_name(d.mType), std::move(out_conn));
            cell_maps.push_back(refine_int64_vector(first_child));
        }
    }

    // --- phase 7: point_data -------------------------------------------------
    for (const std::string& name : rMesh.PointDataNames()) {
        // The two reserved point arrays are this operation's own bookkeeping and
        // are rebuilt below from the output's topology. Letting them through here
        // would interpolate them: a hanging flag would be averaged to 0.5 and
        // truncated, and an entity key to the mean of two unrelated node ids.
        if (name == kRefineHangingName || name == kRefineEntityName)
            continue;
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
        parallel_for_bw(new_entities.size(), [&](std::size_t i) {
            const RefineNodeKey& key = entities[static_cast<std::size_t>(new_entities[i])];
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
                if (d.mType != CellType::Hexahedron || d.mNumCells == 0)
                    continue;
                const NDArray& conn = cb.Conn();
                const std::size_t npc = cb.NodesPerCell();
                const std::size_t ncorners = d.mNumCorners;
                const double inv = 1.0 / static_cast<double>(ncorners);
                parallel_for_bw(d.mNumCells, [&](std::size_t c) {
                    const std::int64_t body = body_of_cell[d.mFirstCell + c];
                    if (body < 0)
                        return;
                    for (std::size_t k = 0; k < ncomp; ++k) {
                        double sum = 0.0;
                        for (std::size_t n = 0; n < ncorners; ++n) {
                            const std::size_t p =
                                static_cast<std::size_t>(detail::read_int(conn, c * npc + n));
                            sum += detail::read_double(a, p * ncomp + k);
                        }
                        detail::write_double(b, static_cast<std::size_t>(body) * ncomp + k,
                                             sum * inv);
                    }
                });
            }
        }
        out.AddPointData(name, std::move(b));
    }

    // --- phase 8: cell_data gathered parent -> children ----------------------
    // base_levels was hoisted above: the Balanced closure compares levels.
    for (const std::string& name : rMesh.CellDataNames()) {
        if (name == kRefineLevelName || name == kRefineCellIdName || name == kRefineParentIdName)
            continue;  // this operation's own bookkeeping, rebuilt separately -- gathering
                       // them here would replicate a split cell's id/parent to every child
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
            const std::vector<std::int64_t>& parents = parent_of_child[b];
            std::vector<std::size_t> shape = a.Shape();
            shape[0] = parents.size();
            const std::size_t row_bytes = a.Nbytes() / in_rows;
            NDArray gathered = NDArray::Uninit(a.Dtype(), std::move(shape));
            const std::byte* src = a.Data();
            std::byte* dst = gathered.Data();
            parallel_for_bw(parents.size(), [&](std::size_t i) {
                std::memcpy(dst + i * row_bytes,
                            src + static_cast<std::size_t>(parents[i]) * row_bytes, row_bytes);
            });
            blocks.push_back(std::move(gathered));
        }
        out.AddCellData(name, std::move(blocks));
    }

    // refine:level. Red children increment; green and untouched cells inherit,
    // because a green split is a closure, not a refinement.
    if (RecordLevels || !base_levels.empty()) {
        std::vector<NDArray> blocks;
        blocks.reserve(nblocks);
        for (std::size_t b = 0; b < nblocks; ++b) {
            const std::vector<std::int64_t>& parents = parent_of_child[b];
            NDArray a = NDArray::Uninit(DType::Int64, {parents.size(), 1});
            std::int64_t* dst = a.As<std::int64_t>();
            for (std::size_t i = 0; i < parents.size(); ++i) {
                const std::size_t parent =
                    descs[b].mFirstCell + static_cast<std::size_t>(parents[i]);
                const std::int64_t base = base_levels.empty() ? 0 : base_levels[parent];
                dst[i] = base + (red_child[b][i] ? 1 : 0);
            }
            blocks.push_back(std::move(a));
        }
        out.AddCellData(kRefineLevelName, std::move(blocks));
    }

    // --- phase 8: the output's own entity keys, and the hanging set ----------
    // Every point's entity, in OUTPUT point ids -- which equal input ids for the
    // points that were already there, since refine never renumbers. Body centres
    // and original points keep the {-1,-1,-1,-1} sentinel: interior to one cell
    // or created by nothing, so no other cell can ever need to look one up.
    std::vector<RefineNodeKey> out_entity_keys(num_points_out, RefineNodeKey{-1, -1, -1, -1});
    for (std::size_t p = 0; p < base_entity_keys.size(); ++p)
        out_entity_keys[p] = base_entity_keys[p];
    for (std::size_t i = 0; i < new_entities.size(); ++i)
        out_entity_keys[num_points + i] = entities[static_cast<std::size_t>(new_entities[i])];

    // Hanging nodes. The rule is "a node sits on an entity of a cell that does
    // not reference it", and it must be evaluated over the EMITTED cells, not
    // over the input cells' entities. Those are not the same question once a
    // pass splits more than one level of the mesh: when the balance rule draws a
    // coarse cell in, its children's own sub-edges are entities that the input
    // cell never had, and a node the neighbouring refinement places on one of
    // them is constrained without any input entity ever naming it. Stating the
    // rule over the input is what left 42 of 84 constrained nodes unreported on
    // a two-level balanced refinement of a 4x4x4 hex block.
    //
    // Serial, and over a map keyed by the entity rather than by coordinates, so
    // there is no tolerance anywhere and the traversal order is a pure function
    // of (block, cell, slot) as everywhere else in this file.
    std::vector<char> hanging(num_points_out, 0);
    std::size_t num_hanging = 0;
    // If every cell was split fully and the input carried no nodes of its own,
    // every entity of every emitted cell is resolved by that cell's own template
    // and nothing can be constrained. That is the uniform path, so it skips the
    // scan below entirely rather than paying a second phase-1-sized pass for an
    // answer that is zero by construction.
    bool any_partial = false;
    for (std::size_t b = 0; b < nblocks && !any_partial; ++b)
        for (std::size_t c = 0; c < descs[b].mNumCells; ++c)
            if (masks[descs[b].mFirstCell + c] != descs[b].mFullMask) {
                any_partial = true;
                break;
            }
    if (any_partial || !base_entity_keys.empty()) {
        std::unordered_map<RefineNodeKey, std::int64_t, RefineNodeKeyHash> node_at_entity;
        node_at_entity.reserve(num_points_out * 2);
        for (std::size_t p = 0; p < num_points_out; ++p) {
            if (out_entity_keys[p][3] < 0)
                continue;
            node_at_entity.emplace(out_entity_keys[p], static_cast<std::int64_t>(p));
        }
        if (!node_at_entity.empty()) {
            std::size_t bi = 0;
            for (const auto cb : out.CellRange()) {
                const RefineBlockDesc& d = descs[bi++];
                const std::vector<detail::CellEdgePair>& edges = *d.mpEdges;
                const std::vector<detail::CellQuadFace>& faces = *d.mpFaces;
                const std::size_t ncells = cb.NumCells();
                if (ncells == 0)
                    continue;
                const NDArray& conn = cb.Conn();
                const std::size_t npc = cb.NodesPerCell();
                for (std::size_t c = 0; c < ncells; ++c) {
                    const std::size_t row = c * npc;
                    const auto mark = [&](const RefineNodeKey& rKey) {
                        auto it = node_at_entity.find(rKey);
                        if (it == node_at_entity.end())
                            return;
                        const std::size_t node = static_cast<std::size_t>(it->second);
                        for (std::size_t n = 0; n < npc; ++n)
                            if (detail::read_int(conn, row + n) == it->second)
                                return;  // this cell references it: not hanging here
                        if (!hanging[node]) {
                            hanging[node] = 1;
                            ++num_hanging;
                        }
                    };
                    for (const detail::CellEdgePair& e : edges)
                        mark(refine_edge_key(detail::read_int(conn, row + e[0]),
                                             detail::read_int(conn, row + e[1])));
                    for (const detail::CellQuadFace& f : faces)
                        mark(refine_face_key(detail::read_int(conn, row + f[0]),
                                             detail::read_int(conn, row + f[1]),
                                             detail::read_int(conn, row + f[2]),
                                             detail::read_int(conn, row + f[3])));
                }
            }
        }
    }

    // refine:hanging -- the constrained nodes a Balanced pass leaves behind. Only
    // attached when there are any, so the conforming closures are unaffected.
    if (num_hanging > 0) {
        NDArray flags = NDArray::Uninit(DType::Int64, {num_points_out, 1});
        std::int64_t* dst = flags.As<std::int64_t>();
        for (std::size_t i = 0; i < num_points_out; ++i)
            dst[i] = hanging[i] ? 1 : 0;
        out.AddPointData(kRefineHangingName, std::move(flags));
        log::info("refine: the balanced closure left {} hanging node(s); see '{}'.", num_hanging,
                  kRefineHangingName);
    }

    // refine:entity -- so a later pass can reuse this pass's nodes instead of
    // allocating coincident duplicates. Attached whenever this pass leaves
    // hanging nodes, and *maintained* thereafter once present, exactly as
    // refine:level is. A pass that leaves none and inherits none writes nothing,
    // which is what keeps every conforming closure's output byte-identical to
    // what it was before the array existed -- UNLESS `ForceEntity`, which
    // attaches it regardless: `mRecordHierarchy` wants it as the multigrid
    // prolongation stencil, and `out_entity_keys` is already fully built above
    // whether or not any node in it is hanging.
    if (num_hanging > 0 || !base_entity_keys.empty() || ForceEntity) {
        NDArray keys_out = NDArray::Uninit(DType::Int64, {num_points_out, std::size_t{4}});
        std::int64_t* dst = keys_out.As<std::int64_t>();
        for (std::size_t p = 0; p < num_points_out; ++p)
            for (std::size_t k = 0; k < 4; ++k)
                dst[p * 4 + k] = out_entity_keys[p][k];
        out.AddPointData(kRefineEntityName, std::move(keys_out));
    }

    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    if (pOutRedChildren != nullptr) {
        pOutRedChildren->clear();
        for (std::size_t b = 0; b < nblocks; ++b)
            pOutRedChildren->insert(pOutRedChildren->end(), red_child[b].begin(),
                                    red_child[b].end());
    }

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
    rResult.mMesh.AddCellData(kRefineParentCellName, std::move(blocks));
}

// Attach refine:cell_id/refine:parent_id from the FINAL composed cell maps --
// one serial pass over the whole call's result, exactly parity with
// refine_attach_parent_ids above and for the same reason: `rResult.mCellMaps`
// already composes through every internal level (refine_compose_maps), so
// there is no need to touch refine_once or thread anything through the levels
// loop. A composed run of length one is a cell no level ever split -- it keeps
// its own id and is its own parent; every other cell in a longer run gets one
// fresh id from a single counter walking the mesh's own block/cell order, and
// carries the run's original (in THIS call's input) ancestor as its parent.
void refine_attach_hierarchy(const Mesh& rIn, RefineResult& rResult, bool RecordHierarchy) {
    const std::vector<std::int64_t> bases = detail::block_bases(rIn);
    std::vector<std::int64_t> input_ids;
    std::int64_t id_base = 0;
    const RefineHierarchyState state = refine_read_hierarchy(rIn, bases, input_ids, id_base);
    const bool maintained = state == RefineHierarchyState::Valid;
    if (!maintained && !RecordHierarchy)
        return;  // nothing valid to carry, and nothing asked for
    if (!maintained) {
        // Absent, or dropped as malformed/non-unique (warned inside
        // refine_read_hierarchy already): start a fresh id space from the
        // mesh's implicit global block-major ids.
        const std::int64_t total = detail::total_cells(bases);
        input_ids.assign(static_cast<std::size_t>(total), 0);
        for (std::int64_t i = 0; i < total; ++i)
            input_ids[static_cast<std::size_t>(i)] = i;
        id_base = total;
    }

    const Mesh& mesh = rResult.mMesh;
    const std::size_t nblocks = mesh.NumCellBlocks();
    if (nblocks != rResult.mCellMaps.size() || nblocks == 0)
        return;
    std::vector<NDArray> id_blocks;
    std::vector<NDArray> parent_blocks;
    id_blocks.reserve(nblocks);
    parent_blocks.reserve(nblocks);
    std::int64_t next_id = id_base;
    std::size_t bi = 0;
    for (const auto cb : mesh.CellRange()) {
        const NDArray& map = rResult.mCellMaps[bi];
        const std::int64_t* first = map.As<std::int64_t>();
        const std::size_t nparents = map.Size();
        const std::size_t ncells = cb.NumCells();
        const std::size_t base = static_cast<std::size_t>(bases[bi]);
        NDArray ids = NDArray::Uninit(DType::Int64, {ncells, 1});
        NDArray parents = NDArray::Uninit(DType::Int64, {ncells, 1});
        std::int64_t* id_dst = ids.As<std::int64_t>();
        std::int64_t* parent_dst = parents.As<std::int64_t>();
        for (std::size_t p = 0; p < nparents; ++p) {
            if (first[p] < 0)
                continue;
            const std::size_t lo = static_cast<std::size_t>(first[p]);
            const std::size_t hi = p + 1 < nparents && first[p + 1] >= 0
                                       ? static_cast<std::size_t>(first[p + 1])
                                       : ncells;
            const std::int64_t original_id = input_ids[base + p];
            if (hi - lo == 1 && lo < ncells) {
                // A run of length one is a cell no level of this call ever
                // split (every admissible non-empty mask yields more than one
                // child -- see detail/refine_templates.hpp), so it keeps its
                // own identity and is its own parent.
                id_dst[lo] = original_id;
                parent_dst[lo] = original_id;
            } else {
                for (std::size_t c = lo; c < hi && c < ncells; ++c) {
                    id_dst[c] = next_id++;
                    parent_dst[c] = original_id;
                }
            }
        }
        id_blocks.push_back(std::move(ids));
        parent_blocks.push_back(std::move(parents));
        ++bi;
    }
    rResult.mMesh.AddCellData(kRefineCellIdName, std::move(id_blocks));
    rResult.mMesh.AddCellData(kRefineParentIdName, std::move(parent_blocks));
}

// Cells the next level would produce, for the pre-flight size warning.
std::size_t refine_projected_cells(const Mesh& rMesh) {
    std::size_t total = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsRagged())
            continue;
        const CellType type = cell_type_from_name(std::string(cb.Type()));
        const std::size_t children =
            detail::refine_mask_template(type, detail::refine_full_mask(type)).NumChildren();
        total += cb.NumCells() * (children == 0 ? 1 : children);
    }
    return total;
}

// --- selection ---------------------------------------------------------------

bool refine_has_selector(const RefineOptions& rOptions) {
    return !rOptions.mCells.empty() || !rOptions.mRegion.empty() ||
           !rOptions.mPredicateArray.empty();
}

std::string refine_region_names(const Mesh& rMesh) {
    std::string names;
    for (const std::string& n : rMesh.RegionNames()) {
        if (!names.empty())
            names += ", ";
        names += "'" + n + "'";
    }
    return names.empty() ? "none" : names;
}

// Resolve the selector into a per-global-cell red flag. Only called when
// `refine_has_selector`, so "no selector" and "a selector that matched nothing"
// stay distinguishable -- the latter is a legitimate request for an unchanged
// mesh, not a silent fall-through to refining everything.
std::vector<char> refine_resolve_selection(const Mesh& rMesh, const RefineOptions& rOptions) {
    const int num_set = (rOptions.mCells.empty() ? 0 : 1) + (rOptions.mRegion.empty() ? 0 : 1) +
                        (rOptions.mPredicateArray.empty() ? 0 : 1);
    if (num_set > 1)
        throw std::invalid_argument(
            "refine: set at most one cell selector -- cells, region or predicate array -- but " +
            std::to_string(num_set) + " were given");

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::int64_t ncells = detail::total_cells(bases);
    std::vector<char> red(static_cast<std::size_t>(ncells), 0);

    if (!rOptions.mCells.empty()) {
        for (std::int64_t g : rOptions.mCells) {
            if (g < 0 || g >= ncells)
                throw std::invalid_argument("refine: cell index " + std::to_string(g) +
                                            " is out of range; the mesh has " +
                                            std::to_string(ncells) + " cells");
            red[static_cast<std::size_t>(g)] = 1;  // duplicates collapse here
        }
        return red;
    }

    if (!rOptions.mRegion.empty()) {
        std::size_t idx = rMesh.FindRegion(rOptions.mRegion, RegionKind::Cell);
        if (idx != Mesh::npos) {
            const Region& region = rMesh.Region(idx);
            const std::int64_t* entries = region.Entries();
            for (std::size_t i = 0; i < region.NumEntries(); ++i) {
                const std::int64_t g = entries[i];
                if (g >= 0 && g < ncells)
                    red[static_cast<std::size_t>(g)] = 1;
            }
            return red;
        }
        idx = rMesh.FindRegion(rOptions.mRegion, RegionKind::Point);
        if (idx != Mesh::npos) {
            const Region& region = rMesh.Region(idx);
            std::vector<char> in_region(rMesh.NumPoints(), 0);
            const std::int64_t* entries = region.Entries();
            for (std::size_t i = 0; i < region.NumEntries(); ++i) {
                const std::int64_t p = entries[i];
                if (p >= 0 && p < static_cast<std::int64_t>(rMesh.NumPoints()))
                    in_region[static_cast<std::size_t>(p)] = 1;
            }
            // A cell is selected when ANY of its nodes is in the region. "All"
            // would select nothing at all for a point set describing a surface.
            std::size_t b = 0;
            for (const auto cb : rMesh.CellRange()) {
                const std::size_t base = static_cast<std::size_t>(bases[b++]);
                if (cb.IsRagged())
                    continue;  // rejected by refine_check_block later, by name
                const NDArray& conn = cb.Conn();
                const std::size_t npc = cb.NodesPerCell();
                for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                    for (std::size_t n = 0; n < npc; ++n) {
                        const std::int64_t p = detail::read_int(conn, c * npc + n);
                        if (p >= 0 && p < static_cast<std::int64_t>(in_region.size()) &&
                            in_region[static_cast<std::size_t>(p)]) {
                            red[base + c] = 1;
                            break;
                        }
                    }
                }
            }
            return red;
        }
        if (rMesh.FindRegion(rOptions.mRegion, RegionKind::Side) != Mesh::npos)
            throw std::invalid_argument(
                "refine: region '" + rOptions.mRegion +
                "' is a side region; a facet is not a cell, so it cannot select what to refine "
                "(use a cell or point region)");
        throw std::invalid_argument("refine: no region named '" + rOptions.mRegion +
                                    "' (available: " + refine_region_names(rMesh) + ")");
    }

    const std::string& name = rOptions.mPredicateArray;
    if (!rMesh.HasCellData(name))
        throw std::invalid_argument("refine: " +
                                    data_unknown_key_message(rMesh, DataLocation::Cell, name));
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (rMesh.CellDataNumBlocks(name) != nblocks)
        throw std::invalid_argument("refine: cell_data '" + name +
                                    "' does not cover every cell block");
    std::size_t b = 0;
    for (const auto cb : rMesh.CellRange()) {
        const NDArray& a = rMesh.CellData(name, b);
        const std::size_t base = static_cast<std::size_t>(bases[b]);
        ++b;
        if (detail::rows(a) != cb.NumCells())
            throw std::invalid_argument("refine: cell_data '" + name + "' has " +
                                        std::to_string(detail::rows(a)) + " rows on a block of " +
                                        std::to_string(cb.NumCells()) + " cells");
        if (cb.NumCells() != 0 && data_num_components(a) != 1)
            throw std::invalid_argument("refine: cell_data '" + name +
                                        "' must be scalar (one value per cell) to be used as a "
                                        "refinement predicate");
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            if (refine_compare_value(detail::read_double(a, c), rOptions.mPredicateOp,
                                     rOptions.mPredicateValue))
                red[base + c] = 1;
        }
    }
    return red;
}

}  // namespace

// --- public API --------------------------------------------------------------

RefineClosure refine_closure_from_name(const std::string& rName) {
    if (rName.empty() || rName == "redgreen" || rName == "red-green" || rName == "green")
        return RefineClosure::RedGreen;
    if (rName == "propagate" || rName == "red")
        return RefineClosure::Propagate;
    if (rName == "balanced" || rName == "2:1")
        return RefineClosure::Balanced;
    throw std::invalid_argument("refine: unknown closure '" + rName +
                                "' (expected 'redgreen'/'green', 'propagate' or 'balanced')");
}

bool refine_compare_value(double Value, RefineCompare Op, double Rhs) {
    // A non-finite value never matches. compute_quality deliberately reports NaN
    // where a metric does not apply, so a predicate over `quality:*` on a mixed
    // mesh is the headline use case -- rejecting the array would break it.
    //
    // Public (declared in refine.hpp) rather than file-private because
    // `crop_predicate` selects cells by this identical rule; see the header.
    if (!std::isfinite(Value))
        return false;
    switch (Op) {
        case RefineCompare::Less:
            return Value < Rhs;
        case RefineCompare::LessEqual:
            return Value <= Rhs;
        case RefineCompare::Greater:
            return Value > Rhs;
        case RefineCompare::GreaterEqual:
            return Value >= Rhs;
        case RefineCompare::Equal:
            return Value == Rhs;
        case RefineCompare::NotEqual:
            return Value != Rhs;
    }
    return false;
}

RefineCompare refine_compare_from_name(const std::string& rName) {
    if (rName == "<" || rName == "lt")
        return RefineCompare::Less;
    if (rName == "<=" || rName == "le")
        return RefineCompare::LessEqual;
    if (rName == ">" || rName == "gt")
        return RefineCompare::Greater;
    if (rName == ">=" || rName == "ge")
        return RefineCompare::GreaterEqual;
    if (rName == "==" || rName == "=" || rName == "eq")
        return RefineCompare::Equal;
    if (rName == "!=" || rName == "ne")
        return RefineCompare::NotEqual;
    throw std::invalid_argument("refine: unknown comparison '" + rName +
                                "' (expected '<', '<=', '>', '>=', '==' or '!=')");
}

namespace {

// Carry the input's named regions onto the output. The cell map is FirstChild:
// a parent's children occupy a contiguous run -- of length 1 for a cell the
// selection left untouched, which is still a run and still non-negative. Side
// regions are dropped: a child cell is a new cell of a subdivided topology, so
// its facets have no correspondence with the parent's. See
// detail/region_remap.hpp.
void refine_carry_regions(const Mesh& rIn, RefineResult& rRes) {
    detail::RegionRemap rmap;
    rmap.pPointMap = &rRes.mPointMap;
    rmap.mCellMapKind = detail::CellMapKind::FirstChild;
    rmap.pCellMaps = &rRes.mCellMaps;
    rmap.mOpName = "refine";
    detail::remap_regions(rIn, rRes.mMesh, rmap);
}

}  // namespace

RefineResult refine(const Mesh& rMesh, const RefineOptions& rOptions) {
    // Resolve (and validate) the selection against the INPUT mesh, before the
    // level loop: a selector names input cells, and level k > 1 refines the
    // children of level k - 1's red cells rather than re-resolving it.
    const bool selective = refine_has_selector(rOptions);
    std::vector<char> seed;
    if (selective)
        seed = refine_resolve_selection(rMesh, rOptions);

    if (rOptions.mLevels <= 0) {
        RefineResult res;
        res.mMesh = detail::clone_mesh(rMesh);
        res.mPointMap = refine_identity_map(rMesh.NumPoints());
        res.mCellMaps.reserve(rMesh.NumCellBlocks());
        for (const auto cb : rMesh.CellRange())
            res.mCellMaps.push_back(refine_identity_map(cb.NumCells()));
        if (rOptions.mRecordParentIds)
            refine_attach_parent_ids(res);
        refine_attach_hierarchy(rMesh, res, rOptions.mRecordHierarchy);
        refine_carry_regions(rMesh, res);
        return res;
    }

    // Refinement is exponential in `levels`; the failure mode at depth is a
    // bad_alloc rather than a wrong answer, so say so before it happens. Only
    // meaningful for the uniform path -- a selection's growth is bounded by the
    // closure, which is exactly what cannot be predicted from the cell count.
    if (!selective) {
        static constexpr std::size_t refine_warn_cells = 20'000'000;
        const std::size_t projected = refine_projected_cells(rMesh);
        if (projected > refine_warn_cells)
            log::warn("refine: one level of this mesh yields ~{} cells; {} level(s) requested.",
                      projected, rOptions.mLevels);
    }

    std::vector<char> next_seed;
    RefineResult acc =
        refine_once(rMesh, selective ? &seed : nullptr, rOptions.mClosure, rOptions.mRecordLevels,
                    rOptions.mRecordHierarchy, selective ? &next_seed : nullptr);
    for (int level = 1; level < rOptions.mLevels; ++level) {
        const std::vector<char> current = std::move(next_seed);
        RefineResult next = refine_once(
            acc.mMesh, selective ? &current : nullptr, rOptions.mClosure, rOptions.mRecordLevels,
            rOptions.mRecordHierarchy, selective ? &next_seed : nullptr);
        refine_compose_maps(acc, next);     // reads next's maps first
        acc.mMesh = std::move(next.mMesh);  // sequenced after
    }
    if (rOptions.mRecordParentIds)
        refine_attach_parent_ids(acc);
    // Reads the ORIGINAL rMesh, never an intermediate: refine_attach_hierarchy
    // resolves entirely through acc.mCellMaps, which already composes through
    // every internal level, so "the mesh passed to this call" is unambiguous
    // regardless of mLevels.
    refine_attach_hierarchy(rMesh, acc, rOptions.mRecordHierarchy);
    refine_carry_regions(rMesh, acc);
    return acc;
}

}  // namespace meshioplusplus
