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
// Mesh partitioning into N balanced parts: Hilbert space-filling-curve cut
// (always available) or KaHIP kaffpa() on the shared-face dual graph (optional,
// -DMESHIOPLUSPLUS_WITH_KAHIP=ON; the KaHIP wrapper is ported from the Kratos
// KaHIPApplication). Built entirely through the uniform mesh API so it compiles
// under every mesh backend. See operations/partition.hpp for the contract.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/detail/cell_adjacency.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/space_filling.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

#ifdef MESHIOPLUSPLUS_HAS_KAHIP
// System includes (KaHIP dual-graph path only)
#include <array>
#include <functional>
#include <unordered_map>
// External includes
#include <kaHIP_interface.h>  // IWYU pragma: keep
// Project includes (KaHIP dual-graph path only)
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#endif

namespace meshioplusplus {

namespace {

// --- shared helpers ----------------------------------------------------------

// Total cell count across all blocks (the global cell numbering size).
std::size_t partition_total_cells(const Mesh& rMesh) {
    std::size_t total = 0;
    for (const auto cb : rMesh.CellRange())
        total += cb.NumCells();
    return total;
}

// Per-cell centroids in global cell order, always 3 doubles per cell (unused
// axes zero). The mean of the cell's nodes; polyhedra use the distinct nodes
// across their faces (a node shared by several faces counts once).
std::vector<double> partition_centroids(const Mesh& rMesh, std::size_t total) {
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    const std::size_t dim = std::min<std::size_t>(pdim, 3);

    std::vector<double> cent(3 * total, 0.0);
    std::size_t base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        if (cb.IsPolyhedron()) {
            parallel_for(ncells, [&, base](std::size_t c) {
                std::vector<std::int64_t> nodes;
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    const auto face = cb.Face(c, f);
                    nodes.insert(nodes.end(), face.first, face.first + face.second);
                }
                std::sort(nodes.begin(), nodes.end());
                nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
                double* out = cent.data() + 3 * (base + c);
                for (std::int64_t node : nodes)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] +=
                            detail::read_double(points, static_cast<std::size_t>(node) * pdim + d);
                if (!nodes.empty())
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] /= static_cast<double>(nodes.size());
            });
        } else if (cb.IsRagged()) {
            parallel_for(ncells, [&, base](std::size_t c) {
                const std::int64_t* row = cb.Row(c);
                const std::size_t nn = cb.RowSize(c);
                double* out = cent.data() + 3 * (base + c);
                for (std::size_t k = 0; k < nn; ++k)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] += detail::read_double(points,
                                                      static_cast<std::size_t>(row[k]) * pdim + d);
                if (nn > 0)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] /= static_cast<double>(nn);
            });
        } else {
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            parallel_for(ncells, [&, base](std::size_t c) {
                double* out = cent.data() + 3 * (base + c);
                for (std::size_t k = 0; k < npc; ++k) {
                    const std::size_t node =
                        static_cast<std::size_t>(detail::read_int(conn, c * npc + k));
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] += detail::read_double(points, node * pdim + d);
                }
                if (npc > 0)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] /= static_cast<double>(npc);
            });
        }
        base += ncells;
    }
    return cent;
}

// The validated per-cell weights (global cell order) from a scalar numeric
// cell_data array.
std::vector<double> partition_weights(const Mesh& rMesh, const std::string& rKey,
                                      std::size_t total) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (!rMesh.HasCellData(rKey))
        throw std::invalid_argument("partition: " +
                                    data_unknown_key_message(rMesh, DataLocation::Cell, rKey));
    if (rMesh.CellDataNumBlocks(rKey) != nblocks)
        throw std::invalid_argument("partition: weights cell_data '" + rKey +
                                    "' does not cover every cell block");

    std::vector<double> weights;
    weights.reserve(total);
    double sum = 0.0;
    std::size_t b = 0;
    for (const auto cb : rMesh.CellRange()) {
        const NDArray& a = rMesh.CellData(rKey, b);
        const std::size_t rows = a.Shape().empty() ? 0 : a.Shape()[0];
        if (rows != cb.NumCells())
            throw std::invalid_argument("partition: weights cell_data '" + rKey + "' has " +
                                        std::to_string(rows) + " rows on a block of " +
                                        std::to_string(cb.NumCells()) + " cells");
        if (data_num_components(a) != 1)
            throw std::invalid_argument("partition: weights cell_data '" + rKey +
                                        "' must be scalar (one component per cell)");
        for (std::size_t c = 0; c < rows; ++c) {
            const double w = detail::read_double(a, c);
            if (!std::isfinite(w))
                throw std::invalid_argument("partition: weights cell_data '" + rKey +
                                            "' contains a non-finite value");
            if (w < 0.0)
                throw std::invalid_argument("partition: weights cell_data '" + rKey +
                                            "' contains a negative value");
            weights.push_back(w);
            sum += w;
        }
        ++b;
    }
    if (!(sum > 0.0))
        throw std::invalid_argument("partition: weights cell_data '" + rKey +
                                    "' has zero total weight");
    return weights;
}

// --- SFC method --------------------------------------------------------------

// Hilbert-curve cut of the cell centroids. Deterministic by construction: the
// keys are filled into disjoint slots in parallel, the argsort is a serial
// stable sort (ties broken by cell index) and the cut is a serial integer /
// prefix-sum rule, so the assignment is byte-identical across mesh backends
// and thread counts.
std::vector<int> partition_sfc_parts(const Mesh& rMesh, const PartitionOptions& rOptions,
                                     const std::vector<double>& rWeights, std::size_t total) {
    const int nparts = rOptions.mNParts;
    std::vector<int> parts(total, 0);
    if (total == 0 || nparts == 1)
        return parts;

    const std::vector<double> cent = partition_centroids(rMesh, total);
    const std::size_t dim = std::min<std::size_t>(rMesh.PointDim(), 3);
    const int bits = detail::sfc_bits;

    // Bounding box of the centroids (serial, deterministic).
    double lo[3] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[3] = {-std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity()};
    for (std::size_t i = 0; i < total; ++i) {
        for (std::size_t d = 0; d < dim; ++d) {
            const double c = cent[3 * i + d];
            lo[d] = std::min(lo[d], c);
            hi[d] = std::max(hi[d], c);
        }
    }
    double scale[3] = {0.0, 0.0, 0.0};
    const double qmax = static_cast<double>((std::int64_t(1) << bits) - 1);
    for (std::size_t d = 0; d < dim; ++d) {
        const double span = hi[d] - lo[d];
        scale[d] = span > 0.0 ? qmax / span : 0.0;
    }

    // Quantize + Hilbert key per centroid (disjoint slots, parallel-safe).
    std::vector<std::uint64_t> keys(total);
    parallel_for(total, [&](std::size_t i) {
        std::uint32_t q[3] = {0, 0, 0};
        for (std::size_t d = 0; d < dim; ++d) {
            const double t = (cent[3 * i + d] - lo[d]) * scale[d];
            std::int64_t qi = static_cast<std::int64_t>(t + 0.5);
            if (qi < 0)
                qi = 0;
            const std::int64_t maxq = (std::int64_t(1) << bits) - 1;
            if (qi > maxq)
                qi = maxq;
            q[d] = static_cast<std::uint32_t>(qi);
        }
        keys[i] = detail::sfc_hilbert_key(q, bits);
    });

    // Serial stable argsort along the curve (ties broken by cell index).
    std::vector<std::int64_t> order(total);
    std::iota(order.begin(), order.end(), std::int64_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::int64_t a, std::int64_t b) {
        return keys[static_cast<std::size_t>(a)] < keys[static_cast<std::size_t>(b)];
    });

    // Serial cut into nparts contiguous ranges.
    if (rWeights.empty()) {
        // Pure integer rule: part sizes are floor(n/k) or ceil(n/k) (<= 1 apart).
        for (std::size_t r = 0; r < total; ++r)
            parts[static_cast<std::size_t>(order[r])] = static_cast<int>(
                (static_cast<std::int64_t>(r) * nparts) / static_cast<std::int64_t>(total));
    } else {
        // Each cell goes to the ideal interval containing its weight midpoint;
        // midpoints are nondecreasing along the curve, so parts stay contiguous.
        double sum = 0.0;
        for (double w : rWeights)
            sum += w;
        double prefix = 0.0;
        for (std::size_t r = 0; r < total; ++r) {
            const double w = rWeights[static_cast<std::size_t>(order[r])];
            int p = static_cast<int>(std::floor((prefix + 0.5 * w) * nparts / sum));
            p = std::min(nparts - 1, std::max(0, p));
            parts[static_cast<std::size_t>(order[r])] = p;
            prefix += w;
        }
    }
    return parts;
}

// --- KaHIP method ------------------------------------------------------------

#ifdef MESHIOPLUSPLUS_HAS_KAHIP

// Sorted corner ids of one facet, padded with -1 up to 4 entries (surface.cpp's
// facet-key idiom; the padding lets triangular and quadrilateral facets share
// one map with no discriminator).
using PartitionFacetKey = std::array<std::int64_t, 4>;

struct PartitionFacetKeyHash {
    std::size_t operator()(const PartitionFacetKey& rKey) const {
        std::size_t h = 0;
        for (std::int64_t v : rKey)
            h ^= std::hash<std::int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

// The corner-node rows of a cell type's facets in the dual-graph sense: faces
// for 3D cells, edges for 2D cells, empty (isolated vertex) otherwise.
struct PartitionFacetDef {
    std::uint8_t mNumCorners;
    std::array<std::uint8_t, 4> mNodes;
};

std::vector<PartitionFacetDef> partition_facets_for(CellType type, int dualDim) {
    std::vector<PartitionFacetDef> defs;
    if (dualDim == 3) {
        for (const detail::CellFaceDef& f : detail::cell_faces(type)) {
            PartitionFacetDef d;
            d.mNumCorners = f.mNumCorners;
            for (std::uint8_t k = 0; k < f.mNumCorners; ++k)
                d.mNodes[k] = f.mNodes[k];
            defs.push_back(d);
        }
    } else if (dualDim == 2) {
        for (const detail::CellEdgeDef& e : detail::cell_edges(type)) {
            PartitionFacetDef d;
            d.mNumCorners = e.mNumCorners;
            for (std::uint8_t k = 0; k < e.mNumCorners; ++k)
                d.mNodes[k] = e.mNodes[k];
            defs.push_back(d);
        }
    }
    return defs;
}

// One facet occurrence recorded in phase 1 (keys built in parallel).
struct PartitionFacetRecord {
    PartitionFacetKey mKey;
    std::int64_t mParent;  // global cell index owning the facet
};

// The dual graph in CSR form: cells are vertices, an edge where two cells
// share a facet (each undirected edge stored twice, 0-indexed).
struct PartitionCsr {
    std::vector<std::int64_t> mXadj;    // size total + 1
    std::vector<std::int64_t> mAdjncy;  // size mXadj.back()
};

PartitionCsr partition_dual_graph(const Mesh& rMesh, std::size_t total) {
    // The dual dimension: faces when any 3D cell exists, else edges. Blocks of
    // lower dimension (or without a facet table) become isolated vertices.
    int dual_dim = 0;
    for (const auto cb : rMesh.CellRange()) {
        const int d = cb.IsPolyhedron() ? 3 : cell_type_dimension(cell_type_from_name(cb.Type()));
        dual_dim = std::max(dual_dim, d);
    }

    // Per-block facet tables + the flat record-buffer layout.
    struct BlockDesc {
        const NDArray* mpConn;
        std::size_t mNpc;
        std::size_t mNumCells;
        std::vector<PartitionFacetDef> mFacets;
        std::size_t mFirstFacet;
        std::int64_t mGlobalCellBase;
    };
    std::vector<BlockDesc> descs;
    std::size_t total_facets = 0;
    std::int64_t global_cell_base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t base = global_cell_base;
        global_cell_base += static_cast<std::int64_t>(cb.NumCells());
        if (cb.IsRagged() || cb.IsPolyhedron())
            continue;  // isolated vertices (no facet table)
        const CellType ct = cell_type_from_name(cb.Type());
        if (cell_type_dimension(ct) != dual_dim)
            continue;
        std::vector<PartitionFacetDef> facets = partition_facets_for(ct, dual_dim);
        if (facets.empty())
            continue;
        BlockDesc d;
        d.mpConn = &cb.Conn();
        d.mNpc = cb.NodesPerCell();
        d.mNumCells = cb.NumCells();
        d.mFacets = std::move(facets);
        d.mFirstFacet = total_facets;
        d.mGlobalCellBase = base;
        total_facets += d.mNumCells * d.mFacets.size();
        descs.push_back(std::move(d));
    }

    // Phase 1: build facet keys into disjoint slots (parallel-safe).
    std::vector<PartitionFacetRecord> recs(total_facets);
    for (const BlockDesc& d : descs) {
        const NDArray& conn = *d.mpConn;
        const std::size_t npc = d.mNpc;
        const std::size_t fpc = d.mFacets.size();
        parallel_for(d.mNumCells * fpc, [&](std::size_t j) {
            const std::size_t cell = j / fpc;
            const PartitionFacetDef& facet = d.mFacets[j % fpc];
            PartitionFacetRecord& r = recs[d.mFirstFacet + j];
            r.mKey = {-1, -1, -1, -1};
            for (std::uint8_t k = 0; k < facet.mNumCorners; ++k)
                r.mKey[k] = detail::read_int(conn, cell * npc + facet.mNodes[k]);
            std::sort(r.mKey.begin(), r.mKey.end());
            r.mParent = d.mGlobalCellBase + static_cast<std::int64_t>(cell);
        });
    }

    // Phase 2: serial first-parent pairing (deterministic). Occurrences beyond
    // the second (a non-manifold facet) all connect to the first owner.
    std::unordered_map<PartitionFacetKey, std::int64_t, PartitionFacetKeyHash> first_owner;
    first_owner.reserve(total_facets * 2);
    std::vector<std::vector<std::int64_t>> adj(total);
    for (const PartitionFacetRecord& r : recs) {
        auto it = first_owner.find(r.mKey);
        if (it == first_owner.end()) {
            first_owner.emplace(r.mKey, r.mParent);
        } else if (it->second != r.mParent) {
            adj[static_cast<std::size_t>(it->second)].push_back(r.mParent);
            adj[static_cast<std::size_t>(r.mParent)].push_back(it->second);
        }
    }

    // Dedupe neighbour lists (cells can share more than one facet) and pack CSR.
    parallel_for(total, [&](std::size_t i) {
        std::vector<std::int64_t>& a = adj[i];
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
    });
    PartitionCsr csr;
    csr.mXadj.resize(total + 1);
    csr.mXadj[0] = 0;
    for (std::size_t i = 0; i < total; ++i)
        csr.mXadj[i + 1] = csr.mXadj[i] + static_cast<std::int64_t>(adj[i].size());
    csr.mAdjncy.reserve(static_cast<std::size_t>(csr.mXadj[total]));
    for (std::size_t i = 0; i < total; ++i)
        csr.mAdjncy.insert(csr.mAdjncy.end(), adj[i].begin(), adj[i].end());
    return csr;
}

// Per-cell part assignment via KaHIP's serial kaffpa() (ported from the Kratos
// KaHIPApplication's KaHIPPartitioner). The CSR buffers are built in the index
// width the *installed* KaHIP library uses -- kahip_sizeof_idx() at runtime,
// never the compile-time kahip_idx typedef, whose default (int32) need not
// match a -D64BITMODE=On build of the library.
std::vector<int> partition_kahip_parts(const Mesh& rMesh, const PartitionOptions& rOptions,
                                       const std::vector<double>& rWeights, std::size_t total) {
    const int nparts = rOptions.mNParts;
    if (total == 0)
        return {};
    if (nparts == 1)
        return std::vector<int>(total, 0);
    if (total > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument(
            "partition: mesh has " + std::to_string(total) +
            " cells; KaHIP's kaffpa interface indexes cells with a 32-bit int");

    PartitionCsr graph = partition_dual_graph(rMesh, total);

    // Vertex weights: kaffpa takes int; scale so the total maps to ~2^30
    // (relative weights preserved, sums safe in a 32-bit KaHIP build), with a
    // floor of 1 so no cell becomes weightless.
    std::vector<int> vwgt;
    if (!rWeights.empty()) {
        double sum = 0.0;
        for (double w : rWeights)
            sum += w;
        const double scale = static_cast<double>(std::int64_t(1) << 30) / sum;
        vwgt.resize(total);
        for (std::size_t i = 0; i < total; ++i)
            vwgt[i] = static_cast<int>(std::max<long long>(1, std::llround(rWeights[i] * scale)));
    }

    int n = static_cast<int>(total);
    int k = nparts;
    double imbalance = rOptions.mImbalance;  // kaffpa takes a non-const pointer
    int mode = ECO;
    switch (rOptions.mMode) {
        case PartitionMode::Fast:
            mode = FAST;
            break;
        case PartitionMode::Eco:
            mode = ECO;
            break;
        case PartitionMode::Strong:
            mode = STRONG;
            break;
    }
    int* p_vwgt = vwgt.empty() ? nullptr : vwgt.data();
    std::vector<int> part(total, 0);

    const int idx_width = kahip_sizeof_idx();
    if (idx_width == 8) {
        // 64-bit KaHIP: our CSR is already int64 -- pass it straight through.
        std::int64_t edgecut = 0;
        kaffpa(&n, p_vwgt, reinterpret_cast<kahip_idx*>(graph.mXadj.data()), nullptr,
               reinterpret_cast<kahip_idx*>(graph.mAdjncy.data()), &k, &imbalance,
               /*suppress_output=*/true, rOptions.mSeed, mode,
               reinterpret_cast<kahip_idx*>(&edgecut), part.data());
    } else {
        // 32-bit KaHIP: narrow the CSR, guarding against overflow.
        if (graph.mXadj[total] > std::numeric_limits<std::int32_t>::max())
            throw std::invalid_argument(
                "partition: the dual graph has " + std::to_string(graph.mXadj[total]) +
                " directed edges, which overflows this KaHIP build's 32-bit indices; "
                "rebuild KaHIP with -D64BITMODE=On");
        std::vector<std::int32_t> xadj(graph.mXadj.begin(), graph.mXadj.end());
        std::vector<std::int32_t> adjncy(graph.mAdjncy.begin(), graph.mAdjncy.end());
        std::int32_t edgecut = 0;
        kaffpa(&n, p_vwgt, reinterpret_cast<kahip_idx*>(xadj.data()), nullptr,
               reinterpret_cast<kahip_idx*>(adjncy.data()), &k, &imbalance,
               /*suppress_output=*/true, rOptions.mSeed, mode,
               reinterpret_cast<kahip_idx*>(&edgecut), part.data());
    }
    return part;
}

#else  // !MESHIOPLUSPLUS_HAS_KAHIP

// Always defined so `method == KaHIP` fails with an error naming the CMake
// option (registry_compiled_out() / vtk_codec_* spirit) -- never a link error,
// never a silent downgrade to SFC.
std::vector<int> partition_kahip_parts(const Mesh& rMesh, const PartitionOptions& rOptions,
                                       const std::vector<double>& rWeights, std::size_t total) {
    (void)rMesh;
    (void)rOptions;
    (void)rWeights;
    (void)total;
    throw std::invalid_argument(
        "partition: method 'kahip' requires a build with -DMESHIOPLUSPLUS_WITH_KAHIP=ON");
}

#endif  // MESHIOPLUSPLUS_HAS_KAHIP

// --- dispatch ----------------------------------------------------------------

// --- ghost (halo) layers -----------------------------------------------------
//
// The cell <-> node incidences live in `detail/cell_adjacency.hpp`, shared with
// `gradient`'s least-squares stencil so the two cannot disagree about what
// "shares a node" means. They used to be private to this file, and read dense
// connectivity as `Conn().As<std::int64_t>()` -- which performs no dtype check,
// so on a MESHIO-backed mesh carrying Int32 connectivity every node id was two
// fused Int32 entries. Most such ids failed the range filter below and simply
// vanished, leaving halos silently too small; the shared helper reads through
// `detail::read_int` instead.

// Validate the options and resolve Auto to a concrete method.
PartitionMethod partition_resolve_method(const PartitionOptions& rOptions) {
    if (rOptions.mNParts < 1)
        throw std::invalid_argument("partition: nparts must be >= 1, got " +
                                    std::to_string(rOptions.mNParts));
    if (rOptions.mGhostLayers < 0)
        throw std::invalid_argument("partition: ghost_layers must be >= 0, got " +
                                    std::to_string(rOptions.mGhostLayers));
    PartitionMethod method = rOptions.mMethod;
    if (method == PartitionMethod::Auto)
        method = partition_has_kahip() ? PartitionMethod::KaHIP : PartitionMethod::SFC;
    if (method == PartitionMethod::KaHIP &&
        !(rOptions.mImbalance > 0.0 && rOptions.mImbalance < 1.0))
        throw std::invalid_argument("partition: imbalance must be in (0, 1), got " +
                                    std::to_string(rOptions.mImbalance));
    return method;
}

// The per-global-cell part assignment (values in [0, nparts)).
std::vector<int> partition_compute_parts(const Mesh& rMesh, const PartitionOptions& rOptions) {
    const PartitionMethod method = partition_resolve_method(rOptions);
    const std::size_t total = partition_total_cells(rMesh);
    std::vector<double> weights;
    if (!rOptions.mWeightsKey.empty())
        weights = partition_weights(rMesh, rOptions.mWeightsKey, total);
    if (method == PartitionMethod::KaHIP)
        return partition_kahip_parts(rMesh, rOptions, weights, total);
    return partition_sfc_parts(rMesh, rOptions, weights, total);
}

}  // namespace

PartitionMethod partition_method_from_name(const std::string& rName) {
    if (rName == "sfc" || rName == "hilbert")
        return PartitionMethod::SFC;
    if (rName == "kahip")
        return PartitionMethod::KaHIP;
    if (rName == "auto" || rName.empty())
        return PartitionMethod::Auto;
    throw std::invalid_argument("partition: unknown method '" + rName +
                                "' (expected 'sfc', 'kahip', or 'auto')");
}

PartitionMode partition_mode_from_name(const std::string& rName) {
    if (rName == "fast")
        return PartitionMode::Fast;
    if (rName == "eco" || rName.empty())
        return PartitionMode::Eco;
    if (rName == "strong")
        return PartitionMode::Strong;
    throw std::invalid_argument("partition: unknown mode '" + rName +
                                "' (expected 'fast', 'eco', or 'strong')");
}

bool partition_has_kahip() noexcept {
#ifdef MESHIOPLUSPLUS_HAS_KAHIP
    return true;
#else
    return false;
#endif
}

std::vector<NDArray> partition_labels(const Mesh& rMesh, const PartitionOptions& rOptions) {
    // One label per input cell IS the ownership map; ghosting is a property of
    // the extracted pieces (a cell can be a ghost of several parts at once), so
    // it cannot be expressed here. Refuse rather than silently ignore it.
    if (rOptions.mGhostLayers != 0)
        throw std::invalid_argument(
            "partition_labels: ghost_layers has no meaning for a flat per-cell label array "
            "(a cell may be a ghost of several parts); use partition() for ghosted pieces");
    const std::vector<int> parts = partition_compute_parts(rMesh, rOptions);
    std::vector<NDArray> out;
    out.reserve(rMesh.NumCellBlocks());
    std::size_t base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {ncells});
        std::int64_t* dst = a.As<std::int64_t>();
        for (std::size_t c = 0; c < ncells; ++c)
            dst[c] = static_cast<std::int64_t>(parts[base + c]);
        out.push_back(std::move(a));
        base += ncells;
    }
    return out;
}

PartitionResult partition(const Mesh& rMesh, const PartitionOptions& rOptions) {
    const std::vector<int> parts = partition_compute_parts(rMesh, rOptions);
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::size_t nparts = static_cast<std::size_t>(rOptions.mNParts);

    const std::size_t total = parts.size();

    // Ghost (halo) layers: cells owned by ANOTHER part that share at least one
    // node with this part's cells, grown breadth-first `mGhostLayers` times.
    // Shared-node (rather than shared-facet) adjacency is the conservative
    // choice -- it is what a node-based assembly actually needs, and it is the
    // same neighbour definition the KaHIP dual graph falls back on.
    const int nghost = rOptions.mGhostLayers;
    detail::CellIncidence cell_nodes, node_cells;
    if (nghost > 0) {
        cell_nodes = detail::build_cell_node_incidence(rMesh, rMesh.NumPoints());
        node_cells = detail::invert_cell_node_incidence(cell_nodes, rMesh.NumPoints(), total);
    }

    // Per part, per block, the kept (ascending) local cell indices, and the
    // matching ghost depth (0 = owned) for the partition:ghost tag.
    std::vector<std::vector<std::vector<std::int64_t>>> kept(nparts);
    std::vector<std::vector<std::vector<std::int64_t>>> depth(nparts);
    for (std::size_t p = 0; p < nparts; ++p) {
        kept[p].assign(nblocks, {});
        depth[p].assign(nblocks, {});
    }

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    // stamp[c] == p means cell c is already in part p's set, so each cell is
    // visited once per part instead of needing a per-part boolean array.
    std::vector<int> stamp(total, -1);
    std::vector<int> cell_depth(total, 0);
    std::vector<std::int64_t> owned, frontier, next;

    for (std::size_t p = 0; p < nparts; ++p) {
        owned.clear();
        for (std::size_t c = 0; c < total; ++c)
            if (static_cast<std::size_t>(parts[c]) == p) {
                stamp[c] = static_cast<int>(p);
                cell_depth[c] = 0;
                owned.push_back(static_cast<std::int64_t>(c));
            }

        frontier = owned;
        for (int layer = 1; layer <= nghost && !frontier.empty(); ++layer) {
            next.clear();
            for (std::int64_t c : frontier)
                for (std::int64_t k = cell_nodes.mOffsets[static_cast<std::size_t>(c)];
                     k < cell_nodes.mOffsets[static_cast<std::size_t>(c) + 1]; ++k) {
                    const std::int64_t nid = cell_nodes.mEntries[static_cast<std::size_t>(k)];
                    if (nid < 0 || static_cast<std::size_t>(nid) >= rMesh.NumPoints())
                        continue;
                    for (std::int64_t j = node_cells.mOffsets[static_cast<std::size_t>(nid)];
                         j < node_cells.mOffsets[static_cast<std::size_t>(nid) + 1]; ++j) {
                        const std::int64_t c2 = node_cells.mEntries[static_cast<std::size_t>(j)];
                        if (stamp[static_cast<std::size_t>(c2)] == static_cast<int>(p))
                            continue;
                        stamp[static_cast<std::size_t>(c2)] = static_cast<int>(p);
                        cell_depth[static_cast<std::size_t>(c2)] = layer;
                        owned.push_back(c2);
                        next.push_back(c2);
                    }
                }
            frontier.swap(next);
        }

        // Ascending global order -> ascending (block, local) order, which is
        // what build_cell_subset expects.
        std::sort(owned.begin(), owned.end());
        for (std::int64_t g : owned) {
            const auto [blk, row] = detail::global_to_block_row(bases, g);
            kept[p][blk].push_back(row);
            depth[p][blk].push_back(cell_depth[static_cast<std::size_t>(g)]);
        }
    }

    // Exactly nparts pieces (possibly empty), input block structure kept 1:1
    // (drop_empty_blocks=false -- deliberately unlike split -- so mCellMaps
    // indexes by input block and the pieces recombine cleanly).
    PartitionResult res;
    res.mPieces.reserve(nparts);
    for (std::size_t p = 0; p < nparts; ++p) {
        detail::SubsetResult sub = detail::build_cell_subset(
            rMesh, kept[p], rOptions.mRecordIds ? "partition:original_point_id" : "",
            rOptions.mRecordIds ? "partition:original_cell_id" : "",
            /*drop_empty_blocks=*/false, "partition");
        PartitionPiece piece;
        piece.mPartId = static_cast<int>(p);
        piece.mMesh = std::move(sub.mMesh);
        piece.mPointMap = std::move(sub.mPointMap);
        piece.mCellMaps = std::move(sub.mCellMaps);
        if (nghost > 0) {
            // 0 = owned by this part, L = reached at BFS layer L. Emitted for
            // every block so the array stays block-aligned (the cell_data
            // invariant), including the empty ones.
            for (std::size_t blk = 0; blk < nblocks; ++blk) {
                NDArray a = NDArray::Uninit(DType::Int64, {depth[p][blk].size()});
                std::int64_t* dst = a.As<std::int64_t>();
                for (std::size_t i = 0; i < depth[p][blk].size(); ++i)
                    dst[i] = depth[p][blk][i];
                piece.mMesh.AppendCellData("partition:ghost", std::move(a));
            }
        }
        res.mPieces.push_back(std::move(piece));
    }
    return res;
}

}  // namespace meshioplusplus
