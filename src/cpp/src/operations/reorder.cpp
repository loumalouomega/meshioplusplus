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
// Mesh renumbering (RCM + Morton/Hilbert space-filling curves). A pure
// permutation of node and element indices — geometry and all data preserved —
// built entirely through the uniform mesh API so it compiles under every mesh
// backend. See operations/reorder.hpp for the contract.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/detail/space_filling.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// --- small helpers ----------------------------------------------------------

// A deep copy of an NDArray that always owns its buffer (the source may be a
// view over foreign memory, e.g. numpy memory on the write path).
NDArray reorder_owned_copy(const NDArray& rArr) {
    NDArray c = rArr;
    c.MakeOwned();
    return c;
}

// --- CSR node adjacency -----------------------------------------------------

// The node graph and its builder live in detail/node_adjacency.hpp, shared with
// operations/smooth.cpp. Reorder wants the *clique* graph — every pair of nodes
// co-occurring in a cell — because that is the sparsity pattern of the FEM
// matrix whose bandwidth this operation exists to reduce. (Smooth wants the
// edge graph instead; see that header.) The Clique path is the code that used
// to live here, unchanged, so RCM's tie-breaking on adjacency row order — and
// therefore every permutation this operation has ever produced — is unaffected.
using ReorderCsr = detail::NodeAdjacency;

// --- Reverse Cuthill-McKee --------------------------------------------------

// Build the rooted level structure from `root` (BFS), returning its eccentricity
// (deepest level) and filling `rLastLevel` with the deepest-level nodes. The
// scratch `rLevel` (size n, all -1 on entry) is reset to all -1 on return.
std::int64_t reorder_rls(const ReorderCsr& rCsr, std::int64_t root,
                         std::vector<std::int64_t>& rLevel, std::vector<std::int64_t>& rLastLevel) {
    std::vector<std::int64_t> touched;
    std::queue<std::int64_t> q;
    rLevel[static_cast<std::size_t>(root)] = 0;
    touched.push_back(root);
    q.push(root);
    std::int64_t ecc = 0;
    while (!q.empty()) {
        const std::int64_t u = q.front();
        q.pop();
        const std::int64_t lu = rLevel[static_cast<std::size_t>(u)];
        if (lu > ecc)
            ecc = lu;
        for (std::int64_t k = rCsr.mXadj[static_cast<std::size_t>(u)];
             k < rCsr.mXadj[static_cast<std::size_t>(u) + 1]; ++k) {
            const std::int64_t v = rCsr.mAdj[static_cast<std::size_t>(k)];
            if (rLevel[static_cast<std::size_t>(v)] == -1) {
                rLevel[static_cast<std::size_t>(v)] = lu + 1;
                touched.push_back(v);
                q.push(v);
            }
        }
    }
    rLastLevel.clear();
    for (std::int64_t t : touched)
        if (rLevel[static_cast<std::size_t>(t)] == ecc)
            rLastLevel.push_back(t);
    std::sort(rLastLevel.begin(), rLastLevel.end());
    for (std::int64_t t : touched)
        rLevel[static_cast<std::size_t>(t)] = -1;  // reset scratch
    return ecc;
}

// George-Liu pseudo-peripheral node heuristic: iterate the rooted level
// structure, hopping to a minimum-degree node of the deepest level until the
// eccentricity stops growing.
std::int64_t reorder_pseudo_peripheral(const ReorderCsr& rCsr, std::int64_t start,
                                       const std::vector<std::int64_t>& rDegree,
                                       std::vector<std::int64_t>& rLevel) {
    std::int64_t root = start;
    std::vector<std::int64_t> last;
    std::int64_t ecc = reorder_rls(rCsr, root, rLevel, last);
    for (int it = 0; it < 50; ++it) {  // guard against pathological non-convergence
        std::int64_t cand = last.empty() ? root : last[0];
        std::int64_t best_deg = rDegree[static_cast<std::size_t>(cand)];
        for (std::int64_t v : last) {
            if (rDegree[static_cast<std::size_t>(v)] < best_deg) {
                best_deg = rDegree[static_cast<std::size_t>(v)];
                cand = v;
            }
        }
        std::vector<std::int64_t> last2;
        const std::int64_t ecc2 = reorder_rls(rCsr, cand, rLevel, last2);
        if (ecc2 > ecc) {
            root = cand;
            ecc = ecc2;
            last = std::move(last2);
        } else {
            break;
        }
    }
    return root;
}

// Reverse Cuthill-McKee permutation (old index -> new index) over the CSR graph.
// Handles disconnected components (restart from the lowest-id unvisited node)
// and isolated nodes (their own singleton component).
std::vector<std::int64_t> reorder_rcm(const ReorderCsr& rCsr, std::size_t n) {
    std::vector<std::int64_t> degree(n);
    for (std::size_t i = 0; i < n; ++i)
        degree[i] = rCsr.mXadj[i + 1] - rCsr.mXadj[i];

    std::vector<char> visited(n, 0);
    std::vector<std::int64_t> order;  // Cuthill-McKee order (new position -> old node)
    order.reserve(n);
    std::vector<std::int64_t> level(n, -1);  // scratch for the RLS
    std::vector<std::int64_t> nbrs;

    for (std::size_t s = 0; s < n; ++s) {
        if (visited[s])
            continue;
        const std::int64_t root =
            reorder_pseudo_peripheral(rCsr, static_cast<std::int64_t>(s), degree, level);
        std::queue<std::int64_t> q;
        visited[static_cast<std::size_t>(root)] = 1;
        q.push(root);
        while (!q.empty()) {
            const std::int64_t u = q.front();
            q.pop();
            order.push_back(u);
            nbrs.clear();
            for (std::int64_t k = rCsr.mXadj[static_cast<std::size_t>(u)];
                 k < rCsr.mXadj[static_cast<std::size_t>(u) + 1]; ++k) {
                const std::int64_t v = rCsr.mAdj[static_cast<std::size_t>(k)];
                if (!visited[static_cast<std::size_t>(v)])
                    nbrs.push_back(v);
            }
            std::sort(nbrs.begin(), nbrs.end(), [&](std::int64_t a, std::int64_t b) {
                if (degree[static_cast<std::size_t>(a)] != degree[static_cast<std::size_t>(b)])
                    return degree[static_cast<std::size_t>(a)] <
                           degree[static_cast<std::size_t>(b)];
                return a < b;
            });
            for (std::int64_t v : nbrs) {
                if (!visited[static_cast<std::size_t>(v)]) {
                    visited[static_cast<std::size_t>(v)] = 1;
                    q.push(v);
                }
            }
        }
    }

    std::reverse(order.begin(), order.end());  // Cuthill-McKee -> reverse
    std::vector<std::int64_t> perm(n);
    for (std::size_t newidx = 0; newidx < order.size(); ++newidx)
        perm[static_cast<std::size_t>(order[newidx])] = static_cast<std::int64_t>(newidx);
    return perm;
}

// --- space-filling curves ---------------------------------------------------

// The key transforms (part1by2 / Morton / Hilbert) live in
// detail/space_filling.hpp, shared with the partition operation. Only the
// Mesh-coupled quantization below is reorder-specific.
constexpr int REORDER_SFC_BITS = detail::sfc_bits;  // 3 * 21 = 63 bits fit in a uint64 key

// Space-filling-curve keys over the node coordinates (bounding box quantized to
// REORDER_SFC_BITS per axis). `hilbert` selects Hilbert vs Morton.
std::vector<std::uint64_t> reorder_sfc_keys(const Mesh& rMesh, std::size_t n, bool hilbert) {
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    const std::size_t dim = std::min<std::size_t>(pdim, 3);
    const int bits = REORDER_SFC_BITS;

    double lo[3] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[3] = {-std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity()};
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t d = 0; d < dim; ++d) {
            const double c = detail::read_double(points, i * pdim + d);
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

    std::vector<std::uint64_t> keys(n);
    parallel_for(n, [&](std::size_t i) {
        std::uint32_t q[3] = {0, 0, 0};
        for (std::size_t d = 0; d < dim; ++d) {
            const double c = detail::read_double(points, i * pdim + d);
            double t = (c - lo[d]) * scale[d];
            std::int64_t qi = static_cast<std::int64_t>(t + 0.5);
            if (qi < 0)
                qi = 0;
            const std::int64_t maxq = (std::int64_t(1) << bits) - 1;
            if (qi > maxq)
                qi = maxq;
            q[d] = static_cast<std::uint32_t>(qi);
        }
        keys[i] = hilbert ? detail::sfc_hilbert_key(q, bits) : detail::sfc_morton_key(q);
    });
    return keys;
}

// Stable argsort of the SFC keys -> node permutation (old index -> new index).
std::vector<std::int64_t> reorder_sfc_perm(const std::vector<std::uint64_t>& rKeys, std::size_t n) {
    std::vector<std::int64_t> order(n);
    std::iota(order.begin(), order.end(), std::int64_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::int64_t a, std::int64_t b) {
        return rKeys[static_cast<std::size_t>(a)] < rKeys[static_cast<std::size_t>(b)];
    });
    std::vector<std::int64_t> perm(n);
    for (std::size_t newidx = 0; newidx < n; ++newidx)
        perm[static_cast<std::size_t>(order[newidx])] = static_cast<std::int64_t>(newidx);
    return perm;
}

// --- applying a permutation -------------------------------------------------

// Scatter the `nrows` rows of `rSrc` so that new row `perm[i]` gets old row `i`,
// preserving dtype and column shape (a bandwidth-bound byte gather).
NDArray reorder_scatter_rows(const NDArray& rSrc, const std::vector<std::int64_t>& rPerm,
                             std::size_t nrows) {
    NDArray dst = NDArray::Uninit(rSrc.Dtype(), rSrc.Shape());
    if (nrows == 0)
        return dst;
    const std::size_t row_bytes = rSrc.Nbytes() / nrows;
    const std::byte* src = rSrc.Data();
    std::byte* out = dst.Data();
    parallel_for_bw(nrows, [&](std::size_t i) {
        std::memcpy(out + static_cast<std::size_t>(rPerm[i]) * row_bytes, src + i * row_bytes,
                    row_bytes);
    });
    return dst;
}

// Per-cell key = the minimum new node index among the cell's nodes.
std::int64_t reorder_cell_key(const std::vector<std::int64_t>& rNodes,
                              const std::vector<std::int64_t>& rNodePerm) {
    std::int64_t m = std::numeric_limits<std::int64_t>::max();
    for (std::int64_t id : rNodes)
        m = std::min(m, rNodePerm[static_cast<std::size_t>(id)]);
    return rNodes.empty() ? 0 : m;
}

ReorderResult reorder_apply(const Mesh& rMesh, std::vector<std::int64_t> node_perm) {
    const std::size_t n = rMesh.NumPoints();
    ReorderResult res;
    Mesh& out = res.mMesh;

    // Node permutation array (old -> new).
    {
        NDArray np = NDArray::Uninit(DType::Int64, {n});
        std::int64_t* d = np.As<std::int64_t>();
        for (std::size_t i = 0; i < n; ++i)
            d[i] = node_perm[i];
        res.mNodePermutation = std::move(np);
    }

    // Points scattered by the node permutation.
    out.AssignPoints(reorder_scatter_rows(rMesh.Points(), node_perm, n));

    // Cell blocks: remap connectivity node ids, then reorder cells by min new
    // node index. Keep each block's cell permutation for the cell_data below.
    std::vector<std::vector<std::int64_t>> block_cell_perms;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t nc = cb.NumCells();

        // Per-cell sort key = min new node index; stable argsort -> cell order.
        std::vector<std::int64_t> key(nc);
        parallel_for_bw(nc, [&](std::size_t c) {
            std::vector<std::int64_t> nodes;
            detail::cell_node_ids(cb, c, n, nodes);
            key[c] = reorder_cell_key(nodes, node_perm);
        });
        std::vector<std::int64_t> cellorder(nc);
        std::iota(cellorder.begin(), cellorder.end(), std::int64_t{0});
        std::stable_sort(cellorder.begin(), cellorder.end(), [&](std::int64_t a, std::int64_t b) {
            return key[static_cast<std::size_t>(a)] < key[static_cast<std::size_t>(b)];
        });
        std::vector<std::int64_t> cellperm(nc);  // old -> new
        for (std::size_t p = 0; p < nc; ++p)
            cellperm[static_cast<std::size_t>(cellorder[p])] = static_cast<std::int64_t>(p);

        auto remap = [&](std::int64_t id) -> std::int64_t {
            return (id >= 0 && static_cast<std::size_t>(id) < n)
                       ? node_perm[static_cast<std::size_t>(id)]
                       : id;
        };
        if (cb.IsPolyhedron()) {
            std::vector<std::vector<std::vector<std::int64_t>>> cells(nc);
            parallel_for_bw(nc, [&](std::size_t p) {
                const std::size_t oc = static_cast<std::size_t>(cellorder[p]);
                cells[p].resize(cb.NumFaces(oc));
                for (std::size_t f = 0; f < cb.NumFaces(oc); ++f) {
                    std::pair<const std::int64_t*, std::size_t> face = cb.Face(oc, f);
                    cells[p][f].reserve(face.second);
                    for (std::size_t k = 0; k < face.second; ++k)
                        cells[p][f].push_back(remap(face.first[k]));
                }
            });
            out.AddPolyhedronBlock(std::string(cb.Type()), std::move(cells));
        } else if (cb.IsRagged()) {
            std::vector<std::vector<std::int64_t>> rows(nc);
            parallel_for_bw(nc, [&](std::size_t p) {
                const std::size_t oc = static_cast<std::size_t>(cellorder[p]);
                const std::int64_t* row = cb.Row(oc);
                const std::size_t sz = cb.RowSize(oc);
                rows[p].reserve(sz);
                for (std::size_t k = 0; k < sz; ++k)
                    rows[p].push_back(remap(row[k]));
            });
            out.AddPolygonBlock(std::string(cb.Type()), std::move(rows));
        } else {
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            NDArray block = NDArray::Uninit(DType::Int64, {nc, npc});
            std::int64_t* dst = block.As<std::int64_t>();
            parallel_for_bw(nc, [&](std::size_t p) {
                const std::size_t oc = static_cast<std::size_t>(cellorder[p]);
                for (std::size_t k = 0; k < npc; ++k)
                    dst[p * npc + k] = remap(detail::read_int(conn, oc * npc + k));
            });
            out.AddCellBlock(std::string(cb.Type()), std::move(block));
        }

        NDArray cp = NDArray::Uninit(DType::Int64, {nc});
        std::memcpy(cp.Data(), cellperm.data(), nc * sizeof(std::int64_t));
        res.mCellPermutations.push_back(std::move(cp));
        block_cell_perms.push_back(std::move(cellperm));
    }

    // point_data: scatter rows that align with the point table; copy the rest.
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (detail::rows(a) == n)
            out.AddPointData(name, reorder_scatter_rows(a, node_perm, n));
        else
            out.AddPointData(name, reorder_owned_copy(a));
    }

    // cell_data: reorder each block's rows by that block's cell permutation.
    for (const std::string& name : rMesh.CellDataNames()) {
        std::vector<NDArray> blocks;
        const std::size_t nb = rMesh.CellDataNumBlocks(name);
        blocks.reserve(nb);
        for (std::size_t b = 0; b < nb; ++b) {
            const NDArray& a = rMesh.CellData(name, b);
            if (b < block_cell_perms.size() && detail::rows(a) == block_cell_perms[b].size())
                blocks.push_back(reorder_scatter_rows(a, block_cell_perms[b], detail::rows(a)));
            else
                blocks.push_back(reorder_owned_copy(a));
        }
        out.AddCellData(name, std::move(blocks));
    }

    // field_data: mesh-global, unaffected by the permutation.
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, reorder_owned_copy(rMesh.FieldData(name)));

    // named regions: a pure permutation, so every entry survives — and because
    // the cell type of each block is unchanged, side-set facets survive too.
    // The maps are Direct (permutations are not monotone, so the FirstChild
    // "next entry" rule would be wrong here).
    {
        detail::RegionRemap rmap;
        rmap.pPointMap = &res.mNodePermutation;
        rmap.mCellMapKind = detail::CellMapKind::Direct;
        rmap.pCellMaps = &res.mCellPermutations;
        rmap.mOpName = "reorder";
        detail::remap_regions(rMesh, out, rmap);
    }

    return res;
}

}  // namespace

// --- public API -------------------------------------------------------------

ReorderMethod reorder_method_from_name(const std::string& rName) {
    if (rName == "rcm" || rName == "RCM" || rName == "Rcm")
        return ReorderMethod::RCM;
    if (rName == "morton" || rName == "z" || rName == "zorder" || rName == "z-order")
        return ReorderMethod::Morton;
    if (rName == "hilbert")
        return ReorderMethod::Hilbert;
    throw std::invalid_argument("reorder: unknown method '" + rName +
                                "' (expected 'rcm', 'morton', or 'hilbert')");
}

std::int64_t compute_bandwidth(const Mesh& rMesh) {
    std::int64_t bw = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t nc = cb.NumCells();
        if (nc == 0)
            continue;
        const bool poly = cb.IsPolyhedron();
        const bool ragged = !poly && cb.IsRagged();
        const NDArray* conn = (!poly && !ragged) ? &cb.Conn() : nullptr;
        const std::size_t npc = cb.NodesPerCell();
        std::vector<std::int64_t> local(nc, 0);
        parallel_for(nc, [&](std::size_t c) {
            std::int64_t lo = std::numeric_limits<std::int64_t>::max();
            std::int64_t hi = std::numeric_limits<std::int64_t>::min();
            auto consider = [&](std::int64_t id) {
                lo = std::min(lo, id);
                hi = std::max(hi, id);
            };
            if (poly) {
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    std::pair<const std::int64_t*, std::size_t> face = cb.Face(c, f);
                    for (std::size_t k = 0; k < face.second; ++k)
                        consider(face.first[k]);
                }
            } else if (ragged) {
                const std::int64_t* row = cb.Row(c);
                for (std::size_t k = 0; k < cb.RowSize(c); ++k)
                    consider(row[k]);
            } else {
                for (std::size_t k = 0; k < npc; ++k)
                    consider(detail::read_int(*conn, c * npc + k));
            }
            local[c] = (hi >= lo) ? (hi - lo) : 0;
        });
        for (std::int64_t v : local)
            bw = std::max(bw, v);
    }
    return bw;
}

ReorderResult reorder(const Mesh& rMesh, ReorderMethod method) {
    const std::size_t n = rMesh.NumPoints();
    std::vector<std::int64_t> perm;
    if (method == ReorderMethod::RCM) {
        ReorderCsr csr = detail::build_node_adjacency(rMesh, n, detail::NodeAdjacencyKind::Clique);
        perm = reorder_rcm(csr, n);
    } else {
        std::vector<std::uint64_t> keys =
            reorder_sfc_keys(rMesh, n, method == ReorderMethod::Hilbert);
        perm = reorder_sfc_perm(keys, n);
    }
    return reorder_apply(rMesh, std::move(perm));
}

}  // namespace meshioplusplus
