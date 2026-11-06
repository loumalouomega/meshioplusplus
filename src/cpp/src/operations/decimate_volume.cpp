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
// Quadric-error-metric tet-edge collapse volume decimation. See
// operations/decimate_volume.hpp for the contract. Built entirely on the
// uniform mesh API plus the placement/error/link-condition machinery hoisted
// into detail/decimate_common.hpp for `decimate.cpp` (2D) to share.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/decimate_volume.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/decimate_common.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/decimate_shared.hpp"

namespace meshioplusplus {

namespace {

using detail::decim_accumulate_quadrics;
using detail::decim_count_common;
using detail::decim_face_planes;
using detail::decim_place;
using detail::decim_quadric_error;
using detail::decim_vertex_faces_csr;
using detail::decim_vertex_ring;
using detail::DecimCsr;
using detail::DecimFaces;
using detail::DecimPlaceCtx;
using detail::DecimPlaced;
using detail::Vec3;

constexpr const char* kDvPrefix = "meshio++: decimate_volume: ";

// --- validation --------------------------------------------------------------

struct DvStop {
    bool mUseTarget = false;
    bool mUseMaxError = false;
    std::int64_t mTargetCells = -1;
    double mMaxError = -1.0;
};

DvStop dv_resolve_stop(const DecimateVolumeOptions& rOptions) {
    const int num_set = (rOptions.mTargetRatio >= 0.0 ? 1 : 0) +
                        (rOptions.mTargetCells >= 0 ? 1 : 0) + (rOptions.mMaxError >= 0.0 ? 1 : 0);
    if (num_set != 1)
        throw std::invalid_argument(
            std::string(kDvPrefix) +
            "exactly one of target_ratio, target_cells, max_error must be set; got " +
            std::to_string(num_set));
    DvStop stop;
    if (rOptions.mTargetRatio >= 0.0) {
        if (!(rOptions.mTargetRatio > 0.0 && rOptions.mTargetRatio <= 1.0))
            throw std::invalid_argument(std::string(kDvPrefix) +
                                        "target_ratio must lie in (0, 1]; got " +
                                        std::to_string(rOptions.mTargetRatio));
        stop.mUseTarget = true;  // filled in once T (tet count) is known
    } else if (rOptions.mTargetCells >= 0) {
        stop.mUseTarget = true;
        stop.mTargetCells = rOptions.mTargetCells;
    } else {
        stop.mUseMaxError = true;
        stop.mMaxError = rOptions.mMaxError;
    }
    return stop;
}

// Reject every construct outside the tet-only scope, by name, before any
// other work.
void dv_check_blocks(const Mesh& rMesh) {
    bool has_tet = false;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                std::string(kDvPrefix) +
                "mesh contains a polyhedron cell block; decimate_volume operates on tet-only "
                "meshes -- run convert_cells(simplexify) first");
        if (cb.IsRagged())
            throw std::invalid_argument(
                std::string(kDvPrefix) + "mesh contains ragged cell block '" +
                std::string(cb.Type()) + "'; decimate_volume operates on tet-only meshes");
        const CellType ct = cell_type_from_name(std::string(cb.Type()));
        if (ct == CellType::Tetra) {
            has_tet = true;
            continue;
        }
        const int dim = cell_type_dimension(ct);
        if (dim == 3)
            throw std::invalid_argument(
                std::string(kDvPrefix) + "mesh contains 3D volume cell block '" +
                std::string(cb.Type()) +
                "' that is not linear tetra; decimate_volume operates on tet-only meshes -- run "
                "convert_cells(simplexify) first");
        throw std::invalid_argument(
            std::string(kDvPrefix) + "mesh contains non-3D cell block '" + std::string(cb.Type()) +
            "' alongside its tets (its nodes would dangle after the collapse; drop it first, "
            "e.g. via split)");
    }
    if (!has_tet)
        throw std::invalid_argument(std::string(kDvPrefix) + "mesh contains no tetra cell block");
}

// --- small utilities ----------------------------------------------------------

// Sorted intersection of two sorted vectors (values present in both).
std::vector<std::int64_t> dv_sorted_intersection(const std::vector<std::int64_t>& rA,
                                                 const std::vector<std::int64_t>& rB) {
    std::vector<std::int64_t> out;
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < rA.size() && j < rB.size()) {
        if (rA[i] < rB[j])
            ++i;
        else if (rB[j] < rA[i])
            ++j;
        else {
            out.push_back(rA[i]);
            ++i;
            ++j;
        }
    }
    return out;
}

// The "opposite triple" of tet `t` after removing corner `v` (the tet's other
// 3 corners, sorted). `t` must contain `v`.
std::array<std::int64_t, 3> dv_opposite_triple(const std::vector<std::int64_t>& rTetConn,
                                               std::int64_t t, std::int64_t v) {
    std::array<std::int64_t, 3> tri{-1, -1, -1};
    const std::int64_t* c = rTetConn.data() + static_cast<std::size_t>(t) * 4;
    std::size_t w = 0;
    for (int k = 0; k < 4; ++k)
        if (c[k] != v)
            tri[w++] = c[k];
    std::sort(tri.begin(), tri.end());
    return tri;
}

// `before` non-zero and `after` on the opposite side of zero (or exactly
// zero): a sign flip or a degenerate result. An already-degenerate `before`
// (exactly zero) imposes no constraint -- `smooth`'s "do no harm" rule.
bool dv_tet_inverts(double before, double after) {
    if (before == 0.0)
        return false;
    if (before > 0.0)
        return after <= 0.0;
    return after >= 0.0;
}

// Signed volume of tet `t`, substituting any corner equal to `A` or `B` with
// `pSub` when given.
double dv_tet_volume(const std::vector<double>& rXyz, const std::vector<std::int64_t>& rTetConn,
                     std::int64_t t, std::int64_t A, std::int64_t B, const double* pSub) {
    const std::int64_t* c = rTetConn.data() + static_cast<std::size_t>(t) * 4;
    Vec3 coords[4];
    for (int k = 0; k < 4; ++k) {
        const std::int64_t id = c[k];
        if (pSub != nullptr && (id == A || id == B)) {
            coords[k] = {pSub[0], pSub[1], pSub[2]};
        } else {
            const double* src = rXyz.data() + static_cast<std::size_t>(id) * 3;
            coords[k] = {src[0], src[1], src[2]};
        }
    }
    return detail::cell_volume_from_corners(coords, CellType::Tetra);
}

// A Float64 working copy of a float-kind point_data array, blended in place
// at each collapse commit (shared with decimate.cpp, detail/decimate_shared.hpp).
using detail::DecimFloatData;

}  // namespace

DecimateVolumeResult decimate_volume(const Mesh& rMesh, const DecimateVolumeOptions& rOptions) {
    DvStop stop = dv_resolve_stop(rOptions);

    const std::size_t n = rMesh.NumPoints();
    if (!rOptions.mFrozen.empty() && rOptions.mFrozen.size() != n)
        throw std::invalid_argument(std::string(kDvPrefix) + "frozen mask has " +
                                    std::to_string(rOptions.mFrozen.size()) +
                                    " entries but the mesh has " + std::to_string(n) + " points");

    dv_check_blocks(rMesh);

    const std::size_t dim = rMesh.PointDim();
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t num_tets = static_cast<std::size_t>(detail::total_cells(bases));

    // --- phase 1: flat buffers (parallel fills) -------------------------------
    std::vector<double> xyz(n * 3, 0.0);
    {
        const NDArray& points = rMesh.Points();
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t d = 0; d < dim && d < 3; ++d)
                xyz[i * 3 + d] = detail::read_double(points, i * dim + d);
        });
    }

    std::vector<std::int64_t> tet_conn(num_tets * 4);
    {
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const NDArray& conn = cb.Conn();
            const std::size_t nc = cb.NumCells();
            const std::size_t base = static_cast<std::size_t>(bases[bi]);
            for (std::size_t c = 0; c < nc; ++c)
                for (std::size_t k = 0; k < 4; ++k)
                    tet_conn[(base + c) * 4 + k] = detail::read_int(conn, c * 4 + k);
            ++bi;
        }
    }

    // --- phase 2: the mesh's own outer skin -----------------------------------
    const detail::GlobalFaces gf = detail::build_global_faces(rMesh);
    if (gf.mNumNonManifold > 0)
        throw std::invalid_argument(
            std::string(kDvPrefix) + "mesh contains " + std::to_string(gf.mNumNonManifold) +
            " face(s) shared by three or more tets (non-manifold); refusing rather than "
            "guessing a boundary classification");

    DecimFaces boundary;
    boundary.mBlockBase.push_back(0);
    for (std::size_t f = 0; f < gf.NumFaces(); ++f) {
        if (gf.mNeighbour[f] >= 0)
            continue;  // interior face
        const std::int64_t* ring = gf.Face(f);
        const std::size_t fs = gf.FaceSize(f);
        // A tet's faces are triangles; a non-triangular boundary face would
        // mean a non-tet block contributed it, which dv_check_blocks already
        // refused.
        for (std::size_t k = 0; k < fs && k < 3; ++k)
            boundary.mCorners.push_back(ring[k]);
        ++boundary.mNumFaces;
    }

    const DecimCsr bcsr = decim_vertex_faces_csr(boundary, n);
    std::vector<std::uint8_t> touches_boundary(n, 0);
    for (std::size_t v = 0; v < n; ++v)
        touches_boundary[v] = bcsr.mXadj[v + 1] > bcsr.mXadj[v] ? 1 : 0;

    std::vector<double> bquad_k;
    std::vector<double> bnormals;
    decim_face_planes(boundary, xyz, bquad_k, bnormals);

    // --- phase 3: the pin mask (caller | boundary | feature) -----------------
    const std::vector<std::uint8_t> pinned = detail::decim_pin_mask(
        rOptions.mFrozen, touches_boundary, rOptions.mPreserveBoundary, rOptions.mPreserveFeatures,
        rOptions.mFeatureAngleDeg, boundary.mCorners, bnormals, n);

    // --- phase 4: per-vertex quadrics (boundary-only, fixed FP order) --------
    std::vector<double> quads = decim_accumulate_quadrics(bcsr, n, bquad_k);
    bquad_k.clear();
    bquad_k.shrink_to_fit();

    // --- phase 5: mutable collapse state ---------------------------------------
    std::vector<std::uint8_t> tet_alive(num_tets, 1);
    std::vector<std::int64_t> successor(n, -1);  // removed vertex -> survivor
    std::vector<std::uint32_t> version(n, 0);

    // vtets[v]: sorted, incident tet ids -- built by hand (no shared CSR
    // builder: decim_vertex_faces_csr is hardcoded to 3 corners per face).
    std::vector<std::vector<std::int64_t>> vtets(n);
    {
        std::vector<std::int64_t> counts(n, 0);
        for (std::size_t t = 0; t < num_tets; ++t)
            for (int k = 0; k < 4; ++k)
                ++counts[static_cast<std::size_t>(tet_conn[t * 4 + k])];
        for (std::size_t v = 0; v < n; ++v)
            vtets[v].reserve(static_cast<std::size_t>(counts[v]));
        for (std::size_t t = 0; t < num_tets; ++t)
            for (int k = 0; k < 4; ++k)
                vtets[static_cast<std::size_t>(tet_conn[t * 4 + k])].push_back(
                    static_cast<std::int64_t>(t));
        // Each tet contributes each of its 4 corners once, in tet-ascending
        // order, so every row is already sorted -- no explicit sort needed.
    }

    // vboundary[v]: sorted, incident ALIVE boundary-face ids into `boundary`
    // (mutated in lockstep with it, exactly `decimate.cpp`'s own `vfaces`
    // applied to the mesh's own skin instead of a triangulated 2D input) --
    // what lets a boundary-touching collapse reuse `decimate`'s own 2D
    // ring/shared-face link condition and normal-flip check.
    std::vector<std::uint8_t> boundary_alive(boundary.mNumFaces, 1);
    std::vector<std::vector<std::int64_t>> vboundary(n);
    for (std::size_t v = 0; v < n; ++v) {
        std::vector<std::int64_t>& row = vboundary[v];
        row.reserve(static_cast<std::size_t>(bcsr.mXadj[v + 1] - bcsr.mXadj[v]));
        for (std::int64_t k = bcsr.mXadj[v]; k < bcsr.mXadj[v + 1]; ++k)
            row.push_back(bcsr.mAdj[static_cast<std::size_t>(k)]);
    }

    std::vector<DecimFloatData> float_data = detail::decim_float_working_copies(rMesh, n);

    DecimPlaceCtx ctx;
    ctx.mpXyz = &xyz;
    ctx.mpQ = &quads;
    ctx.mpPinned = &pinned;
    ctx.mPlacement = rOptions.mPlacement;

    if (stop.mUseTarget && rOptions.mTargetRatio >= 0.0)
        stop.mTargetCells = std::max<std::int64_t>(
            1,
            static_cast<std::int64_t>(rOptions.mTargetRatio * static_cast<double>(num_tets) + 0.5));

    // --- phase 6: unique tet edges (6 per tet, dedup by first occurrence) ----
    using DvEdgeKey = std::array<std::int64_t, 2>;
    struct DvEdgeKeyHash {
        std::size_t operator()(const DvEdgeKey& rKey) const {
            std::size_t h = 0;
            for (std::int64_t v : rKey)
                h ^= std::hash<std::int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            return h;
        }
    };
    static constexpr int kTetEdgeCorners[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    std::vector<DvEdgeKey> edges;
    {
        std::unordered_map<DvEdgeKey, char, DvEdgeKeyHash> seen;
        seen.reserve(num_tets * 6);
        for (std::size_t t = 0; t < num_tets; ++t) {
            const std::int64_t* c = tet_conn.data() + t * 4;
            for (const auto& ec : kTetEdgeCorners) {
                std::int64_t a = c[ec[0]];
                std::int64_t b = c[ec[1]];
                if (a > b)
                    std::swap(a, b);
                DvEdgeKey key{a, b};
                if (seen.emplace(key, 0).second)
                    edges.push_back(key);
            }
        }
    }

    // --- phase 7: seed the queue (parallel scoring, serial pushes) -----------
    // Priority: (regime, score, lo, hi, version_lo, version_hi). regime 0
    // (boundary-touching, real quadric error) sorts ahead of regime 1 (purely
    // interior, squared edge length) -- see the header doc comment.
    using DvHeapEntry =
        std::tuple<int, double, std::int64_t, std::int64_t, std::uint32_t, std::uint32_t>;
    std::priority_queue<DvHeapEntry, std::vector<DvHeapEntry>, std::greater<DvHeapEntry>> heap;

    auto dv_score = [&](std::int64_t a, std::int64_t b) -> std::pair<int, double> {
        if (touches_boundary[static_cast<std::size_t>(a)] ||
            touches_boundary[static_cast<std::size_t>(b)]) {
            return {0, decim_place(ctx, a, b).mErr};
        }
        const double* xa = xyz.data() + static_cast<std::size_t>(a) * 3;
        const double* xb = xyz.data() + static_cast<std::size_t>(b) * 3;
        const double dx = xa[0] - xb[0];
        const double dy = xa[1] - xb[1];
        const double dz = xa[2] - xb[2];
        return {1, dx * dx + dy * dy + dz * dz};
    };

    {
        std::vector<int> seed_regime(edges.size(), 0);
        std::vector<double> seed_score(edges.size(), 0.0);
        parallel_for(edges.size(), [&](std::size_t e) {
            if (pinned[static_cast<std::size_t>(edges[e][0])] &&
                pinned[static_cast<std::size_t>(edges[e][1])])
                return;
            const auto rs = dv_score(edges[e][0], edges[e][1]);
            seed_regime[e] = rs.first;
            seed_score[e] = rs.second;
        });
        for (std::size_t e = 0; e < edges.size(); ++e) {
            const std::int64_t a = edges[e][0];
            const std::int64_t b = edges[e][1];
            if (pinned[static_cast<std::size_t>(a)] && pinned[static_cast<std::size_t>(b)])
                continue;
            if (!std::isfinite(seed_score[e]))
                continue;
            heap.push({seed_regime[e], seed_score[e], a, b, 0, 0});
        }
    }

    // --- phase 8: the SERIAL greedy loop --------------------------------------
    DecimateVolumeResult result;
    std::int64_t alive_tets = static_cast<std::int64_t>(num_tets);
    while (true) {
        if (stop.mUseTarget && alive_tets <= stop.mTargetCells)
            break;
        if (heap.empty()) {
            if (stop.mUseTarget && alive_tets > stop.mTargetCells)
                log::warn(
                    "decimate_volume: queue exhausted at {} tets before reaching the target of "
                    "{} (pinning and the validity guards left no collapsible edge)",
                    alive_tets, stop.mTargetCells);
            break;
        }
        const DvHeapEntry top = heap.top();
        heap.pop();
        const int regime = std::get<0>(top);
        const double score = std::get<1>(top);
        const std::int64_t a = std::get<2>(top);
        const std::int64_t b = std::get<3>(top);
        if (regime == 0 && stop.mUseMaxError && score > stop.mMaxError)
            break;
        if (std::get<4>(top) != version[static_cast<std::size_t>(a)] ||
            std::get<5>(top) != version[static_cast<std::size_t>(b)])
            continue;  // stale
        if (pinned[static_cast<std::size_t>(a)] && pinned[static_cast<std::size_t>(b)])
            continue;

        // Guard 1: T(ab), the shared tets.
        const std::vector<std::int64_t> shared = dv_sorted_intersection(
            vtets[static_cast<std::size_t>(a)], vtets[static_cast<std::size_t>(b)]);
        if (shared.empty() && vtets[static_cast<std::size_t>(a)].empty() &&
            vtets[static_cast<std::size_t>(b)].empty())
            continue;  // defensive: neither endpoint has any tet left

        // Guard 2: the vertex-link condition, exact set equality.
        {
            const std::vector<std::int64_t> link_a =
                detail::decim_corner_ring<4>(vtets[static_cast<std::size_t>(a)], tet_conn, a);
            const std::vector<std::int64_t> link_b =
                detail::decim_corner_ring<4>(vtets[static_cast<std::size_t>(b)], tet_conn, b);
            std::vector<std::int64_t> link_ab;
            link_ab.reserve(shared.size() * 2);
            for (std::int64_t t : shared) {
                const std::int64_t* c = tet_conn.data() + static_cast<std::size_t>(t) * 4;
                for (int k = 0; k < 4; ++k)
                    if (c[k] != a && c[k] != b)
                        link_ab.push_back(c[k]);
            }
            std::sort(link_ab.begin(), link_ab.end());
            link_ab.erase(std::unique(link_ab.begin(), link_ab.end()), link_ab.end());
            const std::vector<std::int64_t> common = dv_sorted_intersection(link_a, link_b);
            if (common != link_ab) {
                ++result.mCollapsesRejected;
                continue;
            }
        }

        // Guard 3: the duplicate-tet guard -- no survivor tet of `a` may share
        // its opposite triple with a survivor tet of `b`.
        bool duplicate = false;
        {
            std::vector<std::array<std::int64_t, 3>> triples_a;
            for (std::int64_t t : vtets[static_cast<std::size_t>(a)])
                if (!std::binary_search(shared.begin(), shared.end(), t))
                    triples_a.push_back(dv_opposite_triple(tet_conn, t, a));
            std::vector<std::array<std::int64_t, 3>> triples_b;
            for (std::int64_t t : vtets[static_cast<std::size_t>(b)])
                if (!std::binary_search(shared.begin(), shared.end(), t))
                    triples_b.push_back(dv_opposite_triple(tet_conn, t, b));
            for (const auto& ta : triples_a)
                for (const auto& tb : triples_b)
                    if (ta == tb) {
                        duplicate = true;
                        break;
                    }
        }
        if (duplicate) {
            ++result.mCollapsesRejected;
            continue;
        }

        // Guard 4: the mesh's own skin's link condition -- `decimate`'s own
        // ring/shared-face check, applied over `vboundary`/`boundary` instead
        // of a triangulated 2D input. `bshared` is the set of alive boundary
        // faces containing BOTH endpoints (empty, and so trivially satisfied,
        // when neither endpoint's boundary-face row has any overlap -- which
        // is always true for a purely interior collapse, since an empty row
        // intersects nothing).
        const std::vector<std::int64_t> bshared = dv_sorted_intersection(
            vboundary[static_cast<std::size_t>(a)], vboundary[static_cast<std::size_t>(b)]);
        if (bshared.size() > 2) {
            ++result.mCollapsesRejected;
            continue;
        }
        {
            const std::vector<std::int64_t> bring_a =
                decim_vertex_ring(vboundary[static_cast<std::size_t>(a)], boundary.mCorners, a);
            const std::vector<std::int64_t> bring_b =
                decim_vertex_ring(vboundary[static_cast<std::size_t>(b)], boundary.mCorners, b);
            if (decim_count_common(bring_a, bring_b) != bshared.size()) {
                ++result.mCollapsesRejected;
                continue;
            }
        }

        // Placement (recomputed -- versions matched, so bit-identical to what
        // was scored at push time).
        const DecimPlaced placed = decim_place(ctx, a, b);

        // Guard 5 (boundary-touching edges only): normal-flip over every
        // surviving alive boundary face incident to `a` or `b`.
        bool bflips = false;
        if (!vboundary[static_cast<std::size_t>(a)].empty() ||
            !vboundary[static_cast<std::size_t>(b)].empty()) {
            std::vector<std::int64_t> bmerged;
            bflips = detail::decim_any_normal_flip(
                xyz, boundary.mCorners, vboundary[static_cast<std::size_t>(a)],
                vboundary[static_cast<std::size_t>(b)], bshared, a, b, placed.mX, bmerged);
        }
        if (bflips) {
            ++result.mCollapsesRejected;
            continue;
        }

        // Guard 6: tet-inversion, over every ALIVE tet incident to `a` or `b`
        // that dies in neither this collapse (i.e. not in `shared`).
        bool flips = false;
        for (std::int64_t t : vtets[static_cast<std::size_t>(a)]) {
            if (std::binary_search(shared.begin(), shared.end(), t))
                continue;
            const double before = dv_tet_volume(xyz, tet_conn, t, -1, -1, nullptr);
            const double after = dv_tet_volume(xyz, tet_conn, t, a, b, placed.mX);
            if (dv_tet_inverts(before, after)) {
                flips = true;
                break;
            }
        }
        if (!flips)
            for (std::int64_t t : vtets[static_cast<std::size_t>(b)]) {
                if (std::binary_search(shared.begin(), shared.end(), t))
                    continue;
                const double before = dv_tet_volume(xyz, tet_conn, t, -1, -1, nullptr);
                const double after = dv_tet_volume(xyz, tet_conn, t, a, b, placed.mX);
                if (dv_tet_inverts(before, after)) {
                    flips = true;
                    break;
                }
            }
        if (flips) {
            ++result.mCollapsesRejected;
            continue;
        }

        // --- commit ------------------------------------------------------------
        const std::int64_t s = pinned[static_cast<std::size_t>(b)] ? b : a;
        const std::int64_t r = s == a ? b : a;

        detail::decim_commit_blend(xyz, quads, float_data, placed.mX, s, r);

        alive_tets -= static_cast<std::int64_t>(shared.size());
        result.mTetsRemoved += static_cast<std::int64_t>(shared.size());
        detail::decim_kill_and_migrate<4>(tet_conn, vtets, tet_alive, shared, s, r);

        // Kill the shared boundary faces, then migrate r's remaining alive
        // boundary faces onto s (decimate.cpp's own vfaces/face_alive commit
        // pattern, applied to the mesh's own skin).
        detail::decim_kill_and_migrate<3>(boundary.mCorners, vboundary, boundary_alive, bshared, s,
                                          r);

        successor[static_cast<std::size_t>(r)] = s;
        ++version[static_cast<std::size_t>(s)];
        ++version[static_cast<std::size_t>(r)];
        touches_boundary[static_cast<std::size_t>(s)] =
            vboundary[static_cast<std::size_t>(s)].empty() ? 0 : 1;
        if (regime == 0)
            result.mMaxErrorApplied = std::max(result.mMaxErrorApplied, score);

        // Re-score the survivor's remaining edges via its updated link.
        const std::vector<std::int64_t> link_s =
            detail::decim_corner_ring<4>(vtets[static_cast<std::size_t>(s)], tet_conn, s);
        for (std::int64_t w : link_s) {
            if (pinned[static_cast<std::size_t>(s)] && pinned[static_cast<std::size_t>(w)])
                continue;
            const std::int64_t lo = s < w ? s : w;
            const std::int64_t hi = s < w ? w : s;
            const auto rs = dv_score(lo, hi);
            if (!std::isfinite(rs.second))
                continue;
            heap.push({rs.first, rs.second, lo, hi, version[static_cast<std::size_t>(lo)],
                       version[static_cast<std::size_t>(hi)]});
        }
    }

    // --- phase 9: emission ------------------------------------------------------
    std::vector<char> used;
    std::vector<std::int64_t> remap;
    const std::size_t num_used =
        detail::decim_compact_points<4>(tet_alive, tet_conn, n, used, remap);

    Mesh& out = result.mMesh;

    out.AssignPoints(detail::decim_emit_points(rMesh.Points(), xyz, used, n, num_used, dim));

    {
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const std::size_t nc = cb.NumCells();
            const std::size_t base = static_cast<std::size_t>(bases[bi]);
            std::size_t alive_in_block = 0;
            for (std::size_t c = 0; c < nc; ++c)
                if (tet_alive[base + c])
                    ++alive_in_block;
            NDArray block = NDArray::Uninit(DType::Int64, {alive_in_block, 4});
            std::int64_t* dst = block.As<std::int64_t>();
            std::vector<std::int64_t> cell_map(nc, -1);
            std::size_t w = 0;
            for (std::size_t c = 0; c < nc; ++c) {
                if (!tet_alive[base + c])
                    continue;
                for (std::size_t k = 0; k < 4; ++k)
                    dst[w * 4 + k] = remap[static_cast<std::size_t>(tet_conn[(base + c) * 4 + k])];
                cell_map[c] = static_cast<std::int64_t>(w);
                ++w;
            }
            out.AddCellBlock(cell_type_name(CellType::Tetra), std::move(block));
            result.mCellMaps.push_back(detail::decim_int64_vector(cell_map));
            ++bi;
        }
    }

    detail::decim_emit_data(rMesh, tet_alive, bases, nullptr, float_data, used, n, num_used, out);

    result.mPointMap = detail::decim_point_map(successor, remap, n);
    result.mPointsRemoved = static_cast<std::int64_t>(n - num_used);

    {
        detail::RegionRemap rmap;
        rmap.pPointMap = &result.mPointMap;
        rmap.mCellMapKind = detail::CellMapKind::Direct;
        rmap.pCellMaps = &result.mCellMaps;
        rmap.mOpName = "decimate_volume";
        detail::remap_regions(rMesh, result.mMesh, rmap);
    }

    return result;
}

}  // namespace meshioplusplus
