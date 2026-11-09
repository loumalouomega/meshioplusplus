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
// Quadric-error-metric (Garland-Heckbert) edge-collapse surface decimation.
// Built entirely through the uniform mesh API so it compiles under every mesh
// backend. See operations/decimate.hpp for the contract.
//
// Determinism: the setup passes are parallel with a fixed FP order (each vertex
// sums its quadric in ascending incident-face index), and the greedy collapse
// loop is SERIAL, driven by a priority queue with the total order (error, lower
// endpoint id, higher endpoint id, versions) and lazy deletion via per-vertex
// version counters. Every floating-point expression that feeds a branch (the
// plane, the quadric, the 3x3 solve, the error form, the blend parameter) is
// transcribed token for token into _decimate.py, and none of them may be
// contracted into FMAs across the C++/numpy boundary -- the same build-flag
// assumption slice/interpolate's parity already relies on. Output is
// byte-identical across mesh backends, thread counts, and the numpy fallback.

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
#include <vector>

// Project includes
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/detail/decimate_common.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/hausdorff_samples.hpp"
#include "../detail/decimate_shared.hpp"
#include "../detail/slot_runs.hpp"
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {

namespace {

// The placement solve, error form, normal-flip check, and the vertex-ring
// link-condition primitives are shared with volume decimation and live in
// detail/decimate_common.{hpp,cpp} -- see that header for why this is a hoist
// rather than a duplication.
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

// --- validation --------------------------------------------------------------

constexpr const char* kDecimPrefix = "meshio++: decimate: ";

// Reject, by name, every construct outside the surface scope -- BEFORE the
// simplexify call, whose own policy is to pass unsupported blocks through with
// only a warning.
void decim_check_blocks(const Mesh& rMesh) {
    bool has_surface = false;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                "meshio++: decimate: mesh contains a polyhedron cell block; decimate operates on "
                "surface meshes -- run extract_surface first (or convert_cells(simplexify) "
                "if you wanted the volume decomposed instead)");
        const CellType ct = cell_type_from_name(std::string(cb.Type()));
        const int dim = cell_type_dimension(ct);
        if (dim == 3)
            throw std::invalid_argument(
                "meshio++: decimate: mesh contains 3D volume cell block '" +
                std::string(cb.Type()) +
                "'; decimate operates on surface meshes (extract_surface first)");
        if (dim == 2)
            has_surface = true;
    }
    if (!has_surface)
        throw std::invalid_argument("meshio++: decimate: mesh contains no surface (2D) cell block");
    for (const auto cb : rMesh.CellRange()) {
        const CellType ct = cell_type_from_name(std::string(cb.Type()));
        const int dim = cell_type_dimension(ct);
        if (dim <= 1)
            throw std::invalid_argument(
                "meshio++: decimate: cannot decimate a mesh containing cell block '" +
                std::string(cb.Type()) +
                "' alongside its surface blocks (its nodes would dangle after the collapse; "
                "drop it first, e.g. via split)");
        if (cb.IsRagged())
            throw std::invalid_argument("meshio++: decimate: cannot decimate ragged cell block '" +
                                        std::string(cb.Type()) + "' (triangulate/clean it first)");
        if (ct != CellType::Triangle && ct != CellType::Quad && ct != CellType::Polygon)
            throw std::invalid_argument(
                "meshio++: decimate: cannot decimate higher-order cell block '" +
                std::string(cb.Type()) + "' (linearize the mesh first)");
    }
}

// --- flat face table ---------------------------------------------------------
// DecimFaces/DecimCsr and the generic per-face/per-vertex quadric machinery
// live in detail/decimate_common.hpp (brought into scope above); this file
// keeps only what ties them to a *simplexified `Mesh`*'s own blocks.

DecimFaces decim_build_faces(const Mesh& rSimp, std::size_t n) {
    DecimFaces t;
    for (const auto cb : rSimp.CellRange()) {
        t.mBlockBase.push_back(static_cast<std::int64_t>(t.mNumFaces));
        if (cb.NumCells() == 0)
            continue;
        if (cell_type_from_name(std::string(cb.Type())) != CellType::Triangle)
            throw std::invalid_argument("meshio++: decimate: internal error: simplexified block '" +
                                        std::string(cb.Type()) + "' is not triangular");
        const NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t nc = cb.NumCells();
        const std::size_t first = t.mCorners.size();
        t.mCorners.resize(first + nc * 3);
        for (std::size_t c = 0; c < nc; ++c) {
            for (std::size_t k = 0; k < 3; ++k) {
                const std::int64_t id = detail::read_int(conn, c * npc + k);
                if (id < 0 || static_cast<std::size_t>(id) >= n)
                    throw std::invalid_argument(
                        "meshio++: decimate: connectivity references point " + std::to_string(id) +
                        " outside [0, " + std::to_string(n) + ")");
                t.mCorners[first + c * 3 + k] = id;
            }
        }
        t.mNumFaces += nc;
    }
    return t;
}

// --- edge table + boundary marks ---------------------------------------------

using DecimEdgeKey = std::array<std::int64_t, 2>;

// Unique edges in first-seen (face-major) order, and the boundary vertex marks
// from the once-used-edge test. Phase 1 fills the per-face-edge keys into
// disjoint slots in parallel; phase 2 is the SERIAL first-seen dedup that pins
// determinism (surface.cpp's phase-split idiom).
std::vector<DecimEdgeKey> decim_build_edges(const DecimFaces& rFaces,
                                            std::vector<std::uint8_t>& rBoundary) {
    const std::size_t nf = rFaces.mNumFaces;
    std::vector<DecimEdgeKey> keys(nf * 3);
    parallel_for(nf, [&](std::size_t f) {
        const std::int64_t* c = rFaces.mCorners.data() + f * 3;
        for (std::size_t e = 0; e < 3; ++e) {
            const std::int64_t a = c[e];
            const std::int64_t b = c[(e + 1) % 3];
            keys[f * 3 + e] = a < b ? DecimEdgeKey{a, b} : DecimEdgeKey{b, a};
        }
    });

    // A degenerate face's self-edge is not an edge: drop those slots, keeping
    // the rest in slot order, then number edges in first-seen order and count
    // their uses (detail/slot_runs.hpp) -- the ids and counts the former
    // serial hash sweep gave.
    keys.erase(std::remove_if(keys.begin(), keys.end(),
                              [](const DecimEdgeKey& k) { return k[0] == k[1]; }),
               keys.end());
    const std::size_t nslots = keys.size();
    const std::size_t npts = rBoundary.size();  // one flag per point
    const detail::SlotRuns runs =
        detail::group_slots(keys, npts + 1, [npts](const DecimEdgeKey& rK) {
            return rK[0] >= 0 && static_cast<std::size_t>(rK[0]) < npts
                       ? static_cast<std::size_t>(rK[0])
                       : npts;
        });
    const detail::FirstSeen seen = detail::number_first_seen(runs, nslots);
    std::vector<DecimEdgeKey> edges(seen.NumIds());
    std::vector<std::uint32_t> use_count(seen.NumIds());
    parallel_for(seen.NumIds(), [&](std::size_t id) {
        const std::size_t r = static_cast<std::size_t>(seen.mRunOfId[id]);
        edges[id] = keys[runs.Head(r)];
        use_count[id] = static_cast<std::uint32_t>(runs.Size(r));
    });
    for (std::size_t e = 0; e < edges.size(); ++e)
        if (use_count[e] == 1) {
            rBoundary[static_cast<std::size_t>(edges[e][0])] = 1;
            rBoundary[static_cast<std::size_t>(edges[e][1])] = 1;
        }
    return edges;
}

// --- placement, error, normal-flip and link-condition primitives -----------
// Hoisted into detail/decimate_common.{hpp,cpp}; brought into scope above.

// A Float64 working copy of a float-kind point_data array, blended in place at
// each collapse commit, and the rest of the commit/emission plumbing shared
// with volume decimation: detail/decimate_shared.hpp.
using detail::DecimFloatData;

}  // namespace

DecimatePlacement decimate_placement_from_name(const std::string& rName) {
    if (rName == "optimal")
        return DecimatePlacement::Optimal;
    if (rName == "midpoint")
        return DecimatePlacement::Midpoint;
    if (rName == "endpoint")
        return DecimatePlacement::Endpoint;
    throw std::invalid_argument("meshio++: decimate: unknown placement '" + rName +
                                "' (expected 'optimal', 'midpoint' or 'endpoint')");
}

namespace {

DecimateResult decim_core(const Mesh& rMesh, const DecimateOptions& rOptions) {
    detail::DecimStop stop =
        detail::decim_resolve_stop(rOptions.mTargetRatio, rOptions.mTargetFaces, rOptions.mMaxError,
                                   kDecimPrefix, "target_faces");

    const std::size_t n = rMesh.NumPoints();
    detail::decim_check_frozen(rOptions.mFrozen.size(), n, kDecimPrefix);

    decim_check_blocks(rMesh);

    // --- phase 0: triangulate (quads/polygons fan into triangles; a pure
    // triangle mesh passes through 1:1). Higher-order input was rejected above,
    // so the points are untouched and the simplexify point map is the identity.
    ConvertCellsResult prep =
        convert_cells(rMesh, {ConvertCellsMode::Simplexify, /*mRecordParentIds=*/true});
    const Mesh& simp = prep.mMesh;
    const std::size_t nblocks = simp.NumCellBlocks();
    const std::size_t dim = simp.PointDim();

    // --- phase 1: flat buffers + incidence (parallel fills, serial dedups) ---
    std::vector<double> xyz = detail::surface_xyz(simp);

    DecimFaces faces = decim_build_faces(simp, n);
    const std::size_t num_faces = faces.mNumFaces;
    const DecimCsr csr = decim_vertex_faces_csr(faces, n);

    std::vector<std::uint8_t> boundary(n, 0);
    const std::vector<DecimEdgeKey> edges = decim_build_edges(faces, boundary);

    std::vector<double> quad_k;
    std::vector<double> fnormals;
    decim_face_planes(faces, xyz, quad_k, fnormals);

    // --- phase 2: the pin mask (caller | boundary | feature) -----------------
    const std::vector<std::uint8_t> pinned = detail::decim_pin_mask(
        rOptions.mFrozen, boundary, rOptions.mPreserveBoundary, rOptions.mPreserveFeatures,
        rOptions.mFeatureAngleDeg, faces.mCorners, fnormals, n);

    // --- phase 3: per-vertex quadrics (fixed ascending-face FP order) --------
    std::vector<double> quads = decim_accumulate_quadrics(csr, n, quad_k);
    quad_k.clear();
    quad_k.shrink_to_fit();

    // --- phase 4: mutable collapse state -------------------------------------
    std::vector<std::uint8_t> face_alive(num_faces, 1);
    std::vector<std::int64_t> successor(n, -1);  // removed vertex -> survivor
    std::vector<std::uint32_t> version(n, 0);
    std::vector<std::vector<std::int64_t>> vfaces(n);
    for (std::size_t v = 0; v < n; ++v) {
        std::vector<std::int64_t>& row = vfaces[v];
        row.reserve(static_cast<std::size_t>(csr.mXadj[v + 1] - csr.mXadj[v]));
        for (std::int64_t k = csr.mXadj[v]; k < csr.mXadj[v + 1]; ++k) {
            const std::int64_t f = csr.mAdj[static_cast<std::size_t>(k)];
            if (!row.empty() && row.back() == f)
                continue;  // a degenerate face repeating this corner
            row.push_back(f);
        }
    }

    // Float64 working copies of the float-kind point_data, blended per commit.
    std::vector<DecimFloatData> float_data = detail::decim_float_working_copies(simp, n);

    DecimPlaceCtx ctx;
    ctx.mpXyz = &xyz;
    ctx.mpQ = &quads;
    ctx.mpPinned = &pinned;
    ctx.mPlacement = rOptions.mPlacement;

    // Resolve a ratio into a face target now that F is known.
    if (stop.mUseTarget && rOptions.mTargetRatio >= 0.0)
        stop.mTarget = detail::decim_ratio_target(rOptions.mTargetRatio, num_faces);

    // --- phase 5: seed the queue (parallel scoring, serial pushes) -----------
    using DecimHeapEntry =
        std::tuple<double, std::int64_t, std::int64_t, std::uint32_t, std::uint32_t>;
    std::priority_queue<DecimHeapEntry, std::vector<DecimHeapEntry>, std::greater<DecimHeapEntry>>
        heap;
    {
        std::vector<double> seed_err(edges.size(), 0.0);
        parallel_for(edges.size(), [&](std::size_t e) {
            if (pinned[static_cast<std::size_t>(edges[e][0])] &&
                pinned[static_cast<std::size_t>(edges[e][1])])
                return;
            seed_err[e] = decim_place(ctx, edges[e][0], edges[e][1]).mErr;
        });
        for (std::size_t e = 0; e < edges.size(); ++e) {
            const std::int64_t a = edges[e][0];
            const std::int64_t b = edges[e][1];
            if (pinned[static_cast<std::size_t>(a)] && pinned[static_cast<std::size_t>(b)])
                continue;  // a pinned-pinned collapse is impossible by definition
            if (!std::isfinite(seed_err[e]))
                continue;
            heap.push({seed_err[e], a, b, 0, 0});
        }
    }

    // --- phase 6: the SERIAL greedy loop -------------------------------------
    // QEM is inherently sequential; serializing it (with the total heap order
    // and lazy version-stamped deletion) is the determinism pin. The guards run
    // in a fixed order and each rejection is counted.
    DecimateResult result;
    std::int64_t alive_faces = static_cast<std::int64_t>(num_faces);
    std::vector<std::int64_t> merged;
    while (true) {
        if (stop.mUseTarget && alive_faces <= stop.mTarget)
            break;
        if (heap.empty()) {
            if (stop.mUseTarget && alive_faces > stop.mTarget)
                log::warn(
                    "decimate: queue exhausted at {} faces before reaching the target of {} "
                    "(pinning and the validity guards left no collapsible edge)",
                    alive_faces, stop.mTarget);
            break;
        }
        const DecimHeapEntry top = heap.top();
        heap.pop();
        const double err = std::get<0>(top);
        const std::int64_t a = std::get<1>(top);
        const std::int64_t b = std::get<2>(top);
        // The heap minimum bounds every remaining candidate (each live edge has
        // a fresh entry), so this stop is valid even on a stale entry.
        if (stop.mUseMaxError && err > stop.mMaxError)
            break;
        if (std::get<3>(top) != version[static_cast<std::size_t>(a)] ||
            std::get<4>(top) != version[static_cast<std::size_t>(b)])
            continue;  // stale: some commit touched an endpoint since the push
        if (pinned[static_cast<std::size_t>(a)] && pinned[static_cast<std::size_t>(b)])
            continue;  // defensive: never pushed in the first place

        // Guard 1: the link condition. `shared` are the alive faces on the
        // edge: 3+ means a non-manifold edge, and two boundary endpoints joined
        // by an interior edge would pinch the surface. The classic condition
        // |ring(a) & ring(b)| == |shared| forbids every other topology change.
        std::vector<std::int64_t> shared;
        {
            const std::vector<std::int64_t>& fa = vfaces[static_cast<std::size_t>(a)];
            const std::vector<std::int64_t>& fb = vfaces[static_cast<std::size_t>(b)];
            std::size_t i = 0;
            std::size_t j = 0;
            while (i < fa.size() && j < fb.size()) {
                if (fa[i] < fb[j])
                    ++i;
                else if (fb[j] < fa[i])
                    ++j;
                else {
                    shared.push_back(fa[i]);
                    ++i;
                    ++j;
                }
            }
        }
        if (shared.empty())
            continue;  // defensive: the edge no longer exists
        if (shared.size() > 2) {
            ++result.mCollapsesRejected;
            continue;
        }
        if (boundary[static_cast<std::size_t>(a)] && boundary[static_cast<std::size_t>(b)] &&
            shared.size() == 2) {
            // Only reachable with mPreserveBoundary off: joining two boundary
            // vertices across the interior would pinch the surface.
            ++result.mCollapsesRejected;
            continue;
        }
        const std::vector<std::int64_t> ring_a =
            decim_vertex_ring(vfaces[static_cast<std::size_t>(a)], faces.mCorners, a);
        const std::vector<std::int64_t> ring_b =
            decim_vertex_ring(vfaces[static_cast<std::size_t>(b)], faces.mCorners, b);
        if (decim_count_common(ring_a, ring_b) != shared.size()) {
            ++result.mCollapsesRejected;
            continue;
        }

        // Guard 2: placement (recomputed -- the versions matched, so this is
        // bit-identical to what was scored at push time) and the normal-flip /
        // degenerate-survivor sweep, in ascending face index. An already
        // degenerate face (zero normal) imposes no constraint -- smooth's "do
        // no harm" rule -- and a survivor collapsing to zero area yields
        // n_after = 0, hence dot = 0, and is rejected by the same test.
        const DecimPlaced placed = decim_place(ctx, a, b);
        const bool flips = detail::decim_any_normal_flip(
            xyz, faces.mCorners, vfaces[static_cast<std::size_t>(a)],
            vfaces[static_cast<std::size_t>(b)], shared, a, b, placed.mX, merged);
        if (flips) {
            ++result.mCollapsesRejected;
            continue;
        }

        // --- commit: v merges into the survivor s ---------------------------
        const std::int64_t s = pinned[static_cast<std::size_t>(b)] ? b : a;
        const std::int64_t r = s == a ? b : a;

        // Blend the data, move the survivor, merge the quadrics.
        detail::decim_commit_blend(xyz, quads, float_data, placed.mX, s, r);

        // Kill the shared faces (removing them from every corner's list), then
        // migrate the removed vertex's remaining faces onto the survivor.
        alive_faces -= static_cast<std::int64_t>(shared.size());
        result.mFacesRemoved += static_cast<std::int64_t>(shared.size());
        detail::decim_kill_and_migrate<3>(faces.mCorners, vfaces, face_alive, shared, s, r);
        successor[static_cast<std::size_t>(r)] = s;
        ++version[static_cast<std::size_t>(s)];
        ++version[static_cast<std::size_t>(r)];
        result.mMaxErrorApplied = std::max(result.mMaxErrorApplied, err);

        // Re-score the survivor's edges with the bumped versions. Rejected
        // edges are deliberately not re-pushed here or at rejection time: an
        // edge re-enters the queue only when a commit touches an endpoint.
        const std::vector<std::int64_t> ring_s =
            decim_vertex_ring(vfaces[static_cast<std::size_t>(s)], faces.mCorners, s);
        for (const std::int64_t w : ring_s) {
            if (pinned[static_cast<std::size_t>(s)] && pinned[static_cast<std::size_t>(w)])
                continue;
            const std::int64_t lo = s < w ? s : w;
            const std::int64_t hi = s < w ? w : s;
            const double e2 = decim_place(ctx, lo, hi).mErr;
            if (!std::isfinite(e2))
                continue;
            heap.push({e2, lo, hi, version[static_cast<std::size_t>(lo)],
                       version[static_cast<std::size_t>(hi)]});
        }
    }

    // --- phase 7: emission ----------------------------------------------------
    // Alive faces in original (block, face) order; points compacted ascending
    // (surface.cpp's used/remap pattern).
    std::vector<char> used;
    std::vector<std::int64_t> remap;
    const std::size_t num_used =
        detail::decim_compact_points<3>(face_alive, faces.mCorners, n, used, remap);

    Mesh& out = result.mMesh;

    // Points: the working coordinates of the used vertices, cast back once to
    // the input dtype.
    out.AssignPoints(detail::decim_emit_points(simp.Points(), xyz, used, n, num_used, dim));

    // One triangle block per input block (1:1, possibly empty), plus the
    // per-input-block cell maps composed through the simplexify child map.
    detail::decim_emit_blocks<3>(rMesh, simp, "convert:parent_cell", face_alive, faces.mCorners,
                                 faces.mBlockBase, remap, CellType::Triangle, out,
                                 result.mCellMaps);

    // cell_data: the surviving faces' own rows (simplexify already replicated a
    // parent's row to its triangles); convert:parent_cell is consumed above,
    // never forwarded. point_data: float-kind arrays from the blended working
    // copies (cast back to the source dtype), everything else as the
    // survivor's own byte rows.
    detail::decim_emit_data(simp, face_alive, faces.mBlockBase, "convert:parent_cell", float_data,
                            used, n, num_used, out);

    // mPointMap: every input point lands on its survivor's output index.
    result.mPointMap = detail::decim_point_map(successor, remap, n);
    result.mPointsRemoved = static_cast<std::int64_t>(n - num_used);

    // Named regions. The cell map is FirstChild — an input quad/polygon becomes
    // a contiguous run of triangles, and a fully-collapsed parent maps to -1,
    // which the run scan steps over. A side edge moves to the triangle edge
    // joining its surviving nodes; one a collapse removed is lost, with a
    // warning. See detail/region_remap.hpp.
    detail::decim_remap_regions(rMesh, result.mMesh, result.mPointMap, result.mCellMaps,
                                detail::CellMapKind::FirstChild, "decimate");

    return result;
}

std::int64_t decim_num_faces(const Mesh& rMesh) {
    std::int64_t n = 0;
    for (const auto cb : rMesh.CellRange())
        n += static_cast<std::int64_t>(cb.NumCells());
    return n;
}

// The input's triangulated face count: a quad is two, a polygon n - 2.
std::int64_t decim_input_faces(const Mesh& rMesh) {
    std::int64_t n = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (!cb.IsRagged()) {
            n += static_cast<std::int64_t>(cb.NumCells()) *
                 std::max<std::int64_t>(static_cast<std::int64_t>(cb.NodesPerCell()) - 2, 0);
            continue;
        }
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            n += std::max<std::int64_t>(static_cast<std::int64_t>(cb.RowSize(c)) - 2, 0);
    }
    return n;
}

// max_deviation and max_faces (roadmap 7.8.5): see decimate.hpp.
DecimateResult decim_bounded(const Mesh& rMesh, const DecimateOptions& rOptions) {
    const DecimateOptions& o = rOptions;
    const bool deviation = o.mMaxDeviation >= 0.0;
    const int others = (o.mTargetRatio >= 0.0 ? 1 : 0) + (o.mTargetFaces >= 0 ? 1 : 0) +
                       (o.mMaxError >= 0.0 ? 1 : 0);
    if (std::isnan(o.mMaxDeviation) || (deviation && !std::isfinite(o.mMaxDeviation)))
        throw std::invalid_argument(std::string(kDecimPrefix) + "max_deviation must be finite");
    if (deviation && others > 0)
        throw std::invalid_argument(std::string(kDecimPrefix) +
                                    "max_deviation replaces target_ratio, target_faces and "
                                    "max_error; set only one");
    if (!(o.mCoverage > 0.0 && o.mCoverage <= 1.0))
        throw std::invalid_argument(std::string(kDecimPrefix) + "coverage must be in (0, 1]");
    const auto at_faces = [&](std::int64_t Faces) {
        DecimateOptions t = o;
        t.mTargetRatio = -1.0;
        t.mMaxError = -1.0;
        t.mMaxDeviation = -1.0;
        t.mMaxFaces = -1;
        t.mTargetFaces = Faces;
        return decim_core(rMesh, t);
    };
    const std::int64_t cap = o.mMaxFaces;
    if (!deviation) {
        if (others > 1)
            throw std::invalid_argument(
                std::string(kDecimPrefix) +
                "exactly one of target_ratio, target_faces, max_error must be set");
        if (others == 0 || o.mTargetFaces >= 0) {
            DecimateResult r = at_faces(o.mTargetFaces >= 0 ? std::min(o.mTargetFaces, cap) : cap);
            r.mRounds = 1;
            return r;
        }
        DecimateOptions t = o;
        t.mMaxFaces = -1;
        DecimateResult r = decim_core(rMesh, t);
        r.mRounds = 1;
        if (decim_num_faces(r.mMesh) > cap) {
            r = at_faces(cap);
            r.mRounds = 2;
        }
        return r;
    }
    // The deviation search.
    const double limit = o.mMaxDeviation;
    std::int64_t rounds = 0;
    // The input at its vertices, the coarser result with 2 x 2 face samples.
    const detail::HausdorffSampler sampler(rMesh, 0, 2);
    const auto measure = [&](const DecimateResult& rR) {
        ++rounds;
        return detail::coverage_quantile(sampler.Distances(rR.mMesh), o.mCoverage);
    };
    const std::int64_t n0 = decim_input_faces(rMesh);
    const std::int64_t hi0 = cap >= 0 ? std::min(cap, n0) : n0;
    DecimateResult best = at_faces(hi0);
    double best_dev = measure(best);
    if (best_dev > limit) {
        best.mDeviation = best_dev;
        best.mDeviationMet = false;
        best.mRounds = rounds;
        return best;
    }
    std::int64_t hi = std::min(hi0, decim_num_faces(best.mMesh));
    std::int64_t lo = 0;
    {
        // The quadric bound: a collapse whose error is at most limit^2 keeps
        // its vertex within limit of every plane it absorbed.
        DecimateOptions t = o;
        t.mMaxDeviation = -1.0;
        t.mMaxFaces = -1;
        t.mMaxError = limit * limit;
        DecimateResult q = decim_core(rMesh, t);
        const std::int64_t nq = decim_num_faces(q.mMesh);
        if (nq < hi) {
            const double d = measure(q);
            if (d <= limit) {
                best = std::move(q);
                best_dev = d;
                hi = nq;
            } else {
                lo = nq;
            }
        }
    }
    while (hi - lo > std::max<std::int64_t>(1, hi / 100) && rounds < 40) {
        const std::int64_t mid = lo + (hi - lo) / 2;
        DecimateResult r = at_faces(mid);
        const std::int64_t got = decim_num_faces(r.mMesh);
        const double d = measure(r);
        if (d <= limit && got < hi) {
            best = std::move(r);
            best_dev = d;
            hi = got;
        } else if (d <= limit) {
            hi = mid;  // the run stopped short of mid: nothing lighter here
        } else {
            lo = mid;
        }
    }
    best.mDeviation = best_dev;
    best.mDeviationMet = true;
    best.mRounds = rounds;
    return best;
}

}  // namespace

DecimateResult decimate(const Mesh& rMesh, const DecimateOptions& rOptions) {
    if (rOptions.mMaxDeviation < 0.0 && rOptions.mMaxFaces < 0)
        return decim_core(rMesh, rOptions);
    return decim_bounded(rMesh, rOptions);
}

}  // namespace meshioplusplus
