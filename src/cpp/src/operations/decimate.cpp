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
#include <cstring>
#include <functional>
#include <iterator>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/detail/decimate_common.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// The placement solve, error form, normal-flip check, and the vertex-ring
// link-condition primitives are shared with volume decimation and live in
// detail/decimate_common.{hpp,cpp} -- see that header for why this is a hoist
// rather than a duplication.
using detail::decim_accumulate_quadrics;
using detail::decim_count_common;
using detail::decim_face_normal;
using detail::decim_face_planes;
using detail::decim_mark_features;
using detail::decim_place;
using detail::decim_quadric_error;
using detail::decim_sorted_erase;
using detail::decim_sorted_insert;
using detail::decim_vertex_faces_csr;
using detail::decim_vertex_ring;
using detail::DecimCsr;
using detail::DecimFaces;
using detail::DecimPlaceCtx;
using detail::DecimPlaced;

// --- validation --------------------------------------------------------------

// The resolved stopping criterion. `mTargetFaces` is only meaningful when
// `mUseTarget` is set (a ratio is turned into a face count once the
// triangulated face count is known).
struct DecimStop {
    bool mUseTarget = false;
    bool mUseMaxError = false;
    std::int64_t mTargetFaces = -1;
    double mMaxError = -1.0;
};

DecimStop decim_resolve_stop(const DecimateOptions& rOptions) {
    const int num_set = (rOptions.mTargetRatio >= 0.0 ? 1 : 0) +
                        (rOptions.mTargetFaces >= 0 ? 1 : 0) + (rOptions.mMaxError >= 0.0 ? 1 : 0);
    if (num_set != 1)
        throw std::invalid_argument(
            "meshio++: decimate: exactly one of target_ratio, target_faces, max_error must be "
            "set; got " +
            std::to_string(num_set));
    DecimStop stop;
    if (rOptions.mTargetRatio >= 0.0) {
        if (!(rOptions.mTargetRatio > 0.0 && rOptions.mTargetRatio <= 1.0))
            throw std::invalid_argument(
                "meshio++: decimate: target_ratio must lie in (0, 1]; got " +
                std::to_string(rOptions.mTargetRatio));
        stop.mUseTarget = true;  // the face count is filled in once F is known
    } else if (rOptions.mTargetFaces >= 0) {
        stop.mUseTarget = true;
        stop.mTargetFaces = rOptions.mTargetFaces;
    } else {
        stop.mUseMaxError = true;
        stop.mMaxError = rOptions.mMaxError;
    }
    return stop;
}

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

struct DecimEdgeKeyHash {
    std::size_t operator()(const DecimEdgeKey& rKey) const {
        std::size_t h = 0;
        for (std::int64_t v : rKey)
            h ^= std::hash<std::int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

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

    std::vector<DecimEdgeKey> edges;
    std::unordered_map<DecimEdgeKey, std::size_t, DecimEdgeKeyHash> edge_id;
    edge_id.reserve(nf * 3 * 2);
    std::vector<std::uint32_t> use_count;
    for (const DecimEdgeKey& key : keys) {
        if (key[0] == key[1])
            continue;  // a degenerate face's self-edge is not an edge
        auto it = edge_id.find(key);
        if (it == edge_id.end()) {
            edge_id.emplace(key, edges.size());
            edges.push_back(key);
            use_count.push_back(1);
        } else {
            ++use_count[it->second];
        }
    }
    for (std::size_t e = 0; e < edges.size(); ++e)
        if (use_count[e] == 1) {
            rBoundary[static_cast<std::size_t>(edges[e][0])] = 1;
            rBoundary[static_cast<std::size_t>(edges[e][1])] = 1;
        }
    return edges;
}

// --- placement, error, normal-flip and link-condition primitives -----------
// Hoisted into detail/decimate_common.{hpp,cpp}; brought into scope above.

// A vector of Int64 as a 1-D NDArray.
NDArray decim_int64_vector(const std::vector<std::int64_t>& rValues) {
    NDArray a = NDArray::Uninit(DType::Int64, {rValues.size()});
    if (!rValues.empty())
        std::memcpy(a.Data(), rValues.data(), rValues.size() * sizeof(std::int64_t));
    return a;
}

// A Float64 working copy of a float-kind point_data array, blended in place at
// each collapse commit.
struct DecimFloatData {
    std::string mName;
    std::size_t mNumComponents = 0;
    std::vector<double> mValues;  // n * mNumComponents
};

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

DecimateResult decimate(const Mesh& rMesh, const DecimateOptions& rOptions) {
    DecimStop stop = decim_resolve_stop(rOptions);

    const std::size_t n = rMesh.NumPoints();
    if (!rOptions.mFrozen.empty() && rOptions.mFrozen.size() != n)
        throw std::invalid_argument("meshio++: decimate: frozen mask has " +
                                    std::to_string(rOptions.mFrozen.size()) +
                                    " entries but the mesh has " + std::to_string(n) + " points");

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
    std::vector<double> xyz(n * 3, 0.0);
    {
        const NDArray& points = simp.Points();
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t d = 0; d < dim && d < 3; ++d)
                xyz[i * 3 + d] = detail::read_double(points, i * dim + d);
        });
    }

    DecimFaces faces = decim_build_faces(simp, n);
    const std::size_t num_faces = faces.mNumFaces;
    const DecimCsr csr = decim_vertex_faces_csr(faces, n);

    std::vector<std::uint8_t> boundary(n, 0);
    const std::vector<DecimEdgeKey> edges = decim_build_edges(faces, boundary);

    std::vector<double> quad_k;
    std::vector<double> fnormals;
    decim_face_planes(faces, xyz, quad_k, fnormals);

    // --- phase 2: the pin mask (caller | boundary | feature) -----------------
    std::vector<std::uint8_t> pinned(n, 0);
    if (!rOptions.mFrozen.empty())
        for (std::size_t i = 0; i < n; ++i)
            pinned[i] = rOptions.mFrozen[i] ? 1 : 0;
    if (rOptions.mPreserveBoundary)
        for (std::size_t i = 0; i < n; ++i)
            if (boundary[i])
                pinned[i] = 1;
    if (rOptions.mPreserveFeatures) {
        const double cos_thr = std::cos(rOptions.mFeatureAngleDeg * 3.14159265358979323846 / 180.0);
        decim_mark_features(csr, n, fnormals, cos_thr, pinned);
    }

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
    std::vector<DecimFloatData> float_data;
    for (const std::string& name : simp.PointDataNames()) {
        const NDArray& a = simp.PointData(name);
        if (!detail::is_float_dtype(a.Dtype()) || detail::rows(a) != n || n == 0)
            continue;
        DecimFloatData fd;
        fd.mName = name;
        fd.mNumComponents = a.Size() / n;
        fd.mValues.resize(a.Size());
        parallel_for_bw(a.Size(),
                        [&](std::size_t i) { fd.mValues[i] = detail::read_double(a, i); });
        float_data.push_back(std::move(fd));
    }

    DecimPlaceCtx ctx;
    ctx.mpXyz = &xyz;
    ctx.mpQ = &quads;
    ctx.mpPinned = &pinned;
    ctx.mPlacement = rOptions.mPlacement;

    // Resolve a ratio into a face target now that F is known. int64(x + 0.5)
    // truncation, never llround, whose tie rule the numpy twin could not
    // reproduce with int(x + 0.5).
    if (stop.mUseTarget && rOptions.mTargetRatio >= 0.0)
        stop.mTargetFaces = std::max<std::int64_t>(
            1, static_cast<std::int64_t>(rOptions.mTargetRatio * static_cast<double>(num_faces) +
                                         0.5));

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
        if (stop.mUseTarget && alive_faces <= stop.mTargetFaces)
            break;
        if (heap.empty()) {
            if (stop.mUseTarget && alive_faces > stop.mTargetFaces)
                log::warn(
                    "decimate: queue exhausted at {} faces before reaching the target of {} "
                    "(pinning and the validity guards left no collapsible edge)",
                    alive_faces, stop.mTargetFaces);
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
        bool flips = false;
        {
            const std::vector<std::int64_t>& fa = vfaces[static_cast<std::size_t>(a)];
            const std::vector<std::int64_t>& fb = vfaces[static_cast<std::size_t>(b)];
            merged.clear();
            std::merge(fa.begin(), fa.end(), fb.begin(), fb.end(), std::back_inserter(merged));
            std::size_t si = 0;
            for (const std::int64_t f : merged) {
                while (si < shared.size() && shared[si] < f)
                    ++si;
                if (si < shared.size() && shared[si] == f)
                    continue;  // dies in the collapse: no constraint
                const std::int64_t* c = faces.mCorners.data() + static_cast<std::size_t>(f) * 3;
                double n0[3];
                double n1[3];
                decim_face_normal(xyz, c, a, b, nullptr, n0);
                if (n0[0] == 0.0 && n0[1] == 0.0 && n0[2] == 0.0)
                    continue;
                decim_face_normal(xyz, c, a, b, placed.mX, n1);
                if (n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] <= 0.0) {
                    flips = true;
                    break;
                }
            }
        }
        if (flips) {
            ++result.mCollapsesRejected;
            continue;
        }

        // --- commit: v merges into the survivor s ---------------------------
        const std::int64_t s = pinned[static_cast<std::size_t>(b)] ? b : a;
        const std::int64_t r = s == a ? b : a;

        // Blend parameter: project the placed point onto the edge, clamped to
        // [0, 1]. Exact for Midpoint/Endpoint, an approximation for Optimal
        // (the optimal point need not lie on the edge). A pinned survivor has
        // placed == its own position, so t == 0 and the data provably keeps the
        // survivor's values.
        {
            const double* xs = xyz.data() + static_cast<std::size_t>(s) * 3;
            const double* xr = xyz.data() + static_cast<std::size_t>(r) * 3;
            const double dx = xr[0] - xs[0];
            const double dy = xr[1] - xs[1];
            const double dz = xr[2] - xs[2];
            const double den = dx * dx + dy * dy + dz * dz;
            double t;
            if (den == 0.0) {
                t = 0.5;
            } else {
                t = ((placed.mX[0] - xs[0]) * dx + (placed.mX[1] - xs[1]) * dy +
                     (placed.mX[2] - xs[2]) * dz) /
                    den;
                if (t < 0.0)
                    t = 0.0;
                if (t > 1.0)
                    t = 1.0;
            }
            for (DecimFloatData& fd : float_data) {
                double* vals = fd.mValues.data();
                const std::size_t so = static_cast<std::size_t>(s) * fd.mNumComponents;
                const std::size_t ro = static_cast<std::size_t>(r) * fd.mNumComponents;
                for (std::size_t k = 0; k < fd.mNumComponents; ++k)
                    vals[so + k] = vals[so + k] + t * (vals[ro + k] - vals[so + k]);
            }
        }

        xyz[static_cast<std::size_t>(s) * 3] = placed.mX[0];
        xyz[static_cast<std::size_t>(s) * 3 + 1] = placed.mX[1];
        xyz[static_cast<std::size_t>(s) * 3 + 2] = placed.mX[2];
        {
            double* qs = quads.data() + static_cast<std::size_t>(s) * 10;
            const double* qr = quads.data() + static_cast<std::size_t>(r) * 10;
            for (int i = 0; i < 10; ++i)
                qs[i] += qr[i];
        }

        // Kill the shared faces (removing them from every corner's list), then
        // migrate the removed vertex's remaining faces onto the survivor.
        for (const std::int64_t f : shared) {
            face_alive[static_cast<std::size_t>(f)] = 0;
            --alive_faces;
            ++result.mFacesRemoved;
            const std::int64_t* c = faces.mCorners.data() + static_cast<std::size_t>(f) * 3;
            for (int k = 0; k < 3; ++k)
                if (c[k] != s && c[k] != r)
                    decim_sorted_erase(vfaces[static_cast<std::size_t>(c[k])], f);
            decim_sorted_erase(vfaces[static_cast<std::size_t>(s)], f);
        }
        for (const std::int64_t f : vfaces[static_cast<std::size_t>(r)]) {
            if (!face_alive[static_cast<std::size_t>(f)])
                continue;  // was in shared
            std::int64_t* c = faces.mCorners.data() + static_cast<std::size_t>(f) * 3;
            for (int k = 0; k < 3; ++k)
                if (c[k] == r)
                    c[k] = s;
            decim_sorted_insert(vfaces[static_cast<std::size_t>(s)], f);
        }
        vfaces[static_cast<std::size_t>(r)].clear();
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
    std::vector<char> used(n, 0);
    for (std::size_t f = 0; f < num_faces; ++f)
        if (face_alive[f])
            for (int k = 0; k < 3; ++k)
                used[static_cast<std::size_t>(faces.mCorners[f * 3 + k])] = 1;
    std::vector<std::int64_t> remap(n, -1);
    std::size_t num_used = 0;
    for (std::size_t i = 0; i < n; ++i)
        if (used[i])
            remap[i] = static_cast<std::int64_t>(num_used++);

    Mesh& out = result.mMesh;

    // Points: the working coordinates of the used vertices, cast back once to
    // the input dtype.
    {
        const NDArray& points = simp.Points();
        NDArray new_points = NDArray::Uninit(points.Dtype(), {num_used, dim});
        detail::dispatch_dtype(points.Dtype(), [&]<class T>() {
            T* dst = new_points.As<T>();
            std::size_t w = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (!used[i])
                    continue;
                for (std::size_t d = 0; d < dim && d < 3; ++d)
                    dst[w * dim + d] = static_cast<T>(xyz[i * 3 + d]);
                ++w;
            }
        });
        out.AssignPoints(std::move(new_points));
    }

    // One triangle block per input block (1:1, possibly empty), plus the
    // per-input-block cell maps composed through the simplexify child map.
    {
        const bool has_parent = simp.HasCellData("convert:parent_cell") &&
                                simp.CellDataNumBlocks("convert:parent_cell") == nblocks;
        std::size_t bi = 0;
        std::vector<std::size_t> in_cells;
        for (const auto cb : rMesh.CellRange())
            in_cells.push_back(cb.NumCells());
        for (const auto cb : simp.CellRange()) {
            const std::size_t nc = cb.NumCells();
            const std::size_t base = static_cast<std::size_t>(faces.mBlockBase[bi]);
            std::size_t alive_in_block = 0;
            for (std::size_t c = 0; c < nc; ++c)
                if (face_alive[base + c])
                    ++alive_in_block;
            NDArray block = NDArray::Uninit(DType::Int64, {alive_in_block, 3});
            std::int64_t* dst = block.As<std::int64_t>();
            const NDArray* p_parent =
                has_parent ? &simp.CellData("convert:parent_cell", bi) : nullptr;
            std::vector<std::int64_t> cell_map(bi < in_cells.size() ? in_cells[bi] : nc, -1);
            std::size_t w = 0;
            for (std::size_t c = 0; c < nc; ++c) {
                if (!face_alive[base + c])
                    continue;
                for (std::size_t k = 0; k < 3; ++k)
                    dst[w * 3 + k] =
                        remap[static_cast<std::size_t>(faces.mCorners[(base + c) * 3 + k])];
                const std::int64_t parent = p_parent != nullptr ? detail::read_int(*p_parent, c)
                                                                : static_cast<std::int64_t>(c);
                if (parent >= 0 && static_cast<std::size_t>(parent) < cell_map.size() &&
                    cell_map[static_cast<std::size_t>(parent)] < 0)
                    cell_map[static_cast<std::size_t>(parent)] = static_cast<std::int64_t>(w);
                ++w;
            }
            out.AddCellBlock(cell_type_name(CellType::Triangle), std::move(block));
            result.mCellMaps.push_back(decim_int64_vector(cell_map));
            ++bi;
        }
    }

    // cell_data: the surviving faces' own rows (simplexify already replicated a
    // parent's row to its triangles); convert:parent_cell is consumed above,
    // never forwarded.
    for (const std::string& name : simp.CellDataNames()) {
        if (name == "convert:parent_cell")
            continue;
        const std::size_t ndata = simp.CellDataNumBlocks(name);
        std::vector<NDArray> blocks;
        blocks.reserve(ndata);
        std::size_t bi = 0;
        for (const auto cb : simp.CellRange()) {
            if (bi >= ndata)
                break;
            const NDArray& a = simp.CellData(name, bi);
            const std::size_t nc = cb.NumCells();
            if (ndata != nblocks || detail::rows(a) != nc || nc == 0) {
                blocks.push_back(detail::data_owned_copy(a));
                ++bi;
                continue;
            }
            const std::size_t base = static_cast<std::size_t>(faces.mBlockBase[bi]);
            std::size_t alive_in_block = 0;
            for (std::size_t c = 0; c < nc; ++c)
                if (face_alive[base + c])
                    ++alive_in_block;
            std::vector<std::size_t> shape = a.Shape();
            shape[0] = alive_in_block;
            const std::size_t row_bytes = a.Nbytes() / nc;
            NDArray sub = NDArray::Uninit(a.Dtype(), std::move(shape));
            const std::byte* src = a.Data();
            std::byte* dst = sub.Data();
            std::size_t w = 0;
            for (std::size_t c = 0; c < nc; ++c)
                if (face_alive[base + c])
                    std::memcpy(dst + (w++) * row_bytes, src + c * row_bytes, row_bytes);
            blocks.push_back(std::move(sub));
            ++bi;
        }
        out.AddCellData(name, std::move(blocks));
    }

    // point_data: float-kind arrays from the blended working copies (cast back
    // to the source dtype), everything else as the survivor's own byte rows.
    for (const std::string& name : simp.PointDataNames()) {
        const NDArray& a = simp.PointData(name);
        if (detail::rows(a) != n) {
            out.AddPointData(name, detail::data_owned_copy(a));
            continue;
        }
        std::vector<std::size_t> shape = a.Shape();
        shape[0] = num_used;
        if (detail::is_float_dtype(a.Dtype())) {
            const DecimFloatData* p_fd = nullptr;
            for (const DecimFloatData& fd : float_data)
                if (fd.mName == name) {
                    p_fd = &fd;
                    break;
                }
            const std::size_t ncomp = p_fd != nullptr ? p_fd->mNumComponents : 0;
            NDArray b = NDArray::Uninit(a.Dtype(), std::move(shape));
            detail::dispatch_dtype(a.Dtype(), [&]<class T>() {
                T* dst = b.As<T>();
                std::size_t w = 0;
                for (std::size_t i = 0; i < n; ++i) {
                    if (!used[i])
                        continue;
                    for (std::size_t k = 0; k < ncomp; ++k)
                        dst[w * ncomp + k] = static_cast<T>(p_fd->mValues[i * ncomp + k]);
                    ++w;
                }
            });
            out.AddPointData(name, std::move(b));
        } else {
            const std::size_t row_bytes = n == 0 ? 0 : a.Nbytes() / n;
            NDArray b = NDArray::Uninit(a.Dtype(), std::move(shape));
            const std::byte* src = a.Data();
            std::byte* dst = b.Data();
            std::size_t w = 0;
            for (std::size_t i = 0; i < n; ++i)
                if (used[i])
                    std::memcpy(dst + (w++) * row_bytes, src + i * row_bytes, row_bytes);
            out.AddPointData(name, std::move(b));
        }
    }
    for (const std::string& name : simp.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(simp.FieldData(name)));

    // mPointMap: every input point lands on its survivor's output index.
    {
        NDArray pm = NDArray::Uninit(DType::Int64, {n});
        std::int64_t* dst = pm.As<std::int64_t>();
        for (std::size_t i = 0; i < n; ++i) {
            std::int64_t v = static_cast<std::int64_t>(i);
            while (successor[static_cast<std::size_t>(v)] >= 0)
                v = successor[static_cast<std::size_t>(v)];
            dst[i] = remap[static_cast<std::size_t>(v)];
        }
        result.mPointMap = std::move(pm);
    }
    result.mPointsRemoved = static_cast<std::int64_t>(n - num_used);

    // Named regions. The cell map is FirstChild — an input quad/polygon becomes
    // a contiguous run of triangles, and a fully-collapsed parent maps to -1,
    // which the run scan steps over. Side regions are dropped: the output is
    // all-triangle, so an input facet has no counterpart. See
    // detail/region_remap.hpp.
    {
        detail::RegionRemap rmap;
        rmap.pPointMap = &result.mPointMap;
        rmap.mCellMapKind = detail::CellMapKind::FirstChild;
        rmap.pCellMaps = &result.mCellMaps;
        rmap.mOpName = "decimate";
        detail::remap_regions(rMesh, result.mMesh, rmap);
    }

    return result;
}

}  // namespace meshioplusplus
