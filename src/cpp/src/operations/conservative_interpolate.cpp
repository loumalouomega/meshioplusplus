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
// Mass-preserving cross-mesh field transfer: both meshes are simplexified,
// overlapping simplex pairs are found via a spatial hash and measured with an
// exact triangle-triangle (2D) or tetra-tetra (3D) clip, and every target
// cell's value is the overlap-measure-weighted mean of the source cells it
// intersects. See operations/conservative_interpolate.hpp for the contract.
// There is no pure-numpy fallback (see the header doc for why); this file is
// therefore the sole reference implementation.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/conservative_interpolate.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kProxySuffix = "__cons_interp_proxy__";

// Effective topological dimension of a block, matching surface.cpp's
// surface_effective_dim: a polyhedron's type name may be parameterized (e.g.
// "polyhedron12") and so miss the cell_type_from_name table.
int cons_effective_dim(const Mesh::CellView& rBlock) {
    if (rBlock.IsPolyhedron())
        return 3;
    return cell_type_dimension(cell_type_from_name(rBlock.Type()));
}

// The mesh's maximum effective cell-block dimension, or -1 if it has none.
int cons_max_dim(const Mesh& rMesh) {
    int max_dim = -1;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.NumCells() == 0)
            continue;
        max_dim = std::max(max_dim, cons_effective_dim(cb));
    }
    return max_dim;
}

// -------------------------------------------------------------------------
// Simplex extraction: flatten every "want"-typed (triangle/tetra) block of a
// simplexified mesh into one flat connectivity array, plus each simplex's
// GLOBAL parent-cell index in the *original* (pre-simplexify) mesh, resolved
// via that mesh's own block_bases + the simplexified mesh's own
// "convert:parent_cell" cell_data (per-output-cell, index within its own
// block -- valid here because convert_cells preserves block structure 1:1).
// -------------------------------------------------------------------------
struct ConsSimplices {
    std::vector<std::int64_t> mConn;    ///< nv per simplex, ascending block-then-cell order
    std::vector<std::int64_t> mParent;  ///< one per simplex: global cell index in the original mesh
    std::size_t mNv = 0;                ///< 3 (triangle) or 4 (tetra)
};

ConsSimplices cons_extract_simplices(const Mesh& rOriginal, const Mesh& rSimplexMesh,
                                     bool wantTet) {
    ConsSimplices out;
    out.mNv = wantTet ? 4 : 3;
    const std::string want = wantTet ? "tetra" : "triangle";
    const std::vector<std::int64_t> bases = detail::block_bases(rOriginal);
    const bool has_parent = rSimplexMesh.HasCellData("convert:parent_cell");

    std::size_t block_idx = 0;
    for (const auto cb : rSimplexMesh.CellRange()) {
        if (cb.Type() == want && cb.NumCells() > 0) {
            const NDArray& conn = cb.Conn();
            const std::size_t n = cb.NumCells() * out.mNv;
            const std::size_t base = out.mConn.size();
            out.mConn.resize(base + n);
            for (std::size_t j = 0; j < n; ++j)
                out.mConn[base + j] = detail::read_int(conn, j);

            const std::size_t pbase = out.mParent.size();
            out.mParent.resize(pbase + cb.NumCells());
            if (has_parent) {
                const NDArray& pc = rSimplexMesh.CellData("convert:parent_cell", block_idx);
                for (std::size_t c = 0; c < cb.NumCells(); ++c)
                    out.mParent[pbase + c] = bases[block_idx] + detail::read_int(pc, c);
            } else {
                // 1:1 pass-through block (no children were ever created for
                // it) -- the local row index is already the parent's own.
                for (std::size_t c = 0; c < cb.NumCells(); ++c)
                    out.mParent[pbase + c] = bases[block_idx] + static_cast<std::int64_t>(c);
            }
        }
        ++block_idx;
    }
    return out;
}

// -------------------------------------------------------------------------
// Small vector helpers on top of detail::Vec3 -- plain doubles, since the
// clip kernel never needs anything the mesh's own dtype provides.
// -------------------------------------------------------------------------

double cons_signed_tet_volume(const detail::Vec3 v[4]) {
    return detail::triple_product(detail::vec3_sub(v[1], v[0]), detail::vec3_sub(v[2], v[0]),
                                  detail::vec3_sub(v[3], v[0])) /
           6.0;
}

// Positively-oriented copy of a tetra's 4 corners (swap the first two if the
// signed volume came out negative). The clip kernel's face table assumes
// this orientation, so this removes any dependence on convert_cells'
// simplexification always producing positive volumes.
std::array<detail::Vec3, 4> cons_oriented_tet(const detail::Vec3 v[4]) {
    if (cons_signed_tet_volume(v) < 0.0)
        return {v[1], v[0], v[2], v[3]};
    return {v[0], v[1], v[2], v[3]};
}

// Sutherland-Hodgman clip of a (planar or non-planar-but-convex) vertex ring
// against the half-space {p : dot(p - planePoint, planeNormal) <= 0}, i.e.
// "inside" is the non-positive side. Every newly interpolated (on-plane)
// vertex is additionally appended to *pNewPts when it is not null -- the
// mechanism the 3D tetra-tetra clip uses to collect the chord bounding its
// new capping face. Works unmodified for the 2D triangle-triangle clip too
// (called with pNewPts == nullptr there, since a planar 2D clip needs no
// capping step).
std::vector<detail::Vec3> cons_clip_ring(const std::vector<detail::Vec3>& rPoly,
                                         const detail::Vec3& rPlanePoint,
                                         const detail::Vec3& rPlaneNormal,
                                         std::vector<detail::Vec3>* pNewPts) {
    if (rPoly.empty())
        return {};
    const std::size_t n = rPoly.size();
    std::vector<detail::Vec3> out;
    out.reserve(n + 1);
    auto s = [&](const detail::Vec3& p) {
        return detail::vec3_dot(detail::vec3_sub(p, rPlanePoint), rPlaneNormal);
    };
    for (std::size_t i = 0; i < n; ++i) {
        const detail::Vec3& cur = rPoly[i];
        const detail::Vec3& prev = rPoly[(i + n - 1) % n];
        const double s_cur = s(cur);
        const double s_prev = s(prev);
        const bool cur_in = s_cur <= 0.0;
        const bool prev_in = s_prev <= 0.0;
        if (cur_in != prev_in) {
            const double t = s_prev / (s_prev - s_cur);
            const detail::Vec3 ip =
                detail::vec3_add(prev, detail::vec3_scale(detail::vec3_sub(cur, prev), t));
            out.push_back(ip);
            if (pNewPts != nullptr)
                pNewPts->push_back(ip);
        }
        if (cur_in)
            out.push_back(cur);
    }
    return out;
}

// Exact overlap area of two triangles, in the xy-plane (z ignored) -- the
// same planar assumption interpolate.cpp's Barycentric triangle path already
// documents. Both inputs are first made CCW so the source triangle's three
// edges define three consistent "inside" half-planes.
double cons_clip_triangle_triangle_area(const detail::Vec3 tgtIn[3], const detail::Vec3 srcIn[3]) {
    auto ccw = [](const detail::Vec3 t[3]) -> std::array<detail::Vec3, 3> {
        const double a2 =
            (t[1][0] - t[0][0]) * (t[2][1] - t[0][1]) - (t[1][1] - t[0][1]) * (t[2][0] - t[0][0]);
        const detail::Vec3 p0{t[0][0], t[0][1], 0.0};
        const detail::Vec3 p1{t[1][0], t[1][1], 0.0};
        const detail::Vec3 p2{t[2][0], t[2][1], 0.0};
        if (a2 < 0.0)
            return {p0, p2, p1};
        return {p0, p1, p2};
    };
    const std::array<detail::Vec3, 3> tgt = ccw(tgtIn);
    const std::array<detail::Vec3, 3> src = ccw(srcIn);

    std::vector<detail::Vec3> poly{tgt[0], tgt[1], tgt[2]};
    for (int e = 0; e < 3 && !poly.empty(); ++e) {
        const detail::Vec3& a = src[static_cast<std::size_t>(e)];
        const detail::Vec3& b = src[static_cast<std::size_t>((e + 1) % 3)];
        // Left of a->b (CCW inside) written as a plane test: s <= 0 iff
        // cross_z(b - a, p - a) >= 0.
        const detail::Vec3 normal{b[1] - a[1], a[0] - b[0], 0.0};
        poly = cons_clip_ring(poly, a, normal, nullptr);
    }
    if (poly.size() < 3)
        return 0.0;
    return detail::polygon_area(poly.data(), poly.size());
}

// Exact overlap volume of two tetrahedra. Both operands are always
// tetrahedra (guaranteed by convert_cells(Simplexify)), which bounds this to
// a genuine but tractable 3D convex-polytope clip: the source tet's own 4
// faces are clipped against the target tet's 4 half-spaces in turn, and each
// clip's new boundary chord is capped with a fan-triangulated, angle-sorted
// polygon -- never a general polyhedron-polyhedron clipper. See the header
// doc comment for the full algorithm description.
double cons_clip_tetra_tetra_volume(const detail::Vec3 tgtIn[4], const detail::Vec3 srcIn[4]) {
    // Recentre on the source tet's own corner average before any arithmetic
    // -- poly_measure's own numerical-stability lesson: this keeps the final
    // divergence-theorem sum (and every plane-classification dot product
    // along the way) accurate for a mesh far from the origin.
    detail::Vec3 c{0.0, 0.0, 0.0};
    for (int i = 0; i < 4; ++i)
        c = detail::vec3_add(c, srcIn[i]);
    c = detail::vec3_scale(c, 0.25);

    detail::Vec3 src_raw[4];
    detail::Vec3 tgt_raw[4];
    for (int i = 0; i < 4; ++i) {
        src_raw[i] = detail::vec3_sub(srcIn[i], c);
        tgt_raw[i] = detail::vec3_sub(tgtIn[i], c);
    }
    const std::array<detail::Vec3, 4> src = cons_oriented_tet(src_raw);
    const std::array<detail::Vec3, 4> tgt = cons_oriented_tet(tgt_raw);

    // Dedup tolerance for chord endpoints, scaled to the source tet's own
    // size (recentring only shifted the origin, not the physical scale).
    double scale = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j)
            scale = std::max(scale,
                             detail::vec3_norm(detail::vec3_sub(src[static_cast<std::size_t>(i)],
                                                                src[static_cast<std::size_t>(j)])));
    const double eps_abs = 1e-9 * scale;
    const double eps_abs2 = eps_abs * eps_abs;

    // The outward-wound tetra face table (detail/cell_faces.hpp's own rows),
    // valid for a positively-oriented tetra.
    static constexpr int kFaces[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};

    std::vector<std::array<detail::Vec3, 3>> tris;
    tris.reserve(4);
    for (const auto& f : kFaces)
        tris.push_back({src[static_cast<std::size_t>(f[0])], src[static_cast<std::size_t>(f[1])],
                        src[static_cast<std::size_t>(f[2])]});

    for (const auto& f : kFaces) {
        if (tris.empty())
            break;
        const detail::Vec3& p0 = tgt[static_cast<std::size_t>(f[0])];
        const detail::Vec3& p1 = tgt[static_cast<std::size_t>(f[1])];
        const detail::Vec3& p2 = tgt[static_cast<std::size_t>(f[2])];
        const detail::Vec3 normal =
            detail::vec3_cross(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0));

        std::vector<std::array<detail::Vec3, 3>> kept;
        std::vector<detail::Vec3> chord;
        for (const auto& t : tris) {
            std::vector<detail::Vec3> new_pts;
            const std::vector<detail::Vec3> ring =
                cons_clip_ring({t[0], t[1], t[2]}, p0, normal, &new_pts);
            for (const detail::Vec3& p : new_pts)
                chord.push_back(p);
            for (std::size_t i = 1; i + 1 < ring.size(); ++i)
                kept.push_back({ring[0], ring[i], ring[i + 1]});
        }
        tris = std::move(kept);

        // Cap the new hole (if this half-space actually cut anything) with
        // the dedup'd chord, angle-sorted around the cutting plane's own 2D
        // basis and fan-triangulated -- always valid, since the cap of a
        // convex polytope clipped by one more plane is itself convex.
        std::vector<detail::Vec3> uniq;
        for (const detail::Vec3& p : chord) {
            bool dup = false;
            for (const detail::Vec3& q : uniq) {
                if (detail::vec3_norm_sq(detail::vec3_sub(p, q)) < eps_abs2) {
                    dup = true;
                    break;
                }
            }
            if (!dup)
                uniq.push_back(p);
        }
        if (uniq.size() >= 3) {
            detail::Vec3 cc{0.0, 0.0, 0.0};
            for (const detail::Vec3& p : uniq)
                cc = detail::vec3_add(cc, p);
            cc = detail::vec3_scale(cc, 1.0 / static_cast<double>(uniq.size()));

            const detail::Vec3 w = detail::vec3_normalize(normal);
            const detail::Vec3 arbitrary =
                std::abs(w[0]) < 0.9 ? detail::Vec3{1.0, 0.0, 0.0} : detail::Vec3{0.0, 1.0, 0.0};
            const detail::Vec3 u = detail::vec3_normalize(
                detail::vec3_sub(arbitrary, detail::vec3_scale(w, detail::vec3_dot(arbitrary, w))));
            const detail::Vec3 v = detail::vec3_cross(w, u);

            std::vector<std::pair<double, std::size_t>> order(uniq.size());
            for (std::size_t i = 0; i < uniq.size(); ++i) {
                const detail::Vec3 d = detail::vec3_sub(uniq[i], cc);
                order[i] = {std::atan2(detail::vec3_dot(d, v), detail::vec3_dot(d, u)), i};
            }
            // Explicit, order-independent tie-break: ascending build index on
            // an exact angle tie (near-coplanar chord points), never
            // insertion order from an unordered container.
            std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) {
                if (a.first != b.first)
                    return a.first < b.first;
                return a.second < b.second;
            });
            for (std::size_t i = 1; i + 1 < order.size(); ++i)
                tris.push_back(
                    {uniq[order[0].second], uniq[order[i].second], uniq[order[i + 1].second]});
        }
    }

    double vol6 = 0.0;
    for (const auto& t : tris)
        vol6 += detail::triple_product(t[0], t[1], t[2]);
    const double vol = vol6 / 6.0;
    return vol > 0.0 ? vol : 0.0;
}

// -------------------------------------------------------------------------
// The shared geometric overlap structure -- purely a function of the two
// meshes' geometry, independent of every array. Built once per call and
// reused for every requested array (direct cell_data and point-composed
// proxies alike).
// -------------------------------------------------------------------------
struct ConsOverlap {
    // Per target simplex: (source simplex id, overlap measure), deduplicated
    // and sorted ascending by source simplex id (broad-phase determinism).
    std::vector<std::vector<std::pair<std::int64_t, double>>> mBySimplex;
    ConsSimplices mSrc;
    ConsSimplices mTgt;
};

// A simplex's own (unclipped) measure, from its own corner coordinates.
double cons_own_measure(const detail::Vec3* pCorners, std::size_t nv) {
    if (nv == 3)
        return detail::polygon_area(pCorners, 3);
    detail::Vec3 v[4] = {pCorners[0], pCorners[1], pCorners[2], pCorners[3]};
    const double vol = detail::cell_volume_from_corners(v, CellType::Tetra);
    return std::abs(vol);
}

ConsOverlap cons_build_overlap(const Mesh& rSource, const Mesh& rTarget, const Mesh& rSourceSimp,
                               const Mesh& rTargetSimp, bool wantTet) {
    ConsOverlap out;
    out.mSrc = cons_extract_simplices(rSource, rSourceSimp, wantTet);
    out.mTgt = cons_extract_simplices(rTarget, rTargetSimp, wantTet);
    const std::size_t nv = out.mSrc.mNv;
    const std::size_t nsrc = out.mSrc.mParent.size();
    const std::size_t ntgt = out.mTgt.mParent.size();
    out.mBySimplex.resize(ntgt);
    if (nsrc == 0 || ntgt == 0)
        return out;

    const NDArray& spts = rSourceSimp.Points();
    const std::size_t sdim = rSourceSimp.PointDim();
    const NDArray& tpts = rTargetSimp.Points();
    const std::size_t tdim = rTargetSimp.PointDim();

    auto scoord = [&](std::size_t simp, std::size_t corner) {
        return detail::read_point(spts, sdim, out.mSrc.mConn[simp * nv + corner]);
    };
    auto tcoord = [&](std::size_t simp, std::size_t corner) {
        return detail::read_point(tpts, tdim, out.mTgt.mConn[simp * nv + corner]);
    };

    // Broad phase: quantize every source simplex's own bounding box into a
    // shared grid (interp_cell_size's cbrt-free sizing rule, kept local since
    // anon-namespace helpers cannot be shared across files).
    double lo[3], hi[3];
    {
        const detail::Vec3 first = scoord(0, 0);
        for (int d = 0; d < 3; ++d)
            lo[d] = hi[d] = first[d];
        for (std::size_t s = 0; s < nsrc; ++s)
            for (std::size_t k = 0; k < nv; ++k) {
                const detail::Vec3 p = scoord(s, k);
                for (int d = 0; d < 3; ++d) {
                    if (p[d] < lo[d])
                        lo[d] = p[d];
                    if (p[d] > hi[d])
                        hi[d] = p[d];
                }
            }
    }
    double max_extent = 0.0;
    for (int d = 0; d < 3; ++d)
        max_extent = std::max(max_extent, hi[d] - lo[d]);
    std::int64_t r = 1;
    while (r * r * r < static_cast<std::int64_t>(nsrc))
        ++r;
    const double cell = max_extent > 0.0 ? max_extent / static_cast<double>(r) : 1.0;

    detail::SpatialGrid grid(cell);
    std::vector<std::pair<detail::GridKey, detail::GridKey>> src_boxes(nsrc);
    parallel_for(nsrc, [&](std::size_t s) {
        detail::Vec3 blo = scoord(s, 0);
        detail::Vec3 bhi = blo;
        for (std::size_t k = 1; k < nv; ++k) {
            const detail::Vec3 p = scoord(s, k);
            for (int d = 0; d < 3; ++d) {
                if (p[d] < blo[d])
                    blo[d] = p[d];
                if (p[d] > bhi[d])
                    bhi[d] = p[d];
            }
        }
        src_boxes[s] = {grid.KeyOf(blo.data()), grid.KeyOf(bhi.data())};
    });
    for (std::size_t s = 0; s < nsrc; ++s)
        grid.InsertBox(src_boxes[s].first, src_boxes[s].second, static_cast<std::int64_t>(s));

    // Narrow phase, independent per target simplex -> parallel-safe.
    parallel_for(ntgt, [&](std::size_t t) {
        detail::Vec3 t_corners[4];
        for (std::size_t k = 0; k < nv; ++k)
            t_corners[k] = tcoord(t, k);
        detail::Vec3 tblo = t_corners[0];
        detail::Vec3 tbhi = t_corners[0];
        for (std::size_t k = 1; k < nv; ++k)
            for (int d = 0; d < 3; ++d) {
                if (t_corners[k][d] < tblo[d])
                    tblo[d] = t_corners[k][d];
                if (t_corners[k][d] > tbhi[d])
                    tbhi[d] = t_corners[k][d];
            }

        std::vector<std::int64_t> cand;
        grid.ForEachInBox(grid.KeyOf(tblo.data()), grid.KeyOf(tbhi.data()),
                          [&](const std::vector<std::int64_t>& rIds) {
                              cand.insert(cand.end(), rIds.begin(), rIds.end());
                          });
        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

        auto& result = out.mBySimplex[t];
        result.reserve(cand.size());
        for (std::int64_t s : cand) {
            detail::Vec3 s_corners[4];
            for (std::size_t k = 0; k < nv; ++k)
                s_corners[k] = scoord(static_cast<std::size_t>(s), k);
            const double measure = nv == 3 ? cons_clip_triangle_triangle_area(t_corners, s_corners)
                                           : cons_clip_tetra_tetra_volume(t_corners, s_corners);
            if (measure > 0.0)
                result.push_back({s, measure});
        }
    });

    return out;
}

// -------------------------------------------------------------------------
// Accumulate one requested array's conservative remap: reads component
// values from `rValueSource.CellData(rSourceName, block)` (resolved through
// the source simplex's own global-parent index, indexing the block the
// *original* source mesh assigns that parent to), scatters
// `value * measure` into per-target-cell accumulators using `rOverlap`, and
// returns one Float64 NDArray per target cell block.
// -------------------------------------------------------------------------
std::vector<NDArray> cons_accumulate_array(const Mesh& rTarget, const ConsOverlap& rOverlap,
                                           const Mesh& rValueSource, const std::string& rSourceName,
                                           const std::vector<double>& rCoveredMeasure,
                                           const std::vector<double>& rOwnMeasure, double eps_rel,
                                           double defaultValue, std::vector<char>& rUncovered) {
    const std::size_t comps = data_num_components(rValueSource.CellData(rSourceName, 0));
    const std::vector<std::int64_t> src_bases = detail::block_bases(rValueSource);
    const std::size_t total_tgt =
        static_cast<std::size_t>(detail::total_cells(detail::block_bases(rTarget)));

    std::vector<double> numerator(total_tgt * comps, 0.0);

    // Serial scatter (target simplex -> its original target cell), ascending
    // target-simplex index -- floating-point addition is not associative, so
    // this many-to-few fold must not be parallelised (data_average.cpp's
    // cell-to-point averaging documents the identical rule).
    for (std::size_t t = 0; t < rOverlap.mBySimplex.size(); ++t) {
        const std::int64_t parent = rOverlap.mTgt.mParent[t];
        double* out_row = numerator.data() + static_cast<std::size_t>(parent) * comps;
        for (const auto& sm : rOverlap.mBySimplex[t]) {
            const std::int64_t src_parent =
                rOverlap.mSrc.mParent[static_cast<std::size_t>(sm.first)];
            const auto br = detail::global_to_block_row(src_bases, src_parent);
            const NDArray& src_arr = rValueSource.CellData(rSourceName, br.first);
            for (std::size_t c = 0; c < comps; ++c)
                out_row[c] +=
                    sm.second *
                    detail::read_double(src_arr, static_cast<std::size_t>(br.second) * comps + c);
        }
    }

    if (rUncovered.empty())
        rUncovered.assign(total_tgt, 0);
    for (std::size_t g = 0; g < total_tgt; ++g) {
        const bool covered = rOwnMeasure[g] > 0.0 && rCoveredMeasure[g] / rOwnMeasure[g] >= eps_rel;
        if (!covered) {
            rUncovered[g] = 1;
            for (std::size_t c = 0; c < comps; ++c)
                numerator[g * comps + c] = defaultValue;
        } else {
            for (std::size_t c = 0; c < comps; ++c)
                numerator[g * comps + c] /= rCoveredMeasure[g];
        }
    }

    std::vector<NDArray> out_blocks;
    out_blocks.reserve(rTarget.NumCellBlocks());
    std::size_t base = 0;
    const NDArray& sample = rValueSource.CellData(rSourceName, 0);
    std::vector<std::size_t> shape = sample.Shape();
    for (const auto cb : rTarget.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        std::vector<std::size_t> block_shape = shape;
        if (block_shape.empty())
            block_shape = {ncells};
        else
            block_shape[0] = ncells;
        NDArray arr = NDArray::Uninit(DType::Float64, block_shape);
        double* o = arr.As<double>();
        for (std::size_t c = 0; c < ncells * comps; ++c)
            o[c] = numerator[base * comps + c];
        out_blocks.push_back(std::move(arr));
        base += ncells;
    }
    return out_blocks;
}

// The standard unknown-key diagnostic, covering both source locations.
std::string cons_unknown_array_message(const Mesh& rSource, const std::string& rName) {
    auto join = [](const std::vector<std::string>& rNames) -> std::string {
        if (rNames.empty())
            return "none";
        std::string s;
        for (std::size_t i = 0; i < rNames.size(); ++i) {
            if (i > 0)
                s += ", ";
            s += rNames[i];
        }
        return s;
    };
    return "meshio++: conservative_interpolate: no source point_data or cell_data array named '" +
           rName + "' (available point_data: " + join(data_names(rSource, DataLocation::Point)) +
           "; cell_data: " + join(data_names(rSource, DataLocation::Cell)) + ")";
}

std::string cons_resolve_name(const Mesh& rTarget, bool isPointData, const std::string& rName,
                              ConservativeInterpolateConflict Conflict) {
    const bool exists = isPointData ? rTarget.HasPointData(rName) : rTarget.HasCellData(rName);
    if (!exists)
        return rName;
    switch (Conflict) {
        case ConservativeInterpolateConflict::Overwrite:
            return rName;
        case ConservativeInterpolateConflict::Suffix: {
            const std::string suffixed = rName + "_interp";
            const bool taken =
                isPointData ? rTarget.HasPointData(suffixed) : rTarget.HasCellData(suffixed);
            if (taken)
                throw std::invalid_argument("meshio++: conservative_interpolate: array '" + rName +
                                            "' already exists on the target and so does '" +
                                            suffixed + "' (on_conflict='suffix')");
            return suffixed;
        }
        case ConservativeInterpolateConflict::Error:
        default:
            throw std::invalid_argument("meshio++: conservative_interpolate: array '" + rName +
                                        "' already exists on the target (on_conflict='error'; pass "
                                        "'overwrite' or 'suffix')");
    }
}

}  // namespace

ConservativeInterpolateConflict conservative_interpolate_conflict_from_name(
    const std::string& rName) {
    if (rName == "error")
        return ConservativeInterpolateConflict::Error;
    if (rName == "overwrite")
        return ConservativeInterpolateConflict::Overwrite;
    if (rName == "suffix")
        return ConservativeInterpolateConflict::Suffix;
    throw std::invalid_argument("meshio++: conservative_interpolate: unknown on_conflict '" +
                                rName + "' (expected 'error', 'overwrite' or 'suffix')");
}

Mesh conservative_interpolate(const Mesh& rSource, const Mesh& rTarget,
                              const ConservativeInterpolateOptions& rOptions) {
    if (rSource.NumPoints() == 0)
        throw std::invalid_argument(
            "meshio++: conservative_interpolate: the source mesh has no points");

    const int src_dim = cons_max_dim(rSource);
    const int tgt_dim = cons_max_dim(rTarget);
    if (src_dim != tgt_dim || (src_dim != 2 && src_dim != 3))
        throw std::invalid_argument(
            "meshio++: conservative_interpolate: source and target must share the same maximum "
            "topological dimension (2 or 3); got source=" +
            std::to_string(src_dim) + ", target=" + std::to_string(tgt_dim));
    const bool want_tet = src_dim == 3;

    // --- resolve which arrays transfer (validation before any work) --------
    std::vector<std::string> point_names;
    std::vector<std::string> cell_names;
    if (rOptions.mArrays.empty()) {
        point_names = rSource.PointDataNames();
        cell_names = rSource.CellDataNames();
    } else {
        for (const std::string& name : rOptions.mArrays) {
            bool any = false;
            if (rSource.HasPointData(name)) {
                point_names.push_back(name);
                any = true;
            }
            if (rSource.HasCellData(name)) {
                cell_names.push_back(name);
                any = true;
            }
            if (!any)
                throw std::invalid_argument(cons_unknown_array_message(rSource, name));
        }
    }
    for (const std::string& name : cell_names) {
        if (rSource.CellDataNumBlocks(name) != rSource.NumCellBlocks())
            throw std::invalid_argument("meshio++: conservative_interpolate: source cell_data '" +
                                        name + "' does not have one array per cell block");
    }

    std::vector<std::string> point_out(point_names.size());
    for (std::size_t i = 0; i < point_names.size(); ++i)
        point_out[i] = cons_resolve_name(rTarget, true, point_names[i], rOptions.mOnConflict);
    std::vector<std::string> cell_out(cell_names.size());
    for (std::size_t i = 0; i < cell_names.size(); ++i)
        cell_out[i] = cons_resolve_name(rTarget, false, cell_names[i], rOptions.mOnConflict);

    if (point_names.empty() && cell_names.empty())
        return detail::clone_mesh(rTarget);

    // --- simplexify both meshes once -----------------------------------------
    ConvertCellsOptions cco;
    cco.mMode = ConvertCellsMode::Simplexify;
    cco.mRecordParentIds = true;
    const ConvertCellsResult src_simp = convert_cells(rSource, cco);
    const ConvertCellsResult tgt_simp = convert_cells(rTarget, cco);
    const std::string want = want_tet ? "tetra" : "triangle";
    auto has_simplex = [&](const Mesh& rMesh) {
        for (const auto cb : rMesh.CellRange())
            if (cb.Type() == want && cb.NumCells() > 0)
                return true;
        return false;
    };
    if (!has_simplex(src_simp.mMesh) || !has_simplex(tgt_simp.mMesh))
        throw std::invalid_argument(
            "meshio++: conservative_interpolate: source and target both need at least one " + want +
            " cell after simplexification");

    // --- build the shared geometric overlap structure once ------------------
    const ConsOverlap overlap =
        cons_build_overlap(rSource, rTarget, src_simp.mMesh, tgt_simp.mMesh, want_tet);

    const std::size_t total_tgt =
        static_cast<std::size_t>(detail::total_cells(detail::block_bases(rTarget)));
    std::vector<double> covered_measure(total_tgt, 0.0);
    std::vector<double> own_measure(total_tgt, 0.0);
    {
        const NDArray& tpts = tgt_simp.mMesh.Points();
        const std::size_t tdim = tgt_simp.mMesh.PointDim();
        const std::size_t nv = overlap.mTgt.mNv;
        for (std::size_t t = 0; t < overlap.mTgt.mParent.size(); ++t) {
            const std::int64_t parent = overlap.mTgt.mParent[t];
            detail::Vec3 corners[4];
            for (std::size_t k = 0; k < nv; ++k)
                corners[k] = detail::read_point(tpts, tdim, overlap.mTgt.mConn[t * nv + k]);
            own_measure[static_cast<std::size_t>(parent)] += cons_own_measure(corners, nv);
            double m = 0.0;
            for (const auto& sm : overlap.mBySimplex[t])
                m += sm.second;
            covered_measure[static_cast<std::size_t>(parent)] += m;
        }
    }
    constexpr double kCoverageEpsRel = 1e-9;

    Mesh out = detail::clone_mesh(rTarget);

    // --- cell_data: direct -----------------------------------------------
    std::vector<char> uncovered;
    for (std::size_t i = 0; i < cell_names.size(); ++i) {
        std::vector<NDArray> blocks =
            cons_accumulate_array(rTarget, overlap, rSource, cell_names[i], covered_measure,
                                  own_measure, kCoverageEpsRel, rOptions.mDefaultValue, uncovered);
        out.AddCellData(cell_out[i], std::move(blocks));
    }

    // --- point_data: composed via a lumped cell proxy -----------------------
    if (!point_names.empty()) {
        DataAverageOptions pd2c;
        pd2c.names = point_names;
        pd2c.suffix = kProxySuffix;
        pd2c.overwrite = true;
        const Mesh src_proxy = point_data_to_cell_data(rSource, pd2c);

        Mesh tgt_proxy_scratch = detail::clone_geometry(rTarget);
        for (std::size_t i = 0; i < point_names.size(); ++i) {
            const std::string proxy_name = point_names[i] + kProxySuffix;
            std::vector<NDArray> blocks = cons_accumulate_array(
                rTarget, overlap, src_proxy, proxy_name, covered_measure, own_measure,
                kCoverageEpsRel, rOptions.mDefaultValue, uncovered);
            tgt_proxy_scratch.AddCellData(proxy_name, std::move(blocks));
        }

        DataAverageOptions c2pd;
        for (const std::string& name : point_names)
            c2pd.names.push_back(name + kProxySuffix);
        c2pd.weight = CellPointWeight::Measure;
        c2pd.suffix = "";
        c2pd.overwrite = true;
        const Mesh tgt_points = cell_data_to_point_data(tgt_proxy_scratch, c2pd);
        for (std::size_t i = 0; i < point_names.size(); ++i)
            out.AddPointData(
                point_out[i],
                detail::data_owned_copy(tgt_points.PointData(point_names[i] + kProxySuffix)));
    }

    const std::size_t n_uncovered =
        static_cast<std::size_t>(std::count(uncovered.begin(), uncovered.end(), char{1}));
    if (n_uncovered > 0)
        log::warn(
            "meshio++: conservative_interpolate: {} of {} target cells had no source overlap; "
            "filled with default_value = {}",
            n_uncovered, total_tgt, rOptions.mDefaultValue);

    return out;
}

}  // namespace meshioplusplus
