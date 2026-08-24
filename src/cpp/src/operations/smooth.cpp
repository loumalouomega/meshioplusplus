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
// Laplacian / Taubin mesh smoothing. A pure coordinate move — connectivity,
// cell_data, field_data and point_data values are carried through unchanged —
// built entirely through the uniform mesh API so it compiles under every mesh
// backend. See operations/smooth.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Project includes
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

// --- small helpers ----------------------------------------------------------

// A deep copy of an NDArray that always owns its buffer (the source may be a
// view over foreign memory, e.g. numpy memory on the write path).
NDArray smooth_owned_copy(const NDArray& rArr) {
    NDArray c = rArr;
    c.MakeOwned();
    return c;
}

// --- resolved iteration parameters ------------------------------------------

struct SmoothParams {
    SmoothMethod mMethod = SmoothMethod::Laplacian;
    double mLambda = 0.0;
    double mMu = 0.0;
    int mNumPasses = 0;
    // Retained alongside mMethod (never a substitute for a switch on it):
    // Taubin's every-other-pass +lambda/-mu alternation is the one piece of
    // per-pass logic still cheaper as a bool than a re-dispatch on mMethod
    // inside the hot loop -- see phase 4 below.
    bool mTaubin = false;
};

// Validate the options and expand the method-dependent lambda sentinel.
//
// Every branch here is an explicit switch (mMethod) with NO default: case --
// the RemeshMetric precedent from last cycle, adopted because the previous
// shape (a single `p.mTaubin = mMethod == Taubin` bool, with "not Taubin"
// standing in for Laplacian everywhere downstream) would have made a third
// enumerator silently take the Laplacian branch: plausible output, wrong
// method, no diagnostic. -Wswitch now catches the next enumerator at compile
// time instead.
SmoothParams smooth_resolve_params(const SmoothOptions& rOptions) {
    SmoothParams p;
    p.mMethod = rOptions.mMethod;
    p.mTaubin = rOptions.mMethod == SmoothMethod::Taubin;
    p.mLambda = rOptions.mLambda;
    if (p.mLambda < 0.0) {
        switch (rOptions.mMethod) {
            case SmoothMethod::Laplacian:
                p.mLambda = 0.5;
                break;
            case SmoothMethod::Taubin:
                p.mLambda = 0.33;
                break;
            case SmoothMethod::Odt:
                p.mLambda = 0.9;
                break;
        }  // the sentinel: each method's own default
    }
    if (!(p.mLambda > 0.0 && p.mLambda < 1.0))
        throw std::invalid_argument("meshio++: smooth: lambda must lie in (0, 1); got " +
                                    std::to_string(p.mLambda));

    if (rOptions.mMethod == SmoothMethod::Taubin) {
        p.mMu = rOptions.mMu;
        // mu < -lambda < 0 is what makes the pass pair a low-pass filter rather
        // than an amplifier; without it Taubin diverges instead of preserving.
        if (!(p.mMu < -p.mLambda))
            throw std::invalid_argument(
                "meshio++: smooth: taubin requires mu < -lambda < 0; got mu=" +
                std::to_string(p.mMu) + ", lambda=" + std::to_string(p.mLambda));
    }
    // Laplacian and Odt silently ignore mu (documented in SmoothOptions::mMu):
    // making only the newest method strict about an already-established
    // silently-ignored field would be a worse inconsistency than either rule
    // applied uniformly.

    const int iterations = rOptions.mIterations > 0 ? rOptions.mIterations : 0;
    switch (rOptions.mMethod) {
        case SmoothMethod::Taubin:
            p.mNumPasses = iterations * 2;  // +lambda then mu, per iteration
            break;
        case SmoothMethod::Laplacian:
        case SmoothMethod::Odt:
            p.mNumPasses = iterations;
            break;
    }
    return p;
}

// --- flat coordinate buffer -------------------------------------------------

// Points as a flat (n, 3) double buffer, z-padded for a 2D mesh. Everything
// downstream reads coordinates from here rather than from the NDArray, so no
// hot loop pays a dtype dispatch or touches a backend accessor.
std::vector<double> smooth_read_coords(const Mesh& rMesh, std::size_t n, std::size_t dim) {
    std::vector<double> xyz(n * 3, 0.0);
    if (n == 0 || dim == 0)
        return xyz;
    const NDArray& points = rMesh.Points();
    parallel_for_bw(n, [&](std::size_t i) {
        for (std::size_t d = 0; d < dim && d < 3; ++d)
            xyz[i * 3 + d] = detail::read_double(points, i * dim + d);
    });
    return xyz;
}

// Write the leading `dim` columns back out, preserving the source dtype.
// Mirrors transform.cpp's transform_apply_points.
NDArray smooth_write_coords(const NDArray& rPoints, const std::vector<double>& rXyz, std::size_t n,
                            std::size_t dim) {
    NDArray out = NDArray::Uninit(rPoints.Dtype(), {n, dim});
    if (n == 0 || dim == 0)
        return out;
    detail::dispatch_dtype(rPoints.Dtype(), [&]<class T>() {
        T* dst = out.As<T>();
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t d = 0; d < dim && d < 3; ++d)
                dst[i * dim + d] = static_cast<T>(rXyz[i * 3 + d]);
        });
    });
    return out;
}

// --- measurable-cell table (for the inversion guard) ------------------------

// How a cell's orientation can be checked. Cells with no signed measure at all
// never enter the table, so `None` is not represented here.
enum class SmoothMeasure : std::uint8_t {
    PolyFan,     ///< polyhedron: the kernel's corner-average fan over the cell's own faces
    FaceFan,     ///< 3D volume cell: signed volume by the outward face fan.
    Shoelace2D,  ///< 2D cell in a 2D mesh: signed area.
    NormalFlip,  ///< 2D cell in a 3D mesh: did the facet normal fold over?
};

// Only cells whose orientation can actually be checked are recorded; the guard
// has nothing to say about the rest, so leaving them out shrinks both this
// table and the incidence CSR built from it.
struct SmoothCellTable {
    std::vector<std::int64_t> mCornerOffset;  // size numMeasurable + 1
    std::vector<std::int64_t> mCornerNodes;   // flat corner ids
    std::vector<SmoothMeasure> mMeasure;      // size numMeasurable
    std::vector<std::int32_t> mFaceTable;     // index into mFaceTables, or -1
    std::vector<const std::vector<detail::CellFaceDef>*> mFaceTables;
    // Polyhedron cells only. A polyhedron's faces belong to the CELL, not to
    // its type, so they cannot share a per-type table like the others: each
    // cell's rings are appended here and located by mPolyStart/mPolyNumFaces.
    std::vector<std::int64_t> mPolyStart;       // per cell: index into mPolyFaceStart, or -1
    std::vector<std::uint32_t> mPolyNumFaces;   // per cell: 0 unless PolyFan
    std::vector<std::uint32_t> mPolyFaceStart;  // concatenated (nfaces + 1) runs
    std::vector<std::uint32_t> mPolyFaceNodes;  // concatenated local ring indices

    std::size_t NumCells() const { return mMeasure.size(); }
};

SmoothCellTable smooth_build_cell_table(const Mesh& rMesh, std::size_t n, bool is2d) {
    SmoothCellTable t;
    t.mCornerOffset.push_back(0);
    std::unordered_map<int, std::int32_t> face_table_ids;

    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron()) {
            // Each cell contributes its own rings. Before v9.16.0 polyhedron
            // blocks were skipped here, which meant smooth MOVED their nodes
            // with no inversion guard at all -- silently, since nothing pins
            // them either.
            detail::CellRings rings;
            std::vector<detail::Vec3> coords;
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                if (!detail::cell_rings(cb, c, rMesh.Points(), rMesh.PointDim(), rings, coords))
                    continue;
                bool ok = true;
                const std::size_t first = t.mCornerNodes.size();
                for (std::int64_t id : rings.mNodes) {
                    if (id < 0 || static_cast<std::size_t>(id) >= n) {
                        ok = false;
                        break;
                    }
                    t.mCornerNodes.push_back(id);
                }
                if (!ok) {
                    t.mCornerNodes.resize(first);
                    continue;
                }
                t.mPolyStart.push_back(static_cast<std::int64_t>(t.mPolyFaceStart.size()));
                t.mPolyNumFaces.push_back(static_cast<std::uint32_t>(rings.NumFaces()));
                const std::uint32_t base = static_cast<std::uint32_t>(t.mPolyFaceNodes.size());
                for (std::uint32_t v : rings.mFaceStart)
                    t.mPolyFaceStart.push_back(base + v);
                t.mPolyFaceNodes.insert(t.mPolyFaceNodes.end(), rings.mFaceNodes.begin(),
                                        rings.mFaceNodes.end());
                t.mMeasure.push_back(SmoothMeasure::PolyFan);
                t.mFaceTable.push_back(-1);
                t.mCornerOffset.push_back(static_cast<std::int64_t>(t.mCornerNodes.size()));
            }
            continue;
        }
        if (cb.IsRagged())
            continue;  // 1-level polygons: no enclosed volume to check
        const CellType ct = cell_type_from_name(cb.Type());
        const int corners = detail::cell_corner_count(ct);
        if (corners <= 0)
            continue;
        const int dim = cell_type_dimension(ct);

        SmoothMeasure measure;
        std::int32_t face_id = -1;
        if (dim == 3 && detail::skin_supported(ct)) {
            measure = SmoothMeasure::FaceFan;
            auto it = face_table_ids.find(static_cast<int>(ct));
            if (it == face_table_ids.end()) {
                face_id = static_cast<std::int32_t>(t.mFaceTables.size());
                t.mFaceTables.push_back(&detail::cell_faces(ct));
                face_table_ids.emplace(static_cast<int>(ct), face_id);
            } else {
                face_id = it->second;
            }
        } else if (dim == 2) {
            measure = is2d ? SmoothMeasure::Shoelace2D : SmoothMeasure::NormalFlip;
        } else {
            continue;  // line/vertex blocks: nothing to flip
        }

        const NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t nc = cb.NumCells();
        for (std::size_t c = 0; c < nc; ++c) {
            bool ok = true;
            const std::size_t first = t.mCornerNodes.size();
            for (int k = 0; k < corners; ++k) {
                const std::int64_t id =
                    detail::read_int(conn, c * npc + static_cast<std::size_t>(k));
                if (id < 0 || static_cast<std::size_t>(id) >= n) {
                    ok = false;
                    break;
                }
                t.mCornerNodes.push_back(id);
            }
            if (!ok) {
                t.mCornerNodes.resize(first);  // drop the partially-written row
                continue;
            }
            t.mMeasure.push_back(measure);
            t.mFaceTable.push_back(face_id);
            t.mPolyStart.push_back(-1);
            t.mPolyNumFaces.push_back(0);
            t.mCornerOffset.push_back(static_cast<std::int64_t>(t.mCornerNodes.size()));
        }
    }
    return t;
}

// --- signed measures --------------------------------------------------------

// Corner coordinate lookup with one node optionally substituted, so the guard
// can evaluate "the cell as it would be if node `SubNode` moved to `*pSub`"
// without materialising a modified coordinate buffer in the hot loop.
struct SmoothCornerReader {
    const std::vector<double>* mpXyz;
    const std::int64_t* mpCorners;
    std::int64_t mSubNode;
    const Vec3* mpSub;

    Vec3 operator()(std::size_t Local) const {
        const std::int64_t id = mpCorners[Local];
        if (mpSub != nullptr && id == mSubNode)
            return *mpSub;
        const std::size_t b = static_cast<std::size_t>(id) * 3;
        return {(*mpXyz)[b], (*mpXyz)[b + 1], (*mpXyz)[b + 2]};
    }
};

// Signed volume via a cell-centroid / face-centroid fan over outward-wound
// faces (robust to non-planar faces). Positive for a well-oriented cell.
//
// Term-for-term the same computation as quality.cpp's quality_facefan_volume,
// deliberately kept as a separate flat-buffer variant: the guard evaluates this
// twice per incident cell per free node per pass, and routing it through
// quality's `const std::vector<Vec3>&` signature would put a heap allocation and
// a gather in the hottest loop in this file.
double smooth_facefan_volume(const SmoothCornerReader& rAt, std::size_t NumCorners,
                             const std::vector<detail::CellFaceDef>& rFaces) {
    Vec3 cc = {0.0, 0.0, 0.0};
    for (std::size_t k = 0; k < NumCorners; ++k)
        cc = detail::vec3_add(cc, rAt(k));
    cc = detail::vec3_scale(cc, 1.0 / static_cast<double>(NumCorners));
    double vol = 0.0;
    for (const detail::CellFaceDef& f : rFaces) {
        Vec3 fc = {0.0, 0.0, 0.0};
        for (std::uint8_t k = 0; k < f.mNumCorners; ++k)
            fc = detail::vec3_add(fc, rAt(f.mNodes[k]));
        fc = detail::vec3_scale(fc, 1.0 / static_cast<double>(f.mNumCorners));
        for (std::uint8_t k = 0; k < f.mNumCorners; ++k) {
            const Vec3 a = rAt(f.mNodes[k]);
            const Vec3 b = rAt(f.mNodes[(k + 1) % f.mNumCorners]);
            vol += detail::triple_product(detail::vec3_sub(a, cc), detail::vec3_sub(b, cc),
                                          detail::vec3_sub(fc, cc)) /
                   6.0;
        }
    }
    return vol;
}

// The polyhedral twin of smooth_facefan_volume: identical arithmetic, but the
// faces come from the CELL's own rings (stored per cell in the table) rather
// than from a per-type `cell_faces` table. Same corner-average fan as
// detail/polyhedron.hpp's poly_measure, so the two agree in sign.
double smooth_polyfan_volume(const SmoothCornerReader& rAt, std::size_t NumCorners,
                             const SmoothCellTable& rTable, std::size_t Cell) {
    const std::int64_t start = rTable.mPolyStart[Cell];
    const std::size_t nfaces = rTable.mPolyNumFaces[Cell];
    if (start < 0 || nfaces == 0)
        return 0.0;
    const std::uint32_t* face_start =
        rTable.mPolyFaceStart.data() + static_cast<std::size_t>(start);

    Vec3 cc = {0.0, 0.0, 0.0};
    for (std::size_t k = 0; k < NumCorners; ++k)
        cc = detail::vec3_add(cc, rAt(k));
    cc = detail::vec3_scale(cc, 1.0 / static_cast<double>(NumCorners));
    double vol = 0.0;
    for (std::size_t f = 0; f < nfaces; ++f) {
        const std::uint32_t* ring = rTable.mPolyFaceNodes.data() + face_start[f];
        const std::size_t m = face_start[f + 1] - face_start[f];
        Vec3 fc = {0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < m; ++k)
            fc = detail::vec3_add(fc, rAt(ring[k]));
        fc = detail::vec3_scale(fc, 1.0 / static_cast<double>(m));
        for (std::size_t k = 0; k < m; ++k) {
            const Vec3 a = rAt(ring[k]);
            const Vec3 b = rAt(ring[(k + 1) % m]);
            vol += detail::triple_product(detail::vec3_sub(a, cc), detail::vec3_sub(b, cc),
                                          detail::vec3_sub(fc, cc)) /
                   6.0;
        }
    }
    return vol;
}

// Signed area of a 2D cell. The triangle and quad expressions are copied
// verbatim from quality.cpp's quality_tri / quality_quad (rather than re-derived
// into an algebraically equal form) so that this guard and `compute_quality`
// cannot disagree about the sign of a cell sitting at area ~ 0 -- which would
// let smooth commit a move that compute_quality then reports as inverted.
double smooth_shoelace_area(const SmoothCornerReader& rAt, std::size_t NumCorners) {
    if (NumCorners == 3) {
        const Vec3 p0 = rAt(0);
        const Vec3 p1 = rAt(1);
        const Vec3 p2 = rAt(2);
        return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
    }
    if (NumCorners == 4) {
        const Vec3 p0 = rAt(0);
        const Vec3 p1 = rAt(1);
        const Vec3 p2 = rAt(2);
        const Vec3 p3 = rAt(3);
        return 0.5 * ((p0[0] * p1[1] - p1[0] * p0[1]) + (p1[0] * p2[1] - p2[0] * p1[1]) +
                      (p2[0] * p3[1] - p3[0] * p2[1]) + (p3[0] * p0[1] - p0[0] * p3[1]));
    }
    double a = 0.0;
    for (std::size_t k = 0; k < NumCorners; ++k) {
        const Vec3 p = rAt(k);
        const Vec3 q = rAt((k + 1) % NumCorners);
        a += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * a;
}

// Newell normal of a corner ring (unnormalized).
Vec3 smooth_newell_normal(const SmoothCornerReader& rAt, std::size_t NumCorners) {
    Vec3 nrm = {0.0, 0.0, 0.0};
    for (std::size_t k = 0; k < NumCorners; ++k)
        nrm = detail::vec3_add(nrm, detail::vec3_cross(rAt(k), rAt((k + 1) % NumCorners)));
    return nrm;
}

// Would moving `Node` to `rCand` turn cell `Cell` from valid into inverted?
//
// The rule is strictly "do no harm", not "preserve the sign":
//
//  - A cell that is currently **valid** may not be made invalid -- that is the
//    guard's whole purpose.
//  - A cell that is **already inverted** imposes no constraint at all. Blocking
//    its sign change would pin the guard's own semantics backwards: smoothing
//    is one of the few things that can *repair* a tangled region, and an
//    early version of this function locked in every pre-existing inversion
//    (measured: 5 inverted cells in, 4 still inverted out with the guard on
//    versus 0 with it off) because un-inverting also changes sign.
//  - A cell that is exactly degenerate (measure 0) likewise imposes nothing;
//    treating it as always-flipping would permanently pin every node of a
//    sliver, the opposite of what a smoother is for.
//
// For the volume and 2D-area modes "valid" means positive, which is the same
// convention `compute_quality` uses for its `inverted` metric. The normal-flip
// mode has no absolute convention available -- a facet in 3D has no intrinsic
// orientation -- so there it stays purely relative: did this facet just fold
// back over itself?
bool smooth_cell_flips(std::size_t Cell, const SmoothCellTable& rTable,
                       const std::vector<double>& rXyz, std::int64_t Node, const Vec3& rCand) {
    const std::int64_t off = rTable.mCornerOffset[Cell];
    const std::size_t ncorner = static_cast<std::size_t>(rTable.mCornerOffset[Cell + 1] - off);
    const std::int64_t* corners = rTable.mCornerNodes.data() + off;

    SmoothCornerReader before{&rXyz, corners, Node, nullptr};
    SmoothCornerReader after{&rXyz, corners, Node, &rCand};

    switch (rTable.mMeasure[Cell]) {
        case SmoothMeasure::PolyFan: {
            // Measured WITHOUT re-orienting: a polyhedron's stored winding is
            // arbitrary but fixed for the whole run, so the sign is consistent
            // between `before` and `after` -- which is all "do no harm" needs.
            // Orienting here would repair the winding on every probe and hide
            // exactly the flip the guard exists to catch.
            const double v0 = smooth_polyfan_volume(before, ncorner, rTable, Cell);
            if (v0 == 0.0)
                return false;  // degenerate on arrival: no constraint
            const double v1 = smooth_polyfan_volume(after, ncorner, rTable, Cell);
            return (v0 > 0.0) ? (v1 <= 0.0) : (v1 >= 0.0);
        }
        case SmoothMeasure::FaceFan: {
            const std::vector<detail::CellFaceDef>& faces =
                *rTable.mFaceTables[static_cast<std::size_t>(rTable.mFaceTable[Cell])];
            const double v0 = smooth_facefan_volume(before, ncorner, faces);
            if (v0 <= 0.0)
                return false;  // already inverted/degenerate: no constraint
            return smooth_facefan_volume(after, ncorner, faces) <= 0.0;
        }
        case SmoothMeasure::Shoelace2D: {
            const double a0 = smooth_shoelace_area(before, ncorner);
            if (a0 <= 0.0)
                return false;  // already inverted/degenerate: no constraint
            return smooth_shoelace_area(after, ncorner) <= 0.0;
        }
        case SmoothMeasure::NormalFlip: {
            const Vec3 n0 = smooth_newell_normal(before, ncorner);
            const Vec3 n1 = smooth_newell_normal(after, ncorner);
            if (detail::vec3_norm_sq(n0) == 0.0)
                return false;
            return detail::vec3_dot(n0, n1) <= 0.0;
        }
    }
    return false;
}

// --- CSR built by counting (node -> measurable cell) ------------------------

struct SmoothCsr {
    std::vector<std::int64_t> mXadj;
    std::vector<std::int64_t> mAdj;
};

// Each (node, cell) incidence is emitted exactly once, so no dedup pass is
// needed; a malformed cell repeating a node merely makes the guard re-check that
// cell, which is harmless.
SmoothCsr smooth_build_incidence(const SmoothCellTable& rTable, std::size_t n) {
    SmoothCsr csr;
    csr.mXadj.assign(n + 1, 0);
    const std::size_t nc = rTable.NumCells();
    for (std::size_t c = 0; c < nc; ++c)
        for (std::int64_t k = rTable.mCornerOffset[c]; k < rTable.mCornerOffset[c + 1]; ++k)
            ++csr.mXadj[static_cast<std::size_t>(rTable.mCornerNodes[static_cast<std::size_t>(k)]) +
                        1];
    for (std::size_t i = 0; i < n; ++i)
        csr.mXadj[i + 1] += csr.mXadj[i];
    csr.mAdj.resize(static_cast<std::size_t>(csr.mXadj[n]));
    std::vector<std::int64_t> cursor(csr.mXadj.begin(), csr.mXadj.end() - 1);
    for (std::size_t c = 0; c < nc; ++c) {
        for (std::int64_t k = rTable.mCornerOffset[c]; k < rTable.mCornerOffset[c + 1]; ++k) {
            const std::size_t node =
                static_cast<std::size_t>(rTable.mCornerNodes[static_cast<std::size_t>(k)]);
            csr.mAdj[static_cast<std::size_t>(cursor[node]++)] = static_cast<std::int64_t>(c);
        }
    }
    return csr;
}

// --- boundary + feature detection -------------------------------------------

// Sorted corner ids of one facet, of any arity. detail::FacetKey rather than
// the fixed array<int64_t,4> this used before v9.16.1, for the reason
// surface.cpp switched: a polyhedron's face can have any number of corners, and
// ONE shared key type is what lets a hexahedron and a polyhedron meeting on a
// face cancel each other out instead of both reporting it as boundary.
using SmoothFacetKey = detail::FacetKey;
using SmoothFacetKeyHash = detail::FacetKeyHash;

// One facet of a cell, corners only: unifies CellFaceDef (3D) and CellEdgeDef
// (2D) so the two-phase extractor is dimension-agnostic.
struct SmoothFacetDef {
    std::uint8_t mNumCorners = 0;
    std::array<std::uint8_t, 4> mNodes = {};
};

struct SmoothFacetBlock {
    const NDArray* mpConn = nullptr;
    std::size_t mNpc = 0;
    std::size_t mNumCells = 0;
    std::vector<SmoothFacetDef> mFacets;
    std::size_t mFirstFacet = 0;
    // Polyhedron blocks: faces per cell vary, so record offsets are tabulated.
    bool mPolyhedron = false;
    std::size_t mBlock = 0;
    std::vector<std::size_t> mFaceStart;  // mNumCells + 1
};

struct SmoothFacetRecord {
    SmoothFacetKey mKey;
    std::uint32_t mBlock = 0;
    std::uint32_t mCell = 0;
    std::uint32_t mSlot = 0;
};

// A boundary facet, kept only when feature detection needs its normal.
struct SmoothBoundaryFacet {
    std::array<std::int64_t, 4> mNodes = {-1, -1, -1, -1};
    std::uint8_t mNumCorners = 0;
    Vec3 mNormal = {0.0, 0.0, 0.0};
};

std::vector<SmoothFacetDef> smooth_facets_for(CellType Type, bool FaceMode) {
    std::vector<SmoothFacetDef> out;
    if (FaceMode) {
        for (const detail::CellFaceDef& fd : detail::cell_faces(Type)) {
            SmoothFacetDef d;
            d.mNumCorners = fd.mNumCorners;
            for (std::uint8_t k = 0; k < fd.mNumCorners && k < 4; ++k)
                d.mNodes[k] = fd.mNodes[k];
            out.push_back(d);
        }
    } else {
        for (const detail::CellEdgeDef& ed : detail::cell_edges(Type)) {
            SmoothFacetDef d;
            d.mNumCorners = 2;
            d.mNodes[0] = ed.mNodes[0];
            d.mNodes[1] = ed.mNodes[1];
            out.push_back(d);
        }
    }
    return out;
}

// Marks boundary nodes, and (when rpFacets is non-null) collects the boundary
// facets with their normals for the feature pass.
//
// This is surface.cpp's phase-split idiom re-implemented locally with smooth_
// prefixes, following the v7.6.0 partition precedent: surface.cpp's
// anon-namespace machinery stays untouched. The two serial passes are the
// determinism pin and must never become concurrent hash inserts.
void smooth_mark_boundary(const Mesh& rMesh, std::size_t n, bool FaceMode,
                          const std::vector<double>& rXyz, std::vector<std::uint8_t>& rBoundary,
                          std::vector<SmoothBoundaryFacet>* pFacets) {
    std::vector<SmoothFacetBlock> blocks;
    std::size_t total_facets = 0;
    std::size_t block_index = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t this_block = block_index++;
        if (cb.IsPolyhedron()) {
            if (!FaceMode)
                continue;  // a polyhedron has no 2D boundary edges to mark
            SmoothFacetBlock b;
            b.mPolyhedron = true;
            b.mBlock = this_block;
            b.mNumCells = cb.NumCells();
            b.mFirstFacet = total_facets;
            b.mFaceStart.reserve(b.mNumCells + 1);
            std::size_t at = 0;
            b.mFaceStart.push_back(0);
            for (std::size_t c = 0; c < b.mNumCells; ++c) {
                at += cb.NumFaces(c);
                b.mFaceStart.push_back(at);
            }
            total_facets += at;
            blocks.push_back(std::move(b));
            continue;
        }
        if (cb.IsRagged())
            continue;
        const CellType ct = cell_type_from_name(cb.Type());
        if (cell_type_dimension(ct) != (FaceMode ? 3 : 2))
            continue;
        SmoothFacetBlock b;
        b.mFacets = smooth_facets_for(ct, FaceMode);
        if (b.mFacets.empty())
            continue;
        b.mpConn = &cb.Conn();
        b.mNpc = cb.NodesPerCell();
        b.mNumCells = cb.NumCells();
        b.mFirstFacet = total_facets;
        total_facets += b.mNumCells * b.mFacets.size();
        blocks.push_back(std::move(b));
    }
    if (total_facets == 0)
        return;

    // --- phase 1: build facet keys into disjoint slots (parallel-safe) ---
    std::vector<SmoothFacetRecord> recs(total_facets);
    for (std::uint32_t bi = 0; bi < blocks.size(); ++bi) {
        const SmoothFacetBlock& b = blocks[bi];
        if (b.mPolyhedron) {
            const auto cb = rMesh.Cells(b.mBlock);
            parallel_for(b.mNumCells, [&, bi](std::size_t cell) {
                const std::size_t nf = cb.NumFaces(cell);
                for (std::size_t f = 0; f < nf; ++f) {
                    const auto face = cb.Face(cell, f);
                    SmoothFacetRecord& r = recs[b.mFirstFacet + b.mFaceStart[cell] + f];
                    r.mKey = SmoothFacetKey(face.first, face.second);
                    r.mBlock = bi;
                    r.mCell = static_cast<std::uint32_t>(cell);
                    r.mSlot = static_cast<std::uint32_t>(f);
                }
            });
            continue;
        }
        const NDArray& conn = *b.mpConn;
        const std::size_t npc = b.mNpc;
        const std::size_t fpc = b.mFacets.size();
        parallel_for(b.mNumCells * fpc, [&, bi](std::size_t j) {
            const std::size_t cell = j / fpc;
            const std::size_t slot = j % fpc;
            const SmoothFacetDef& fd = b.mFacets[slot];
            SmoothFacetRecord& r = recs[b.mFirstFacet + j];
            std::array<std::int64_t, 4> ids{};
            const std::uint8_t nk = fd.mNumCorners < 4 ? fd.mNumCorners : 4;
            for (std::uint8_t k = 0; k < nk; ++k)
                ids[k] = detail::read_int(conn, cell * npc + fd.mNodes[k]);
            r.mKey = SmoothFacetKey(ids.data(), nk);
            r.mBlock = bi;
            r.mCell = static_cast<std::uint32_t>(cell);
            r.mSlot = static_cast<std::uint32_t>(slot);
        });
    }

    // --- phase 2, pass A: count key occurrences (serial -> deterministic) ---
    std::unordered_map<SmoothFacetKey, std::uint32_t, SmoothFacetKeyHash> counts;
    counts.reserve(total_facets * 2);
    for (const SmoothFacetRecord& r : recs)
        ++counts[r.mKey];

    // --- phase 2, pass B: mark once-used facets (serial, stored order) ---
    for (const SmoothFacetRecord& r : recs) {
        if (counts[r.mKey] != 1)
            continue;
        const SmoothFacetBlock& b = blocks[r.mBlock];
        if (b.mPolyhedron) {
            const auto cb = rMesh.Cells(b.mBlock);
            const auto face =
                cb.Face(static_cast<std::size_t>(r.mCell), static_cast<std::size_t>(r.mSlot));
            bool ok = true;
            for (std::size_t k = 0; k < face.second; ++k) {
                const std::int64_t id = face.first[k];
                if (id < 0 || static_cast<std::size_t>(id) >= n) {
                    ok = false;
                    break;
                }
                rBoundary[static_cast<std::size_t>(id)] = 1;
            }
            if (!ok || pFacets == nullptr || face.second < 3)
                continue;
            // One normal for the whole face, computed over ALL its corners.
            // SmoothBoundaryFacet holds at most four node ids, so an n-gon is
            // emitted as several records sharing that normal -- every corner
            // then takes part in the feature test, which a single truncated
            // record would silently deny to corners 5+.
            std::vector<std::int64_t> ids(face.first, face.first + face.second);
            const SmoothCornerReader at{&rXyz, ids.data(), -1, nullptr};
            const Vec3 nrm = detail::vec3_normalize(smooth_newell_normal(at, ids.size()));
            for (std::size_t base = 0; base < ids.size(); base += 4) {
                SmoothBoundaryFacet bf;
                bf.mNormal = nrm;
                const std::size_t take = std::min<std::size_t>(4, ids.size() - base);
                bf.mNumCorners = static_cast<std::uint8_t>(take);
                for (std::size_t k = 0; k < take; ++k)
                    bf.mNodes[k] = ids[base + k];
                pFacets->push_back(bf);
            }
            continue;
        }
        const SmoothFacetDef& fd = b.mFacets[r.mSlot];
        const std::size_t row = static_cast<std::size_t>(r.mCell) * b.mNpc;

        SmoothBoundaryFacet bf;
        bf.mNumCorners = fd.mNumCorners;
        bool ok = true;
        for (std::uint8_t k = 0; k < fd.mNumCorners && k < 4; ++k) {
            const std::int64_t id = detail::read_int(*b.mpConn, row + fd.mNodes[k]);
            if (id < 0 || static_cast<std::size_t>(id) >= n) {
                ok = false;
                break;
            }
            bf.mNodes[k] = id;
            rBoundary[static_cast<std::size_t>(id)] = 1;
        }
        if (!ok || pFacets == nullptr)
            continue;

        // Facet normal: Newell for a face, the in-plane perpendicular for a 2D
        // boundary edge (so a polyline's corners read as features too).
        const SmoothCornerReader at{&rXyz, bf.mNodes.data(), -1, nullptr};
        if (fd.mNumCorners >= 3) {
            bf.mNormal = detail::vec3_normalize(smooth_newell_normal(at, fd.mNumCorners));
        } else {
            const Vec3 p0 = at(0);
            const Vec3 p1 = at(1);
            bf.mNormal = detail::vec3_normalize(Vec3{p1[1] - p0[1], p0[0] - p1[0], 0.0});
        }
        pFacets->push_back(bf);
    }
}

// Pin boundary nodes whose incident boundary facets disagree in orientation by
// more than the feature angle. O(d^2) in the boundary valence d, which is 4-8 in
// practice; each iteration writes only its own slot, so it parallelizes cleanly.
void smooth_mark_features(const std::vector<SmoothBoundaryFacet>& rFacets, std::size_t n,
                          double CosThreshold, std::vector<std::uint8_t>& rFrozen) {
    if (rFacets.empty())
        return;
    SmoothCsr inc;
    inc.mXadj.assign(n + 1, 0);
    for (const SmoothBoundaryFacet& f : rFacets)
        for (std::uint8_t k = 0; k < f.mNumCorners && k < 4; ++k)
            ++inc.mXadj[static_cast<std::size_t>(f.mNodes[k]) + 1];
    for (std::size_t i = 0; i < n; ++i)
        inc.mXadj[i + 1] += inc.mXadj[i];
    inc.mAdj.resize(static_cast<std::size_t>(inc.mXadj[n]));
    std::vector<std::int64_t> cursor(inc.mXadj.begin(), inc.mXadj.end() - 1);
    for (std::size_t fi = 0; fi < rFacets.size(); ++fi) {
        const SmoothBoundaryFacet& f = rFacets[fi];
        for (std::uint8_t k = 0; k < f.mNumCorners && k < 4; ++k) {
            const std::size_t node = static_cast<std::size_t>(f.mNodes[k]);
            inc.mAdj[static_cast<std::size_t>(cursor[node]++)] = static_cast<std::int64_t>(fi);
        }
    }

    parallel_for(n, [&](std::size_t i) {
        const std::int64_t b = inc.mXadj[i];
        const std::int64_t e = inc.mXadj[i + 1];
        for (std::int64_t p = b; p < e; ++p) {
            const Vec3& na =
                rFacets[static_cast<std::size_t>(inc.mAdj[static_cast<std::size_t>(p)])].mNormal;
            for (std::int64_t q = p + 1; q < e; ++q) {
                const Vec3& nb =
                    rFacets[static_cast<std::size_t>(inc.mAdj[static_cast<std::size_t>(q)])]
                        .mNormal;
                if (detail::vec3_dot(na, nb) < CosThreshold) {
                    rFrozen[i] = 1;
                    return;
                }
            }
        }
    });
}

// Does the mesh contain any 3D block? Selects face mode vs edge mode, mirroring
// extract_surface's automatic dimension pick.
bool smooth_has_volume_cells(const Mesh& rMesh) {
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            return true;
        if (!cb.IsRagged() && cell_type_dimension(cell_type_from_name(cb.Type())) == 3)
            return true;
    }
    return false;
}

// Pin every node referenced by a block whose edge topology is unknown (the
// higher-order family, VTK-Lagrange, custom). One warn per distinct type.
void smooth_pin_unknown_topology(const Mesh& rMesh, std::size_t n,
                                 std::vector<std::uint8_t>& rFrozen) {
    std::unordered_set<std::string> warned;
    std::vector<std::int64_t> nodes;
    for (const auto cb : rMesh.CellRange()) {
        if (detail::node_edge_topology_known(cb))
            continue;
        const std::string type(cb.Type());
        if (warned.insert(type).second)
            log::warn(
                "smooth: cell type '{}' has no known edge topology; its nodes are pinned "
                "(an unknown neighbourhood gives no defined smoothing target)",
                type);
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            detail::cell_node_ids(cb, c, n, nodes);
            for (std::int64_t id : nodes)
                rFrozen[static_cast<std::size_t>(id)] = 1;
        }
    }
}

// --- ODT: tet-only scope check + circumcenter --------------------------------

// ODT's closed-form vertex update needs a genuine circumsphere per incident
// tet, which only a linear tetra has -- no other 3D cell type (hex/wedge/
// pyramid, any quadratic family, a polyhedron) has one at all. Reject every
// construct outside that scope by name before any node moves; the same shape
// as decimate_volume.cpp's own dv_check_blocks for its identical tet-only
// restriction (a separate, private copy here rather than a shared helper --
// the two operations' error prefixes and follow-up advice genuinely differ).
void smooth_check_odt_blocks(const Mesh& rMesh) {
    bool has_tet = false;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                "meshio++: smooth: method 'odt' operates on tet-only meshes, but the mesh "
                "contains a polyhedron cell block -- run convert_cells(mode='simplexify') first");
        if (cb.IsRagged())
            throw std::invalid_argument(
                "meshio++: smooth: method 'odt' operates on tet-only meshes, but the mesh "
                "contains ragged cell block '" +
                std::string(cb.Type()) + "'");
        const CellType ct = cell_type_from_name(cb.Type());
        if (ct == CellType::Tetra) {
            has_tet = true;
            continue;
        }
        const int dim = cell_type_dimension(ct);
        if (dim == 3)
            throw std::invalid_argument(
                "meshio++: smooth: method 'odt' operates on tet-only meshes, but the mesh "
                "contains 3D cell block '" +
                std::string(cb.Type()) +
                "' that is not linear tetra -- run convert_cells(mode='simplexify') first");
        throw std::invalid_argument(
            "meshio++: smooth: method 'odt' operates on tet-only meshes, but the mesh contains "
            "non-3D cell block '" +
            std::string(cb.Type()) +
            "' alongside its tets (its nodes would dangle after smoothing; drop it first, e.g. "
            "via split)");
    }
    if (!has_tet)
        throw std::invalid_argument(
            "meshio++: smooth: method 'odt' requires at least one tetra cell block");
}

// Circumcenter of tetrahedron (p0, p1, p2, p3): the point equidistant from all
// four corners, found from the three equal-distance planes |x-p1|^2=|x-p0|^2,
// |x-p2|^2=|x-p0|^2, |x-p3|^2=|x-p0|^2, which linearize into A*o=b with rows
// `2*(p_k - p0)` and `b_k = |p_k|^2 - |p0|^2`. Solved by Cramer's rule via
// `detail::det3` (this file's own vocabulary, not `quality.cpp`'s unrelated
// -- and file-private, hence unreachable from here -- `quality_solve3`; the
// two are independently derived from the same textbook construction rather
// than one transcribing the other). Pinned against a regular tetrahedron's
// closed-form circumcenter in `Smooth.OdtCircumcenterMatchesClosedFormOnARegularTetra`.
bool smooth_tet_circumcenter(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3,
                             double eps, Vec3& rOut) {
    const Vec3 a0 = detail::vec3_scale(detail::vec3_sub(p1, p0), 2.0);
    const Vec3 a1 = detail::vec3_scale(detail::vec3_sub(p2, p0), 2.0);
    const Vec3 a2 = detail::vec3_scale(detail::vec3_sub(p3, p0), 2.0);
    const Vec3 rhs = {detail::vec3_norm_sq(p1) - detail::vec3_norm_sq(p0),
                      detail::vec3_norm_sq(p2) - detail::vec3_norm_sq(p0),
                      detail::vec3_norm_sq(p3) - detail::vec3_norm_sq(p0)};
    const double d = detail::det3(a0, a1, a2);
    if (std::fabs(d) < eps)
        return false;
    for (int k = 0; k < 3; ++k) {
        Vec3 c0 = a0, c1 = a1, c2 = a2;
        c0[k] = rhs[0];
        c1[k] = rhs[1];
        c2[k] = rhs[2];
        rOut[k] = detail::det3(c0, c1, c2) / d;
    }
    return true;
}

}  // namespace

SmoothMethod smooth_method_from_name(const std::string& rName) {
    if (rName == "laplacian")
        return SmoothMethod::Laplacian;
    if (rName == "taubin")
        return SmoothMethod::Taubin;
    if (rName == "odt")
        return SmoothMethod::Odt;
    throw std::invalid_argument("meshio++: smooth: unknown method '" + rName +
                                "' (expected 'laplacian', 'taubin' or 'odt')");
}

SmoothResult smooth(const Mesh& rMesh, const SmoothOptions& rOptions) {
    const SmoothParams params = smooth_resolve_params(rOptions);
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    if (!rOptions.mFrozen.empty() && rOptions.mFrozen.size() != n)
        throw std::invalid_argument("meshio++: smooth: frozen mask has " +
                                    std::to_string(rOptions.mFrozen.size()) +
                                    " entries but the mesh has " + std::to_string(n) + " points");

    const bool is_odt = params.mMethod == SmoothMethod::Odt;
    if (is_odt)
        smooth_check_odt_blocks(rMesh);
    // Degenerate-tet / degenerate-circumsphere threshold, matching
    // quality.cpp's own `eps = 1e-14` for the analogous cofactor solve.
    constexpr double kOdtEps = 1e-14;

    // --- phase 0: coordinates as a flat double buffer ---
    const std::vector<double> original = smooth_read_coords(rMesh, n, dim);
    std::vector<double> prev = original;
    std::vector<double> cur(prev.size(), 0.0);

    // --- phase 1: edge adjacency ---
    const detail::NodeAdjacency csr =
        detail::build_node_adjacency(rMesh, n, detail::NodeAdjacencyKind::Edge);

    // --- phase 2: the pin mask (boundary | feature | unknown | caller) ---
    std::vector<std::uint8_t> frozen(n, 0);
    if (!rOptions.mFrozen.empty())
        for (std::size_t i = 0; i < n; ++i)
            frozen[i] = rOptions.mFrozen[i] ? 1 : 0;
    smooth_pin_unknown_topology(rMesh, n, frozen);

    if (rOptions.mFixBoundary) {
        const bool face_mode = smooth_has_volume_cells(rMesh);
        std::vector<std::uint8_t> boundary(n, 0);
        std::vector<SmoothBoundaryFacet> facets;
        smooth_mark_boundary(rMesh, n, face_mode, prev, boundary,
                             rOptions.mPreserveFeatures ? &facets : nullptr);
        for (std::size_t i = 0; i < n; ++i)
            if (boundary[i])
                frozen[i] = 1;
        if (rOptions.mPreserveFeatures) {
            const double cos_thr =
                std::cos(rOptions.mFeatureAngleDeg * 3.14159265358979323846 / 180.0);
            smooth_mark_features(facets, n, cos_thr, frozen);
        }
    } else if (rOptions.mPreserveFeatures) {
        // Features are a subset of the boundary, so asking to preserve them
        // while explicitly freeing the boundary is contradictory rather than
        // merely redundant -- say so instead of silently doing nothing.
        log::warn(
            "smooth: preserve_features has no effect when fix_boundary is off (feature nodes "
            "are boundary nodes)");
    }

    // --- phase 3: the inversion guard's tables ---
    // Also built, unconditionally, when method == Odt: under Odt's tet-only
    // scope this table is exactly the tet corner list, so its node -> cell
    // incidence serves double duty -- ODT's own target computation (phase 4)
    // AND the inversion guard, with no second CSR anywhere. This is the
    // reason ODT was scoped tet-only rather than left general: a general
    // vertex -> incident-cell structure already existed here for the guard,
    // and restricting the scope is what let it be reused rather than
    // duplicated (e.g. via detail/cell_adjacency.hpp's differently-keyed
    // node incidence, which this file has no other use for).
    SmoothCellTable cells;
    SmoothCsr incidence;
    const bool guard = rOptions.mGuardInversion;
    if (guard || is_odt) {
        cells = smooth_build_cell_table(rMesh, n, dim == 2);
        incidence = smooth_build_incidence(cells, n);
    }

    // --- phase 4: the Jacobi iteration ---
    std::vector<std::uint8_t> skipped(n, 0);
    std::int64_t num_skipped = 0;
    for (int pass = 0; pass < params.mNumPasses; ++pass) {
        // Taubin alternates the shrinking (+lambda) and un-shrinking (mu) pass.
        const double factor = (!params.mTaubin || (pass % 2 == 0)) ? params.mLambda : params.mMu;

        std::fill(skipped.begin(), skipped.end(), 0);
        parallel_for(n, [&](std::size_t i) {
            const std::size_t o = i * 3;

            if (frozen[i]) {
                cur[o] = prev[o];
                cur[o + 1] = prev[o + 1];
                cur[o + 2] = prev[o + 2];
                return;
            }

            // The move target: an edge-neighbour mean for Laplacian/Taubin,
            // the volume-weighted circumcenter average of incident tets for
            // Odt. `has_target` is each method's own "isolated node" case
            // (no edge neighbours; no tet contributed a valid circumcenter),
            // which both leave the node exactly where it was.
            Vec3 target = {0.0, 0.0, 0.0};
            bool has_target = false;
            if (is_odt) {
                const std::int64_t ib = incidence.mXadj[i];
                const std::int64_t ie = incidence.mXadj[i + 1];
                // Accumulated over incident tets in ascending cell index
                // (the incidence rows are built by ascending cell, per
                // smooth_build_incidence), which pins the FP order the same
                // way the edge-mean branch pins ascending neighbour id.
                double sum_v = 0.0;
                Vec3 sum_wc = {0.0, 0.0, 0.0};
                for (std::int64_t k = ib; k < ie; ++k) {
                    const std::size_t c =
                        static_cast<std::size_t>(incidence.mAdj[static_cast<std::size_t>(k)]);
                    const std::int64_t cb = cells.mCornerOffset[c];
                    const SmoothCornerReader at{&prev, &cells.mCornerNodes[static_cast<std::size_t>(cb)],
                                                -1, nullptr};
                    const Vec3 p0 = at(0), p1 = at(1), p2 = at(2), p3 = at(3);
                    Vec3 cc;
                    if (!smooth_tet_circumcenter(p0, p1, p2, p3, kOdtEps, cc))
                        continue;  // degenerate tet: no circumsphere, no contribution
                    const Vec3 coords[4] = {p0, p1, p2, p3};
                    const double vol = std::fabs(detail::cell_volume_from_corners(coords, CellType::Tetra));
                    if (vol < kOdtEps)
                        continue;
                    sum_v += vol;
                    sum_wc[0] += vol * cc[0];
                    sum_wc[1] += vol * cc[1];
                    sum_wc[2] += vol * cc[2];
                }
                has_target = sum_v > kOdtEps;
                if (has_target) {
                    const double inv = 1.0 / sum_v;
                    target = {sum_wc[0] * inv, sum_wc[1] * inv, sum_wc[2] * inv};
                }
            } else {
                const std::int64_t b = csr.mXadj[i];
                const std::int64_t e = csr.mXadj[i + 1];
                has_target = b != e;
                if (has_target) {
                    // Summed in ascending neighbour id (the adjacency rows are
                    // sorted), which is what pins the FP accumulation order
                    // across backends and thread counts.
                    Vec3 sum = {0.0, 0.0, 0.0};
                    for (std::int64_t k = b; k < e; ++k) {
                        const std::size_t p =
                            static_cast<std::size_t>(csr.mAdj[static_cast<std::size_t>(k)]) * 3;
                        sum[0] += prev[p];
                        sum[1] += prev[p + 1];
                        sum[2] += prev[p + 2];
                    }
                    const double inv = 1.0 / static_cast<double>(e - b);
                    target = {sum[0] * inv, sum[1] * inv, sum[2] * inv};
                }
            }

            if (!has_target) {
                cur[o] = prev[o];
                cur[o + 1] = prev[o + 1];
                cur[o + 2] = prev[o + 2];
                return;
            }

            // For a 2D Laplacian/Taubin mesh every z is exactly +0.0, so this
            // reduces to 0.0 + factor * (0.0 - 0.0) and the z column stays
            // bit-exactly +0.0 through arbitrarily many passes -- no masking
            // needed. Odt is tet-only, so dim == 3 always and this note does
            // not apply there.
            const Vec3 cand = {
                prev[o] + factor * (target[0] - prev[o]),
                prev[o + 1] + factor * (target[1] - prev[o + 1]),
                prev[o + 2] + factor * (target[2] - prev[o + 2]),
            };

            if (guard) {
                const std::int64_t ib = incidence.mXadj[i];
                const std::int64_t ie = incidence.mXadj[i + 1];
                for (std::int64_t k = ib; k < ie; ++k) {
                    const std::size_t c =
                        static_cast<std::size_t>(incidence.mAdj[static_cast<std::size_t>(k)]);
                    if (smooth_cell_flips(c, cells, prev, static_cast<std::int64_t>(i), cand)) {
                        skipped[i] = 1;
                        cur[o] = prev[o];
                        cur[o + 1] = prev[o + 1];
                        cur[o + 2] = prev[o + 2];
                        return;
                    }
                }
            }
            cur[o] = cand[0];
            cur[o + 1] = cand[1];
            cur[o + 2] = cand[2];
        });

        // Serial fold out of disjoint per-node slots -- never a reduction inside
        // the parallel region (surface.cpp's phase-split philosophy).
        for (std::size_t i = 0; i < n; ++i)
            num_skipped += skipped[i];
        prev.swap(cur);
    }
    if (guard && num_skipped > 0)
        log::warn(
            "smooth: the inversion guard rejected {} node moves; those nodes held still to "
            "keep every incident cell correctly oriented",
            num_skipped);

    // --- phase 5: summary, measured against the input ---
    SmoothResult result;
    result.mNumSkippedInversion = num_skipped;
    if (n > 0) {
        Vec3 lo = {original[0], original[1], original[2]};
        Vec3 hi = lo;
        for (std::size_t i = 1; i < n; ++i)
            for (std::size_t d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], original[i * 3 + d]);
                hi[d] = std::max(hi[d], original[i * 3 + d]);
            }
        const double diag = detail::vec3_norm(detail::vec3_sub(hi, lo));
        const double tol = rOptions.mMoveTolerance * (diag > 0.0 ? diag : 1.0);
        for (std::size_t i = 0; i < n; ++i) {
            const Vec3 d = {prev[i * 3] - original[i * 3], prev[i * 3 + 1] - original[i * 3 + 1],
                            prev[i * 3 + 2] - original[i * 3 + 2]};
            const double len = detail::vec3_norm(d);
            if (len > tol)
                ++result.mNumNodesMoved;
            result.mMaxDisplacement = std::max(result.mMaxDisplacement, len);
        }
    }

    // --- phase 6: write-back (only the coordinates change) ---
    Mesh& out = result.mMesh;
    out.AssignPoints(smooth_write_coords(rMesh.Points(), prev, n, dim));

    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron()) {
            std::vector<std::vector<std::vector<std::int64_t>>> blocks(cb.NumCells());
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                blocks[c].resize(cb.NumFaces(c));
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    auto face = cb.Face(c, f);
                    blocks[c][f].assign(face.first, face.first + face.second);
                }
            }
            out.AddPolyhedronBlock(std::string(cb.Type()), std::move(blocks));
        } else if (cb.IsRagged()) {
            std::vector<std::vector<std::int64_t>> rows(cb.NumCells());
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                rows[c].assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
            out.AddPolygonBlock(std::string(cb.Type()), std::move(rows));
        } else {
            out.AddCellBlock(std::string(cb.Type()), smooth_owned_copy(cb.Conn()));
        }
    }

    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, smooth_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.CellDataNames()) {
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b)
            blocks.push_back(smooth_owned_copy(rMesh.CellData(name, b)));
        out.AddCellData(name, std::move(blocks));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, smooth_owned_copy(rMesh.FieldData(name)));

    // Named regions pass through verbatim: smoothing is a pure coordinate move,
    // so no point, cell or facet is renumbered.
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        out.AddRegion(rMesh.Region(i));
    // Property sets ride along too: these operations preserve the mesh's shape,
    // so dropping a deck's material data here would be new lossiness. They are
    // keyed by id, not by entity index, so there is nothing to remap.
    for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i)
        out.AddPropertySet(rMesh.GetPropertySet(i));

    return result;
}

}  // namespace meshioplusplus
