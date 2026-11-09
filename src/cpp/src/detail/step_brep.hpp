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
#pragma once

/**
 * @file detail/step_brep.hpp
 * @brief The kernel-free STEP B-rep writer behind `StepFaces::Planar` and
 * `StepFaces::Fitted`: a closed skin as one face per surface instead of one
 * per facet, with no CAD kernel in the loop.
 *
 *  - **Faces.** With fitted patches (`segment:patch` and the `fit:patch_*`
 *    table of `fit_primitives`), every accepted plane, sphere, cylinder, cone
 *    or torus patch whose points all lie within the tolerance of its surface
 *    is one `ADVANCED_FACE` on a `PLANE`, `SPHERICAL_SURFACE`,
 *    `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE` or `TOROIDAL_SURFACE`. Every
 *    other facet joins an edge-connected coplanar group (largest facet
 *    first, every corner within `1e-6` of the diagonal of the seed's plane),
 *    one `PLANE` face each, holes included; a group whose points all lie
 *    within the tolerance of an analytic neighbour's surface joins that face.
 *  - **Edges.** The mesh edges between two faces form chains from corner to
 *    corner (a point three faces meet at, or where a boundary branches). A
 *    chain between two planes is split at its corners (Douglas-Peucker) into
 *    `LINE`s; one beside a curved analytic face is a `CIRCLE` when its points
 *    lie within the tolerance of a circle about a shared axis (a bore's rim),
 *    a `LINE` when they lie on a segment, and otherwise a degree-1
 *    `B_SPLINE_CURVE_WITH_KNOTS` through them; one between a curved face and
 *    a coplanar group keeps the mesh's straight edges. A closed chain with no
 *    corner is one closed edge.
 *  - **Vertices.** A point among analytic faces only is moved onto all their
 *    surfaces (alternating projections), when that moves it no further than
 *    the tolerance; any other keeps the mesh position.
 *  - **Loops.** Each face's boundary is walked through the facet fan around
 *    every point, so a loop touching itself stays two loops; a planar face's
 *    loop of largest area is its outer bound. A periodic face has no seam
 *    edge: STEP allows its bounds without one, and CAD readers add it.
 *
 * A patch whose loops do not close is written as coplanar groups instead,
 * and counted demoted. Deterministic and serial.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/mesh.hpp"

// Project includes (private, not installed)
#include "bspline_curve.hpp"
#include "bspline_lsq.hpp"
#include "closed_facets.hpp"
#include "step_stream.hpp"

namespace meshioplusplus {
namespace detail {

/// A freeform B-spline face for a patch (v17.12.0, `fit_cad`'s native
/// backend outside the legacy mode): written as a `B_SPLINE_SURFACE_WITH_KNOTS`
/// face, its edges polylines through the skin's points, which the fit holds
/// within the patch's allowance.
struct StepFreeformFace {
    /// The `segment:patch` label; its row of the fit table becomes the face.
    std::int64_t mLabel = -1;
    BsplineSurface mSurface;
    /// Whether the surface's normal (`S_u x S_v`) agrees with the facets'.
    bool mSense = true;
    /// A patch of a network: its boundary rows are pieces of `StepEdgeCurve`s.
    bool mNetwork = false;
};

/// A patch network's shared edge curve (v17.12.0): the mesh points of the
/// chain it stands for, their parameters on it, and the curve. An edge
/// between two network faces that runs along one of these is written as the
/// exact piece of it both faces' boundaries are.
struct StepEdgeCurve {
    std::vector<std::int64_t> mVerts;
    std::vector<double> mParams;
    BsplineCurve mCurve;
};

/// How the B-rep writer builds its faces.
struct StepBrepOptions {
    /// Use the fitted patches of `rSurface` (else coplanar groups only).
    bool mFitted = false;
    /// The largest distance, as a length, a point may lie from its face's
    /// surface. At most 0 takes three times the largest accepted
    /// `fit:patch_tolerance`, else `1e-3` of the bounding-box diagonal.
    double mTolerance = 0.0;
    /// Freeform faces for patches the fit table leaves unexplained (only with
    /// `mFitted`); null or empty writes none, as before the option existed.
    const std::vector<StepFreeformFace>* mFreeform = nullptr;
    /// The networks' shared edge curves; null writes polyline edges, as before.
    const std::vector<StepEdgeCurve>* mEdgeCurves = nullptr;
    /// Also build the solid as trimmed NURBS faces (`StepBrepResult::mBrep`,
    /// `detail/brep_nurbs`), the file unchanged.
    bool mBrep = false;
    /// The length unit the file declares, metres per unit; at most 0 is
    /// millimetres (`StepStream`).
    double mLengthUnit = -1.0;
    /// More items written into the same file after the solids (roadmap 7.6.5,
    /// `write_step_elements`' balls and bounded planes): given the stream, it
    /// appends `MANIFOLD_SOLID_BREP` ids to the first list and
    /// `SHELL_BASED_SURFACE_MODEL` ids to the second. Empty writes none, as before.
    std::function<void(StepStream&, std::vector<std::int64_t>&, std::vector<std::int64_t>&)>
        mExtraItems;
};

/// What the B-rep writer wrote.
struct StepBrepResult {
    std::vector<std::int64_t> mFaceOfFacet;   ///< per facet: its face (writing order)
    std::vector<std::int64_t> mTypeOfFacet;   ///< per facet: 0 a coplanar group, else 1..5, 6
    std::vector<std::int64_t> mPatchOfFacet;  ///< per facet: its `segment:patch`, or -1
    double mTolerance = 0.0;
    double mMaxDeviation = 0.0;  ///< the largest distance of a point from its face's surface
    double mVolume = 0.0;        ///< the skin's enclosed volume
    bool mValid = false;         ///< every face's loops closed, every edge used twice
    /// the faces by type, the fit table's rows and the accepted patches written
    /// as coplanar groups (`mNumDemoted`)
    BrepFaceCounts mCounts;
    std::int64_t mNumSolids = 0;
    std::int64_t mNumLineEdges = 0;
    std::int64_t mNumCircleEdges = 0;
    std::int64_t mNumSplineEdges = 0;
    /// Of those, the networks' shared cubic edge curves.
    std::int64_t mNumSharedEdges = 0;
    /// Of those, the ones between two network faces written as a
    /// SURFACE_CURVE with a p-curve on each.
    std::int64_t mNumPcurveEdges = 0;
    /// Faces written from `StepBrepOptions::mFreeform` (`fit_cad:type` 6).
    std::int64_t mNumBspline = 0;
    /// Per row of the fit table (its order): the `segment:patch` label, and
    /// why an accepted patch was written as coplanar groups (a
    /// `FitCadDemotionReason` value: 2 a point beyond the tolerance, 3 loops
    /// that do not close or no boundary at all), else 0.
    std::vector<std::int64_t> mRowLabel;
    std::vector<std::int32_t> mRowReason;
    /// With `StepBrepOptions::mBrep`: the solid as NURBS, its faces in the
    /// file's order (id `i + 1`), complete or empty with `mBrepWhy`.
    BrepInfo mBrep;
    std::string mBrepWhy;
};

/**
 * @brief Write @p rShell (the facets of @p rSurface) as a STEP B-rep.
 * @param rSurface the surface the facets came from; its `segment:patch` and
 *        `fit:patch_*` data are read when `mFitted`.
 * @throws WriteError on a skin whose faces cannot be closed into loops, or
 *         when the file cannot be written.
 */
StepBrepResult step_brep_write(const std::string& rPath, const Mesh& rSurface,
                               const ClosedFacets& rShell, const StepBrepOptions& rOptions);

/// `write_step` with options, also handing back the B-rep writer's report
/// (left as constructed for `StepFaces::Facet`). Implemented in
/// `formats/step.cpp`; `fit_cad`'s native backend reads the report.
void write_step_report(const std::string& rPath, const Mesh& rMesh,
                       const StepWriteOptions& rOptions, StepBrepResult* pReport,
                       const std::vector<StepFreeformFace>* pFreeform = nullptr,
                       const std::vector<StepEdgeCurve>* pEdgeCurves = nullptr,
                       bool BuildBrep = false);

/**
 * @brief Label each facet with its edge-connected piece and sum each piece's
 * signed volume (positive when wound outwards).
 * @param nf Number of facets.
 * @param FindRoot The union-find root of a facet, after the facets sharing
 *        an edge have been united.
 * @param rStart, rNodes The facets as CSR rows into @p rXyz (x, y, z each).
 * @param rPieceOf Out: each facet's piece, numbered by its root's order.
 * @param rVolume Out: each piece's signed volume, summed facet by facet.
 * @return The number of pieces.
 */
template <class Find>
std::size_t step_piece_volumes(std::size_t nf, Find&& FindRoot,
                               const std::vector<std::int64_t>& rStart,
                               const std::vector<std::int64_t>& rNodes,
                               const std::vector<double>& rXyz, std::vector<std::size_t>& rPieceOf,
                               std::vector<double>& rVolume) {
    rPieceOf.assign(nf, 0);
    std::vector<std::size_t> roots;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t r = FindRoot(f);
        if (r == f)
            roots.push_back(f);
        rPieceOf[f] = static_cast<std::size_t>(std::lower_bound(roots.begin(), roots.end(), r) -
                                               roots.begin());
    }
    const std::size_t npieces = roots.size();
    rVolume.assign(npieces, 0.0);
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t b = rStart[f], s = rStart[f + 1] - b;
        const double* p0 = &rXyz[static_cast<std::size_t>(rNodes[static_cast<std::size_t>(b)]) * 3];
        for (std::int64_t k = 1; k + 1 < s; ++k) {
            const double* p1 =
                &rXyz[static_cast<std::size_t>(rNodes[static_cast<std::size_t>(b + k)]) * 3];
            const double* p2 =
                &rXyz[static_cast<std::size_t>(rNodes[static_cast<std::size_t>(b + k + 1)]) * 3];
            rVolume[rPieceOf[f]] +=
                (p0[0] * (p1[1] * p2[2] - p1[2] * p2[1]) - p0[1] * (p1[0] * p2[2] - p1[2] * p2[0]) +
                 p0[2] * (p1[0] * p2[1] - p1[1] * p2[0])) /
                6.0;
        }
    }
    return npieces;
}

/// A hash of an undirected mesh edge `(lo, hi)`. Buckets only: the maps it
/// serves are looked up, never iterated, so it cannot reach any output.
struct StepEdgeHash {
    std::size_t operator()(const std::pair<std::int64_t, std::int64_t>& rKey) const {
        const std::uint64_t lo = static_cast<std::uint64_t>(rKey.first);
        const std::uint64_t hi = static_cast<std::uint64_t>(rKey.second);
        return static_cast<std::size_t>(((hi << 32) ^ lo ^ (hi >> 32)) * 0x9E3779B97F4A7C15ull);
    }
};

/// The undirected key of the mesh edge `U`-`V`.
inline std::pair<std::int64_t, std::int64_t> step_edge(std::int64_t U, std::int64_t V) {
    return {std::min(U, V), std::max(U, V)};
}

/**
 * @brief The chains of boundary edges between face groups, as the kernel-free
 * writer and the OpenCASCADE backend both cut them: from every junction (in
 * ascending point order) along the unvisited edges, each chain ending at the
 * next junction, then the closed loops no junction lies on (from their lowest
 * edge).
 *
 * @param NumPoints the skin's point count.
 * @param rEdges the boundary edges, distinct `(lo, hi)` pairs, sorted ascending.
 * @param rJunction `rJunction(V, NumNeighbours)`: whether point `V`, with that
 *        many boundary neighbours, ends a chain.
 * @param rChain `rChain(std::vector<std::int64_t> Verts)`, called once per
 *        chain in walking order; a loop's `Verts` ends where it starts.
 */
template <class TJunction, class TChain>
void step_boundary_chains(std::size_t NumPoints,
                          const std::vector<std::pair<std::int64_t, std::int64_t>>& rEdges,
                          TJunction&& rJunction, TChain&& rChain) {
    std::vector<std::vector<std::int64_t>> nbr(NumPoints);
    for (const auto& e : rEdges) {
        nbr[static_cast<std::size_t>(e.first)].push_back(e.second);
        nbr[static_cast<std::size_t>(e.second)].push_back(e.first);
    }
    for (auto& n : nbr)
        std::sort(n.begin(), n.end());
    auto junction = [&](std::int64_t V) {
        return rJunction(V, nbr[static_cast<std::size_t>(V)].size());
    };
    std::unordered_map<std::pair<std::int64_t, std::int64_t>, std::uint8_t, StepEdgeHash> taken;
    taken.reserve(rEdges.size());
    for (const auto& e : rEdges)
        taken.emplace(e, 0);
    auto take = [&](std::int64_t U, std::int64_t V) {
        auto it = taken.find(step_edge(U, V));
        if (it == taken.end() || it->second)
            return false;
        it->second = 1;
        return true;
    };
    auto walk = [&](std::int64_t Start, std::int64_t First, bool Loop) {
        std::vector<std::int64_t> verts{Start, First};
        std::int64_t prev = Start, cur = First;
        while (!(Loop ? cur == Start : junction(cur))) {
            std::int64_t next = -1;
            for (std::int64_t n : nbr[static_cast<std::size_t>(cur)])
                if (n != prev && take(cur, n)) {
                    next = n;
                    break;
                }
            if (next < 0) {
                // Back to the start of a loop through its last edge.
                if (Loop && take(cur, Start))
                    verts.push_back(Start);
                break;
            }
            verts.push_back(next);
            prev = cur;
            cur = next;
        }
        return verts;
    };
    for (std::int64_t v = 0; v < static_cast<std::int64_t>(NumPoints); ++v) {
        if (nbr[static_cast<std::size_t>(v)].empty() || !junction(v))
            continue;
        for (std::int64_t n : nbr[static_cast<std::size_t>(v)])
            if (take(v, n))
                rChain(walk(v, n, false));
    }
    for (const auto& e : rEdges)
        if (take(e.first, e.second))
            rChain(walk(e.first, e.second, true));
}

/**
 * @brief The patch networks' shared curves (`StepEdgeCurve`), indexed by the
 * mesh edges they run along: which curve, and where on it, a chain of mesh
 * edges follows.
 */
class StepCurveIndex {
public:
    /// The index of @p pCurves (none when null).
    explicit StepCurveIndex(const std::vector<StepEdgeCurve>* pCurves) : mpCurves(pCurves) {
        if (pCurves == nullptr)
            return;
        for (std::size_t k = 0; k < pCurves->size(); ++k) {
            const StepEdgeCurve& ec = (*pCurves)[k];
            for (std::size_t i = 0; i + 1 < ec.mVerts.size(); ++i)
                mAt.emplace(step_edge(ec.mVerts[i], ec.mVerts[i + 1]), std::make_pair(k, i));
        }
    }

    /// Whether no curve runs along any mesh edge.
    bool Empty() const { return mAt.empty(); }

    /**
     * @brief Whether every edge of the chain @p rVs lies on one curve, with the
     * parameters strictly monotone along the chain and its two ends apart;
     * then @p rCurve is that curve and `[rT0, rT1]` the chain's piece of it, in
     * the chain's order.
     */
    bool ChainOnCurve(const std::vector<std::int64_t>& rVs, std::size_t& rCurve, double& rT0,
                      double& rT1) const {
        bool one = true;
        std::size_t curve = 0;
        std::vector<double> ts;
        for (std::size_t i = 0; i + 1 < rVs.size() && one; ++i) {
            const auto it = mAt.find(step_edge(rVs[i], rVs[i + 1]));
            if (it == mAt.end() || (i > 0 && it->second.first != curve)) {
                one = false;
                break;
            }
            curve = it->second.first;
            const StepEdgeCurve& ec = (*mpCurves)[curve];
            const std::size_t at = it->second.second;
            // This mesh edge's two ends' parameters, in the chain's order.
            const bool fwd = ec.mVerts[at] == rVs[i];
            if (i == 0)
                ts.push_back(ec.mParams[fwd ? at : at + 1]);
            ts.push_back(ec.mParams[fwd ? at + 1 : at]);
        }
        for (std::size_t i = 2; one && i < ts.size(); ++i)
            one = (ts[i] - ts[i - 1]) * (ts[1] - ts[0]) > 0.0;
        if (!(one && ts.size() >= 2 && ts.front() != ts.back()))
            return false;
        rCurve = curve;
        rT0 = ts.front();
        rT1 = ts.back();
        return true;
    }

private:
    const std::vector<StepEdgeCurve>* mpCurves;
    std::unordered_map<std::pair<std::int64_t, std::int64_t>, std::pair<std::size_t, std::size_t>,
                       StepEdgeHash>
        mAt;
};

/// The point after @p V in facet @p F's ring (CSR rows @p rStart, @p rNodes), or
/// -1 when @p V is not on it.
inline std::int64_t step_ring_next(const std::vector<std::int64_t>& rStart,
                                   const std::vector<std::int64_t>& rNodes, std::size_t F,
                                   std::int64_t V) {
    const std::int64_t b = rStart[F], n = rStart[F + 1] - b;
    for (std::int64_t k = 0; k < n; ++k)
        if (rNodes[static_cast<std::size_t>(b + k)] == V)
            return rNodes[static_cast<std::size_t>(b + (k + 1) % n)];
    return -1;
}

/**
 * @brief Around point @p V, after the boundary edge @p X -> @p V of group @p G
 * (walked that way by a facet of G): the point W of G's next boundary edge
 * V -> W, found by turning through G's facets around V; -1 when the fan breaks.
 *
 * @param rWalker `rWalker(U, V)`: the facet walking the edge U -> V, or -1.
 * @param rGroupOf `rGroupOf(F)`: facet F's group.
 * @param nf the facet count, which bounds the turn.
 */
template <class TWalker, class TGroupOf>
std::int64_t step_fan_next(const std::vector<std::int64_t>& rStart,
                           const std::vector<std::int64_t>& rNodes, std::size_t nf,
                           TWalker&& rWalker, TGroupOf&& rGroupOf, std::int64_t G, std::int64_t X,
                           std::int64_t V) {
    std::int64_t f = rWalker(X, V);
    for (std::size_t guard = 0; f >= 0 && guard <= nf; ++guard) {
        const std::int64_t w = step_ring_next(rStart, rNodes, static_cast<std::size_t>(f), V);
        if (w < 0)
            return -1;
        const std::int64_t across = rWalker(w, V);
        if (across < 0 || rGroupOf(static_cast<std::size_t>(across)) != G)
            return w;
        f = across;
    }
    return -1;
}

}  // namespace detail
}  // namespace meshioplusplus
