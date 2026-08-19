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
//  Attribution:     The isotropic clustering algorithm is derived from
//                   pyacvd (https://github.com/pyvista/pyacvd), MIT License,
//                   Copyright (c) 2017-2024 The PyVista Developers, which is
//                   itself an independent implementation of the research of
//                   S. Valette and J.-M. Chassery. RemeshMetric::Quadric has
//                   no pyacvd counterpart and is a fresh synthesis built on
//                   meshio++'s own pre-existing quadric machinery
//                   (detail/decimate_common.hpp). See doc/remesh.md.
//
// Built entirely through the uniform mesh API so it compiles under every mesh
// backend. See operations/remesh.hpp for the contract.
//
// Determinism: seeding is RNG-free (fronts grow from the lowest unassigned
// vertex), the energy sweep visits the edge list in a fixed order and is
// SERIAL -- every accepted move mutates the accumulators the next test reads,
// so the objective is inherently sequential in exactly the way decimate's
// greedy queue is -- and the setup passes fill disjoint slots with a fixed FP
// order. Output is byte-identical across mesh backends and thread counts.
//
// There is deliberately NO numpy twin: each sweep branches on a comparison
// between two candidate badness values and an accepted move mutates the
// accumulators the next comparison reads, so a single near-tie decided
// differently yields a different clustering rather than a last-ulp difference
// (subdivide / agglomerate / decimate_volume / conservative_interpolate
// precedent).

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
#include "meshioplusplus/operations/remesh.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/decimate_common.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// Anonymous-namespace helpers are uniquely prefixed `remesh_` because the
// amalgamation concatenates every operations TU into one translation unit.

/// The triangulated working surface: flat corner ids plus Float64 coordinates.
struct RemeshSurface {
    std::vector<std::int64_t> mCorners;  ///< `3 * mNumFaces` corner ids.
    std::vector<double> mXyz;            ///< `3 * mNumPoints` coordinates.
    std::size_t mNumFaces = 0;
    std::size_t mNumPoints = 0;
};

/// Rejects everything outside the surface scope, by name. Mirrors
/// `decim_check_blocks` -- the two operations have deliberately identical
/// scope, so their diagnostics should read the same way.
void remesh_check_blocks(const Mesh& rMesh) {
    int max_dim = -1;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                "meshio++: remesh: cannot remesh polyhedron cell block '" + cb.Type() +
                "'; this operation is surface-only (use extract_surface first)");
        const CellType ct = cell_type_from_name(cb.Type());
        max_dim = std::max(max_dim, cell_type_dimension(ct));
    }
    if (max_dim == 3)
        throw std::invalid_argument(
            "meshio++: remesh: mesh contains 3D volume cells; this operation is surface-only "
            "(use extract_surface first, or decimate_volume for a tetrahedral mesh)");
    if (max_dim < 2)
        throw std::invalid_argument("meshio++: remesh: mesh contains no surface (2D) cell block");

    for (const auto cb : rMesh.CellRange()) {
        const CellType ct = cell_type_from_name(cb.Type());
        if (cell_type_dimension(ct) != 2)
            throw std::invalid_argument(
                "meshio++: remesh: cannot remesh a mesh mixing cell block '" + cb.Type() +
                "' with its surface blocks; split the mesh first");
        if (cb.IsRagged())
            throw std::invalid_argument("meshio++: remesh: cannot remesh ragged cell block '" +
                                        cb.Type() + "'");
        if (ct != CellType::Triangle && ct != CellType::Quad && ct != CellType::Polygon)
            throw std::invalid_argument(
                "meshio++: remesh: cannot remesh higher-order cell block '" + cb.Type() +
                "'; call convert_cells(Linearize) first");
    }
}

/// Flattens an all-triangle mesh into the working surface.
RemeshSurface remesh_build_surface(const Mesh& rMesh) {
    RemeshSurface s;
    s.mNumPoints = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    s.mXyz.assign(s.mNumPoints * 3, 0.0);
    const NDArray& points = rMesh.Points();
    for (std::size_t i = 0; i < s.mNumPoints; ++i)
        for (std::size_t d = 0; d < dim && d < 3; ++d)
            s.mXyz[i * 3 + d] = detail::read_double(points, i * dim + d);

    for (const auto cb : rMesh.CellRange()) {
        if (cell_type_from_name(cb.Type()) != CellType::Triangle)
            throw std::invalid_argument(
                "meshio++: remesh: internal error: simplexified block '" + cb.Type() +
                "' is not a triangle block");
        const std::size_t nc = cb.NumCells();
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < nc; ++c)
            for (std::size_t k = 0; k < 3; ++k) {
                const std::int64_t v = detail::read_int(conn, c * 3 + k);
                if (v < 0 || static_cast<std::size_t>(v) >= s.mNumPoints)
                    throw std::invalid_argument(
                        "meshio++: remesh: connectivity references point " + std::to_string(v) +
                        ", which is outside the mesh's " + std::to_string(s.mNumPoints) +
                        " points");
                s.mCorners.push_back(v);
            }
        s.mNumFaces += nc;
    }
    return s;
}

/// Per-item weight (a third of the incident triangle area, so the weights sum
/// to the surface area) and the weighted position each cluster accumulates.
void remesh_item_weights(const RemeshSurface& rS, std::vector<double>& rArea,
                         std::vector<double>& rWeighted) {
    rArea.assign(rS.mNumPoints, 0.0);
    rWeighted.assign(rS.mNumPoints * 3, 0.0);

    for (std::size_t f = 0; f < rS.mNumFaces; ++f) {
        const std::int64_t i0 = rS.mCorners[f * 3];
        const std::int64_t i1 = rS.mCorners[f * 3 + 1];
        const std::int64_t i2 = rS.mCorners[f * 3 + 2];
        const double* v0 = rS.mXyz.data() + static_cast<std::size_t>(i0) * 3;
        const double* v1 = rS.mXyz.data() + static_cast<std::size_t>(i1) * 3;
        const double* v2 = rS.mXyz.data() + static_cast<std::size_t>(i2) * 3;

        const double e0[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        const double e1[3] = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
        const double c[3] = {e0[1] * e1[2] - e0[2] * e1[1], e0[2] * e1[0] - e0[0] * e1[2],
                             e0[0] * e1[1] - e0[1] * e1[0]};
        const double area = 0.5 * std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) / 3.0;

        rArea[static_cast<std::size_t>(i0)] += area;
        rArea[static_cast<std::size_t>(i1)] += area;
        rArea[static_cast<std::size_t>(i2)] += area;
    }
    for (std::size_t i = 0; i < rS.mNumPoints; ++i)
        for (int d = 0; d < 3; ++d)
            rWeighted[i * 3 + d] = rArea[i] * rS.mXyz[i * 3 + d];
}

/// Each vertex's own accumulated quadric (sum of its incident face planes),
/// via `detail/decimate_common.hpp`'s existing, already-tested machinery.
/// Empty (never computed) when `Metric == Isotropic`, which never reads it.
std::vector<double> remesh_vertex_quadrics(const RemeshSurface& rS, RemeshMetric Metric) {
    if (Metric != RemeshMetric::Quadric)
        return {};
    detail::DecimFaces faces;
    faces.mCorners = rS.mCorners;
    faces.mBlockBase = {0};
    faces.mNumFaces = rS.mNumFaces;

    const detail::DecimCsr csr = detail::decim_vertex_faces_csr(faces, rS.mNumPoints);
    std::vector<double> quad_k;
    std::vector<double> normals;
    detail::decim_face_planes(faces, rS.mXyz, quad_k, normals);
    return detail::decim_accumulate_quadrics(csr, rS.mNumPoints, quad_k);
}

/// The undirected edge list of the adjacency graph, each edge once with
/// `lo < hi`, in ascending `(lo, hi)` order -- the sweep order, hence part of
/// the determinism contract.
std::vector<std::int64_t> remesh_unique_edges(const detail::NodeAdjacency& rAdj,
                                              std::size_t NumPoints) {
    std::vector<std::int64_t> edges;
    for (std::size_t i = 0; i < NumPoints; ++i)
        for (std::int64_t k = rAdj.mXadj[i]; k < rAdj.mXadj[i + 1]; ++k) {
            const std::int64_t j = rAdj.mAdj[static_cast<std::size_t>(k)];
            if (j > static_cast<std::int64_t>(i)) {
                edges.push_back(static_cast<std::int64_t>(i));
                edges.push_back(j);
            }
        }
    return edges;
}

/// Grows `NumClusters` seeds breadth-first from the lowest unassigned vertex,
/// each absorbing neighbours until it holds its share of the remaining area.
/// RNG-free by construction, unlike pyacvd.
void remesh_seed_clusters(std::vector<std::int64_t>& rClusters, const detail::NodeAdjacency& rAdj,
                          const std::vector<double>& rArea, std::int64_t NumClusters,
                          std::size_t NumPoints) {
    double area_remain = 0.0;
    for (std::size_t i = 0; i < NumPoints; ++i)
        area_remain += rArea[i];

    const double target = area_remain / static_cast<double>(NumClusters);
    std::size_t next_free = 0;
    std::vector<std::int64_t> ring;
    std::vector<std::int64_t> next;

    for (std::int64_t c = 0; c < NumClusters; ++c) {
        const double t_area = area_remain - target * static_cast<double>(NumClusters - c - 1);
        double c_area = 0.0;

        while (next_free < NumPoints && rClusters[next_free] != -1)
            ++next_free;
        if (next_free >= NumPoints)
            break;

        c_area += rArea[next_free];
        rClusters[next_free] = c;
        ring.clear();
        ring.push_back(static_cast<std::int64_t>(next_free));

        while (!ring.empty()) {
            next.clear();
            for (const std::int64_t p : ring)
                for (std::int64_t k = rAdj.mXadj[static_cast<std::size_t>(p)];
                     k < rAdj.mXadj[static_cast<std::size_t>(p) + 1]; ++k) {
                    const std::int64_t nbr = rAdj.mAdj[static_cast<std::size_t>(k)];
                    const std::size_t un = static_cast<std::size_t>(nbr);
                    if (rClusters[un] == -1 && rArea[un] + c_area < t_area) {
                        c_area += rArea[un];
                        rClusters[un] = c;
                        next.push_back(nbr);
                    }
                }
            ring.swap(next);
        }
        area_remain -= c_area;
    }
}

/// Propagates cluster ids across edges into vertices the seeding never
/// reached. Returns true when some vertex is still unassigned afterwards.
bool remesh_grow_null(const std::vector<std::int64_t>& rEdges,
                      std::vector<std::int64_t>& rClusters) {
    const std::size_t n_edges = rEdges.size() / 2;
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::size_t e = 0; e < n_edges; ++e) {
            const std::size_t a = static_cast<std::size_t>(rEdges[2 * e]);
            const std::size_t b = static_cast<std::size_t>(rEdges[2 * e + 1]);
            if (rClusters[a] == -1 && rClusters[b] != -1) {
                rClusters[a] = rClusters[b];
                changed = true;
            } else if (rClusters[b] == -1 && rClusters[a] != -1) {
                rClusters[b] = rClusters[a];
                changed = true;
            }
        }
    }
    for (const std::int64_t c : rClusters)
        if (c == -1)
            return true;
    return false;
}

/// Isotropic centroid of a cluster from its raw accumulators -- the shared
/// ill-conditioning fallback the quadric solve always has available.
inline void remesh_centroid_of(const double* pSgamma, double Srho, double pOut[3]) {
    if (Srho > 0.0) {
        pOut[0] = pSgamma[0] / Srho;
        pOut[1] = pSgamma[1] / Srho;
        pOut[2] = pSgamma[2] / Srho;
    } else {
        pOut[0] = pOut[1] = pOut[2] = 0.0;
    }
}

/// Quadric error of a cluster whose combined quadric is `pQ`, evaluated at
/// its own optimal point (`pCentroidFallback` used when that solve is
/// ill-conditioned -- an exactly planar cluster, the expected case).
inline double remesh_quadric_badness(const double* pQ, const double* pCentroidFallback) {
    double pt[3];
    detail::decim_quadric_optimal_point(pQ, pCentroidFallback, pt);
    return detail::decim_quadric_error(pQ, pt[0], pt[1], pt[2]);
}

/// The "badness to minimize" of a cluster GAINING (`Adding = true`) or LOSING
/// (`Adding = false`) one item, evaluated from its CURRENT accumulators
/// without mutating them -- an O(1) trial, mode-dispatched behind one uniform
/// scalar: isotropic badness is `-|sgamma|^2/srho` (a sign flip of the
/// maximized energy form, so accepting the lower badness is exactly
/// accepting the higher energy -- no behaviour change from that framing);
/// quadric badness is the summed quadric error at the candidate cluster's own
/// optimal point.
double remesh_trial_badness(RemeshMetric Metric, const double* pSgamma, double Srho,
                            const double* pQuad, const double* pItemWeighted, double ItemArea,
                            const double* pItemQuad, bool Adding) {
    const double sign = Adding ? 1.0 : -1.0;
    double new_sgamma[3] = {pSgamma[0] + sign * pItemWeighted[0],
                            pSgamma[1] + sign * pItemWeighted[1],
                            pSgamma[2] + sign * pItemWeighted[2]};
    const double new_srho = Srho + sign * ItemArea;
    const double iso_badness = -(new_sgamma[0] * new_sgamma[0] + new_sgamma[1] * new_sgamma[1] +
                                 new_sgamma[2] * new_sgamma[2]) /
                               new_srho;

    if (Metric == RemeshMetric::Isotropic)
        return iso_badness;

    // Quadric mode ADDS the quadric flatness error to the isotropic
    // compactness term rather than using the quadric error alone -- measured
    // (see remesh.hpp's doc comment): pure quadric-error badness has no
    // compactness pull at all, since a cluster confined to a common plane can
    // grow arbitrarily thin and snake-like along a low-curvature direction
    // without its quadric error rising, which repeatedly produced
    // disconnected, non-manifold clusters even after every repair pass (17 of
    // 150 clusters still split on a smooth closed sphere in gtest). The
    // isotropic term is the SAME `-|sgamma|^2/srho` form `Isotropic` mode
    // uses on its own -- both terms are area-weighted squared-distance-like
    // quantities (dimensionally consistent, no free blend weight to tune),
    // and their sum keeps the compactness pull that prevents pathological
    // elongation while still rewarding flatness/feature alignment where a
    // real crease or corner makes the quadric term the deciding factor.
    double new_quad[10];
    for (int k = 0; k < 10; ++k)
        new_quad[k] = pQuad[k] + sign * pItemQuad[k];
    double centroid[3];
    remesh_centroid_of(new_sgamma, new_srho, centroid);
    return iso_badness + remesh_quadric_badness(new_quad, centroid);
}

/// The accumulators one cluster carries. `mQuad` is filled and read only
/// under `RemeshMetric::Quadric`; `Isotropic` mode never touches it.
struct RemeshState {
    std::vector<double> mSgamma;  ///< `3 * n_clus`, sum of weighted positions.
    std::vector<double> mSrho;    ///< `n_clus`, sum of weights.
    std::vector<double> mBadness; ///< `n_clus`, the metric's own badness.
    std::vector<double> mQuad;    ///< `10 * n_clus`, quadric mode only.
    std::vector<std::int64_t> mCount;
};

RemeshState remesh_init_state(const std::vector<std::int64_t>& rClusters,
                              const std::vector<double>& rArea,
                              const std::vector<double>& rWeighted,
                              const std::vector<double>& rVertexQuad, RemeshMetric Metric,
                              std::int64_t NumClusters, std::size_t NumPoints) {
    RemeshState st;
    const std::size_t n_clus = static_cast<std::size_t>(NumClusters);
    st.mSgamma.assign(n_clus * 3, 0.0);
    st.mSrho.assign(n_clus, 0.0);
    st.mBadness.assign(n_clus, 0.0);
    st.mCount.assign(n_clus, 0);
    if (Metric == RemeshMetric::Quadric)
        st.mQuad.assign(n_clus * 10, 0.0);

    for (std::size_t i = 0; i < NumPoints; ++i) {
        const std::size_t c = static_cast<std::size_t>(rClusters[i]);
        st.mSrho[c] += rArea[i];
        st.mSgamma[c * 3] += rWeighted[i * 3];
        st.mSgamma[c * 3 + 1] += rWeighted[i * 3 + 1];
        st.mSgamma[c * 3 + 2] += rWeighted[i * 3 + 2];
        ++st.mCount[c];
        if (Metric == RemeshMetric::Quadric)
            for (int k = 0; k < 10; ++k)
                st.mQuad[c * 10 + static_cast<std::size_t>(k)] +=
                    rVertexQuad[i * 10 + static_cast<std::size_t>(k)];
    }
    for (std::size_t c = 0; c < n_clus; ++c) {
        if (Metric == RemeshMetric::Isotropic) {
            st.mBadness[c] = -(st.mSgamma[c * 3] * st.mSgamma[c * 3] +
                               st.mSgamma[c * 3 + 1] * st.mSgamma[c * 3 + 1] +
                               st.mSgamma[c * 3 + 2] * st.mSgamma[c * 3 + 2]) /
                             st.mSrho[c];
        } else {
            const double iso = -(st.mSgamma[c * 3] * st.mSgamma[c * 3] +
                                 st.mSgamma[c * 3 + 1] * st.mSgamma[c * 3 + 1] +
                                 st.mSgamma[c * 3 + 2] * st.mSgamma[c * 3 + 2]) /
                               st.mSrho[c];
            double centroid[3];
            remesh_centroid_of(st.mSgamma.data() + c * 3, st.mSrho[c], centroid);
            st.mBadness[c] = iso + remesh_quadric_badness(st.mQuad.data() + c * 10, centroid);
        }
    }
    return st;
}

/// Moves item `Item` from `From` to `To`: updates every accumulator, then
/// recomputes both clusters' badness under `Metric` from the updated state
/// (not the trial value already computed by the caller, to keep this
/// function correct on its own).
void remesh_commit(RemeshState& rSt, std::vector<std::int64_t>& rClusters,
                   const std::vector<double>& rArea, const std::vector<double>& rWeighted,
                   const std::vector<double>& rVertexQuad, RemeshMetric Metric, std::size_t Item,
                   std::int64_t From, std::int64_t To) {
    const std::size_t uf = static_cast<std::size_t>(From);
    const std::size_t ut = static_cast<std::size_t>(To);
    rClusters[Item] = To;
    --rSt.mCount[uf];
    ++rSt.mCount[ut];
    rSt.mSrho[ut] += rArea[Item];
    rSt.mSrho[uf] -= rArea[Item];
    for (int d = 0; d < 3; ++d) {
        rSt.mSgamma[ut * 3 + d] += rWeighted[Item * 3 + d];
        rSt.mSgamma[uf * 3 + d] -= rWeighted[Item * 3 + d];
    }
    if (Metric == RemeshMetric::Isotropic) {
        rSt.mBadness[uf] = -(rSt.mSgamma[uf * 3] * rSt.mSgamma[uf * 3] +
                             rSt.mSgamma[uf * 3 + 1] * rSt.mSgamma[uf * 3 + 1] +
                             rSt.mSgamma[uf * 3 + 2] * rSt.mSgamma[uf * 3 + 2]) /
                           rSt.mSrho[uf];
        rSt.mBadness[ut] = -(rSt.mSgamma[ut * 3] * rSt.mSgamma[ut * 3] +
                             rSt.mSgamma[ut * 3 + 1] * rSt.mSgamma[ut * 3 + 1] +
                             rSt.mSgamma[ut * 3 + 2] * rSt.mSgamma[ut * 3 + 2]) /
                           rSt.mSrho[ut];
    } else {
        for (int k = 0; k < 10; ++k) {
            rSt.mQuad[ut * 10 + static_cast<std::size_t>(k)] +=
                rVertexQuad[Item * 10 + static_cast<std::size_t>(k)];
            rSt.mQuad[uf * 10 + static_cast<std::size_t>(k)] -=
                rVertexQuad[Item * 10 + static_cast<std::size_t>(k)];
        }
        double cen_f[3], cen_t[3];
        remesh_centroid_of(rSt.mSgamma.data() + uf * 3, rSt.mSrho[uf], cen_f);
        remesh_centroid_of(rSt.mSgamma.data() + ut * 3, rSt.mSrho[ut], cen_t);
        const double iso_f = -(rSt.mSgamma[uf * 3] * rSt.mSgamma[uf * 3] +
                               rSt.mSgamma[uf * 3 + 1] * rSt.mSgamma[uf * 3 + 1] +
                               rSt.mSgamma[uf * 3 + 2] * rSt.mSgamma[uf * 3 + 2]) /
                             rSt.mSrho[uf];
        const double iso_t = -(rSt.mSgamma[ut * 3] * rSt.mSgamma[ut * 3] +
                               rSt.mSgamma[ut * 3 + 1] * rSt.mSgamma[ut * 3 + 1] +
                               rSt.mSgamma[ut * 3 + 2] * rSt.mSgamma[ut * 3 + 2]) /
                             rSt.mSrho[ut];
        rSt.mBadness[uf] = iso_f + remesh_quadric_badness(rSt.mQuad.data() + uf * 10, cen_f);
        rSt.mBadness[ut] = iso_t + remesh_quadric_badness(rSt.mQuad.data() + ut * 10, cen_t);
    }
}

/// Sweeps the boundary edges until nothing moves or `MaxIterations` is
/// reached. Serial by necessity: each accepted move changes the accumulators
/// the next test reads. Clusters unmodified in the previous sweep are
/// skipped.
std::int64_t remesh_minimize(const std::vector<std::int64_t>& rEdges,
                             std::vector<std::int64_t>& rClusters,
                             const std::vector<double>& rArea,
                             const std::vector<double>& rWeighted,
                             const std::vector<double>& rVertexQuad, RemeshMetric Metric,
                             RemeshState& rSt, std::int64_t NumClusters, int MaxIterations) {
    const std::size_t n_edges = rEdges.size() / 2;
    const std::size_t n_clus = static_cast<std::size_t>(NumClusters);
    std::vector<char> mod_cur(n_clus, 1);
    std::vector<char> mod_next(n_clus, 0);

    std::int64_t iterations = 0;
    bool changed = true;
    while (changed && iterations < MaxIterations) {
        changed = false;
        std::fill(mod_next.begin(), mod_next.end(), 0);

        for (std::size_t e = 0; e < n_edges; ++e) {
            const std::size_t pa = static_cast<std::size_t>(rEdges[2 * e]);
            const std::size_t pb = static_cast<std::size_t>(rEdges[2 * e + 1]);
            const std::int64_t ca = rClusters[pa];
            const std::int64_t cb = rClusters[pb];
            if (ca == cb)
                continue;
            const std::size_t uca = static_cast<std::size_t>(ca);
            const std::size_t ucb = static_cast<std::size_t>(cb);
            if (!mod_cur[uca] && !mod_cur[ucb])
                continue;

            // A cluster of one item may not be emptied -- losing its last
            // item would divide by a zero srho (or leave a rank-deficient
            // empty quadric) and delete the cluster outright.
            const bool can_move_a = rSt.mCount[uca] > 1;
            const bool can_move_b = rSt.mCount[ucb] > 1;
            const double badness_orig = rSt.mBadness[uca] + rSt.mBadness[ucb];

            const double* q_a = Metric == RemeshMetric::Quadric ? rSt.mQuad.data() + uca * 10
                                                                 : nullptr;
            const double* q_b = Metric == RemeshMetric::Quadric ? rSt.mQuad.data() + ucb * 10
                                                                 : nullptr;
            const double* item_q_b =
                Metric == RemeshMetric::Quadric ? rVertexQuad.data() + pb * 10 : nullptr;
            const double* item_q_a =
                Metric == RemeshMetric::Quadric ? rVertexQuad.data() + pa * 10 : nullptr;

            // Candidate 1: move b into a (b leaves its own cluster).
            double badness_move_b = badness_orig;  // sentinel: "no better than orig"
            if (can_move_b)
                badness_move_b =
                    remesh_trial_badness(Metric, rSt.mSgamma.data() + uca * 3, rSt.mSrho[uca],
                                         q_a, rWeighted.data() + pb * 3, rArea[pb], item_q_b,
                                         /*Adding=*/true) +
                    remesh_trial_badness(Metric, rSt.mSgamma.data() + ucb * 3, rSt.mSrho[ucb],
                                         q_b, rWeighted.data() + pb * 3, rArea[pb], item_q_b,
                                         /*Adding=*/false);

            // Candidate 2: move a into b (a leaves its own cluster).
            double badness_move_a = badness_orig;
            if (can_move_a)
                badness_move_a =
                    remesh_trial_badness(Metric, rSt.mSgamma.data() + uca * 3, rSt.mSrho[uca],
                                         q_a, rWeighted.data() + pa * 3, rArea[pa], item_q_a,
                                         /*Adding=*/false) +
                    remesh_trial_badness(Metric, rSt.mSgamma.data() + ucb * 3, rSt.mSrho[ucb],
                                         q_b, rWeighted.data() + pa * 3, rArea[pa], item_q_a,
                                         /*Adding=*/true);

            // Accept the strictly better candidate; ties resolve toward moving
            // b -- a fixed, documented rule (this repo's determinism contract
            // requires one), not an accident of evaluation order.
            if (can_move_b && badness_move_b < badness_orig && badness_move_b <= badness_move_a) {
                remesh_commit(rSt, rClusters, rArea, rWeighted, rVertexQuad, Metric, pb, cb, ca);
                mod_next[uca] = mod_next[ucb] = 1;
                changed = true;
            } else if (can_move_a && badness_move_a < badness_orig) {
                remesh_commit(rSt, rClusters, rArea, rWeighted, rVertexQuad, Metric, pa, ca, cb);
                mod_next[uca] = mod_next[ucb] = 1;
                changed = true;
            }
        }
        mod_cur.swap(mod_next);
        ++iterations;
    }
    return iterations;
}

/// Finds clusters split into several connected components, keeps the largest
/// (lowest starting vertex on a tie) and unassigns the rest. Returns how many
/// clusters were split.
std::int64_t remesh_split_disconnected(const detail::NodeAdjacency& rAdj,
                                       std::vector<std::int64_t>& rClusters,
                                       std::int64_t NumClusters, std::size_t NumPoints) {
    std::vector<char> visited(NumPoints, 0);
    std::vector<std::vector<std::int64_t>> components;
    std::vector<std::int64_t> best_size(static_cast<std::size_t>(NumClusters), -1);
    std::vector<std::int64_t> best_comp(static_cast<std::size_t>(NumClusters), -1);
    std::vector<std::int64_t> owner;
    std::vector<std::int64_t> stack;

    for (std::size_t seed = 0; seed < NumPoints; ++seed) {
        if (visited[seed] || rClusters[seed] < 0)
            continue;
        const std::int64_t c = rClusters[seed];
        std::vector<std::int64_t> comp;
        stack.clear();
        stack.push_back(static_cast<std::int64_t>(seed));
        visited[seed] = 1;
        while (!stack.empty()) {
            const std::int64_t p = stack.back();
            stack.pop_back();
            comp.push_back(p);
            for (std::int64_t k = rAdj.mXadj[static_cast<std::size_t>(p)];
                 k < rAdj.mXadj[static_cast<std::size_t>(p) + 1]; ++k) {
                const std::int64_t nbr = rAdj.mAdj[static_cast<std::size_t>(k)];
                const std::size_t un = static_cast<std::size_t>(nbr);
                if (!visited[un] && rClusters[un] == c) {
                    visited[un] = 1;
                    stack.push_back(nbr);
                }
            }
        }
        const std::int64_t size = static_cast<std::int64_t>(comp.size());
        const std::size_t uc = static_cast<std::size_t>(c);
        if (size > best_size[uc]) {
            best_size[uc] = size;
            best_comp[uc] = static_cast<std::int64_t>(components.size());
        }
        owner.push_back(c);
        components.push_back(std::move(comp));
    }

    std::int64_t num_split = 0;
    std::vector<char> counted(static_cast<std::size_t>(NumClusters), 0);
    for (std::size_t ci = 0; ci < components.size(); ++ci) {
        const std::size_t uc = static_cast<std::size_t>(owner[ci]);
        if (best_comp[uc] == static_cast<std::int64_t>(ci))
            continue;
        for (const std::int64_t p : components[ci])
            rClusters[static_cast<std::size_t>(p)] = -1;
        if (!counted[uc]) {
            counted[uc] = 1;
            ++num_split;
        }
    }
    return num_split;
}

/// Compacts cluster ids to `[0, out)` preserving ascending order, returning the
/// number that survived.
std::int64_t remesh_renumber(std::vector<std::int64_t>& rClusters, std::int64_t NumClusters,
                             std::size_t NumPoints) {
    std::vector<std::int64_t> remap(static_cast<std::size_t>(NumClusters), -1);
    for (std::size_t i = 0; i < NumPoints; ++i)
        if (rClusters[i] >= 0)
            remap[static_cast<std::size_t>(rClusters[i])] = 1;
    std::int64_t out = 0;
    for (std::size_t c = 0; c < static_cast<std::size_t>(NumClusters); ++c)
        if (remap[c] > 0)
            remap[c] = out++;
    for (std::size_t i = 0; i < NumPoints; ++i)
        if (rClusters[i] >= 0)
            rClusters[i] = remap[static_cast<std::size_t>(rClusters[i])];
    return out;
}

/// The output vertex position of each FINAL (post-renumber) cluster: a fresh
/// accumulation pass keyed on the renumbered `rClusters` array (the old
/// `RemeshState`'s accumulators are indexed by the pre-renumber id space and
/// cannot be reused here). Isotropic mode's position is the area-weighted
/// centroid; quadric mode's is `decim_quadric_optimal_point` of the summed
/// per-cluster quadric, falling back to that same centroid.
std::vector<double> remesh_final_positions(const std::vector<std::int64_t>& rClusters,
                                           const std::vector<double>& rArea,
                                           const std::vector<double>& rWeighted,
                                           const std::vector<double>& rVertexQuad,
                                           RemeshMetric Metric, std::int64_t NumClusters,
                                           std::size_t NumPoints) {
    const std::size_t n_clus = static_cast<std::size_t>(NumClusters);
    std::vector<double> sgamma(n_clus * 3, 0.0);
    std::vector<double> srho(n_clus, 0.0);
    std::vector<double> quad;
    if (Metric == RemeshMetric::Quadric)
        quad.assign(n_clus * 10, 0.0);

    for (std::size_t i = 0; i < NumPoints; ++i) {
        const std::size_t c = static_cast<std::size_t>(rClusters[i]);
        srho[c] += rArea[i];
        for (int d = 0; d < 3; ++d)
            sgamma[c * 3 + d] += rWeighted[i * 3 + d];
        if (Metric == RemeshMetric::Quadric)
            for (int k = 0; k < 10; ++k)
                quad[c * 10 + static_cast<std::size_t>(k)] +=
                    rVertexQuad[i * 10 + static_cast<std::size_t>(k)];
    }

    std::vector<double> out(n_clus * 3, 0.0);
    for (std::size_t c = 0; c < n_clus; ++c) {
        double centroid[3];
        remesh_centroid_of(sgamma.data() + c * 3, srho[c], centroid);
        if (Metric == RemeshMetric::Isotropic) {
            out[c * 3] = centroid[0];
            out[c * 3 + 1] = centroid[1];
            out[c * 3 + 2] = centroid[2];
        } else {
            double pt[3];
            detail::decim_quadric_optimal_point(quad.data() + c * 10, centroid, pt);
            out[c * 3] = pt[0];
            out[c * 3 + 1] = pt[1];
            out[c * 3 + 2] = pt[2];
        }
    }
    return out;
}

/// Mean unit normal of each cluster, used only to orient the dual triangles.
/// Metric-agnostic: it reads the input surface's own triangle normals, never
/// a cluster's quadric.
std::vector<double> remesh_cluster_normals(const RemeshSurface& rS,
                                           const std::vector<std::int64_t>& rClusters,
                                           std::int64_t NumClusters) {
    std::vector<double> nrm(static_cast<std::size_t>(NumClusters) * 3, 0.0);
    for (std::size_t f = 0; f < rS.mNumFaces; ++f) {
        const double* v0 = rS.mXyz.data() + static_cast<std::size_t>(rS.mCorners[f * 3]) * 3;
        const double* v1 = rS.mXyz.data() + static_cast<std::size_t>(rS.mCorners[f * 3 + 1]) * 3;
        const double* v2 = rS.mXyz.data() + static_cast<std::size_t>(rS.mCorners[f * 3 + 2]) * 3;
        const double e0[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        const double e1[3] = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
        const double c[3] = {e0[1] * e1[2] - e0[2] * e1[1], e0[2] * e1[0] - e0[0] * e1[2],
                             e0[0] * e1[1] - e0[1] * e1[0]};
        for (int k = 0; k < 3; ++k) {
            const std::size_t uc = static_cast<std::size_t>(rClusters[static_cast<std::size_t>(
                rS.mCorners[f * 3 + k])]);
            for (int d = 0; d < 3; ++d)
                nrm[uc * 3 + d] += c[d];
        }
    }
    for (std::size_t c = 0; c < static_cast<std::size_t>(NumClusters); ++c) {
        const double len = std::sqrt(nrm[c * 3] * nrm[c * 3] + nrm[c * 3 + 1] * nrm[c * 3 + 1] +
                                     nrm[c * 3 + 2] * nrm[c * 3 + 2]);
        if (len > 0.0)
            for (int d = 0; d < 3; ++d)
                nrm[c * 3 + d] /= len;
    }
    return nrm;
}

/// Builds the dual: one triangle per input face spanning three distinct
/// clusters, deduplicated on the sorted id triple and emitted in ascending
/// order, then wound to agree with the mean normal of its three clusters.
std::vector<std::int64_t> remesh_dual_faces(const RemeshSurface& rS,
                                            const std::vector<std::int64_t>& rClusters,
                                            const std::vector<double>& rPositions,
                                            const std::vector<double>& rNormals) {
    std::vector<std::array<std::int64_t, 3>> keys;
    keys.reserve(rS.mNumFaces);
    for (std::size_t f = 0; f < rS.mNumFaces; ++f) {
        std::array<std::int64_t, 3> t = {
            rClusters[static_cast<std::size_t>(rS.mCorners[f * 3])],
            rClusters[static_cast<std::size_t>(rS.mCorners[f * 3 + 1])],
            rClusters[static_cast<std::size_t>(rS.mCorners[f * 3 + 2])]};
        std::sort(t.begin(), t.end());
        if (t[0] == t[1] || t[1] == t[2])
            continue;  // the face does not join three distinct clusters
        keys.push_back(t);
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

    std::vector<std::int64_t> out;
    out.reserve(keys.size() * 3);
    for (const auto& t : keys) {
        double ref[3] = {0.0, 0.0, 0.0};
        for (int k = 0; k < 3; ++k)
            for (int d = 0; d < 3; ++d)
                ref[d] += rNormals[static_cast<std::size_t>(t[k]) * 3 + d];

        const double* a = rPositions.data() + static_cast<std::size_t>(t[0]) * 3;
        const double* b = rPositions.data() + static_cast<std::size_t>(t[1]) * 3;
        const double* c = rPositions.data() + static_cast<std::size_t>(t[2]) * 3;
        const double e0[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double e1[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double nx = e0[1] * e1[2] - e0[2] * e1[1];
        const double ny = e0[2] * e1[0] - e0[0] * e1[2];
        const double nz = e0[0] * e1[1] - e0[1] * e1[0];
        const double dot = nx * ref[0] + ny * ref[1] + nz * ref[2];

        if (dot < 0.0) {
            out.push_back(t[2]);
            out.push_back(t[1]);
            out.push_back(t[0]);
        } else {
            out.push_back(t[0]);
            out.push_back(t[1]);
            out.push_back(t[2]);
        }
    }
    return out;
}

}  // namespace

RemeshMetric remesh_metric_from_name(const std::string& rName) {
    if (rName == "isotropic")
        return RemeshMetric::Isotropic;
    if (rName == "quadric")
        return RemeshMetric::Quadric;
    throw std::invalid_argument("meshio++: remesh: unknown metric '" + rName +
                                "' (expected 'isotropic' or 'quadric')");
}

int remesh_suggest_subdivide(std::int64_t NumPoints, std::int64_t NumClusters,
                             double SubsampleRatio, int MaxSubdivide) {
    if (NumPoints <= 0 || NumClusters <= 0 || MaxSubdivide <= 0)
        return 0;
    const double target = SubsampleRatio * static_cast<double>(NumClusters);
    double have = static_cast<double>(NumPoints);
    int n = 0;
    while (have < target && n < MaxSubdivide) {
        // A uniform triangle refinement quadruples the face count; the vertex
        // count grows to roughly 4x too on a large mesh (V - E + F is fixed,
        // and each pass adds one vertex per edge).
        have *= 4.0;
        ++n;
    }
    return n;
}

RemeshResult remesh(const Mesh& rMesh, const RemeshOptions& rOptions) {
    if (rOptions.mNumClusters < 4)
        throw std::invalid_argument(
            "meshio++: remesh: num_clusters must be at least 4, got " +
            std::to_string(rOptions.mNumClusters) +
            " (fewer vertices cannot bound a closed surface)");
    if (rOptions.mMaxIterations <= 0)
        throw std::invalid_argument("meshio++: remesh: max_iterations must be positive, got " +
                                    std::to_string(rOptions.mMaxIterations));
    if (rOptions.mMaxRepairPasses < 0)
        throw std::invalid_argument(
            "meshio++: remesh: max_repair_passes must not be negative, got " +
            std::to_string(rOptions.mMaxRepairPasses));
    if (rOptions.mSubsampleRatio <= 0.0)
        throw std::invalid_argument("meshio++: remesh: subsample_ratio must be positive");

    remesh_check_blocks(rMesh);
    detail::warn_regions_dropped(rMesh, "remesh");

    RemeshResult result;
    const RemeshMetric metric = rOptions.mMetric;

    // Triangulate, then subdivide. Both are existing operations: a quad or a
    // rectangular polygon block is fanned by convert_cells, and the
    // subsampling step is exactly a uniform refine pass.
    Mesh work =
        convert_cells(rMesh, {ConvertCellsMode::Simplexify, /*mRecordParentIds=*/false}).mMesh;

    int subdivide = rOptions.mSubdivide;
    if (subdivide < 0)
        subdivide = remesh_suggest_subdivide(static_cast<std::int64_t>(work.NumPoints()),
                                             rOptions.mNumClusters, rOptions.mSubsampleRatio,
                                             rOptions.mMaxSubdivide);
    if (subdivide > 0) {
        RefineOptions ropts;
        ropts.mLevels = subdivide;
        work = refine(work, ropts).mMesh;
    }
    result.mSubdivideApplied = subdivide;

    RemeshSurface surf = remesh_build_surface(work);
    if (surf.mNumFaces == 0)
        throw std::invalid_argument("meshio++: remesh: mesh has no triangles to remesh");
    if (rOptions.mNumClusters > static_cast<std::int64_t>(surf.mNumPoints))
        throw std::invalid_argument(
            "meshio++: remesh: num_clusters (" + std::to_string(rOptions.mNumClusters) +
            ") exceeds the " + std::to_string(surf.mNumPoints) +
            " vertices available after subdivision; raise subdivide or lower num_clusters");

    std::vector<double> area;
    std::vector<double> weighted;
    remesh_item_weights(surf, area, weighted);
    const std::vector<double> vertex_quad = remesh_vertex_quadrics(surf, metric);

    const detail::NodeAdjacency adj =
        detail::build_node_adjacency(work, surf.mNumPoints, detail::NodeAdjacencyKind::Edge);
    const std::vector<std::int64_t> edges = remesh_unique_edges(adj, surf.mNumPoints);

    std::vector<std::int64_t> clusters(surf.mNumPoints, -1);
    remesh_seed_clusters(clusters, adj, area, rOptions.mNumClusters, surf.mNumPoints);
    if (remesh_grow_null(edges, clusters))
        for (std::size_t i = 0; i < surf.mNumPoints; ++i)
            if (clusters[i] == -1)
                clusters[i] = 0;

    RemeshState state = remesh_init_state(clusters, area, weighted, vertex_quad, metric,
                                          rOptions.mNumClusters, surf.mNumPoints);
    result.mNumIterations =
        remesh_minimize(edges, clusters, area, weighted, vertex_quad, metric, state,
                        rOptions.mNumClusters, rOptions.mMaxIterations);

    // Repair: a minimisation sweep can leave a cluster in two pieces, whose
    // dual is where non-manifold output comes from. Split, re-grow, minimise
    // again, and report whatever is still split at the end.
    std::int64_t disconnected = remesh_split_disconnected(adj, clusters, rOptions.mNumClusters,
                                                          surf.mNumPoints);
    for (int pass = 0; pass < rOptions.mMaxRepairPasses && disconnected > 0; ++pass) {
        if (remesh_grow_null(edges, clusters))
            for (std::size_t i = 0; i < surf.mNumPoints; ++i)
                if (clusters[i] == -1)
                    clusters[i] = 0;
        state = remesh_init_state(clusters, area, weighted, vertex_quad, metric,
                                  rOptions.mNumClusters, surf.mNumPoints);
        result.mNumIterations +=
            remesh_minimize(edges, clusters, area, weighted, vertex_quad, metric, state,
                            rOptions.mNumClusters, rOptions.mMaxIterations);
        disconnected = remesh_split_disconnected(adj, clusters, rOptions.mNumClusters,
                                                 surf.mNumPoints);
    }
    result.mNumIsolatedClusters = disconnected;
    if (remesh_grow_null(edges, clusters))
        for (std::size_t i = 0; i < surf.mNumPoints; ++i)
            if (clusters[i] == -1)
                clusters[i] = 0;

    const std::int64_t n_clus =
        remesh_renumber(clusters, rOptions.mNumClusters, surf.mNumPoints);
    result.mNumClusters = n_clus;
    if (n_clus < rOptions.mNumClusters)
        log::warn("remesh: produced {} clusters of the {} requested (the surface could not be "
                  "split more finely)",
                  n_clus, rOptions.mNumClusters);

    const std::vector<double> positions = remesh_final_positions(
        clusters, area, weighted, vertex_quad, metric, n_clus, surf.mNumPoints);
    const std::vector<double> normals = remesh_cluster_normals(surf, clusters, n_clus);
    const std::vector<std::int64_t> faces = remesh_dual_faces(surf, clusters, positions, normals);

    Mesh& out = result.mMesh;
    {
        NDArray pts = NDArray::Uninit(DType::Float64, {static_cast<std::size_t>(n_clus), 3});
        double* dst = pts.As<double>();
        for (std::size_t i = 0; i < static_cast<std::size_t>(n_clus) * 3; ++i)
            dst[i] = positions[i];
        out.AssignPoints(std::move(pts));
    }
    {
        const std::size_t nf = faces.size() / 3;
        NDArray conn = NDArray::Uninit(DType::Int64, {nf, 3});
        std::int64_t* dst = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < faces.size(); ++i)
            dst[i] = faces[i];
        out.AddCellBlock(cell_type_name(CellType::Triangle), std::move(conn));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, rMesh.FieldData(name));

    return result;
}

}  // namespace meshioplusplus
