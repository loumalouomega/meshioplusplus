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
#include <map>
#include <set>
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
/// to the surface area -- unless curvature gradation rescales it, see below)
/// and the weighted position each cluster accumulates.
///
/// When `Gradation != 0.0` and `rCurvature` is non-empty, `rArea` (which
/// becomes the general "item weight" from this point on, not literal area,
/// everywhere it is used downstream -- seeding, badness, final placement)
/// is rescaled to `area * max(kappa, kappa_floor)^Gradation`, where
/// `kappa_floor` is a small internal fraction of the mesh's own maximum
/// curvature -- a defensive floor, not a tunable, guarding against an
/// exact-zero weight degenerating a flat vertex's cluster membership (the
/// same weight-clamping idea the real ACVD reference documents for its own
/// isotropic weighting). A perfectly flat mesh (`kappa` identically zero
/// everywhere) leaves `rArea` untouched -- there is nothing to grade.
void remesh_item_weights(const RemeshSurface& rS, const std::vector<double>& rCurvature,
                         double Gradation, std::vector<double>& rArea,
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

    if (Gradation != 0.0 && !rCurvature.empty()) {
        double kappa_max = 0.0;
        for (const double k : rCurvature)
            kappa_max = std::max(kappa_max, k);
        if (kappa_max > 0.0) {
            const double kappa_floor = kappa_max * 1e-3;
            for (std::size_t i = 0; i < rS.mNumPoints; ++i)
                rArea[i] *= std::pow(std::max(rCurvature[i], kappa_floor), Gradation);
        }
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

/// Each vertex's own curvature-tensor quadric (`RemeshMetric::Anisotropic`
/// only), packed into the SAME 10-entry layout `remesh_vertex_quadrics`
/// builds for `Quadric` mode -- the sibling that construction, not a second
/// accumulator shape. `E(v) = sum(w * (x-v)^T M (x-v))` expands to
/// `x^T A x + 2 b.x + c` with `A = sum(w*M)`, `b = -sum(w*M*x)`, `c =
/// sum(w * x^T M x)` -- exactly `decim_quadric_error`'s own convention, so
/// the same `decim_quadric_optimal_point` solves it with no new solver.
/// `rArea` must already be computed (the item weight): folded in HERE
/// rather than in `remesh_vertex_curvature`, which has no notion of
/// clustering weights and is also used, tensor-free, under `mGradation`
/// alone.
std::vector<double> remesh_vertex_metrics(const RemeshSurface& rS,
                                          const std::vector<double>& rTensors,
                                          const std::vector<double>& rArea) {
    std::vector<double> quad(rS.mNumPoints * 10, 0.0);
    for (std::size_t i = 0; i < rS.mNumPoints; ++i) {
        const double* m = rTensors.data() + i * 6;  // [xx,xy,xz,yy,yz,zz]
        const double* x = rS.mXyz.data() + i * 3;
        const double w = rArea[i];

        const double mx0 = m[0] * x[0] + m[1] * x[1] + m[2] * x[2];
        const double mx1 = m[1] * x[0] + m[3] * x[1] + m[4] * x[2];
        const double mx2 = m[2] * x[0] + m[4] * x[1] + m[5] * x[2];
        const double c_term = w * (x[0] * mx0 + x[1] * mx1 + x[2] * mx2);

        double* q = quad.data() + i * 10;
        q[0] = w * m[0];
        q[1] = w * m[1];
        q[2] = w * m[2];
        q[3] = -w * mx0;
        q[4] = w * m[3];
        q[5] = w * m[4];
        q[6] = -w * mx1;
        q[7] = w * m[5];
        q[8] = -w * mx2;
        q[9] = c_term;
    }
    return quad;
}

/// Area-weighted unit normal at each vertex, from its incident face normals
/// (the per-vertex analogue of `remesh_cluster_normals` below). Only computed
/// when curvature gradation is requested.
std::vector<double> remesh_vertex_normals(const RemeshSurface& rS) {
    std::vector<double> nrm(rS.mNumPoints * 3, 0.0);
    for (std::size_t f = 0; f < rS.mNumFaces; ++f) {
        const double* v0 = rS.mXyz.data() + static_cast<std::size_t>(rS.mCorners[f * 3]) * 3;
        const double* v1 = rS.mXyz.data() + static_cast<std::size_t>(rS.mCorners[f * 3 + 1]) * 3;
        const double* v2 = rS.mXyz.data() + static_cast<std::size_t>(rS.mCorners[f * 3 + 2]) * 3;
        const double e0[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        const double e1[3] = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
        const double c[3] = {e0[1] * e1[2] - e0[2] * e1[1], e0[2] * e1[0] - e0[0] * e1[2],
                             e0[0] * e1[1] - e0[1] * e1[0]};
        for (int k = 0; k < 3; ++k) {
            const std::size_t vi = static_cast<std::size_t>(rS.mCorners[f * 3 + k]);
            for (int d = 0; d < 3; ++d)
                nrm[vi * 3 + d] += c[d];
        }
    }
    for (std::size_t i = 0; i < rS.mNumPoints; ++i) {
        const double len = std::sqrt(nrm[i * 3] * nrm[i * 3] + nrm[i * 3 + 1] * nrm[i * 3 + 1] +
                                     nrm[i * 3 + 2] * nrm[i * 3 + 2]);
        if (len > 0.0)
            for (int d = 0; d < 3; ++d)
                nrm[i * 3 + d] /= len;
    }
    return nrm;
}

/// Solves the symmetric 3x3 system `M x = b` via Cramer's rule with explicit
/// 3x3 determinants (never the cofactor/adjugate shortcut -- this system has
/// no special structure to exploit, unlike `decim_quadric_optimal_point`'s).
/// Returns `false` (leaving `pOut` untouched) when the system is
/// ill-conditioned relative to its own scale.
bool remesh_solve_sym3(double m00, double m01, double m02, double m11, double m12, double m22,
                       double b0, double b1, double b2, double pOut[3]) {
    const double det = m00 * (m11 * m22 - m12 * m12) - m01 * (m01 * m22 - m12 * m02) +
                       m02 * (m01 * m12 - m11 * m02);
    double scale = std::abs(m00);
    scale = std::max(scale, std::abs(m01));
    scale = std::max(scale, std::abs(m02));
    scale = std::max(scale, std::abs(m11));
    scale = std::max(scale, std::abs(m12));
    scale = std::max(scale, std::abs(m22));
    if (std::abs(det) <= 1e-12 * (scale * scale * scale))
        return false;

    const double det0 = b0 * (m11 * m22 - m12 * m12) - m01 * (b1 * m22 - m12 * b2) +
                        m02 * (b1 * m12 - m11 * b2);
    const double det1 = m00 * (b1 * m22 - m12 * b2) - b0 * (m01 * m22 - m12 * m02) +
                        m02 * (m01 * b2 - b1 * m02);
    const double det2 = m00 * (m11 * b2 - b1 * m12) - m01 * (m01 * b2 - b1 * m02) +
                        b0 * (m01 * m12 - m11 * m02);
    pOut[0] = det0 / det;
    pOut[1] = det1 / det;
    pOut[2] = det2 / det;
    return true;
}

/// Per-vertex curvature magnitude via a local osculating-paraboloid fit over
/// the 1-ring (the same `detail::NodeAdjacency` graph the edge sweep already
/// builds -- no separate neighbourhood machinery). At vertex `v` with unit
/// normal `n` and tangent basis `(u, w)`, each 1-ring neighbour is projected
/// into local `(u, w, h)` coordinates (`h` = height out of the tangent
/// plane) and `h = a*u^2 + b*u*w + c*w^2` is fit by least squares; `kappa` is
/// the larger-magnitude eigenvalue of `[[2a, b], [b, 2c]]` (the dominant
/// bend). A low-valence (<3) or near-degenerate neighbourhood falls back to
/// `kappa = 0` -- the same ill-conditioning-means-flat convention
/// `decim_quadric_optimal_point` already established.
///
/// `pTensors`, when non-null, is filled with the SAME fit's principal
/// directions instead of discarding them: `RemeshMetric::Anisotropic`'s only
/// consumer. Six doubles per vertex, the upper triangle `[xx,xy,xz,yy,yz,zz]`
/// of a world-space symmetric 3x3 -- see `remesh_vertex_metrics` for what it
/// is built into. Every early-return path below (low valence, zero normal,
/// degenerate tangent, ill-conditioned fit) writes the identity, matching
/// `kappa = 0`'s "flat" convention: no curvature signal, no shape signal.
std::vector<double> remesh_vertex_curvature(const RemeshSurface& rS,
                                            const detail::NodeAdjacency& rAdj,
                                            const std::vector<double>& rNormals,
                                            std::vector<double>* pTensors = nullptr,
                                            double MaxAnisotropy = kRemeshDefaultMaxAnisotropy) {
    if (pTensors) {
        pTensors->assign(rS.mNumPoints * 6, 0.0);
        for (std::size_t i = 0; i < rS.mNumPoints; ++i) {
            (*pTensors)[i * 6] = 1.0;
            (*pTensors)[i * 6 + 3] = 1.0;
            (*pTensors)[i * 6 + 5] = 1.0;
        }
    }
    std::vector<double> kappa(rS.mNumPoints, 0.0);
    for (std::size_t i = 0; i < rS.mNumPoints; ++i) {
        const std::int64_t lo = rAdj.mXadj[i];
        const std::int64_t hi = rAdj.mXadj[i + 1];
        if (hi - lo < 3)
            continue;

        const double* p0 = rS.mXyz.data() + i * 3;
        const double* n = rNormals.data() + i * 3;
        if (n[0] == 0.0 && n[1] == 0.0 && n[2] == 0.0)
            continue;  // an unreferenced or degenerate-fan vertex

        // An orthonormal tangent basis (u, w) perpendicular to n: seed u from
        // whichever coordinate axis is least aligned with n (so the
        // Gram-Schmidt subtraction below never nearly cancels), then
        // Gram-Schmidt it against n and take w = n x u.
        const double ax[3] = {1.0, 0.0, 0.0};
        const double ay[3] = {0.0, 1.0, 0.0};
        const double* seed = (std::abs(n[0]) < 0.9) ? ax : ay;
        const double dp = seed[0] * n[0] + seed[1] * n[1] + seed[2] * n[2];
        double u[3] = {seed[0] - dp * n[0], seed[1] - dp * n[1], seed[2] - dp * n[2]};
        const double ulen = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        if (ulen <= 0.0)
            continue;
        u[0] /= ulen;
        u[1] /= ulen;
        u[2] /= ulen;
        const double w[3] = {n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2],
                             n[0] * u[1] - n[1] * u[0]};

        double m00 = 0.0, m01 = 0.0, m02 = 0.0, m11 = 0.0, m12 = 0.0, m22 = 0.0;
        double b0 = 0.0, b1 = 0.0, b2 = 0.0;
        for (std::int64_t k = lo; k < hi; ++k) {
            const std::size_t j = static_cast<std::size_t>(rAdj.mAdj[static_cast<std::size_t>(k)]);
            const double d[3] = {rS.mXyz[j * 3] - p0[0], rS.mXyz[j * 3 + 1] - p0[1],
                                 rS.mXyz[j * 3 + 2] - p0[2]};
            const double du = d[0] * u[0] + d[1] * u[1] + d[2] * u[2];
            const double dw = d[0] * w[0] + d[1] * w[1] + d[2] * w[2];
            const double dh = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
            const double r0 = du * du, r1 = du * dw, r2 = dw * dw;
            m00 += r0 * r0;
            m01 += r0 * r1;
            m02 += r0 * r2;
            m11 += r1 * r1;
            m12 += r1 * r2;
            m22 += r2 * r2;
            b0 += r0 * dh;
            b1 += r1 * dh;
            b2 += r2 * dh;
        }

        double abc[3];
        if (!remesh_solve_sym3(m00, m01, m02, m11, m12, m22, b0, b1, b2, abc))
            continue;  // ill-conditioned -> kappa stays 0

        // Principal curvatures = eigenvalues of [[2a, b], [b, 2c]].
        const double a = abc[0], b = abc[1], c = abc[2];
        const double tr = 2.0 * a + 2.0 * c;
        const double det2x2 = 4.0 * a * c - b * b;
        const double disc = std::max(0.0, tr * tr - 4.0 * det2x2);
        const double sq = std::sqrt(disc);
        const double k1 = 0.5 * (tr + sq);
        const double k2 = 0.5 * (tr - sq);
        kappa[i] = std::max(std::abs(k1), std::abs(k2));

        if (pTensors) {
            // The 2x2 eigenvector for k1 in the local (u, w) frame is
            // proportional to (b, k1 - 2a) -- or (k1 - 2c, b) when that is
            // better conditioned. k2's eigenvector is its exact in-plane
            // perpendicular: symmetric matrices with distinct eigenvalues
            // have orthogonal eigenvectors, so this needs no second solve.
            // At an umbilic (k1 == k2) any orthogonal pair is equally valid,
            // since only the OUTER PRODUCTS e*e^T enter the tensor below --
            // e and -e give the same e*e^T, and an umbilic's ratio is 1
            // regardless of which pair is picked.
            const double tensor_eps = 1e-20 * (a * a + b * b + c * c + 1.0);
            double alpha = b, beta = k1 - 2.0 * a;
            double len2 = alpha * alpha + beta * beta;
            if (len2 < tensor_eps) {
                alpha = k1 - 2.0 * c;
                beta = b;
                len2 = alpha * alpha + beta * beta;
            }
            if (len2 < tensor_eps) {
                alpha = 1.0;
                beta = 0.0;
                len2 = 1.0;
            }
            const double inv_len = 1.0 / std::sqrt(len2);
            alpha *= inv_len;
            beta *= inv_len;

            // World-space unit vectors: e1 for the k1 direction, e2 its
            // in-plane perpendicular (the k2 direction).
            const double e1[3] = {alpha * u[0] + beta * w[0], alpha * u[1] + beta * w[1],
                                  alpha * u[2] + beta * w[2]};
            const double e2[3] = {-beta * u[0] + alpha * w[0], -beta * u[1] + alpha * w[1],
                                  -beta * u[2] + alpha * w[2]};

            const double k1_abs = std::abs(k1), k2_abs = std::abs(k2);
            const bool k1_is_sharp = k1_abs >= k2_abs;
            const double hi_kappa = k1_is_sharp ? k1_abs : k2_abs;
            const double lo_kappa = k1_is_sharp ? k2_abs : k1_abs;
            const double* e_hi = k1_is_sharp ? e1 : e2;
            const double* e_lo = k1_is_sharp ? e2 : e1;

            // Length ratio (long axis / short axis) the tensor targets,
            // clamped to MaxAnisotropy. hi_kappa == 0 means no curvature
            // signal at all (an isotropic point), so ratio == 1; hi_kappa >
            // 0 with lo_kappa == 0 (the canonical cylinder) goes straight to
            // the clamp rather than dividing by zero -- the same formula's
            // limit, not a special case grafted on.
            double ratio = 1.0;
            if (hi_kappa > 0.0)
                ratio = (lo_kappa > 0.0)
                            ? std::min(std::sqrt(hi_kappa / lo_kappa), MaxAnisotropy)
                            : MaxAnisotropy;

            // Unit-determinant SPD tensor: in-plane eigenvalues (lam_hi,
            // lam_lo), normal eigenvalue pinned to lam_hi (keeps clusters on
            // the surface, no sharper than the sharpest in-plane direction
            // already demands) -- the pre-scale determinant is ratio^4, so
            // scaling every eigenvalue by ratio^(-4/3) makes det == 1 while
            // preserving the ratio itself (see remesh.hpp's doc comment).
            // ratio == 1 (flat, or MaxAnisotropy == 1.0) makes both
            // eigenvalues exactly 1.0, i.e. M == I.
            const double lam_hi = std::pow(ratio, 2.0 / 3.0);
            const double lam_lo = std::pow(ratio, -4.0 / 3.0);

            for (int r = 0; r < 3; ++r)
                for (int cc = r; cc < 3; ++cc) {
                    const double m_rc = lam_hi * (e_hi[r] * e_hi[cc] + n[r] * n[cc]) +
                                        lam_lo * e_lo[r] * e_lo[cc];
                    const std::size_t idx =
                        static_cast<std::size_t>(r == 0 ? cc : (r == 1 ? cc + 2 : 5));
                    (*pTensors)[i * 6 + idx] = m_rc;
                }
        }
    }
    return kappa;
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

/// The input's open boundary: an edge used by exactly one triangle is a
/// boundary edge; a vertex is a boundary vertex iff it is an endpoint of at
/// least one. A local, self-contained edge-use-count pass -- deliberately
/// NOT shared with `decimate.cpp`'s own `decim_build_edges`, whose header
/// comment already scopes that helper to decimate's two current callers
/// only ("stays private to decimate.cpp ... neither concept applies" to the
/// other one). `std::map` (not a hash map) is deliberate: it costs nothing
/// here and gives a canonical, deterministic edge-visitation order for free.
struct RemeshBoundaryInfo {
    std::vector<std::uint8_t> mIsBoundaryVertex;  ///< size NumPoints.
    std::vector<std::pair<std::int64_t, std::int64_t>> mBoundaryEdges;  ///< lo < hi, sorted.
};

RemeshBoundaryInfo remesh_boundary_info(const RemeshSurface& rS) {
    RemeshBoundaryInfo info;
    info.mIsBoundaryVertex.assign(rS.mNumPoints, 0);
    std::map<std::pair<std::int64_t, std::int64_t>, int> use_count;
    for (std::size_t f = 0; f < rS.mNumFaces; ++f) {
        const std::int64_t v[3] = {rS.mCorners[f * 3], rS.mCorners[f * 3 + 1],
                                   rS.mCorners[f * 3 + 2]};
        for (int k = 0; k < 3; ++k) {
            std::int64_t a = v[k], b = v[(k + 1) % 3];
            if (a > b)
                std::swap(a, b);
            ++use_count[{a, b}];
        }
    }
    for (const auto& [edge, count] : use_count) {
        if (count == 1) {
            info.mBoundaryEdges.push_back(edge);
            info.mIsBoundaryVertex[static_cast<std::size_t>(edge.first)] = 1;
            info.mIsBoundaryVertex[static_cast<std::size_t>(edge.second)] = 1;
        }
    }
    return info;
}

/// The vertex visitation order `remesh_seed_clusters` picks its next seed
/// from: every boundary vertex (ascending id) before every interior vertex
/// (ascending id) when `rIsBoundaryVertex` names any, else plain ascending
/// id -- reproducing today's order exactly. This is the "pinned" half of
/// boundary preservation: a boundary vertex is always available to become a
/// seed before an interior one, so a cluster's boundary segment is anchored
/// early rather than being absorbed piecemeal by whichever interior cluster
/// reaches it first during the general BFS growth that follows unchanged.
std::vector<std::int64_t> remesh_seed_order(std::size_t NumPoints,
                                            const std::vector<std::uint8_t>& rIsBoundaryVertex) {
    std::vector<std::int64_t> order;
    order.reserve(NumPoints);
    if (rIsBoundaryVertex.empty()) {
        for (std::size_t i = 0; i < NumPoints; ++i)
            order.push_back(static_cast<std::int64_t>(i));
        return order;
    }
    for (std::size_t i = 0; i < NumPoints; ++i)
        if (rIsBoundaryVertex[i])
            order.push_back(static_cast<std::int64_t>(i));
    for (std::size_t i = 0; i < NumPoints; ++i)
        if (!rIsBoundaryVertex[i])
            order.push_back(static_cast<std::int64_t>(i));
    return order;
}

/// Grows `NumClusters` seeds breadth-first from the next unassigned vertex in
/// `rOrder`, each absorbing neighbours until it holds its share of the
/// remaining area. RNG-free by construction, unlike pyacvd. `rOrder` is
/// plain ascending-id order unless boundary preservation reorders it (see
/// `remesh_seed_order`); this function itself has no boundary-specific logic
/// at all -- it only ever asks "what's the next free vertex in this order."
void remesh_seed_clusters(std::vector<std::int64_t>& rClusters, const detail::NodeAdjacency& rAdj,
                          const std::vector<double>& rArea,
                          const std::vector<std::int64_t>& rOrder, std::int64_t NumClusters,
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

        while (next_free < NumPoints &&
               rClusters[static_cast<std::size_t>(rOrder[next_free])] != -1)
            ++next_free;
        if (next_free >= NumPoints)
            break;
        const std::size_t seed_v = static_cast<std::size_t>(rOrder[next_free]);

        c_area += rArea[seed_v];
        rClusters[seed_v] = c;
        ring.clear();
        ring.push_back(static_cast<std::int64_t>(seed_v));

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

    switch (Metric) {
    case RemeshMetric::Isotropic:
        return iso_badness;
    case RemeshMetric::Quadric: {
        // Quadric mode ADDS the quadric flatness error to the isotropic
        // compactness term rather than using the quadric error alone --
        // measured (see remesh.hpp's doc comment): pure quadric-error
        // badness has no compactness pull at all, since a cluster confined
        // to a common plane can grow arbitrarily thin and snake-like along a
        // low-curvature direction without its quadric error rising, which
        // repeatedly produced disconnected, non-manifold clusters even after
        // every repair pass (17 of 150 clusters still split on a smooth
        // closed sphere in gtest). The isotropic term is the SAME
        // `-|sgamma|^2/srho` form `Isotropic` mode uses on its own -- both
        // terms are area-weighted squared-distance-like quantities
        // (dimensionally consistent, no free blend weight to tune), and
        // their sum keeps the compactness pull that prevents pathological
        // elongation while still rewarding flatness/feature alignment where
        // a real crease or corner makes the quadric term the deciding
        // factor.
        double new_quad[10];
        for (int k = 0; k < 10; ++k)
            new_quad[k] = pQuad[k] + sign * pItemQuad[k];
        double centroid[3];
        remesh_centroid_of(new_sgamma, new_srho, centroid);
        return iso_badness + remesh_quadric_badness(new_quad, centroid);
    }
    case RemeshMetric::Anisotropic: {
        // The per-vertex curvature tensor is SPD (never rank-deficient the
        // way a flat quadric is), so unlike `Quadric` this metric's own
        // error term already penalises spread in every direction -- just
        // anisotropically -- and is used PURE, not additively stabilised.
        // See remesh.hpp's doc comment and doc/remesh.md for the measurement
        // that justified this (mNumIsolatedClusters/mNumNonManifoldVertices
        // swept across mMaxAnisotropy on the fixtures below).
        double new_quad[10];
        for (int k = 0; k < 10; ++k)
            new_quad[k] = pQuad[k] + sign * pItemQuad[k];
        double centroid[3];
        remesh_centroid_of(new_sgamma, new_srho, centroid);
        return remesh_quadric_badness(new_quad, centroid);
    }
    }
    // Unreachable: every RemeshMetric enumerator is handled above, and
    // -Wswitch (no `default:`) makes a future enumerator a compile error
    // here rather than a silent quadric/anisotropic misclassification.
    return iso_badness;
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
    // The 10-double accumulator is shared by both Quadric and Anisotropic
    // (the latter packs a per-vertex curvature-tensor quadric into the exact
    // same layout, see remesh_vertex_metrics) -- Isotropic is the one metric
    // that never reads it.
    const bool needs_quad = Metric != RemeshMetric::Isotropic;
    if (needs_quad)
        st.mQuad.assign(n_clus * 10, 0.0);

    for (std::size_t i = 0; i < NumPoints; ++i) {
        const std::size_t c = static_cast<std::size_t>(rClusters[i]);
        st.mSrho[c] += rArea[i];
        st.mSgamma[c * 3] += rWeighted[i * 3];
        st.mSgamma[c * 3 + 1] += rWeighted[i * 3 + 1];
        st.mSgamma[c * 3 + 2] += rWeighted[i * 3 + 2];
        ++st.mCount[c];
        if (needs_quad)
            for (int k = 0; k < 10; ++k)
                st.mQuad[c * 10 + static_cast<std::size_t>(k)] +=
                    rVertexQuad[i * 10 + static_cast<std::size_t>(k)];
    }
    for (std::size_t c = 0; c < n_clus; ++c) {
        const double iso = -(st.mSgamma[c * 3] * st.mSgamma[c * 3] +
                             st.mSgamma[c * 3 + 1] * st.mSgamma[c * 3 + 1] +
                             st.mSgamma[c * 3 + 2] * st.mSgamma[c * 3 + 2]) /
                           st.mSrho[c];
        switch (Metric) {
        case RemeshMetric::Isotropic:
            st.mBadness[c] = iso;
            break;
        case RemeshMetric::Quadric: {
            double centroid[3];
            remesh_centroid_of(st.mSgamma.data() + c * 3, st.mSrho[c], centroid);
            st.mBadness[c] = iso + remesh_quadric_badness(st.mQuad.data() + c * 10, centroid);
            break;
        }
        case RemeshMetric::Anisotropic: {
            double centroid[3];
            remesh_centroid_of(st.mSgamma.data() + c * 3, st.mSrho[c], centroid);
            st.mBadness[c] = remesh_quadric_badness(st.mQuad.data() + c * 10, centroid);
            break;
        }
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
    switch (Metric) {
    case RemeshMetric::Isotropic:
        rSt.mBadness[uf] = -(rSt.mSgamma[uf * 3] * rSt.mSgamma[uf * 3] +
                             rSt.mSgamma[uf * 3 + 1] * rSt.mSgamma[uf * 3 + 1] +
                             rSt.mSgamma[uf * 3 + 2] * rSt.mSgamma[uf * 3 + 2]) /
                           rSt.mSrho[uf];
        rSt.mBadness[ut] = -(rSt.mSgamma[ut * 3] * rSt.mSgamma[ut * 3] +
                             rSt.mSgamma[ut * 3 + 1] * rSt.mSgamma[ut * 3 + 1] +
                             rSt.mSgamma[ut * 3 + 2] * rSt.mSgamma[ut * 3 + 2]) /
                           rSt.mSrho[ut];
        break;
    case RemeshMetric::Quadric:
    case RemeshMetric::Anisotropic: {
        for (int k = 0; k < 10; ++k) {
            rSt.mQuad[ut * 10 + static_cast<std::size_t>(k)] +=
                rVertexQuad[Item * 10 + static_cast<std::size_t>(k)];
            rSt.mQuad[uf * 10 + static_cast<std::size_t>(k)] -=
                rVertexQuad[Item * 10 + static_cast<std::size_t>(k)];
        }
        double cen_f[3], cen_t[3];
        remesh_centroid_of(rSt.mSgamma.data() + uf * 3, rSt.mSrho[uf], cen_f);
        remesh_centroid_of(rSt.mSgamma.data() + ut * 3, rSt.mSrho[ut], cen_t);
        // Quadric adds the isotropic compactness term back in (see
        // remesh_trial_badness); Anisotropic's own SPD tensor already
        // penalises spread in every direction and uses the quadric error
        // pure.
        double iso_f = 0.0, iso_t = 0.0;
        if (Metric == RemeshMetric::Quadric) {
            iso_f = -(rSt.mSgamma[uf * 3] * rSt.mSgamma[uf * 3] +
                     rSt.mSgamma[uf * 3 + 1] * rSt.mSgamma[uf * 3 + 1] +
                     rSt.mSgamma[uf * 3 + 2] * rSt.mSgamma[uf * 3 + 2]) /
                   rSt.mSrho[uf];
            iso_t = -(rSt.mSgamma[ut * 3] * rSt.mSgamma[ut * 3] +
                     rSt.mSgamma[ut * 3 + 1] * rSt.mSgamma[ut * 3 + 1] +
                     rSt.mSgamma[ut * 3 + 2] * rSt.mSgamma[ut * 3 + 2]) /
                   rSt.mSrho[ut];
        }
        rSt.mBadness[uf] = iso_f + remesh_quadric_badness(rSt.mQuad.data() + uf * 10, cen_f);
        rSt.mBadness[ut] = iso_t + remesh_quadric_badness(rSt.mQuad.data() + ut * 10, cen_t);
        break;
    }
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

            // Both Quadric and Anisotropic read the 10-double accumulator;
            // only Isotropic never touches it (see remesh_init_state).
            const bool needs_quad = Metric != RemeshMetric::Isotropic;
            const double* q_a = needs_quad ? rSt.mQuad.data() + uca * 10 : nullptr;
            const double* q_b = needs_quad ? rSt.mQuad.data() + ucb * 10 : nullptr;
            const double* item_q_b = needs_quad ? rVertexQuad.data() + pb * 10 : nullptr;
            const double* item_q_a = needs_quad ? rVertexQuad.data() + pa * 10 : nullptr;

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
    const bool needs_quad = Metric != RemeshMetric::Isotropic;
    if (needs_quad)
        quad.assign(n_clus * 10, 0.0);

    for (std::size_t i = 0; i < NumPoints; ++i) {
        const std::size_t c = static_cast<std::size_t>(rClusters[i]);
        srho[c] += rArea[i];
        for (int d = 0; d < 3; ++d)
            sgamma[c * 3 + d] += rWeighted[i * 3 + d];
        if (needs_quad)
            for (int k = 0; k < 10; ++k)
                quad[c * 10 + static_cast<std::size_t>(k)] +=
                    rVertexQuad[i * 10 + static_cast<std::size_t>(k)];
    }

    std::vector<double> out(n_clus * 3, 0.0);
    for (std::size_t c = 0; c < n_clus; ++c) {
        double centroid[3];
        remesh_centroid_of(sgamma.data() + c * 3, srho[c], centroid);
        switch (Metric) {
        case RemeshMetric::Isotropic:
            out[c * 3] = centroid[0];
            out[c * 3 + 1] = centroid[1];
            out[c * 3 + 2] = centroid[2];
            break;
        case RemeshMetric::Quadric:
        case RemeshMetric::Anisotropic: {
            // Both metrics place the dual vertex at their own quadric's
            // optimal point -- Anisotropic's quadric already encodes the
            // curvature-tensor shape, so no separate placement rule exists.
            double pt[3];
            detail::decim_quadric_optimal_point(quad.data() + c * 10, centroid, pt);
            out[c * 3] = pt[0];
            out[c * 3 + 1] = pt[1];
            out[c * 3 + 2] = pt[2];
            break;
        }
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

/// The deduplicated set of (sorted) cluster-id triples the current
/// clustering's dual would have as triangles -- topology only, no
/// positions/winding. Shared by `remesh_dual_faces` (which adds those on
/// top, and is only ever called once every vertex is resolved) and the
/// repair loop's non-manifold check (which may call this on a NOT-yet-fully
/// resolved cluster array -- `remesh_split_disconnected` unassigns vertices
/// back to `-1` and this may be read before the next `remesh_grow_null`
/// re-fills them, so a face touching any still-unassigned corner is skipped
/// outright: with no valid output vertex for that corner yet, it cannot
/// contribute to even a provisional dual).
std::vector<std::array<std::int64_t, 3>> remesh_dual_triangle_keys(
    const RemeshSurface& rS, const std::vector<std::int64_t>& rClusters) {
    std::vector<std::array<std::int64_t, 3>> keys;
    keys.reserve(rS.mNumFaces);
    for (std::size_t f = 0; f < rS.mNumFaces; ++f) {
        std::array<std::int64_t, 3> t = {
            rClusters[static_cast<std::size_t>(rS.mCorners[f * 3])],
            rClusters[static_cast<std::size_t>(rS.mCorners[f * 3 + 1])],
            rClusters[static_cast<std::size_t>(rS.mCorners[f * 3 + 2])]};
        if (t[0] < 0 || t[1] < 0 || t[2] < 0)
            continue;  // a corner has no cluster yet
        std::sort(t.begin(), t.end());
        if (t[0] == t[1] || t[1] == t[2])
            continue;  // the face does not join three distinct clusters
        keys.push_back(t);
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    return keys;
}

/// Builds the dual: one triangle per deduplicated key from
/// `remesh_dual_triangle_keys`, wound to agree with the mean normal of its
/// three clusters.
std::vector<std::int64_t> remesh_dual_faces(const RemeshSurface& rS,
                                            const std::vector<std::int64_t>& rClusters,
                                            const std::vector<double>& rPositions,
                                            const std::vector<double>& rNormals) {
    const std::vector<std::array<std::int64_t, 3>> keys =
        remesh_dual_triangle_keys(rS, rClusters);

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

/// One `line` cell per genuine input boundary edge whose two endpoints land
/// in different clusters, connecting those two clusters' representative
/// points -- deduplicated by unordered cluster-id pair, the same idiom
/// `remesh_dual_triangle_keys` already uses for the triangle dual. Without
/// this, a boundary-adjacent input triangle simply has no "third neighbour"
/// across the missing side and is silently dropped from the triangle dual,
/// leaving no coherent output boundary at all -- this is the "extra dual
/// elements... inserted along the boundary" half of boundary preservation.
std::vector<std::int64_t> remesh_dual_boundary_edges(
    const std::vector<std::pair<std::int64_t, std::int64_t>>& rBoundaryEdges,
    const std::vector<std::int64_t>& rClusters) {
    std::vector<std::pair<std::int64_t, std::int64_t>> keys;
    keys.reserve(rBoundaryEdges.size());
    for (const auto& [a, b] : rBoundaryEdges) {
        std::int64_t ca = rClusters[static_cast<std::size_t>(a)];
        std::int64_t cb = rClusters[static_cast<std::size_t>(b)];
        if (ca == cb)
            continue;
        if (ca > cb)
            std::swap(ca, cb);
        keys.push_back({ca, cb});
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

    std::vector<std::int64_t> out;
    out.reserve(keys.size() * 2);
    for (const auto& [ca, cb] : keys) {
        out.push_back(ca);
        out.push_back(cb);
    }
    return out;
}

/// Output (cluster) vertices whose incident dual-triangle fan is NOT a
/// single loop (interior vertex, closed fan) or open chain (boundary
/// vertex, open fan) -- a "bowtie." `rKeys` is `remesh_dual_triangle_keys`'
/// own output, so this works identically whether called on the FINAL
/// (post-renumber, contiguous) cluster space or the repair loop's CURRENT
/// (possibly-gapped) one -- it only ever asks about ids that actually
/// appear in `rKeys`.
///
/// For vertex `v`, each incident triangle `(v, a, b)` contributes the
/// "opposite edge" `(a, b)` (the edge not touching `v`); manifold means
/// those opposite edges, viewed as a graph on `v`'s neighbours, are a
/// single simple cycle (every node degree 2) or a single simple path
/// (exactly two nodes of degree 1, the open chain's ends, the rest degree
/// 2) -- checked as: no node of degree > 2, either 0 or 2 nodes of degree
/// 1, and one connected component spanning every opposite edge.
std::vector<std::uint8_t> remesh_detect_nonmanifold_vertices(
    const std::vector<std::array<std::int64_t, 3>>& rKeys, std::int64_t NumClusters) {
    const std::size_t n_clus = static_cast<std::size_t>(NumClusters);
    std::vector<std::vector<std::pair<std::int64_t, std::int64_t>>> incident(n_clus);
    for (const auto& t : rKeys) {
        for (int k = 0; k < 3; ++k) {
            const std::int64_t v = t[k], a = t[(k + 1) % 3], b = t[(k + 2) % 3];
            incident[static_cast<std::size_t>(v)].push_back({a, b});
        }
    }

    std::vector<std::uint8_t> bad(n_clus, 0);
    for (std::size_t v = 0; v < n_clus; ++v) {
        const auto& edges = incident[v];
        if (edges.size() < 2)
            continue;  // 0 or 1 incident dual triangle is trivially fine

        std::map<std::int64_t, std::vector<std::int64_t>> adj;
        for (const auto& [a, b] : edges) {
            adj[a].push_back(b);
            adj[b].push_back(a);
        }

        bool ok = true;
        int deg1_count = 0;
        for (const auto& [node, nbrs] : adj) {
            if (nbrs.size() > 2) {
                ok = false;
                break;
            }
            if (nbrs.size() == 1)
                ++deg1_count;
        }
        if (ok && deg1_count != 0 && deg1_count != 2)
            ok = false;
        if (ok) {
            std::set<std::int64_t> visited;
            std::vector<std::int64_t> stack = {edges[0].first};
            visited.insert(edges[0].first);
            while (!stack.empty()) {
                const std::int64_t cur = stack.back();
                stack.pop_back();
                for (const std::int64_t nb : adj.at(cur))
                    if (visited.insert(nb).second)
                        stack.push_back(nb);
            }
            if (visited.size() != adj.size())
                ok = false;
        }
        if (!ok)
            bad[v] = 1;
    }
    return bad;
}

}  // namespace

RemeshMetric remesh_metric_from_name(const std::string& rName) {
    if (rName == "isotropic")
        return RemeshMetric::Isotropic;
    if (rName == "quadric")
        return RemeshMetric::Quadric;
    if (rName == "anisotropic")
        return RemeshMetric::Anisotropic;
    throw std::invalid_argument("meshio++: remesh: unknown metric '" + rName +
                                "' (expected 'isotropic', 'quadric' or 'anisotropic')");
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
    if (rOptions.mMetric == RemeshMetric::Anisotropic && rOptions.mMaxAnisotropy < 1.0)
        throw std::invalid_argument(
            "meshio++: remesh: max_anisotropy must be at least 1.0, got " +
            std::to_string(rOptions.mMaxAnisotropy));
    if (rOptions.mMetric != RemeshMetric::Anisotropic &&
        rOptions.mMaxAnisotropy != kRemeshDefaultMaxAnisotropy)
        throw std::invalid_argument(
            "meshio++: remesh: max_anisotropy is only meaningful under "
            "metric=\"anisotropic\" (there is nothing for it to shape under '" +
            std::string(rOptions.mMetric == RemeshMetric::Isotropic ? "isotropic" : "quadric") +
            "')");

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

    // Adjacency is needed by curvature (1-ring) as well as the edge sweep and
    // seeding below, so it is built once, up front, with no dependency on the
    // item weights that follow.
    const detail::NodeAdjacency adj =
        detail::build_node_adjacency(work, surf.mNumPoints, detail::NodeAdjacencyKind::Edge);
    const std::vector<std::int64_t> edges = remesh_unique_edges(adj, surf.mNumPoints);

    std::vector<double> curvature;
    std::vector<double> curvature_tensors;  // Anisotropic only, 6 doubles/vertex
    if (rOptions.mGradation != 0.0 || metric == RemeshMetric::Anisotropic) {
        const std::vector<double> vertex_normals = remesh_vertex_normals(surf);
        curvature = remesh_vertex_curvature(
            surf, adj, vertex_normals,
            metric == RemeshMetric::Anisotropic ? &curvature_tensors : nullptr,
            rOptions.mMaxAnisotropy);
    }
    std::vector<double> area;
    std::vector<double> weighted;
    remesh_item_weights(surf, curvature, rOptions.mGradation, area, weighted);
    // Anisotropic packs its per-vertex curvature tensor into the accumulator
    // AFTER item weights exist (it needs rArea); Quadric's plain
    // face-quadric accumulation needs no weight at all.
    const std::vector<double> vertex_quad = metric == RemeshMetric::Anisotropic
                                                ? remesh_vertex_metrics(surf, curvature_tensors, area)
                                                : remesh_vertex_quadrics(surf, metric);

    const RemeshBoundaryInfo boundary =
        rOptions.mPreserveBoundary ? remesh_boundary_info(surf) : RemeshBoundaryInfo{};
    const std::vector<std::int64_t> seed_order =
        remesh_seed_order(surf.mNumPoints, boundary.mIsBoundaryVertex);

    std::vector<std::int64_t> clusters(surf.mNumPoints, -1);
    remesh_seed_clusters(clusters, adj, area, seed_order, rOptions.mNumClusters, surf.mNumPoints);
    if (remesh_grow_null(edges, clusters))
        for (std::size_t i = 0; i < surf.mNumPoints; ++i)
            if (clusters[i] == -1)
                clusters[i] = 0;

    RemeshState state = remesh_init_state(clusters, area, weighted, vertex_quad, metric,
                                          rOptions.mNumClusters, surf.mNumPoints);
    result.mNumIterations =
        remesh_minimize(edges, clusters, area, weighted, vertex_quad, metric, state,
                        rOptions.mNumClusters, rOptions.mMaxIterations);

    // Repair: a minimisation sweep can leave a cluster in two pieces, or the
    // dual can come out non-manifold at some vertex even with every cluster
    // itself connected -- two distinct causes, both folded into the same
    // split/re-grow/re-minimise loop and both reported separately at the end.
    // The non-manifold check runs on the CURRENT (possibly-gapped,
    // pre-renumber) cluster id space via remesh_dual_triangle_keys, which
    // makes no assumption of contiguity.
    auto unassign_bad_clusters = [&](const std::vector<std::uint8_t>& rBad) {
        for (std::size_t i = 0; i < surf.mNumPoints; ++i)
            if (clusters[i] >= 0 && rBad[static_cast<std::size_t>(clusters[i])])
                clusters[i] = -1;
    };

    std::int64_t disconnected = remesh_split_disconnected(adj, clusters, rOptions.mNumClusters,
                                                          surf.mNumPoints);
    std::vector<std::uint8_t> bad_vertices = remesh_detect_nonmanifold_vertices(
        remesh_dual_triangle_keys(surf, clusters), rOptions.mNumClusters);
    std::int64_t nonmanifold =
        std::count(bad_vertices.begin(), bad_vertices.end(), std::uint8_t{1});
    if (nonmanifold > 0)
        unassign_bad_clusters(bad_vertices);

    for (int pass = 0; pass < rOptions.mMaxRepairPasses && (disconnected > 0 || nonmanifold > 0);
        ++pass) {
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
        bad_vertices = remesh_detect_nonmanifold_vertices(
            remesh_dual_triangle_keys(surf, clusters), rOptions.mNumClusters);
        nonmanifold = std::count(bad_vertices.begin(), bad_vertices.end(), std::uint8_t{1});
        if (nonmanifold > 0)
            unassign_bad_clusters(bad_vertices);
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

    // The honest, final, vertex-level manifoldness count: whatever repair
    // could not fix, checked once more on the true final (renumbered)
    // clustering -- distinct in kind from mNumIsolatedClusters, which the
    // loop above already tracks at cluster granularity.
    const std::vector<std::uint8_t> final_bad = remesh_detect_nonmanifold_vertices(
        remesh_dual_triangle_keys(surf, clusters), n_clus);
    result.mNumNonManifoldVertices =
        std::count(final_bad.begin(), final_bad.end(), std::uint8_t{1});

    const std::vector<double> positions = remesh_final_positions(
        clusters, area, weighted, vertex_quad, metric, n_clus, surf.mNumPoints);
    const std::vector<double> normals = remesh_cluster_normals(surf, clusters, n_clus);
    const std::vector<std::int64_t> faces = remesh_dual_faces(surf, clusters, positions, normals);
    const std::vector<std::int64_t> boundary_lines =
        rOptions.mPreserveBoundary
            ? remesh_dual_boundary_edges(boundary.mBoundaryEdges, clusters)
            : std::vector<std::int64_t>{};

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
    if (!boundary_lines.empty()) {
        const std::size_t nl = boundary_lines.size() / 2;
        NDArray lconn = NDArray::Uninit(DType::Int64, {nl, 2});
        std::int64_t* dst = lconn.As<std::int64_t>();
        for (std::size_t i = 0; i < boundary_lines.size(); ++i)
            dst[i] = boundary_lines[i];
        out.AddCellBlock(cell_type_name(CellType::Line), std::move(lconn));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, rMesh.FieldData(name));

    return result;
}

}  // namespace meshioplusplus
