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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#include "meshioplusplus/operations/transfer_operator.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/parallel.hpp"

#include "../detail/point_knn.hpp"
#include "../detail/point_locator.hpp"
#include "../detail/sparse_solve.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

std::string to_lower_copy(const std::string& rS) {
    std::string out = rS;
    for (char& c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

std::vector<Vec3> read_xyz_vec3(const Mesh& rMesh) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& pts = rMesh.Points();
    std::vector<Vec3> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = {detail::read_double(pts, i * pdim + 0),
                  pdim > 1 ? detail::read_double(pts, i * pdim + 1) : 0.0,
                  pdim > 2 ? detail::read_double(pts, i * pdim + 2) : 0.0};
    }
    return out;
}

// Recursively resolve a fine node's coarse weights from `refine:entity`.
// Returns false if a node cannot be resolved strictly to [0, n_coarse).
bool to_resolve_nested(std::int64_t Id, std::size_t NCoarse, const std::int64_t* pEnt,
                       std::size_t NFine, double Scale,
                       std::vector<std::pair<std::size_t, double>>& rWeights, int Depth) {
    if (Depth > 16 || Id < 0 || static_cast<std::size_t>(Id) >= NFine)
        return false;
    if (static_cast<std::size_t>(Id) < NCoarse) {
        rWeights.emplace_back(static_cast<std::size_t>(Id), Scale);
        return true;
    }
    const std::int64_t* e = pEnt + static_cast<std::size_t>(Id) * 4;
    if (e[0] < 0 && e[1] < 0 && e[2] >= 0 && e[3] >= 0) {
        return to_resolve_nested(e[2], NCoarse, pEnt, NFine, Scale * 0.5, rWeights, Depth + 1) &&
               to_resolve_nested(e[3], NCoarse, pEnt, NFine, Scale * 0.5, rWeights, Depth + 1);
    }
    if (e[0] >= 0 && e[1] >= 0 && e[2] >= 0 && e[3] >= 0) {
        for (int k = 0; k < 4; ++k) {
            if (!to_resolve_nested(e[k], NCoarse, pEnt, NFine, Scale * 0.25, rWeights, Depth + 1))
                return false;
        }
        return true;
    }
    return false;
}

bool to_is_surface_in_3d(const Mesh& rMesh) {
    bool has_2d = false;
    bool has_3d = false;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.NumCells() == 0)
            continue;
        const int dim = cell_type_dimension(cell_type_from_name(cb.Type()));
        if (dim == 3)
            has_3d = true;
        else if (dim == 2)
            has_2d = true;
    }
    if (has_3d || !has_2d)
        return false;
    // Check if z varies across points; if planar z == const, 2-D PointLocator works directly.
    if (rMesh.PointDim() < 3 || rMesh.NumPoints() == 0)
        return false;
    const NDArray& pts = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    const double z0 = detail::read_double(pts, 2);
    for (std::size_t i = 1; i < rMesh.NumPoints(); ++i) {
        if (std::abs(detail::read_double(pts, i * pdim + 2) - z0) > 1e-12)
            return true;
    }
    return false;
}

}  // namespace

TransferOperatorMode transfer_operator_mode_from_name(const std::string& rName) {
    const std::string low = to_lower_copy(rName);
    if (low == "auto")
        return TransferOperatorMode::Auto;
    if (low == "nested")
        return TransferOperatorMode::Nested;
    if (low == "barycentric")
        return TransferOperatorMode::Barycentric;
    throw std::invalid_argument("transfer_operator: unknown mode '" + rName +
                                "' (expected 'auto', 'nested', or 'barycentric')");
}

MeshOperatorResult transfer_operator(const Mesh& rCoarse, const Mesh& rFine,
                                     const TransferOperatorOptions& rOptions) {
    if (!(rOptions.mTolerance >= 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument("transfer_operator: tolerance must be finite and >= 0");
    const std::size_t n_coarse = rCoarse.NumPoints();
    const std::size_t n_fine = rFine.NumPoints();
    if (n_coarse == 0 || n_fine == 0)
        throw std::invalid_argument("transfer_operator: meshes must be non-empty");

    const std::vector<Vec3> coarse_xyz = read_xyz_vec3(rCoarse);
    const std::vector<Vec3> fine_xyz = read_xyz_vec3(rFine);

    // Per-fine-point list of (coarse_col, weight).
    std::vector<std::vector<std::pair<std::size_t, double>>> rows_w(n_fine);
    std::vector<std::uint8_t> resolved(n_fine, 0);

    // 1. Try nested refine:entity path if Auto or Nested.
    bool has_entity = false;
    if (rOptions.mMode != TransferOperatorMode::Barycentric &&
        rFine.HasPointData(kRefineEntityName)) {
        const NDArray& ent = rFine.PointData(kRefineEntityName);
        if (ent.Dtype() == DType::Int64 && ent.Shape().size() == 2 && ent.Shape()[0] == n_fine &&
            ent.Shape()[1] == 4 && n_fine >= n_coarse) {
            // Verify leading coarse points match coordinates within tolerance.
            bool coords_match = true;
            const double tol2 = std::max(rOptions.mTolerance, 1e-9) *
                                std::max(rOptions.mTolerance, 1e-9);
            for (std::size_t i = 0; i < n_coarse; ++i) {
                const Vec3 d = detail::vec3_sub(fine_xyz[i], coarse_xyz[i]);
                if (detail::vec3_norm_sq(d) > tol2) {
                    coords_match = false;
                    break;
                }
            }
            if (coords_match) {
                has_entity = true;
                const std::int64_t* p_ent = ent.As<std::int64_t>();
                for (std::size_t i = 0; i < n_fine; ++i) {
                    std::vector<std::pair<std::size_t, double>> w;
                    if (to_resolve_nested(static_cast<std::int64_t>(i), n_coarse, p_ent, n_fine,
                                          1.0, w, 0)) {
                        rows_w[i] = std::move(w);
                        resolved[i] = 1;
                    }
                }
            }
        }
    }

    if (rOptions.mMode == TransferOperatorMode::Nested && !has_entity)
        throw std::invalid_argument(
            "transfer_operator: mode='nested' requires a compatible 'refine:entity' array on fine");

    // 2. Locate any remaining unresolved points geometrically.
    bool any_unresolved = false;
    for (std::uint8_t r : resolved)
        if (!r) {
            any_unresolved = true;
            break;
        }

    std::int64_t num_unlocated = 0;
    if (any_unresolved) {
        if (to_is_surface_in_3d(rCoarse)) {
            // 3-D curved surface: project onto nearest triangle via DistanceQuery.
            const detail::TriangleSoup soup = detail::build_triangle_soup(rCoarse, "");
            const detail::DistanceQuery dq =
                detail::build_distance_query(soup, SurfaceDistanceOptions{});
            std::vector<std::size_t> todo;
            std::vector<Vec3> qpts;
            for (std::size_t i = 0; i < n_fine; ++i) {
                if (!resolved[i]) {
                    todo.push_back(i);
                    qpts.push_back(fine_xyz[i]);
                }
            }
            const std::vector<detail::SurfaceProjection> hits =
                detail::query_surface_projections(dq, qpts);
            parallel_for(todo.size(), [&](std::size_t idx) {
                const std::size_t fi = todo[idx];
                const detail::SurfaceProjection& hit = hits[idx];
                if (!hit.mFound || hit.mTriangle < 0)
                    return;
                const auto& tri = soup.mVertices[static_cast<std::size_t>(hit.mTriangle)];
                const Vec3 p0 = soup.mPoints[static_cast<std::size_t>(tri[0])];
                const Vec3 p1 = soup.mPoints[static_cast<std::size_t>(tri[1])];
                const Vec3 p2 = soup.mPoints[static_cast<std::size_t>(tri[2])];
                const Vec3& q = hit.mPoint;
                const Vec3 v0 = detail::vec3_sub(p1, p0);
                const Vec3 v1 = detail::vec3_sub(p2, p0);
                const Vec3 v2 = detail::vec3_sub(q, p0);
                const double d00 = detail::vec3_dot(v0, v0);
                const double d01 = detail::vec3_dot(v0, v1);
                const double d11 = detail::vec3_dot(v1, v1);
                const double d20 = detail::vec3_dot(v2, v0);
                const double d21 = detail::vec3_dot(v2, v1);
                const double denom = d00 * d11 - d01 * d01;
                double b1 = 0.0, b2 = 0.0;
                if (std::abs(denom) > 1e-30) {
                    b1 = (d11 * d20 - d01 * d21) / denom;
                    b2 = (d00 * d21 - d01 * d20) / denom;
                }
                const double b0 = 1.0 - b1 - b2;
                rows_w[fi] = {{static_cast<std::size_t>(tri[0]), b0},
                              {static_cast<std::size_t>(tri[1]), b1},
                              {static_cast<std::size_t>(tri[2]), b2}};
                resolved[fi] = 1;
            });
        } else {
            // Volume or 2-D planar mesh: simplexify + PointLocator.
            ConvertCellsOptions copt;
            copt.mMode = ConvertCellsMode::Simplexify;
            const auto simp_res = convert_cells(rCoarse, copt);
            const detail::PointLocator locator(simp_res.mMesh, rOptions.mTolerance, "transfer_operator");
            const std::size_t k = locator.NodesPerSimplex();
            const std::vector<std::int64_t>& conn = locator.Connectivity();

            parallel_for(n_fine, [&](std::size_t i) {
                if (resolved[i])
                    return;
                double w[4] = {0.0, 0.0, 0.0, 0.0};
                const std::int64_t sid = locator.Locate(fine_xyz[i], w);
                if (sid >= 0) {
                    rows_w[i].reserve(k);
                    for (std::size_t a = 0; a < k; ++a) {
                        rows_w[i].emplace_back(
                            static_cast<std::size_t>(conn[static_cast<std::size_t>(sid) * k + a]),
                            w[a]);
                    }
                    resolved[i] = 1;
                }
            });
        }
    }

    // 3. Optional nearest-vertex fallback for any still-unlocated points.
    std::vector<std::size_t> missing;
    for (std::size_t i = 0; i < n_fine; ++i)
        if (!resolved[i])
            missing.push_back(i);

    if (!missing.empty()) {
        if (rOptions.mFallbackNearest) {
            Vec3 lo = coarse_xyz[0], hi = coarse_xyz[0];
            for (const Vec3& p : coarse_xyz)
                for (int d = 0; d < 3; ++d) {
                    lo[d] = std::min(lo[d], p[d]);
                    hi[d] = std::max(hi[d], p[d]);
                }
            const double diag = detail::vec3_norm(detail::vec3_sub(hi, lo));
            const double h =
                diag > 0.0 ? diag / std::cbrt(static_cast<double>(n_coarse) + 1.0) : 1.0;
            const detail::PointKnn knn(coarse_xyz, std::max(h, 1e-12));
            parallel_for(missing.size(), [&](std::size_t idx) {
                const std::size_t fi = missing[idx];
                const std::int64_t nid = knn.Nearest(fine_xyz[fi]);
                if (nid >= 0) {
                    rows_w[fi] = {{static_cast<std::size_t>(nid), 1.0}};
                    resolved[fi] = 1;
                }
            });
        } else {
            num_unlocated = static_cast<std::int64_t>(missing.size());
        }
    }

    // 4. Assemble CSR (prolongation P or transposed restriction R = P^T).
    std::vector<detail::SparseTriplet> triplets;
    triplets.reserve(n_fine * 4);
    for (std::size_t i = 0; i < n_fine; ++i) {
        for (const auto& [col, val] : rows_w[i]) {
            if (std::abs(val) <= 1e-15)
                continue;
            if (rOptions.mTranspose)
                triplets.push_back({col, i, val});
            else
                triplets.push_back({i, col, val});
        }
    }

    const std::size_t out_rows = rOptions.mTranspose ? n_coarse : n_fine;
    const std::size_t out_cols = rOptions.mTranspose ? n_fine : n_coarse;
    const detail::SparseCsr csr =
        detail::csr_from_triplets(out_rows, out_cols, std::move(triplets));

    MeshOperatorResult out;
    out.mRows = static_cast<std::int64_t>(csr.mRows);
    out.mCols = static_cast<std::int64_t>(csr.mCols);
    out.mIndptr.assign(csr.mStart.begin(), csr.mStart.end());
    out.mIndices.assign(csr.mCol.begin(), csr.mCol.end());
    out.mValues = csr.mValue;
    out.mNumDegenerate = num_unlocated;
    return out;
}

}  // namespace meshioplusplus
