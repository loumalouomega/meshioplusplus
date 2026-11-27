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

#include "meshioplusplus/operations/denoise_surface.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

#include "../detail/crease_edges.hpp"
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

std::string dn_lower(const std::string& rS) {
    std::string out = rS;
    for (char& c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

void dn_face_geometry(const detail::SurfaceRings& rRings, const std::vector<Vec3>& rXyz,
                      std::vector<Vec3>& rCentroids, std::vector<Vec3>& rNormals,
                      std::vector<double>& rAreas) {
    const std::size_t nf = rRings.NumFaces();
    rCentroids.resize(nf);
    rNormals.resize(nf);
    rAreas.resize(nf);
    parallel_for(nf, [&](std::size_t f) {
        const std::int64_t s = rRings.mStart[f];
        const std::int64_t e = rRings.mStart[f + 1];
        const std::size_t k = static_cast<std::size_t>(e - s);
        Vec3 c{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < k; ++i) {
            const Vec3& p = rXyz[static_cast<std::size_t>(rRings.mNodes[static_cast<std::size_t>(s) + i])];
            c[0] += p[0];
            c[1] += p[1];
            c[2] += p[2];
        }
        const double inv_k = k > 0 ? 1.0 / static_cast<double>(k) : 0.0;
        c[0] *= inv_k;
        c[1] *= inv_k;
        c[2] *= inv_k;
        rCentroids[f] = c;

        // Newell normal / area vector.
        Vec3 n{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < k; ++i) {
            const Vec3& p0 = rXyz[static_cast<std::size_t>(
                rRings.mNodes[static_cast<std::size_t>(s) + i])];
            const Vec3& p1 = rXyz[static_cast<std::size_t>(
                rRings.mNodes[static_cast<std::size_t>(s) + ((i + 1) % k)])];
            const Vec3 d0 = detail::vec3_sub(p0, c);
            const Vec3 d1 = detail::vec3_sub(p1, c);
            const Vec3 cr = detail::vec3_cross(d0, d1);
            n[0] += cr[0];
            n[1] += cr[1];
            n[2] += cr[2];
        }
        const double len = detail::vec3_norm(n);
        rAreas[f] = 0.5 * len;
        if (len > 1e-30)
            rNormals[f] = detail::vec3_scale(n, 1.0 / len);
        else
            rNormals[f] = {0.0, 0.0, 0.0};
    });
}

Vec3 dn_single_face_normal(const detail::SurfaceRings& rRings, std::size_t f,
                           const std::vector<Vec3>& rXyz, std::size_t OverrideV,
                           const Vec3& rOverridePos) {
    const std::int64_t s = rRings.mStart[f];
    const std::int64_t e = rRings.mStart[f + 1];
    const std::size_t k = static_cast<std::size_t>(e - s);
    auto get_p = [&](std::size_t i) -> Vec3 {
        const auto vid =
            static_cast<std::size_t>(rRings.mNodes[static_cast<std::size_t>(s) + i]);
        return vid == OverrideV ? rOverridePos : rXyz[vid];
    };
    const Vec3 p0 = get_p(0);
    Vec3 n{0.0, 0.0, 0.0};
    for (std::size_t i = 1; i + 1 < k; ++i) {
        const Vec3 cr =
            detail::vec3_cross(detail::vec3_sub(get_p(i), p0), detail::vec3_sub(get_p(i + 1), p0));
        n[0] += cr[0];
        n[1] += cr[1];
        n[2] += cr[2];
    }
    return n;
}

}  // namespace

DenoiseMethod denoise_method_from_name(const std::string& rName) {
    const std::string low = dn_lower(rName);
    if (low == "bilateral")
        return DenoiseMethod::Bilateral;
    if (low == "guided")
        return DenoiseMethod::Guided;
    throw std::invalid_argument("denoise_surface: unknown method '" + rName +
                                "' (expected 'bilateral' or 'guided')");
}

DenoiseSurfaceResult denoise_surface(const Mesh& rMesh, const DenoiseSurfaceOptions& rOptions) {
    if (rOptions.mNormalIterations < 0 || rOptions.mVertexIterations < 0)
        throw std::invalid_argument("denoise_surface: iteration counts must be >= 0");
    if (!(rOptions.mSigmaSpatial >= 0.0) || !std::isfinite(rOptions.mSigmaSpatial))
        throw std::invalid_argument("denoise_surface: sigma_spatial must be >= 0 and finite");
    if (!(rOptions.mSigmaNormal > 0.0) || !std::isfinite(rOptions.mSigmaNormal))
        throw std::invalid_argument("denoise_surface: sigma_normal must be > 0 and finite");

    const std::size_t npts = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    if (npts == 0)
        throw std::invalid_argument("denoise_surface: mesh has no points");

    // Resolve optional Cell region selection.
    const auto bases = detail::block_bases(rMesh);
    const std::size_t total_cells = bases.empty() ? 0 : static_cast<std::size_t>(bases.back());
    std::vector<std::uint8_t> selected(total_cells, rOptions.mRegion.empty() ? 1 : 0);
    if (!rOptions.mRegion.empty()) {
        const std::size_t ridx = rMesh.FindRegion(rOptions.mRegion, RegionKind::Cell);
        if (ridx == Mesh::npos)
            throw std::invalid_argument("denoise_surface: unknown Cell region '" + rOptions.mRegion +
                                        "'");
        const Region& reg = rMesh.Region(ridx);
        const std::int64_t* ent = reg.Entries();
        for (std::size_t i = 0; i < reg.NumEntries(); ++i)
            if (ent[i] >= 0 && static_cast<std::size_t>(ent[i]) < total_cells)
                selected[static_cast<std::size_t>(ent[i])] = 1;
    }

    const detail::SurfaceRings rings = detail::surface_rings(rMesh, selected, "denoise_surface");
    const std::size_t nf = rings.NumFaces();
    if (nf == 0)
        throw std::invalid_argument("denoise_surface: mesh has no surface faces");

    // Read coordinates.
    const NDArray& pts_in = rMesh.Points();
    std::vector<Vec3> xyz(npts);
    for (std::size_t i = 0; i < npts; ++i) {
        xyz[i] = {detail::read_double(pts_in, i * pdim + 0),
                  pdim > 1 ? detail::read_double(pts_in, i * pdim + 1) : 0.0,
                  pdim > 2 ? detail::read_double(pts_in, i * pdim + 2) : 0.0};
    }
    const std::vector<Vec3> xyz_orig = xyz;

    // Build vertex -> incident faces CSR.
    std::vector<std::vector<std::size_t>> v_faces(npts);
    for (std::size_t f = 0; f < nf; ++f) {
        for (std::int64_t p = rings.mStart[f]; p < rings.mStart[f + 1]; ++p) {
            v_faces[static_cast<std::size_t>(rings.mNodes[static_cast<std::size_t>(p)])].push_back(
                f);
        }
    }

    // Build face -> vertex-adjacent faces (sorted, deduplicated, excluding self)
    // and identify boundary vertices via crease_edge_uses.
    std::vector<std::uint8_t> is_boundary(npts, 0);
    std::vector<std::vector<std::size_t>> edge_adj(nf);
    {
        const std::vector<detail::CreaseEdgeUse> uses =
            detail::crease_edge_uses(rings.mStart, rings.mNodes);
        std::size_t i = 0;
        while (i < uses.size()) {
            std::size_t j = i + 1;
            while (j < uses.size() && uses[j].mLo == uses[i].mLo && uses[j].mHi == uses[i].mHi)
                ++j;
            if (j - i == 1) {
                is_boundary[static_cast<std::size_t>(uses[i].mLo)] = 1;
                is_boundary[static_cast<std::size_t>(uses[i].mHi)] = 1;
            } else if (j - i == 2) {
                const auto f0 = static_cast<std::size_t>(uses[i].mFace);
                const auto f1 = static_cast<std::size_t>(uses[i + 1].mFace);
                edge_adj[f0].push_back(f1);
                edge_adj[f1].push_back(f0);
            }
            i = j;
        }
    }

    std::vector<std::vector<std::size_t>> face_nbrs(nf);
    parallel_for(nf, [&](std::size_t f) {
        auto& nbrs = face_nbrs[f];
        for (std::int64_t p = rings.mStart[f]; p < rings.mStart[f + 1]; ++p) {
            const auto v = static_cast<std::size_t>(rings.mNodes[static_cast<std::size_t>(p)]);
            for (std::size_t g : v_faces[v])
                if (g != f)
                    nbrs.push_back(g);
        }
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    });

    // Initial face geometry.
    std::vector<Vec3> centroids, normals;
    std::vector<double> areas;
    dn_face_geometry(rings, xyz, centroids, normals, areas);

    // Determine sigma_spatial automatically if 0.
    double sigma_c = rOptions.mSigmaSpatial;
    if (sigma_c <= 0.0) {
        double sum_d = 0.0;
        std::size_t cnt = 0;
        for (std::size_t f = 0; f < nf; ++f) {
            for (std::size_t g : edge_adj[f]) {
                if (g > f) {
                    sum_d += detail::vec3_norm(detail::vec3_sub(centroids[f], centroids[g]));
                    ++cnt;
                }
            }
        }
        sigma_c = (cnt > 0 && sum_d > 0.0) ? (sum_d / static_cast<double>(cnt)) : 1.0;
    }
    const double inv_2sc2 = 1.0 / (2.0 * sigma_c * sigma_c);
    const double inv_2ss2 = 1.0 / (2.0 * rOptions.mSigmaNormal * rOptions.mSigmaNormal);

    // Stage 1: Iterative Jacobi normal filtering.
    std::vector<Vec3> next_normals(nf);
    std::vector<Vec3> guide_normals;
    if (rOptions.mMethod == DenoiseMethod::Guided)
        guide_normals.resize(nf);

    for (std::int64_t iter = 0; iter < rOptions.mNormalIterations; ++iter) {
        if (rOptions.mMethod == DenoiseMethod::Guided) {
            // Zhang et al. (2015): compute per-face patch consistency C(P_k) and guidance normal.
            std::vector<double> patch_cons(nf, 0.0);
            std::vector<Vec3> patch_avg(nf, {0.0, 0.0, 0.0});
            parallel_for(nf, [&](std::size_t f) {
                const auto& nbrs = face_nbrs[f];
                Vec3 avg = detail::vec3_scale(normals[f], areas[f]);
                double max_diff = 0.0;
                for (std::size_t g : nbrs) {
                    avg = detail::vec3_add(avg, detail::vec3_scale(normals[g], areas[g]));
                    const double d = detail::vec3_norm(detail::vec3_sub(normals[f], normals[g]));
                    max_diff = std::max(max_diff, d);
                }
                const double len = detail::vec3_norm(avg);
                patch_avg[f] = len > 1e-30 ? detail::vec3_scale(avg, 1.0 / len) : normals[f];
                // Consistency = max normal spread in patch * max adjacent edge salience.
                double max_sal = 0.0;
                for (std::size_t g : edge_adj[f])
                    max_sal =
                        std::max(max_sal, detail::vec3_norm(detail::vec3_sub(normals[f], normals[g])));
                patch_cons[f] = max_diff * (max_sal + 1e-6);
            });
            parallel_for(nf, [&](std::size_t f) {
                std::size_t best = f;
                double best_c = patch_cons[f];
                for (std::size_t g : face_nbrs[f]) {
                    if (patch_cons[g] < best_c) {
                        best_c = patch_cons[g];
                        best = g;
                    }
                }
                guide_normals[f] = patch_avg[best];
            });
        }

        const std::vector<Vec3>& ref_normals =
            (rOptions.mMethod == DenoiseMethod::Guided) ? guide_normals : normals;

        parallel_for(nf, [&](std::size_t f) {
            Vec3 acc = detail::vec3_scale(normals[f], areas[f]);
            const Vec3& ci = centroids[f];
            const Vec3& ri = ref_normals[f];
            for (std::size_t g : face_nbrs[f]) {
                const double dc2 = detail::vec3_norm_sq(detail::vec3_sub(ci, centroids[g]));
                const double dn2 = detail::vec3_norm_sq(detail::vec3_sub(ri, ref_normals[g]));
                const double w = areas[g] * std::exp(-dc2 * inv_2sc2 - dn2 * inv_2ss2);
                acc[0] += w * normals[g][0];
                acc[1] += w * normals[g][1];
                acc[2] += w * normals[g][2];
            }
            const double len = detail::vec3_norm(acc);
            next_normals[f] = len > 1e-30 ? detail::vec3_scale(acc, 1.0 / len) : normals[f];
        });
        normals.swap(next_normals);
    }

    // Stage 2: Vertex position update (Sun et al. 2007) with optional fold-over guard.
    std::int64_t skipped_foldover = 0;
    std::vector<Vec3> next_xyz(npts);
    std::vector<std::uint8_t> skip_flags(npts, 0);

    for (std::int64_t iter = 0; iter < rOptions.mVertexIterations; ++iter) {
        // Recompute centroids from current vertex positions.
        parallel_for(nf, [&](std::size_t f) {
            const std::int64_t s = rings.mStart[f];
            const std::int64_t e = rings.mStart[f + 1];
            const std::size_t k = static_cast<std::size_t>(e - s);
            Vec3 c{0.0, 0.0, 0.0};
            for (std::size_t i = 0; i < k; ++i) {
                const Vec3& p =
                    xyz[static_cast<std::size_t>(rings.mNodes[static_cast<std::size_t>(s) + i])];
                c[0] += p[0];
                c[1] += p[1];
                c[2] += p[2];
            }
            const double inv_k = k > 0 ? 1.0 / static_cast<double>(k) : 0.0;
            centroids[f] = {c[0] * inv_k, c[1] * inv_k, c[2] * inv_k};
        });

        std::fill(skip_flags.begin(), skip_flags.end(), 0);
        parallel_for(npts, [&](std::size_t v) {
            next_xyz[v] = xyz[v];
            if (v_faces[v].empty())
                return;
            if (rOptions.mFixBoundary && is_boundary[v])
                return;
            Vec3 delta{0.0, 0.0, 0.0};
            for (std::size_t f : v_faces[v]) {
                const double proj =
                    detail::vec3_dot(normals[f], detail::vec3_sub(centroids[f], xyz[v]));
                delta[0] += normals[f][0] * proj;
                delta[1] += normals[f][1] * proj;
                delta[2] += normals[f][2] * proj;
            }
            const double inv_deg = 1.0 / static_cast<double>(v_faces[v].size());
            const Vec3 cand{xyz[v][0] + delta[0] * inv_deg, xyz[v][1] + delta[1] * inv_deg,
                            xyz[v][2] + delta[2] * inv_deg};
            if (rOptions.mGuardFoldover) {
                for (std::size_t f : v_faces[v]) {
                    const Vec3 nf_new = dn_single_face_normal(rings, f, xyz, v, cand);
                    if (detail::vec3_dot(nf_new, normals[f]) <= 0.0) {
                        skip_flags[v] = 1;
                        return;
                    }
                }
            }
            next_xyz[v] = cand;
        });

        for (std::uint8_t fl : skip_flags)
            skipped_foldover += fl;
        xyz.swap(next_xyz);
    }

    // Write back to cloned mesh.
    DenoiseSurfaceResult result;
    result.mMesh = detail::clone_mesh(rMesh);
    result.mSigmaSpatialUsed = sigma_c;
    result.mNumSkippedFoldover = skipped_foldover;

    NDArray out_pts(DType::Float64, {npts, pdim});
    double* po = out_pts.As<double>();
    for (std::size_t i = 0; i < npts; ++i) {
        const double d = detail::vec3_norm(detail::vec3_sub(xyz[i], xyz_orig[i]));
        if (d > 0.0) {
            ++result.mNumNodesMoved;
            result.mMaxDisplacement = std::max(result.mMaxDisplacement, d);
        }
        po[i * pdim + 0] = xyz[i][0];
        if (pdim > 1)
            po[i * pdim + 1] = xyz[i][1];
        if (pdim > 2)
            po[i * pdim + 2] = xyz[i][2];
    }
    result.mMesh.AssignPoints(std::move(out_pts));
    return result;
}

}  // namespace meshioplusplus
