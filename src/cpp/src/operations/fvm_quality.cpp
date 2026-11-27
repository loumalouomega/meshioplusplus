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

#include "meshioplusplus/operations/fvm_quality.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/parallel.hpp"

#include "../detail/cell_centroids.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

QualityMetricSummary fvm_summarize(const std::vector<double>& rValues) {
    QualityMetricSummary s{};
    s.mMin = std::numeric_limits<double>::infinity();
    s.mMax = -std::numeric_limits<double>::infinity();
    s.mMean = 0.0;
    s.mCount = 0;
    s.mHistogram.fill(0);

    double sum = 0.0;
    for (double v : rValues) {
        if (!std::isfinite(v))
            continue;
        s.mMin = std::min(s.mMin, v);
        s.mMax = std::max(s.mMax, v);
        sum += v;
        ++s.mCount;
    }
    if (s.mCount == 0) {
        s.mMin = std::numeric_limits<double>::quiet_NaN();
        s.mMax = std::numeric_limits<double>::quiet_NaN();
        s.mMean = std::numeric_limits<double>::quiet_NaN();
        return s;
    }
    s.mMean = sum / static_cast<double>(s.mCount);
    const double span = s.mMax - s.mMin;
    for (double v : rValues) {
        if (!std::isfinite(v))
            continue;
        int bin = 0;
        if (span > 0.0) {
            bin = static_cast<int>(std::floor((v - s.mMin) / span * QualityMetricSummary::K));
            bin = std::clamp(bin, 0, QualityMetricSummary::K - 1);
        }
        ++s.mHistogram[static_cast<std::size_t>(bin)];
    }
    return s;
}

}  // namespace

FvmQualityReport compute_fvm_quality(const Mesh& rMesh) {
    FvmQualityReport report;
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t total_cells = bases.empty() ? 0 : static_cast<std::size_t>(bases.back());

    // Allocate NaN-initialized per-block output cell arrays for all 4 metrics.
    const double qnan = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> cell_nonortho(nblocks);
    std::vector<std::vector<double>> cell_skew(nblocks);
    std::vector<std::vector<double>> cell_vratio(nblocks);
    std::vector<std::vector<double>> cell_pyr(nblocks);
    for (std::size_t b = 0; b < nblocks; ++b) {
        const std::size_t nc = rMesh.Cells(b).NumCells();
        cell_nonortho[b].assign(nc, qnan);
        cell_skew[b].assign(nc, qnan);
        cell_vratio[b].assign(nc, qnan);
        cell_pyr[b].assign(nc, qnan);
    }

    const detail::GlobalFaces gf = detail::build_global_faces(rMesh);
    report.mNumCells = static_cast<std::int64_t>(gf.NumCells());
    report.mNumFaces = static_cast<std::int64_t>(gf.NumFaces());

    const std::size_t npts = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& pts_arr = rMesh.Points();
    std::vector<Vec3> xyz(npts);
    for (std::size_t i = 0; i < npts; ++i) {
        xyz[i] = {detail::read_double(pts_arr, i * pdim + 0),
                  pdim > 1 ? detail::read_double(pts_arr, i * pdim + 1) : 0.0,
                  pdim > 2 ? detail::read_double(pts_arr, i * pdim + 2) : 0.0};
    }

    // Cell centroids and volumes (global block-major indexing).
    const std::vector<double> centroids = detail::cell_centroids(rMesh, total_cells);
    std::vector<double> volumes(total_cells, 0.0);
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t base = static_cast<std::size_t>(bases[b]);
        if (cell_type_dimension(cell_type_from_name(cb.Type())) != 3)
            continue;
        parallel_for(cb.NumCells(), [&](std::size_t c) {
            detail::CellRings rings;
            std::vector<Vec3> local_pts;
            if (detail::cell_rings(cb, c, pts_arr, pdim, rings, local_pts)) {
                const detail::PolyMeasure pm = detail::poly_measure(rings, local_pts.data());
                volumes[base + c] = std::abs(pm.mVolume);
            }
        });
    }

    // Collect interior face indices.
    std::vector<std::size_t> interior_faces;
    interior_faces.reserve(gf.NumFaces());
    for (std::size_t f = 0; f < gf.NumFaces(); ++f) {
        if (gf.mNeighbour[f] >= 0)
            interior_faces.push_back(f);
    }
    const std::size_t n_int = interior_faces.size();
    report.mNumInteriorFaces = static_cast<std::int64_t>(n_int);

    std::vector<double> face_nonortho(n_int, qnan);
    std::vector<double> face_skew(n_int, qnan);
    std::vector<double> face_vratio(n_int, qnan);
    std::vector<double> face_pyr(n_int, qnan);

    constexpr double kRadToDeg = 180.0 / 3.14159265358979323846;

    parallel_for(n_int, [&](std::size_t idx) {
        const std::size_t f = interior_faces[idx];
        const std::int64_t o_compact = gf.mOwner[f];
        const std::int64_t n_compact = gf.mNeighbour[f];
        const std::size_t o_glob =
            static_cast<std::size_t>(gf.mCellToGlobal[static_cast<std::size_t>(o_compact)]);
        const std::size_t n_glob =
            static_cast<std::size_t>(gf.mCellToGlobal[static_cast<std::size_t>(n_compact)]);

        const Vec3 c_o{centroids[o_glob * 3 + 0], centroids[o_glob * 3 + 1],
                       centroids[o_glob * 3 + 2]};
        const Vec3 c_n{centroids[n_glob * 3 + 0], centroids[n_glob * 3 + 1],
                       centroids[n_glob * 3 + 2]};

        const std::size_t fn = gf.FaceSize(f);
        const std::int64_t* fids = gf.Face(f);
        std::vector<Vec3> fpts(fn);
        std::vector<std::uint32_t> ring(fn);
        for (std::size_t i = 0; i < fn; ++i) {
            fpts[i] = xyz[static_cast<std::size_t>(fids[i])];
            ring[i] = static_cast<std::uint32_t>(i);
        }
        const Vec3 area_vec = detail::polygon_area_vector(fpts.data(), fn);
        const double area = detail::vec3_norm(area_vec);
        const Vec3 c_f = detail::polygon_centroid(fpts.data(), ring.data(), fn);

        // Face diameter: max pairwise distance among face vertices.
        double diam2 = 0.0;
        for (std::size_t i = 0; i < fn; ++i) {
            const Vec3& pi = xyz[static_cast<std::size_t>(fids[i])];
            for (std::size_t j = i + 1; j < fn; ++j) {
                const Vec3& pj = xyz[static_cast<std::size_t>(fids[j])];
                diam2 = std::max(diam2, detail::vec3_norm_sq(detail::vec3_sub(pj, pi)));
            }
        }
        const double diam = std::sqrt(diam2);

        const Vec3 d = detail::vec3_sub(c_n, c_o);
        const double d_len = detail::vec3_norm(d);

        if (area > 1e-30 && d_len > 1e-30) {
            const Vec3 n_hat = detail::vec3_scale(area_vec, 1.0 / area);
            const Vec3 d_hat = detail::vec3_scale(d, 1.0 / d_len);
            const double cos_a = std::clamp(detail::vec3_dot(n_hat, d_hat), -1.0, 1.0);
            face_nonortho[idx] = std::acos(cos_a) * kRadToDeg;

            // Face skewness: intersection of C_O -> C_N with face plane.
            const double denom = detail::vec3_dot(d, n_hat);
            const double numer = detail::vec3_dot(detail::vec3_sub(c_f, c_o), n_hat);
            double skew = 0.0;
            if (std::abs(denom) > 1e-30 && diam > 1e-30) {
                const double t = numer / denom;
                const Vec3 p_i = detail::vec3_add(c_o, detail::vec3_scale(d, t));
                skew = detail::vec3_norm(detail::vec3_sub(p_i, c_f)) / diam;
            } else if (diam > 1e-30) {
                skew = std::numeric_limits<double>::infinity();
            }
            face_skew[idx] = skew;

            // Signed pyramid heights from owner and neighbour.
            const double h_o = detail::vec3_dot(detail::vec3_sub(c_f, c_o), n_hat);
            const double h_n = detail::vec3_dot(detail::vec3_sub(c_n, c_f), n_hat);
            face_pyr[idx] = std::min(h_o, h_n);
        }

        // Adjacent volume growth ratio.
        const double vo = volumes[o_glob];
        const double vn = volumes[n_glob];
        if (vo > 1e-30 && vn > 1e-30)
            face_vratio[idx] = std::max(vo, vn) / std::min(vo, vn);
    });

    // Serial deterministic projection onto owner and neighbour cells + severity counts.
    auto update_max = [](double& rDst, double Val) {
        if (!std::isfinite(Val))
            return;
        if (!std::isfinite(rDst) || Val > rDst)
            rDst = Val;
    };
    auto update_min = [](double& rDst, double Val) {
        if (!std::isfinite(Val))
            return;
        if (!std::isfinite(rDst) || Val < rDst)
            rDst = Val;
    };

    for (std::size_t idx = 0; idx < n_int; ++idx) {
        if (std::isfinite(face_nonortho[idx]) && face_nonortho[idx] > 70.0)
            ++report.mNumSevereNonOrtho;
        if (std::isfinite(face_pyr[idx]) && face_pyr[idx] < 0.0)
            ++report.mNumNegativePyramids;

        const std::size_t f = interior_faces[idx];
        const std::int64_t g_pair[2] = {
            gf.mCellToGlobal[static_cast<std::size_t>(gf.mOwner[f])],
            gf.mCellToGlobal[static_cast<std::size_t>(gf.mNeighbour[f])],
        };
        for (std::int64_t g : g_pair) {
            const auto [b, row] = detail::global_to_block_row(bases, g);
            if (b == static_cast<std::size_t>(-1))
                continue;
            const std::size_t c = static_cast<std::size_t>(row);
            update_max(cell_nonortho[b][c], face_nonortho[idx]);
            update_max(cell_skew[b][c], face_skew[idx]);
            update_max(cell_vratio[b][c], face_vratio[idx]);
            update_min(cell_pyr[b][c], face_pyr[idx]);
        }
    }

    report.mMetrics.emplace_back(kFvmNonOrthogonalityName, fvm_summarize(face_nonortho));
    report.mMetrics.emplace_back(kFvmFaceSkewnessName, fvm_summarize(face_skew));
    report.mMetrics.emplace_back(kFvmVolumeRatioName, fvm_summarize(face_vratio));
    report.mMetrics.emplace_back(kFvmMinPyramidHeightName, fvm_summarize(face_pyr));

    auto to_arrays = [&](const std::vector<std::vector<double>>& rBlocks) {
        std::vector<NDArray> out;
        out.reserve(nblocks);
        for (std::size_t b = 0; b < nblocks; ++b) {
            NDArray a(DType::Float64, {rBlocks[b].size()});
            std::copy(rBlocks[b].begin(), rBlocks[b].end(), a.As<double>());
            out.push_back(std::move(a));
        }
        return out;
    };

    report.mCellArrays.emplace_back(kFvmNonOrthogonalityName, to_arrays(cell_nonortho));
    report.mCellArrays.emplace_back(kFvmFaceSkewnessName, to_arrays(cell_skew));
    report.mCellArrays.emplace_back(kFvmVolumeRatioName, to_arrays(cell_vratio));
    report.mCellArrays.emplace_back(kFvmMinPyramidHeightName, to_arrays(cell_pyr));

    return report;
}

Mesh attach_fvm_quality(const Mesh& rMesh) {
    FvmQualityReport rep = compute_fvm_quality(rMesh);
    Mesh out = detail::clone_mesh(rMesh);
    for (auto& [name, blocks] : rep.mCellArrays)
        out.AddCellData(name, std::move(blocks));
    return out;
}

}  // namespace meshioplusplus
