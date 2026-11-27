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

#include "meshioplusplus/operations/critical_points.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/parallel.hpp"

#include "../detail/jacobi_eigen.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/reference_element.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

bool cp_in_ref_cell(const detail::ReferenceElementInfo& rInfo, const std::array<double, 3>& rXi,
                    double Eps = 1e-6) {
    const double u = rXi[0], v = rXi[1], w = rXi[2];
    if (u < -Eps || u > 1.0 + Eps)
        return false;
    if (rInfo.mDim >= 2 && (v < -Eps || v > 1.0 + Eps))
        return false;
    if (rInfo.mDim >= 3 && (w < -Eps || w > 1.0 + Eps))
        return false;
    if (rInfo.mShape == lagrange::Shape::Triangle || rInfo.mShape == lagrange::Shape::Wedge) {
        if (u + v > 1.0 + Eps)
            return false;
    } else if (rInfo.mShape == lagrange::Shape::Tetra) {
        if (u + v + w > 1.0 + Eps)
            return false;
    }
    return true;
}

std::vector<std::array<double, 3>> cp_seeds(const detail::ReferenceElementInfo& rInfo, int N) {
    std::vector<std::array<double, 3>> out;
    const int n = std::max(N, 2);
    auto coord = [n](int i) { return (static_cast<double>(i) + 0.5) / static_cast<double>(n); };
    if (rInfo.mDim == 1) {
        for (int i = 0; i < n; ++i)
            out.push_back({coord(i), 0.0, 0.0});
    } else if (rInfo.mDim == 2) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                const std::array<double, 3> xi{coord(i), coord(j), 0.0};
                if (cp_in_ref_cell(rInfo, xi, 0.0))
                    out.push_back(xi);
            }
    } else {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    const std::array<double, 3> xi{coord(i), coord(j), coord(k)};
                    if (cp_in_ref_cell(rInfo, xi, 0.0))
                        out.push_back(xi);
                }
    }
    return out;
}

struct CpCandidate {
    Vec3 mPos;
    double mValue = 0.0;
    std::int32_t mType = 3;
    std::int64_t mCell = -1;
};

}  // namespace

CriticalPointsResult critical_points(const Mesh& rMesh, const CriticalPointsOptions& rOptions) {
    if (rOptions.mArray.empty())
        throw std::invalid_argument("critical_points: array name is required");
    if (!rMesh.HasPointData(rOptions.mArray))
        throw std::invalid_argument(
            data_unknown_key_message(rMesh, DataLocation::Point, rOptions.mArray));
    const NDArray& pdata = rMesh.PointData(rOptions.mArray);
    if (data_num_components(pdata) != 1)
        throw std::invalid_argument("critical_points: array must be scalar (1 component)");
    if (!(rOptions.mTolerance > 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument("critical_points: tolerance must be > 0");

    const NDArray& pts_in = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    const auto bases = detail::block_bases(rMesh);

    std::vector<CpCandidate> all_cands;

    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t nc = cb.NumCells();
        if (nc == 0 || cb.IsRagged())
            continue;
        const std::size_t k = cb.NodesPerCell();
        const CellType ct = cell_type_from_name(cb.Type());
        const NDArray& conn = cb.Conn();
        detail::ReferenceElementInfo info;
        if (!detail::reference_element_info(ct, k, info))
            continue;
        if (info.mOrder < 2)
            continue;  // Linear fields on linear cells have constant gradient (no isolated interior
                       // extrema).

        const int dim = info.mDim;
        const int n_sub = rOptions.mSeedSubdivisions > 0
                              ? static_cast<int>(rOptions.mSeedSubdivisions)
                              : std::max(info.mOrder + 1, 2);
        const auto seeds = cp_seeds(info, n_sub);
        const std::int64_t base = bases[b];

        std::vector<std::vector<CpCandidate>> per_cell(nc);
        parallel_for(nc, [&](std::size_t c) {
            std::vector<Vec3> cpts(k);
            std::vector<double> fvals(k);
            for (std::size_t i = 0; i < k; ++i) {
                const auto u = static_cast<std::size_t>(detail::read_int(conn, c * k + i));
                cpts[i] = {detail::read_double(pts_in, u * pdim + 0),
                           pdim > 1 ? detail::read_double(pts_in, u * pdim + 1) : 0.0,
                           pdim > 2 ? detail::read_double(pts_in, u * pdim + 2) : 0.0};
                fvals[i] = detail::read_double(pdata, u);
            }

            auto eval_ref_grad = [&](const std::array<double, 3>& rXi, double& rVal,
                                     std::array<double, 3>& rGrad) -> bool {
                std::vector<double> nv;
                std::vector<std::array<double, 3>> ng;
                if (!detail::eval_reference_shape(ct, k, rXi, nv, ng))
                    return false;
                rVal = 0.0;
                rGrad = {0.0, 0.0, 0.0};
                for (std::size_t i = 0; i < k; ++i) {
                    rVal += nv[i] * fvals[i];
                    for (int d = 0; d < dim; ++d)
                        rGrad[static_cast<std::size_t>(d)] +=
                            ng[i][static_cast<std::size_t>(d)] * fvals[i];
                }
                return true;
            };

            for (const auto& seed : seeds) {
                std::array<double, 3> xi = seed;
                double val = 0.0;
                std::array<double, 3> g{0.0, 0.0, 0.0};
                bool converged = false;

                for (std::int64_t iter = 0; iter < rOptions.mMaxIterations; ++iter) {
                    if (!eval_ref_grad(xi, val, g))
                        break;
                    double gnorm = 0.0;
                    for (int d = 0; d < dim; ++d)
                        gnorm += g[static_cast<std::size_t>(d)] * g[static_cast<std::size_t>(d)];
                    gnorm = std::sqrt(gnorm);
                    if (gnorm <= rOptions.mTolerance) {
                        converged = true;
                        break;
                    }

                    // Reference Hessian by central differences of exact reference gradient.
                    constexpr double h = 1e-6;
                    std::vector<double> H(static_cast<std::size_t>(dim * dim), 0.0);
                    for (int j = 0; j < dim; ++j) {
                        std::array<double, 3> xp = xi, xm = xi;
                        xp[static_cast<std::size_t>(j)] += h;
                        xm[static_cast<std::size_t>(j)] -= h;
                        double vp = 0.0, vm = 0.0;
                        std::array<double, 3> gp{}, gm{};
                        if (!eval_ref_grad(xp, vp, gp) || !eval_ref_grad(xm, vm, gm))
                            break;
                        for (int i = 0; i < dim; ++i)
                            H[static_cast<std::size_t>(i * dim + j)] =
                                (gp[static_cast<std::size_t>(i)] -
                                 gm[static_cast<std::size_t>(i)]) /
                                (2.0 * h);
                    }
                    // Symmetrize H and solve H * step = -g.
                    for (int i = 0; i < dim; ++i)
                        for (int j = i + 1; j < dim; ++j) {
                            const double avg = 0.5 * (H[static_cast<std::size_t>(i * dim + j)] +
                                                      H[static_cast<std::size_t>(j * dim + i)]);
                            H[static_cast<std::size_t>(i * dim + j)] = avg;
                            H[static_cast<std::size_t>(j * dim + i)] = avg;
                        }
                    std::vector<double> rhs(static_cast<std::size_t>(dim));
                    for (int i = 0; i < dim; ++i)
                        rhs[static_cast<std::size_t>(i)] = -g[static_cast<std::size_t>(i)];
                    if (!detail::lsq_solve(H, rhs, static_cast<std::size_t>(dim)))
                        break;

                    // Damped step to stay well-behaved.
                    double step_norm = 0.0;
                    for (int i = 0; i < dim; ++i)
                        step_norm +=
                            rhs[static_cast<std::size_t>(i)] * rhs[static_cast<std::size_t>(i)];
                    step_norm = std::sqrt(step_norm);
                    const double scale = step_norm > 0.5 ? (0.5 / step_norm) : 1.0;
                    for (int i = 0; i < dim; ++i)
                        xi[static_cast<std::size_t>(i)] += scale * rhs[static_cast<std::size_t>(i)];
                    if (!cp_in_ref_cell(info, xi, 0.25))
                        break;
                }

                if (!converged || !cp_in_ref_cell(info, xi, 1e-5))
                    continue;

                // Evaluate physical position & classify by reference Hessian eigenvalues
                // (Sylvester's law of inertia: at a critical point nabla f = 0, the signature of
                // the physical Hessian J^{-T} H_xi J^{-1} is identical to H_xi's!).
                std::vector<double> nv;
                std::vector<std::array<double, 3>> ng;
                if (!detail::eval_reference_shape(ct, k, xi, nv, ng))
                    continue;
                Vec3 pos{0.0, 0.0, 0.0};
                for (std::size_t i = 0; i < k; ++i) {
                    pos[0] += nv[i] * cpts[i][0];
                    pos[1] += nv[i] * cpts[i][1];
                    pos[2] += nv[i] * cpts[i][2];
                }

                constexpr double h = 1e-6;
                std::vector<double> H(static_cast<std::size_t>(dim * dim), 0.0);
                for (int j = 0; j < dim; ++j) {
                    std::array<double, 3> xp = xi, xm = xi;
                    xp[static_cast<std::size_t>(j)] += h;
                    xm[static_cast<std::size_t>(j)] -= h;
                    double vp = 0.0, vm = 0.0;
                    std::array<double, 3> gp{}, gm{};
                    eval_ref_grad(xp, vp, gp);
                    eval_ref_grad(xm, vm, gm);
                    for (int i = 0; i < dim; ++i)
                        H[static_cast<std::size_t>(i * dim + j)] =
                            (gp[static_cast<std::size_t>(i)] - gm[static_cast<std::size_t>(i)]) /
                            (2.0 * h);
                }
                for (int i = 0; i < dim; ++i)
                    for (int j = i + 1; j < dim; ++j) {
                        const double avg = 0.5 * (H[static_cast<std::size_t>(i * dim + j)] +
                                                  H[static_cast<std::size_t>(j * dim + i)]);
                        H[static_cast<std::size_t>(i * dim + j)] = avg;
                        H[static_cast<std::size_t>(j * dim + i)] = avg;
                    }
                std::vector<double> evals(static_cast<std::size_t>(dim), 0.0);
                std::vector<double> evecs(static_cast<std::size_t>(dim * dim), 0.0);
                detail::jacobi_eigen_sym(H.data(), static_cast<std::size_t>(dim), evals.data(),
                                         evecs.data(), 1e-12);

                int n_pos = 0, n_neg = 0;
                const double eig_tol = 1e-8;
                for (int i = 0; i < dim; ++i) {
                    if (evals[static_cast<std::size_t>(i)] > eig_tol)
                        ++n_pos;
                    else if (evals[static_cast<std::size_t>(i)] < -eig_tol)
                        ++n_neg;
                }
                CriticalPointType ctype = CriticalPointType::Degenerate;
                if (n_pos == dim)
                    ctype = CriticalPointType::Minimum;
                else if (n_neg == dim)
                    ctype = CriticalPointType::Maximum;
                else if (n_pos > 0 && n_neg > 0)
                    ctype = CriticalPointType::Saddle;

                per_cell[c].push_back({pos, val, static_cast<std::int32_t>(ctype),
                                       base + static_cast<std::int64_t>(c)});
            }
        });

        for (auto& vec : per_cell)
            for (auto& cand : vec)
                all_cands.push_back(cand);
    }

    // Spatial deduplication in deterministic cell-order.
    const double tol2 = rOptions.mMergeTolerance * rOptions.mMergeTolerance;
    std::vector<CpCandidate> uniq;
    for (const auto& c : all_cands) {
        bool dup = false;
        for (const auto& u : uniq) {
            if (detail::vec3_norm_sq(detail::vec3_sub(c.mPos, u.mPos)) <= tol2) {
                dup = true;
                break;
            }
        }
        if (!dup)
            uniq.push_back(c);
    }

    CriticalPointsResult res;
    const std::size_t nu = uniq.size();
    Mesh out;
    NDArray pts(DType::Float64, {nu, std::size_t{3}});
    NDArray vals(DType::Float64, {nu});
    NDArray types(DType::Int32, {nu});
    NDArray cells(DType::Int64, {nu});
    double* pp = pts.As<double>();
    double* pv = vals.As<double>();
    std::int32_t* pt = types.As<std::int32_t>();
    std::int64_t* pc = cells.As<std::int64_t>();
    for (std::size_t i = 0; i < nu; ++i) {
        pp[i * 3 + 0] = uniq[i].mPos[0];
        pp[i * 3 + 1] = uniq[i].mPos[1];
        pp[i * 3 + 2] = uniq[i].mPos[2];
        pv[i] = uniq[i].mValue;
        pt[i] = uniq[i].mType;
        pc[i] = uniq[i].mCell;
        switch (static_cast<CriticalPointType>(uniq[i].mType)) {
            case CriticalPointType::Minimum:
                ++res.mNumMinima;
                break;
            case CriticalPointType::Maximum:
                ++res.mNumMaxima;
                break;
            case CriticalPointType::Saddle:
                ++res.mNumSaddles;
                break;
            default:
                ++res.mNumDegenerate;
                break;
        }
    }
    out.AssignPoints(std::move(pts));
    if (nu > 0) {
        NDArray vconn(DType::Int64, {nu, std::size_t{1}});
        std::int64_t* pvc = vconn.As<std::int64_t>();
        for (std::size_t i = 0; i < nu; ++i)
            pvc[i] = static_cast<std::int64_t>(i);
        out.AddCellBlock(cell_type_name(CellType::Vertex), std::move(vconn));
    }
    out.AddPointData(kCriticalValueName, std::move(vals));
    out.AddPointData(kCriticalTypeName, std::move(types));
    out.AddPointData(kCriticalCellName, std::move(cells));
    res.mMesh = std::move(out);
    return res;
}

}  // namespace meshioplusplus
