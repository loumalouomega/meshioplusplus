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

#include "meshioplusplus/operations/deform.hpp"

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

#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

#include "../detail/bounding_box_utils.hpp"
#include "../detail/iga_kernels.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/point_knn.hpp"
#include "../detail/sparse_solve.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

std::string df_lower(const std::string& rS) {
    std::string out = rS;
    for (char& c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

std::vector<Vec3> df_read_xyz(const Mesh& rMesh) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& pts = rMesh.Points();
    std::vector<Vec3> xyz(n);
    for (std::size_t i = 0; i < n; ++i) {
        xyz[i] = {detail::read_double(pts, i * pdim + 0),
                  pdim > 1 ? detail::read_double(pts, i * pdim + 1) : 0.0,
                  pdim > 2 ? detail::read_double(pts, i * pdim + 2) : 0.0};
    }
    return xyz;
}

std::vector<std::uint8_t> df_region_mask(const Mesh& rMesh, const std::string& rRegion,
                                         const char* pOp) {
    const std::size_t n = rMesh.NumPoints();
    if (rRegion.empty())
        return std::vector<std::uint8_t>(n, 1);
    std::vector<std::uint8_t> mask(n, 0);
    const std::size_t pidx = rMesh.FindRegion(rRegion, RegionKind::Point);
    if (pidx != Mesh::npos) {
        const Region& reg = rMesh.Region(pidx);
        const std::int64_t* e = reg.Entries();
        for (std::size_t i = 0; i < reg.NumEntries(); ++i)
            if (e[i] >= 0 && static_cast<std::size_t>(e[i]) < n)
                mask[static_cast<std::size_t>(e[i])] = 1;
        return mask;
    }
    const std::size_t cidx = rMesh.FindRegion(rRegion, RegionKind::Cell);
    if (cidx != Mesh::npos) {
        const auto bases = detail::block_bases(rMesh);
        const Region& reg = rMesh.Region(cidx);
        const std::int64_t* e = reg.Entries();
        for (std::size_t i = 0; i < reg.NumEntries(); ++i) {
            const auto [b, row] = detail::global_to_block_row(bases, e[i]);
            if (b == static_cast<std::size_t>(-1))
                continue;
            const auto cb = rMesh.Cells(b);
            if (cb.IsRagged()) {
                const auto r = static_cast<std::size_t>(row);
                const std::int64_t* p = cb.Row(r);
                const std::size_t sz = cb.RowSize(r);
                for (std::size_t j = 0; j < sz; ++j)
                    if (p[j] >= 0 && static_cast<std::size_t>(p[j]) < n)
                        mask[static_cast<std::size_t>(p[j])] = 1;
            } else {
                const std::size_t k = cb.NodesPerCell();
                const NDArray& conn = cb.Conn();
                for (std::size_t j = 0; j < k; ++j) {
                    const std::int64_t v =
                        detail::read_int(conn, static_cast<std::size_t>(row) * k + j);
                    if (v >= 0 && static_cast<std::size_t>(v) < n)
                        mask[static_cast<std::size_t>(v)] = 1;
                }
            }
        }
        return mask;
    }
    throw std::invalid_argument(std::string(pOp) + ": unknown region '" + rRegion + "'");
}

std::vector<double> df_clamped_uniform_knots(std::size_t N, int P) {
    const std::size_t m = N + static_cast<std::size_t>(P) + 1;
    std::vector<double> knots(m, 0.0);
    for (std::size_t i = 0; i <= static_cast<std::size_t>(P); ++i)
        knots[i] = 0.0;
    const std::size_t n_int = N - static_cast<std::size_t>(P) - 1;
    for (std::size_t i = 1; i <= n_int; ++i)
        knots[static_cast<std::size_t>(P) + i] =
            static_cast<double>(i) / static_cast<double>(n_int + 1);
    for (std::size_t i = N; i < m; ++i)
        knots[i] = 1.0;
    return knots;
}

std::int64_t df_count_inverted(const Mesh& rMesh) {
    std::int64_t inverted = 0;
    const NDArray& pts = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cell_type_dimension(cell_type_from_name(cb.Type())) != 3)
            continue;
        detail::CellRings rings;
        std::vector<Vec3> local_pts;
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            if (detail::cell_rings(cb, c, pts, pdim, rings, local_pts)) {
                const detail::PolyMeasure pm = detail::poly_measure(rings, local_pts.data());
                if (pm.mVolume < 0.0)
                    ++inverted;
            }
        }
    }
    return inverted;
}

DeformResult df_finalize(const Mesh& rMesh, const std::vector<Vec3>& rOrig,
                         const std::vector<Vec3>& rNew) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    DeformResult res;
    res.mMesh = detail::clone_mesh(rMesh);
    NDArray out_pts(DType::Float64, {n, pdim});
    double* po = out_pts.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        const double d = detail::vec3_norm(detail::vec3_sub(rNew[i], rOrig[i]));
        if (d > 0.0) {
            ++res.mNumNodesMoved;
            res.mMaxDisplacement = std::max(res.mMaxDisplacement, d);
        }
        po[i * pdim + 0] = rNew[i][0];
        if (pdim > 1)
            po[i * pdim + 1] = rNew[i][1];
        if (pdim > 2)
            po[i * pdim + 2] = rNew[i][2];
    }
    res.mMesh.AssignPoints(std::move(out_pts));
    res.mNumInverted = df_count_inverted(res.mMesh);
    return res;
}

inline double df_wendland_c2(double R, double Radius) {
    if (R >= Radius)
        return 0.0;
    const double q = R / Radius;
    const double om = 1.0 - q;
    const double om2 = om * om;
    return om2 * om2 * (4.0 * q + 1.0);
}

inline double df_rbf_eval(RbfKernel Kernel, double R, double Radius) {
    if (Kernel == RbfKernel::Tps)
        return R;  // 3-D biharmonic RBF
    return df_wendland_c2(R, Radius);
}

}  // namespace

RbfKernel rbf_kernel_from_name(const std::string& rName) {
    const std::string low = df_lower(rName);
    if (low == "tps")
        return RbfKernel::Tps;
    if (low == "wendland_c2" || low == "wendland-c2" || low == "wendland")
        return RbfKernel::WendlandC2;
    throw std::invalid_argument("deform_rbf: unknown kernel '" + rName +
                                "' (expected 'tps' or 'wendland_c2')");
}

DeformResult deform_ffd(const Mesh& rMesh, const FfdOptions& rOptions) {
    const std::size_t npts = rMesh.NumPoints();
    if (npts == 0)
        throw std::invalid_argument("deform_ffd: mesh has no points");

    std::array<std::size_t, 3> counts{};
    std::array<int, 3> degs{};
    for (int d = 0; d < 3; ++d) {
        if (rOptions.mGridDims[static_cast<std::size_t>(d)] < 2)
            throw std::invalid_argument("deform_ffd: grid_dims must be >= 2 in each direction");
        counts[static_cast<std::size_t>(d)] =
            static_cast<std::size_t>(rOptions.mGridDims[static_cast<std::size_t>(d)]);
        const std::int64_t req = rOptions.mDegrees[static_cast<std::size_t>(d)];
        const int p = req < 0 ? static_cast<int>(counts[static_cast<std::size_t>(d)] - 1)
                              : static_cast<int>(req);
        if (p < 1 || static_cast<std::size_t>(p) >= counts[static_cast<std::size_t>(d)])
            throw std::invalid_argument("deform_ffd: degree must satisfy 1 <= p <= grid_dims - 1");
        degs[static_cast<std::size_t>(d)] = p;
    }

    const std::size_t n_ctrl = counts[0] * counts[1] * counts[2];
    if (rOptions.mDisplacements.size() != n_ctrl * 3)
        throw std::invalid_argument("deform_ffd: displacements length must be 3 * nx * ny * nz");

    const std::vector<Vec3> xyz = df_read_xyz(rMesh);
    const std::vector<std::uint8_t> mask = df_region_mask(rMesh, rOptions.mRegion, "deform_ffd");

    Vec3 lo = rOptions.mBoxMin;
    Vec3 hi = rOptions.mBoxMax;
    if (!rOptions.mHasBounds) {
        const auto bb = detail::BoundingBoxUtils::FromPoints(xyz);
        lo = bb.mMin;
        hi = bb.mMax;
    }
    Vec3 inv_span{};
    for (int d = 0; d < 3; ++d) {
        const double span = hi[static_cast<std::size_t>(d)] - lo[static_cast<std::size_t>(d)];
        inv_span[static_cast<std::size_t>(d)] = span > 1e-15 ? (1.0 / span) : 0.0;
    }

    const std::vector<double> kx = df_clamped_uniform_knots(counts[0], degs[0]);
    const std::vector<double> ky = df_clamped_uniform_knots(counts[1], degs[1]);
    const std::vector<double> kz = df_clamped_uniform_knots(counts[2], degs[2]);

    detail::IgaEvalNet net;
    net.mDim = 3;
    net.mDegrees = degs;
    net.mKnots = {&kx, &ky, &kz};
    net.mCounts = counts;
    net.mValues = rOptions.mDisplacements.data();
    net.mColumns = 3;
    net.mWeights = nullptr;

    std::vector<Vec3> out_xyz = xyz;
    parallel_for(npts, [&](std::size_t i) {
        if (!mask[i])
            return;
        double u[3];
        std::size_t span[3];
        for (int d = 0; d < 3; ++d) {
            const double t =
                (xyz[i][static_cast<std::size_t>(d)] - lo[static_cast<std::size_t>(d)]) *
                inv_span[static_cast<std::size_t>(d)];
            u[d] = std::clamp(t, 0.0, 1.0);
            span[d] = detail::nurbs_find_span(degs[static_cast<std::size_t>(d)],
                                              *net.mKnots[static_cast<std::size_t>(d)],
                                              counts[static_cast<std::size_t>(d)], u[d]);
        }
        double disp[3] = {0.0, 0.0, 0.0};
        detail::iga_eval_point(net, span, u, 0, disp, nullptr);
        out_xyz[i] = {xyz[i][0] + disp[0], xyz[i][1] + disp[1], xyz[i][2] + disp[2]};
    });

    return df_finalize(rMesh, xyz, out_xyz);
}

DeformResult deform_rbf(const Mesh& rMesh, const RbfOptions& rOptions) {
    const std::size_t npts = rMesh.NumPoints();
    if (npts == 0)
        throw std::invalid_argument("deform_rbf: mesh has no points");
    if (rOptions.mControlPoints.empty() || rOptions.mControlPoints.size() % 3 != 0)
        throw std::invalid_argument("deform_rbf: control_points must be non-empty flat 3*M");
    if (rOptions.mDisplacements.size() != rOptions.mControlPoints.size())
        throw std::invalid_argument("deform_rbf: displacements must match control_points length");
    if (!(rOptions.mRegularization >= 0.0) || !std::isfinite(rOptions.mRegularization))
        throw std::invalid_argument("deform_rbf: regularization must be >= 0 and finite");

    const std::size_t m = rOptions.mControlPoints.size() / 3;
    std::vector<Vec3> ctrl(m);
    for (std::size_t i = 0; i < m; ++i) {
        ctrl[i] = {rOptions.mControlPoints[i * 3 + 0], rOptions.mControlPoints[i * 3 + 1],
                   rOptions.mControlPoints[i * 3 + 2]};
    }

    const std::vector<Vec3> xyz = df_read_xyz(rMesh);
    const std::vector<std::uint8_t> mask = df_region_mask(rMesh, rOptions.mRegion, "deform_rbf");

    double radius = rOptions.mRadius;
    if (rOptions.mKernel == RbfKernel::WendlandC2 && radius <= 0.0) {
        const auto bb = detail::BoundingBoxUtils::FromPoints(xyz);
        radius = 0.5 * detail::BoundingBoxUtils::Diagonal(bb);
        if (!(radius > 0.0))
            radius = 1.0;
    }

    // Solve for RBF weights gamma (m x 3) and optional linear polynomial coeffs poly (4 x 3).
    std::array<std::vector<double>, 3> gamma;
    std::array<std::array<double, 4>, 3> poly = {{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}};
    for (int d = 0; d < 3; ++d)
        gamma[static_cast<std::size_t>(d)].assign(m, 0.0);

    const bool use_sparse = (rOptions.mKernel == RbfKernel::WendlandC2 && m > 512);
    const bool use_poly = rOptions.mPolynomial && m >= 4;

    if (!use_sparse) {
        const std::size_t k = use_poly ? (m + 4) : m;
        std::vector<double> mat(k * k, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            mat[i * k + i] = df_rbf_eval(rOptions.mKernel, 0.0, radius) + rOptions.mRegularization;
            for (std::size_t j = i + 1; j < m; ++j) {
                const double r = detail::vec3_norm(detail::vec3_sub(ctrl[i], ctrl[j]));
                const double val = df_rbf_eval(rOptions.mKernel, r, radius);
                mat[i * k + j] = val;
                mat[j * k + i] = val;
            }
            if (use_poly) {
                mat[i * k + m + 0] = 1.0;
                mat[i * k + m + 1] = ctrl[i][0];
                mat[i * k + m + 2] = ctrl[i][1];
                mat[i * k + m + 3] = ctrl[i][2];
                mat[(m + 0) * k + i] = 1.0;
                mat[(m + 1) * k + i] = ctrl[i][0];
                mat[(m + 2) * k + i] = ctrl[i][1];
                mat[(m + 3) * k + i] = ctrl[i][2];
            }
        }
        for (int d = 0; d < 3; ++d) {
            std::vector<double> a = mat;
            std::vector<double> b(k, 0.0);
            for (std::size_t i = 0; i < m; ++i)
                b[i] = rOptions.mDisplacements[i * 3 + static_cast<std::size_t>(d)];
            if (!detail::lsq_solve(a, b, k))
                throw std::invalid_argument(
                    "deform_rbf: singular RBF system (coincident or coplanar control points; try "
                    "adding regularization or disabling polynomial reproduction)");
            for (std::size_t i = 0; i < m; ++i)
                gamma[static_cast<std::size_t>(d)][i] = b[i];
            if (use_poly) {
                for (std::size_t p = 0; p < 4; ++p)
                    poly[static_cast<std::size_t>(d)][p] = b[m + p];
            }
        }
    } else {
        // Large Wendland C2: subtract least-squares linear trend if requested, then solve SPD Phi.
        std::array<std::vector<double>, 3> rhs;
        for (int d = 0; d < 3; ++d) {
            rhs[static_cast<std::size_t>(d)].resize(m);
            for (std::size_t i = 0; i < m; ++i)
                rhs[static_cast<std::size_t>(d)][i] =
                    rOptions.mDisplacements[i * 3 + static_cast<std::size_t>(d)];
        }
        if (use_poly) {
            std::vector<double> ata(16, 0.0);
            for (std::size_t i = 0; i < m; ++i) {
                const double row[4] = {1.0, ctrl[i][0], ctrl[i][1], ctrl[i][2]};
                for (std::size_t a = 0; a < 4; ++a)
                    for (std::size_t b = 0; b < 4; ++b)
                        ata[a * 4 + b] += row[a] * row[b];
            }
            for (int d = 0; d < 3; ++d) {
                std::vector<double> a = ata;
                std::vector<double> b(4, 0.0);
                for (std::size_t i = 0; i < m; ++i) {
                    const double row[4] = {1.0, ctrl[i][0], ctrl[i][1], ctrl[i][2]};
                    const double val = rhs[static_cast<std::size_t>(d)][i];
                    for (std::size_t k = 0; k < 4; ++k)
                        b[k] += row[k] * val;
                }
                if (detail::lsq_solve(a, b, 4)) {
                    for (std::size_t k = 0; k < 4; ++k)
                        poly[static_cast<std::size_t>(d)][k] = b[k];
                    for (std::size_t i = 0; i < m; ++i) {
                        rhs[static_cast<std::size_t>(d)][i] -=
                            b[0] + b[1] * ctrl[i][0] + b[2] * ctrl[i][1] + b[3] * ctrl[i][2];
                    }
                }
            }
        }
        const detail::PointKnn knn(ctrl, std::max(radius, 1e-12));
        std::vector<detail::SparseTriplet> triplets;
        triplets.reserve(m * 16);
        std::vector<std::pair<double, std::int64_t>> nbrs;
        for (std::size_t i = 0; i < m; ++i) {
            knn.WithinRadius(ctrl[i], radius, nbrs);
            for (const auto& [d2, j] : nbrs) {
                const double val =
                    df_wendland_c2(std::sqrt(d2), radius) +
                    (static_cast<std::size_t>(j) == i ? rOptions.mRegularization : 0.0);
                if (val != 0.0)
                    triplets.push_back({i, static_cast<std::size_t>(j), val});
            }
        }
        const detail::SparseCsr phi = detail::csr_from_triplets(m, m, std::move(triplets));
        detail::CsrIc0 precond;
        precond.Build(phi);
        detail::CsrPcgWorkspace ws;
        const double* rhs_ptrs[3] = {rhs[0].data(), rhs[1].data(), rhs[2].data()};
        double* x_ptrs[3] = {gamma[0].data(), gamma[1].data(), gamma[2].data()};
        detail::csr_ic0_pcg_block(phi, precond, 3, rhs_ptrs, x_ptrs, ws, 1e-12, m + 200);
    }

    // Evaluate displacements in parallel_for.
    std::vector<Vec3> out_xyz = xyz;
    if (rOptions.mKernel == RbfKernel::WendlandC2) {
        const detail::PointKnn knn(ctrl, std::max(radius, 1e-12));
        parallel_for(npts, [&](std::size_t i) {
            if (!mask[i])
                return;
            std::vector<std::pair<double, std::int64_t>> nbrs;
            knn.WithinRadius(xyz[i], radius, nbrs);
            Vec3 disp{0.0, 0.0, 0.0};
            for (const auto& [d2, j] : nbrs) {
                const double w = df_wendland_c2(std::sqrt(d2), radius);
                const auto ju = static_cast<std::size_t>(j);
                disp[0] += w * gamma[0][ju];
                disp[1] += w * gamma[1][ju];
                disp[2] += w * gamma[2][ju];
            }
            if (use_poly) {
                for (int d = 0; d < 3; ++d) {
                    const auto& p = poly[static_cast<std::size_t>(d)];
                    disp[static_cast<std::size_t>(d)] +=
                        p[0] + p[1] * xyz[i][0] + p[2] * xyz[i][1] + p[3] * xyz[i][2];
                }
            }
            out_xyz[i] = {xyz[i][0] + disp[0], xyz[i][1] + disp[1], xyz[i][2] + disp[2]};
        });
    } else {
        parallel_for(npts, [&](std::size_t i) {
            if (!mask[i])
                return;
            Vec3 disp{0.0, 0.0, 0.0};
            for (std::size_t j = 0; j < m; ++j) {
                const double r = detail::vec3_norm(detail::vec3_sub(xyz[i], ctrl[j]));
                disp[0] += r * gamma[0][j];
                disp[1] += r * gamma[1][j];
                disp[2] += r * gamma[2][j];
            }
            if (use_poly) {
                for (int d = 0; d < 3; ++d) {
                    const auto& p = poly[static_cast<std::size_t>(d)];
                    disp[static_cast<std::size_t>(d)] +=
                        p[0] + p[1] * xyz[i][0] + p[2] * xyz[i][1] + p[3] * xyz[i][2];
                }
            }
            out_xyz[i] = {xyz[i][0] + disp[0], xyz[i][1] + disp[1], xyz[i][2] + disp[2]};
        });
    }

    return df_finalize(rMesh, xyz, out_xyz);
}

}  // namespace meshioplusplus
