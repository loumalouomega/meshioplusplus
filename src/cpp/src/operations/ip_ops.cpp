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

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <vector>

// Project includes
#include "../detail/integration_points.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/quadrature.hpp"
#include "../detail/reference_element.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/ip_ops.hpp"

namespace meshioplusplus {
namespace {

Mesh ip_clone_mesh(const Mesh& rMesh, const std::set<std::string>& rSkipPointData = {},
                   const std::set<std::string>& rSkipCellData = {},
                   const std::set<std::string>& rSkipFieldData = {}) {
    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rMesh.Points()));
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.IsRagged()) {
            if (cb.IsPolyhedron()) {
                std::vector<std::vector<std::vector<std::int64_t>>> cells(cb.NumCells());
                for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                    cells[r].resize(cb.NumFaces(r));
                    for (std::size_t f = 0; f < cb.NumFaces(r); ++f) {
                        const auto face = cb.Face(r, f);
                        cells[r][f].assign(face.first, face.first + face.second);
                    }
                }
                out.AddPolyhedronBlock(cb.Type(), cells);
            } else {
                std::vector<std::vector<std::int64_t>> rows(cb.NumCells());
                for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                    const auto* pRow = cb.Row(r);
                    rows[r].assign(pRow, pRow + cb.RowSize(r));
                }
                out.AddPolygonBlock(cb.Type(), rows);
            }
        } else {
            out.AddCellBlock(cb.Type(), detail::data_owned_copy(cb.Conn()));
        }
    }
    for (const std::string& name : rMesh.PointDataNames()) {
        if (!rSkipPointData.count(name))
            out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    }
    for (const std::string& name : rMesh.CellDataNames()) {
        if (rSkipCellData.count(name))
            continue;
        std::vector<NDArray> blocks;
        blocks.reserve(rMesh.CellDataNumBlocks(name));
        for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b)
            blocks.push_back(detail::data_owned_copy(rMesh.CellData(name, b)));
        out.AddCellData(name, std::move(blocks));
    }
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (!rSkipFieldData.count(name))
            out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));
    }
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
        out.AddRegion(rMesh.Region(r));
    return out;
}

void ip_collect_layout_keys(const std::string& rName, std::set<std::string>& rSkipField) {
    const std::string prefixes[] = {
        detail::kIpLayoutPrefix,         "abaqus:layout:",
        "lsdyna_d3plot:layout:",         "gid:gauss_points:",
        "gmsh:element_node_components:",
    };
    for (const std::string& pre : prefixes)
        rSkipField.insert(pre + rName);
}

std::vector<std::string> ip_select_cell_arrays(const Mesh& rMesh,
                                               const std::vector<std::string>& rRequested,
                                               std::vector<std::string>& rWarnings) {
    std::vector<std::string> out;
    if (rRequested.empty()) {
        for (const std::string& name : rMesh.CellDataNames()) {
            if (detail::get_ip_layout(rMesh, name).mValid)
                out.push_back(name);
        }
        return out;
    }
    for (const std::string& name : rRequested) {
        if (!rMesh.HasCellData(name)) {
            rWarnings.push_back("cell_data array '" + name + "' not found");
            continue;
        }
        if (!detail::get_ip_layout(rMesh, name).mValid) {
            rWarnings.push_back("cell_data array '" + name + "' has no integration-point layout");
            continue;
        }
        out.push_back(name);
    }
    return out;
}

/// Build the `npc x n_ip` extrapolation matrix `rE` (row-major) from the
/// `n_ip x npc` shape-function evaluation matrix `M` (`M[q * npc + j] = N_j(xi_q)`).
/// Returns false if singular.
bool ip_build_extrapolation_matrix(const std::vector<double>& rM, std::size_t NIp, std::size_t Npc,
                                   std::vector<double>& rE) {
    rE.assign(Npc * NIp, 0.0);
    if (NIp == 0 || Npc == 0)
        return false;

    if (NIp >= Npc) {
        // Solve (M^T M) E = M^T column-by-column for each of the NIp columns.
        std::vector<double> mtm(Npc * Npc, 0.0);
        for (std::size_t i = 0; i < Npc; ++i) {
            for (std::size_t j = 0; j < Npc; ++j) {
                double s = 0.0;
                for (std::size_t q = 0; q < NIp; ++q)
                    s += rM[q * Npc + i] * rM[q * Npc + j];
                mtm[i * Npc + j] = s;
            }
        }
        for (std::size_t q = 0; q < NIp; ++q) {
            std::vector<double> a = mtm;
            std::vector<double> rhs(Npc, 0.0);
            for (std::size_t i = 0; i < Npc; ++i)
                rhs[i] = rM[q * Npc + i];
            if (!detail::lsq_solve(a, rhs, Npc))
                return false;
            for (std::size_t i = 0; i < Npc; ++i)
                rE[i * NIp + q] = rhs[i];
        }
        return true;
    }

    // Underdetermined (NIp < Npc): minimum-norm least-squares E = M^T (M M^T)^{-1}.
    // Since sum_j N_j(xi_q) = 1, constant fields are reproduced identically.
    std::vector<double> mmt(NIp * NIp, 0.0);
    for (std::size_t p = 0; p < NIp; ++p) {
        for (std::size_t q = 0; q < NIp; ++q) {
            double s = 0.0;
            for (std::size_t j = 0; j < Npc; ++j)
                s += rM[p * Npc + j] * rM[q * Npc + j];
            mmt[p * NIp + q] = s;
        }
    }
    std::vector<double> inv(NIp * NIp, 0.0);
    for (std::size_t col = 0; col < NIp; ++col) {
        std::vector<double> a = mmt;
        std::vector<double> rhs(NIp, 0.0);
        rhs[col] = 1.0;
        if (!detail::lsq_solve(a, rhs, NIp))
            return false;
        for (std::size_t row = 0; row < NIp; ++row)
            inv[row * NIp + col] = rhs[row];
    }
    for (std::size_t j = 0; j < Npc; ++j) {
        for (std::size_t col = 0; col < NIp; ++col) {
            double s = 0.0;
            for (std::size_t p = 0; p < NIp; ++p)
                s += rM[p * Npc + j] * inv[p * NIp + col];
            rE[j * NIp + col] = s;
        }
    }
    return true;
}

}  // namespace

IpOpResult ip_to_points(const Mesh& rMesh, const IpToPointsOptions& rOptions) {
    IpOpResult res;
    const std::vector<std::string> targets =
        ip_select_cell_arrays(rMesh, rOptions.mArrays, res.mWarnings);

    std::set<std::string> skip_cell;
    std::set<std::string> skip_field;
    if (rOptions.mDropSource) {
        for (const std::string& name : targets) {
            skip_cell.insert(name);
            ip_collect_layout_keys(name, skip_field);
        }
    }
    res.mMesh = ip_clone_mesh(rMesh, {}, skip_cell, skip_field);

    const std::size_t n_pts = rMesh.NumPoints();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (const std::string& name : targets) {
        const detail::IpLayout layout = detail::get_ip_layout(rMesh, name);
        const std::size_t n_ip = static_cast<std::size_t>(layout.mPoints);
        const std::size_t width = static_cast<std::size_t>(layout.mWidth);

        std::vector<double> accum(n_pts * width, 0.0);
        std::vector<std::size_t> count(n_pts, 0);

        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            if (b >= rMesh.CellDataNumBlocks(name))
                continue;
            const auto cb = rMesh.Cells(b);
            if (cb.IsRagged() || cb.NumCells() == 0)
                continue;
            const NDArray& arr = rMesh.CellData(name, b);
            const std::size_t n_cells = cb.NumCells();
            const std::size_t total_cols = detail::cols(arr);
            if (detail::rows(arr) != n_cells || total_cols < n_ip * width)
                continue;

            const std::size_t npc = cb.NodesPerCell();
            const CellType ct = cell_type_from_name(cb.Type());

            // Determine extrapolation strategy for this block.
            bool is_elno = (layout.mRuleId == 0 && n_ip == npc);
            bool has_extrap = false;
            std::vector<double> extrap;

            if (!is_elno) {
                detail::QuadratureRule rule;
                if (detail::resolve_ip_quadrature(rMesh, name, ct, layout, rule)) {
                    std::vector<double> M(n_ip * npc, 0.0);
                    bool ok = true;
                    std::vector<double> N;
                    std::vector<std::array<double, 3>> dN;
                    for (std::size_t q = 0; q < n_ip; ++q) {
                        if (!detail::eval_reference_shape(ct, npc, rule.mPoints[q], N, dN) ||
                            N.size() != npc) {
                            ok = false;
                            break;
                        }
                        for (std::size_t j = 0; j < npc; ++j)
                            M[q * npc + j] = N[j];
                    }
                    if (ok)
                        has_extrap = ip_build_extrapolation_matrix(M, n_ip, npc, extrap);
                }
            }

            const NDArray& conn = cb.Conn();
            std::vector<double> ip_vals(n_ip, 0.0);
            for (std::size_t r = 0; r < n_cells; ++r) {
                for (std::size_t c = 0; c < width; ++c) {
                    bool all_finite = true;
                    double mean = 0.0;
                    std::size_t n_fin = 0;
                    for (std::size_t q = 0; q < n_ip; ++q) {
                        const double v = detail::read_double(arr, r * total_cols + q * width + c);
                        ip_vals[q] = v;
                        if (std::isfinite(v)) {
                            mean += v;
                            ++n_fin;
                        } else {
                            all_finite = false;
                        }
                    }
                    if (n_fin == 0)
                        continue;
                    mean /= static_cast<double>(n_fin);

                    for (std::size_t j = 0; j < npc; ++j) {
                        const std::int64_t node = detail::read_int(conn, r * npc + j);
                        if (node < 0 || static_cast<std::size_t>(node) >= n_pts)
                            continue;
                        double v_node = mean;
                        if (is_elno) {
                            v_node = std::isfinite(ip_vals[j]) ? ip_vals[j] : mean;
                        } else if (has_extrap && all_finite) {
                            double s = 0.0;
                            for (std::size_t q = 0; q < n_ip; ++q)
                                s += extrap[j * n_ip + q] * ip_vals[q];
                            v_node = s;
                        }
                        accum[static_cast<std::size_t>(node) * width + c] += v_node;
                        if (c == 0)
                            count[static_cast<std::size_t>(node)] += 1;
                    }
                }
            }
        }

        std::vector<std::size_t> shape =
            width == 1 ? std::vector<std::size_t>{n_pts} : std::vector<std::size_t>{n_pts, width};
        NDArray pd(DType::Float64, shape);
        double* dst = pd.As<double>();
        for (std::size_t i = 0; i < n_pts; ++i) {
            if (count[i] == 0) {
                for (std::size_t c = 0; c < width; ++c)
                    dst[i * width + c] = nan;
            } else {
                const double inv = 1.0 / static_cast<double>(count[i]);
                for (std::size_t c = 0; c < width; ++c)
                    dst[i * width + c] = accum[i * width + c] * inv;
            }
        }
        res.mMesh.AddPointData(name, std::move(pd));
        res.mConverted.push_back(name);
    }

    return res;
}

IpOpResult ip_to_cells(const Mesh& rMesh, const IpToCellsOptions& rOptions) {
    IpOpResult res;
    const std::vector<std::string> targets =
        ip_select_cell_arrays(rMesh, rOptions.mArrays, res.mWarnings);

    std::set<std::string> skip_field;
    if (rOptions.mDropMetadata) {
        for (const std::string& name : targets)
            ip_collect_layout_keys(name, skip_field);
    }
    res.mMesh = ip_clone_mesh(rMesh, {}, {}, skip_field);

    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (const std::string& name : targets) {
        const detail::IpLayout layout = detail::get_ip_layout(rMesh, name);
        const std::size_t n_ip = static_cast<std::size_t>(layout.mPoints);
        const std::size_t width = static_cast<std::size_t>(layout.mWidth);

        std::vector<NDArray> out_blocks;
        out_blocks.reserve(rMesh.NumCellBlocks());

        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const auto cb = rMesh.Cells(b);
            const std::size_t n_cells = cb.NumCells();
            std::vector<std::size_t> shape = width == 1 ? std::vector<std::size_t>{n_cells}
                                                        : std::vector<std::size_t>{n_cells, width};
            NDArray cd(DType::Float64, shape);
            double* dst = cd.As<double>();
            std::fill_n(dst, cd.Size(), nan);

            if (b < rMesh.CellDataNumBlocks(name) && n_cells > 0) {
                const NDArray& arr = rMesh.CellData(name, b);
                const std::size_t total_cols = detail::cols(arr);
                if (detail::rows(arr) == n_cells && total_cols >= n_ip * width) {
                    const CellType ct = cell_type_from_name(cb.Type());
                    std::vector<double> weights(n_ip, 1.0);
                    detail::QuadratureRule rule;
                    if (detail::resolve_ip_quadrature(rMesh, name, ct, layout, rule) &&
                        rule.mWeights.size() == n_ip) {
                        weights = rule.mWeights;
                    }

                    for (std::size_t r = 0; r < n_cells; ++r) {
                        for (std::size_t c = 0; c < width; ++c) {
                            double wsum = 0.0;
                            double vsum = 0.0;
                            for (std::size_t q = 0; q < n_ip; ++q) {
                                const double v =
                                    detail::read_double(arr, r * total_cols + q * width + c);
                                if (std::isfinite(v)) {
                                    vsum += weights[q] * v;
                                    wsum += weights[q];
                                }
                            }
                            if (wsum > 0.0)
                                dst[r * width + c] = vsum / wsum;
                        }
                    }
                }
            }
            out_blocks.push_back(std::move(cd));
        }

        res.mMesh.AddCellData(name, std::move(out_blocks));
        res.mConverted.push_back(name);
    }

    return res;
}

IpOpResult points_to_ip(const Mesh& rMesh, const PointsToIpOptions& rOptions) {
    IpOpResult res;

    if (rOptions.mRuleId < 0)
        throw std::invalid_argument(
            "points_to_ip: rule_id must be >= 0 (0 = element-nodal, >= 1 = quadrature order)");

    std::vector<std::string> targets;
    if (rOptions.mArrays.empty()) {
        for (const std::string& name : rMesh.PointDataNames()) {
            const NDArray& pd = rMesh.PointData(name);
            if (pd.Size() > 0)
                targets.push_back(name);
        }
    } else {
        for (const std::string& name : rOptions.mArrays) {
            if (!rMesh.HasPointData(name)) {
                res.mWarnings.push_back("point_data array '" + name + "' not found");
                continue;
            }
            targets.push_back(name);
        }
    }

    // Precompute shape matrices per block.
    const std::size_t nb = rMesh.NumCellBlocks();
    std::vector<std::size_t> block_nip(nb, 0);
    std::vector<std::vector<double>> block_M(nb);
    std::size_t max_nip = 0;

    for (std::size_t b = 0; b < nb; ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.IsRagged() || cb.NumCells() == 0)
            continue;
        const std::size_t npc = cb.NodesPerCell();
        const CellType ct = cell_type_from_name(cb.Type());
        if (rOptions.mRuleId == 0) {
            block_nip[b] = npc;
            max_nip = std::max(max_nip, npc);
        } else {
            const detail::QuadratureRule rule =
                detail::quadrature_rule(ct, static_cast<int>(rOptions.mRuleId));
            if (rule.mPoints.empty())
                continue;
            const std::size_t n_ip = rule.mPoints.size();
            std::vector<double> M(n_ip * npc, 0.0);
            std::vector<double> N;
            std::vector<std::array<double, 3>> dN;
            bool ok = true;
            for (std::size_t q = 0; q < n_ip; ++q) {
                if (!detail::eval_reference_shape(ct, npc, rule.mPoints[q], N, dN) ||
                    N.size() != npc) {
                    ok = false;
                    break;
                }
                for (std::size_t j = 0; j < npc; ++j)
                    M[q * npc + j] = N[j];
            }
            if (ok) {
                block_nip[b] = n_ip;
                block_M[b] = std::move(M);
                max_nip = std::max(max_nip, n_ip);
            }
        }
    }

    if (max_nip == 0) {
        res.mMesh = ip_clone_mesh(rMesh);
        res.mWarnings.push_back("points_to_ip: no supported cell blocks for rule_id " +
                                std::to_string(rOptions.mRuleId));
        return res;
    }

    const std::size_t n_pts = rMesh.NumPoints();
    std::set<std::string> skip_point;
    if (rOptions.mDropSource) {
        for (const std::string& name : targets) {
            if (detail::rows(rMesh.PointData(name)) == n_pts)
                skip_point.insert(name);
        }
    }
    res.mMesh = ip_clone_mesh(rMesh, skip_point, {}, {});

    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (const std::string& name : targets) {
        const NDArray& pd = rMesh.PointData(name);
        if (detail::rows(pd) != n_pts) {
            res.mWarnings.push_back("point_data array '" + name + "' row count mismatch");
            continue;
        }
        const std::size_t width = detail::cols(pd);
        const std::size_t total_cols = max_nip * width;

        std::vector<NDArray> out_blocks;
        out_blocks.reserve(nb);
        for (std::size_t b = 0; b < nb; ++b) {
            const auto cb = rMesh.Cells(b);
            const std::size_t n_cells = cb.NumCells();
            NDArray cd(DType::Float64, {n_cells, total_cols});
            double* dst = cd.As<double>();
            std::fill_n(dst, cd.Size(), nan);

            const std::size_t n_ip = block_nip[b];
            if (n_ip > 0 && n_cells > 0) {
                const std::size_t npc = cb.NodesPerCell();
                const NDArray& conn = cb.Conn();
                for (std::size_t r = 0; r < n_cells; ++r) {
                    if (rOptions.mRuleId == 0) {
                        for (std::size_t j = 0; j < npc; ++j) {
                            const std::int64_t node = detail::read_int(conn, r * npc + j);
                            if (node < 0 || static_cast<std::size_t>(node) >= n_pts)
                                continue;
                            for (std::size_t c = 0; c < width; ++c) {
                                dst[r * total_cols + j * width + c] = detail::read_double(
                                    pd, static_cast<std::size_t>(node) * width + c);
                            }
                        }
                    } else {
                        const std::vector<double>& M = block_M[b];
                        for (std::size_t q = 0; q < n_ip; ++q) {
                            for (std::size_t c = 0; c < width; ++c) {
                                double s = 0.0;
                                for (std::size_t j = 0; j < npc; ++j) {
                                    const std::int64_t node = detail::read_int(conn, r * npc + j);
                                    if (node < 0 || static_cast<std::size_t>(node) >= n_pts)
                                        continue;
                                    s += M[q * npc + j] *
                                         detail::read_double(
                                             pd, static_cast<std::size_t>(node) * width + c);
                                }
                                dst[r * total_cols + q * width + c] = s;
                            }
                        }
                    }
                }
            }
            out_blocks.push_back(std::move(cd));
        }

        res.mMesh.AddCellData(name, std::move(out_blocks));
        detail::set_ip_layout(res.mMesh, name, static_cast<std::int64_t>(max_nip),
                              static_cast<std::int64_t>(width), rOptions.mRuleId);
        res.mConverted.push_back(name);
    }

    return res;
}

}  // namespace meshioplusplus
