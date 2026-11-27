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

#include "meshioplusplus/operations/evaluate_cells.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/parallel.hpp"

#include "../detail/reference_element.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

std::string ec_lower(const std::string& rS) {
    std::string out = rS;
    for (char& c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

CellType ec_target_cell_type(lagrange::Shape S, BasisTarget Target) {
    if (Target == BasisTarget::Bezier) {
        switch (S) {
            case lagrange::Shape::Line:
                return CellType::VtkBezierCurve;
            case lagrange::Shape::Triangle:
                return CellType::VtkBezierTriangle;
            case lagrange::Shape::Quad:
                return CellType::VtkBezierQuadrilateral;
            case lagrange::Shape::Tetra:
                return CellType::VtkBezierTetrahedron;
            case lagrange::Shape::Hexahedron:
                return CellType::VtkBezierHexahedron;
            case lagrange::Shape::Wedge:
                return CellType::VtkBezierWedge;
        }
    }
    switch (S) {
        case lagrange::Shape::Line:
            return CellType::VtkLagrangeCurve;
        case lagrange::Shape::Triangle:
            return CellType::VtkLagrangeTriangle;
        case lagrange::Shape::Quad:
            return CellType::VtkLagrangeQuadrilateral;
        case lagrange::Shape::Tetra:
            return CellType::VtkLagrangeTetrahedron;
        case lagrange::Shape::Hexahedron:
            return CellType::VtkLagrangeHexahedron;
        case lagrange::Shape::Wedge:
            return CellType::VtkLagrangeWedge;
    }
    return CellType::VtkLagrangeCurve;
}

}  // namespace

BasisTarget basis_target_from_name(const std::string& rName) {
    const std::string low = ec_lower(rName);
    if (low == "bezier" || low == "bézier")
        return BasisTarget::Bezier;
    if (low == "lagrange")
        return BasisTarget::Lagrange;
    throw std::invalid_argument("change_basis: unknown target '" + rName +
                                "' (expected 'bezier' or 'lagrange')");
}

EvaluateCellsResult evaluate_cells(const Mesh& rMesh, const EvaluateCellsOptions& rOptions) {
    const std::size_t nq = rOptions.mCells.size();
    if (rOptions.mReferencePoints.size() != nq * 3)
        throw std::invalid_argument("evaluate_cells: reference_points length must be 3 * len(cells)");

    std::vector<std::string> fields = rOptions.mFields;
    if (fields.empty()) {
        fields = rMesh.PointDataNames();
    } else {
        for (const std::string& f : fields)
            if (!rMesh.HasPointData(f))
                throw std::invalid_argument(data_unknown_key_message(rMesh, DataLocation::Point, f));
    }

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const NDArray& pts_in = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    const double qnan = std::numeric_limits<double>::quiet_NaN();

    EvaluateCellsResult result;
    result.mPoints = NDArray(DType::Float64, {nq, std::size_t{3}});
    double* p_out = result.mPoints.As<double>();
    std::fill(p_out, p_out + nq * 3, qnan);

    std::vector<std::size_t> f_comps(fields.size());
    std::vector<const NDArray*> f_ptrs(fields.size());
    result.mValues.reserve(fields.size());
    if (rOptions.mDerivatives)
        result.mGradients.reserve(fields.size());

    for (std::size_t fi = 0; fi < fields.size(); ++fi) {
        const NDArray& arr = rMesh.PointData(fields[fi]);
        f_ptrs[fi] = &arr;
        const std::size_t c = data_num_components(arr);
        f_comps[fi] = c;
        NDArray v_arr =
            c == 1 ? NDArray(DType::Float64, {nq}) : NDArray(DType::Float64, {nq, c});
        std::fill(v_arr.As<double>(), v_arr.As<double>() + v_arr.Size(), qnan);
        result.mValues.emplace_back(fields[fi], std::move(v_arr));
        if (rOptions.mDerivatives) {
            NDArray g_arr = c == 1 ? NDArray(DType::Float64, {nq, std::size_t{3}})
                                   : NDArray(DType::Float64, {nq, c, std::size_t{3}});
            std::fill(g_arr.As<double>(), g_arr.As<double>() + g_arr.Size(), qnan);
            result.mGradients.emplace_back(fields[fi], std::move(g_arr));
        }
    }

    std::vector<std::uint8_t> ok_flags(nq, 0);
    parallel_for(nq, [&](std::size_t q) {
        const auto [b, row] = detail::global_to_block_row(bases, rOptions.mCells[q]);
        if (b == static_cast<std::size_t>(-1))
            return;
        const auto cb = rMesh.Cells(b);
        const CellType ct = cell_type_from_name(cb.Type());
        const std::size_t c = static_cast<std::size_t>(row);
        std::vector<std::int64_t> nids;
        if (cb.IsRagged()) {
            const std::int64_t* r = cb.Row(c);
            nids.assign(r, r + cb.RowSize(c));
        } else {
            const std::size_t k = cb.NodesPerCell();
            const NDArray& conn = cb.Conn();
            nids.resize(k);
            for (std::size_t i = 0; i < k; ++i)
                nids[i] = detail::read_int(conn, c * k + i);
        }
        const std::size_t k = nids.size();
        std::vector<Vec3> cpts(k);
        for (std::size_t i = 0; i < k; ++i) {
            const auto u = static_cast<std::size_t>(nids[i]);
            cpts[i] = {detail::read_double(pts_in, u * pdim + 0),
                       pdim > 1 ? detail::read_double(pts_in, u * pdim + 1) : 0.0,
                       pdim > 2 ? detail::read_double(pts_in, u * pdim + 2) : 0.0};
        }
        const std::array<double, 3> xi = {rOptions.mReferencePoints[q * 3 + 0],
                                          rOptions.mReferencePoints[q * 3 + 1],
                                          rOptions.mReferencePoints[q * 3 + 2]};
        std::vector<double> vals;
        std::vector<Vec3> grads;
        double meas = 0.0;
        if (rOptions.mDerivatives) {
            if (!detail::shape_eval_and_jacobian(ct, cpts, xi, vals, grads, meas))
                return;
        } else {
            std::vector<std::array<double, 3>> dummy;
            if (!detail::eval_reference_shape(ct, k, xi, vals, dummy))
                return;
        }

        // Physical position.
        Vec3 pos{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < k; ++i) {
            pos[0] += vals[i] * cpts[i][0];
            pos[1] += vals[i] * cpts[i][1];
            pos[2] += vals[i] * cpts[i][2];
        }
        p_out[q * 3 + 0] = pos[0];
        p_out[q * 3 + 1] = pos[1];
        p_out[q * 3 + 2] = pos[2];

        // Field values & optional physical gradients.
        for (std::size_t fi = 0; fi < fields.size(); ++fi) {
            const NDArray& src = *f_ptrs[fi];
            const std::size_t nc = f_comps[fi];
            double* pv = result.mValues[fi].second.As<double>() + q * nc;
            for (std::size_t comp = 0; comp < nc; ++comp) {
                double acc = 0.0;
                for (std::size_t i = 0; i < k; ++i)
                    acc += vals[i] *
                           detail::read_double(src, static_cast<std::size_t>(nids[i]) * nc + comp);
                pv[comp] = acc;
            }
            if (rOptions.mDerivatives) {
                double* pg = result.mGradients[fi].second.As<double>() + q * nc * 3;
                for (std::size_t comp = 0; comp < nc; ++comp) {
                    double gx = 0.0, gy = 0.0, gz = 0.0;
                    for (std::size_t i = 0; i < k; ++i) {
                        const double fv = detail::read_double(
                            src, static_cast<std::size_t>(nids[i]) * nc + comp);
                        gx += grads[i][0] * fv;
                        gy += grads[i][1] * fv;
                        gz += grads[i][2] * fv;
                    }
                    pg[comp * 3 + 0] = gx;
                    pg[comp * 3 + 1] = gy;
                    pg[comp * 3 + 2] = gz;
                }
            }
        }
        ok_flags[q] = 1;
    });

    for (std::uint8_t fl : ok_flags)
        result.mNumEvaluated += fl;
    return result;
}

ChangeBasisResult change_basis(const Mesh& rMesh, const ChangeBasisOptions& rOptions) {
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& pts_in = rMesh.Points();

    std::vector<std::string> field_names;
    std::vector<std::size_t> field_comps;
    if (rOptions.mTransformFields) {
        field_names = rMesh.PointDataNames();
        field_comps.resize(field_names.size());
        for (std::size_t i = 0; i < field_names.size(); ++i)
            field_comps[i] = data_num_components(rMesh.PointData(field_names[i]));
    }

    // Per-cell discontinuous control-point transformation (exact and general for any order),
    // preserving exact polynomial representation on every cell.
    std::vector<Vec3> new_pts;
    std::vector<std::vector<double>> new_fields(field_names.size());
    Mesh out;
    std::vector<NDArray> degrees_blocks;
    std::int64_t num_converted = 0;

    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t nc = cb.NumCells();
        if (nc == 0)
            continue;
        if (cb.IsRagged())
            throw std::invalid_argument("change_basis: ragged blocks are not supported");
        const std::size_t k = cb.NodesPerCell();
        const CellType ct = cell_type_from_name(cb.Type());
        const NDArray& in_conn = cb.Conn();
        detail::ReferenceElementInfo info;
        if (!detail::reference_element_info(ct, k, info) || info.mIsSerendipity ||
            info.mIsPyramid)
            throw std::invalid_argument(
                "change_basis: unsupported cell type '" + cb.Type() +
                "' (expected full Lagrange or Bezier cells)");

        const bool same_basis =
            (info.mIsBezier && rOptions.mTarget == BasisTarget::Bezier) ||
            (!info.mIsBezier && rOptions.mTarget == BasisTarget::Lagrange);

        std::vector<double> mat;
        if (!same_basis) {
            mat = rOptions.mTarget == BasisTarget::Bezier
                      ? lagrange::lagrange_to_bezier_matrix(info.mShape, info.mOrder)
                      : lagrange::bezier_to_lagrange_matrix(info.mShape, info.mOrder);
        }

        const CellType out_type = ec_target_cell_type(info.mShape, rOptions.mTarget);
        NDArray conn(DType::Int64, {nc, k});
        std::int64_t* pc = conn.As<std::int64_t>();
        NDArray deg_arr(DType::Int32, {nc, std::size_t{3}});
        std::int32_t* pdeg = deg_arr.As<std::int32_t>();

        for (std::size_t c = 0; c < nc; ++c) {
            const std::size_t base_pt = new_pts.size();
            std::vector<std::int64_t> nids(k);
            for (std::size_t i = 0; i < k; ++i)
                nids[i] = detail::read_int(in_conn, c * k + i);
            // Fix Hexahedron27 vertical-edge swap (18 <-> 19) to VTK Lagrange lattice order.
            if (ct == CellType::Hexahedron27)
                std::swap(nids[18], nids[19]);

            std::vector<Vec3> cpts(k);
            for (std::size_t i = 0; i < k; ++i) {
                const auto u = static_cast<std::size_t>(nids[i]);
                cpts[i] = {detail::read_double(pts_in, u * pdim + 0),
                           pdim > 1 ? detail::read_double(pts_in, u * pdim + 1) : 0.0,
                           pdim > 2 ? detail::read_double(pts_in, u * pdim + 2) : 0.0};
            }

            for (std::size_t row = 0; row < k; ++row) {
                Vec3 p{0.0, 0.0, 0.0};
                if (same_basis) {
                    p = cpts[row];
                } else {
                    for (std::size_t col = 0; col < k; ++col) {
                        const double w = mat[row * k + col];
                        p[0] += w * cpts[col][0];
                        p[1] += w * cpts[col][1];
                        p[2] += w * cpts[col][2];
                    }
                }
                new_pts.push_back(p);
                pc[c * k + row] = static_cast<std::int64_t>(base_pt + row);
            }

            for (std::size_t fi = 0; fi < field_names.size(); ++fi) {
                const NDArray& src = rMesh.PointData(field_names[fi]);
                const std::size_t ncomp = field_comps[fi];
                auto& dst = new_fields[fi];
                for (std::size_t row = 0; row < k; ++row) {
                    for (std::size_t comp = 0; comp < ncomp; ++comp) {
                        double val = 0.0;
                        if (same_basis) {
                            val = detail::read_double(
                                src, static_cast<std::size_t>(nids[row]) * ncomp + comp);
                        } else {
                            for (std::size_t col = 0; col < k; ++col)
                                val += mat[row * k + col] *
                                       detail::read_double(
                                           src, static_cast<std::size_t>(nids[col]) * ncomp + comp);
                        }
                        dst.push_back(val);
                    }
                }
            }

            const std::int32_t ord = static_cast<std::int32_t>(info.mOrder);
            pdeg[c * 3 + 0] = ord;
            pdeg[c * 3 + 1] = info.mDim >= 2 ? ord : 0;
            pdeg[c * 3 + 2] = info.mDim >= 3 ? ord : 0;
            ++num_converted;
        }
        out.AddCellBlock(cell_type_name(out_type), std::move(conn));
        degrees_blocks.push_back(std::move(deg_arr));
    }

    NDArray out_pts(DType::Float64, {new_pts.size(), pdim});
    double* po = out_pts.As<double>();
    for (std::size_t i = 0; i < new_pts.size(); ++i) {
        po[i * pdim + 0] = new_pts[i][0];
        if (pdim > 1)
            po[i * pdim + 1] = new_pts[i][1];
        if (pdim > 2)
            po[i * pdim + 2] = new_pts[i][2];
    }
    out.AssignPoints(std::move(out_pts));

    for (std::size_t fi = 0; fi < field_names.size(); ++fi) {
        const std::size_t ncomp = field_comps[fi];
        NDArray fa = ncomp == 1 ? NDArray(DType::Float64, {new_pts.size()})
                                : NDArray(DType::Float64, {new_pts.size(), ncomp});
        std::copy(new_fields[fi].begin(), new_fields[fi].end(), fa.As<double>());
        out.AddPointData(field_names[fi], std::move(fa));
    }

    if (!degrees_blocks.empty())
        out.AddCellData("HigherOrderDegrees", std::move(degrees_blocks));

    ChangeBasisResult res;
    res.mMesh = std::move(out);
    res.mNumConverted = num_converted;
    return res;
}

}  // namespace meshioplusplus
