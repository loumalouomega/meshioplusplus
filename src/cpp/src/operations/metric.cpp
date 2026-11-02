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
// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/metric.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {

namespace {

using MtrMat = std::array<double, 9>;  // row-major 3x3 (2-D metrics fill the top-left 2x2)

/// Symmetric eigen-decomposition of the top-left D x D block by cyclic
/// Jacobi: eigenvalues `rW`, eigenvectors the columns of `rV`.
void mtr_eigen(MtrMat A, int D, std::array<double, 3>& rW, MtrMat& rV) {
    rV = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for (int sweep = 0; sweep < 50; ++sweep) {
        double off = 0.0, diag = 0.0;
        for (int i = 0; i < D; ++i) {
            diag += A[3 * i + i] * A[3 * i + i];
            for (int j = i + 1; j < D; ++j)
                off += A[3 * i + j] * A[3 * i + j];
        }
        if (off == 0.0 || off <= 1e-32 * diag)
            break;
        for (int p = 0; p < D; ++p)
            for (int q = p + 1; q < D; ++q) {
                const double apq = A[3 * p + q];
                if (apq == 0.0)
                    continue;
                const double theta = (A[3 * q + q] - A[3 * p + p]) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < D; ++k) {
                    const double akp = A[3 * k + p], akq = A[3 * k + q];
                    A[3 * k + p] = c * akp - s * akq;
                    A[3 * k + q] = s * akp + c * akq;
                }
                for (int k = 0; k < D; ++k) {
                    const double apk = A[3 * p + k], aqk = A[3 * q + k];
                    A[3 * p + k] = c * apk - s * aqk;
                    A[3 * q + k] = s * apk + c * aqk;
                }
                for (int k = 0; k < D; ++k) {
                    const double vkp = rV[3 * k + p], vkq = rV[3 * k + q];
                    rV[3 * k + p] = c * vkp - s * vkq;
                    rV[3 * k + q] = s * vkp + c * vkq;
                }
            }
    }
    for (int i = 0; i < 3; ++i)
        rW[static_cast<std::size_t>(i)] = i < D ? A[3 * i + i] : 0.0;
}

/// V diag(w) V^T over the top-left D x D block.
MtrMat mtr_compose(const MtrMat& rV, const std::array<double, 3>& rW, int D) {
    MtrMat m{};
    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j) {
            double s = 0.0;
            for (int k = 0; k < D; ++k)
                s += rV[3 * i + k] * rW[static_cast<std::size_t>(k)] * rV[3 * j + k];
            m[3 * i + j] = s;
        }
    return m;
}

MtrMat mtr_mul(const MtrMat& rA, const MtrMat& rB, int D) {
    MtrMat m{};
    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j) {
            double s = 0.0;
            for (int k = 0; k < D; ++k)
                s += rA[3 * i + k] * rB[3 * k + j];
            m[3 * i + j] = s;
        }
    return m;
}

/// The intersection A cap B by simultaneous reduction.
MtrMat mtr_intersect(const MtrMat& rA, const MtrMat& rB, int D) {
    std::array<double, 3> w;
    MtrMat v;
    mtr_eigen(rA, D, w, v);
    double top = 0.0;
    for (int k = 0; k < D; ++k)
        top = std::max(top, std::abs(w[static_cast<std::size_t>(k)]));
    if (!(top > 0.0))
        return rB;  // A asks for nothing
    std::array<double, 3> root{}, inv_root{};
    for (int k = 0; k < D; ++k) {
        const double lam = std::max(w[static_cast<std::size_t>(k)], 1e-300 + 1e-14 * top);
        root[static_cast<std::size_t>(k)] = std::sqrt(lam);
        inv_root[static_cast<std::size_t>(k)] = 1.0 / std::sqrt(lam);
    }
    const MtrMat half = mtr_compose(v, root, D);
    const MtrMat inv_half = mtr_compose(v, inv_root, D);
    const MtrMat c = mtr_mul(mtr_mul(inv_half, rB, D), inv_half, D);
    std::array<double, 3> mu;
    MtrMat q;
    mtr_eigen(c, D, mu, q);
    for (int k = 0; k < D; ++k)
        mu[static_cast<std::size_t>(k)] = std::max(1.0, mu[static_cast<std::size_t>(k)]);
    return mtr_mul(mtr_mul(half, mtr_compose(q, mu, D), D), half, D);
}

double mtr_det(const MtrMat& rM, int D) {
    if (D == 2)
        return rM[0] * rM[4] - rM[1] * rM[3];
    return rM[0] * (rM[4] * rM[8] - rM[5] * rM[7]) - rM[1] * (rM[3] * rM[8] - rM[5] * rM[6]) +
           rM[2] * (rM[3] * rM[7] - rM[4] * rM[6]);
}

/// Clamp the eigenvalues to [Lo, Hi] (either 0: that side open).
MtrMat mtr_clamp(const MtrMat& rM, int D, double Lo, double Hi) {
    std::array<double, 3> w;
    MtrMat v;
    mtr_eigen(rM, D, w, v);
    for (int k = 0; k < D; ++k) {
        double& x = w[static_cast<std::size_t>(k)];
        if (Lo > 0.0)
            x = std::max(x, Lo);
        if (Hi > 0.0)
            x = std::min(x, Hi);
    }
    return mtr_compose(v, w, D);
}

/// Medit lower-triangle order: m11 m21 m22 (m31 m32 m33).
MtrMat mtr_read(const NDArray& rA, std::size_t Row, int D) {
    MtrMat m{};
    const std::size_t w = D == 3 ? 6 : 3;
    const auto at = [&](std::size_t k) { return detail::read_double(rA, Row * w + k); };
    m[0] = at(0);
    m[1] = m[3] = at(1);
    m[4] = at(2);
    if (D == 3) {
        m[2] = m[6] = at(3);
        m[5] = m[7] = at(4);
        m[8] = at(5);
    }
    return m;
}

void mtr_write(const MtrMat& rM, int D, double* pOut) {
    pOut[0] = rM[0];
    pOut[1] = 0.5 * (rM[1] + rM[3]);
    pOut[2] = rM[4];
    if (D == 3) {
        pOut[3] = 0.5 * (rM[2] + rM[6]);
        pOut[4] = 0.5 * (rM[5] + rM[7]);
        pOut[5] = rM[8];
    }
}

double mtr_change(const MtrMat& rA, const MtrMat& rB, int D) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j) {
            num += (rA[3 * i + j] - rB[3 * i + j]) * (rA[3 * i + j] - rB[3 * i + j]);
            den += rB[3 * i + j] * rB[3 * i + j];
        }
    return den > 0.0 ? std::sqrt(num / den) : std::sqrt(num);
}

}  // namespace

PrepareMetricResult prepare_metric(const Mesh& rMesh, const PrepareMetricOptions& rOptions) {
    const char* op = "prepare_metric";
    if (rOptions.mArrays.empty())
        throw std::invalid_argument(std::string(op) + ": name at least one metric array");
    if (rOptions.mOutput.empty())
        throw std::invalid_argument(std::string(op) + ": the output name is empty");
    if (!(rOptions.mTargetVertices >= 0.0))
        throw std::invalid_argument(std::string(op) + ": target_vertices must not be negative");
    if (!(rOptions.mNorm == 0.0 || rOptions.mNorm >= 1.0))
        throw std::invalid_argument(std::string(op) + ": norm must be at least 1 (0: infinity)");
    if (!(rOptions.mGradation == 0.0 || rOptions.mGradation > 1.0))
        throw std::invalid_argument(std::string(op) + ": gradation must exceed 1 (0: none)");
    if (!(rOptions.mHmin >= 0.0) || !(rOptions.mHmax >= 0.0) ||
        (rOptions.mHmin > 0.0 && rOptions.mHmax > 0.0 && rOptions.mHmin > rOptions.mHmax))
        throw std::invalid_argument(std::string(op) + ": need 0 <= hmin <= hmax");
    const std::size_t n = rMesh.NumPoints();
    int d = 0;
    for (const std::string& name : rOptions.mArrays) {
        if (!rMesh.HasPointData(name))
            throw std::invalid_argument(std::string(op) + ": no point data named '" + name + "'");
        const std::size_t size = rMesh.PointData(name).Size();
        const int width = n > 0 && size == 6 * n ? 3 : (n > 0 && size == 3 * n ? 2 : 0);
        if (width == 0 || (d != 0 && width != d))
            throw std::invalid_argument(std::string(op) + ": '" + name +
                                        "' must be (n, 6) (3-D) or (n, 3) (2-D), like the others");
        d = width;
    }

    // 1. Intersection.
    std::vector<MtrMat> m(n);
    for (std::size_t i = 0; i < n; ++i) {
        m[i] = mtr_read(rMesh.PointData(rOptions.mArrays[0]), i, d);
        for (std::size_t a = 1; a < rOptions.mArrays.size(); ++a)
            m[i] = mtr_intersect(m[i], mtr_read(rMesh.PointData(rOptions.mArrays[a]), i, d), d);
    }

    // The dimension the measures are taken in: the metric's, or, for a 3-D
    // metric on a planar mesh (2-D points, or every z equal) with no volume
    // cell, the plane's (the metric's x-y block measures it).
    const auto dim_of = [](const Mesh::CellView& rCb) {
        return rCb.IsRagged() ? (rCb.IsPolyhedron() ? 3 : 2)
                              : cell_type_dimension(cell_type_from_name(rCb.Type()));
    };
    int deff = d;
    {
        bool has_d = false, has_2 = false;
        for (const auto cb : rMesh.CellRange()) {
            has_d = has_d || dim_of(cb) == d;
            has_2 = has_2 || dim_of(cb) == 2;
        }
        bool planar = rMesh.PointDim() < 3;
        if (!planar && n > 0) {
            planar = true;
            const double z0 = detail::read_double(rMesh.Points(), 2);
            for (std::size_t i = 1; i < n && planar; ++i)
                planar = detail::read_double(rMesh.Points(), i * rMesh.PointDim() + 2) == z0;
        }
        if (d == 3 && !has_d && has_2 && planar)
            deff = 2;
    }
    // The lumped point measures of the cells of that dimension.
    std::vector<double> weight(n, 0.0);
    for (const auto cb : rMesh.CellRange()) {
        const int cd = dim_of(cb);
        if (cd != deff)
            continue;
        std::vector<std::int64_t> ids;
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const double meas = detail::cell_measure(rMesh.Points(), rMesh.PointDim(), cb, c);
            if (!std::isfinite(meas))
                continue;
            detail::cell_node_ids(cb, c, n, ids);
            std::sort(ids.begin(), ids.end());
            ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
            for (const std::int64_t v : ids)
                weight[static_cast<std::size_t>(v)] += meas / static_cast<double>(ids.size());
        }
    }
    const auto complexity = [&]() {
        double c = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            c += weight[i] * std::sqrt(std::max(0.0, mtr_det(m[i], deff)));
        return c;
    };
    PrepareMetricResult result;
    result.mComplexityBefore = complexity();

    // 2. Normalisation.
    if (rOptions.mTargetVertices > 0.0) {
        double total = 0.0;
        for (const double w : weight)
            total += w;
        if (!(total > 0.0))
            throw std::invalid_argument(std::string(op) +
                                        ": normalising needs cells of the metric's dimension");
        const double dd = static_cast<double>(deff);
        const double p = rOptions.mNorm;
        const double e_int = p == 0.0 ? 0.5 : p / (2.0 * p + dd);
        const double e_loc = p == 0.0 ? 0.0 : -1.0 / (2.0 * p + dd);
        double integral = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            integral += weight[i] * std::pow(std::max(0.0, mtr_det(m[i], deff)), e_int);
        if (!(integral > 0.0))
            throw std::invalid_argument(std::string(op) + ": the metric is zero everywhere");
        const double global =
            std::pow(rOptions.mTargetVertices, 2.0 / dd) * std::pow(integral, -2.0 / dd);
        for (std::size_t i = 0; i < n; ++i) {
            const double det = mtr_det(m[i], deff);
            const double s = global * (det > 0.0 ? std::pow(det, e_loc) : 1.0);
            for (double& x : m[i])
                x *= s;
        }
    }
    // 3. Clamping.
    const double lo = rOptions.mHmax > 0.0 ? 1.0 / (rOptions.mHmax * rOptions.mHmax) : 0.0;
    const double hi = rOptions.mHmin > 0.0 ? 1.0 / (rOptions.mHmin * rOptions.mHmin) : 0.0;
    if (lo > 0.0 || hi > 0.0)
        for (MtrMat& x : m)
            x = mtr_clamp(x, d, lo, hi);

    // 4. Gradation over the edges, in order.
    if (rOptions.mGradation > 1.0) {
        const detail::NodeAdjacency adj =
            detail::build_node_adjacency(rMesh, n, detail::NodeAdjacencyKind::Edge);
        const double log_beta = std::log(rOptions.mGradation);
        const std::size_t pdim = rMesh.PointDim();
        std::vector<char> changed(n, 0);
        const auto edge_length = [&](const MtrMat& rM, std::size_t P, std::size_t Q) {
            double e[3] = {0.0, 0.0, 0.0};
            for (std::size_t k = 0; k < pdim && k < 3; ++k)
                e[k] = detail::read_double(rMesh.Points(), Q * pdim + k) -
                       detail::read_double(rMesh.Points(), P * pdim + k);
            double l2 = 0.0;
            for (int i = 0; i < d; ++i)
                for (int j = 0; j < d; ++j)
                    l2 += e[i] * rM[3 * i + j] * e[j];
            return std::sqrt(std::max(0.0, l2));
        };
        const auto spread = [&](std::size_t From, std::size_t To) {
            const double l = edge_length(m[From], From, To);
            const double eta = 1.0 / ((1.0 + l * log_beta) * (1.0 + l * log_beta));
            MtrMat grown = m[From];
            for (double& x : grown)
                x *= eta;
            const MtrMat next = mtr_intersect(m[To], grown, d);
            if (mtr_change(next, m[To], d) > 1e-3) {
                m[To] = next;
                changed[To] = 1;
                return true;
            }
            return false;
        };
        for (std::int64_t sweep = 1; sweep <= rOptions.mMaxSweeps; ++sweep) {
            bool any = false;
            for (std::size_t p = 0; p < n; ++p)
                for (std::int64_t k = adj.mXadj[p]; k < adj.mXadj[p + 1]; ++k) {
                    const auto q = static_cast<std::size_t>(adj.mAdj[static_cast<std::size_t>(k)]);
                    if (q <= p)
                        continue;
                    any = spread(p, q) || any;
                    any = spread(q, p) || any;
                }
            result.mSweeps = sweep;
            if (!any)
                break;
        }
        for (const char c : changed)
            result.mNumGraded += c;
        if (lo > 0.0 || hi > 0.0)
            for (MtrMat& x : m)
                x = mtr_clamp(x, d, lo, hi);
    }
    result.mComplexityAfter = complexity();

    const std::size_t width = d == 3 ? 6 : 3;
    NDArray out(DType::Float64, {n, width});
    for (std::size_t i = 0; i < n; ++i)
        mtr_write(m[i], d, out.As<double>() + i * width);
    result.mMesh = detail::clone_mesh(rMesh);
    result.mMesh.AddPointData(rOptions.mOutput, std::move(out));
    return result;
}

}  // namespace meshioplusplus
