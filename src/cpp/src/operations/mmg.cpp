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
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/mmg.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/providers.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/lsq_fit.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

bool mmg_available() noexcept {
    return remesher() != nullptr;
}

bool mmg_has_elas() noexcept {
    const auto p_remesher = remesher();
    return p_remesher && p_remesher->mElas;
}

int mmg_index_bits() noexcept {
    const auto p_remesher = remesher();
    return p_remesher ? p_remesher->mIndexBits : 0;
}

namespace {

[[noreturn]] void mmg_absent(const char* pOp) {
    throw Unsupported(std::string(pOp) +
                      ": the MMG backend is not available: install the meshioplusplus-mmg "
                      "package (pip install meshioplusplus[mmg]), or rebuild with "
                      "-DMESHIOPLUSPLUS_WITH_MMG=ON and a shared MMG install "
                      "(build/build-mmg.sh; see doc/mmg.md)");
}

}  // namespace
// --- metric_from_hessian: no MMG needed ---------------------------------------------------

Mesh metric_from_hessian(const Mesh& rMesh, const MetricFromHessianOptions& rOptions) {
    const char* op = "metric_from_hessian";
    if (!rMesh.HasPointData(rOptions.mHessian))
        throw std::invalid_argument(std::string(op) + ": no point_data '" + rOptions.mHessian +
                                    "'");
    if (!(rOptions.mError > 0.0) || !(rOptions.mHmin > 0.0) || !(rOptions.mHmax >= rOptions.mHmin))
        throw std::invalid_argument(std::string(op) +
                                    ": error and hmin must be positive and hmax >= hmin");
    if (rOptions.mOutput.empty())
        throw std::invalid_argument(std::string(op) + ": the output name is empty");
    const NDArray& h = rMesh.PointData(rOptions.mHessian);
    const std::size_t n = rMesh.NumPoints();
    const std::size_t width = n ? h.Size() / n : 0;
    if (h.Size() != n * width || (width != 9 && width != 4))
        throw std::invalid_argument(std::string(op) + ": '" + rOptions.mHessian +
                                    "' needs 9 components (4 in 2-D) per point, a row-major "
                                    "Hessian such as hessian() writes");
    const std::size_t d = width == 9 ? 3 : 2;
    const double c = 0.5 * std::pow(static_cast<double>(d) / static_cast<double>(d + 1), 2.0);
    const double lo = 1.0 / (rOptions.mHmax * rOptions.mHmax);
    const double hi = 1.0 / (rOptions.mHmin * rOptions.mHmin);
    const std::size_t out_width = d == 3 ? 6 : 3;
    NDArray metric(DType::Float64, {n, out_width});
    double* m = metric.As<double>();
    const detail::DoubleView hv(h);
    for (std::size_t i = 0; i < n; ++i) {
        // The symmetrised Hessian as a 3x3 (a 2-D one padded with an inert axis).
        double a[9] = {0, 0, 0, 0, 0, 0, 0, 0, 1};
        for (std::size_t r = 0; r < d; ++r)
            for (std::size_t q = 0; q < d; ++q)
                a[r * 3 + q] = 0.5 * (hv[i * width + r * d + q] + hv[i * width + q * d + r]);
        double w[3], v[9];
        detail::lsq_sym3_eigen(a, w, v);
        double full[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        for (std::size_t k = 0; k < 3; ++k) {
            // The padding axis has eigenvector e3 and no say in a 2-D metric.
            if (d == 2 && std::abs(v[k * 3 + 2]) > 0.5)
                continue;
            double lambda = c * std::abs(w[k]) / rOptions.mError;
            if (!std::isfinite(lambda))
                lambda = hi;
            lambda = std::clamp(lambda, lo, hi);
            for (std::size_t r = 0; r < 3; ++r)
                for (std::size_t q = 0; q < 3; ++q)
                    full[r * 3 + q] += lambda * v[k * 3 + r] * v[k * 3 + q];
        }
        double* row = m + i * out_width;
        // The lower triangle, row by row: the order of a Medit .sol.
        row[0] = full[0];
        row[1] = full[3];
        row[2] = full[4];
        if (d == 3) {
            row[3] = full[6];
            row[4] = full[7];
            row[5] = full[8];
        }
    }
    Mesh out = detail::clone_mesh(rMesh);
    out.AddPointData(rOptions.mOutput, std::move(metric));
    return out;
}

// --- the provider's operations -------------------------------------------------------------

AdaptResult adapt(const Mesh& rMesh, const AdaptOptions& rOptions) {
    const auto p_remesher = remesher();
    if (!p_remesher || !p_remesher->mAdapt)
        mmg_absent("adapt");
    return p_remesher->mAdapt(rMesh, rOptions);
}

DiscretizeLevelsetResult discretize_levelset(const Mesh& rMesh,
                                             const DiscretizeLevelsetOptions& rOptions) {
    const auto p_remesher = remesher();
    if (!p_remesher || !p_remesher->mDiscretizeLevelset)
        mmg_absent("discretize_levelset");
    return p_remesher->mDiscretizeLevelset(rMesh, rOptions);
}

MeshDomainResult mesh_domain(const Mesh& rBoundary, const MeshDomainOptions& rOptions) {
    const auto p_remesher = remesher();
    if (!p_remesher || !p_remesher->mMeshDomain)
        mmg_absent("mesh_domain");
    return p_remesher->mMeshDomain(rBoundary, rOptions);
}

AdaptResult move(const Mesh& rMesh, const MoveOptions& rOptions) {
    const auto p_remesher = remesher();
    if (!p_remesher || !p_remesher->mMove)
        mmg_absent("move");
    return p_remesher->mMove(rMesh, rOptions);
}

}  // namespace meshioplusplus
