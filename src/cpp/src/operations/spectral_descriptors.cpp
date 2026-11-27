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

#include "meshioplusplus/operations/spectral_descriptors.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/eigenmodes.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

std::string sd_lower(const std::string& rS) {
    std::string out = rS;
    for (char& c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

}  // namespace

SpectralDescriptorKind spectral_descriptor_kind_from_name(const std::string& rName) {
    const std::string low = sd_lower(rName);
    if (low == "hks")
        return SpectralDescriptorKind::Hks;
    if (low == "wks")
        return SpectralDescriptorKind::Wks;
    throw std::invalid_argument("spectral_descriptors: unknown kind '" + rName +
                                "' (expected 'hks' or 'wks')");
}

SpectralDescriptorsResult spectral_descriptors(const Mesh& rMesh,
                                               const SpectralDescriptorsOptions& rOptions) {
    if (rOptions.mNumScales < 1)
        throw std::invalid_argument("spectral_descriptors: num_scales must be >= 1");
    if (rOptions.mWksVariance < 0.0 || !std::isfinite(rOptions.mWksVariance))
        throw std::invalid_argument("spectral_descriptors: wks_variance must be >= 0 and finite");

    const std::size_t n_pts = rMesh.NumPoints();
    SpectralDescriptorsResult result;

    // 1. Obtain eigenmodes (reuse precomputed if present and requested).
    std::vector<double> eigvals;
    std::vector<double> modes_flat;
    std::size_t n_modes = 0;

    if (rOptions.mReuseEigenmodes && rMesh.HasPointData(kLaplacianModesName) &&
        rMesh.HasFieldData(kLaplacianEigenvaluesName)) {
        const NDArray& modes_arr = rMesh.PointData(kLaplacianModesName);
        const NDArray& vals_arr = rMesh.FieldData(kLaplacianEigenvaluesName);
        if (modes_arr.Shape().size() == 2 && modes_arr.Shape()[0] == n_pts &&
            vals_arr.Size() == modes_arr.Shape()[1] && vals_arr.Size() >= 2 &&
            modes_arr.Dtype() == DType::Float64 && vals_arr.Dtype() == DType::Float64) {
            n_modes = vals_arr.Size();
            const double* pv = vals_arr.As<double>();
            const double* pm = modes_arr.As<double>();
            eigvals.assign(pv, pv + n_modes);
            modes_flat.assign(pm, pm + n_pts * n_modes);
            result.mMesh = detail::clone_mesh(rMesh);
            result.mReusedEigenmodes = true;
        }
    }

    if (!result.mReusedEigenmodes) {
        if (rOptions.mNumModes < 2)
            throw std::invalid_argument(
                "spectral_descriptors: num_modes must be >= 2 to span positive eigenvalues");
        EigenmodeOptions eopt;
        eopt.mNumModes = rOptions.mNumModes;
        eopt.mTolerance = rOptions.mTolerance;
        eopt.mMaxIterations = rOptions.mMaxIterations;
        EigenmodeResult eres = laplacian_eigenmodes(rMesh, eopt);
        eigvals = std::move(eres.mEigenvalues);
        n_modes = eigvals.size();
        const NDArray& modes_arr = eres.mMesh.PointData(kLaplacianModesName);
        const double* pm = modes_arr.As<double>();
        modes_flat.assign(pm, pm + n_pts * n_modes);
        result.mMesh = std::move(eres.mMesh);
    }

    // 2. Filter strictly positive eigenvalues (skip constant nullspace mode(s)).
    double lam_max = 0.0;
    for (double l : eigvals)
        if (std::isfinite(l) && l > lam_max)
            lam_max = l;
    const double cutoff = 1e-10 * lam_max;
    std::vector<std::size_t> pos_idx;
    pos_idx.reserve(n_modes);
    for (std::size_t k = 0; k < n_modes; ++k) {
        if (std::isfinite(eigvals[k]) && eigvals[k] > cutoff)
            pos_idx.push_back(k);
    }
    if (pos_idx.empty())
        throw std::invalid_argument(
            "spectral_descriptors: no positive eigenvalues available (increase num_modes)");

    const double lam_min = eigvals[pos_idx.front()];
    lam_max = eigvals[pos_idx.back()];
    result.mNumModesUsed = static_cast<std::int64_t>(pos_idx.size());

    // 3. Build log-spaced scales and per-(scale, mode) kernel weights W[s * m + j].
    const std::size_t n_scales = static_cast<std::size_t>(rOptions.mNumScales);
    const std::size_t m_used = pos_idx.size();
    result.mScales.resize(n_scales);
    std::vector<double> weights(n_scales * m_used, 0.0);

    if (rOptions.mKind == SpectralDescriptorKind::Hks) {
        // Sun, Ovsjanikov & Guibas (2009): t in [4 ln 10 / lam_max, 4 ln 10 / lam_min]
        const double c = 4.0 * std::log(10.0);
        const double t_min = c / lam_max;
        const double t_max = c / lam_min;
        const double log_min = std::log(t_min);
        const double log_max = std::log(t_max);
        for (std::size_t s = 0; s < n_scales; ++s) {
            const double alpha =
                n_scales > 1 ? static_cast<double>(s) / static_cast<double>(n_scales - 1) : 0.5;
            const double t = std::exp(log_min + alpha * (log_max - log_min));
            result.mScales[s] = t;
            double trace = 0.0;
            for (std::size_t j = 0; j < m_used; ++j) {
                const double w = std::exp(-eigvals[pos_idx[j]] * t);
                weights[s * m_used + j] = w;
                trace += w;
            }
            if (rOptions.mScaleInvariant && trace > 0.0) {
                const double inv = 1.0 / trace;
                for (std::size_t j = 0; j < m_used; ++j)
                    weights[s * m_used + j] *= inv;
            }
        }
    } else {
        // Aubry, Schlickewei & Cremers (2011): WKS on log-energy e in [log(lam_min), log(lam_max)]
        const double e_min = std::log(lam_min);
        const double e_max = std::log(lam_max);
        const double span = std::max(e_max - e_min, 1e-6);
        const double sigma =
            rOptions.mWksVariance > 0.0
                ? rOptions.mWksVariance
                : 7.0 * span / static_cast<double>(std::max<std::size_t>(n_scales, 2));
        const double inv_2s2 = 1.0 / (2.0 * sigma * sigma);
        const double lo = e_min + 2.0 * sigma * (span > 4.0 * sigma ? 1.0 : 0.0);
        const double hi = e_max - 2.0 * sigma * (span > 4.0 * sigma ? 1.0 : 0.0);
        const double e_lo = (hi > lo) ? lo : e_min;
        const double e_hi = (hi > lo) ? hi : e_max;

        std::vector<double> log_lam(m_used);
        for (std::size_t j = 0; j < m_used; ++j)
            log_lam[j] = std::log(eigvals[pos_idx[j]]);

        for (std::size_t s = 0; s < n_scales; ++s) {
            const double alpha =
                n_scales > 1 ? static_cast<double>(s) / static_cast<double>(n_scales - 1) : 0.5;
            const double e = e_lo + alpha * (e_hi - e_lo);
            result.mScales[s] = e;
            double c_e = 0.0;
            for (std::size_t j = 0; j < m_used; ++j) {
                const double diff = e - log_lam[j];
                const double w = std::exp(-diff * diff * inv_2s2);
                weights[s * m_used + j] = w;
                c_e += w;
            }
            if (c_e > 0.0) {
                const double inv = 1.0 / c_e;
                for (std::size_t j = 0; j < m_used; ++j)
                    weights[s * m_used + j] *= inv;
            }
        }
    }

    // 4. Evaluate descriptors per vertex in parallel_for (deterministic across backends).
    NDArray out(DType::Float64, {n_pts, n_scales});
    double* p_out = out.As<double>();
    parallel_for(n_pts, [&](std::size_t v) {
        const double* row = modes_flat.data() + v * n_modes;
        // Propagate NaN for points not belonging to any triangle.
        if (!std::isfinite(row[pos_idx[0]])) {
            for (std::size_t s = 0; s < n_scales; ++s)
                p_out[v * n_scales + s] = std::numeric_limits<double>::quiet_NaN();
            return;
        }
        for (std::size_t s = 0; s < n_scales; ++s) {
            const double* w_row = weights.data() + s * m_used;
            double acc = 0.0;
            for (std::size_t j = 0; j < m_used; ++j) {
                const double phi = row[pos_idx[j]];
                acc += w_row[j] * phi * phi;
            }
            p_out[v * n_scales + s] = acc;
        }
    });

    const std::string out_name =
        !rOptions.mOutput.empty()
            ? rOptions.mOutput
            : (rOptions.mKind == SpectralDescriptorKind::Hks ? kSpectralHksName : kSpectralWksName);
    result.mMesh.AddPointData(out_name, std::move(out));

    NDArray scales_arr(DType::Float64, {n_scales});
    std::copy(result.mScales.begin(), result.mScales.end(), scales_arr.As<double>());
    result.mMesh.AddFieldData(kSpectralScalesName, std::move(scales_arr));

    return result;
}

}  // namespace meshioplusplus
