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
#pragma once

/**
 * @file operations/spectral_descriptors.hpp
 * @brief Intrinsic, isometry-invariant spectral shape descriptors (roadmap 15.6):
 * Heat Kernel Signature (HKS; Sun, Ovsjanikov and Guibas 2009) and Wave Kernel
 * Signature (WKS; Aubry, Schlickewei and Cremers 2011), built directly on the
 * Laplace-Beltrami eigenpairs of `laplacian_eigenmodes`.
 *
 * When `mReuseEigenmodes` is true and the mesh already carries
 * `point_data["laplacian:modes"]` and `field_data["laplacian:eigenvalues"]`,
 * those eigenpairs are reused without re-running the eigensolver.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Default point-data array name for HKS (`(n, num_scales)` Float64).
inline constexpr const char* kSpectralHksName = "spectral:hks";
/// Default point-data array name for WKS (`(n, num_scales)` Float64).
inline constexpr const char* kSpectralWksName = "spectral:wks";
/// Field-data array name for the log-spaced scales (`(num_scales,)` Float64).
inline constexpr const char* kSpectralScalesName = "spectral:scales";

/// Which spectral shape descriptor to compute.
enum class SpectralDescriptorKind : std::int32_t {
    Hks = 0,
    Wks = 1,
};

/// Parse `"hks"` or `"wks"` (case-insensitive).
MESHIOPLUSPLUS_API SpectralDescriptorKind spectral_descriptor_kind_from_name(
    const std::string& rName);

/// Options for `spectral_descriptors`.
struct SpectralDescriptorsOptions {
    /// Descriptor family (`Hks` or `Wks`).
    SpectralDescriptorKind mKind = SpectralDescriptorKind::Hks;
    /// Number of Laplace-Beltrami eigenmodes to compute when not reusing existing ones.
    std::int64_t mNumModes = 32;
    /// Number of log-spaced time (HKS) or log-energy (WKS) scales.
    std::int64_t mNumScales = 16;
    /// Normalize each scale by its spectral trace (HKS) or scale-invariant log-derivative.
    bool mScaleInvariant = false;
    /// Reuse `laplacian:modes` / `laplacian:eigenvalues` if already attached to the mesh.
    bool mReuseEigenmodes = true;
    /// WKS Gaussian standard deviation in log-energy; 0 selects `7 * (e_max - e_min) / num_scales`.
    double mWksVariance = 0.0;
    /// Tolerance forwarded to `laplacian_eigenmodes` when eigenmodes are computed.
    double mTolerance = 1e-8;
    /// Max iterations forwarded to `laplacian_eigenmodes`; 0 uses its default.
    std::int64_t mMaxIterations = 0;
    /// Output point-data array name; empty defaults to `spectral:hks` or `spectral:wks`.
    std::string mOutput;
};

/// Result of `spectral_descriptors`.
struct SpectralDescriptorsResult {
    /// Output mesh carrying the descriptor point array and `spectral:scales` field array.
    Mesh mMesh;
    /// The log-spaced time (HKS) or log-energy (WKS) scales evaluated.
    std::vector<double> mScales;
    /// Number of positive eigenmodes used in the descriptor sum.
    std::int64_t mNumModesUsed = 0;
    /// True when precomputed `laplacian:modes` / `laplacian:eigenvalues` were reused.
    bool mReusedEigenmodes = false;
};

/**
 * @brief Compute HKS or WKS spectral shape descriptors on a triangle surface.
 * @throws std::invalid_argument on invalid options or unsupported mesh topology.
 */
MESHIOPLUSPLUS_API SpectralDescriptorsResult
spectral_descriptors(const Mesh& rMesh, const SpectralDescriptorsOptions& rOptions = {});

}  // namespace meshioplusplus
