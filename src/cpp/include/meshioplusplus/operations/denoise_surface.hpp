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
 * @file operations/denoise_surface.hpp
 * @brief Feature-preserving surface denoising via two-stage bilateral and
 * guided normal filtering (roadmap 9.12; Zheng et al. 2011; Zhang et al. 2015;
 * Sun et al. 2007).
 *
 * Stage 1 filters face normals iteratively in deterministic Jacobi
 * `parallel_for` passes using spatial and normal-similarity Gaussian weights.
 * Stage 2 updates vertex positions to align with the filtered normals, with an
 * optional facet fold-over guard.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Normal-filtering method for `denoise_surface`.
enum class DenoiseMethod : std::int32_t {
    /// Local bilateral normal filter (Zheng et al. 2011).
    Bilateral = 0,
    /// Patch-guided bilateral normal filter (Zhang et al. 2015).
    Guided = 1,
};

/// Parse `"bilateral"` or `"guided"`.
MESHIOPLUSPLUS_API DenoiseMethod denoise_method_from_name(const std::string& rName);

/// Options for `denoise_surface`.
struct DenoiseSurfaceOptions {
    /// Filtering method (`Bilateral` or `Guided`).
    DenoiseMethod mMethod = DenoiseMethod::Bilateral;
    /// Number of Jacobi normal-filtering iterations.
    std::int64_t mNormalIterations = 10;
    /// Number of vertex-position update iterations.
    std::int64_t mVertexIterations = 10;
    /// Spatial Gaussian sigma over face centroid distances; 0 selects mean adjacent centroid
    /// distance.
    double mSigmaSpatial = 0.0;
    /// Normal-similarity Gaussian sigma (in unit-normal Euclidean distance, typically 0.2..0.5).
    double mSigmaNormal = 0.35;
    /// Keep boundary vertices fixed.
    bool mFixBoundary = true;
    /// Skip vertex moves that would flip an incident face normal against its filtered normal.
    bool mGuardFoldover = true;
    /// Optional named Cell region to denoise (empty = whole surface).
    std::string mRegion;
};

/// Result of `denoise_surface`.
struct DenoiseSurfaceResult {
    /// The denoised mesh.
    Mesh mMesh;
    /// Number of vertices whose position changed.
    std::int64_t mNumNodesMoved = 0;
    /// Maximum Euclidean displacement of any vertex.
    double mMaxDisplacement = 0.0;
    /// Number of per-iteration vertex updates skipped by the fold-over guard.
    std::int64_t mNumSkippedFoldover = 0;
    /// Spatial sigma actually used (helpful when `mSigmaSpatial == 0`).
    double mSigmaSpatialUsed = 0.0;
};

/**
 * @brief Denoise a surface mesh while preserving sharp dihedral features.
 */
MESHIOPLUSPLUS_API DenoiseSurfaceResult denoise_surface(const Mesh& rMesh,
                                                        const DenoiseSurfaceOptions& rOptions = {});

}  // namespace meshioplusplus
