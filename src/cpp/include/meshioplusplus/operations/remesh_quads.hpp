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
#pragma once

/**
 * @file operations/remesh_quads.hpp
 * @brief Field-aligned quad remeshing of a triangle surface (roadmap 7.7.3).
 *
 * `remesh_quads` hands a triangle surface to a quad remesher and returns its
 * quads: QuadriFlow (Huang et al. 2018, MIT) or Instant Meshes (Jakob et al.
 * 2015, BSD-3), installed programs this build runs as separate processes and
 * registers by name (`quad_remesher_names()`, providers.hpp); a plugin can add
 * more. They align the quads with a smooth cross field, so a sphere comes out
 * all quads with no pole of slivers, which recombining a mesher's own
 * triangles (`mesh_cad_surface(quads=all)`) leaves. The result is a new
 * surface near the input, not on its points: `mHausdorffToSurface` says how
 * near. No program runs in WebAssembly.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

struct RemeshQuadsOptions {
    /// `auto` (QuadriFlow, else Instant Meshes) or a quad remesher's name.
    std::string mBackend = "auto";
    /// The number of faces asked for; 0 takes half the input's triangles.
    std::int64_t mTargetFaces = 0;
    /// Align with and keep the sharp edges.
    bool mSharp = false;
    /// Keep the boundary of an open surface.
    bool mBoundary = true;
    /// Let the face size follow the curvature (QuadriFlow).
    bool mAdaptive = false;
    /// Allow triangles where quads fit badly (Instant Meshes' quad-dominant mode).
    bool mDominant = false;
    /// The program's random seed, so a run repeats.
    std::int32_t mSeed = 0;
    /// Refuse a target of more faces.
    std::int64_t mMaxCells = 50000000;
};

struct RemeshQuadsResult {
    /// A `quad` block, and a `triangle` block when the remesher left any.
    Mesh mMesh;
    std::int64_t mNumQuads = 0;
    std::int64_t mNumTriangles = 0;
    std::int64_t mNumPoints = 0;
    /// The least and mean scaled Jacobian of the quads (Verdict's, in [-1, 1]).
    double mMinScaledJacobian = 0.0;
    double mMeanScaledJacobian = 0.0;
    /// The two-sided distance between the result and the input surface.
    double mHausdorffToSurface = 0.0;
    std::string mBackend;  ///< the remesher that ran
    std::vector<std::string> mWarnings;
};

/// The backend names `remesh_quads` takes besides `auto`: every registered quad
/// remesher (the installed `quadriflow` and `instant_meshes`), sorted.
MESHIOPLUSPLUS_API std::vector<std::string> remesh_quads_backends();

/**
 * @brief Remesh the triangle surface `rSurface` into field-aligned quads.
 * @throws std::invalid_argument for a negative target or one past `mMaxCells`,
 *         or a surface that is not linear triangles (lines and vertices are
 *         skipped).
 * @throws Unsupported when no quad remesher is installed, or none of that name.
 * @throws std::runtime_error when the remesher fails or writes a bad mesh.
 */
MESHIOPLUSPLUS_API RemeshQuadsResult remesh_quads(const Mesh& rSurface,
                                                  const RemeshQuadsOptions& rOptions = {});

}  // namespace meshioplusplus
