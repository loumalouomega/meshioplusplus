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
 * @file operations/streamlines.hpp
 * @brief Streamlines of a point vector field, as a mesh of line segments.
 *
 * The field is the piecewise-linear interpolation of the 2- or 3-component
 * point array `mArray` over the mesh's cells (each split into simplices: a
 * volume mesh is traced in its volume, a surface mesh on its surface, the
 * field projected onto the triangle it is on). Each line is traced from its
 * seed by the **field's own velocity** (`dx/dt = v`), forward, backward or both,
 * with fixed-step RK4 or adaptive Dormand-Prince 5(4), so `IntegrationTime` is
 * the field's time. The result is a new mesh: the points of every line, one
 * `line` segment between consecutive points, and the point data
 * `IntegrationTime` (negative backward of the seed), `Speed`, `SeedId` and the
 * field itself under its own name; the cell data `SeedId` per segment.
 *
 * Seeds are `mSeeds` simplex centroids ranked at equal measure along a
 * golden-ratio sequence (the choice `render`'s streamlines make), or explicit
 * points. A line stops in each direction because it left the mesh, the speed
 * fell to `mTerminalSpeed`, the path reached `mMaxLength` or the step limit was
 * used up; the reasons are in the result. Lines are traced in parallel into
 * per-seed buffers joined in seed order, so the result does not depend on the
 * backend or the thread count. `render`'s streamlines keep their own fixed
 * arc-length tracer.
 *
 * C++ core only, with no numpy twin: the adaptive step control branches on
 * rounded error estimates, which a second implementation could decide
 * differently.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Which way from each seed.
enum class StreamlineDirection : std::int32_t { Both = 0, Forward = 1, Backward = 2 };

/// The integrator.
enum class StreamlineIntegrator : std::int32_t { RK4 = 0, RK45 = 1 };

/// Why a line stopped growing in one direction (`StreamlinesResult::mEnd*`).
enum class StreamlineStop : std::int32_t {
    LeftDomain = 0,
    Stagnation = 1,
    MaxLength = 2,
    MaxSteps = 3,
};

/// Point data written by `streamlines` (the field keeps its own name).
inline constexpr const char* kStreamlineTimeName = "IntegrationTime";
inline constexpr const char* kStreamlineSpeedName = "Speed";
inline constexpr const char* kStreamlineSeedIdName = "SeedId";

/// Options for `streamlines`.
struct StreamlinesOptions {
    /// The point array to follow: 2 or 3 components (a 2-component field is
    /// `(u, v, 0)`). Required.
    std::string mArray;
    /// About this many seeds, used when `mSeedPoints` is empty.
    std::int64_t mSeeds = 40;
    /// Explicit seeds, 3 coordinates each.
    std::vector<double> mSeedPoints;
    StreamlineDirection mDirection = StreamlineDirection::Both;
    StreamlineIntegrator mIntegrator = StreamlineIntegrator::RK4;
    /// The step as a length at the fastest speed; 0 is `0.005` of the diagonal.
    double mStep = 0.0;
    /// The longest path in each direction, as a fraction of the diagonal.
    double mMaxLength = 0.5;
    /// The step limit in each direction; 0 is 20000.
    std::int64_t mMaxSteps = 0;
    /// The speed at or below which a line stops; 0 is `1e-9` of the fastest.
    double mTerminalSpeed = 0.0;
    /// The adaptive integrator's tolerance as a fraction of the diagonal.
    double mTolerance = 1e-6;
};

/// The lines, and how each one ended.
struct StreamlinesResult {
    /// Points and `line` segments, with the data described above.
    Mesh mMesh;
    std::int64_t mNumLines = 0;
    /// Explicit seeds that were not inside the mesh (no line is made for them).
    std::int64_t mNumSeedsOutside = 0;
    /// Per line: the seed (an index into the seeds as listed) and the
    /// `StreamlineStop` of the backward and the forward end.
    std::vector<std::int64_t> mSeedIds;
    std::vector<std::int32_t> mEndBackward;
    std::vector<std::int32_t> mEndForward;
    /// 2 (traced on a surface) or 3 (in a volume).
    std::int32_t mDim = 0;
};

/**
 * @brief Traces streamlines of `mArray` over @p rMesh.
 * @throws std::invalid_argument on an unknown array or one that is not 2- or
 *         3-component, on a mesh with no triangle, quad, tetra, hexahedron,
 *         wedge or pyramid cell, or on invalid options.
 */
MESHIOPLUSPLUS_API StreamlinesResult streamlines(const Mesh& rMesh,
                                                 const StreamlinesOptions& rOptions);

}  // namespace meshioplusplus
