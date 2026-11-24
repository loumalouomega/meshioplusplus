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
 * @file detail/streamlines.hpp
 * @brief Streamlines of a point vector field over a mesh, for `render`.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * The field is the piecewise-linear interpolation of the vector point array
 * over the mesh's cells, each split into simplices: triangles and quads
 * (two triangles) for a surface, tetrahedra, hexahedra (six), wedges (three)
 * and pyramids (two) for a volume; a higher-order cell contributes its corner
 * nodes. A volume mesh is traced in its volume, a surface mesh on its surface
 * (the field is projected onto the triangle it is on, and the point stays on
 * it), and the dimension of the highest cells found decides which.
 *
 * Everything is a pure function of the mesh and the arguments: the point
 * locator is a uniform bucket grid filled in simplex order, the seeds are the
 * centroids of simplices ranked at equal measure along a golden-ratio
 * sequence, the integrator is fixed-step RK4 in arc length, and the lines are
 * traced in parallel into per-seed buffers joined in seed order, so the
 * result does not depend on the parallel backend or the thread count.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

struct StreamlineOptions {
    /// About this many lines (at most the number of simplices).
    std::size_t mSeeds = 40;
    /// The longest a line may grow in each direction from its seed, as a
    /// fraction of the diagonal of the box around the traced cells.
    double mLength = 0.5;
};

/// Polylines, concatenated: line `l` is the vertices `[mStart[l], mStart[l+1])`.
struct Streamlines {
    std::vector<double> mXyz;  ///< 3 per vertex
    std::vector<std::size_t> mStart = {0};
    std::size_t mNumSimplices = 0;
    int mDim = 0;  ///< 2 (on a surface) or 3 (in a volume)

    std::size_t NumLines() const { return mStart.size() - 1; }
};

/// Trace the streamlines of `pVec` (3 doubles per point of `rMesh`) over the
/// cells of `rMesh`, whose points are `pXyz` (3 per point; a warp may have
/// moved them). Throws std::invalid_argument, naming the cell types, when the
/// mesh has no cell of a supported type.
Streamlines trace_streamlines(const Mesh& rMesh, const double* pXyz, const double* pVec,
                              const StreamlineOptions& rOptions);

/// Why a streamline stopped growing in one direction.
enum class StreamlineEnd : int {
    LeftDomain = 0,  ///< the next point is outside the mesh
    Stagnation = 1,  ///< the speed fell to the terminal speed or below
    MaxLength = 2,   ///< the path reached the length limit
    MaxSteps = 3,    ///< the step limit was used up
};

/// The options of the time-domain tracer behind the public `streamlines`.
struct StreamlineTraceOptions {
    /// About this many seeds when `mSeedPoints` is empty (at most the number of
    /// simplices), chosen as `trace_streamlines` chooses them.
    std::size_t mSeeds = 40;
    /// Explicit seeds, 3 doubles each; used instead of the count when not empty.
    std::vector<double> mSeedPoints;
    bool mBackward = true;
    bool mForward = true;
    /// Dormand-Prince 5(4) with step control, instead of fixed-step RK4.
    bool mAdaptive = false;
    /// The step as a length, at the fastest speed; 0 is `0.005` of the diagonal.
    double mStep = 0.0;
    /// The longest path in each direction, as a fraction of the diagonal.
    double mMaxLength = 0.5;
    /// The step limit in each direction; 0 is 20000.
    std::size_t mMaxSteps = 0;
    /// The speed at or below which a line stops; 0 is `1e-9` of the fastest.
    double mTerminalSpeed = 0.0;
    /// The adaptive integrator's tolerance as a fraction of the diagonal.
    double mTolerance = 1e-6;
};

/// One line: the points in order along the path (backward end first), the
/// velocity of the field there, and the time since the seed (negative backward).
struct StreamlineTrace {
    std::vector<double> mXyz;       ///< 3 per point
    std::vector<double> mVelocity;  ///< 3 per point
    std::vector<double> mTime;      ///< 1 per point
    std::int64_t mSeedId = 0;
    StreamlineEnd mBackwardEnd = StreamlineEnd::LeftDomain;
    StreamlineEnd mForwardEnd = StreamlineEnd::LeftDomain;
};

struct StreamlineTraces {
    std::vector<StreamlineTrace> mLines;
    std::int64_t mNumSeedsOutside = 0;  ///< explicit seeds that were not in the mesh
    int mDim = 0;                       ///< 2 (on a surface) or 3 (in a volume)
};

/// Trace by true velocity (`dx/dt = v`), so the time is the field's own time.
/// Throws std::invalid_argument, naming the cell types, as `trace_streamlines`.
StreamlineTraces trace_streamlines_timed(const Mesh& rMesh, const double* pXyz, const double* pVec,
                                         const StreamlineTraceOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus
