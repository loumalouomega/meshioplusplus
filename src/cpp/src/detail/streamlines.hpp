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

}  // namespace detail
}  // namespace meshioplusplus
