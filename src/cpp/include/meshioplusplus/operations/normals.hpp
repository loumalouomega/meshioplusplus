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
 * @file operations/normals.hpp
 * @brief Point and cell normals of a surface mesh, optionally splitting
 * vertices at creases.
 *
 * A *vertex* normal is a property of a smooth patch, not of a position: at the
 * edge of a cube the one position has three normals. With `mSplit` off the
 * operation returns one normal per point, the angle- (or area-) weighted mean of
 * the incident faces, and smooths the crease over. With `mSplit` on it does what
 * a renderer needs: the corners around a vertex are grouped into smooth fans
 * (see `detail/surface_normals.hpp` for how) and every fan beyond the first gets
 * its own copy of the point, so each point carries exactly one normal.
 *
 * The split layout is the one `repair` uses for its bowtie copies: the original
 * points keep their indices and the copies are appended, so
 * `point_data` gathers by row and `normals:parent_point` (opt-in) names each
 * copy's source. A copy joins its source's Point regions. Cell numbering is
 * untouched, so Cell and Side regions and cell data ride through verbatim.
 *
 * ### What it will not do
 *
 * It **never reorients**. Two triangles that disagree about which side is out
 * cannot both be right, and averaging them gives a wrong normal that looks
 * plausible. `mQuality.mInconsistentPairs` reports the count, a split always
 * cuts at such an edge, and the fix is `repair(mesh, {.mFixOrientation = true})`
 * first -- the same contract `compute_curvature` has.
 *
 * It works on a *surface*. A volume block is refused by name, pointing at
 * `extract_surface`, rather than being silently reduced to its skin.
 *
 * ### The output name
 *
 * The arrays are called `normals` -- the name the PCD and XYZ readers and
 * writers already use -- so `compute_normals` followed by a write to `.pcd` or
 * `.xyz` emits the normal columns with no further step. An existing array of
 * that name is replaced.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/sdf.hpp"

namespace meshioplusplus {

/// Point data `(n, 3)` and cell data `(cells, 3)`: the unit normals, Float64.
inline constexpr const char* kNormalsName = "normals";
/// Point data: for each output point, the input point it came from. Opt-in.
inline constexpr const char* kNormalsParentPointName = "normals:parent_point";

/// What `compute_normals` should compute and attach.
struct NormalsOptions {
    /// Attach point normals.
    bool mPointNormals = true;
    /// Attach cell normals (the unit vector area of each cell).
    bool mCellNormals = false;
    /// How incident faces are weighted into a point normal. Cell normals do not
    /// depend on it.
    SdfPseudonormalWeight mWeight = SdfPseudonormalWeight::Angle;
    /// Duplicate points where the surface creases by more than `mSplitAngle`.
    bool mSplit = false;
    /// The largest dihedral angle, in degrees, still treated as smooth. Used
    /// only when `mSplit` is set; must lie in `[0, 180]`.
    double mSplitAngle = 30.0;
    /// Attach `normals:parent_point`.
    bool mRecordParentIds = false;
    /// Restrict to this named `Cell` region; empty takes every surface cell.
    std::string mRegion;
};

/// What `compute_normals` computed, and what it found on the way.
struct NormalsResult {
    /// The input with the requested arrays attached (and, with `mSplit`, the
    /// split points appended).
    Mesh mMesh;
    /// The INPUT surface's defect counts. `mInconsistentPairs != 0` means some
    /// normals are averaged across faces that disagree about "out".
    SurfaceQuality mQuality;
    /// Points no selected surface triangle touches; their normal is NaN.
    std::int64_t mNumIsolated = 0;
    /// Touched points whose incident faces sum to nothing; their normal is NaN.
    std::int64_t mNumUndefined = 0;
    /// Triangles with no area, which contribute no direction.
    std::int64_t mNumDegenerate = 0;
    /// Input points that received at least one copy.
    std::int64_t mNumSplitPoints = 0;
    /// Points appended to the mesh.
    std::int64_t mNumAddedPoints = 0;
};

/**
 * @brief Point and/or cell normals of a surface mesh.
 *
 * Triangles come from `detail::build_triangle_soup`, so quads and polygons are
 * fanned on the diagonal `convert_cells(Simplexify)` uses; a polygon's cell
 * normal is the sum of its fan (Newell's normal), which does not depend on which
 * corner the fan starts from. Lines and vertices are ignored and get NaN.
 *
 * @param rMesh a surface mesh (2-D cells in 2-D or 3-D space, or a mix with
 *        lines and vertices).
 * @param rOptions what to compute; see `NormalsOptions`.
 * @throws std::invalid_argument on a volume or higher-order block, an unknown
 *         region name, or a split angle outside `[0, 180]`.
 */
MESHIOPLUSPLUS_API NormalsResult compute_normals(const Mesh& rMesh,
                                                 const NormalsOptions& rOptions = {});

}  // namespace meshioplusplus
