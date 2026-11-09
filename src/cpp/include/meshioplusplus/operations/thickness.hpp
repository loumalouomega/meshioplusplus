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
 * @file operations/thickness.hpp
 * @brief Wall thickness and clearance per facet, from a ray along each
 * facet's normal to the opposite surface (roadmap 3.1.4).
 *
 * ### The ray
 *
 * From each 2-D cell's area-weighted centroid, along its unit normal (the
 * sum of its fan triangles' normals, so a quad or polygon has one), a ray is
 * cast through a bounding-volume hierarchy of the target surface's triangles
 * with the core's one ray–triangle kernel (Möller and Trumbore, JGT 2(1),
 * 1997). The surface must be consistently oriented with outward normals
 * (`repair` orients a closed shell):
 *
 *  - **Inward** (`ThicknessDirection::Inward`, the default for one surface):
 *    the ray runs against the normal and stops at the first triangle it
 *    leaves through, one whose normal points along the ray -- the far side
 *    of the wall. The distance is the wall thickness.
 *  - **Outward** (`Outward`, the default for two surfaces): the ray runs
 *    along the normal and stops at the first triangle it enters, one whose
 *    normal faces it -- the next wall across a gap. The distance is the
 *    clearance, to the same surface or to a second one.
 *
 * A cell's own triangles are never hit. A ray that finds nothing within
 * `mMaxDistance` is a miss: its distance is NaN and its hit cell -1.
 *
 * ### Beyond one ray (roadmap 7.5)
 *
 *  - **`ThicknessMethod::Sphere`**: from the ray's midpoint, the distance to
 *    the nearest target triangle outside the facet's one-ring (the triangles
 *    sharing a vertex with it), doubled, and never more than the ray: exact
 *    between parallel walls like the ray, and where the faces diverge the
 *    perpendicular distance to the opposite face rather than the ray's
 *    1/cos longer one. A miss of the ray is a miss here too.
 *  - **`ThicknessMethod::Cone`**: the shortest hit of a fixed fan of
 *    `mConeRays` rays within `mConeAngleDeg` of the normal (the central ray,
 *    then a Fibonacci spiral in the facet's frame), which finds the narrow
 *    spot of a chamfer or a tapered rib.
 *  - **`mTwoSided`**: scan winding is unreliable, so rays run both ways along
 *    the normal line and either face of every triangle counts; the shorter
 *    hit wins. With `mMaxDistance` 0 the ray is then a fifth of the diagonal.
 *  - **`mMaxNormalDeviationDeg`**: a triangle more edge-on to the ray than
 *    this is stepped over (a cone ray grazing a convex edge); 90, the
 *    default, accepts any triangle facing the right way.
 *
 * The report adds the 2nd and 95th percentiles of the hit distances, the
 * range a colour map is set from.
 *
 * ### What comes back
 *
 * The input, bit-identical, with the distance and the cell it reached as
 * cell data (cells that are not tested, such as lines, read as misses), and
 * a report: the least, largest and mean distance over the hits, where the
 * least was measured, the misses and, against a threshold, how many cells
 * fall below it.
 *
 * Deterministic: each ray's hit is the nearest accepted one with ties to the
 * lower triangle, whatever order the tree is walked in, and the report is
 * summed in cell order, so the result is identical across parallel backends
 * and thread counts.
 *
 * @see doc/thickness.md
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Cell data (Float64): the distance along the ray; NaN for a miss.
inline constexpr const char* kThicknessDistanceName = "thickness:distance";
/// Cell data (Int64): the global cell of the target surface the ray reached;
/// -1 for a miss.
inline constexpr const char* kThicknessHitCellName = "thickness:hit_cell";

/// Which way the ray runs from a facet.
enum class ThicknessDirection : std::int32_t {
    Auto = 0,     ///< inward for one surface, outward for two
    Inward = 1,   ///< against the normal, to the far side of the wall
    Outward = 2,  ///< along the normal, across a gap
};

/// How `thickness` measures from a facet.
enum class ThicknessMethod : std::int32_t {
    Ray = 0,     ///< one ray along the normal
    Sphere = 1,  ///< twice the distance from the ray's midpoint to the nearest other surface
    Cone = 2,    ///< the shortest of a fan of rays within the cone angle
};

/// Options for `thickness`.
struct ThicknessOptions {
    ThicknessDirection mDirection = ThicknessDirection::Auto;
    /// Longest ray; 0 takes the diagonal of both surfaces' bounding box. Must
    /// be finite and non-negative.
    double mMaxDistance = 0.0;
    /// When positive, count the cells whose distance is below it.
    double mThreshold = 0.0;
    /// Restrict the measured cells to this named `Cell` region.
    std::string mRegion;
    /// Restrict the target surface (the second mesh, else the mesh itself) to
    /// this named `Cell` region; empty takes all of it.
    std::string mRegionTarget;
    /// Attach `thickness:distance` and `thickness:hit_cell`.
    bool mAttach = true;
    /// Rays both ways along the normal, either face of a triangle counting.
    bool mTwoSided = false;
    ThicknessMethod mMethod = ThicknessMethod::Ray;
    /// `Cone`: the half-angle in degrees, in (0, 90).
    double mConeAngleDeg = 15.0;
    /// `Cone`: the rays of the fan, the central one included; at least 1.
    std::int64_t mConeRays = 16;
    /// Step over triangles more edge-on to the ray than this, in degrees in
    /// (0, 90]; 90 accepts any triangle facing the right way.
    double mMaxNormalDeviationDeg = 90.0;
};

/// The measured mesh and the summary.
struct ThicknessResult {
    /// The input, bit-identical, with the two arrays when asked for.
    Mesh mMesh;
    /// The direction the rays ran, `Inward` or `Outward`.
    ThicknessDirection mDirection = ThicknessDirection::Inward;
    std::int64_t mNumCells = 0;   ///< cells measured
    std::int64_t mNumHit = 0;     ///< of those, rays that reached a wall
    std::int64_t mNumMissed = 0;  ///< rays that did not
    std::int64_t mNumBelow = 0;   ///< hits below `mThreshold` (when positive)
    double mMin = 0.0;            ///< least distance over the hits (0 when none)
    double mMax = 0.0;            ///< largest distance over the hits
    double mMean = 0.0;           ///< mean distance over the hits
    std::int64_t mMinCell = -1;   ///< the global cell of the least distance
    std::array<double, 3> mMinPoint{{0.0, 0.0, 0.0}};  ///< where its ray started
    double mMaxDistance = 0.0;                         ///< the ray length used
    std::vector<std::string> mWarnings;
    ThicknessMethod mMethod = ThicknessMethod::Ray;  ///< the method used
    /// The 2nd and 95th percentiles of the hit distances (0 when none).
    double mP2 = 0.0;
    double mP95 = 0.0;
};

/** @brief `"ray"`, `"sphere"` or `"cone"`; anything else throws naming them. */
MESHIOPLUSPLUS_API ThicknessMethod thickness_method_from_name(const std::string& rName);

/**
 * @brief Measure the wall thickness (or, outward, the clearance) of the
 * surface @p rMesh against itself.
 * @param rMesh a consistently oriented surface mesh (triangles, quads,
 *        polygons; other cells are not measured).
 * @throws std::invalid_argument on a volume, polyhedron or higher-order block
 *         (run `extract_surface` or `linearize` first), an unknown region or a
 *         negative or non-finite distance.
 */
MESHIOPLUSPLUS_API ThicknessResult thickness(const Mesh& rMesh,
                                             const ThicknessOptions& rOptions = {});

/**
 * @brief Measure the clearance (or, inward, the depth) from each facet of
 * @p rMesh to the surface @p rTarget.
 * @throws std::invalid_argument as the one-surface form, for either input.
 */
MESHIOPLUSPLUS_API ThicknessResult thickness(const Mesh& rMesh, const Mesh& rTarget,
                                             const ThicknessOptions& rOptions = {});

}  // namespace meshioplusplus
