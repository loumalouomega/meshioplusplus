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
 * @file operations/sdf.hpp
 * @brief Distance to a surface: how far a point is from a triangle skin, and
 * which side of it the point is on.
 *
 * This is the primitive collision detection, offsetting, inside/outside queries
 * and voxel-style preprocessing all reduce to, and the one spatial query
 * meshio++ did not have. `extract_surface` gives you the skin, `slice` cuts it
 * and `isosurface` contours a field on it -- but nothing answered "how far is
 * this point from the surface".
 *
 * ### What lives here
 *
 * - `sample_distance` -- distances at arbitrary query points, the batch form.
 * - `distance_to_surface` -- the same, attached to a query mesh as `point_data`.
 * - `surface_watertight_check` -- what is wrong with this skin, in numbers.
 * - `compute_sdf` -- the umbrella that generates its own grid, dense or adaptive.
 *
 * `compute_sdf`'s types were **complete in this header from v9.24.0**, one
 * release before its body, deliberately:
 * `SdfOptions` embeds `SurfaceDistanceOptions` by value, so adding a member to
 * the inner struct later would shift the outer struct's tail under a consumer
 * compiled against the older header -- a silent Tier A ABI break, and exactly the
 * `Pipeline`-embeds-`PipelineInput` trap `doc/abi_reviews.md` records for
 * v9.12.0. Shipping the full layout up front means every later release changes
 * function declarations only, which is unambiguously additive.
 *
 * ### The sign is the hard part
 *
 * An unsigned distance is a minimisation and has no interesting failure mode. A
 * *signed* distance needs to know which side of the surface the query is on, and
 * the standard mistake -- taking the sign from the nearest **triangle's** normal
 * -- is exactly right on convex geometry and exactly wrong on the concave side
 * of every crease. The fix is the angle-weighted **pseudonormal** (Baerentzen &
 * Aanaes): the sign comes from the normal of the nearest *feature*, which is the
 * triangle only when the closest point lies in its interior, and is a blend of
 * the incident faces' normals when it lies on an edge or at a vertex.
 *
 * For a skin that is not watertight -- which real STLs frequently are not -- the
 * pseudonormal is not merely inaccurate, it is undefined, so `SdfSign` also
 * offers the generalized winding number, which is robust to holes,
 * self-intersection and inconsistent winding at the cost of being O(number of
 * triangles) per query with no acceleration structure that helps.
 *
 * @see doc/sdf.md
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// The Float64 signed distance array. Negative inside, by the usual convention.
inline constexpr const char* kSdfDistanceName = "sdf:distance";

/// The Int64 0/1 array marking values that were actually computed rather than
/// clamped to the band limit. Attached whenever a band is in force, and
/// **not optional**: a clamped value is byte-indistinguishable from a computed
/// one, so without this array a caller cannot tell the difference.
inline constexpr const char* kSdfBandName = "sdf:band";

/// The Int64 0/1 array marking points inside the surface (opt-in).
inline constexpr const char* kSdfInsideName = "sdf:inside";

/// The Int64 array naming the nearest **input** cell per query (opt-in).
inline constexpr const char* kSdfClosestCellName = "sdf:closest_cell";

/// How to decide which side of the surface a query point is on.
enum class SdfSign {
    /// No sign at all: the distance is the unsigned nearest-point distance.
    /// The only mode that is meaningful on an open sheet.
    Unsigned = 0,
    /// The angle-weighted pseudonormal of the nearest *feature* (face, edge or
    /// vertex). Exact for a watertight, consistently wound surface, and costs
    /// nothing beyond the nearest-triangle query that already ran.
    Pseudonormal = 1,
    /// The generalized winding number. Robust to holes, self-intersections and
    /// per-component orientation flips, but **O(number of triangles) per query**
    /// -- there is no acceleration structure short of a fast-multipole
    /// expansion, which is a project in its own right. Pair it with a band.
    WindingNumber = 2,
};

/// How the incident faces are weighted when building a vertex pseudonormal.
enum class SdfPseudonormalWeight {
    /// Weight each incident face by the angle it subtends at the vertex. The
    /// geometrically correct choice (Baerentzen & Aanaes) and the default.
    Angle = 0,
    /// Weight by face area. Not exact at a *vertex*, but exact at edges and
    /// faces, which is where nearly every nearest point on a fine mesh lands.
    /// It exists because it is free of `acos`, and therefore has a bit-exact
    /// numpy twin -- see `doc/sdf.md` on why that is worth a second mode.
    Area = 1,
};

/// Where a distance is evaluated on a mesh.
enum class SdfLocation {
    Corner = 0,  ///< At the points, as `point_data`.
    Center = 1,  ///< At the cell centroids, as `cell_data`.
};

/// What `compute_sdf` generates to carry the field.
enum class SdfStructure {
    Voxel = 0,  ///< A dense uniform lattice over the (padded) bounding box.
    /// Adaptive: a coarse root lattice, refined only near the surface. The
    /// output is **1-irregular** (it has hanging nodes) -- see `refine.hpp`'s
    /// `RefineClosure::Balanced`, which is the mechanism.
    Octree = 1,
};

/// What to do when the surface turns out not to be watertight.
enum class SdfWatertightCheck {
    Off = 0,    ///< Do not look.
    Warn = 1,   ///< Log the counts by name and carry on. The default.
    Error = 2,  ///< Throw, naming the counts.
};

/// Parse an `SdfSign` name (`unsigned`, `pseudonormal`, `winding-number`).
MESHIOPLUSPLUS_API SdfSign sdf_sign_from_name(const std::string& rName);
/// Parse an `SdfPseudonormalWeight` name (`angle`, `area`).
MESHIOPLUSPLUS_API SdfPseudonormalWeight sdf_weight_from_name(const std::string& rName);
/// Parse an `SdfLocation` name (`corner`/`point`, `center`/`cell`).
MESHIOPLUSPLUS_API SdfLocation sdf_location_from_name(const std::string& rName);
/// Parse an `SdfStructure` name (`voxel`, `octree`).
MESHIOPLUSPLUS_API SdfStructure sdf_structure_from_name(const std::string& rName);
/// Parse an `SdfWatertightCheck` name (`off`, `warn`, `error`).
MESHIOPLUSPLUS_API SdfWatertightCheck sdf_watertight_check_from_name(const std::string& rName);

/**
 * @brief What is wrong with a surface, in numbers rather than a bare bool.
 *
 * "Your STL has 12 boundary edges" is actionable; "not watertight" is not. All
 * four edge counts are reported separately because they need different fixes.
 */
struct SurfaceQuality {
    /// Edges used by exactly one triangle: the surface has a hole or is a sheet.
    std::int64_t mBoundaryEdges = 0;
    /// Edges used by three or more triangles: the surface is non-manifold.
    std::int64_t mNonManifoldEdges = 0;
    /// Edges whose two triangles wind the same way, i.e. disagree about which
    /// side is out. A sign computed from these is meaningless locally.
    std::int64_t mInconsistentPairs = 0;
    /// Triangles with zero area. They carry no normal and are skipped.
    std::int64_t mDegenerateTriangles = 0;
    /// True when all four counts are zero.
    bool mWatertight = false;
};

/// Options shared by every distance query.
struct SurfaceDistanceOptions {
    /// How to sign the distance.
    SdfSign mSign = SdfSign::Pseudonormal;
    /// How to weight incident faces at a vertex (`Pseudonormal` only).
    SdfPseudonormalWeight mWeight = SdfPseudonormalWeight::Angle;
    /// Where to evaluate, for the mesh-taking entry points.
    SdfLocation mLocation = SdfLocation::Corner;
    /// Distances beyond this are clamped to it and marked in `sdf:band`.
    /// Non-positive means the full field, computed everywhere.
    double mBand = 0.0;
    /// Attach `sdf:closest_cell`.
    bool mRecordClosestCell = false;
    /// Attach `sdf:inside`.
    bool mRecordInside = false;
    /// What to do about a surface that is not watertight.
    SdfWatertightCheck mWatertightCheck = SdfWatertightCheck::Warn;
    /// Restrict the surface to a named `Cell` region; empty means all of it.
    std::string mSurfaceRegion;
    /// Bucket size of the search grid; 0 derives one from the triangle sizes.
    ///
    /// This is a tuning knob, but it is public mainly so that the accelerator
    /// can be *proved* not to affect the answer: the same query at several cell
    /// sizes must give byte-identical distances and nearest-cell ids. See
    /// `doc/sdf.md`.
    double mGridCellSize = 0.0;
    /// Refuse a `WindingNumber` query whose `n_query * n_tri` exceeds this,
    /// naming the option, rather than silently running for an hour.
    double mMaxWindingWork = 2.0e9;
};

/// The result of `distance_to_surface`: the query mesh plus what was found.
struct SurfaceDistanceResult {
    /// The query mesh with the distance arrays attached.
    Mesh mMesh;
    /// What was wrong with the surface (all zeros when `mWatertightCheck` is Off).
    SurfaceQuality mQuality;
    /// How many queries were clamped to the band rather than computed.
    std::int64_t mNumBanded = 0;
};

/**
 * @brief What is wrong with @p rSurface, without computing any distances.
 * @param rSurface a triangle (or triangulatable) surface mesh.
 * @return the four edge counts and the resulting verdict.
 * @throws std::invalid_argument on a volume or polyhedron block.
 */
MESHIOPLUSPLUS_API SurfaceQuality surface_watertight_check(const Mesh& rSurface);

/**
 * @brief Distances from arbitrary points to a surface.
 * @param rSurface the surface to measure against.
 * @param rPoints Float64 `(n, 2)` or `(n, 3)` query coordinates; a 2-D array is
 *        z-padded, matching `detail::read_point`.
 * @param rOptions sign, band and accelerator settings.
 * @return a Float64 `(n,)` array of signed distances.
 * @throws std::invalid_argument on a mis-shaped query array, a surface with no
 *         triangles, or a `WindingNumber` request exceeding `mMaxWindingWork`.
 */
MESHIOPLUSPLUS_API NDArray sample_distance(const Mesh& rSurface, const NDArray& rPoints,
                                           const SurfaceDistanceOptions& rOptions = {});

/**
 * @brief Attach distances from a query mesh's points (or cell centres) to a
 *        surface.
 * @param rQuery the mesh to annotate; its geometry is copied unchanged.
 * @param rSurface the surface to measure against.
 * @param rOptions sign, location, band and accelerator settings.
 * @return the annotated mesh plus the surface verdict and the banded count.
 */
MESHIOPLUSPLUS_API SurfaceDistanceResult distance_to_surface(
    const Mesh& rQuery, const Mesh& rSurface, const SurfaceDistanceOptions& rOptions = {});

// --- the umbrella operation (types complete; implementation lands later) -----

/// `field_data` keys describing the generated grid. They are numeric, so they
/// cross every binding -- but note **no file format persists arbitrary
/// `field_data`**, so a written-and-reread grid loses them and they are
/// recovered from the geometry instead. See `doc/sdf.md`.
inline constexpr const char* kSdfOriginName = "sdf:origin";
inline constexpr const char* kSdfSpacingName = "sdf:spacing";
inline constexpr const char* kSdfDimsName = "sdf:dims";
inline constexpr const char* kSdfBoundsName = "sdf:bounds";
inline constexpr const char* kSdfMaxDepthName = "sdf:max_depth";
inline constexpr const char* kSdfStructureName = "sdf:structure";

/// Options for `compute_sdf`: generate a grid, then fill it with distances.
struct SdfOptions {
    /// Which structure to generate.
    SdfStructure mStructure = SdfStructure::Voxel;
    /// Cell counts per axis; unset derives them from `mCellSize`.
    std::optional<std::array<std::int64_t, 3>> mResolution;
    /// Cell size; unset derives it from `mResolution`.
    std::optional<double> mCellSize;
    /// Explicit bounds `{lo[3], hi[3]}`; unset uses the surface's bounding box.
    std::optional<std::array<double, 6>> mBounds;
    /// Padding added to every side, in world units.
    double mPadding = 0.0;
    /// Padding added to every side, as a fraction of the bounding-box diagonal.
    /// An SDF that stops at the surface is not much use, hence the non-zero
    /// default.
    double mPaddingRelative = 0.1;
    /// Octree: cell count per axis of the root lattice. `mResolution` and
    /// `mCellSize` size a `Voxel` grid and are an **error** with `Octree`, whose
    /// finest cell is `root cell / 2^depth` and is therefore already determined.
    std::int64_t mRootResolution = 8;
    /// Octree: how many refinement passes. Each halves the cell size in the band.
    std::int64_t mMaxDepth = 4;
    /// Octree: refine a cell while `|distance| <= mBandCells * cell diagonal`.
    /// The diagonal is the cell's own, so the band narrows as the tree deepens --
    /// which is what makes the cost bounded rather than cubic.
    double mBandCells = 1.0;
    /// Octree: attach `refine:level`.
    bool mRecordLevels = true;
    /// Refuse to generate more cells than this, naming the option.
    std::int64_t mMaxCells = 20000000;
    /// How the distances themselves are computed.
    ///
    /// **Embedded by value on purpose, and the reason this whole header ships
    /// complete**: growing `SurfaceDistanceOptions` later would move everything
    /// after it here.
    SurfaceDistanceOptions mDistance;
};

/// The result of `compute_sdf`: the generated grid and what was found.
struct SdfResult {
    Mesh mMesh;                                       ///< The grid, carrying `sdf:distance`.
    SurfaceQuality mQuality;                          ///< The surface verdict.
    std::array<std::int64_t, 3> mDims{{0, 0, 0}};     ///< Root cell counts.
    std::array<double, 3> mOrigin{{0.0, 0.0, 0.0}};   ///< Lattice lo corner.
    std::array<double, 3> mSpacing{{0.0, 0.0, 0.0}};  ///< Finest cell size.
    std::int64_t mMaxDepth = 0;                       ///< 0 for a dense grid.
    std::int64_t mNumBanded = 0;                      ///< Queries clamped to the band.
};

/**
 * @brief Generate a grid over @p rSurface and fill it with signed distances.
 *
 * The one call that turns a surface into a field. `Voxel` builds the dense
 * lattice the options describe; `Octree` builds a coarse root lattice and
 * refines the cells within `mBandCells` diagonals of the surface, `mMaxDepth`
 * times, through `refine`'s `Balanced` closure.
 *
 * **The field is computed once, on the final mesh.** An octree pass cannot
 * carry it: `refine` interpolates `point_data`, so a field attached mid-way
 * would come out a smooth, plausible interpolation of the coarse values rather
 * than the distance. The generated grid also carries the `sdf:*` `field_data`
 * header describing itself.
 *
 * @param rSurface the surface to measure against.
 * @param rOptions the grid to generate and how to fill it.
 * @return the grid carrying `sdf:distance`, plus the surface verdict and the
 *         lattice geometry.
 * @throws std::invalid_argument when neither or both of `mResolution` and
 *         `mCellSize` are set for `Voxel`, when either is set for `Octree`, on a
 *         non-positive `mRootResolution`/`mBandCells`, a negative `mMaxDepth`,
 *         and when the grid exceeds `mMaxCells` (checked again after every
 *         octree pass, so a runaway band fails by name rather than by OOM).
 */
MESHIOPLUSPLUS_API SdfResult compute_sdf(const Mesh& rSurface, const SdfOptions& rOptions = {});

}  // namespace meshioplusplus
