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
 * @file operations/extrude.hpp
 * @brief Sweep a surface into a volume, or a curve into a surface, along a
 * vector (`extrude`, roadmap 3.1.2) or around an axis (`revolve`, 3.1.3).
 *
 * The input's top-dimensional cells are the **profile**: `triangle` and `quad`
 * cells sweep into a volume, or, when there are none, `line` cells into a
 * surface. Each profile cell gives one cell per layer -- a triangle a `wedge`,
 * a quad a `hexahedron`, a line a `quad` -- between consecutive levels of
 * copies of its points. Cells of lower dimension (a gmsh mesh's boundary
 * lines, the vertices at a curve's ends) are not swept: what they name is
 * carried as described below. Higher-order, polygon, polyhedron and volume
 * cells are refused (`convert_cells` linearizes a quadratic profile).
 *
 * **extrude** copies the points to `L` levels along `mDirection`, the whole
 * displacement: level `k` is `x + t_k d`, with `t_k = k / L`, a geometric
 * grading `t_k = (r^k - 1) / (r^L - 1)` (`mRatio = r`, the ratio of each layer
 * to the one below), or explicit increasing fractions ending at 1 (`mOffsets`).
 * A face is turned to point along the direction first (counted in
 * `mNumFlipped`), so every cell has a positive volume; a face or line parallel
 * to the direction is refused by its index. A curve's quads face `(b - a) x d`
 * for a line from `a` to `b`.
 *
 * **revolve** turns the points about the axis through `mAxisPoint` along
 * `mAxis` by `mAngle` degrees (the right-hand rule; negative turns the other
 * way) in `L` equal steps. A full turn (|angle| = 360) closes on itself: the
 * last level is the profile again. A point on the axis (within `mTolerance` of
 * the points' bounding-box diagonal) stays where it is, so the cells around it
 * collapse: a triangle with a corner on the axis gives pyramids, one with two
 * tetrahedra, a quad with an edge on the axis wedges, a quad with one corner
 * on it pyramids and tetrahedra (split along its first diagonal), and a line
 * with one end on it triangles. A face is turned to point along its motion;
 * one that crosses the axis or moves within its own plane is refused, and so
 * is a turn of more than 180 degrees of faces on both sides of the axis,
 * which would sweep through each other. A full turn needs at least 3 steps,
 * and no step may reach 180 degrees.
 *
 * **Output**: the swept cells only (wedge, hexahedron, pyramid, tetra blocks;
 * quad, triangle for a curve), layer by layer, so extruding `rectangle` gives
 * `box`'s points and cells exactly. Points are the profile's, then each
 * level's own copies. Named `Side` regions, as the parametric primitives name
 * theirs: `bottom` (the profile) and `top` (the last level), or `start` and
 * `end` for a partial turn, and `side` (what the profile's open border sweeps
 * to); a full turn has neither `start` nor `end`. An input region of the same
 * name wins, with a warning.
 *
 * **Carried**: point data per copy, cell data per layer (an array whose dtype
 * or width differs between the profile's blocks is dropped, with a warning),
 * field data. A `Cell` region of profile cells names their columns; a `Point`
 * region every copy of its points; a `Side` region of profile facets the faces
 * their facets sweep to; a `Cell` region of lower-dimensional cells the faces
 * (or edges) those cells sweep to, as a `Side` region. Dimensions grow by one.
 * With `mRecordLayers`, an Int32 `extrude:layer` (`revolve:layer`) cell array
 * holds each cell's layer, 1 to `L`.
 *
 * Deterministic: every index is a function of the input's order.
 *
 * @see doc/extrude.md
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

/// Cell data (Int32) naming each extruded cell's layer, with `mRecordLayers`.
inline constexpr const char* kExtrudeLayerName = "extrude:layer";
/// Cell data (Int32) naming each revolved cell's layer, with `mRecordLayers`.
inline constexpr const char* kRevolveLayerName = "revolve:layer";

/// Options for `extrude`.
struct ExtrudeOptions {
    /// The whole displacement of the last level; not zero.
    std::array<double, 3> mDirection{{0.0, 0.0, 1.0}};
    /// Layers, 1 or more; with `mOffsets`, 1 (the default) or their count.
    std::int32_t mNumLayers = 1;
    /// Each layer's thickness over the one below's; positive, 1 is uniform.
    double mRatio = 1.0;
    /// Explicit level fractions: increasing, positive, the last 1; empty for
    /// `mNumLayers` levels graded by `mRatio`.
    std::vector<double> mOffsets;
    /// Attach `extrude:layer`.
    bool mRecordLayers = false;
    /// Refuse a sweep of more cells (profile cells times layers); non-positive
    /// disables the budget.
    std::int64_t mMaxCells = 20000000;
};

/// Options for `revolve`.
struct RevolveOptions {
    /// A point of the axis.
    std::array<double, 3> mAxisPoint{{0.0, 0.0, 0.0}};
    /// The axis direction; not zero. The default turns an `(r, z)` profile in
    /// the xy plane about the y axis.
    std::array<double, 3> mAxis{{0.0, 1.0, 0.0}};
    /// Degrees, 0 < |angle| <= 360; the sign follows the right-hand rule.
    double mAngle = 360.0;
    /// Steps of the turn; at least 3 for a full turn.
    std::int32_t mNumLayers = 16;
    /// "On the axis": within this fraction of the bounding-box diagonal.
    double mTolerance = 1e-9;
    /// Attach `revolve:layer`.
    bool mRecordLayers = false;
    /// As `ExtrudeOptions::mMaxCells`.
    std::int64_t mMaxCells = 20000000;
};

/// The swept mesh and what went into it.
struct SweepResult {
    Mesh mMesh;
    std::int64_t mNumCells = 0;
    std::int64_t mNumHexahedra = 0;
    std::int64_t mNumWedges = 0;
    std::int64_t mNumPyramids = 0;
    std::int64_t mNumTetrahedra = 0;
    std::int64_t mNumQuads = 0;      ///< a swept curve's
    std::int64_t mNumTriangles = 0;  ///< a swept curve's, where an end is on the axis
    /// Profile cells turned to face the sweep before it.
    std::int64_t mNumFlipped = 0;
    /// Profile points on the axis (revolve).
    std::int64_t mNumOnAxis = 0;
    /// The least scaled Jacobian of a volume cell (1 for a right prism); 0 for
    /// a swept curve.
    double mMinScaledJacobian = 0.0;
    /// A full turn, closed on itself (revolve).
    bool mFull = false;
    std::vector<std::string> mWarnings;
};

/**
 * @brief Sweep @p rMesh's surface (or curve) along a vector.
 * @throws std::invalid_argument on a bad option, a face or line parallel to
 *         the direction, or an exceeded budget.
 * @throws Unsupported for higher-order, polygon, polyhedron or volume cells.
 */
MESHIOPLUSPLUS_API SweepResult extrude(const Mesh& rMesh, const ExtrudeOptions& rOptions = {});

/**
 * @brief Sweep @p rMesh's surface (or curve) around an axis.
 * @throws std::invalid_argument on a bad option, a face that crosses the axis
 *         or turns within its own plane, or an exceeded budget.
 * @throws Unsupported as `extrude`.
 */
MESHIOPLUSPLUS_API SweepResult revolve(const Mesh& rMesh, const RevolveOptions& rOptions = {});

}  // namespace meshioplusplus
