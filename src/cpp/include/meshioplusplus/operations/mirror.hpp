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
 * @file operations/mirror.hpp
 * @brief Reflect a mesh across a plane, keeping every cell right side out
 * (roadmap 3.2.5).
 *
 * A reflecting `transform` (a matrix of negative determinant) moves the points
 * and leaves the cells as they were, so every one of them comes out inverted:
 * a tetrahedron of negative volume, a triangle whose normal points inward.
 * `mirror` reflects the points across the plane through `mOrigin` with normal
 * `mNormal` and permutes the nodes of every cell so its orientation is valid
 * again (`detail/cell_flip`): mid-edge, mid-face and centre nodes stay in their
 * slots, a polygon's ring and a polyhedron's face rings are reversed. `Side`
 * region entries follow their facets to their new local numbers; cell and
 * point regions keep their indices.
 *
 * With `mKeepOriginal`, the result is the input followed by its mirror image:
 * each block holds its own cells, then their reflected copies, and the points
 * are the input's, then the copies'. With `mWeld` as well (the default), a
 * point on the plane -- within `mTolerance` of the bounding-box diagonal -- is
 * its own image, so the two halves share it: a half model becomes a whole one,
 * conforming across the plane. A cell lying entirely on the plane would be its
 * own image too, so its copy is dropped; a `Side` facet lying on the plane is
 * interior once both halves meet, so its entries (both copies) are dropped and
 * counted, and its region is kept even when that leaves it empty.
 *
 * Data is carried row by row (a copy's row is its source's). With
 * `mReflectVectorData`, float arrays of three columns (vectors) are reflected
 * `v' = H v` and of nine (row-major 3x3 tensors) `A' = H A H`, with
 * `H = I - 2 n n^T`, on the reflected copies only -- as `transform`'s
 * `rotate_vector_data` treats them; a pseudo-vector (a vorticity, a normal
 * stored as data) would need the opposite sign and is not singled out. Field
 * data and property sets are copied.
 *
 * Two-dimensional points stay two-dimensional: the plane must then be a line
 * of the xy plane (`mNormal[2] == 0`), likewise for one-dimensional points.
 * Cell types without a flip (`triangle10` and beyond, the other Lagrange
 * orders, the VTK Lagrange and Bezier families) throw `Unsupported` by name.
 *
 * Deterministic: every output index is a function of the input's order.
 *
 * @see doc/mirror.md
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

/// Options for `mirror`.
struct MirrorOptions {
    /// A point of the mirror plane.
    std::array<double, 3> mOrigin{{0.0, 0.0, 0.0}};
    /// The plane's normal; any non-zero length.
    std::array<double, 3> mNormal{{1.0, 0.0, 0.0}};
    /// Return the input followed by its mirror image instead of the image alone.
    bool mKeepOriginal = false;
    /// With `mKeepOriginal`, share the points on the plane between both halves.
    bool mWeld = true;
    /// "On the plane": a distance at or below this fraction of the points'
    /// bounding-box diagonal; finite and non-negative.
    double mTolerance = 1e-9;
    /// Reflect float vector (3-column) and tensor (9-column) point and cell data.
    bool mReflectVectorData = false;
};

/// The mirrored mesh and what was done to it.
struct MirrorResult {
    Mesh mMesh;
    std::int64_t mNumFlippedCells = 0;   ///< reflected cells written, after drops
    std::int64_t mNumPointsOnPlane = 0;  ///< input points within the tolerance of the plane
    std::int64_t mNumDroppedCells = 0;   ///< copies dropped as lying on the plane (weld only)
    std::int64_t mNumDroppedSides = 0;   ///< Side entries dropped as interior (weld only)
    /// One entry per region a drop emptied or shrank, and per region not carried.
    std::vector<std::string> mWarnings;
};

/**
 * @brief Reflect @p rMesh across a plane, cells right side out.
 * @param rMesh any mesh.
 * @param rOptions the plane and what to keep.
 * @throws std::invalid_argument on a zero or non-finite normal, a non-finite
 *         origin, a bad tolerance, or a plane that leaves the points' dimension.
 * @throws Unsupported for a cell type with no orientation flip.
 */
MESHIOPLUSPLUS_API MirrorResult mirror(const Mesh& rMesh, const MirrorOptions& rOptions = {});

}  // namespace meshioplusplus
