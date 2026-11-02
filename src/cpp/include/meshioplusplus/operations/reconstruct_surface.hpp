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
 * @file operations/reconstruct_surface.hpp
 * @brief A triangle surface from an oriented point cloud (roadmap 7.5.2):
 * scans read from PLY, PCD or XYZ carry points and normals but no cells.
 *
 * `reconstruct_surface(cloud, options)` contours the signed distance to the
 * samples' tangent planes (H. Hoppe, T. DeRose, T. Duchamp, J. McDonald and
 * W. Stuetzle, "Surface reconstruction from unorganized points", SIGGRAPH
 * 1992):
 *
 * - **The samples** are the points with a finite, non-zero normal in
 *   `point_data[mNormalsName]` (normalised); with the default name and no
 *   such array, the `nx`, `ny` and `nz` scalars a PLY scan reads as.
 *   Without either the normals are estimated first (`estimate_normals`,
 *   `mNormalsK` neighbours, consistent orientation), when `mEstimateNormals`
 *   allows it.
 * - **The lattice** has cells of `mCellSize`; 0 takes `extent / mResolution`
 *   when a resolution is given, else the sample spacing (the median distance
 *   from a sample to its fourth-nearest other sample: about the
 *   nearest-neighbour distance on a regular sampling, and the gaps' size on a
 *   random one). Only the cells within
 *   `ceil(s / h) + 1` cells of a sample are evaluated (a band about the
 *   samples, `s` and `h` below), and only those the surface crosses are cut,
 *   so the cost grows with the surface's area, not the box's volume.
 * - **The field** at a lattice node `x` blends the signed distances
 *   `(x - p_i) . n_i` to the tangent planes of its `mK` nearest samples,
 *   weighted by `exp(-|x - p_i|^2 / s^2)` (`s` the larger of the spacing and
 *   the cell size); `mK = 1` is Hoppe's field exactly. A node is
 *   **undefined** when its projection onto the nearest sample's plane,
 *   `x - ((x - p_0) . n_0) n_0`, lies farther than `mMaxDistance` from that
 *   sample (Hoppe's test, 0 taking twice the larger of the spacing and the
 *   cell size): that is what keeps an open scan open instead of capping it.
 * - **The surface** is the field's zero level set by marching tetrahedra
 *   (`detail/marching.hpp`, the cutter `isosurface` uses) over the band's
 *   hexahedra split into tetrahedra; a cell with an undefined corner is not
 *   cut. A node value closer to zero than the larger of a thousandth of the
 *   cell and 1e-4 of the cloud's extent is moved to that distance (its sign
 *   kept), so no crossing collapses into a speck the cutter would drop.
 *   Quads are split along their 0-2 diagonal and coincident crossings
 *   welded, so the result is one `triangle` block,
 *   watertight wherever the field is defined, wound so its normals point the
 *   way the samples' normals do (outward, for outward normals).
 *
 * The result has no data arrays; it states the input's length unit when the
 * input does. The nearest-sample queries (`detail/point_knn.hpp`) are exact
 * and the cut serial, so the output is identical on every backend and thread
 * count.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// What `reconstruct_surface` contours, and on what lattice.
struct ReconstructSurfaceOptions {
    /// The point data holding the normals, Float `(n, 3)`.
    std::string mNormalsName = "normals";
    /// Tangent planes blended per lattice node; at least 1 (clamped to the
    /// sample count).
    std::int64_t mK = 8;
    /// The lattice cell size; 0 for automatic (see the file comment).
    double mCellSize = 0.0;
    /// Cells along the cloud's longest axis, when `mCellSize` is 0; 0 for
    /// the sample spacing.
    std::int64_t mResolution = 0;
    /// Hoppe's trim distance; 0 for automatic.
    double mMaxDistance = 0.0;
    /// Estimate the normals when the array is absent (else throw).
    bool mEstimateNormals = true;
    /// Neighbours per normal fit when they are estimated.
    std::int64_t mNormalsK = 16;
    /// The most lattice cells the band may hold; a larger band throws,
    /// naming a cell size that fits.
    std::int64_t mMaxCells = 20'000'000;
};

/// What `reconstruct_surface` built.
struct ReconstructSurfaceResult {
    /// The surface: points and one `triangle` block.
    Mesh mMesh;
    /// The samples used (points with a finite, non-zero normal).
    std::int64_t mNumSamples = 0;
    /// The median distance from a sample to its fourth-nearest other sample.
    double mSpacing = 0.0;
    /// The lattice cell size used.
    double mCellSize = 0.0;
    /// The trim distance used.
    double mMaxDistance = 0.0;
    /// Lattice nodes evaluated (the band's), and how many were undefined.
    std::int64_t mNumLatticePoints = 0;
    std::int64_t mNumUndefined = 0;
    /// Triangles, and their edge-connected pieces.
    std::int64_t mNumTriangles = 0;
    std::int64_t mNumComponents = 0;
    /// Whether the normals were estimated (the array was absent).
    bool mEstimatedNormals = false;
};

/**
 * @brief A triangle surface through an oriented point cloud (see the file
 *        comment).
 * @throws std::invalid_argument on points without three coordinates, a
 *         non-finite coordinate, fewer than three samples, samples that all
 *         coincide, a missing normals array with `mEstimateNormals` off, a
 *         malformed one, a negative size or distance, both `mCellSize` and
 *         `mResolution`, `mK` below 1, or a band over `mMaxCells`.
 */
MESHIOPLUSPLUS_API ReconstructSurfaceResult
reconstruct_surface(const Mesh& rMesh, const ReconstructSurfaceOptions& rOptions = {});

}  // namespace meshioplusplus
