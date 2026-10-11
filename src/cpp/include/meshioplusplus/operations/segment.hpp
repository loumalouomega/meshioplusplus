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
 * @file operations/segment.hpp
 * @brief Split a surface into smooth patches by region growing, bounded by its
 * feature edges -- the first step of reverse engineering a triangle skin into
 * a CAD model (`fit_primitives`, then the `step` writer).
 *
 * The surface is the mesh's 2-D cells (triangles, quads, polygons and their
 * quadratic variants, which contribute their corners). Two cells are
 * neighbours when they share an edge that exactly those two cells use. A patch
 * grows from its seed -- the lowest-numbered cell not yet in a patch -- across
 * every neighbour edge that is
 *
 *  - **smooth**: the dihedral angle between the two cells' unit normals, after
 *    reorienting an inconsistently wound pair, is at most `mFeatureAngleDeg`
 *    (the crease test `feature_edges` reports); and, when enabled,
 *  - **within the normal tolerance**: the neighbour's (orientation-corrected)
 *    normal lies within `mNormalToleranceDeg` of the seed's; and
 *  - **within the curvature tolerance**: the neighbour's curvature estimate
 *    differs by at most `mCurvatureTolerance` from the mean estimate of the
 *    cells the patch holds so far (taken in walk order, so it is
 *    deterministic, and robust to one noisy seed).
 *
 * Open, non-manifold and sharp edges are never crossed. A patch smaller than
 * `mMinPatchCells` is then merged, in one ascending pass, into the neighbouring
 * patch it shares the most edges with (ties to the lowest patch id); a patch
 * with no neighbour stays.
 *
 * The curvature estimate of a cell is the largest, over its smooth paired
 * edges, of the edge's dihedral angle in radians divided by the distance
 * between the two cells' centroids: an estimate of the larger principal
 * curvature -- `1/r` on a cylinder or sphere of radius `r`, `0` on a plane.
 * It is in inverse length units, so `mCurvatureTolerance` is too.
 *
 * ### Output
 *
 * The input, unchanged -- points, connectivity, block order, data and regions
 * are bit-identical -- plus:
 *
 *  - cell data `segment:patch` (Int64): the patch of each cell, numbered from 0
 *    in the order of each patch's lowest cell; `-1` on cells that are not
 *    surface cells (lines, vertices) or lie outside `mRegion`;
 *  - one `Cell` region per patch, named `<prefix>:<k>`. Cell regions already
 *    named `<prefix>:<integer>` (an earlier run's) are replaced;
 *  - field data `segment:adjacency` (Int64, `(m, 3)`): one row `(a, b, edges)`
 *    per pair of patches `a < b` sharing at least one paired edge, ascending.
 *
 * ### Determinism
 *
 * Normals and the sorted edge table are built in parallel; the growth itself is
 * a serial breadth-first walk (FIFO queue, neighbours in ascending cell order)
 * because which patch claims a cell depends on the order the walk reaches it.
 * The output is byte-identical on every backend and thread count, and to the
 * numpy twin `meshioplusplus.segment`.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Cell data (Int64): the patch of each cell, or -1.
inline constexpr const char* kSegmentPatchName = "segment:patch";
/// Field data (Int64, `(m, 3)`): `(a, b, shared edges)` per adjacent patch pair.
inline constexpr const char* kSegmentAdjacencyName = "segment:adjacency";

/// How `segment` grows its patches.
struct SegmentOptions {
    /// The largest dihedral angle, in degrees, still crossed; in `[0, 180]`.
    double mFeatureAngleDeg = 30.0;
    /// The largest angle, in degrees, between a cell's normal and its patch
    /// seed's; `0` disables the test. In `[0, 180]`.
    double mNormalToleranceDeg = 0.0;
    /// The largest difference between a cell's curvature estimate and its
    /// patch's mean so far, in inverse length units; `0` disables the test.
    double mCurvatureTolerance = 0.0;
    /// Patches with fewer cells are merged into a neighbour; at least 1.
    std::int64_t mMinPatchCells = 1;
    /// Restrict to the cells of this named `Cell` region; empty takes every cell.
    std::string mRegion;
    /// The name prefix of the patch regions: `<prefix>:<k>`. Must not be empty.
    std::string mPrefix = "patch";
};

/// The segmented mesh and what the segmentation found.
struct SegmentResult {
    /// The input with `segment:patch`, `segment:adjacency` and the patch regions.
    Mesh mMesh;
    /// Patches after merging.
    std::int64_t mNumPatches = 0;
    /// Patches merged into a neighbour for being smaller than `mMinPatchCells`.
    std::int64_t mNumMerged = 0;
    /// Paired edges sharper than the feature angle.
    std::int64_t mNumFeatureEdges = 0;
    /// Edges used by one cell.
    std::int64_t mNumBoundaryEdges = 0;
    /// Edges used by three or more cells.
    std::int64_t mNumNonManifoldEdges = 0;
};

/**
 * @brief Segment a surface into smooth patches.
 *
 * @param rMesh a surface mesh; lines and vertices are carried and labelled -1.
 * @param rOptions see `SegmentOptions`.
 * @throws std::invalid_argument on an option out of range, an unknown region
 *         name, or a volume cell block (take its skin first with `extract_skin`).
 */
MESHIOPLUSPLUS_API SegmentResult segment(const Mesh& rMesh, const SegmentOptions& rOptions = {});

}  // namespace meshioplusplus
