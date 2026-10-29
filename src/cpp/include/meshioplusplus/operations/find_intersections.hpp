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
 * @file operations/find_intersections.hpp
 * @brief Where a surface intersects itself, or where two surfaces
 * interpenetrate: the cell pairs, the intersection segments as a line mesh
 * and a count per cell (roadmap 3.1.2).
 *
 * ### What is tested
 *
 * The surfaces are reduced to triangles the way `build_triangle_soup` does
 * (quads and polygons fanned from their first corner), candidate pairs are
 * the triangles whose bounding boxes meet over a bucket grid, and each pair
 * is decided by an exact triangle–triangle test built on Shewchuk's
 * `orient3d` and `orient2d`. Every reported hit is therefore exact for the
 * input's double coordinates; there is no tolerance and no "possible" class.
 *
 *  - **One surface.** Two triangles of one cell (a quad's fan) are never
 *    paired. Two triangles that share corners -- bitwise-equal coordinates,
 *    so an unwelded seam is not mistaken for a crossing -- are an
 *    intersection only when they meet beyond what they share: on a shared
 *    edge, when they fold flat onto each other; on a shared corner, when an
 *    edge of one crosses the other or runs into it from the corner. Any other
 *    contact, touching included, counts.
 *  - **Two surfaces.** Every contact counts, touching included: two parts
 *    resting on each other are reported, so a contact interface that must not
 *    count is better checked for clearance with `thickness`.
 *
 * ### What comes back
 *
 * The cell pairs as global (block-major) cell ids, sorted; the segments
 * where triangles cross, one `line` cell each with two points of its own and
 * the cells it came from as cell data; and the inputs, bit-identical, with
 * the number of cells each cell intersects as cell data. Two triangles in one
 * plane that overlap meet over an area, not a segment: they are counted in
 * `mNumCoplanar` and flagged, but give no segment.
 *
 * Deterministic: candidates are collected per triangle and sorted, the exact
 * test does not depend on order, and the results are reduced in triangle
 * order, so the output is identical across parallel backends and thread
 * counts.
 *
 * @see doc/find_intersections.md
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

/// Cell data (Int64) on the returned meshes: how many cells each cell
/// intersects, 0 for a clean cell.
inline constexpr const char* kIntersectionsCountName = "intersections:count";
/// Cell data (Int64) on the segment mesh: the cell of the first surface.
inline constexpr const char* kIntersectionsCellAName = "intersections:cell_a";
/// Cell data (Int64) on the segment mesh: the cell of the second surface (of
/// the first, for a self-intersection).
inline constexpr const char* kIntersectionsCellBName = "intersections:cell_b";

/// Options for `find_intersections`.
struct FindIntersectionsOptions {
    /// Restrict the (first) surface to this named `Cell` region; empty takes
    /// every 2-D cell.
    std::string mRegion;
    /// Restrict the second surface to this named `Cell` region.
    std::string mRegionB;
    /// Attach `intersections:count` to the returned meshes.
    bool mAttachFlags = true;
    /// Build the segment mesh.
    bool mSegments = true;
    /// Keep at most this many pairs in `mPairs` (0 keeps all); the counts,
    /// flags and segments still cover every pair.
    std::int64_t mMaxPairs = 0;
};

/// The pairs, the segments and the counts.
struct FindIntersectionsResult {
    /// The (first) input, bit-identical, with `intersections:count`.
    Mesh mMesh;
    /// The second input, likewise; empty for a self-intersection check.
    Mesh mMeshB;
    /// One `line` cell per crossing segment, with `intersections:cell_a` and
    /// `intersections:cell_b`; two points per segment, not welded.
    Mesh mSegments;
    /// The intersecting cell pairs, `(cell of A, cell of B)` (for one surface,
    /// the lower cell first), sorted.
    std::vector<std::array<std::int64_t, 2>> mPairs;
    std::int64_t mNumPairs = 0;          ///< cell pairs, every one (not capped)
    std::int64_t mNumTrianglePairs = 0;  ///< triangle pairs that meet
    std::int64_t mNumCoplanar = 0;       ///< of those, overlapping in one plane
    std::int64_t mNumSegments = 0;       ///< cells of `mSegments`
    std::int64_t mNumCellsA = 0;         ///< cells of A that intersect anything
    std::int64_t mNumCellsB = 0;         ///< cells of B that do (0 for one surface)
    /// `mPairs` was capped by `mMaxPairs`.
    bool mTruncated = false;
    /// No pair at all.
    bool mClean = true;
    std::vector<std::string> mWarnings;
};

/**
 * @brief Find where the surface @p rMesh intersects itself.
 * @param rMesh a surface mesh (triangles, quads, polygons; lines and vertices
 *        are ignored).
 * @param rOptions the region, what to build and the pair cap.
 * @throws std::invalid_argument on a volume, polyhedron or higher-order block
 *         (run `extract_surface` or `linearize` first), an unknown region or a
 *         negative cap.
 */
MESHIOPLUSPLUS_API FindIntersectionsResult
find_intersections(const Mesh& rMesh, const FindIntersectionsOptions& rOptions = {});

/**
 * @brief Find where the surfaces @p rA and @p rB interpenetrate or touch.
 * @throws std::invalid_argument as the one-surface form, for either input.
 */
MESHIOPLUSPLUS_API FindIntersectionsResult
find_intersections(const Mesh& rA, const Mesh& rB, const FindIntersectionsOptions& rOptions = {});

}  // namespace meshioplusplus
