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
 * @file slice.hpp
 * @brief Planar cross-section of a mesh: the intersection of a mesh with a
 * plane, one topological dimension below the cut cells. This is an operation,
 * not a file format — deliberately not in the registry — built entirely
 * through the uniform mesh API so it compiles under every mesh backend.
 *
 * `slice(mesh, options)` returns a **new mesh** that is the cross-section:
 *
 * - a **3D volume** mesh cut by the plane yields a **surface** mesh
 *   (`triangle`/`quad` section faces);
 * - a **2D surface** mesh cut by the plane yields a **line** mesh (segments).
 *
 * Unlike `crop` (plane mode), which keeps *whole cells* on one side, `slice`
 * computes the actual intersection and lowers the dimension.
 *
 * Robust construction (marching tetrahedra): the input is first simplexified
 * via `convert_cells` (so every 3D cell is a tetrahedron and every 2D cell a
 * triangle), which makes each cell's cross-section a well-defined convex
 * primitive with no angular point-sorting and no non-convex-hex ambiguity.
 * Consequence: a hex/wedge section is the union of its simplices' sections
 * (correctly located, just triangulated).
 *
 * The signed distance of node `i` to the plane is `d_i = dot(x_i - origin,
 * normal)`. An edge `(i,j)` is crossed iff `d_i` and `d_j` fall on opposite
 * sides; the crossing point is `x_i + t*(x_j - x_i)` with `t = d_i/(d_i - d_j)`,
 * and every `point_data` array is interpolated at the same `t` (output promoted
 * to Float64). Crossing points are computed from the **sorted** endpoint pair
 * and deduped by that key, so an edge shared by two simplices yields a single
 * output node — the section is watertight — and the result is byte-identical
 * across the three mesh backends, thread counts, and the C++/numpy boundary
 * (`tests/python/test_slice.py::test_cpp_matches_python`).
 *
 * Degeneracy rule (uniform): a node exactly on the plane (`d_i == 0`) is
 * classified on the **positive** side (`d_i >= 0`), which makes the sign mask
 * total; a plane grazing a shared face is therefore emitted by exactly one
 * incident simplex — no double emission. Any section primitive whose crossings
 * collapse to a point/line (area/length below a bbox-relative tolerance) is
 * dropped. Section faces are wound so their Newell normal points toward the
 * `+normal` side — a consistent global orientation, pinned by a gtest.
 *
 * Each section cell inherits its parent cell's `cell_data` row (replicated,
 * like simplexify); with `mRecordParentIds` an Int64 `slice:parent_cell`
 * `cell_data` array records the originating **input** cell (global, block-major
 * over the input mesh). Section points are all new, so `point_sets`/`cell_sets`,
 * `mesh.info` and `gmsh_periodic` are **not** carried.
 *
 * ### Chained sections (roadmap 7.6.1)
 *
 * With `mChain`, a surface's line section comes back as polylines: its
 * segments are walked into chains (`detail/polyline_chain.hpp`: open chains
 * from their lowest-id end, then loops from their lowest point, each walk
 * taking the lowest-id neighbour), emitted in chain order and walking
 * direction, the points renumbered in order of first visit, with the Int64
 * cell data `section:chain` (the chain's index) and `section:closed` (1 on a
 * loop). Points at exactly one position are one point of a chain (a plane
 * through a vertex cuts each of its edges there; the lowest id stands in), so
 * a segment between two of them is left out. Chains shorter than
 * `mMinLength` are dropped. A volume's section is
 * a surface and is refused by name. `sections` slices along a direction at
 * several offsets and merges the chained sections, with the Float64 cell
 * data `section:offset` and the chains numbered across them.
 */

// System includes
#include <array>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Options for `slice`.
struct SliceOptions {
    /// A point on the cutting plane.
    std::array<double, 3> mOrigin{{0.0, 0.0, 0.0}};

    /// The plane normal (need not be unit length; must be non-zero).
    std::array<double, 3> mNormal{{0.0, 0.0, 1.0}};

    /// Attach an Int64 `slice:parent_cell` `cell_data` array recording, per
    /// section cell, the global (block-major) index of the input cell it was
    /// cut from.
    bool mRecordParentIds = false;

    /// Walk a line section into chains (see the file comment).
    bool mChain = false;

    /// With `mChain`, drop chains shorter than this.
    double mMinLength = 0.0;
};

/// Cell data (Int64) of a chained section: each segment's chain.
inline constexpr const char* kSectionChainName = "section:chain";
/// Cell data (Int64) of a chained section: 1 when the segment's chain is a loop.
inline constexpr const char* kSectionClosedName = "section:closed";
/// Cell data (Float64) of `sections`: the offset each segment was cut at.
inline constexpr const char* kSectionOffsetName = "section:offset";

/// Options for `sections`.
struct SectionsOptions {
    /// The planes' normal; any non-zero length. Plane `k` is `n . x = mOffsets[k]`
    /// with `n` the unit direction.
    std::array<double, 3> mDirection{{0.0, 0.0, 1.0}};
    /// The planes' signed distances from the origin along the unit direction.
    std::vector<double> mOffsets;
    /// Drop chains shorter than this.
    double mMinLength = 0.0;
    /// As `SliceOptions::mRecordParentIds`.
    bool mRecordParentIds = false;
};

/**
 * @brief Computes the planar cross-section of @p rMesh.
 * @param rMesh the mesh to cut (never modified).
 * @param rOptions the plane (origin/normal) and the provenance flag.
 * @return a new mesh one topological dimension below the cut cells: a
 *         `triangle`/`quad` surface for a volume mesh, a `line` mesh for a 2D
 *         surface mesh, or an empty mesh when the plane misses the geometry (or
 *         the mesh has no 2D/3D cells).
 * @throws std::invalid_argument on a zero-length normal, or when a cell block
 *         cannot be simplexified (a polyhedron block).
 */
MESHIOPLUSPLUS_API Mesh slice(const Mesh& rMesh, const SliceOptions& rOptions = {});

/**
 * @brief The chained sections of the surface @p rMesh at several parallel
 * planes, merged into one line mesh (see the file comment).
 * @throws std::invalid_argument on a zero direction, no offsets, a
 *         non-finite offset, or a mesh whose sections are not lines.
 */
MESHIOPLUSPLUS_API Mesh sections(const Mesh& rMesh, const SectionsOptions& rOptions);

}  // namespace meshioplusplus
