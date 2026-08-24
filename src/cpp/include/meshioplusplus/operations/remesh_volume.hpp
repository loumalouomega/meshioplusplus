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
//  Attribution:     Implemented from the PUBLISHED DESCRIPTION of Labelle &
//                   Shewchuk, "Isosurface Stuffing: Fast Tetrahedral Meshes
//                   with Good Dihedral Angles" (SIGGRAPH 2007) only. Neither
//                   the paper's reference implementation (Stellar,
//                   non-commercial) nor TetGen (AGPL-3.0, also by Shewchuk)
//                   is read or vendored here. See doc/remesh_volume.md.
//
#pragma once

/**
 * @file remesh_volume.hpp
 * @brief Retetrahedralize a volume (or a closed surface) at a caller-chosen
 * resolution by isosurface stuffing: a body-centered cubic (BCC) lattice,
 * warped onto the surface, then cut.
 *
 * The volumetric sibling of `remesh` -- the same relationship `decimate_volume`
 * has to `decimate`, and the roadmap's own framing of this bullet. Nothing
 * else in this repo can *raise* the quality of a tetrahedral mesh at a chosen
 * target resolution: `refine` subdivides the input's own (possibly bad) cells,
 * `decimate_volume` can only remove elements, `smooth` (see
 * `SmoothMethod::Odt`) only moves points. `remesh_volume` instead discards the
 * input's own tetrahedra entirely and generates a fresh lattice-based mesh,
 * so output element quality is a property of the lattice construction, not of
 * the input.
 *
 * ### Why a lattice rather than Delaunay/CVD/ODT remeshing
 *
 * The roadmap bullet this closes names "optimal Delaunay triangulation / CVD
 * in 3D" -- the literal 3D counterpart of `remesh.hpp`'s own surface
 * clustering. That method needs robust orientation/in-sphere predicates,
 * which is exactly where this project's dependency-free posture stops paying
 * (the same trade `doc/roadmap.md`'s Delaunay-meshing bullet already
 * records): the permissive candidates are large new dependencies (Geogram),
 * and the smaller ones (TetGen, CGAL's `Mesh_3`) are AGPL/GPL, a poor fit for
 * an MIT project even as an optional dependency. Isosurface stuffing needs no
 * predicate at all -- only a signed-distance field, which this repo already
 * has (`sample_distance`/`distance_to_surface`, `operations/sdf.hpp`) -- and
 * gives every uncut lattice tet a dihedral angle from a small FIXED set
 * (below) rather than an unbounded, mesh-dependent one, so it closes the
 * bullet's *capability* (retetrahedralize at a chosen resolution with good
 * element quality) even though it is not the named method.
 * `SmoothMethod::Odt` (`operations/smooth.hpp`) closes the "ODT" half of the
 * bullet's name honestly and separately: ODT *smoothing* on an existing tet
 * mesh's fixed connectivity, not ODT *remeshing*.
 *
 * ### The algorithm
 *
 * 1. **Lattice.** A body-centered cubic lattice -- two interleaved simple
 *    cubic lattices, the second offset by half a cell on every axis -- is
 *    generated over the (padded) bounding box, sized exactly like
 *    `SdfOptions`' own `mResolution`/`mCellSize`/`mBounds`/`mPadding`/
 *    `mPaddingRelative` vocabulary (reused verbatim, not reinvented, since
 *    `remesh_volume` is a lattice generator exactly as `compute_sdf` is);
 *    `mResolution` is validated to resolve to genuinely CUBIC cells (a
 *    non-cube bounding box combined with a per-axis resolution can otherwise
 *    resolve to a rectangular lattice, which the construction below assumes
 *    away) -- named by error, never silently distorted.
 *
 *    Its natural tetrahedralization -- 12 congruent tets per cell, each
 *    cell's body-center connected across one fixed diagonal of each of the
 *    cell's 6 faces -- gives every UNCUT tet dihedral angles from the fixed
 *    set {45, 60, 90, 120} degrees, independent of cell size (every cell is a
 *    translate of the same unit construction). This is this project's OWN
 *    derivation, not a transcription of Labelle & Shewchuk's own
 *    tetrahedralization (which is not read): the paper reports a tighter
 *    {60, 90, 120} set for its own construction, and this one is a
 *    deliberately simpler, independently-verified alternative that still
 *    delivers a small, mesh-size-independent bound rather than an unbounded
 *    one -- honestly reported as what it is, not as a reproduction of the
 *    paper's exact angles. The 12-tet table is DERIVED and algebraically
 *    pinned in `remesh_volume.cpp` (tiling with no gap or overlap,
 *    orientation, volume-sum, and the dihedral set itself), never
 *    transcribed from anywhere.
 * 2. **Classify.** Every lattice vertex gets a signed distance to the surface
 *    (`sample_distance`'s own machinery) and an inside/outside LABEL from
 *    that sign, fixed for the rest of the algorithm.
 * 3. **Warp.** A lattice vertex within `mWarpFraction * h` of the surface
 *    (`h` the lattice's own cell size -- every incident edge's length is a
 *    fixed multiple of `h` regardless of vertex type, so one uniform
 *    threshold covers both the A and B sublattices) is moved (its LABEL is
 *    unchanged) onto the nearest point on the surface -- this is what
 *    removes the slivers a naive cut would otherwise produce near the
 *    boundary. Warping is load-bearing,
 *    not an optimization: disabling it (`mWarpFraction = 0`) still produces a
 *    watertight, correctly-oriented mesh, but with arbitrarily bad dihedral
 *    angles on boundary-adjacent tets -- plain marching tetrahedra's usual
 *    failure mode.
 * 4. **Cut.** A tet whose 4 vertices share one label is kept whole (using
 *    each vertex's current, possibly-warped position) or discarded. A
 *    straddling tet is cut by the same sign-mask case table marching
 *    tetrahedra uses (`detail/marching.hpp`'s sibling for VOLUME rather than
 *    surface output): each crossing edge gets one cut point, reusing an
 *    incident vertex's warped position when one of the edge's two endpoints
 *    was warped (which is the common case near a well-resolved surface,
 *    and is what carries the quality benefit of warping into the cut
 *    region), falling back to a fresh linear-interpolation point along the
 *    edge otherwise (`detail/marching.hpp`'s own convention) -- so a cut
 *    tet's boundary is always an honest approximation of the true crossing,
 *    never expanded by a warp that moved a vertex further than its own edge.
 *    This is a considered simplification of the paper's own two-threshold,
 *    priority-ordered warp/cut interaction (see doc/remesh_volume.md), kept
 *    honestly distinct from a claim of literal fidelity to it; the achieved
 *    dihedral-angle improvement over an unwarped cut is MEASURED, not merely
 *    asserted from the paper, exactly as `RemeshOptions::mMaxAnisotropy`'s own
 *    default was. Reusing a warped vertex's position for every cut edge it
 *    touches is also what can leave a SMALL number of the output's own
 *    boundary edges non-manifold on some inputs, at `mWarpFraction > 0` only
 *    -- reported as `RemeshVolumeResult::mNumNonManifoldEdges` rather than
 *    hidden, and always exactly 0 at `mWarpFraction = 0` (a plain cut has no
 *    such reuse to create the defect). See that field's own doc comment and
 *    doc/remesh_volume.md for the measured magnitude.
 *
 * Determinism follows the repo's standard phase split: signed distances and
 * warps are computed in parallel into disjoint per-vertex slots; cut points
 * are deduplicated by a SERIAL, ascending-edge-key pass (the
 * `refine.cpp`/`surface.cpp` idiom), so two lattice tets sharing a face agree
 * bit-for-bit on that face's cut points and the output is watertight by
 * construction, not by luck.
 *
 * ### What does NOT survive
 *
 * Identical to `remesh`: the output has entirely new points and new
 * connectivity, so there is no point or cell map, and `point_data`/
 * `cell_data`/regions are dropped (`warn_regions_dropped`, never silent);
 * `field_data` carries through. Compose with `interpolate`/
 * `conservative_interpolate` to transfer a field onto the result.
 *
 * ### Scope
 *
 * Accepts a volume mesh (its boundary is extracted internally via
 * `extract_surface`, which is closed by construction) OR a closed surface
 * directly (triangle/quad/rectangular-polygon, exactly `build_triangle_soup`'s
 * own scope) -- unlike plain `remesh`, which only ever works on a surface and
 * throws on a volume input. The BOUNDARY of the result is therefore only ever
 * as good as the lattice resolution and the warp step can make it; a caller
 * wanting a genuinely high-quality boundary composes with `extract_surface`
 * + `remesh` on the result rather than expecting this operation to build one
 * in.
 *
 * Standard C++ and the uniform mesh API only, so it compiles under every mesh
 * backend. This is an operation, not a file format -- it is deliberately not
 * in the format registry.
 */

// System includes
#include <cstdint>
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/sdf.hpp"

namespace meshioplusplus {

/// Options for `remesh_volume`. Mirrors `SdfOptions`' lattice-sizing
/// vocabulary field-for-field (this operation, like `compute_sdf`, is a
/// lattice generator) but never embeds `SdfOptions` itself: `SdfOptions`
/// carries octree-specific fields (`mRootResolution`, `mMaxDepth`,
/// `mBandCells`, `mRecordLevels`) this operation has no use for -- there is
/// no adaptive refinement here, only one fixed-resolution BCC lattice.
struct RemeshVolumeOptions {
    /// Cell counts per axis of the ROOT (pre-cut) BCC lattice; unset derives
    /// them from `mCellSize`. Exactly one of `mResolution`/`mCellSize` must be
    /// set, `SdfOptions`' own rule.
    std::optional<std::array<std::int64_t, 3>> mResolution;
    /// Cubic cell size of the root lattice; unset derives it from `mResolution`.
    std::optional<double> mCellSize;
    /// Explicit bounds `{lo[3], hi[3]}`; unset uses the surface's bounding box.
    std::optional<std::array<double, 6>> mBounds;
    /// Padding added to every side, in world units.
    double mPadding = 0.0;
    /// Padding added to every side, as a fraction of the bounding-box
    /// diagonal. A lattice that stops exactly at the surface leaves nowhere
    /// for a fully-outside cell to exist, hence the non-zero default -- the
    /// same reasoning `SdfOptions::mPaddingRelative` states.
    double mPaddingRelative = 0.1;
    /// Refuse to generate a root lattice with more cells than this, naming
    /// the option -- checked BEFORE any lattice is built, `SdfOptions::mMaxCells`'s
    /// own contract.
    std::int64_t mMaxCells = 20000000;
    /// Refuse an output with more tets than this, naming the option --
    /// checked AFTER cutting, since the true output count (unlike the root
    /// lattice's) depends on how much of the lattice the surface excludes and
    /// cannot be bounded in advance without doing the work.
    std::int64_t mMaxTets = 20000000;

    /// Fraction of a lattice vertex's own shortest incident edge length within
    /// which it may be warped onto the surface. `0` disables warping entirely
    /// (a plain, unwarped marching-tetrahedra cut of the BCC lattice --
    /// watertight and correctly oriented, but with arbitrarily bad boundary
    /// dihedral angles); see `doc/remesh_volume.md` for the measured
    /// sweep this default was chosen from.
    double mWarpFraction = 0.35;

    /// Sign, region restriction and watertight-check policy for every surface
    /// query this operation makes. Embedded by value on purpose, and the
    /// reason this whole header ships complete: growing
    /// `SurfaceDistanceOptions` later would move everything after it here
    /// (`SdfOptions::mDistance`'s own documented ABI caveat).
    SurfaceDistanceOptions mDistance;
};

/// The result of `remesh_volume`: the generated tet mesh and what was found.
struct RemeshVolumeResult {
    Mesh mMesh;                ///< The output mesh, a single `tetra` block.
    SurfaceQuality mQuality;   ///< The surface verdict (see `surface_watertight_check`).
    std::int64_t mNumTets = 0;         ///< Tets emitted.
    std::int64_t mNumVerticesWarped = 0;  ///< Lattice vertices moved onto the surface.
    std::int64_t mNumTetsRejected = 0;    ///< Cut sub-tets discarded as degenerate.
    /// Non-manifold edges of the OUTPUT tet mesh's own boundary (via
    /// `extract_surface` + `surface_watertight_check` on `mMesh` itself, not
    /// the input surface `mQuality` reports on). `mWarpFraction = 0` gives a
    /// mathematically watertight cut and this is always 0 there -- pinned by
    /// `AWarplessCutIsExactlyWatertight`. A nonzero `mWarpFraction` reuses a
    /// warped lattice vertex's position for every cut edge it touches
    /// (`rvol_cut_tet`'s doc comment), and on some inputs this can leave a
    /// small number of coincident-but-topologically-distinct boundary edges
    /// near the warped region -- MEASURED here rather than silently accepted:
    /// `doc/remesh_volume.md` records the swept magnitude (small and bounded
    /// relative to the output's own boundary edge count on every fixture
    /// tried) and a caller needing an exactly watertight result should check
    /// this counter, or fall back to `mWarpFraction = 0` / post-process with
    /// `clean(weld=true)`.
    std::int64_t mNumNonManifoldEdges = 0;
};

/**
 * @brief Retetrahedralize @p rMesh at a caller-chosen resolution by isosurface
 *        stuffing.
 * @param rMesh a volume mesh (its boundary is extracted via `extract_surface`)
 *        or a closed surface (`build_triangle_soup`'s own scope: triangle,
 *        quad, rectangular polygon).
 * @param rOptions the lattice to generate, the warp fraction, and the surface
 *        query policy.
 * @return the output tet mesh plus the surface verdict and counters.
 * @throws std::invalid_argument when neither or both of `mResolution` and
 *         `mCellSize` are set, on a non-positive resolution/cell size, on a
 *         negative `mWarpFraction`, when the root lattice would exceed
 *         `mMaxCells`, when the cut output would exceed `mMaxTets`, and on
 *         any input `build_triangle_soup`/`extract_surface` itself refuses
 *         (a higher-order or ragged surface block, an empty surface).
 */
MESHIOPLUSPLUS_API RemeshVolumeResult remesh_volume(const Mesh& rMesh,
                                                    const RemeshVolumeOptions& rOptions = {});

}  // namespace meshioplusplus
