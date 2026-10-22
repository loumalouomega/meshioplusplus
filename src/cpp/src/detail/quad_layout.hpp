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
 * @file detail/quad_layout.hpp
 * @brief A quad layout of a freeform region: a partition of its triangles
 * into four-sided discs, the curves they share and the extraordinary
 * vertices where other than four meet -- the skeleton `patch_network` fits
 * B-spline patches on (`fit_cad`'s network mode, v17.12.0).
 *
 *  1. **Features.** The region starts as the components its feature edges
 *     (dihedral angle above `mFeatureAngleDeg`, when `mSnapToFeatures`) and
 *     its initial labels separate: a crease is always a layout curve.
 *  2. **Discs.** A piece that is not a disc (or is pinched at a vertex) is
 *     cut by `patch_split_to_discs` (genus 0) or halved across its longest
 *     extent (`patch_bisect`, genus 1 and up) until every piece is one; a
 *     piece under `mMinTriangles` then joins the neighbour it shares most
 *     edges with, where the union stays a disc.
 *  3. **Singularities.** The region's smoothest cross field (`cross_field`,
 *     aligned with its boundary and features) places extraordinary
 *     vertices: a piece holding exactly one singularity of index `i` well
 *     inside it is cut along the `4 - i` separatrices traced from it (greedy
 *     walks along the mesh edges best aligned with the field, from
 *     directions spread evenly round the vertex and turned to the field) to
 *     its boundary, into `4 - i` pieces meeting there. A piece holding
 *     several stays whole: a B-spline patch follows a smooth surface across
 *     a singularity of its field.
 *  4. **Corners.** Each piece takes four corners on its boundary, in this
 *     order: the ones its construction fixed (a singularity, a separatrix's
 *     end), its sharp turns (the chord angle over a few boundary vertices,
 *     under 150 degrees, locally sharpest), the corners its neighbours
 *     already chose there (so a shared curve breaks at the same vertices on
 *     both sides), and the vertices that best split its longest gaps. A piece
 *     with five or more sharp turns (at most `4 + mMaxValenceDeviation`) is
 *     split into as many quads by paths from its innermost vertex to the
 *     middles of its sides instead.
 *  5. **Size.** With `mTargetSize`, every piece is cut into strips along
 *     each of its two directions, conformingly: a *sheet* is the set of
 *     directions (a piece's sides 0 and 2, or 1 and 3) linked across shared
 *     sides, each sheet takes the strip count its longest side needs (at
 *     most 16), and each piece cuts its sides at the same arc-length
 *     fractions its neighbours cut them at, joined across it by shortest
 *     paths. A piece that does not cut cleanly stays whole, and its
 *     neighbours' cuts end on its side as T-junctions.
 *  6. **Curves.** Every piece's sides are broken at every corner and junction
 *     into arcs, each shared by two pieces or on the region's boundary; the
 *     arcs a side runs over form a chain, the curve fitted once and shared
 *     (a corner of one piece in the middle of a neighbour's side is a
 *     T-junction: both sides are stretches of the same chain). Where the
 *     arcs so linked do not make a simple path, each is a chain of its own
 *     and a side runs over several, joined end to end.
 *
 * The turns are measured in the tangent plane (the chords projected off the
 * piece's normal there), so a boundary that curves only out of the surface,
 * along its normal curvature, has no corner.
 *
 * The layout is checked (every piece a disc with four distinct corners,
 * every interior edge shared by two pieces in opposite directions, every
 * chain a simple path, every side the stretches of chains between its
 * corners) and comes back with `mOk` false otherwise: the caller
 * falls back, never uses half a layout. Serial and deterministic: fixed loop
 * orders, ties by index.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/// What `quad_layout` aims for.
struct QuadLayoutOptions {
    /// The pieces' target extent (absolute); 0: no size split.
    double mTargetSize = 0.0;
    /// The valences allowed at an extraordinary vertex: `4 +- this`.
    int mMaxValenceDeviation = 2;
    /// Cut along feature edges (dihedral angle above `mFeatureAngleDeg`).
    bool mSnapToFeatures = true;
    double mFeatureAngleDeg = 30.0;
    /// Per triangle an initial label whose changes are cuts too; empty: none.
    std::vector<std::int64_t> mLabels;
    /// Place extraordinary vertices at the cross field's singularities.
    bool mUseCrossField = true;
    /// A piece with fewer triangles is not split further.
    std::size_t mMinTriangles = 8;
    /// Halvings before the layout gives up on a piece.
    int mMaxDepth = 12;
};

/// One four-sided piece.
struct QuadPatch {
    /// Its triangles (indices into the region's triangles), ascending.
    std::vector<std::size_t> mTris;
    /// Its boundary, a vertex cycle in the triangles' orientation, starting
    /// at its first corner.
    std::vector<std::size_t> mBoundary;
    /// The corners' positions in `mBoundary`, ascending, the first 0.
    std::array<std::size_t, 4> mCornerPos{};
};

/// A curve the layout fits once: a simple vertex path.
struct QuadChain {
    std::vector<std::size_t> mVerts;
    /// Whether every arc of it lies on the region's boundary.
    bool mOnBoundary = false;
    /// Whether it is a cut of the layout's own (not the region's boundary,
    /// not a crease): its zigzag of mesh edges stands for a smooth curve.
    bool mSmooth = false;
};

/// A stretch of a chain.
struct QuadSegment {
    std::size_t mChain = 0;
    /// Positions in the chain's vertices of the stretch's start and end, in
    /// the piece's boundary direction (`mBegin > mEnd` when it runs backwards).
    std::size_t mBegin = 0, mEnd = 0;
};

/// A piece's side: one stretch of a chain, or, where the T-junctions on it
/// leave no single curve both its pieces can share, stretches of several,
/// joined end to end.
struct QuadSide {
    std::vector<QuadSegment> mSegments;
};

/// Why a layout failed.
enum class QuadLayoutFailure : std::int32_t {
    None = 0,
    NonManifold,
    CutFailed,      ///< a piece no cut made a disc
    CornersFailed,  ///< a piece without four distinct corners
    Inconsistent,   ///< an interior edge not shared by two pieces in opposite directions
    CyclicChain,    ///< a chain that is not a simple path
};

/// A quad layout.
struct QuadLayout {
    bool mOk = false;
    QuadLayoutFailure mFailure = QuadLayoutFailure::None;
    /// Per region triangle, its piece.
    std::vector<std::int64_t> mPatchOfTri;
    std::vector<QuadPatch> mPatches;
    /// Per piece, its sides from corner 0 (side k from corner k to k + 1).
    std::vector<std::array<QuadSide, 4>> mSides;
    std::vector<QuadChain> mChains;
    /// Interior vertices where other than four pieces meet at a corner (and
    /// no side passes through), ascending, and how many meet there.
    std::vector<std::size_t> mExtraordinary;
    std::vector<int> mValence;
    /// Arcs shared by two pieces, and corners in the middle of a side.
    std::size_t mNumInteriorArcs = 0;
    std::size_t mNumTJunctions = 0;
    /// Singularities the separatrices split pieces at.
    std::size_t mNumSeparatrixSplits = 0;
};

/// The quad layout of a region (see the file comment). `rTris` must be
/// consistently oriented; `rPts` may hold points no triangle uses.
QuadLayout quad_layout(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       const QuadLayoutOptions& rOptions);

/// Check a layout's invariants (see the file comment) against its region:
/// what `quad_layout` checks before it sets `mOk`.
QuadLayoutFailure quad_layout_check(const QuadLayout& rLayout, const std::vector<PatchTri>& rTris);

}  // namespace detail
}  // namespace meshioplusplus
