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
 * @file operations/fit_curves.hpp
 * @brief Lines, arcs and circles from chained sections (roadmap 7.6.2): the
 * sketch a drawing, a DXF file or a CAD profile is made of.
 *
 * The input is a chained line mesh (`slice(chain=True)` or `sections`): its
 * `section:chain` cell data groups the segments, in walking order. Each
 * chain is planar; its plane (the frame) is `mNormal` through the chain's
 * points when given, else the least-spread direction of the points of all
 * the chains cut at one `section:offset` (all chains when there is none).
 * The frame's in-plane axes complete the normal (`lsq_frame`), its origin is
 * the normal times the plane's offset, so the sketches of parallel sections
 * share one 2-D frame.
 *
 * Per chain, in its frame:
 *
 * 1. **A whole circle.** A closed chain whose points all lie within
 *    `mTolerance` of one circle (Kasa's start, then the geometric fit) is
 *    one circle.
 * 2. **Split.** Otherwise the run is split recursively (a closed chain first
 *    at its first point and the point farthest from it): a run that a line
 *    (total least squares) explains within the tolerance is a line; failing
 *    that, one a circle explains, with at least `mMinArcPoints` points and a
 *    sweep of at least `mMinArcSweepDeg`, is an arc; failing both it is
 *    split at its point farthest from its chord.
 * 3. **Merge.** Neighbouring pieces whose union one line or one arc explains
 *    within the tolerance are merged, left to right, until none is (across
 *    the seam of a closed chain too).
 * 4. **Snap.** With `mSnapAngleDeg` > 0, a line within that angle of the
 *    frame's axes or their 45-degree diagonals is turned onto it about its
 *    midpoint.
 * 5. **Corners.** Consecutive pieces' shared end moves to their intersection
 *    (a line-line, line-arc or arc-arc one, the nearest to the data point)
 *    when one lies within ten tolerances of it.
 *
 * Each piece reports its type, its ends, centre, radius and signed sweep in
 * the frame (counter-clockwise positive), its largest residual over the
 * chain points it spans, and its chain and frame. The input comes back with
 * the Int64 cell data `sketch:piece` (each segment's piece) and the field
 * data `sketch:type`, `sketch:params` (ten numbers a piece: start x, y, end
 * x, y, centre x, y, radius, sweep in degrees, residual, points),
 * `sketch:chain`, `sketch:frame` (its frame index) and `sketch:frames` (a
 * row-major 4x4 per frame whose columns are the in-plane axes, the normal
 * and the origin), so a format or a later step reads the sketch. Residuals
 * are measured at the chain's points, not between them. Core only, no numpy
 * twin. Deterministic: no sampling, every tie to the lower index.
 * @see doc/fit_curves.md
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

/// Cell data (Int64): each segment's sketch piece.
inline constexpr const char* kSketchPieceName = "sketch:piece";
/// Field data: per piece its type (`SketchPieceType`), parameters, chain and frame.
inline constexpr const char* kSketchTypeName = "sketch:type";
inline constexpr const char* kSketchParamsName = "sketch:params";
inline constexpr const char* kSketchChainName = "sketch:chain";
inline constexpr const char* kSketchFrameName = "sketch:frame";
/// Field data: the frames, a row-major 4x4 each.
inline constexpr const char* kSketchFramesName = "sketch:frames";

/// What a sketch piece is.
enum class SketchPieceType : std::int32_t {
    Line = 1,
    Arc = 2,
    Circle = 3,
};

/// Options for `fit_curves`.
struct FitCurvesOptions {
    /// The largest residual a piece may have; 0 takes 1e-3 of the input's
    /// bounding-box diagonal.
    double mTolerance = 0.0;
    /// Snap lines within this many degrees of the frame's axes or their
    /// diagonals; 0 does not.
    double mSnapAngleDeg = 0.0;
    /// The sketch plane's normal; zero takes it from the points.
    std::array<double, 3> mNormal{{0.0, 0.0, 0.0}};
    /// An arc spans at least this many chain points and degrees.
    std::int64_t mMinArcPoints = 4;
    double mMinArcSweepDeg = 5.0;
};

/// One sketch piece, in its frame's 2-D coordinates.
struct SketchPiece {
    SketchPieceType mType = SketchPieceType::Line;
    std::array<double, 2> mStart{{0.0, 0.0}};
    std::array<double, 2> mEnd{{0.0, 0.0}};
    std::array<double, 2> mCentre{{0.0, 0.0}};  ///< arcs and circles
    double mRadius = 0.0;                       ///< arcs and circles
    double mSweepDeg = 0.0;  ///< arcs: signed, counter-clockwise positive; circles 360
    double mResidual = 0.0;  ///< the largest distance of its chain points
    std::int64_t mNumPoints = 0;
    std::int64_t mChain = 0;
    std::int64_t mFrame = 0;
};

/// The fitted sketch.
struct FitCurvesResult {
    /// The input with `sketch:piece` and the `sketch:*` field data.
    Mesh mMesh;
    std::vector<SketchPiece> mPieces;
    /// Per frame, the row-major 4x4 whose columns are e1, e2, the normal and
    /// the origin: `x = F (u, v, 0, 1)`.
    std::vector<std::array<double, 16>> mFrames;
    std::int64_t mNumLines = 0;
    std::int64_t mNumArcs = 0;
    std::int64_t mNumCircles = 0;
    double mMaxResidual = 0.0;
    double mTolerance = 0.0;  ///< the tolerance applied
    std::vector<std::string> mWarnings;
};

/**
 * @brief Lines, arcs and circles from the chains of @p rSections (see the
 * file comment).
 * @throws std::invalid_argument on a bad option, a mesh with cells other
 *         than lines, or one without `section:chain` cell data.
 */
MESHIOPLUSPLUS_API FitCurvesResult fit_curves(const Mesh& rSections,
                                              const FitCurvesOptions& rOptions = {});

/** @brief `"line"`, `"arc"` or `"circle"`. */
MESHIOPLUSPLUS_API const char* sketch_piece_type_name(SketchPieceType Type);

}  // namespace meshioplusplus
