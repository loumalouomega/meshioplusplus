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
 * @file detail/primitive_grow.hpp
 * @brief Regions of a surface grown from a seed facet while one primitive
 * (plane, sphere, cylinder, cone, torus) explains them: the growth
 * `segment(method="primitive")` runs from its own seeds and `fit_primitives`
 * runs from a seed point a user names (roadmap 7.1.1), one implementation for
 * both.
 *
 * `build_primitive_grow_surface` collects the 2-D cells as face rings with
 * their unit normals and centroids, the edges two faces share (smooth below
 * the feature angle) and each face's neighbours across them;
 * `prepare_primitive_grow` adds the areas, the bounding-box diagonal and the
 * tolerances; `primitive_smooth_areas` numbers the smooth areas.
 * `primitive_grow_from_seed` then grows a plane and a curved surface through a
 * seed, trims each to what its own refit explains exactly and returns the
 * larger.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI. Every step is
 * serial and ordered (faces by index, neighbours ascending), so a region is
 * the same on every backend and thread count.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
// Project includes (private, not installed)
#include "primitive_geom.hpp"
#include "surface_rings.hpp"

namespace meshioplusplus {
namespace detail {

/// The largest RMS of a grown region, as a share of the tolerance.
inline constexpr double kPrimitiveGrowExactShare = 0.02;

/// An edge shared by exactly two distinct faces.
struct PrimitiveGrowPair {
    std::int64_t mA;  ///< the lower face
    std::int64_t mB;  ///< the higher face
    bool mInconsistent;
    bool mSmooth;
    double mAngleRad;  ///< NaN when either face is degenerate
};

/// The growth's view of the surface: per face its ring, area, normal and
/// centroid, its smooth neighbours, and the smooth area it lies in.
struct PrimitiveGrowInput {
    const std::vector<std::int64_t>* mStart;
    const std::vector<std::int64_t>* mNodes;
    const std::vector<double>* mXyz;
    const std::vector<double>* mNormals;
    const std::vector<double>* mCentroid;
    const std::vector<std::int64_t>* mAdjStart;
    const std::vector<std::pair<std::int64_t, std::int64_t>>* mAdj;
    const std::vector<PrimitiveGrowPair>* mPairs;
    std::vector<double> mArea;
    std::vector<std::int64_t> mSmoothArea;
    double mTol = 0.0;
    double mGrowTol = 0.0;  ///< how far a corner may lie from a growing region's primitive
    double mDiag = 0.0;
};

/// A surface prepared for growth, owning what `PrimitiveGrowInput` points at.
struct PrimitiveGrowSurface {
    SurfaceRings mRings;
    std::vector<double> mXyz;       ///< every point, padded to 3-D
    std::vector<double> mNormals;   ///< unit, three per face (zero when degenerate)
    std::vector<double> mCentroid;  ///< three per face
    std::vector<PrimitiveGrowPair> mPairs;
    std::vector<std::int64_t> mAdjStart;  ///< per face, its first entry in `mAdj`
    std::vector<std::pair<std::int64_t, std::int64_t>> mAdj;  ///< (neighbour, pair), ascending
    std::int64_t mNumBoundaryEdges = 0;
    std::int64_t mNumNonManifoldEdges = 0;
    std::int64_t mNumFeatureEdges = 0;

    std::size_t NumFaces() const { return mRings.NumFaces(); }

    /// The growth's input over this surface (areas, tolerances unset).
    PrimitiveGrowInput Input() const;
};

/**
 * @brief The selected 2-D cells of @p rMesh as a growth surface, the edges
 * between faces whose normals turn by more than @p FeatureAngleDeg (an
 * inconsistently wound pair is turned first) marked not smooth.
 * @throws std::invalid_argument (prefixed) on a volume or polyhedron block.
 */
PrimitiveGrowSurface build_primitive_grow_surface(const Mesh& rMesh,
                                                  const std::vector<std::uint8_t>& rSelected,
                                                  double FeatureAngleDeg, const char* pPrefix);

/// Fill `mArea` and `mDiag`, then `mTol = RelTolerance * mDiag` and
/// `mGrowTol = kPrimitiveGrowExactShare * mTol`.
void prepare_primitive_grow(PrimitiveGrowInput& rIn, std::size_t NumFaces, double RelTolerance);

/// Number the smooth areas (faces joined across smooth edges) in the order of
/// their lowest face into `mSmoothArea`; returns each area's faces, ascending.
std::vector<std::vector<std::int64_t>> primitive_smooth_areas(PrimitiveGrowInput& rIn,
                                                              std::size_t NumFaces);

/// The samples a fit takes from faces: their distinct corners, normals and centroids.
void primitive_grow_samples(const PrimitiveGrowInput& rIn, const std::vector<std::int64_t>& rFaces,
                            std::vector<double>& rPts, std::vector<double>& rNormals,
                            std::vector<double>& rCentres);

/// A fit accepted with every point within the tolerance, and not a curved
/// surface so flat it stands in for a plane.
PrimitiveFit primitive_grow_fit(const PrimitiveGrowInput& rIn,
                                const std::vector<std::int64_t>& rFaces, std::uint32_t Types);

/// Whether every corner of face @p F lies within @p Within of the fit.
bool primitive_grow_face_fits(const PrimitiveGrowInput& rIn, std::int64_t F,
                              const PrimitiveFit& rFit, double Within);

/// Grow a region from @p Faces across smooth edges, within the faces' smooth
/// area, over faces with @p rPatch -1, while every corner lies within the growth
/// tolerance; the primitive is refitted each time the region doubles past
/// @p FittedOn faces. @p rMark is scratch, all zero on entry and exit.
std::vector<std::int64_t> primitive_grow(const PrimitiveGrowInput& rIn,
                                         std::vector<std::int64_t> Faces, PrimitiveFit Fit,
                                         std::size_t FittedOn,
                                         const std::vector<std::int64_t>& rPatch,
                                         std::vector<std::uint8_t>& rMark);

/// Trim a grown region to the cells its own refit explains within @p Within,
/// keeping the part connected to its first cell, three times over.
void primitive_grow_trim(const PrimitiveGrowInput& rIn, std::vector<std::int64_t>& rRegion,
                         std::uint32_t Types, double Within, std::vector<std::uint8_t>& rMark);

/// The faces' summed area.
double primitive_grow_area_of(const PrimitiveGrowInput& rIn,
                              const std::vector<std::int64_t>& rFaces);

/// What `primitive_grow_from_seed` grew.
struct PrimitiveGrowRegion {
    std::vector<std::int64_t> mFaces;  ///< empty when neither surface holds
    bool mCurved = false;
    std::vector<std::int64_t> mHood;  ///< the seed's neighbourhood the curved fit started on
};

/**
 * @brief A plane through the seed's own facet and a curved surface (of
 * @p CurvedTypes) through its two-ring neighbourhood of up to 32 faces, each
 * grown, trimmed to what its refit explains within @p ExactShare of the
 * tolerance (`segment` passes `kPrimitiveGrowExactShare`) and kept only when
 * that refit holds and a curved one has
 * three points per parameter; the curved region wins when its area is over
 * 1.5 times the plane's, the plane when it has two faces and at least
 * @p MinPlaneShare of @p AreaTotal. The step `segment`'s primitive method runs
 * per seed.
 */
PrimitiveGrowRegion primitive_grow_from_seed(const PrimitiveGrowInput& rIn, std::int64_t Seed,
                                             bool TryPlane, std::uint32_t CurvedTypes,
                                             double AreaTotal, double MinPlaneShare,
                                             double ExactShare,
                                             const std::vector<std::int64_t>& rPatch,
                                             std::vector<std::uint8_t>& rMark);

/// The face whose centroid lies nearest @p pPoint (ties to the lowest
/// index), or -1 for none: `select`'s seed lookup.
std::int64_t nearest_centroid_face(const std::vector<double>& rCentroid, const double* pPoint);

}  // namespace detail
}  // namespace meshioplusplus
