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
 * @file patch_flatten.hpp
 * @brief Flattening a triangulated surface patch over its best-fit plane or an
 * unrolled cylinder, the fold test that says whether that flattening is a
 * parametrisation, and splitting a patch no flattening parametrises into
 * pieces that one does -- `fit_cad`'s freeform faces (roadmap 6.3.2), kept free
 * of OpenCASCADE so it is testable in every build.
 *
 * Not installed: a core-private helper of `operations/fit_cad.cpp`.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace meshioplusplus {
namespace detail {

using PatchVec = std::array<double, 3>;
using PatchTri = std::array<std::size_t, 3>;

/// A patch flattened: per vertex its `(s, t)` parameters and its offset `h`
/// from the plane or cylinder it is flattened over, with the map back.
struct PatchFlat {
    bool mCylinder = false;
    PatchVec mOrigin{}, mE1{}, mE2{}, mAxis{};  // plane: origin, in-plane axes, normal
    double mRadius = 0.0;
    double mTheta0 = 0.0;
    std::vector<double> mS, mT, mH;

    /// The point at parameters `(S, T)` and offset `H`.
    PatchVec At(double S, double T, double H) const;
};

/**
 * @brief Flatten over the best-fit plane (`Kind` 0) or unroll around the axis
 * the normals are most perpendicular to (`Kind` 1), the cylinder cut in the
 * widest gap between the points' angles.
 * @param rPts the patch's points.
 * @param rNormals its triangles' unit normals.
 * @param rAreas its triangles' areas.
 * @return false when that flattening is degenerate.
 */
bool patch_flatten(int Kind, const std::vector<PatchVec>& rPts,
                   const std::vector<PatchVec>& rNormals, const std::vector<double>& rAreas,
                   PatchFlat& rOut);

/// Whether a flattening folds: triangles of both orientations (beyond
/// `1e-9` of the largest), or no area at all.
bool patch_flat_folds(const PatchFlat& rFlat, const std::vector<PatchTri>& rTris);

/// Unit normals and areas of the triangles.
void patch_normals(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                   std::vector<PatchVec>& rNormals, std::vector<double>& rAreas);

/// Whether the plane or the cylinder flattens the patch without folding.
bool patch_flattenable(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris);

/// How `patch_split_until_flat` splits.
struct PatchSplitOptions {
    /// Halvings before a piece is given up on: at most `2^depth` pieces per patch.
    int mMaxDepth = 4;
    /// A piece with fewer triangles stays facets.
    std::size_t mMinTriangles = 16;
};

/**
 * @brief Split a patch no flattening parametrises into pieces one does.
 *
 * The triangles are halved at the area-weighted median of `n . e`, `e` the
 * direction their normals spread along most (the principal axis of the
 * area-weighted normal covariance): a cut along an iso-line of the normal, a
 * smooth curve on a smooth surface. Each half is split into its edge-connected
 * components, and each component that does not flatten is halved again, up to
 * `mMaxDepth` times. Deterministic: ties go by triangle index, and pieces are
 * numbered in the order of their lowest triangle.
 *
 * @return per triangle, its piece (from 0), or -1 for one left as facets
 *         (a piece under `mMinTriangles`, or still folding at the depth cap).
 */
std::vector<std::int64_t> patch_split_until_flat(const std::vector<PatchVec>& rPts,
                                                 const std::vector<PatchTri>& rTris,
                                                 const PatchSplitOptions& rOptions,
                                                 std::size_t& rNumPieces);

}  // namespace detail
}  // namespace meshioplusplus
