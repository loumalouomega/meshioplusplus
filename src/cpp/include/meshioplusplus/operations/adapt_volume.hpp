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
//  Attribution:     Implemented from published descriptions only:
//                   Dapogny, Dobrzynski and Frey, "Three-dimensional adaptive
//                   domain remeshing, implicit domain meshing, and
//                   applications to free and moving boundary problems",
//                   J. Comput. Phys. 262 (2014) (metric edge length, the
//                   split/collapse/swap/move cycle, Hausdorff-bounded surface
//                   changes); Klingner and Shewchuk, "Aggressive Tetrahedral
//                   Mesh Improvement", IMR 2007, and Klingner's thesis
//                   UCB/EECS-2008-145 (schedule and quality measures);
//                   Freitag and Ollivier-Gooch, IJNME 40(21), 1997; Hu et al.,
//                   "Fast Tetrahedral Meshing in the Wild", ACM TOG 39(4),
//                   2020 (conformal AMIPS); Jiang et al., "Declarative
//                   Specification for Unstructured Mesh Editing Algorithms",
//                   ACM TOG 41(6), 2022 (precondition / edit / postcondition
//                   / rollback); Botsch and Kobbelt 2004 (the 4/5 and 4/3
//                   length thresholds); Alliez et al., "Variational
//                   Tetrahedral Meshing", SIGGRAPH 2005 (ODT). No code of a
//                   copyleft or non-commercial mesher (MMG, CGAL, TetGen,
//                   Gmsh, Netgen, Mesquite, TetWild, Stellar) is read or
//                   vendored.
//

#pragma once

/**
 * @file operations/adapt_volume.hpp
 * @brief Isotropic, size-field-driven remeshing of a tetrahedral mesh:
 *        edge collapse, edge split, swaps and smoothing under a target size,
 *        a gradation and a boundary Hausdorff bound (doc/adapt_volume.md).
 *
 * The input is a valid tetrahedral mesh, typically `octree_mesh`'s graded
 * output, which has far more tetrahedra than its boundary resolution needs.
 * Every vertex gets a target size `h`: the boundary keeps its own resolution
 * (or the caller's `size_array`, or a curvature bound), the interior is free
 * to grow to `hmax`, and the field is graded so neighbouring sizes differ by
 * at most `hgrad - 1` per unit length. Each iteration then collapses the
 * edges shorter than 4/5 in that metric, splits those longer than 4/3, swaps
 * (2-3, 3-2, edge removal and multi-face removal, `optimize_volume`'s kernels)
 * and smooths (ODT or AMIPS in the interior, tangential moves projected back
 * onto the input skin on the boundary).
 *
 * Every local edit runs its preconditions, applies the edit, checks its
 * postconditions and is rolled back when one fails: no tetrahedron is ever
 * inverted or flat (exact `orient3d`), the least quality of the tetrahedra an
 * edit touches never drops below `min(old least, quality_min)`, the boundary
 * stays within `hausd` of the input boundary with its normals turning by at
 * most `max_normal_deviation`, corners never move, ridge vertices move only
 * along their chain and surface vertices only within their patch, and the
 * interfaces between materials and regions are kept as internal boundaries.
 * With `preserve_boundary` the boundary vertices and triangles are left
 * bit for bit.
 *
 * Every pass is serial in a fixed order (only per-element pure computations
 * run in parallel, into disjoint slots), so the output is byte-identical
 * across mesh backends, parallel backends and thread counts.
 */

// System includes
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The quality `adapt_volume` guards and reports.
enum class AdaptVolumeQuality : std::int32_t {
    /// The least dihedral angle, in degrees (the default).
    MinDihedral = 0,
    /// The mean ratio `12 (3 V)^(2/3) / sum(l^2)`, 1 for the regular tetrahedron.
    MeanRatio = 1,
    /// The conformal AMIPS energy, guarded as its reciprocal (the mean ratio);
    /// the interior smoothing then minimises the AMIPS energy of each star
    /// instead of moving to the ODT point.
    Amips = 2,
};

/// Options for `adapt_volume`.
struct AdaptVolumeOptions {
    /// The smallest target edge length; <= 0 takes the smallest boundary size.
    double mHmin = 0.0;
    /// The largest target edge length; <= 0 takes an eighth of the bounding-box diagonal.
    double mHmax = 0.0;
    /// The gradation: neighbouring sizes differ by at most `(hgrad - 1) |uv|`.
    double mHgrad = 1.3;
    /// How far the boundary may move from the input boundary. <= 0 takes 1e-3
    /// of the bounding-box diagonal; with `mHausdRelative`, a fraction of it.
    double mHausd = 0.0;
    /// Read `mHausd` as a fraction of the bounding-box diagonal.
    bool mHausdRelative = false;
    /// A float `point_data` array of target sizes (one component); empty for none.
    std::string mSizeArray;
    /// Bound the boundary size by `sqrt(8 hausd R)`, R the local radius of curvature.
    bool mCurvatureSizing = false;
    /// Boundary edges whose faces turn by more than this many degrees are ridges.
    double mFeatureAngle = 30.0;
    /// A ridge vertex where its chain turns by more than this many degrees is a corner.
    double mCornerAngle = 30.0;
    /// The most a boundary face's normal may turn in one edit, in degrees.
    double mMaxNormalDeviation = 30.0;
    /// Freeze the boundary: its vertices and triangles are left bit for bit.
    bool mPreserveBoundary = false;
    /// Keep sharp edges (ridges) and corners of the boundary.
    bool mPreserveFeatures = true;
    /// The quality the edits guard and the report measures.
    AdaptVolumeQuality mQualityMetric = AdaptVolumeQuality::MinDihedral;
    /// The quality floor of an edit: it may lower its star's least quality
    /// only down to this. < 0 takes 10 degrees (MinDihedral) or 0.2.
    double mQualityMin = -1.0;
    /// Iterations of the collapse / split / swap / smooth cycle.
    int mMaxIterations = 10;
    /// Stop early once the least dihedral angle reaches this; 0 disables.
    double mTargetMinDihedral = 0.0;
    /// Per point, non-zero for a vertex that never moves and is never removed.
    std::vector<std::uint8_t> mFrozen;
    /// Integer cell_data arrays that name materials, kept by parent; empty
    /// takes the known names (`octree_mesh:material`, `medit:ref`,
    /// `gmsh:physical`, `gmsh:geometrical`, `cell_tags`).
    std::vector<std::string> mMaterialArrays;
    /// The four passes, for debugging and ablation.
    bool mEnableCollapse = true;
    bool mEnableSplit = true;
    bool mEnableSwap = true;
    bool mEnableSmooth = true;
};

/// The output of `adapt_volume`.
struct AdaptVolumeResult {
    /// The adapted mesh: one `tetra` block per input block.
    Mesh mMesh;
    std::int64_t mNumTetsBefore = 0;
    std::int64_t mNumTetsAfter = 0;
    std::int64_t mNumPointsBefore = 0;
    std::int64_t mNumPointsAfter = 0;
    /// Edits applied.
    std::int64_t mNumCollapsed = 0;
    std::int64_t mNumSplit = 0;
    std::int64_t mNumSwapped = 0;
    std::int64_t mNumSmoothed = 0;
    /// Candidate edits refused, by the first guard that refused them.
    std::int64_t mNumRejectedByLink = 0;
    std::int64_t mNumRejectedByInversion = 0;
    std::int64_t mNumRejectedByQuality = 0;
    std::int64_t mNumRejectedByLength = 0;
    std::int64_t mNumRejectedByClass = 0;
    std::int64_t mNumRejectedByHausdorff = 0;
    std::int64_t mNumRejectedByNormal = 0;
    /// Dihedral angles, in degrees.
    double mMinDihedralBefore = 0.0;
    double mMinDihedralAfter = 0.0;
    double mMaxDihedralBefore = 0.0;
    double mMaxDihedralAfter = 0.0;
    /// The mean of the chosen quality (the mean ratio for AMIPS).
    double mMeanQualityBefore = 0.0;
    double mMeanQualityAfter = 0.0;
    /// The symmetric Hausdorff distance between the output and the input boundary.
    double mMaxBoundaryDeviation = 0.0;
    /// The sizes and the Hausdorff bound the run used.
    double mHmin = 0.0;
    double mHmax = 0.0;
    double mHausd = 0.0;
    /// Total volumes, and the change of each material class (input order).
    double mVolumeBefore = 0.0;
    double mVolumeAfter = 0.0;
    std::vector<std::pair<std::string, double>> mRegionVolumeChange;
    /// The classification of the input.
    std::int64_t mNumCorners = 0;
    std::int64_t mNumRidgeVertices = 0;
    std::int64_t mNumChains = 0;
    std::int64_t mNumPatches = 0;
    /// Iterations run.
    int mIterations = 0;
    /// Seconds per stage (`setup`, `collapse`, `split`, `swap`, `smooth`,
    /// `polish`, `output`), each summed over the iterations.
    std::vector<std::pair<std::string, double>> mTimings;
    /// What was dropped or could not be kept.
    std::vector<std::string> mWarnings;
};

/**
 * @brief Parse an `AdaptVolumeQuality` from `"min_dihedral"`, `"mean_ratio"`
 *        or `"amips"`.
 * @throws std::invalid_argument naming the accepted values.
 */
MESHIOPLUSPLUS_API AdaptVolumeQuality adapt_volume_quality_from_name(const std::string& rName);

/// The name `adapt_volume_quality_from_name` parses back.
MESHIOPLUSPLUS_API const char* adapt_volume_quality_name(AdaptVolumeQuality Quality);

/**
 * @brief Remesh a tetrahedral mesh to a graded isotropic size field.
 * @param rMesh a tetra-only volume mesh (one or more `tetra` blocks).
 * @param rOptions sizes, gradation, Hausdorff bound, feature policy, quality
 *        floor and pass switches.
 * @return the adapted mesh and the report.
 * @throws std::invalid_argument on a non-`tetra` or ragged block, a mesh with
 *         no tetrahedra, an inverted or flat input tetrahedron, a face shared
 *         by more than two tetrahedra, a bad option, or an unknown
 *         `size_array` / material array.
 */
MESHIOPLUSPLUS_API AdaptVolumeResult adapt_volume(const Mesh& rMesh,
                                                  const AdaptVolumeOptions& rOptions = {});

}  // namespace meshioplusplus
