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
// octree_mesh: option checks, the shared preparation (skin, winding-number
// tree, distance query, size field, balanced octree) and the dispatch to the
// outputs. See operations/octree_mesh.hpp for the contract. Anonymous-namespace
// helpers are prefixed `octm_` (the amalgamation concatenates every
// translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "octree_mesh_internal.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

[[noreturn]] void octm_fail(const std::string& rWhat) {
    throw std::invalid_argument(std::string(detail::kOctreeMeshPrefix) + rWhat);
}

void octm_check_options(const OctreeMeshOptions& rO) {
    auto finite_nonneg = [](double v) { return std::isfinite(v) && v >= 0.0; };
    if (!finite_nonneg(rO.mMinSize))
        octm_fail("min_size must be finite and not negative");
    if (!finite_nonneg(rO.mMaxSize))
        octm_fail("max_size must be finite and not negative");
    if (rO.mMinSize > 0.0 && rO.mMaxSize > 0.0 && rO.mMinSize > rO.mMaxSize)
        octm_fail("min_size must not exceed max_size");
    if (rO.mMaxLevel < -1 || rO.mMaxLevel > 19)
        octm_fail("max_level must be -1 (automatic) or in [0, 19]");
    if (!finite_nonneg(rO.mPadding) || !finite_nonneg(rO.mPaddingRelative))
        octm_fail("padding and padding_relative must be finite and not negative");
    if (!(rO.mCurvatureAngle >= 0.0 && rO.mCurvatureAngle < 180.0))
        octm_fail("curvature_angle must be in [0, 180) degrees");
    if (!(rO.mFeatureAngle >= 0.0 && rO.mFeatureAngle <= 180.0))
        octm_fail("feature_angle must be in [0, 180] degrees");
    if (!(rO.mCornerAngle >= 0.0 && rO.mCornerAngle <= 180.0))
        octm_fail("corner_angle must be in [0, 180] degrees");
    if (!finite_nonneg(rO.mFeatureSize))
        octm_fail("feature_size must be finite and not negative");
    if (!finite_nonneg(rO.mGapCells))
        octm_fail("gap_cells must be finite and not negative");
    if (rO.mGradingLayers < 0 || rO.mGradingLayers > 8)
        octm_fail("grading_layers must be in [0, 8]");
    if (!finite_nonneg(rO.mSurfaceSize))
        octm_fail("surface_size must be finite and not negative");
    if (!(rO.mWindingThreshold > 0.0 && rO.mWindingThreshold < 1.0))
        octm_fail("winding_threshold must be in (0, 1)");
    if (!finite_nonneg(rO.mWindingBeta))
        octm_fail("winding_beta must be finite and not negative");
    if (!(rO.mAlphaLong > 0.0 && rO.mAlphaLong < 0.5) ||
        !(rO.mAlphaShort > 0.0 && rO.mAlphaShort < 0.5))
        octm_fail("alpha_long and alpha_short must be in (0, 0.5)");
    if (rO.mFitIterations < 0)
        octm_fail("fit_iterations must not be negative");
    if (rO.mUntangleIterations < 0)
        octm_fail("untangle_iterations must not be negative");
    if (rO.mMaxCells < 1 || rO.mMaxTets < 1)
        octm_fail("max_cells and max_tets must be positive");
    if (rO.mBounds) {
        const auto& b = *rO.mBounds;
        for (std::size_t k = 0; k < 3; ++k)
            if (!(std::isfinite(b[k]) && std::isfinite(b[k + 3]) && b[k] < b[k + 3]))
                octm_fail("bounds must be finite with lo < hi on every axis");
    }
}

// Prepare everything the outputs share, in place (the trees point into the skin).
void octm_prepare(const Mesh& rMesh, const OctreeMeshOptions& rO, detail::OctreeMeshSetup& rS,
                  OctreeMeshResult& rResult, detail::OctreeMeshClock& rClock) {
    using namespace detail;
    rS.mSkin = octree_skin(rMesh, rO.mSurfaceRegion, kOctreeMeshPrefix);
    const OctreeSkin& skin = rS.mSkin;
    rS.mDiagonal = vec3_norm(vec3_sub(skin.mHi, skin.mLo));
    if (!(rS.mDiagonal > 0.0))
        octm_fail("the skin has no extent");
    rS.mOnSkin = 1e-12 * rS.mDiagonal;
    rS.mBvh = build_tri_bvh(skin.mSoup);
    SurfaceDistanceOptions dopt;
    dopt.mWatertightCheck = SdfWatertightCheck::Off;
    rS.mQuery = build_distance_query(skin.mSoup, dopt);
    // A skin wound inside out has negative volume: its winding numbers are -1
    // inside, so they are negated.
    double volume = 0.0;
    for (std::size_t t = 0; t < skin.mSoup.NumTriangles(); ++t) {
        const Vec3* c = &skin.mSoup.mCorners[3 * t];
        volume += vec3_dot(c[0], vec3_cross(c[1], c[2])) / 6.0;
    }
    if (volume < 0.0) {
        rS.mWindingSign = -1.0;
        octree_mesh_warn(rResult, std::string(kOctreeMeshPrefix) +
                                      "the skin is wound inside out (negative volume); its "
                                      "winding numbers are negated");
    }
    rClock.Lap(rResult, "skin");

    // Sizes.
    const double hmax = rO.mMaxSize > 0.0 ? rO.mMaxSize : rS.mDiagonal / 8.0;
    const double hmin = rO.mMinSize > 0.0 ? rO.mMinSize : hmax / 32.0;
    if (hmin > hmax)
        octm_fail(
            "min_size must not exceed max_size (max_size defaults to an eighth of the "
            "skin's diagonal)");
    rS.mFeatures = octree_skin_features(
        skin, rO.mPreserveFeatures || rO.mFeatureSize > 0.0 ? rO.mFeatureAngle : 0.0);
    SizeFieldSpec spec;
    spec.mMinSize = hmin;
    spec.mMaxSize = hmax;
    spec.mCurvatureAngleDeg = rO.mCurvatureAngle;
    spec.mFeatureSize =
        rO.mFeatureAngle > 0.0 ? (rO.mFeatureSize > 0.0 ? rO.mFeatureSize : hmin) : 0.0;
    spec.mGapCells = rO.mGapCells;
    spec.mGapOutside = rO.mGapOutside;
    std::vector<double> user;
    if (!rO.mSizeArray.empty()) {
        const Mesh& surface = skin.Surface(rMesh);
        if (!surface.HasPointData(rO.mSizeArray))
            octm_fail("size_array '" + rO.mSizeArray + "' is not a point_data array");
        const NDArray& a = surface.PointData(rO.mSizeArray);
        if (a.Size() != surface.NumPoints())
            octm_fail("size_array '" + rO.mSizeArray + "' must have one value per point");
        user.resize(a.Size());
        for (std::size_t i = 0; i < a.Size(); ++i)
            user[i] = read_double(a, i);
        spec.pUserSize = user.data();
    }
    rS.mField = octree_size_field(skin, rS.mBvh, rS.mFeatures, spec);
    rClock.Lap(rResult, "size_field");

    // The frame, over the padded box.
    const bool uniform = rO.mSkinGrading == OctreeMeshSkinGrading::Uniform;
    const double hs =
        rO.mSurfaceSize > 0.0 ? std::clamp(rO.mSurfaceSize, hmin, hmax) : rS.mField.mSkinMinSize;
    double pad = rO.mPadding + rO.mPaddingRelative * rS.mDiagonal;
    Vec3 lo = skin.mLo, hi = skin.mHi;
    if (rO.mBounds) {
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = (*rO.mBounds)[k];
            hi[k] = (*rO.mBounds)[k + 3];
        }
    } else {
        if (uniform && pad < 3.0 * hs) {
            log::info("{}padding raised to {} (three skin sizes) for uniform skin grading",
                      kOctreeMeshPrefix, 3.0 * hs);
            pad = 3.0 * hs;
        }
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] -= pad;
            hi[k] += pad;
        }
    }
    rS.mFrame = octree_frame(lo, hi, hmin, hmax, rO.mMaxLevel, kOctreeMeshPrefix);
    if (rS.mFrame.mMinSize > hmin * (1.0 + 1e-12))
        octree_mesh_warn(rResult, std::string(kOctreeMeshPrefix) + "max_level " +
                                      std::to_string(rO.mMaxLevel) + " coarsens min_size to " +
                                      std::to_string(rS.mFrame.mMinSize));
    const std::int32_t surface_level = uniform ? rS.mFrame.LevelFor(hs) : -1;
    if (uniform && rS.mFrame.mRootSize / static_cast<double>(std::int64_t{1} << surface_level) >
                       hs * (1.0 + 1e-9))
        octree_mesh_warn(rResult, std::string(kOctreeMeshPrefix) +
                                      "uniform skin grading: the skin size is coarser than asked");

    // The octree: refined by the size field, the cut leaves dilated, balanced.
    rS.mOctree = octree_adaptive_build(rS.mFrame, rS.mField.mPrims, surface_level, rO.mMaxCells);
    const std::int32_t rings =
        uniform ? std::max<std::int32_t>(2, static_cast<std::int32_t>(rO.mGradingLayers))
                : static_cast<std::int32_t>(rO.mGradingLayers);
    if (rings > 0) {
        std::vector<std::int64_t> origin;
        linear_octree_dilate(rS.mOctree.mTree, rS.mOctree.mCut, rings, rO.mMaxCells, &origin);
        octree_adaptive_refresh(rS.mOctree, rS.mFrame, rS.mField.mPrims, origin,
                                std::move(rS.mOctree.mCand));
    }
    std::vector<std::int64_t> origin;
    linear_octree_balance(rS.mOctree.mTree, OctreeAdjacency::Full, rO.mMaxCells, &origin);
    octree_adaptive_refresh(rS.mOctree, rS.mFrame, rS.mField.mPrims, origin,
                            std::move(rS.mOctree.mCand));
    rResult.mOrigin = rS.mFrame.mOrigin;
    rResult.mRootSize = rS.mFrame.mRootSize;
    rResult.mRootCells = rS.mFrame.mRoots;
    rResult.mMaxLevel = rS.mFrame.mMaxLevel;
    rResult.mSurfaceLevel = surface_level;
    rResult.mMinSize = rS.mFrame.mMinSize;
    rResult.mMaxSize = hmax;
    rClock.Lap(rResult, "octree");
}

// The output-independent part of the report.
void octm_report_mesh(const OctreeMeshOptions& rO, OctreeMeshResult& rResult) {
    const Mesh& m = rResult.mMesh;
    rResult.mNumPoints = static_cast<std::int64_t>(m.NumPoints());
    rResult.mNumCells = 0;
    rResult.mNumHexahedra = 0;
    rResult.mNumTetrahedra = 0;
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        const auto block = m.Cells(b);
        const auto n = static_cast<std::int64_t>(block.NumCells());
        rResult.mNumCells += n;
        if (block.Type() == "hexahedron")
            rResult.mNumHexahedra += n;
        else if (block.Type() == "tetra")
            rResult.mNumTetrahedra += n;
    }
    rResult.mQuality = compute_quality(m);
    rResult.mNumInverted = rResult.mQuality.mNumInverted;
    for (const auto& [name, summary] : rResult.mQuality.mMetrics) {
        if (summary.mCount == 0)
            continue;
        if (name == "quality:scaled_jacobian")
            rResult.mMinScaledJacobian = summary.mMin;
        else if (name == "quality:min_dihedral")
            rResult.mMinDihedral = summary.mMin;
        else if (name == "quality:max_dihedral")
            rResult.mMaxDihedral = summary.mMax;
    }
    if (rResult.mNumTetrahedra > 0 && rResult.mMinDihedral < 5.0)
        detail::octree_mesh_warn(rResult, std::string(detail::kOctreeMeshPrefix) +
                                              "the smallest dihedral angle is " +
                                              std::to_string(rResult.mMinDihedral) + " degrees");
    if (rResult.mNumBoundaryEdgesOpen > 0 || rResult.mNumNonManifoldBoundaryEdges > 0)
        detail::octree_mesh_warn(rResult,
                                 std::string(detail::kOctreeMeshPrefix) + "the boundary has " +
                                     std::to_string(rResult.mNumBoundaryEdgesOpen) + " open and " +
                                     std::to_string(rResult.mNumNonManifoldBoundaryEdges) +
                                     " non-manifold edges");
    if (rResult.mNumInverted > 0)
        detail::octree_mesh_warn(rResult, std::string(detail::kOctreeMeshPrefix) +
                                              std::to_string(rResult.mNumInverted) +
                                              " cells have a non-positive volume");
    if (rO.mAttachQuality)
        rResult.mMesh = attach_quality(rResult.mMesh);
}

}  // namespace

namespace detail {

void OctreeMeshClock::Lap(OctreeMeshResult& rResult, const char* pStage) {
    const auto now = std::chrono::steady_clock::now();
    rResult.mTimings.emplace_back(pStage, std::chrono::duration<double>(now - mStart).count());
    mStart = now;
}

void octree_mesh_warn(OctreeMeshResult& rResult, const std::string& rMessage) {
    log::warn("{}", rMessage);
    rResult.mWarnings.push_back(rMessage);
}

std::vector<double> octree_mesh_winding(const OctreeMeshSetup& rSetup,
                                        const std::vector<Vec3>& rPoints,
                                        const OctreeMeshOptions& rOptions,
                                        OctreeMeshResult& rResult) {
    std::vector<double> w = fast_windings(rSetup.mBvh, rPoints, rOptions.mWindingBeta);
    for (double& v : w)
        v *= rSetup.mWindingSign;
    if (rOptions.mWindingBeta > 0.0) {
        std::vector<std::size_t> near;
        for (std::size_t i = 0; i < w.size(); ++i)
            if (std::abs(w[i] - rOptions.mWindingThreshold) < 0.1)
                near.push_back(i);
        const double beta = std::max(8.0, 4.0 * rOptions.mWindingBeta);
        parallel_for(
            near.size(),
            [&](std::size_t k) {
                w[near[k]] =
                    rSetup.mWindingSign * fast_winding(rSetup.mBvh, rPoints[near[k]], beta);
            },
            64);
        rResult.mNumWindingRechecks += static_cast<std::int64_t>(near.size());
    }
    return w;
}

std::vector<Vec3> octree_mesh_leaf_centres(const OctreeFrame& rFrame, const LinearOctree& rTree) {
    std::vector<Vec3> c(rTree.Size());
    const double half_fine = 0.5 * rFrame.Fine();
    parallel_for(rTree.Size(), [&](std::size_t i) {
        const auto lo = rTree.Lo(i);
        const std::int64_t s = rTree.Span(i);
        for (std::size_t k = 0; k < 3; ++k)
            c[i][k] = rFrame.mOrigin[k] + static_cast<double>(2 * lo[k] + s) * half_fine;
    });
    return c;
}

void octree_mesh_signed_distance(const OctreeMeshSetup& rSetup, const std::vector<Vec3>& rPoints,
                                 const OctreeMeshOptions& rOptions, OctreeMeshResult& rResult,
                                 std::vector<double>& rDistance, std::vector<Vec3>& rClosest) {
    const std::vector<ClosestPointHit> hits = query_closest_points(rSetup.mQuery, rPoints);
    const std::vector<double> w = octree_mesh_winding(rSetup, rPoints, rOptions, rResult);
    rDistance.resize(rPoints.size());
    rClosest.resize(rPoints.size());
    for (std::size_t i = 0; i < rPoints.size(); ++i) {
        rClosest[i] = hits[i].mPoint;
        const double d = hits[i].mDistance;
        rDistance[i] = d <= rSetup.mOnSkin ? 0.0 : (w[i] > rOptions.mWindingThreshold ? -d : d);
    }
}

void octree_mesh_report_tree(const LinearOctree& rTree, OctreeMeshResult& rResult) {
    rResult.mNumLeaves = static_cast<std::int64_t>(rTree.Size());
    rResult.mLeavesPerLevel.assign(static_cast<std::size_t>(rTree.mMaxLevel) + 1, 0);
    rResult.mFinestLevel = 0;
    for (std::uint8_t l : rTree.mLevel) {
        ++rResult.mLeavesPerLevel[l];
        rResult.mFinestLevel = std::max<std::int64_t>(rResult.mFinestLevel, l);
    }
}

}  // namespace detail

OctreeMeshOutput octree_mesh_output_from_name(const std::string& rName) {
    if (rName == "castellated")
        return OctreeMeshOutput::Castellated;
    if (rName == "tet" || rName == "tetra")
        return OctreeMeshOutput::Tetra;
    if (rName == "hex" || rName == "hexa")
        return OctreeMeshOutput::Hexa;
    octm_fail("unknown output '" + rName + "' (expected 'castellated', 'tet' or 'hex')");
}

const char* octree_mesh_output_name(OctreeMeshOutput Output) {
    switch (Output) {
        case OctreeMeshOutput::Castellated:
            return "castellated";
        case OctreeMeshOutput::Tetra:
            return "tet";
        case OctreeMeshOutput::Hexa:
            return "hex";
    }
    return "tet";
}

OctreeMeshConformity octree_mesh_conformity_from_name(const std::string& rName) {
    if (rName == "stuffing")
        return OctreeMeshConformity::Stuffing;
    if (rName == "cleaving")
        return OctreeMeshConformity::Cleaving;
    octm_fail("unknown conformity '" + rName + "' (expected 'stuffing' or 'cleaving')");
}

OctreeMeshSkinGrading octree_mesh_skin_grading_from_name(const std::string& rName) {
    if (rName == "graded")
        return OctreeMeshSkinGrading::Graded;
    if (rName == "uniform")
        return OctreeMeshSkinGrading::Uniform;
    octm_fail("unknown skin_grading '" + rName + "' (expected 'graded' or 'uniform')");
}

OctreeMeshFill octree_mesh_fill_from_name(const std::string& rName) {
    if (rName == "all")
        return OctreeMeshFill::All;
    if (rName == "inside")
        return OctreeMeshFill::Inside;
    if (rName == "cover")
        return OctreeMeshFill::Cover;
    octm_fail("unknown fill '" + rName + "' (expected 'all', 'inside' or 'cover')");
}

OctreeMeshResult octree_mesh(const Mesh& rSkin, const OctreeMeshOptions& rOptions) {
    octm_check_options(rOptions);
    OctreeMeshResult result;
    result.mOutput = rOptions.mOutput;
    detail::OctreeMeshClock clock;
    auto setup = std::make_unique<detail::OctreeMeshSetup>();
    octm_prepare(rSkin, rOptions, *setup, result, clock);
    switch (rOptions.mOutput) {
        case OctreeMeshOutput::Castellated:
            detail::octree_mesh_castellated(*setup, rOptions, result);
            break;
        case OctreeMeshOutput::Tetra: {
            detail::OctreeMeshTets tets;
            detail::octree_mesh_stuffed_lattice(*setup, rOptions, result, tets);
            clock.Lap(result, "stuffing");
            detail::octree_mesh_snap_features(*setup, rOptions, result, tets);
            clock.Lap(result, "features");
            detail::octree_mesh_optimize(*setup, rOptions, result, tets);
            clock.Lap(result, "optimize");
            detail::octree_mesh_emit_tets(*setup, rOptions, result, tets);
            result.mGuaranteeApplies = rOptions.mSkinGrading == OctreeMeshSkinGrading::Uniform &&
                                       result.mNumWarpRefused == 0 &&
                                       result.mNumFeatureMoves == 0 &&
                                       (result.mNumFlips == 0 || rOptions.mDihedralFloorGuard);
            break;
        }
        case OctreeMeshOutput::Hexa:
            detail::octree_mesh_hexahedra(*setup, rOptions, result);
            break;
    }
    clock.Lap(result, "output");
    octm_report_mesh(rOptions, result);
    clock.Lap(result, "quality");
    return result;
}

}  // namespace meshioplusplus
