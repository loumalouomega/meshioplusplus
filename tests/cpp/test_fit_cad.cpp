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
// fit_cad: analytic CAD faces through the optional OpenCASCADE backend.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "analytic_surfaces.hpp"
#include "fit_cad_fixtures.hpp"
#include "mesh_fixtures.hpp"

#include "../../src/cpp/src/detail/nurbs_eval.hpp"

using namespace meshioplusplus;

namespace {

using fcfix::fitted;
using fcfix::saddle_box;

#ifdef MESHIOPLUSPLUS_HAS_OCCT
// The NURBS side channel holds every face of the solid, each trimmed, and every
// trimming curve, sampled over its active range, inside its face's surface
// parameter range: a consumer evaluates the surface there.
void expect_brep_consistent(const FitCadResult& rR) {
    ASSERT_EQ(rR.mBrep.mFaces.size(), static_cast<std::size_t>(rR.mNumFaces));
    for (const BrepFace& f : rR.mBrep.mFaces) {
        EXPECT_TRUE(f.mTrimmed) << "face " << f.mId;
        const NurbsPatch& patch = rR.mBrep.mPatches.at(static_cast<std::size_t>(f.mPatch));
        double lo[2], hi[2];
        for (std::size_t d = 0; d < 2; ++d) {
            const auto p = static_cast<std::size_t>(patch.mDegrees[d]);
            const std::vector<double>& k = patch.mKnots[d];
            lo[d] = k[p];
            hi[d] = k[k.size() - p - 1];
        }
        for (const BrepLoop& loop : f.mLoops)
            for (const BrepTrim& t : loop.mTrims) {
                const BrepCurve& c = t.mCurve;
                const std::size_t n = c.mKnots.size() - static_cast<std::size_t>(c.mDegree) - 1;
                for (int i = 0; i <= 8; ++i) {
                    const double u =
                        c.mActiveRange[0] + (c.mActiveRange[1] - c.mActiveRange[0]) * i / 8.0;
                    const std::size_t span = detail::nurbs_find_span(c.mDegree, c.mKnots, n, u);
                    double uv[3];
                    detail::nurbs_curve_point_in_span(c.mDegree, c.mKnots,
                                                      c.mControlPoints.As<double>(),
                                                      c.mWeights.As<double>(), 3, span, u, uv);
                    for (std::size_t d = 0; d < 2; ++d) {
                        const double slack = 1e-7 * (hi[d] - lo[d]);
                        EXPECT_GE(uv[d], lo[d] - slack) << "face " << f.mId;
                        EXPECT_LE(uv[d], hi[d] + slack) << "face " << f.mId;
                    }
                }
            }
    }
}
#endif

std::string read_bytes(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

}  // namespace

TEST(FitCad, CapabilityMatchesTheBuild) {
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    EXPECT_TRUE(fit_cad_has_occt());
#else
    EXPECT_FALSE(fit_cad_has_occt());
#endif
}

// The kernel-free backend runs in every build.
TEST(FitCad, NativeBackendWritesAnalyticFaces) {
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    const std::string path =
        (std::filesystem::temp_directory_path() / "native_plate.step").string();
    // A raw skin: segmented by the primitive method and fitted first.
    const FitCadResult r = fit_cad(analytic::filleted_plate(), path, o);
    EXPECT_EQ(r.mBackend, FitCadBackend::Native);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumSolids, 1);
    EXPECT_EQ(r.mNumFaces, 8);
    EXPECT_EQ(r.mNumPlane, 6);
    EXPECT_EQ(r.mNumCylinder, 2);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_EQ(r.mNumDemoted, 0);
    EXPECT_TRUE(std::isnan(r.mVolume));
    EXPECT_GT(r.mSkinVolume, 0.0);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    // The skin itself comes back, each cell with its face, type and patch.
    ASSERT_TRUE(r.mMesh.HasCellData(kFitCadFaceName));
    const NDArray& t = r.mMesh.CellData(kFitCadTypeName, 0);
    EXPECT_EQ(t.Size(), r.mMesh.Cells(0).NumCells());
    for (std::size_t i = 0; i < t.Size(); ++i) {
        const std::int64_t v = t.As<std::int64_t>()[i];
        EXPECT_TRUE(v == 1 || v == 3) << v;
    }
    const std::string text = read_bytes(path);
    std::remove(path.c_str());
    EXPECT_NE(text.find("CYLINDRICAL_SURFACE"), std::string::npos);
}

// The native backend is a write entry point: a provenance note an earlier,
// unrelated operation left on this thread must not reach its STEP header.
TEST(FitCad, NativeStepIgnoresEarlierProvenanceNotes) {
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    const Mesh skin = fitted(analytic::box_skin(3));
    const std::string path = ::testing::TempDir() + "native_provenance.step";
    detail::provenance_begin_write();
    fit_cad(skin, path, o);
    const std::string clean = read_bytes(path);
    detail::provenance_note("left-behind", "a note from an earlier operation");
    fit_cad(skin, path, o);
    const std::string after = read_bytes(path);
    std::remove(path.c_str());
    EXPECT_EQ(after.find("left-behind"), std::string::npos);
    EXPECT_EQ(clean, after);
}

// Refinement before building: an angle segmentation of the filleted plate
// leaves the fillet in one patch with a side it is tangent to, which no
// primitive explains; refined, both are analytic faces (v17.11.0).
TEST(FitCad, NativeBackendRefines) {
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    const Mesh skin = fitted(analytic::filleted_plate());
    const FitCadResult before = fit_cad(skin, "", o);
    o.mRefineSegmentation = true;
    const FitCadResult after = fit_cad(skin, "", o);
    EXPECT_FALSE(before.mRefined);
    EXPECT_TRUE(after.mRefined);
    EXPECT_GT(after.mNumRefineSplits, 0);
    EXPECT_GT(after.mExplainedAreaFraction, after.mInitialExplainedAreaFraction);
    EXPECT_NEAR(after.mExplainedAreaFraction, 1.0, 1e-12);
    EXPECT_GT(after.mNumPlane + after.mNumCylinder, before.mNumPlane + before.mNumCylinder);
    EXPECT_EQ(after.mNumCylinder, 2);
    EXPECT_TRUE(after.mValid);
    EXPECT_TRUE(after.mMesh.HasFieldData(kRefinePatchProxyName));
    EXPECT_FALSE(before.mMesh.HasFieldData(kRefinePatchProxyName));
}

TEST(FitCad, NativeBackendLeavesTypesOutAsPolygons) {
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTypes = 1u << 1;  // planes only: the bore and the fillet become polygons
    const FitCadResult r = fit_cad(analytic::filleted_plate(), "", o);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumPlane, 6);
    EXPECT_EQ(r.mNumCylinder, 0);
    EXPECT_GT(r.mNumFacetedFaces, 0);
}

TEST(FitCad, NativeBackendReportsAClosedSolid) {
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    const FitCadResult r = fit_cad(fitted(analytic::box_skin(3)), "", o);
    EXPECT_TRUE(r.mClosed);
    EXPECT_TRUE(std::isnan(r.mHausdorffMax));
    EXPECT_EQ(r.mNumExtendedEdges, 0);
    EXPECT_EQ(r.mNumSnappedJunctions, 0);
}

TEST(FitCad, BadFreeformGeneralOptionsAreRefused) {
    const auto bad = [](auto&& rSet) {
        FitCadOptions o;
        o.mBackend = FitCadBackend::Native;
        rSet(o);
        EXPECT_THROW(fit_cad(analytic::box_skin(2), "", o), std::invalid_argument);
    };
    bad([](FitCadOptions& o) { o.mFreeformTolerance = -1.0; });
    bad([](FitCadOptions& o) { o.mFreeformSmoothing = -1.0; });
    bad([](FitCadOptions& o) { o.mMaxFreeformPoles = 3; });
    bad([](FitCadOptions& o) { o.mMaxFreeformIterations = 0; });
    bad([](FitCadOptions& o) { o.mMaxFreeformSplits = -1; });
    bad([](FitCadOptions& o) { o.mSewTolerance = -1.0; });
    bad([](FitCadOptions& o) { o.mSewTolerance = std::numeric_limits<double>::infinity(); });
}

TEST(FitCad, BadBackendIsRefused) {
    FitCadOptions o;
    o.mBackend = static_cast<FitCadBackend>(5);
    EXPECT_THROW(fit_cad(analytic::box_skin(2), "", o), std::invalid_argument);
}

#ifndef MESHIOPLUSPLUS_HAS_OCCT

TEST(FitCad, AbsentBackendIsNamed) {
    FitCadOptions o;
    o.mBackend = FitCadBackend::Occt;
    try {
        fit_cad(analytic::box_skin(2), "", o);
        FAIL() << "expected Unsupported";
    } catch (const Unsupported& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("MESHIOPLUSPLUS_WITH_OCCT"), std::string::npos);
    }
    // Without OpenCASCADE, the automatic backend is the native one.
    EXPECT_EQ(fit_cad(analytic::box_skin(2)).mBackend, FitCadBackend::Native);
}

// The automatic fallback says so once per process (an earlier test may have
// used it already); asking for the native backend by name never warns.
TEST(FitCad, FallbackIsAnnouncedOnce) {
    const std::string notice = "kernel-free native backend is used";
    auto count = [&](const std::string& rText) {
        std::size_t n = 0;
        for (std::size_t at = rText.find(notice); at != std::string::npos;
             at = rText.find(notice, at + 1))
            ++n;
        return n;
    };
    testing::internal::CaptureStderr();
    fit_cad(analytic::box_skin(2));
    fit_cad(analytic::box_skin(2));
    EXPECT_LE(count(testing::internal::GetCapturedStderr()), 1u);
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    testing::internal::CaptureStderr();
    fit_cad(analytic::box_skin(2), "", o);
    EXPECT_EQ(count(testing::internal::GetCapturedStderr()), 0u);
}

#else

TEST(FitCad, FilletedPlateIsEightExactFaces) {
    // A raw skin with a fillet tangent to its sides and a bore: segmented by
    // the surfaces the facets lie on, so neither is lost to facets.
    const FitCadResult r = fit_cad(analytic::filleted_plate());
    EXPECT_EQ(r.mBackend, FitCadBackend::Occt);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumFaces, 8);
    EXPECT_EQ(r.mNumPlane, 6);
    EXPECT_EQ(r.mNumCylinder, 2);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_EQ(r.mNumDemoted, 0);
    // The exact solid: the plate less the fillet's corner and the bore.
    const double rc = 0.8, rb = 0.6;
    const double exact = (16.0 - (rc * rc - analytic::kPi * rc * rc / 4.0) - analytic::kPi * rb * rb) * 0.5;
    EXPECT_NEAR(r.mVolume, exact, 1e-6);
}

TEST(FitCad, FacetsLeftOutMergeIntoCoplanarFaces) {
    // With every type left out the box is all facets, written as one planar
    // face per coplanar group rather than one per triangle.
    FitCadOptions o;
    o.mTypes = 0;
    const FitCadResult r = fit_cad(fitted(analytic::box_skin(3)), "", o);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumFaces, 6);
    EXPECT_EQ(r.mNumFacetedFaces, 6);
    EXPECT_NEAR(r.mVolume, 1.0, 1e-9);
    o.mMergeFaces = false;
    EXPECT_EQ(fit_cad(fitted(analytic::box_skin(3)), "", o).mNumFaces, 6);
}

TEST(FitCad, BoxIsSixExactPlanes) {
    const FitCadResult r = fit_cad(fitted(analytic::box_skin(3)));
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumSolids, 1);
    EXPECT_EQ(r.mNumFaces, 6);
    EXPECT_EQ(r.mNumPlane, 6);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_EQ(r.mNumDemoted, 0);
    EXPECT_EQ(r.mNumExactEdges, 12);
    EXPECT_EQ(r.mNumApproxEdges, 0);
    EXPECT_NEAR(r.mVolume, 1.0, 1e-9);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-9);
    EXPECT_LT(r.mMaxDeviation, 1e-9);
    ASSERT_TRUE(r.mMesh.HasCellData(kFitCadFaceName));
    ASSERT_TRUE(r.mMesh.HasCellData(kFitCadTypeName));
    ASSERT_TRUE(r.mMesh.HasCellData(kFitCadPatchName));
    const NDArray& t = r.mMesh.CellData(kFitCadTypeName, 0);
    for (std::size_t i = 0; i < t.Size(); ++i)
        EXPECT_EQ(t.As<std::int64_t>()[i], 1);
}

TEST(FitCad, CappedCylinderGetsASeamedWall) {
    const FitCadResult r = fit_cad(fitted(analytic::capped_cylinder()));
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumExactEdges, 4);  // two circles, each closed boundary in two arcs
    EXPECT_EQ(r.mNumApproxEdges, 0);
    EXPECT_EQ(r.mNumPlane, 2);
    EXPECT_EQ(r.mNumCylinder, 1);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_EQ(r.mNumDemoted, 0);
    // The exact cylinder holds the polygon's volume within the chord error.
    EXPECT_NEAR(r.mVolume, analytic::kPi * 2.0, 0.02);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 0.01);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
}

TEST(FitCad, TurnedPartKeepsEveryPrimitive) {
    const Mesh fit = fitted(analytic::turned_part(), 1.0, 2e-3);
    const FitPrimitivesResult counts = fit_primitives(segment(analytic::turned_part(),
                                                              [] {
                                                                  SegmentOptions o;
                                                                  o.mCurvatureTolerance = 1.0;
                                                                  return o;
                                                              }())
                                                          .mMesh,
                                                      [] {
                                                          FitPrimitivesOptions o;
                                                          o.mTolerance = 2e-3;
                                                          return o;
                                                      }());
    const FitCadResult r = fit_cad(fit);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumDemoted, 0);
    EXPECT_EQ(r.mNumPlane, counts.mNumPlane);
    EXPECT_EQ(r.mNumCylinder, counts.mNumCylinder);
    EXPECT_EQ(r.mNumCone, counts.mNumCone);
    EXPECT_EQ(r.mNumSphere, counts.mNumSphere);
    EXPECT_EQ(r.mNumTorus, counts.mNumTorus);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 0.01);
    // A cone whose converted surface would stop short of its trims.
    expect_brep_consistent(r);
}

TEST(FitCad, ClosedTorusAndSphereAreNaturalFaces) {
    const FitCadResult t = fit_cad(fitted(analytic::torus()));
    EXPECT_TRUE(t.mValid);
    EXPECT_EQ(t.mNumTorus, 1);
    EXPECT_EQ(t.mNumFaces, 1);
    EXPECT_NEAR(t.mVolumeRatio, 1.0, 0.02);

    std::vector<std::pair<double, double> > prof;
    for (int i = 1; i < 16; ++i) {
        const double th = analytic::kPi * i / 16;
        prof.push_back({std::sin(th), -std::cos(th)});
    }
    const FitCadResult s = fit_cad(fitted(analytic::closed_revolve(prof, -1.0, 1.0)));
    EXPECT_TRUE(s.mValid);
    EXPECT_EQ(s.mNumSphere, 1);
    EXPECT_EQ(s.mNumFaces, 1);
    EXPECT_NEAR(s.mVolumeRatio, 1.0, 0.03);
}

TEST(FitCad, UnexplainedPatchComesBackFacetedWithoutFreeform) {
    FitCadOptions o;
    o.mFreeform = false;
    const FitCadResult r = fit_cad(fitted(saddle_box(6)), "", o);
    EXPECT_TRUE(r.mValid);
    EXPECT_GT(r.mNumFacetedFaces, 0);
    EXPECT_EQ(r.mNumBspline, 0);
    EXPECT_EQ(r.mNumPlane, 5);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-9);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
}

TEST(FitCad, UnexplainedPatchBecomesAFreeformFace) {
    // The saddle top no primitive explains is a B-spline face (roadmap 6.3.2).
    const FitCadResult r = fit_cad(fitted(saddle_box(24)));
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumBspline, 1);
    EXPECT_EQ(r.mNumPlane, 5);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_EQ(r.mNumFaces, 6);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-3);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    // The solid comes back as NURBS too, its faces numbered as fit_cad:face.
    ASSERT_EQ(r.mBrep.mFaces.size(), 6u);
    expect_brep_consistent(r);
    EXPECT_EQ(r.mBrep.mFaces[0].mId, 1);
    bool any_bspline_type = false;
    const NDArray& types = r.mMesh.CellData(kFitCadTypeName, 0);
    for (std::size_t i = 0; i < types.Size(); ++i)
        any_bspline_type = any_bspline_type || types.As<std::int64_t>()[i] == kFitCadBsplineType;
    EXPECT_TRUE(any_bspline_type);
}

TEST(FitCad, CurledPatchSplitsIntoFreeformFaces) {
    // The wavy wall folds over its plane and, closed around the axis, over its
    // unrolled cylinder too: it is split into pieces that each unroll, and
    // each becomes a B-spline face (roadmap 6.3.2).
    const Mesh skin = fitted(analytic::wavy_capped());
    const FitCadResult r = fit_cad(skin);
    EXPECT_TRUE(r.mValid);
    EXPECT_GE(r.mNumSplitPatches, 2);
    EXPECT_EQ(r.mNumBspline, r.mNumSplitPatches);
    EXPECT_EQ(r.mNumPlane, 2);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-2);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    // The pieces carry new patch labels, after the largest existing one.
    std::int64_t max_in = -1, max_out = -1;
    for (std::size_t b = 0; b < skin.NumCellBlocks(); ++b) {
        const NDArray& a = skin.CellData(kSegmentPatchName, b);
        for (std::size_t i = 0; i < a.Size(); ++i)
            max_in = std::max(max_in, a.As<std::int64_t>()[i]);
    }
    const NDArray& out = r.mMesh.CellData(kFitCadPatchName, 0);
    for (std::size_t i = 0; i < out.Size(); ++i)
        max_out = std::max(max_out, out.As<std::int64_t>()[i]);
    EXPECT_EQ(max_out, max_in + r.mNumSplitPatches);
}

TEST(FitCad, CurledPatchStaysFacetedWithoutSplitting) {
    FitCadOptions o;
    o.mFreeformSplit = false;
    const FitCadResult r = fit_cad(fitted(analytic::wavy_capped()), "", o);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumSplitPatches, 0);
    EXPECT_EQ(r.mNumBspline, 0);
    EXPECT_GT(r.mNumFacetedFaces, 0);
}

TEST(FitCad, InclinedBossFilletIsARollingBallFace) {
    // The foot of a leaning boss rounded with a constant radius: a canal
    // surface along an ellipse, no primitive (roadmap 6.3.1).
    const FitCadResult r = fit_cad(fitted(analytic::inclined_boss(), 3.0));
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumFillet, 1);
    EXPECT_EQ(r.mNumVariableFillet, 0);
    EXPECT_EQ(r.mNumBlend, 0);
    EXPECT_GE(r.mNumCylinder, 1);  // the boss wall (the short rim may fit a sphere)
    // OpenCASCADE 7.9's NURBS conversion throws on the trimmed cylinder: the
    // side channel is then approximated in each face's own parameters.
    expect_brep_consistent(r);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-2);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    bool any_fillet_type = false;
    const NDArray& types = r.mMesh.CellData(kFitCadTypeName, 0);
    for (std::size_t i = 0; i < types.Size(); ++i)
        any_fillet_type = any_fillet_type || types.As<std::int64_t>()[i] == kFitCadFilletType;
    EXPECT_TRUE(any_fillet_type);
}

TEST(FitCad, FilletsOffLeaveTheFilletToTheFreeformFit) {
    FitCadOptions o;
    o.mFillets = false;
    const FitCadResult r = fit_cad(fitted(analytic::inclined_boss(), 3.0), "", o);
    EXPECT_EQ(r.mNumFillet, 0);
    EXPECT_TRUE(r.mValid);
}

TEST(FitCad, ThreadFlanksAreHelicalFaces) {
    // A trapezoidal thread: root and crest are cylinder strips, both flanks
    // ruled helicoids about the same axis (roadmap 6.3.1).
    const FitCadResult r = fit_cad(fitted(analytic::threaded_rod()));
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumHelical, 2);
    EXPECT_EQ(r.mNumRoundedThread, 0);
    EXPECT_EQ(r.mNumCylinder, 2);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-2);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    // The root and crest are cylinders wound several turns round.
    expect_brep_consistent(r);
}

TEST(FitCad, ThreadsOffLeaveTheFlanksToTheFreeformFit) {
    FitCadOptions o;
    o.mThreads = false;
    const FitCadResult r = fit_cad(fitted(analytic::threaded_rod()), "", o);
    EXPECT_EQ(r.mNumHelical, 0);
}

TEST(FitCad, FilletOfVaryingRadiusIsOneFace) {
    // A block's edge rounded by a ball whose radius goes from 0.3 to 0.2 along
    // it: no one radius explains the patch, a radius law along it does
    // (roadmap 6.3.1).
    const Mesh fit = fitted(analytic::variable_fillet_block(), 3.0);
    const FitCadResult r = fit_cad(fit);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumFillet, 1);
    EXPECT_EQ(r.mNumVariableFillet, 1);
    EXPECT_EQ(r.mNumBlend, 0);
    EXPECT_EQ(r.mNumPlane, 6);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-2);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    FitCadOptions o;
    o.mFillets = false;
    const FitCadResult off = fit_cad(fit, "", o);
    EXPECT_EQ(off.mNumFillet, 0);
    EXPECT_EQ(off.mNumVariableFillet, 0);
}

TEST(FitCad, BlendBetweenTwoFreeformFacesIsOneFace) {
    // A bar whose wavy top and leaning front are freeform B-spline faces, the
    // edge between them rounded by a rolling ball: the blend is fitted after
    // the freeform faces, against them (roadmap 6.3.1).
    const Mesh fit = fitted(analytic::wavy_bar(), 3.0);
    const FitCadResult r = fit_cad(fit);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumFillet, 1);
    EXPECT_EQ(r.mNumBlend, 1);
    EXPECT_EQ(r.mNumBspline, 2);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-2);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    FitCadOptions o;
    o.mFillets = false;
    const FitCadResult off = fit_cad(fit, "", o);
    EXPECT_EQ(off.mNumFillet, 0);
    EXPECT_EQ(off.mNumBlend, 0);
}

TEST(FitCad, RoundedThreadWithNoCylinderIsHelical) {
    // Root and crest rounded into arcs and no cylinder anywhere: every face
    // of the thread finds the screw from its own facets, the flanks straight
    // and the root and crest rounded (roadmap 6.3.1).
    const Mesh fit = fitted(analytic::rounded_thread_rod(0.5, 0.8, 1.0, 0.2, 3));
    const FitCadResult r = fit_cad(fit);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumHelical, 4);
    EXPECT_EQ(r.mNumRoundedThread, 2);
    EXPECT_EQ(r.mNumCylinder, 0);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-2);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    expect_brep_consistent(r);
}

// Roadmap 6.3 (v17.2.0): tapered and variable-pitch threads, chamfers, conic
// blends and vertex blends.

TEST(FitCad, TaperedThreadFlanksAreHelicalFaces) {
    // threaded_rod cut on a cone: root and crest are cone strips, and the
    // flanks beside them ruled helicoids reported tapered.
    const FitCadResult r = fit_cad(fitted(analytic::tapered_thread_rod()));
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumHelical, 2);
    EXPECT_EQ(r.mNumTaperedThread, 2);
    EXPECT_EQ(r.mNumVariablePitchThread, 0);
    EXPECT_EQ(r.mNumCone, 2);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
}

TEST(FitCad, ThreadOfVaryingPitchIsHelical) {
    const FitCadResult r = fit_cad(fitted(analytic::variable_pitch_thread_rod()));
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumHelical, 2);
    EXPECT_EQ(r.mNumVariablePitchThread, 2);
    EXPECT_EQ(r.mNumTaperedThread, 0);
    EXPECT_EQ(r.mNumCylinder, 2);
    EXPECT_EQ(r.mNumDemoted, 0);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    FitCadOptions o;
    o.mThreads = false;
    EXPECT_EQ(fit_cad(fitted(analytic::variable_pitch_thread_rod()), "", o).mNumHelical, 0);
}

TEST(FitCad, TaperedRoundedThreadFindsItsScrewAndTaper) {
    // No cylinder or cone anywhere: the rounded root and crest find the
    // screw, its taper and their profile themselves.
    const FitCadResult r = fit_cad(fitted(analytic::tapered_rounded_thread_rod()));
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumHelical, 4);
    EXPECT_EQ(r.mNumRoundedThread, 2);
    EXPECT_EQ(r.mNumTaperedThread, 2);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
}

TEST(FitCad, ChamferOfVaryingWidthIsOneFace) {
    // A chamfer whose width grows along the edge is a twisted ruled surface,
    // which no plane explains.
    const Mesh skin = analytic::section_blend_block(0.0);
    for (const Mesh& fit : {fitted(skin), fit_primitives(skin).mMesh}) {
        const FitCadResult r = fit_cad(fit);
        EXPECT_TRUE(r.mValid);
        EXPECT_EQ(r.mNumChamfer, 1);
        EXPECT_EQ(r.mNumPlane, 6);
        EXPECT_EQ(r.mNumFacetedFaces, 0);
        EXPECT_EQ(r.mNumFaces, 7);
        EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-3);
        EXPECT_LE(r.mMaxDeviation, r.mTolerance);
        expect_brep_consistent(r);
    }
    FitCadOptions o;
    o.mFillets = false;
    EXPECT_EQ(fit_cad(fitted(skin), "", o).mNumChamfer, 0);
}

TEST(FitCad, ConicBlendIsAFilletFace) {
    // Elliptic cross-sections tangent to both faces, legs of unequal and
    // varying length: neither a rolling ball nor a chamfer.
    const FitCadResult r = fit_cad(fit_primitives(analytic::section_blend_block(0.5)).mMesh);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumFillet, 1);
    EXPECT_EQ(r.mNumConicBlend, 1);
    EXPECT_EQ(r.mNumChamfer, 0);
    EXPECT_EQ(r.mNumPlane, 6);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    expect_brep_consistent(r);
}

TEST(FitCad, RoundedBoxCornersAreVertexBlends) {
    // A rounded box scaled unequally: its edge blends are conic blends and
    // each corner a vertex blend tangent to the three meeting there.
    const FitCadResult r = fit_cad(fit_primitives(analytic::rounded_box()).mMesh);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumPlane, 6);
    EXPECT_EQ(r.mNumConicBlend, 12);
    EXPECT_EQ(r.mNumVertexBlend, 8);
    EXPECT_EQ(r.mNumFaces, 26);
    EXPECT_EQ(r.mNumFacetedFaces, 0);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    expect_brep_consistent(r);
    // Unscaled, the blends are cylinders and the corners spheres.
    const FitCadResult u = fit_cad(fit_primitives(analytic::rounded_box(0.3, 1.0, 1.0, 1.0)).mMesh);
    EXPECT_EQ(u.mNumVertexBlend, 0);
    EXPECT_EQ(u.mNumSphere, 8);
    EXPECT_EQ(u.mNumCylinder, 12);
    FitCadOptions o;
    o.mFillets = false;
    const FitCadResult off = fit_cad(fit_primitives(analytic::rounded_box()).mMesh, "", o);
    EXPECT_EQ(off.mNumVertexBlend, 0);
    EXPECT_EQ(off.mNumConicBlend, 0);
}

TEST(FitCad, TypesMaskFacetsTheRest) {
    FitCadOptions o;
    o.mTypes = 1u << static_cast<unsigned>(FitPrimitiveType::Plane);
    const FitCadResult r = fit_cad(fitted(analytic::capped_cylinder()), "", o);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumPlane, 2);
    EXPECT_EQ(r.mNumCylinder, 0);
    EXPECT_GT(r.mNumFacetedFaces, 0);
    EXPECT_EQ(r.mNumDemoted, 0);  // left out by the mask, not demoted
}

TEST(FitCad, ApproximateEdgesStillClose) {
    FitCadOptions o;
    o.mExactEdges = false;
    const FitCadResult r = fit_cad(fitted(analytic::capped_cylinder()), "", o);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumExactEdges, 0);
    EXPECT_GT(r.mNumApproxEdges, 0);
    EXPECT_EQ(r.mNumCylinder, 1);
}

TEST(FitCad, RunsSegmentAndFitWhenMissing) {
    const FitCadResult r = fit_cad(analytic::box_skin(2));
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumPlane, 6);
}

// The same with OpenCASCADE: the turned part segmented by angles alone keeps
// smooth junctions inside patches; refined, every surface is an analytic face
// of a valid solid (v17.11.0).
TEST(FitCad, RefinementExplainsTangentJunctions) {
    const Mesh skin = fitted(analytic::turned_part());
    const FitCadResult before = fit_cad(skin);
    FitCadOptions o;
    o.mRefineSegmentation = true;
    const FitCadResult after = fit_cad(skin, "", o);
    EXPECT_TRUE(after.mValid);
    EXPECT_TRUE(after.mRefined);
    EXPECT_NEAR(after.mExplainedAreaFraction, 1.0, 1e-12);
    const auto analytic_faces = [](const FitCadResult& rR) {
        return rR.mNumPlane + rR.mNumSphere + rR.mNumCylinder + rR.mNumCone + rR.mNumTorus;
    };
    EXPECT_GT(analytic_faces(after), analytic_faces(before));
    EXPECT_EQ(after.mNumFacetedFaces, 0);
    EXPECT_EQ(after.mNumDemoted, 0);
    EXPECT_NEAR(after.mVolumeRatio, 1.0, 1e-2);
    EXPECT_LE(after.mMaxDeviation, after.mTolerance);
}

// The general freeform face (v17.11.0): a wavy top the grid fit cannot reach
// with at most four grid points a side comes back faceted; with the general
// fit it is one B-spline face of a valid solid.
TEST(FitCad, GeneralFreeformFitsWhatTheGridCannot) {
    const Mesh skin = fitted(fcfix::wavy_top_box(24));
    FitCadOptions o;
    o.mFreeformMaxControl = 4;
    const FitCadResult before = fit_cad(skin, "", o);
    EXPECT_EQ(before.mNumBspline, 0);
    EXPECT_GT(before.mNumFacetedFaces, 0);
    EXPECT_EQ(before.mNumFreeformGeneral, 0);
    o.mFreeformGeneral = true;
    const FitCadResult after = fit_cad(skin, "", o);
    EXPECT_TRUE(after.mValid);
    EXPECT_EQ(after.mNumFreeformGeneral, 1);
    EXPECT_EQ(after.mNumBspline, 1);
    EXPECT_EQ(after.mNumPlane, 5);
    EXPECT_EQ(after.mNumFacetedFaces, 0);
    EXPECT_GT(after.mAvgFreeformPoles, 0.0);
    EXPECT_NEAR(after.mVolumeRatio, 1.0, 1e-3);
    EXPECT_LE(after.mMaxDeviation, after.mTolerance);
    expect_brep_consistent(after);
    // The general fit runs only where the grid fit fails: with the default
    // grid the same top is the grid's face, as before.
    FitCadOptions d;
    d.mFreeformGeneral = true;
    const FitCadResult dflt = fit_cad(skin, "", d);
    EXPECT_EQ(dflt.mNumFreeformGeneral, 0);
    EXPECT_EQ(dflt.mNumBspline, 1);
}

// Never silently wrong: with the grid fit starved, the freeform fixtures go
// through the general path (the closed band through its halving too), and
// every solid stays valid and within the tolerance; a face that fails comes
// back as facets.
TEST(FitCad, GeneralFreeformKeepsEverySolidValid) {
    FitCadOptions o;
    o.mFreeformMaxControl = 4;
    o.mFreeformGeneral = true;
    for (const fcfix::Fixture& f : fcfix::golden_fixtures()) {
        if (f.mName != "saddle_box" && f.mName != "wavy_capped")
            continue;  // the freeform ones; the rest have no freeform patch
        SCOPED_TRACE(f.mName);
        const FitCadResult r = fit_cad(f.mMesh, "", o);
        EXPECT_TRUE(r.mValid);
        EXPECT_LE(r.mMaxDeviation, r.mTolerance);
        EXPECT_NEAR(r.mVolumeRatio, 1.0, 2e-2);
        if (f.mName == "wavy_capped") {
            // Starved, the grid fit leaves the wavy wall's halves faceted; the
            // general fit makes them two B-spline faces.
            EXPECT_EQ(r.mNumFreeformGeneral, 2);
            EXPECT_EQ(r.mNumFacetedFaces, 0);
        }
    }
}

TEST(FitCad, WritesReproducibleStep) {
    const std::string a = ::testing::TempDir() + "fit_cad_a.step";
    const Mesh fit = fitted(analytic::capped_cylinder());
    fit_cad(fit, a);
    const std::string bytes = read_bytes(a);
    fit_cad(fit, a);
    EXPECT_NE(bytes.find("ISO-10303-21"), std::string::npos);
    EXPECT_NE(bytes.find("CYLINDRICAL_SURFACE"), std::string::npos);
    EXPECT_NE(bytes.find("1970-01-01T00:00:00"), std::string::npos);
    EXPECT_NE(bytes.find("meshio++ fit_cad"), std::string::npos);
    EXPECT_EQ(bytes, read_bytes(a));
    std::remove(a.c_str());
}

TEST(FitCad, PipelineStepReplacesTheMesh) {
    const std::string path = ::testing::TempDir() + "fit_cad_pipeline.step";
    PipelineStep fit;
    fit.mOp = "FitCad";
    fit.mParams["Path"] = std::string(path);
    fit.mParams["Types"] = std::vector<std::string>{"plane"};
    PipelineReport report;
    const Mesh out = run_pipeline_steps(analytic::box_skin(2), {fit}, report);
    EXPECT_TRUE(out.HasCellData(kFitCadFaceName));
    EXPECT_FALSE(read_bytes(path).empty());
    std::remove(path.c_str());
    fit.mParams["Types"] = std::vector<std::string>{"blob"};
    EXPECT_THROW(run_pipeline_steps(analytic::box_skin(2), {fit}, report), std::invalid_argument);
}

TEST(FitCad, NoisyTurnedPartStaysWithinTolerance) {
    // Every point moved by up to 1e-3, as in doc/cad_fit.md.
    Mesh skin = analytic::turned_part();
    const NDArray& p = skin.Points();
    NDArray q = NDArray::Uninit(DType::Float64, {skin.NumPoints(), 3});
    std::uint64_t state = 7;
    for (std::size_t i = 0; i < skin.NumPoints() * 3; ++i) {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        const double u = static_cast<double>(state >> 11) / 9007199254740992.0;
        q.As<double>()[i] = p.As<double>()[i] + 1e-3 * (2.0 * u - 1.0);
    }
    skin.AssignPoints(std::move(q));
    const FitCadResult r = fit_cad(fitted(skin, 1.0, 2e-3, 20));
    EXPECT_TRUE(r.mValid);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 0.01);
    // Every fitted patch is kept; the facets are a sliver no primitive explains.
    EXPECT_EQ(r.mNumDemoted, 0);
    EXPECT_EQ(r.mNumPlane + r.mNumCylinder + r.mNumCone + r.mNumSphere + r.mNumTorus, 6);
}

TEST(FitCad, RefusesOpenSkinsAndBadOptions) {
    Mesh open = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}});
    EXPECT_THROW(fit_cad(open), std::invalid_argument);
    FitCadOptions o;
    o.mMaxDemotionRounds = -1;
    EXPECT_THROW(fit_cad(analytic::box_skin(2), "", o), std::invalid_argument);
    o = FitCadOptions{};
    o.mTypes = 1u;
    EXPECT_THROW(fit_cad(analytic::box_skin(2), "", o), std::invalid_argument);
    o = FitCadOptions{};
    o.mTolerance = std::nan("");
    EXPECT_THROW(fit_cad(analytic::box_skin(2), "", o), std::invalid_argument);
    o = FitCadOptions{};
    o.mFreeformDegree = 2;
    EXPECT_THROW(fit_cad(analytic::box_skin(2), "", o), std::invalid_argument);
    o = FitCadOptions{};
    o.mFreeformMaxControl = 3;
    EXPECT_THROW(fit_cad(analytic::box_skin(2), "", o), std::invalid_argument);
}

TEST(FitCad, RefinementVerdictOverridesALooserPrimitive) {
    // refine_patches finds no primitive for the saddle top within the
    // tolerance, while fit_primitives, by its own measure, accepts a torus's
    // inner side for it; that torus's edges fail and its demotion took the
    // box's planes with it. fit_cad takes the refinement's verdict.
    FitCadOptions o;
    o.mRefineSegmentation = true;
    for (const FitCadBackend b : {FitCadBackend::Occt, FitCadBackend::Native}) {
        o.mBackend = b;
        const FitCadResult r = fit_cad(fitted(saddle_box(24)), "", o);
        EXPECT_TRUE(r.mValid);
        EXPECT_TRUE(r.mRefined);
        EXPECT_EQ(r.mNumPlane, 5);
        EXPECT_EQ(r.mNumTorus, 0);
        EXPECT_EQ(r.mNumDemoted, 0);
        if (b == FitCadBackend::Occt) {
            EXPECT_EQ(r.mNumFaces, 6);
            EXPECT_EQ(r.mNumBspline, 1);
        }
    }
}

TEST(FitCad, SewToleranceAtTheToleranceIsTheDefault) {
    const Mesh skin = fitted(analytic::capped_cylinder());
    const std::string a = ::testing::TempDir() + "fit_cad_sew.step";
    const FitCadResult r = fit_cad(skin, a);
    const std::string bytes = read_bytes(a);
    FitCadOptions o;
    o.mSewTolerance = r.mTolerance;
    fit_cad(skin, a, o);
    EXPECT_EQ(bytes, read_bytes(a));
    // A tighter sewing still closes a solid whose edges are shared already.
    o.mSewTolerance = 1e-3 * r.mTolerance;
    const FitCadResult t = fit_cad(skin, a, o);
    EXPECT_TRUE(t.mValid);
    EXPECT_TRUE(t.mClosed);
    std::remove(a.c_str());
}

TEST(FitCad, ClosedSolidReported) {
    const FitCadResult r = fit_cad(fitted(analytic::box_skin(3)));
    EXPECT_TRUE(r.mClosed);
    EXPECT_TRUE(std::isnan(r.mHausdorffMax));  // not measured unless asked
    EXPECT_FALSE(r.mMesh.HasFieldData(kFitCadFaceHausdorffName));
}

TEST(FitCad, JunctionSnapPlacesTheCorners) {
    FitCadOptions o;
    o.mSnapJunctions = true;
    const FitCadResult r = fit_cad(fitted(analytic::box_skin(3)), "", o);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mNumSnappedJunctions, 8);  // the box's corners, three planes each
    EXPECT_NEAR(r.mVolume, 1.0, 1e-9);
    EXPECT_LT(r.mMaxDeviation, 1e-9);
    const FitCadResult t = fit_cad(fitted(analytic::turned_part()), "", o);
    EXPECT_TRUE(t.mValid);
    EXPECT_LE(t.mMaxDeviation, t.mTolerance);
}

TEST(FitCad, FaceHausdorffBoundsGlobal) {
    FitCadOptions o;
    o.mFaceDeviation = true;
    const FitCadResult r = fit_cad(fitted(analytic::turned_part()), "", o);
    ASSERT_TRUE(r.mMesh.HasFieldData(kFitCadFaceHausdorffName));
    ASSERT_TRUE(r.mMesh.HasFieldData(kFitCadFaceRmsName));
    const NDArray& h = r.mMesh.FieldData(kFitCadFaceHausdorffName);
    const NDArray& rms = r.mMesh.FieldData(kFitCadFaceRmsName);
    ASSERT_EQ(h.Shape().size(), 2u);
    EXPECT_EQ(h.Shape()[0], static_cast<std::size_t>(r.mNumFaces));
    EXPECT_EQ(h.Shape()[1], 2u);
    EXPECT_EQ(rms.Size(), static_cast<std::size_t>(r.mNumFaces));
    double worst = 0.0;
    for (std::size_t i = 0; i < h.Size(); ++i)
        worst = std::max(worst, h.As<double>()[i]);
    EXPECT_EQ(worst, r.mHausdorffMax);
    EXPECT_LE(r.mHausdorffRms, r.mHausdorffMax);
    for (std::size_t f = 0; f < rms.Size(); ++f)
        EXPECT_LE(rms.As<double>()[f], std::max(h.As<double>()[f * 2], h.As<double>()[f * 2 + 1]));
    // Every face stands within the tolerance of the skin, both ways.
    EXPECT_LE(r.mHausdorffMax, r.mTolerance);
    // The measure does not touch the solid.
    const FitCadResult plain = fit_cad(fitted(analytic::turned_part()));
    EXPECT_EQ(plain.mNumFaces, r.mNumFaces);
    EXPECT_EQ(plain.mMaxDeviation, r.mMaxDeviation);
}

TEST(FitCad, ExtendedSurfacesLeaveFoundEdgesAlone) {
    // Where the surfaces' own intersection has a branch along every chain, the
    // extensions are never intersected: the solid is the one built without.
    FitCadOptions o;
    o.mFreeformMaxControl = 4;
    o.mFreeformGeneral = true;
    for (const fcfix::Fixture& f : fcfix::golden_fixtures()) {
        if (f.mName != "saddle_box" && f.mName != "wavy_capped")
            continue;
        SCOPED_TRACE(f.mName);
        o.mExtendSurfaces = false;
        const FitCadResult plain = fit_cad(f.mMesh, "", o);
        o.mExtendSurfaces = true;
        const FitCadResult r = fit_cad(f.mMesh, "", o);
        EXPECT_TRUE(r.mValid);
        EXPECT_TRUE(r.mClosed);
        EXPECT_EQ(r.mNumExtendedEdges, 0);
        EXPECT_EQ(r.mNumFaces, plain.mNumFaces);
        EXPECT_EQ(r.mNumExactEdges, plain.mNumExactEdges);
        EXPECT_EQ(r.mVolume, plain.mVolume);
    }
}

#endif
