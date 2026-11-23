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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// Small versions of the native fit_cad tests (`test_fit_cad*.cpp`). The full suites fit the
// freeform fixtures at 24 divisions, the closed cages and the tori with the coupled G1 solve,
// which takes minutes each in the instrumented coverage build, so a pull request skips them
// there; these fit the same families at 10 to 16 divisions in every freeform mode (adaptive,
// network, THB refinement, resegmentation, the coupled G1 solve, the kinematic surfaces) and
// check the same properties: a valid fit, the facet area that is left, the written STEP
// surfaces, determinism. Their suite name starts with `Quick`, which the job's exclusion does
// not match. The OpenCASCADE routes stay with the full tests.

// System includes
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "analytic_surfaces.hpp"
#include "fit_cad_fixtures.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"

using namespace meshioplusplus;

namespace {

std::string quick_slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string quick_step(const char* pName) {
    return (std::filesystem::temp_directory_path() / pName).string();
}

// A skin segmented and fitted at the tolerance, and the options of the native backend at it.
struct QuickCase {
    Mesh mFitted;
    FitCadOptions mOptions;
    double mTolerance = 0.0;
};

QuickCase quick_case(const Mesh& rSkin, FitCadFreeformMode Mode) {
    QuickCase c;
    c.mTolerance = 2e-3 * fcfix::diagonal(rSkin);
    c.mFitted = fcfix::fitted_abs(rSkin, c.mTolerance);
    c.mOptions.mBackend = FitCadBackend::Native;
    c.mOptions.mTolerance = c.mTolerance;
    c.mOptions.mFreeformMode = Mode;
    return c;
}

}  // namespace

// The native backend records every fitted patch of a box, writes the per-patch table outside
// the legacy mode, faces everything when no primitive is asked, and demotes a patch past the
// tolerance with a reason.
TEST(QuickFitCad, TheNativeBackendRecordsFacesAndDemotions) {
    const Mesh skin = fcfix::fitted(analytic::box_skin(2));
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    const FitCadResult legacy = fit_cad(skin, "", o);
    EXPECT_FALSE(legacy.mBrepChecked);
    EXPECT_EQ(legacy.mFreeformMode, FitCadFreeformMode::Legacy);
    ASSERT_EQ(legacy.mFaces.size(), 6u);
    for (const FitCadFaceRecord& f : legacy.mFaces) {
        EXPECT_EQ(f.mKind, 1);
        EXPECT_EQ(f.mDemotionReason, FitCadDemotionReason::None);
    }
    EXPECT_FALSE(legacy.mMesh.HasFieldData(kFitCadFaceTableName));
    o.mFreeformMode = FitCadFreeformMode::Adaptive;
    const FitCadResult adaptive = fit_cad(skin, "", o);
    ASSERT_TRUE(adaptive.mMesh.HasFieldData(kFitCadFaceTableName));
    EXPECT_EQ(adaptive.mMesh.FieldData(kFitCadFaceTableName).Shape()[0], 6u);

    FitCadOptions facets;
    facets.mBackend = FitCadBackend::Native;
    facets.mTypes = 0;
    const FitCadResult all = fit_cad(skin, "", facets);
    EXPECT_DOUBLE_EQ(all.mFacetAreaFraction, 1.0);
    EXPECT_TRUE(all.mFaces.empty());

    FitCadOptions tight;
    tight.mBackend = FitCadBackend::Native;
    tight.mTolerance = 1e-6;
    const FitCadResult d =
        fit_cad(fcfix::fitted(fcfix::jittered(analytic::capped_cylinder(1.0, 2.0, 3, 16), 1e-4)),
                "", tight);
    EXPECT_GT(d.mCounts.mNumDemoted, 0);
    bool found = false;
    for (const FitCadFaceRecord& f : d.mFaces)
        found = found || f.mDemotionReason == FitCadDemotionReason::HausdorffExceeded;
    EXPECT_TRUE(found);
}

// The adaptive mode turns the freeform top of a box into one B-spline face, written as
// B_SPLINE_SURFACE_WITH_KNOTS, where the legacy mode leaves facets; THB refinement fits the
// bump with no more facet area and the same faces twice.
TEST(QuickFitCad, AdaptiveFacesAreBsplinesAndThbRefinesLocally) {
    QuickCase c = quick_case(analytic::gaussian_bump_box(8), FitCadFreeformMode::Legacy);
    const FitCadResult legacy = fit_cad(c.mFitted, "", c.mOptions);
    c.mOptions.mFreeformMode = FitCadFreeformMode::Adaptive;
    const std::string path = quick_step("quick_adaptive.step");
    const FitCadResult adaptive = fit_cad(c.mFitted, path, c.mOptions);
    EXPECT_GE(adaptive.mNumBspline, 1);
    EXPECT_LT(adaptive.mFacetAreaFraction, legacy.mFacetAreaFraction);
    EXPECT_TRUE(adaptive.mValid);
    EXPECT_LE(adaptive.mMaxFreeformDistance, c.mTolerance);
    EXPECT_NE(quick_slurp(path).find("B_SPLINE_SURFACE_WITH_KNOTS"), std::string::npos);
    c.mOptions.mLocalRefine = FitCadLocalRefine::Thb;
    const FitCadResult thb = fit_cad(c.mFitted, "", c.mOptions);
    EXPECT_TRUE(thb.mValid);
    EXPECT_GT(thb.mNumFreeformAdaptive, 0);
    const FitCadResult again = fit_cad(c.mFitted, "", c.mOptions);
    EXPECT_EQ(again.mFacetAreaFraction, thb.mFacetAreaFraction);
    ASSERT_EQ(again.mFaces.size(), thb.mFaces.size());
}

// A side merged over a missed crease is cut again by the resegmentation.
TEST(QuickFitCad, ResegmentationCutsAMissedCrease) {
    const Mesh skin = analytic::hood_box(8);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    const Mesh merged = fcfix::merged_side(fcfix::fitted_abs(skin, tol));
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTolerance = tol;
    o.mFreeformMode = FitCadFreeformMode::Adaptive;
    o.mResegmentRounds = 0;
    const FitCadResult off = fit_cad(merged, "", o);
    o.mResegmentRounds = 2;
    const FitCadResult on = fit_cad(merged, "", o);
    EXPECT_EQ(off.mNumResegmented, 0);
    EXPECT_GT(on.mNumResegmented, 0);
    EXPECT_LE(on.mFacetAreaFraction, off.mFacetAreaFraction);
    EXPECT_TRUE(on.mValid);
}

// The network mode writes B-spline patches that share their edges, G1 across them; the coupled
// solve is no worse; the same bytes again.
TEST(QuickFitCad, TheNetworkSharesEdgesAndTheCoupledSolveIsNoWorse) {
    QuickCase c = quick_case(analytic::monkey_saddle_box(16), FitCadFreeformMode::Network);
    const std::string path = quick_step("quick_network.step");
    const FitCadResult sw = fit_cad(c.mFitted, path, c.mOptions);
    EXPECT_GE(sw.mNumNetworkPatches, 2);
    EXPECT_GT(sw.mNumNetworkEdges, 0);
    EXPECT_TRUE(sw.mValid);
    EXPECT_NE(quick_slurp(path).find("B_SPLINE_CURVE_WITH_KNOTS('',3,"), std::string::npos);
    c.mOptions.mFreeformG1Solve = FitCadNetworkG1Solve::Coupled;
    c.mOptions.mValidationReport = true;
    const FitCadResult cp = fit_cad(c.mFitted, quick_step("quick_coupled.step"), c.mOptions);
    EXPECT_TRUE(cp.mValid);
    EXPECT_LE(cp.mFacetAreaFraction, sw.mFacetAreaFraction);
    EXPECT_GE(cp.mNumCoupledRegions, 1);
    EXPECT_LE(cp.mMaxG1AngleAwayDeg, sw.mMaxG1AngleAwayDeg + 1e-12);
    EXPECT_EQ(sw.mNumCoupledRegions, 0);
    const FitCadResult again = fit_cad(c.mFitted, "", c.mOptions);
    EXPECT_EQ(again.mMaxG1AngleDeg, cp.mMaxG1AngleDeg);
}

// A closed skin of genus one is fitted in the network mode, the whole of it by faces or facets.
TEST(QuickFitCad, AClosedTorusIsFittedByTheNetwork) {
    QuickCase c = quick_case(analytic::organic_torus(24, 12), FitCadFreeformMode::Network);
    const FitCadResult r = fit_cad(c.mFitted, "", c.mOptions);
    EXPECT_TRUE(r.mValid);
    EXPECT_GE(r.mFacetAreaFraction, 0.0);
    EXPECT_LE(r.mFacetAreaFraction, 1.0);
    EXPECT_FALSE(r.mFaces.empty());
}

// The kinematic surfaces: an extruded lobe is one extrusion, a wavy cap one revolution,
// written as the exact STEP surfaces, within the tolerance; off by default.
TEST(QuickFitCad, KinematicSurfacesAreExtrusionsAndRevolutions) {
    struct Case {
        const char* mName;
        Mesh (*pMake)();
        std::int64_t FitCadResult::*pCount;
        const char* mStep;
    };
    for (Case& k : std::vector<Case>{{"extruded_lobe", [] { return analytic::extruded_lobe(); },
                                      &FitCadResult::mNumExtrusion, "SURFACE_OF_LINEAR_EXTRUSION"},
                                     {"wavy_capped", [] { return analytic::wavy_capped(); },
                                      &FitCadResult::mNumRevolution, "SURFACE_OF_REVOLUTION"}}) {
        SCOPED_TRACE(k.mName);
        QuickCase c = quick_case(k.pMake(), FitCadFreeformMode::Legacy);
        const FitCadResult off = fit_cad(c.mFitted, "", c.mOptions);
        EXPECT_EQ(off.*k.pCount, 0);
        c.mOptions.mKinematicSurfaces = true;
        const std::string path = quick_step("quick_kinematic.step");
        const FitCadResult r = fit_cad(c.mFitted, path, c.mOptions);
        EXPECT_GT(r.*k.pCount, 0);
        EXPECT_TRUE(r.mValid);
        EXPECT_NE(quick_slurp(path).find(k.mStep), std::string::npos);
    }
}

// Options out of range are refused naming the option.
TEST(QuickFitCad, OptionsAreValidated) {
    const Mesh skin = fcfix::fitted(analytic::box_skin(2));
    const auto refused = [&](auto Set, const char* pName) {
        FitCadOptions o;
        o.mBackend = FitCadBackend::Native;
        Set(o);
        try {
            fit_cad(skin, "", o);
            ADD_FAILURE() << pName << " was accepted";
        } catch (const std::invalid_argument& rErr) {
            EXPECT_NE(std::string(rErr.what()).find(pName), std::string::npos) << rErr.what();
        }
    };
    refused([](FitCadOptions& o) { o.mFreeformMode = static_cast<FitCadFreeformMode>(4); },
            "freeform_mode");
    refused([](FitCadOptions& o) { o.mLocalRefine = static_cast<FitCadLocalRefine>(2); },
            "local_refine");
    refused([](FitCadOptions& o) { o.mFreeformG1ToleranceDeg = 0.0; }, "freeform_g1_tolerance_deg");
    refused([](FitCadOptions& o) { o.mResegmentRounds = -1; }, "resegment_rounds");
    refused([](FitCadOptions& o) { o.mFreeformG1Solve = static_cast<FitCadNetworkG1Solve>(2); },
            "freeform_g1_solve");
}
