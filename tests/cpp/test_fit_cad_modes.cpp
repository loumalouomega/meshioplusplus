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
// fit_cad's freeform modes (v17.12.0): the option surface, the per-patch
// records and demotion reasons, the facet-area measure, and the facet-area
// table of the freeform fixtures.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

#include "meshioplusplus/operations/fit_cad.hpp"
#include "analytic_surfaces.hpp"
#include "fit_cad_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// Whether every point of a skin's facets is used by an edge two triangles share
// in opposite senses: the closed, consistently wound input fit_cad needs.
bool closed_and_wound(const Mesh& rSkin) {
    std::vector<std::pair<std::int64_t, std::int64_t>> half;
    for (std::size_t b = 0; b < rSkin.NumCellBlocks(); ++b) {
        const NDArray& c = rSkin.Cells(b).Conn();
        const std::size_t k = rSkin.Cells(b).NodesPerCell();
        for (std::size_t i = 0; i < rSkin.Cells(b).NumCells(); ++i)
            for (std::size_t j = 0; j < k; ++j)
                half.emplace_back(detail::read_int(c, i * k + j),
                                  detail::read_int(c, i * k + (j + 1) % k));
    }
    std::sort(half.begin(), half.end());
    std::size_t degenerate = 0, unpaired = 0;
    for (const auto& e : half) {
        degenerate += e.first == e.second;
        unpaired +=
            !std::binary_search(half.begin(), half.end(), std::make_pair(e.second, e.first));
    }
    const bool repeated = std::adjacent_find(half.begin(), half.end()) != half.end();
    if (degenerate || unpaired || repeated)
        std::printf("closed_and_wound: %zu degenerate, %zu unpaired half-edges, repeated %d\n",
                    degenerate, unpaired, repeated ? 1 : 0);
    return !degenerate && !unpaired && !repeated;
}

bool kpi_full() {
    const char* p = std::getenv("MESHIOPLUSPLUS_FIT_CAD_KPI");
    return p && std::string(p) == "1";
}

}  // namespace

// Every new option is refused out of range, naming it, in every build.
TEST(FitCadModes, OptionsAreValidated) {
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
    refused(
        [](FitCadOptions& o) {
            o.mFreeformParameterization = static_cast<FitCadParameterization>(-1);
        },
        "freeform_parameterization");
    refused(
        [](FitCadOptions& o) {
            o.mExtraordinaryVertex = static_cast<FitCadExtraordinaryVertex>(7);
        },
        "extraordinary_vertex");
    refused([](FitCadOptions& o) { o.mLocalRefine = static_cast<FitCadLocalRefine>(2); },
            "local_refine");
    refused([](FitCadOptions& o) { o.mFreeformG1ToleranceDeg = 0.0; }, "freeform_g1_tolerance_deg");
    refused([](FitCadOptions& o) { o.mFreeformMaxControlPoints = 15; },
            "freeform_max_control_points");
    refused([](FitCadOptions& o) { o.mLayoutTargetPatchSize = -1.0; }, "layout_target_patch_size");
    refused([](FitCadOptions& o) { o.mLayoutMaxValenceDeviation = 4; },
            "layout_max_valence_deviation");
    refused([](FitCadOptions& o) { o.mResegmentRounds = -1; }, "resegment_rounds");
}

// The native backend records every fitted patch, measures the facet area and
// never claims a kernel check; the per-patch table appears only outside the
// legacy mode, so a default run's outputs are what they were.
TEST(FitCadModes, NativeRecordsAndFacetArea) {
    const Mesh skin = fcfix::fitted(analytic::box_skin(3));
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    const FitCadResult legacy = fit_cad(skin, "", o);
    EXPECT_FALSE(legacy.mBrepChecked);
    EXPECT_EQ(legacy.mFreeformMode, FitCadFreeformMode::Legacy);
    EXPECT_DOUBLE_EQ(legacy.mFacetAreaFraction, 0.0);
    ASSERT_EQ(legacy.mFaces.size(), 6u);
    for (const FitCadFaceRecord& f : legacy.mFaces) {
        EXPECT_EQ(f.mKind, 1);
        EXPECT_GE(f.mFace, 0);
        EXPECT_EQ(f.mDemotionReason, FitCadDemotionReason::None);
    }
    EXPECT_FALSE(legacy.mMesh.HasFieldData(kFitCadFaceTableName));

    o.mFreeformMode = FitCadFreeformMode::Adaptive;
    const FitCadResult adaptive = fit_cad(skin, "", o);
    EXPECT_EQ(adaptive.mFreeformMode, FitCadFreeformMode::Adaptive);
    ASSERT_TRUE(adaptive.mMesh.HasFieldData(kFitCadFaceTableName));
    const NDArray& table = adaptive.mMesh.FieldData(kFitCadFaceTableName);
    ASSERT_EQ(table.Shape().size(), 2u);
    EXPECT_EQ(table.Shape()[0], 6u);
    EXPECT_EQ(table.Shape()[1], 13u);
}

// Everything faceted: the whole area is on facet faces.
TEST(FitCadModes, NativeAllFacetsIsTheWholeArea) {
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTypes = 0;
    const FitCadResult r = fit_cad(fcfix::fitted(analytic::box_skin(3)), "", o);
    EXPECT_DOUBLE_EQ(r.mFacetAreaFraction, 1.0);
    EXPECT_TRUE(r.mFaces.empty());
}

// A patch whose points lie beyond the tolerance of its surface is written as
// facets, with the reason recorded.
TEST(FitCadModes, NativeDemotionHasAReason) {
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTolerance = 1e-6;
    const FitCadResult r =
        fit_cad(fcfix::fitted(fcfix::jittered(analytic::capped_cylinder(), 1e-4)), "", o);
    ASSERT_GT(r.mNumDemoted, 0);
    bool found = false;
    for (const FitCadFaceRecord& f : r.mFaces)
        if (f.mDemotionReason == FitCadDemotionReason::HausdorffExceeded) {
            EXPECT_EQ(f.mFace, -1);
            EXPECT_EQ(f.mKind, 0);
            found = true;
        }
    EXPECT_TRUE(found);
    EXPECT_GT(r.mFacetAreaFraction, 0.0);
}

// The freeform fixtures are closed, consistently wound skins.
TEST(FitCadModes, FreeformFixturesAreClosedSkins) {
    for (const Mesh& m :
         {analytic::monkey_saddle_box(8), analytic::gaussian_bump_box(8), analytic::s_band_box(8),
          analytic::hood_box(8, true), analytic::organic_torus(16, 8),
          analytic::sphere_with_handle(16), analytic::sphere_with_handle(24),
          analytic::sphere_with_handle(40), analytic::double_torus(20), analytic::double_torus(40),
          analytic::cc_cage(3, 2), analytic::cc_cage(5, 2), analytic::cc_cage(6, 2),
          analytic::g1_loaf(8, 6)})
        EXPECT_TRUE(closed_and_wound(m)) << m.NumPoints() << " points";
}

// The same noise on every platform: the Irwin-Hall draw is integer-built.
TEST(FitCadModes, NoiseIsDeterministic) {
    const Mesh a = analytic::gaussian_noise(analytic::hood_box(6), 1e-3, 11);
    const Mesh b = analytic::gaussian_noise(analytic::hood_box(6), 1e-3, 11);
    ASSERT_EQ(a.NumPoints(), b.NumPoints());
    for (std::size_t i = 0; i < a.NumPoints() * 3; ++i)
        ASSERT_EQ(a.Points().As<double>()[i], b.Points().As<double>()[i]);
}

// The trig-free freeform fixtures' points, bit for bit the NumPy-free Python
// ones (tests/python/test_fit_cad_modes.py pins the same digests): FNV-1a of
// the little-endian float64 coordinates.
TEST(FitCadModes, FixturesMatchPython) {
    const auto digest = [](const Mesh& rM) {
        std::uint64_t h = 0xcbf29ce484222325ull;
        const auto* p = reinterpret_cast<const unsigned char*>(rM.Points().As<double>());
        for (std::size_t i = 0; i < rM.NumPoints() * 3 * sizeof(double); ++i)
            h = (h ^ p[i]) * 0x100000001b3ull;
        return h;
    };
    EXPECT_EQ(digest(analytic::monkey_saddle_box(24)), 0x1b4f3a52c6d6f26cull);
    EXPECT_EQ(digest(analytic::s_band_box(24)), 0xe77136aca0de5eafull);
    EXPECT_EQ(digest(analytic::hood_box(24, true)), 0x4b1fe1e9066260d6ull);
    EXPECT_EQ(digest(analytic::cc_cage(3)), 0x85f0b2dd4c5350e1ull);
    EXPECT_EQ(digest(analytic::gaussian_noise(analytic::hood_box(12), 1e-3, 11)),
              0x6be481cc0438d8c3ull);
}

namespace {

std::string slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string temp_step(const char* pName) {
    return (std::filesystem::temp_directory_path() / pName).string();
}

}  // namespace

// The adaptive mode in the native backend (v17.12.0): the freeform top of a
// box becomes one B-spline face, written as B_SPLINE_SURFACE_WITH_KNOTS, where
// the legacy native backend writes it as facets.
TEST(FitCadModes, NativeAdaptiveBuildsBsplineFaces) {
    int index = 0;
    for (const Mesh& skin : {analytic::hood_box(16), analytic::monkey_saddle_box(16),
                             analytic::s_band_box(16), analytic::gaussian_bump_box(16)}) {
        SCOPED_TRACE(index++);
        const double tol = 2e-3 * fcfix::diagonal(skin);
        const Mesh fitted = fcfix::fitted_abs(skin, tol);
        FitCadOptions o;
        o.mBackend = FitCadBackend::Native;
        o.mTolerance = tol;
        const FitCadResult legacy = fit_cad(fitted, "", o);
        o.mFreeformMode = FitCadFreeformMode::Adaptive;
        const std::string path = temp_step("adaptive_native.step");
        const FitCadResult adaptive = fit_cad(fitted, path, o);
        EXPECT_GT(legacy.mFacetAreaFraction, 0.0);
        EXPECT_GE(adaptive.mNumBspline, 1);
        EXPECT_EQ(adaptive.mNumFreeformAdaptive, adaptive.mNumBspline);
        EXPECT_LT(adaptive.mFacetAreaFraction, legacy.mFacetAreaFraction);
        EXPECT_TRUE(adaptive.mValid);
        EXPECT_FALSE(adaptive.mBrepChecked);
        EXPECT_LE(adaptive.mMaxFreeformDistance, tol);
        const std::string bytes = slurp(path);
        EXPECT_NE(bytes.find("B_SPLINE_SURFACE_WITH_KNOTS"), std::string::npos);
        bool bspline = false;
        for (const FitCadFaceRecord& f : adaptive.mFaces)
            if (f.mKind == kFitCadBsplineType) {
                bspline = true;
                EXPECT_NE(f.mParamMethod, FitCadParamMethod::None);
                EXPECT_LE(f.mNumControl, o.mFreeformMaxControlPoints);
                EXPECT_LE(f.mMaxDist, tol);
            }
        EXPECT_TRUE(bspline);
    }
}

// The native network: B-spline patches whose shared edges are written as
// the cubic curves both patches' boundaries are, and the same bytes again.
TEST(FitCadModes, NativeNetworkWritesSharedEdges) {
    const Mesh skin = analytic::monkey_saddle_box(24);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    const Mesh fitted = fcfix::fitted_abs(skin, tol);
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTolerance = tol;
    o.mFreeformMode = FitCadFreeformMode::Network;
    const std::string path = temp_step("network_native.step");
    const FitCadResult r = fit_cad(fitted, path, o);
    EXPECT_GE(r.mNumNetworkPatches, 4);
    EXPECT_GT(r.mNumNetworkEdges, 0);
    EXPECT_EQ(r.mFacetAreaFraction, 0.0);
    EXPECT_TRUE(r.mValid);
    EXPECT_LE(r.mMaxG1AngleDeg, 5.0);
    const std::string bytes = slurp(path);
    std::size_t cubic = 0;
    for (std::size_t at = bytes.find("B_SPLINE_CURVE_WITH_KNOTS('',3,"); at != std::string::npos;
         at = bytes.find("B_SPLINE_CURVE_WITH_KNOTS('',3,", at + 1))
        ++cubic;
    EXPECT_GT(cubic, 0U);
    const std::string again = temp_step("network_native_again.step");
    fit_cad(fitted, again, o);
    std::string b2 = slurp(again);
    // Only the file name, which the header and the product carry, differs.
    for (std::size_t at = b2.find("network_native_again.step"); at != std::string::npos;
         at = b2.find("network_native_again.step"))
        b2.replace(at, std::string("network_native_again.step").size(), "network_native.step");
    EXPECT_EQ(bytes, b2);
}

// The same bytes at every thread count, and run after run.
TEST(FitCadModes, NativeAdaptiveIsDeterministic) {
    const Mesh skin = analytic::monkey_saddle_box(16, true);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    const Mesh fitted = fcfix::fitted_abs(skin, tol);
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTolerance = tol;
    o.mFreeformMode = FitCadFreeformMode::Adaptive;
    std::string ref;
    for (int threads : {1, 2, 8, 1}) {
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
        omp_set_num_threads(threads);
#endif
        const std::string path = temp_step("adaptive_det.step");
        fit_cad(fitted, path, o);
        const std::string bytes = slurp(path);
        if (ref.empty())
            ref = bytes;
        else
            EXPECT_EQ(bytes, ref) << threads;
    }
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(omp_get_num_procs());
#endif
}

// A noisy skin: the face follows the surface, not the noise -- its poles stay
// within the budget and its distance within the noise allowance.
TEST(FitCadModes, NativeAdaptiveOnNoisySkins) {
    const Mesh clean = analytic::hood_box(20);
    const double tol = 2e-3 * fcfix::diagonal(clean);
    for (double k : {0.1, 0.5}) {
        const Mesh fitted = fcfix::fitted_abs(analytic::gaussian_noise(clean, k * tol), tol, 20);
        FitCadOptions o;
        o.mBackend = FitCadBackend::Native;
        o.mTolerance = tol;
        const FitCadResult legacy = fit_cad(fitted, "", o);
        o.mFreeformMode = FitCadFreeformMode::Adaptive;
        const FitCadResult r = fit_cad(fitted, "", o);
        EXPECT_LT(r.mFacetAreaFraction, legacy.mFacetAreaFraction) << k;
        for (const FitCadFaceRecord& f : r.mFaces)
            if (f.mKind == kFitCadBsplineType) {
                EXPECT_LE(f.mNumControl, o.mFreeformMaxControlPoints) << k;
                EXPECT_LE(f.mMaxDist, std::max(tol, 4.5 * k * tol * 1.6)) << k;
            }
    }
}

// A degenerate patch never crashes: it is recorded with a reason.
TEST(FitCadModes, DegeneratePatchesDemoteWithAReason) {
    // A box whose top is one patch of two triangles: too small to fit.
    Mesh skin = fcfix::fitted_abs(analytic::monkey_saddle_box(1), 2e-3);
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mFreeformMode = FitCadFreeformMode::Auto;
    const FitCadResult r = fit_cad(skin, "", o);
    for (const FitCadFaceRecord& f : r.mFaces)
        if (f.mFace < 0)
            EXPECT_NE(f.mDemotionReason, FitCadDemotionReason::None);
}

#ifdef MESHIOPLUSPLUS_HAS_OCCT

// The adaptive mode through OpenCASCADE: checked B-spline faces, no more facet
// area than the legacy mode, and strictly less where the legacy mode left any.
TEST(FitCadModes, OcctAdaptiveReducesFacetArea) {
    for (const Mesh& skin :
         {analytic::gaussian_noise(analytic::hood_box(16), 2e-3 * 1.7320508 * 0.5),
          analytic::cc_cage(3, 3), analytic::monkey_saddle_box(16, true)}) {
        const double tol = 2e-3 * fcfix::diagonal(skin);
        const Mesh fitted = fcfix::fitted_abs(skin, tol, 20);
        FitCadOptions o;
        o.mTolerance = tol;
        const FitCadResult legacy = fit_cad(fitted, "", o);
        o.mFreeformMode = FitCadFreeformMode::Adaptive;
        const FitCadResult adaptive = fit_cad(fitted, "", o);
        EXPECT_TRUE(adaptive.mBrepChecked);
        EXPECT_LE(adaptive.mFacetAreaFraction, legacy.mFacetAreaFraction + 1e-12);
        if (legacy.mFacetAreaFraction > 0.0)
            EXPECT_LT(adaptive.mFacetAreaFraction, legacy.mFacetAreaFraction);
        EXPECT_TRUE(adaptive.mValid);
    }
}

// The OpenCASCADE backend records each patch it tried as a surface: a freeform
// face with its parametrisation and pole count, and a kernel check.
TEST(FitCadModes, OcctRecordsFreeformFaces) {
    const FitCadResult r = fit_cad(fcfix::fitted(fcfix::saddle_box(12)));
    EXPECT_TRUE(r.mBrepChecked);
    ASSERT_GT(r.mNumBspline, 0);
    bool found = false;
    for (const FitCadFaceRecord& f : r.mFaces)
        if (f.mKind == kFitCadBsplineType) {
            EXPECT_GE(f.mFace, 0);
            EXPECT_TRUE(f.mParamMethod == FitCadParamMethod::Plane ||
                        f.mParamMethod == FitCadParamMethod::Cylinder);
            EXPECT_GT(f.mNumControl, 0);
            found = true;
        }
    EXPECT_TRUE(found);
    EXPECT_GE(r.mFacetAreaFraction, 0.0);
    EXPECT_LE(r.mFacetAreaFraction, 1.0);
    EXPECT_FALSE(r.mMesh.HasFieldData(kFitCadFaceTableName));
}

// A patch every fit fails is recorded as demoted, with a reason.
TEST(FitCadModes, OcctUnfittedPatchHasAReason) {
    FitCadOptions o;
    o.mFreeformMaxControl = 4;  // too coarse a grid for the saddle's top
    o.mFreeformSplit = false;
    const FitCadResult r = fit_cad(fcfix::fitted(fcfix::saddle_box(12)), "", o);
    bool demoted = false;
    for (const FitCadFaceRecord& f : r.mFaces)
        if (f.mFace < 0) {
            EXPECT_NE(f.mDemotionReason, FitCadDemotionReason::None);
            demoted = true;
        }
    if (r.mNumBspline == 0) {
        EXPECT_TRUE(demoted);
        EXPECT_GT(r.mFacetAreaFraction, 0.0);
    }
}

#endif

// The facet-area table of the freeform fixtures, by mode. By default three
// fixtures run; MESHIOPLUSPLUS_FIT_CAD_KPI=1 runs the whole set with its
// variants and prints the table doc/fit_cad.md records.
TEST(FitCadModes, FacetAreaTable) {
    std::vector<fcfix::FreeformFixture> set = fcfix::freeform_fixtures(kpi_full());
    if (!kpi_full()) {
        std::vector<fcfix::FreeformFixture> some;
        for (auto& f : set)
            if (f.mName == "monkey_saddle" || f.mName == "cc_valence3" || f.mName == "hood")
                some.push_back(std::move(f));
        set = std::move(some);
    }
    std::printf("%-24s %10s %10s %10s %10s\n", "fixture", "legacy", "adaptive", "network", "auto");
    for (const fcfix::FreeformFixture& f : set) {
        FitCadOptions o;
        o.mTolerance = f.mTolerance;
        const FitCadResult legacy = fit_cad(f.mMesh, "", o);
        EXPECT_GE(legacy.mFacetAreaFraction, 0.0) << f.mName;
        EXPECT_LE(legacy.mFacetAreaFraction, 1.0) << f.mName;
        double frac[3] = {0.0, 0.0, 0.0};
        int k = 0;
        for (FitCadFreeformMode mode :
             {FitCadFreeformMode::Adaptive, FitCadFreeformMode::Network, FitCadFreeformMode::Auto}) {
            o.mFreeformMode = mode;
            const FitCadResult r = fit_cad(f.mMesh, "", o);
            EXPECT_LE(r.mFacetAreaFraction, legacy.mFacetAreaFraction + 1e-12) << f.mName << " " << k;
            // Every record of a patch that is not a face says why.
            for (const FitCadFaceRecord& rec : r.mFaces)
                if (rec.mFace < 0)
                    EXPECT_NE(rec.mDemotionReason, FitCadDemotionReason::None) << f.mName;
            frac[k++] = r.mFacetAreaFraction;
            if (mode == FitCadFreeformMode::Network)
                std::printf("  (network: %lld patches, %lld edges, %lld EVs, G1 %.2f deg, %lld adaptive)\n",
                            static_cast<long long>(r.mNumNetworkPatches),
                            static_cast<long long>(r.mNumNetworkEdges),
                            static_cast<long long>(r.mNumExtraordinaryVertices), r.mMaxG1AngleDeg,
                            static_cast<long long>(r.mNumFreeformAdaptive));
        }
        std::printf("%-24s %10.4f %10.4f %10.4f %10.4f\n", f.mName.c_str(), legacy.mFacetAreaFraction,
                    frac[0], frac[1], frac[2]);
    }
}
