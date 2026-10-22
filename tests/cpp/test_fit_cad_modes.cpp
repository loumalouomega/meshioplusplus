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

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/primitives.hpp"
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
    EXPECT_EQ(digest(analytic::monkey_saddle_box(24)), 0x06c59d1fc5a203c4ull);
    EXPECT_EQ(digest(analytic::s_band_box(24)), 0x04808850ab51bf6dull);
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
    // Between two network faces, each edge carries a p-curve on both.
    EXPECT_NE(bytes.find("SURFACE_CURVE("), std::string::npos);
    EXPECT_NE(bytes.find("PCURVE("), std::string::npos);
    EXPECT_NE(bytes.find("PARAMETRIC_REPRESENTATION_CONTEXT"), std::string::npos);
    const std::string again = temp_step("network_native_again.step");
    fit_cad(fitted, again, o);
    std::string b2 = slurp(again);
    // Only the file name, which the header and the product carry, differs.
    for (std::size_t at = b2.find("network_native_again.step"); at != std::string::npos;
         at = b2.find("network_native_again.step"))
        b2.replace(at, std::string("network_native_again.step").size(), "network_native.step");
    EXPECT_EQ(bytes, b2);
}

namespace {

// The area of a mesh's polygons (fanned), by blocks of fixed width.
double polygon_area(const Mesh& rMesh) {
    const double* xyz = rMesh.Points().As<double>();
    double area = 0.0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cells = rMesh.Cells(b);
        const NDArray& conn = cells.Conn();
        const std::size_t w = cells.NodesPerCell();
        // Lagrange quadrilaterals: their four corners.
        const std::size_t n = w >= 4 && w != 6 ? 4 : w;
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            const double* p0 = xyz + 3 * detail::read_int(conn, c * w);
            for (std::size_t k = 1; k + 1 < n; ++k) {
                const double* p1 = xyz + 3 * detail::read_int(conn, c * w + k);
                const double* p2 = xyz + 3 * detail::read_int(conn, c * w + k + 1);
                const double u[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
                const double v[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
                const double x[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                     u[0] * v[1] - u[1] * v[0]};
                area += 0.5 * std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
            }
        }
    }
    return area;
}

}  // namespace

// `local_refine="thb"` (experimental): the adaptive faces fitted with
// truncated hierarchical B-splines and flattened, their records counting the
// hierarchical degrees of freedom.
TEST(FitCadModes, NativeThbRefinementFitsTheBump) {
    const Mesh skin = analytic::gaussian_bump_box(24);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    const Mesh fitted = fcfix::fitted_abs(skin, tol);
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTolerance = tol;
    o.mFreeformMode = FitCadFreeformMode::Adaptive;
    const FitCadResult knots = fit_cad(fitted, "", o);
    o.mLocalRefine = FitCadLocalRefine::Thb;
    const FitCadResult thb = fit_cad(fitted, "", o);
    EXPECT_TRUE(thb.mValid);
    EXPECT_GT(thb.mNumFreeformAdaptive, 0);
    EXPECT_LE(thb.mFacetAreaFraction, knots.mFacetAreaFraction + 1e-12);
    for (const FitCadFaceRecord& rec : thb.mFaces)
        if (rec.mKind == kFitCadBsplineType) {
            EXPECT_GT(rec.mNumControl, 0);
            EXPECT_LE(rec.mMaxDist, 0.5 * tol);
        }
    // Deterministic, as every freeform face.
    const FitCadResult again = fit_cad(fitted, "", o);
    EXPECT_EQ(again.mFacetAreaFraction, thb.mFacetAreaFraction);
    ASSERT_EQ(again.mFaces.size(), thb.mFaces.size());
    for (std::size_t i = 0; i < thb.mFaces.size(); ++i)
        EXPECT_EQ(again.mFaces[i].mNumControl, thb.mFaces[i].mNumControl);
}

TEST(FitCadModes, NativeModesReturnTheSolidAsNurbs) {
    std::vector<fcfix::Fixture> set = fcfix::golden_fixtures();
    {
        const Mesh skin = analytic::monkey_saddle_box(24);
        set.push_back(
            {"monkey_saddle_network", fcfix::fitted_abs(skin, 2e-3 * fcfix::diagonal(skin))});
    }
    for (const fcfix::Fixture& f : set) {
        FitCadOptions o;
        o.mBackend = FitCadBackend::Native;
        o.mFreeformMode = f.mName == "monkey_saddle_network" ? FitCadFreeformMode::Network
                                                             : FitCadFreeformMode::Adaptive;
        if (f.mName == "monkey_saddle_network")
            o.mTolerance = 2e-3 * fcfix::diagonal(f.mMesh);
        const FitCadResult r = fit_cad(f.mMesh, "", o);
        ASSERT_EQ(static_cast<std::int64_t>(r.mBrep.mFaces.size()), r.mNumFaces) << f.mName;
        ASSERT_EQ(r.mBrep.mPatches.size(), r.mBrep.mFaces.size()) << f.mName;
        // Every edge between two faces, or a seam of one, is used twice.
        for (const BrepEdge& e : r.mBrep.mEdges)
            EXPECT_TRUE(e.mTopology.size() == 2 || (!e.mHasCurve && e.mTopology.size() == 1))
                << f.mName << " edge " << e.mId;
        // Sampled, the faces cover the skin's area.
        CadJsonReadOptions so;
        so.mSubdivisions = 32;
        so.mEdges = false;
        const double sampled = polygon_area(sample_brep(r.mBrep, so));
        const double skin = polygon_area(f.mMesh);
        EXPECT_NEAR(sampled, skin, 0.03 * skin) << f.mName;
        std::printf("%-24s faces %zu edges %zu area ratio %.5f\n", f.mName.c_str(),
                    r.mBrep.mFaces.size(), r.mBrep.mEdges.size(), sampled / skin);
        // And Kratos CAD JSON round-trips it.
        const std::string path = temp_step((f.mName + ".cad.json").c_str());
        write_cad_json(path, r.mMesh, r.mBrep);
        const BrepInfo back = read_cad_json_info(path);
        EXPECT_EQ(back.mFaces.size(), r.mBrep.mFaces.size()) << f.mName;
        EXPECT_EQ(back.mEdges.size(), r.mBrep.mEdges.size()) << f.mName;
    }
    // The legacy mode leaves it empty, as before.
    FitCadOptions legacy;
    legacy.mBackend = FitCadBackend::Native;
    EXPECT_TRUE(fit_cad(set.front().mMesh, "", legacy).mBrep.mFaces.empty());
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

// The re-segmentation rung: a hood whose segmentation missed the creases
// around one side (merged into the freeform top) is cut along them and its
// normals, and the pieces become faces, where without the rung the region is
// halved and partly left as facets.
TEST(FitCadModes, ResegmentationCutsAMissedCrease) {
    const Mesh skin = analytic::hood_box(16);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    const Mesh merged = fcfix::merged_side(fcfix::fitted_abs(skin, tol));
    for (const auto mode : {FitCadFreeformMode::Adaptive, FitCadFreeformMode::Auto}) {
        SCOPED_TRACE(static_cast<int>(mode));
        FitCadOptions o;
        o.mBackend = FitCadBackend::Native;
        o.mTolerance = tol;
        o.mFreeformMode = mode;
        o.mResegmentRounds = 0;
        const FitCadResult off = fit_cad(merged, "", o);
        o.mResegmentRounds = 2;
        const FitCadResult on = fit_cad(merged, "", o);
        EXPECT_EQ(off.mNumResegmented, 0);
        EXPECT_GT(on.mNumResegmented, 0);
        EXPECT_GT(off.mFacetAreaFraction, 0.04);
        EXPECT_LT(on.mFacetAreaFraction, 0.01);
        EXPECT_LT(on.mNumFaces, off.mNumFaces);
        EXPECT_TRUE(on.mValid);
        // And the same result again.
        const FitCadResult again = fit_cad(merged, "", o);
        EXPECT_EQ(again.mNumFaces, on.mNumFaces);
        EXPECT_EQ(again.mFacetAreaFraction, on.mFacetAreaFraction);
    }
}

// A torus resolved finely enough for the tolerance keeps network patches:
// the default organic torus's facets sag from any smooth surface by about the
// tolerance in its tube's tight bends, so no patch through its vertices also
// holds its facet centres, and it falls back to adaptive faces.
TEST(FitCadModes, NetworkOnAResolvedTorus) {
    const Mesh skin = analytic::organic_torus(72, 36);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTolerance = tol;
    o.mFreeformMode = FitCadFreeformMode::Network;
    const FitCadResult r = fit_cad(fcfix::fitted_abs(skin, tol), "", o);
    EXPECT_GT(r.mNumNetworkPatches, 0);
    EXPECT_LT(r.mFacetAreaFraction, 0.01);
    EXPECT_TRUE(r.mValid);
}

// A rounded cube (|x|^4 + |y|^4 + |z|^4 = 1): next to its layout's valence-3
// vertices the tangent-plane treatment once bent the edges between them
// (G1 angles near 180 degrees); the network keeps G0 there now, so every
// mode gives a network whose angles stay those of G0.
TEST(FitCadModes, NetworkKeepsG0WhereTheTangentPlaneBends) {
    Mesh skin = icosa_sphere(1.0, 4);
    NDArray p = skin.Points();
    double* q = p.As<double>();
    for (std::size_t i = 0; i < skin.NumPoints(); ++i) {
        double s = 0.0;
        for (int j = 0; j < 3; ++j)
            s += q[3 * i + j] * q[3 * i + j] * q[3 * i + j] * q[3 * i + j];
        const double k = 1.0 / std::sqrt(std::sqrt(s));
        for (int j = 0; j < 3; ++j)
            q[3 * i + j] *= k;
    }
    skin.AssignPoints(std::move(p));
    const double tol = 2e-3 * fcfix::diagonal(skin);
    const Mesh fitted = fcfix::fitted_abs(skin, tol);
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTolerance = tol;
    o.mFreeformMode = FitCadFreeformMode::Network;
    o.mExtraordinaryVertex = FitCadExtraordinaryVertex::G0;
    const FitCadResult g0 = fit_cad(fitted, "", o);
    ASSERT_GT(g0.mNumExtraordinaryVertices, 0);
    for (const auto ev : {FitCadExtraordinaryVertex::G1, FitCadExtraordinaryVertex::G1x,
                          FitCadExtraordinaryVertex::G1xx}) {
        SCOPED_TRACE(static_cast<int>(ev));
        o.mExtraordinaryVertex = ev;
        const FitCadResult r = fit_cad(fitted, "", o);
        EXPECT_LT(r.mMaxG1AngleDeg, 90.0);
        EXPECT_LE(r.mMaxG1AngleDeg, g0.mMaxG1AngleDeg + 1e-9);
        EXPECT_EQ(r.mNumNetworkPatches, g0.mNumNetworkPatches);
    }
}

// A genus-2 skin in the network mode: its layout (cut along its handles)
// keeps network patches, the patches that miss become adaptive faces, and it
// leaves less area as facets than the adaptive faces alone.
TEST(FitCadModes, NetworkOnAGenusTwoSkin) {
    const Mesh skin = analytic::double_torus(20);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    const Mesh fitted = fcfix::fitted_abs(skin, tol);
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTolerance = tol;
    o.mFreeformMode = FitCadFreeformMode::Network;
    const FitCadResult net = fit_cad(fitted, "", o);
    o.mFreeformMode = FitCadFreeformMode::Adaptive;
    const FitCadResult adaptive = fit_cad(fitted, "", o);
    EXPECT_GT(net.mNumNetworkPatches, 0);
    EXPECT_TRUE(net.mValid);
    EXPECT_LT(net.mFacetAreaFraction, adaptive.mFacetAreaFraction);
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

// The saddle box through OpenCASCADE: its flat sides survive the checks and
// the saddle is a B-spline face in every freeform mode (its sides once folded
// over where the saddle dips more than a row, and failed with it).
TEST(FitCadModes, OcctSaddleBoxKeepsItsSides) {
    const Mesh skin = analytic::monkey_saddle_box(24);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    const Mesh fitted = fcfix::fitted_abs(skin, tol);
    FitCadOptions o;
    o.mTolerance = tol;
    for (const auto mode :
         {FitCadFreeformMode::Adaptive, FitCadFreeformMode::Network, FitCadFreeformMode::Auto}) {
        SCOPED_TRACE(static_cast<int>(mode));
        o.mFreeformMode = mode;
        const FitCadResult r = fit_cad(fitted, "", o);
        EXPECT_TRUE(r.mValid);
        EXPECT_LT(r.mFacetAreaFraction, 1e-6);
    }
}

// The re-segmentation rung through OpenCASCADE: the missed crease is cut and
// the faceted share falls.
TEST(FitCadModes, OcctResegmentationCutsAMissedCrease) {
    const Mesh skin = analytic::hood_box(16);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    const Mesh merged = fcfix::merged_side(fcfix::fitted_abs(skin, tol));
    FitCadOptions o;
    o.mTolerance = tol;
    o.mFreeformMode = FitCadFreeformMode::Auto;
    o.mResegmentRounds = 0;
    const FitCadResult off = fit_cad(merged, "", o);
    o.mResegmentRounds = 2;
    const FitCadResult on = fit_cad(merged, "", o);
    EXPECT_EQ(off.mNumResegmented, 0);
    EXPECT_GT(on.mNumResegmented, 0);
    EXPECT_LT(on.mFacetAreaFraction, off.mFacetAreaFraction);
    EXPECT_TRUE(on.mValid);
}

// The patch network in the OpenCASCADE backend: the Catmull-Clark cage's
// limit surface, which the adaptive faces leave as facets, becomes a valid
// solid of network faces, and the auto mode rebuilds its demoted adaptive
// faces in the network.
TEST(FitCadModes, OcctNetworkBuildsTheCage) {
    for (const fcfix::FreeformFixture& f : fcfix::freeform_fixtures(false)) {
        if (f.mName != "cc_valence3")
            continue;
        FitCadOptions o;
        o.mTolerance = f.mTolerance;
        const FitCadResult legacy = fit_cad(f.mMesh, "", o);
        for (FitCadFreeformMode mode : {FitCadFreeformMode::Network, FitCadFreeformMode::Auto}) {
            o.mFreeformMode = mode;
            const FitCadResult r = fit_cad(f.mMesh, "", o);
            EXPECT_TRUE(r.mValid);
            EXPECT_TRUE(r.mBrepChecked);
            EXPECT_LT(r.mFacetAreaFraction, legacy.mFacetAreaFraction);
            EXPECT_EQ(r.mFacetAreaFraction, 0.0);
            EXPECT_EQ(r.mNumNetworkPatches, r.mNumFaces);
            EXPECT_EQ(r.mNumExtraordinaryVertices, 8);
            EXPECT_GT(r.mNumNetworkEdges, 0);
            for (const FitCadFaceRecord& rec : r.mFaces) {
                EXPECT_GE(rec.mFace, 0);
                EXPECT_EQ(rec.mParamMethod, FitCadParamMethod::Square);
                EXPECT_TRUE(std::isfinite(rec.mMaxG1AngleDeg));
            }
        }
    }
}

// The native patch network's STEP file read back by OpenCASCADE: every face,
// its B-spline faces as B-splines, and the triangulation within the tolerance
// of the skin both ways.
TEST(FitCadModes, NativeNetworkStepReadsBack) {
    const Mesh skin = analytic::monkey_saddle_box(24);
    const double tol = 2e-3 * fcfix::diagonal(skin);
    FitCadOptions o;
    o.mBackend = FitCadBackend::Native;
    o.mTolerance = tol;
    o.mFreeformMode = FitCadFreeformMode::Network;
    const std::string path = temp_step("network_readback.step");
    const FitCadResult r = fit_cad(fcfix::fitted_abs(skin, tol), path, o);
    ASSERT_GT(r.mNumNetworkPatches, 0);
    BrepInfo info;
    CadReadOptions ro;
    ro.mDeflection = 0.1 * tol;
    const Mesh back = read_step(path, info, ro);
    EXPECT_EQ(static_cast<std::int64_t>(info.mFaces.size()), r.mNumFaces);
    std::int64_t cubic = 0;
    for (const BrepFace& f : info.mFaces) {
        const NurbsPatch& p = info.mPatches[static_cast<std::size_t>(f.mPatch)];
        cubic += p.mDegrees[0] == 3 && p.mDegrees[1] == 3 ? 1 : 0;
    }
    EXPECT_GE(cubic, r.mNumNetworkPatches);
    EXPECT_LE(hausdorff_distance(back, skin).mDistance, 1.5 * tol);
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
