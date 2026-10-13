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

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "analytic_surfaces.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

Mesh fitted(const Mesh& rSkin, double CurvatureTolerance = 0.0, double Tolerance = 1e-3,
            std::int64_t MinPatchCells = 1) {
    SegmentOptions so;
    so.mCurvatureTolerance = CurvatureTolerance;
    so.mMinPatchCells = MinPatchCells;
    FitPrimitivesOptions fo;
    fo.mTolerance = Tolerance;
    return fit_primitives(segment(rSkin, so).mMesh, fo).mMesh;
}

std::string read_bytes(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// A closed box whose top is the saddle z = 1 + 0.2 (x - 1/2)(y - 1/2) * 4: no
// primitive explains the top, so it must come back faceted.
Mesh saddle_box(int n) {
    Mesh box = analytic::box_skin(n);
    const NDArray& p = box.Points();
    NDArray q = NDArray::Uninit(DType::Float64, {box.NumPoints(), 3});
    for (std::size_t i = 0; i < box.NumPoints(); ++i) {
        const double x = p.As<double>()[i * 3], y = p.As<double>()[i * 3 + 1];
        double z = p.As<double>()[i * 3 + 2];
        z *= 1.0 + 0.8 * (x - 0.5) * (y - 0.5);
        q.As<double>()[i * 3] = x;
        q.As<double>()[i * 3 + 1] = y;
        q.As<double>()[i * 3 + 2] = z;
    }
    box.AssignPoints(std::move(q));
    return box;
}

}  // namespace

TEST(FitCad, CapabilityMatchesTheBuild) {
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    EXPECT_TRUE(fit_cad_has_occt());
#else
    EXPECT_FALSE(fit_cad_has_occt());
#endif
}

#ifndef MESHIOPLUSPLUS_HAS_OCCT

TEST(FitCad, AbsentBackendIsNamed) {
    try {
        fit_cad(analytic::box_skin(2));
        FAIL() << "expected Unsupported";
    } catch (const Unsupported& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("MESHIOPLUSPLUS_WITH_OCCT"), std::string::npos);
    }
}

#else

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

TEST(FitCad, UnexplainedPatchComesBackFaceted) {
    const FitCadResult r = fit_cad(fitted(saddle_box(6)));
    EXPECT_TRUE(r.mValid);
    EXPECT_GT(r.mNumFacetedFaces, 0);
    EXPECT_EQ(r.mNumPlane, 5);
    EXPECT_NEAR(r.mVolumeRatio, 1.0, 1e-9);
    EXPECT_LE(r.mMaxDeviation, r.mTolerance);
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
}

#endif
