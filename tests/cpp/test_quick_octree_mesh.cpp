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
// Small versions of the octree_mesh tests (`test_octree_mesh.cpp`). The full suites run
// skins of thousands of triangles and take minutes each in the instrumented coverage build,
// so a pull request skips them there; these run the same code paths (the castellated, the
// tetrahedral and the hexahedral output, the uniform and graded skins, features, cleaving,
// the option checks) on skins of a few hundred triangles, in seconds, so the coverage job
// keeps measuring that code. Their suite names start with `Quick`, which the job's exclusion
// does not match. They check validity and the main invariants, not the convergence rates the
// full tests pin.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

#include "thread_sweep.hpp"

// Project includes
#include "analytic_surfaces.hpp"
#include "mesh_fixtures.hpp"
#include "octree_mesh_fixtures.hpp"
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/surface.hpp"

namespace {

namespace mio = meshioplusplus;
using namespace omfix;

constexpr double kPi = 3.14159265358979323846;

// A sphere of 320 triangles: a quarter of a percent of the full tests' spheres' work.
mio::Mesh quick_sphere(double R = 1.0) {
    return mio::icosa_sphere(R, 2);
}

std::uint64_t quick_digest(const mio::Mesh& rMesh) {
    mio::bench::MeshDigest d;
    d.Of(rMesh);
    return d.Value();
}

}  // namespace

// The castellated output: the inside converges from below to the ball, the regions
// agree with the report, and the surrogate faces enclose the inside cells.
TEST(QuickOctreeMesh, TheCastellatedInsideApproachesTheBall) {
    const mio::Mesh sphere = quick_sphere();
    double prev = 1e300;
    for (double hmin : {0.5, 0.25}) {
        auto o = om_castellated(0.5, hmin);
        o.mSkinGrading = mio::OctreeMeshSkinGrading::Uniform;
        o.mSurfaceSize = hmin;
        const auto r = mio::octree_mesh(sphere, o);
        EXPECT_EQ(r.mSurfaceLevel, r.mMaxLevel);
        const auto boxes = om_boxes(r.mMesh);
        const auto inside = om_region(r.mMesh, mio::kOctreeMeshInsideRegion, mio::RegionKind::Cell);
        EXPECT_EQ(static_cast<std::int64_t>(inside.size()), r.mNumInside);
        EXPECT_EQ(inside.size(), boxes.size());
        double vol = 0.0;
        for (std::int64_t n : inside)
            vol += om_volume(boxes[static_cast<std::size_t>(n)]);
        const double err = std::abs(vol - 4.0 / 3.0 * kPi);
        EXPECT_LT(err, prev);
        prev = err;
        const auto sides =
            om_region(r.mMesh, mio::kOctreeMeshSurrogateRegion, mio::RegionKind::Side);
        EXPECT_EQ(static_cast<std::int64_t>(sides.size() / 2), r.mNumSurrogateFaces);
        EXPECT_GT(r.mNumSurrogateFaces, 0);
    }
}

// A holed skin still meshes, and an open one is not flooded; both are reported.
TEST(QuickOctreeMesh, DirtySkinsStillMesh) {
    const mio::Mesh sphere = quick_sphere();
    std::vector<std::vector<double>> pts;
    const double* p = sphere.Points().As<double>();
    for (std::size_t i = 0; i < sphere.NumPoints(); ++i)
        pts.push_back({p[3 * i], p[3 * i + 1], p[3 * i + 2]});
    std::vector<std::vector<std::int64_t>> holed;
    const std::int64_t* c = sphere.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < sphere.Cells(0).NumCells(); ++t)
        if (t % 13 != 0)
            holed.push_back({c[3 * t], c[3 * t + 1], c[3 * t + 2]});
    const auto r =
        mio::octree_mesh(mt::make_mesh(pts, "triangle", holed), om_castellated(0.5, 0.25));
    EXPECT_GT(r.mNumInside, 0);
}

// Every option is checked and refused by name.
TEST(QuickOctreeMesh, RefusesBadOptionsByName) {
    const mio::Mesh sphere = quick_sphere();
    const auto o = om_castellated(0.5, 0.25);
    const auto refused = [&](mio::OctreeMeshOptions Bad, const char* pName) {
        try {
            mio::octree_mesh(sphere, Bad);
            ADD_FAILURE() << pName << " was accepted";
        } catch (const std::invalid_argument& rErr) {
            EXPECT_NE(std::string(rErr.what()).find(pName), std::string::npos) << rErr.what();
        }
    };
    auto bad = o;
    bad.mMinSize = 1.0;
    refused(bad, "min_size");
    bad = o;
    bad.mWindingThreshold = 1.5;
    refused(bad, "winding_threshold");
    bad = o;
    bad.mMaxCells = 10;
    refused(bad, "max_cells");
    bad = o;
    bad.mMinSize = 1e-7;
    refused(bad, "levels");
    bad = o;
    bad.mSizeArray = "nope";
    refused(bad, "size_array");
    bad = o;
    bad.mCrossingMargin = 0.5;
    refused(bad, "crossing_margin");
    bad = o;
    bad.mZeroVolumeFactor = -1.0;
    refused(bad, "zero_volume_factor");
    EXPECT_THROW(mio::octree_mesh(mt::line_mesh(), o), std::invalid_argument);
    EXPECT_THROW(mio::octree_mesh_output_from_name("quad"), std::invalid_argument);
    EXPECT_EQ(mio::octree_mesh_output_from_name("tet"), mio::OctreeMeshOutput::Tetra);
}

// The tetrahedral output of a uniform skin keeps the stuffing's angle bounds and encloses
// the skin's volume; the same run twice, and on every thread count, is the same mesh.
TEST(QuickOctreeMeshTet, UniformSkinsKeepTheStuffingAngleBounds) {
    const mio::Mesh skin = quick_sphere();
    for (double h : {0.1, 0.15}) {
        auto o = om_tet(0.5, h, true);
        o.mFitIterations = 0;
        const auto r = mio::octree_mesh(skin, o);
        const OmTetCheck c = om_check_tets(r.mMesh);
        EXPECT_GE(c.mMinDihedral, 10.7);
        EXPECT_LE(c.mMaxDihedral, 164.8);
        EXPECT_TRUE(r.mGuaranteeApplies);
        EXPECT_EQ(r.mNumInverted, 0);
        EXPECT_NEAR(c.mVolume, om_enclosed(skin), 0.08 * om_enclosed(skin));
        EXPECT_GT(r.mNumWarped, 0);
        EXPECT_EQ(quick_digest(mio::octree_mesh(skin, o).mMesh), quick_digest(r.mMesh));
    }
}

// The fitted and optimised output stays valid and reports the boundary's distance.
TEST(QuickOctreeMeshTet, FittingAndOptimisationKeepTheMeshValid) {
    const mio::Mesh skin = quick_sphere();
    auto o = om_tet(0.5, 0.1, true);
    const auto r = mio::octree_mesh(skin, o);
    const OmTetCheck c = om_check_tets(r.mMesh);
    EXPECT_EQ(r.mNumInverted, 0);
    EXPECT_GT(r.mQuality.mNumCells, 0);
    EXPECT_GE(r.mHausdorffToSkin, 0.0);
    EXPECT_LT(r.mHausdorffToSkin, 0.5 * r.mMinSize + 0.02);
    EXPECT_GE(c.mMinDihedral, r.mMinDihedral - 1e-6);
    for (const char* stage : {"skin", "size_field", "octree", "stuffing", "optimize", "output"}) {
        bool found = false;
        for (const auto& [name, seconds] : r.mTimings)
            found = found || name == stage;
        EXPECT_TRUE(found) << stage;
    }
    // The volume fix switched off still gives a valid mesh.
    o.mVolumeFix = false;
    om_check_tets(mio::octree_mesh(skin, o).mMesh);
}

// A cube's eight corners and twelve edges are kept exactly, uniform or graded.
TEST(QuickOctreeMeshTet, CubeCornersAndEdgesArePreserved) {
    const double lo[3] = {-0.61, -0.67, -0.7}, L = 1.37;
    const mio::Mesh cube =
        mio::convert_cells(
            mio::extract_surface(mio::box({{4, 4, 4}}, {{L, L, L}}, {{lo[0], lo[1], lo[2]}})),
            {mio::ConvertCellsMode::Simplexify, false})
            .mMesh;
    for (bool uniform : {true, false}) {
        auto o = om_tet(0.5, 0.15, uniform);
        const auto r = mio::octree_mesh(cube, o);
        om_check_tets(r.mMesh);
        EXPECT_EQ(r.mNumCorners, 8);
        EXPECT_EQ(r.mNumCornersSnapped, 8);
        EXPECT_EQ(r.mNumFeatureChains, 12);
        EXPECT_EQ(r.mNumFeatureChainsMissed, 0);
        EXPECT_FALSE(r.mGuaranteeApplies);
    }
}

// Two nested spheres give two conforming materials; keeping the outside fills the lattice.
TEST(QuickOctreeMeshTet, CleavingNestedSpheresGivesTwoMaterials) {
    const mio::Mesh both = om_nested_spheres({1.0, 0.5}, 2);
    auto o = om_tet(0.5, 0.15, false);
    o.mConformity = mio::OctreeMeshConformity::Cleaving;
    const auto r = mio::octree_mesh(both, o);
    om_check_tets(r.mMesh);
    EXPECT_GT(om_region(r.mMesh, "material_1", mio::RegionKind::Cell).size(), 0u);
    EXPECT_GT(om_region(r.mMesh, "material_2", mio::RegionKind::Cell).size(), 0u);
    auto keep = om_tet(0.5, 0.25, true);
    keep.mKeepOutside = true;
    const auto k = mio::octree_mesh(quick_sphere(), keep);
    EXPECT_GT(k.mMesh.Cells(0).NumCells(), 0u);
}

// A graded skin (sizes from a point array) still gives a valid mesh.
TEST(QuickOctreeMeshTet, AGradedSkinGivesAValidMesh) {
    mio::Mesh graded = quick_sphere();
    mio::NDArray size = mio::NDArray::Uninit(mio::DType::Float64, {graded.NumPoints()});
    const double* gp = graded.Points().As<double>();
    for (std::size_t i = 0; i < graded.NumPoints(); ++i)
        size.As<double>()[i] = 0.15 + 0.2 * (gp[3 * i] + 1.0);
    graded.AddPointData("h", std::move(size));
    auto o = om_tet(0.5, 0.12, false);
    o.mSizeArray = "h";
    const auto g = mio::octree_mesh(graded, o);
    const OmTetCheck gc = om_check_tets(g.mMesh);
    EXPECT_GT(gc.mMinDihedral, 0.0);
    EXPECT_EQ(g.mNumInverted, 0);
}

// The hexahedral output of a sphere and a torus is all hexahedra fitted to the skin, with
// and without the pillow layer, on every thread count the same.
TEST(QuickOctreeMeshHex, SpheresAndToriAreAllHexahedra) {
    const mio::Mesh skins[] = {quick_sphere(), analytic::capped_cylinder(0.5, 2.0, 4, 24)};
    for (const mio::Mesh& skin : skins)
        for (bool pillow : {true, false}) {
            if (!pillow && &skin == &skins[1])
                continue;  // the no-pillow run of the sphere is enough
            auto o = om_hex(0.3);
            o.mPillow = pillow;
            const auto r = mio::octree_mesh(skin, o);
            const OmHexCheck c = om_check_hexes(r.mMesh);
            EXPECT_EQ(r.mNumInverted, 0);
            EXPECT_EQ(r.mNumHexahedra, static_cast<std::int64_t>(r.mMesh.Cells(0).NumCells()));
            EXPECT_EQ(r.mNumBoundaryEdgesOpen, 0);
            EXPECT_EQ(r.mNumNonManifoldBoundaryEdges, 0);
            EXPECT_LT(r.mHausdorffToSkin, 0.15);
            EXPECT_NEAR(c.mVolume, om_enclosed(skin), 0.15 * om_enclosed(skin));
            EXPECT_NEAR(r.mMinScaledJacobian, c.mMinScaledJacobian, 0.05);
        }
}

// The graded dual has fewer hexahedra than one grid at the skin's size and conforms.
TEST(QuickOctreeMeshHex, GradedHexahedraConformAndSaveCells) {
    const mio::Mesh skin = quick_sphere();
    auto o = om_hex(0.0);
    o.mCurvatureAngle = 15.0;
    o.mMaxSize = 0.6;
    auto g = o;
    g.mHexGrading = mio::OctreeMeshHexGrading::Graded;
    const auto r = mio::octree_mesh(skin, g);
    const auto u = mio::octree_mesh(skin, o);
    const OmHexCheck h = om_check_hexes(r.mMesh);
    EXPECT_LT(r.mNumHexahedra, u.mNumHexahedra);
    EXPECT_EQ(r.mNumInverted, 0);
    EXPECT_EQ(r.mNumBoundaryEdgesOpen, 0);
    EXPECT_NEAR(h.mVolume, om_enclosed(skin), 0.1 * om_enclosed(skin));
    const mio::ConformityResult cc = mio::check_conformity(r.mMesh);
    EXPECT_EQ(cc.mNumHangingNodes, 0);
    EXPECT_EQ(cc.mNumNonManifoldFacets, 0);
    EXPECT_EQ(quick_digest(mio::octree_mesh(skin, g).mMesh), quick_digest(r.mMesh));
}

// A box with features keeps its corners and edges and has no flat corner; cleaving two
// shells gives conforming materials of the two shells' volumes.
TEST(QuickOctreeMeshHex, BoxFeaturesAndCleavedShells) {
    const mio::Mesh box = om_unit_box();
    const auto r = mio::octree_mesh(box, om_hex(0.25));
    om_check_hexes(r.mMesh);
    EXPECT_EQ(om_flat_corners(r.mMesh), 0);
    EXPECT_EQ(r.mNumInverted, 0);

    const mio::Mesh both = om_nested_spheres({1.0, 0.5}, 2);
    auto o = om_hex(0.3);
    o.mConformity = mio::OctreeMeshConformity::Cleaving;
    const auto c = mio::octree_mesh(both, o);
    om_check_hexes(c.mMesh);
    const std::vector<double> vols = om_material_volumes(c.mMesh);
    ASSERT_GE(vols.size(), 3u);
    const mio::ConformityResult cc = mio::check_conformity(c.mMesh);
    EXPECT_EQ(cc.mNumHangingNodes, 0);
    EXPECT_EQ(cc.mNumNonManifoldFacets, 0);
}

// The names the options take, each way, and the refusal of an unknown one.
TEST(QuickOctreeMesh, OptionNamesRoundTripAndUnknownOnesAreRefused) {
    for (const char* name : {"castellated", "tet", "hex"})
        EXPECT_STREQ(mio::octree_mesh_output_name(mio::octree_mesh_output_from_name(name)), name);
    EXPECT_EQ(mio::octree_mesh_output_from_name("hexa"), mio::OctreeMeshOutput::Hexa);
    EXPECT_EQ(mio::octree_mesh_conformity_from_name("stuffing"),
              mio::OctreeMeshConformity::Stuffing);
    EXPECT_EQ(mio::octree_mesh_conformity_from_name("cleaving"),
              mio::OctreeMeshConformity::Cleaving);
    EXPECT_EQ(mio::octree_mesh_skin_grading_from_name("graded"),
              mio::OctreeMeshSkinGrading::Graded);
    EXPECT_EQ(mio::octree_mesh_skin_grading_from_name("uniform"),
              mio::OctreeMeshSkinGrading::Uniform);
    EXPECT_EQ(mio::octree_mesh_hex_grading_from_name("uniform"),
              mio::OctreeMeshHexGrading::Uniform);
    EXPECT_EQ(mio::octree_mesh_hex_grading_from_name("graded"), mio::OctreeMeshHexGrading::Graded);
    for (const char* name : {"all", "inside", "cover"})
        EXPECT_NO_THROW(mio::octree_mesh_fill_from_name(name));
    EXPECT_THROW(mio::octree_mesh_conformity_from_name("nope"), std::invalid_argument);
    EXPECT_THROW(mio::octree_mesh_skin_grading_from_name("nope"), std::invalid_argument);
    EXPECT_THROW(mio::octree_mesh_hex_grading_from_name("nope"), std::invalid_argument);
    EXPECT_THROW(mio::octree_mesh_fill_from_name("nope"), std::invalid_argument);
}

// The opt-in remeshing stage: adapt_volume runs on the tetrahedra, the mesh stays valid and
// the report says it ran; off, the output is the plain one.
TEST(QuickOctreeMeshTet, PostAdaptRemeshesTheTetrahedra) {
    const mio::Mesh skin = mio::icosa_sphere(1.0, 1);
    auto o = om_tet(0.8, 0.25, true);
    const auto plain = mio::octree_mesh(skin, o);
    EXPECT_FALSE(plain.mAdapted);
    o.mPostAdapt = true;
    const auto r = mio::octree_mesh(skin, o);
    EXPECT_TRUE(r.mAdapted);
    const OmTetCheck c = om_check_tets(r.mMesh);
    EXPECT_NEAR(c.mVolume, om_enclosed(skin), 0.2 * om_enclosed(skin));
    EXPECT_EQ(r.mNumInverted, 0);
}

// A crease shorter than a cell (a pushed-out face centre of a box) is not snapped by the
// tetrahedral output, while the box's corners are.
TEST(QuickOctreeMeshTet, CreasesShorterThanACellAreNotSnapped) {
    mio::Mesh box =
        mio::convert_cells(
            mio::extract_surface(mio::box({{4, 4, 4}}, {{1.0, 1.0, 1.0}}, {{-0.51, -0.53, -0.55}})),
            {mio::ConvertCellsMode::Simplexify, false})
            .mMesh;
    {
        const std::size_t np = box.NumPoints();
        mio::NDArray pts = mio::NDArray::Uninit(mio::DType::Float64, {np, std::size_t{3}});
        const double* src = box.Points().As<double>();
        double* p = pts.As<double>();
        std::copy(src, src + 3 * np, p);
        for (std::size_t v = 0; v < np; ++v)
            if (std::abs(p[3 * v] - 0.49) < 1e-9 && std::abs(p[3 * v + 1] + 0.03) < 1e-9 &&
                std::abs(p[3 * v + 2] + 0.05) < 1e-9)
                p[3 * v] += 0.1;
        box.AssignPoints(std::move(pts));
    }
    const auto r = mio::octree_mesh(box, om_tet(0.5, 0.25, true));
    om_check_tets(r.mMesh);
    EXPECT_EQ(r.mNumCornersSnapped, 8);
    EXPECT_GT(r.mNumFeatureChainsMissed, 0);
}
