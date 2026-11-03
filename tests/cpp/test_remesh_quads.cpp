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
/**
 * @file test_remesh_quads.cpp
 * @brief `remesh_quads` (roadmap 7.7.3): the quad remesher registry through a
 *        fake provider, the refusals, and the roadmap's probe through the
 *        installed QuadriFlow (skipped without it).
 */

// System includes
#include <algorithm>
#include <cstdio>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "brep_solids.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/remesh_quads.hpp"
#include "meshioplusplus/providers.hpp"

namespace {

namespace mio = meshioplusplus;

/// The unit cube's skin as 12 outward triangles.
mio::Mesh rq_cube_triangles() {
    return mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "triangle",
        {{0, 2, 1},
         {0, 3, 2},
         {4, 5, 6},
         {4, 6, 7},
         {0, 1, 5},
         {0, 5, 4},
         {1, 2, 6},
         {1, 6, 5},
         {2, 3, 7},
         {2, 7, 6},
         {3, 0, 4},
         {3, 4, 7}});
}

/// A remesher returning the cube as its six quads (and, when asked to be
/// dominant, one face as two triangles instead), recording the request.
mio::QuadRemesherProvider rq_fake(mio::QuadRemeshRequest& rSeen) {
    mio::QuadRemesherProvider p;
    p.mName = "fake-quads";
    p.mVersion = "0.1";
    p.mRemesh = [&rSeen](const mio::QuadRemeshRequest& rRequest) {
        rSeen = rRequest;
        mio::QuadRemeshResult out;
        out.mPoints = rRequest.mPoints;
        out.mQuads = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}};
        if (rRequest.mDominant)
            out.mTriangles = {{3, 0, 4}, {3, 4, 7}};
        else
            out.mQuads.push_back({3, 0, 4, 7});
        return out;
    };
    return p;
}

}  // namespace

TEST(RemeshQuads, AFakeRemesherThroughTheRegistry) {
    mio::QuadRemeshRequest seen;
    mio::add_quad_remesher(rq_fake(seen));
    const auto names = mio::remesh_quads_backends();
    EXPECT_NE(std::find(names.begin(), names.end(), "fake-quads"), names.end());

    mio::RemeshQuadsOptions o;
    o.mBackend = "fake-quads";
    o.mSeed = 7;
    const auto r = mio::remesh_quads(rq_cube_triangles(), o);
    EXPECT_EQ(r.mBackend, "fake-quads");
    EXPECT_EQ(r.mNumQuads, 6);
    EXPECT_EQ(r.mNumTriangles, 0);
    EXPECT_EQ(r.mNumPoints, 8);
    EXPECT_DOUBLE_EQ(r.mMinScaledJacobian, 1.0);
    EXPECT_DOUBLE_EQ(r.mMeanScaledJacobian, 1.0);
    EXPECT_NEAR(r.mHausdorffToSurface, 0.0, 1e-12);
    EXPECT_TRUE(r.mWarnings.empty());
    // The automatic target is half the triangles, and the options go through.
    EXPECT_EQ(seen.mTargetFaces, 6);
    EXPECT_EQ(seen.mSeed, 7);
    EXPECT_TRUE(seen.mBoundary);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "quad");

    o.mDominant = true;
    o.mTargetFaces = 9;
    const auto d = mio::remesh_quads(rq_cube_triangles(), o);
    EXPECT_EQ(seen.mTargetFaces, 9);
    EXPECT_EQ(d.mNumQuads, 5);
    EXPECT_EQ(d.mNumTriangles, 2);
    EXPECT_EQ(d.mWarnings.size(), 1u);
    mio::remove_quad_remesher("fake-quads");
}

TEST(RemeshQuads, RefusalsAreByName) {
    mio::RemeshQuadsOptions o;
    o.mTargetFaces = -1;
    EXPECT_THROW(mio::remesh_quads(rq_cube_triangles(), o), std::invalid_argument);
    o.mTargetFaces = 100;
    o.mMaxCells = 10;
    EXPECT_THROW(mio::remesh_quads(rq_cube_triangles(), o), std::invalid_argument);
    const mio::Mesh quads =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "quad", {{0, 1, 2, 3}});
    try {
        mio::remesh_quads(quads, {});
        FAIL() << "expected invalid_argument";
    } catch (const std::invalid_argument& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("'quad' block"), std::string::npos);
    }
    mio::RemeshQuadsOptions named;
    named.mBackend = "no-such-remesher";
    try {
        mio::remesh_quads(rq_cube_triangles(), named);
        FAIL() << "expected Unsupported";
    } catch (const mio::Unsupported& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("'no-such-remesher'"), std::string::npos);
    }
    // A remesher that writes a face past its points is refused, naming it.
    mio::QuadRemesherProvider bad;
    bad.mName = "fake-bad";
    bad.mRemesh = [](const mio::QuadRemeshRequest& rRequest) {
        mio::QuadRemeshResult out;
        out.mPoints = rRequest.mPoints;
        out.mQuads = {{0, 1, 2, 99}};
        return out;
    };
    mio::add_quad_remesher(std::move(bad));
    named.mBackend = "fake-bad";
    try {
        mio::remesh_quads(rq_cube_triangles(), named);
        FAIL() << "expected runtime_error";
    } catch (const std::runtime_error& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("fake-bad"), std::string::npos);
    }
    mio::remove_quad_remesher("fake-bad");
}

// The roadmap's probe: QuadriFlow makes the sphere all quads, square nearly
// everywhere (mean scaled Jacobian ~0.99), its worst quad far above what
// recombining mesh_cad_surface's own triangles leaves at the poles (0.12).
TEST(RemeshQuads, QuadriFlowMakesAnAllQuadSphere) {
    if (!mio::quad_remesher("quadriflow"))
        GTEST_SKIP() << "quadriflow is not installed";
    mio::MeshCadSurfaceOptions mo;
    mo.mTargetSize = 0.1;
    mo.mEdges = false;
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    const mio::Mesh skin = mio::mesh_cad_surface(sphere, mo).mMesh;
    mio::MeshCadSurfaceExtras all;
    all.mQuads = mio::MeshCadQuads::All;
    mio::MeshCadSurfaceExtrasReport xq;
    mio::mesh_cad_surface(sphere, mo, mio::MeshCadSizeSources{}, all, nullptr, &xq);
    // The least scaled Jacobian depends on QuadriFlow's seed; the mean, and
    // the margin over recombination, do not.
    mio::RemeshQuadsOptions o;
    o.mBackend = "quadriflow";
    mio::RemeshQuadsResult r;
    for (const int seed : {1, 0}) {
        o.mSeed = seed;
        r = mio::remesh_quads(skin, o);
        std::printf(
            "[7.7.3] sphere at h = 0.1, seed %d: quadriflow %lld quads, scaled Jacobian "
            "%.3f (mean %.3f), Hausdorff %.4f; recombination %lld quads, %.3f\n",
            seed, static_cast<long long>(r.mNumQuads), r.mMinScaledJacobian, r.mMeanScaledJacobian,
            r.mHausdorffToSurface, static_cast<long long>(xq.mNumQuads), xq.mMinQuadScaledJacobian);
        EXPECT_EQ(r.mBackend, "quadriflow");
        EXPECT_GT(r.mNumQuads, 0);
        EXPECT_EQ(r.mNumTriangles, 0);
        EXPECT_GE(r.mMinScaledJacobian, 0.2) << seed;
        EXPECT_GE(r.mMeanScaledJacobian, 0.95) << seed;
        EXPECT_GT(r.mMinScaledJacobian, xq.mMinQuadScaledJacobian) << seed;
        EXPECT_LT(r.mHausdorffToSurface, 0.05) << seed;
    }
    // The same seed repeats the same mesh.
    const auto again = mio::remesh_quads(skin, o);
    EXPECT_EQ(again.mNumQuads, r.mNumQuads);
    EXPECT_EQ(again.mMinScaledJacobian, r.mMinScaledJacobian);
}
