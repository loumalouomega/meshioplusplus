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
 * @file test_octree_goldens.cpp
 * @brief Pins the output of `octree`, in 3-D and 2-D, before `octree_mesh`
 *        lands beside it (doc/octree_mesh_3d_plan.md).
 *
 * `octree_mesh` is a separate operation and does not touch `octree.cpp`, but
 * the two share the distance kernel, the triangle soup and the curvature and
 * feature passes; these digests prove that the existing octree and quadtree
 * output (every balance, transition and fill, snapping, the refinement
 * criteria, forests, padding and planes) stays byte for byte what it was. The
 * skins are built from integer arithmetic and square roots only. The
 * repeated-run checks hold everywhere, the golden digests on x86-64 GCC/Clang
 * only (no FMA contraction there), as in `test_op_goldens.cpp`.
 */

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "golden_fixtures.hpp"
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/surface.hpp"

namespace {

using namespace golden;
namespace mio = meshioplusplus;

// A jittered cube's skin (sharp edges and corners, no libm).
Mesh ogt_cube() {
    return mio::extract_surface(ft_tet_cube(6, 0.2));
}

// An icosphere (square roots only), shifted off the lattice planes.
Mesh ogt_sphere() {
    return mio::icosa_sphere(1.0, 2, {{0.03125, -0.0625, 0.015625}});
}

// A closed, slightly irregular octagon of `line` cells in the plane `Axis`
// (2: z = const, 1: y = const), coordinates dyadic.
Mesh ogt_octagon(int Axis) {
    const double xy[8][2] = {{1.0, 0.0},       {0.6875, 0.75},   {0.0, 1.0625},  {-0.75, 0.6875},
                             {-1.0, -0.03125}, {-0.6875, -0.75}, {0.0625, -1.0}, {0.75, -0.6875}};
    std::vector<std::vector<double>> pts;
    for (const auto& p : xy) {
        if (Axis == 2)
            pts.push_back({p[0], p[1], 0.25});
        else
            pts.push_back({p[0], -0.5, p[1]});
    }
    std::vector<std::vector<std::int64_t>> lines;
    for (std::int64_t i = 0; i < 8; ++i)
        lines.push_back({i, (i + 1) % 8});
    return mt::make_mesh(std::move(pts), "line", std::move(lines));
}

std::uint64_t ogt_digest(const mio::OctreeResult& rR) {
    MeshDigest d;
    d.Of(rR.mMesh);
    d.U64(static_cast<std::uint64_t>(rR.mDimension));
    d.Bytes(rR.mOrigin.data(), sizeof(double) * 3);
    d.Bytes(&rR.mRootSize, sizeof(double));
    d.U64(static_cast<std::uint64_t>(rR.mNumLeaves));
    d.U64(static_cast<std::uint64_t>(rR.mNumHanging));
    d.U64(static_cast<std::uint64_t>(rR.mFinestLevel));
    d.U64(static_cast<std::uint64_t>(rR.mNumSnapped));
    d.U64(static_cast<std::uint64_t>(rR.mNumRefused));
    return d.Value();
}

mio::OctreeOptions ogt_options(std::int64_t MaxLevel) {
    mio::OctreeOptions o;
    o.mMaxLevel = MaxLevel;
    o.mPaddingRelative = 0.0625;
    return o;
}

std::uint64_t ogt_hanging_face() {
    return ogt_digest(mio::octree(ogt_cube(), ogt_options(4)));
}

std::uint64_t ogt_inside_full_distance() {
    auto o = ogt_options(4);
    o.mBalance = mio::OctreeBalance::Full;
    o.mFill = mio::VoxelFill::Inside;
    o.mAttachDistance = true;
    return ogt_digest(mio::octree(ogt_sphere(), o));
}

std::uint64_t ogt_winding_surface() {
    auto o = ogt_options(4);
    o.mFill = mio::VoxelFill::Surface;
    o.mDistance.mSign = mio::SdfSign::WindingNumber;
    o.mAttachDistance = true;
    return ogt_digest(mio::octree(ogt_sphere(), o));
}

std::uint64_t ogt_templates() {
    auto o = ogt_options(4);
    o.mTransition = mio::OctreeTransition::Templates;
    return ogt_digest(mio::octree(ogt_sphere(), o));
}

std::uint64_t ogt_polyhedra_full() {
    auto o = ogt_options(4);
    o.mBalance = mio::OctreeBalance::Full;
    o.mTransition = mio::OctreeTransition::Polyhedra;
    o.mFill = mio::VoxelFill::Surface;
    return ogt_digest(mio::octree(ogt_cube(), o));
}

std::uint64_t ogt_snap() {
    auto o = ogt_options(3);
    o.mFill = mio::VoxelFill::Inside;
    o.mSnap = true;
    return ogt_digest(mio::octree(ogt_sphere(), o));
}

std::uint64_t ogt_criteria() {
    auto o = ogt_options(5);
    o.mSurfaceLevel = 3;
    o.mFeatureAngle = 30.0;
    return ogt_digest(mio::octree(ogt_cube(), o));
}

std::uint64_t ogt_curvature() {
    auto o = ogt_options(5);
    o.mSurfaceLevel = 3;
    o.mCurvatureThreshold = 1.5;
    return ogt_digest(mio::octree(ogt_sphere(), o));
}

std::uint64_t ogt_forest_band() {
    auto o = ogt_options(3);
    o.mRootCells = {{2, 1, 1}};
    o.mPadding = 0.125;
    o.mSampleDistance = 0.0625;
    o.mBalance = mio::OctreeBalance::None;
    return ogt_digest(mio::octree(ogt_cube(), o));
}

std::uint64_t ogt_quad_inside() {
    auto o = ogt_options(5);
    o.mFill = mio::VoxelFill::Inside;
    o.mAttachDistance = true;
    return ogt_digest(mio::octree(ogt_octagon(2), o));
}

std::uint64_t ogt_quad_templates() {
    auto o = ogt_options(5);
    o.mTransition = mio::OctreeTransition::Templates;
    return ogt_digest(mio::octree(ogt_octagon(2), o));
}

std::uint64_t ogt_quad_polygons_full() {
    auto o = ogt_options(5);
    o.mBalance = mio::OctreeBalance::Full;
    o.mTransition = mio::OctreeTransition::Polyhedra;
    o.mFill = mio::VoxelFill::Surface;
    return ogt_digest(mio::octree(ogt_octagon(2), o));
}

std::uint64_t ogt_quad_snap_even_odd() {
    auto o = ogt_options(4);
    o.mFill = mio::VoxelFill::Inside;
    o.mEvenOdd = true;
    o.mSnap = true;
    return ogt_digest(mio::octree(ogt_octagon(2), o));
}

std::uint64_t ogt_quad_plane_forest() {
    auto o = ogt_options(4);
    o.mRootCells = {{2, 3, 1}};
    o.mPadding = 0.125;
    o.mPlane = mio::detail::PlanarPlane::XZ;
    return ogt_digest(mio::octree(ogt_octagon(1), o));
}

struct OgtCase {
    const char* mName;
    std::uint64_t (*mRun)();
    std::uint64_t mGolden;
};

const OgtCase kOgtCases[] = {
    {"hanging_face", ogt_hanging_face, 0x1b6f3e935694539cull},
    {"inside_full_distance", ogt_inside_full_distance, 0x60064f3c8cdcfc6bull},
    {"winding_surface", ogt_winding_surface, 0x8a20e18335396c84ull},
    {"templates", ogt_templates, 0xea9ed4cfe94411a6ull},
    {"polyhedra_full", ogt_polyhedra_full, 0x3b782b6cd8dcc541ull},
    {"snap", ogt_snap, 0x230ec29df5e3f93eull},
    {"criteria", ogt_criteria, 0xe683ddcb780b5e0eull},
    {"curvature", ogt_curvature, 0xd561989e683c3769ull},
    {"forest_band", ogt_forest_band, 0xb8dd61337e1fda03ull},
    {"quad_inside", ogt_quad_inside, 0xa39bafe4b7843a4full},
    {"quad_templates", ogt_quad_templates, 0x0c68ddc958d943c6ull},
    {"quad_polygons_full", ogt_quad_polygons_full, 0xc24abff76f3b175cull},
    {"quad_snap_even_odd", ogt_quad_snap_even_odd, 0xf7e4e38e43e526a3ull},
    {"quad_plane_forest", ogt_quad_plane_forest, 0x8c89cc346c71370eull},
};

TEST(OctreeGoldens, ResultsAreStableAcrossRepeatedRuns) {
    for (const OgtCase& c : kOgtCases) {
        const std::uint64_t first = c.mRun();
        EXPECT_EQ(c.mRun(), first) << c.mName;
    }
}

// The fixtures reach the paths they pin.
TEST(OctreeGoldens, FixturesExerciseTheirPaths) {
    auto o = ogt_options(3);
    o.mFill = mio::VoxelFill::Inside;
    o.mSnap = true;
    EXPECT_GT(mio::octree(ogt_sphere(), o).mNumSnapped, 0);
    EXPECT_GT(mio::octree(ogt_cube(), ogt_options(4)).mNumHanging, 0);
    EXPECT_EQ(mio::octree(ogt_octagon(2), ogt_options(5)).mDimension, 2);
    auto t = ogt_options(4);
    t.mTransition = mio::OctreeTransition::Templates;
    EXPECT_GE(mio::octree(ogt_sphere(), t).mMesh.NumCellBlocks(), 2u);
}

#ifdef GOLDEN_PINNED
TEST(OctreeGoldens, ResultsMatchThePreviousImplementation) {
    for (const OgtCase& c : kOgtCases) {
        const std::uint64_t got = c.mRun();
        EXPECT_EQ(got, c.mGolden) << c.mName << ": 0x" << std::hex << got;
    }
}
#endif

}  // namespace
