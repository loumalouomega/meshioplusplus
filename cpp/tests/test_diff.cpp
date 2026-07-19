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
// GoogleTest suite for the mesh comparison operation (operations/diff.hpp).
// Runs under every mesh backend via the cpp-tests matrix.

// System includes
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/operations/diff.hpp"
#include "mesh_fixtures.hpp"

using meshioplusplus::diff;
using meshioplusplus::DiffOptions;
using meshioplusplus::DiffReport;
using meshioplusplus::DiffVerdict;
using meshioplusplus::Mesh;
using meshioplusplus::meshes_equal;
using meshioplusplus::NDArray;

namespace {

// A small triangle mesh whose points can be perturbed.
Mesh tri_with_points(const std::vector<std::vector<double>>& rPts) {
    return mt::make_mesh(rPts, "triangle", {{0, 1, 2}, {0, 2, 3}});
}

const std::vector<std::vector<double>> kSquare = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};

TEST(Diff, IdenticalMeshes) {
    Mesh a = tri_with_points(kSquare);
    Mesh b = tri_with_points(kSquare);
    DiffReport rep = diff(a, b);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::Identical);
    EXPECT_TRUE(meshes_equal(a, b));
    EXPECT_EQ(meshioplusplus::diff_verdict_name(rep.mVerdict), "identical");
    EXPECT_FALSE(rep.mPointCountMismatch);
    EXPECT_FALSE(rep.mBlockCountMismatch);
}

TEST(Diff, PerturbBelowTolerance) {
    std::vector<std::vector<double>> pts = kSquare;
    pts[2][1] += 1e-10;  // below atol default? use explicit tolerance
    Mesh a = tri_with_points(kSquare);
    Mesh b = tri_with_points(pts);
    DiffOptions opts;
    opts.atol = 1e-8;
    opts.rtol = 0.0;
    DiffReport rep = diff(a, b, opts);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::EqualWithinTolerance);
    EXPECT_TRUE(meshes_equal(a, b, 1e-8, 0.0));
    EXPECT_GT(rep.mPoints.mMaxAbsError, 0.0);
    EXPECT_EQ(rep.mPoints.mNumExceeding, 0);
}

TEST(Diff, PerturbAboveTolerance) {
    std::vector<std::vector<double>> pts = kSquare;
    pts[2][1] += 1.0;  // point index 2, coordinate y
    Mesh a = tri_with_points(kSquare);
    Mesh b = tri_with_points(pts);
    DiffOptions opts;
    opts.atol = 1e-8;
    opts.rtol = 1e-6;
    DiffReport rep = diff(a, b, opts);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::Different);
    EXPECT_FALSE(meshes_equal(a, b, 1e-8, 1e-6));
    EXPECT_GE(rep.mPoints.mNumExceeding, 1);
    // Worst offender is the y-coordinate of point 2 -> flat index 2*3 + 1 = 7.
    EXPECT_EQ(rep.mPoints.mWorstIndex, 7);
}

TEST(Diff, DifferentCellCount) {
    Mesh a = tri_with_points(kSquare);
    Mesh b = mt::make_mesh(kSquare, "triangle", {{0, 1, 2}});
    DiffReport rep = diff(a, b);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::Different);
    ASSERT_FALSE(rep.mBlocks.empty());
    EXPECT_TRUE(rep.mBlocks[0].mCountMismatch);
    EXPECT_EQ(rep.mBlocks[0].mCountA, 2);
    EXPECT_EQ(rep.mBlocks[0].mCountB, 1);
}

TEST(Diff, DifferentCellType) {
    Mesh a = mt::make_mesh(kSquare, "quad", {{0, 1, 2, 3}});
    Mesh b = tri_with_points(kSquare);
    DiffReport rep = diff(a, b);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::Different);
    ASSERT_FALSE(rep.mBlocks.empty());
    EXPECT_TRUE(rep.mBlocks[0].mTypeMismatch);
    EXPECT_EQ(rep.mBlocks[0].mTypeA, "quad");
    EXPECT_EQ(rep.mBlocks[0].mTypeB, "triangle");
}

TEST(Diff, DifferentConnectivity) {
    Mesh a = tri_with_points(kSquare);
    Mesh b = mt::make_mesh(kSquare, "triangle", {{0, 1, 2}, {1, 2, 3}});
    DiffReport rep = diff(a, b);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::Different);
    ASSERT_FALSE(rep.mBlocks.empty());
    EXPECT_EQ(rep.mBlocks[0].mConnMismatchCount, 1);
    ASSERT_FALSE(rep.mBlocks[0].mFirstMismatches.empty());
    EXPECT_EQ(rep.mBlocks[0].mFirstMismatches[0], 1);
}

TEST(Diff, PointDataKeyOnlyInOne) {
    Mesh a = tri_with_points(kSquare);
    Mesh b = tri_with_points(kSquare);
    NDArray d = NDArray::Uninit(meshioplusplus::DType::Float64, {4, 1});
    for (int i = 0; i < 4; ++i)
        d.As<double>()[i] = static_cast<double>(i);
    a.AddPointData("temperature", std::move(d));
    DiffReport rep = diff(a, b);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::Different);
    ASSERT_EQ(rep.mPointData.mOnlyInA.size(), 1u);
    EXPECT_EQ(rep.mPointData.mOnlyInA[0], "temperature");
    EXPECT_TRUE(rep.mPointData.mOnlyInB.empty());
}

TEST(Diff, PointDataValuesWithinTolerance) {
    Mesh a = tri_with_points(kSquare);
    Mesh b = tri_with_points(kSquare);
    NDArray da = NDArray::Uninit(meshioplusplus::DType::Float64, {4, 1});
    NDArray db = NDArray::Uninit(meshioplusplus::DType::Float64, {4, 1});
    for (int i = 0; i < 4; ++i) {
        da.As<double>()[i] = static_cast<double>(i);
        db.As<double>()[i] = static_cast<double>(i) + 1e-10;
    }
    a.AddPointData("t", std::move(da));
    b.AddPointData("t", std::move(db));
    DiffOptions opts;
    opts.atol = 1e-8;
    opts.rtol = 0.0;
    DiffReport rep = diff(a, b, opts);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::EqualWithinTolerance);
    ASSERT_EQ(rep.mPointData.mShared.size(), 1u);
    EXPECT_EQ(rep.mPointData.mShared[0].mName, "t");
    EXPECT_GT(rep.mPointData.mShared[0].mMaxAbsError, 0.0);
}

TEST(Diff, UnorderedShuffledNodes) {
    // Same geometry + connectivity as `a`, but node order permuted.
    Mesh a = tri_with_points(kSquare);
    // Permutation old->new: 0->2, 1->0, 2->3, 3->1
    std::vector<std::vector<double>> pts_b = {kSquare[1], kSquare[3], kSquare[0], kSquare[2]};
    // Cells in a reference points {0,1,2} and {0,2,3}; remap through inverse.
    // new index of old i: inv[0]=2, inv[1]=0, inv[2]=3, inv[3]=1
    Mesh b = mt::make_mesh(pts_b, "triangle", {{2, 0, 3}, {2, 3, 1}});
    DiffOptions opts;
    opts.unordered = true;
    DiffReport rep = diff(a, b, opts);
    EXPECT_TRUE(rep.mUnordered);
    EXPECT_FALSE(rep.mCorrespondenceFailed);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::Identical);
    ASSERT_FALSE(rep.mBlocks.empty());
    EXPECT_EQ(rep.mBlocks[0].mConnMismatchCount, 0);
}

TEST(Diff, UnorderedCorrespondenceFails) {
    Mesh a = tri_with_points(kSquare);
    // Move every point far away -> no proximity match.
    std::vector<std::vector<double>> far = {{10, 10, 10}, {11, 10, 10}, {11, 11, 10}, {10, 11, 10}};
    Mesh b = tri_with_points(far);
    DiffOptions opts;
    opts.unordered = true;
    opts.atol = 1e-6;
    opts.rtol = 1e-6;
    DiffReport rep = diff(a, b, opts);
    EXPECT_TRUE(rep.mCorrespondenceFailed);
    EXPECT_EQ(rep.mVerdict, DiffVerdict::Different);
}

}  // namespace
