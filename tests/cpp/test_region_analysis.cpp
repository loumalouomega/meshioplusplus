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

// hausdorff_distance, edit_regions and match_periodic_nodes (roadmap §5,
// v16.23.0).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/periodic.hpp"
#include "meshioplusplus/operations/region_ops.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/region.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

using mt::make_mesh;

// An n x n grid of unit quads in z = 0.
Mesh quad_grid(int N, double Z = 0.0) {
    std::vector<std::vector<double>> pts;
    for (int j = 0; j <= N; ++j)
        for (int i = 0; i <= N; ++i)
            pts.push_back({double(i), double(j), Z});
    std::vector<std::vector<std::int64_t>> quads;
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            const std::int64_t p = j * (N + 1) + i;
            quads.push_back({p, p + 1, p + N + 2, p + N + 1});
        }
    return make_mesh(pts, "quad", quads);
}

Region region(const std::string& rName, RegionKind Kind, std::vector<std::int64_t> Entries,
              int Dim = -1, std::int64_t Tag = -1) {
    Region r;
    r.mName = rName;
    r.mKind = Kind;
    r.mDim = Dim;
    r.mTag = Tag;
    const bool side = Kind == RegionKind::Side;
    r.mEntries = side ? NDArray::Uninit(DType::Int64, {Entries.size() / 2, std::size_t{2}})
                      : NDArray::Uninit(DType::Int64, {Entries.size()});
    for (std::size_t i = 0; i < Entries.size(); ++i)
        r.mEntries.As<std::int64_t>()[i] = Entries[i];
    return r;
}

std::vector<std::int64_t> entries(const Mesh& rMesh, const std::string& rName) {
    RegionSelector s;
    s.mName = rName;
    const Region& r = rMesh.Region(find_region(rMesh, s));
    return std::vector<std::int64_t>(r.mEntries.As<std::int64_t>(),
                                     r.mEntries.As<std::int64_t>() + r.mEntries.Size());
}

RegionEdit edit(RegionOp Op, std::vector<std::string> Inputs, std::string Output = "") {
    RegionEdit e;
    e.mOp = Op;
    for (auto& n : Inputs) {
        RegionSelector s;
        s.mName = n;
        e.mInputs.push_back(s);
    }
    e.mOutputName = std::move(Output);
    return e;
}

}  // namespace

// ---------------------------------------------------------------- hausdorff

TEST(Hausdorff, AMeshIsAtDistanceZeroFromItself) {
    const HausdorffResult r = hausdorff_distance(quad_grid(3), quad_grid(3));
    EXPECT_EQ(r.mDistance, 0.0);
    EXPECT_EQ(r.mNumSamplesA, 16);
}

TEST(Hausdorff, ATranslatedCopyIsAtTheTranslation) {
    const HausdorffResult r = hausdorff_distance(quad_grid(3), quad_grid(3, 0.25));
    EXPECT_NEAR(r.mDistance, 0.25, 1e-12);
    EXPECT_NEAR(r.mAtoB, 0.25, 1e-12);
    EXPECT_NEAR(r.mBtoA, 0.25, 1e-12);
    EXPECT_NEAR(r.mMeanAtoB, 0.25, 1e-12);
    EXPECT_NEAR(r.mRmsBtoA, 0.25, 1e-12);
}

TEST(Hausdorff, IsAsymmetricPerSideAndSymmetricOverall) {
    // B covers only half of A: A's far half is up to 1.5 away from B, but all
    // of B lies on A.
    const Mesh a = quad_grid(3);
    const Mesh b =
        make_mesh({{0, 0, 0}, {1.5, 0, 0}, {1.5, 3, 0}, {0, 3, 0}}, "quad", {{0, 1, 2, 3}});
    const HausdorffResult ab = hausdorff_distance(a, b);
    const HausdorffResult ba = hausdorff_distance(b, a);
    EXPECT_NEAR(ab.mAtoB, 1.5, 1e-12);
    EXPECT_EQ(ab.mBtoA, 0.0);
    EXPECT_EQ(ab.mDistance, ba.mDistance);
    EXPECT_EQ(ab.mAtoB, ba.mBtoA);
    EXPECT_NEAR(ab.mWorstPointA[0], 3.0, 1e-12);
}

TEST(Hausdorff, FaceSamplesSeeWhatVerticesMiss) {
    // A flat square against a tent over it: every vertex of the square lies on
    // the tent's rim, but its centre is 1 below the apex.
    const Mesh flat =
        make_mesh({{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}}, "quad", {{0, 1, 2, 3}});
    const Mesh tent = make_mesh({{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {1, 1, 1}}, "triangle",
                                {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}});
    const HausdorffResult coarse = hausdorff_distance(flat, tent);
    HausdorffOptions o;
    o.mFaceSamples = 4;
    const HausdorffResult fine = hausdorff_distance(flat, tent, o);
    EXPECT_EQ(coarse.mAtoB, 0.0);
    EXPECT_GT(fine.mAtoB, 0.3);
    EXPECT_EQ(fine.mNumSamplesA, 4 + 2 * 16);
}

TEST(Hausdorff, VolumeMeshesCompareTheirSkin) {
    const HausdorffResult r = hausdorff_distance(mt::hex_mesh(), mt::hex_mesh());
    EXPECT_EQ(r.mDistance, 0.0);
    HausdorffOptions o;
    o.mRegionA = "any";
    EXPECT_THROW(hausdorff_distance(mt::hex_mesh(), mt::hex_mesh(), o), std::invalid_argument);
}

TEST(Hausdorff, RefusesAMeshWithoutTriangles) {
    EXPECT_THROW(hausdorff_distance(mt::line_mesh(), quad_grid(1)), std::invalid_argument);
}

// ------------------------------------------------------------- edit_regions

TEST(RegionOps, SetAlgebraOnCellRegions) {
    Mesh m = quad_grid(3);
    m.AddRegion(region("a", RegionKind::Cell, {0, 1, 2, 3}));
    m.AddRegion(region("b", RegionKind::Cell, {2, 3, 4, 5}));
    const Mesh out = edit_regions(
        m, {edit(RegionOp::Union, {"a", "b"}, "u"), edit(RegionOp::Intersection, {"a", "b"}, "i"),
            edit(RegionOp::Difference, {"a", "b"}, "d")});
    EXPECT_EQ(entries(out, "u"), (std::vector<std::int64_t>{0, 1, 2, 3, 4, 5}));
    EXPECT_EQ(entries(out, "i"), (std::vector<std::int64_t>{2, 3}));
    EXPECT_EQ(entries(out, "d"), (std::vector<std::int64_t>{0, 1}));
    EXPECT_EQ(out.NumRegions(), 5u);
    EXPECT_EQ(m.NumRegions(), 2u);  // the input is untouched
}

TEST(RegionOps, SideRegionsCombinePairs) {
    Mesh m = quad_grid(2);
    m.AddRegion(region("s1", RegionKind::Side, {0, 0, 1, 0}));
    m.AddRegion(region("s2", RegionKind::Side, {1, 0, 1, 1}));
    const Mesh out = edit_regions(m, {edit(RegionOp::Union, {"s1", "s2"}, "s")});
    EXPECT_EQ(entries(out, "s"), (std::vector<std::int64_t>{0, 0, 1, 0, 1, 1}));
}

TEST(RegionOps, EmptyResultIsKept) {
    Mesh m = quad_grid(2);
    m.AddRegion(region("a", RegionKind::Cell, {0}));
    m.AddRegion(region("b", RegionKind::Cell, {1}));
    const Mesh out = edit_regions(m, {edit(RegionOp::Intersection, {"a", "b"}, "none")});
    EXPECT_TRUE(entries(out, "none").empty());
}

TEST(RegionOps, KeepInputsOffRemovesThem) {
    Mesh m = quad_grid(2);
    m.AddRegion(region("a", RegionKind::Cell, {0}));
    m.AddRegion(region("b", RegionKind::Cell, {1}));
    RegionEdit e = edit(RegionOp::Union, {"a", "b"}, "ab");
    e.mKeepInputs = false;
    const Mesh out = edit_regions(m, {e});
    EXPECT_EQ(out.NumRegions(), 1u);
    EXPECT_EQ(out.Region(0).mName, "ab");
}

TEST(RegionOps, RenameAndRetag) {
    Mesh m = quad_grid(2);
    m.AddRegion(region("Surface 3", RegionKind::Cell, {0, 1}, 2, 3));
    RegionEdit rt = edit(RegionOp::Retag, {"Outlet"});
    rt.mOutputTag = 10;
    const Mesh out = edit_regions(m, {edit(RegionOp::Rename, {"Surface 3"}, "Outlet"), rt});
    ASSERT_EQ(out.NumRegions(), 1u);
    EXPECT_EQ(out.Region(0).mName, "Outlet");
    EXPECT_EQ(out.Region(0).mTag, 10);
    EXPECT_EQ(out.Region(0).mDim, 2);
    EXPECT_EQ(entries(out, "Outlet"), (std::vector<std::int64_t>{0, 1}));
}

TEST(RegionOps, AmbiguousAndMissingSelectorsThrow) {
    Mesh m = quad_grid(2);
    m.AddRegion(region("wall", RegionKind::Cell, {0}, 2, 1));
    m.AddRegion(region("wall", RegionKind::Point, {0}, 0, 1));
    EXPECT_THROW(edit_regions(m, {edit(RegionOp::Delete, {"wall"})}), std::invalid_argument);
    EXPECT_THROW(edit_regions(m, {edit(RegionOp::Delete, {"nope"})}), std::invalid_argument);
    RegionEdit d = edit(RegionOp::Delete, {"wall"});
    d.mInputs[0].mKind = static_cast<std::int32_t>(RegionKind::Point);
    const Mesh out = edit_regions(m, {d});
    ASSERT_EQ(out.NumRegions(), 1u);
    EXPECT_EQ(out.Region(0).mKind, RegionKind::Cell);
}

TEST(RegionOps, MixedKindsAndCollisionsThrow) {
    Mesh m = quad_grid(2);
    m.AddRegion(region("a", RegionKind::Cell, {0}));
    m.AddRegion(region("p", RegionKind::Point, {0}));
    m.AddRegion(region("c", RegionKind::Cell, {1}));
    EXPECT_THROW(edit_regions(m, {edit(RegionOp::Union, {"a", "p"}, "x")}), std::invalid_argument);
    EXPECT_THROW(edit_regions(m, {edit(RegionOp::Rename, {"a"}, "c")}), std::invalid_argument);
    EXPECT_THROW(edit_regions(m, {edit(RegionOp::Union, {"a"}, "x")}), std::invalid_argument);
    EXPECT_THROW(region_op_from_name("xor"), std::invalid_argument);
    EXPECT_EQ(region_op_from_name("intersect"), RegionOp::Intersection);
}

TEST(RegionOps, RemoveRegionKeepsTheOthersInOrder) {
    Mesh m = quad_grid(2);
    m.AddRegion(region("a", RegionKind::Cell, {0}));
    m.AddRegion(region("b", RegionKind::Cell, {1}));
    m.AddRegion(region("c", RegionKind::Cell, {2}));
    m.RemoveRegion(1);
    ASSERT_EQ(m.NumRegions(), 2u);
    EXPECT_EQ(m.Region(0).mName, "a");
    EXPECT_EQ(m.Region(1).mName, "c");
    m.RemoveRegion(7);  // out of range: a no-op
    EXPECT_EQ(m.NumRegions(), 2u);
}

// ------------------------------------------------------------------ periodic

TEST(Periodic, TranslationPairsOppositeEdges) {
    Mesh m = quad_grid(3);
    // Left (x = 0) and right (x = 3) columns of nodes, as Point regions.
    m.AddRegion(region("left", RegionKind::Point, {0, 4, 8, 12}));
    m.AddRegion(region("right", RegionKind::Point, {3, 7, 11, 15}));
    PeriodicOptions o;
    o.mTransform = transform_translation(3.0, 0.0, 0.0);
    RegionSelector s, t;
    s.mName = "left";
    t.mName = "right";
    const PeriodicPairs p = match_periodic_nodes(m, s, t, o);
    ASSERT_EQ(p.mSlave.Size(), 4u);
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_EQ(p.mSlave.As<std::int64_t>()[i], static_cast<std::int64_t>(i * 4));
        EXPECT_EQ(p.mMaster.As<std::int64_t>()[i], static_cast<std::int64_t>(i * 4 + 3));
    }
    EXPECT_EQ(p.mMaxResidual, 0.0);
    EXPECT_EQ(p.mNumFixed, 0);
}

TEST(Periodic, CellAndSideRegionsContributeTheirNodes) {
    Mesh m = quad_grid(3);
    // Cells 0, 3, 6 are the left column; side (c, 3) of a quad is its edge 3-0,
    // the x = 0 edge.
    m.AddRegion(region("left", RegionKind::Side, {0, 3, 3, 3, 6, 3}));
    m.AddRegion(region("right", RegionKind::Cell, {2, 5, 8}));
    PeriodicOptions o;
    o.mTransform = transform_translation(3.0, 0.0, 0.0);
    o.mRequireComplete = false;
    RegionSelector s, t;
    s.mName = "left";
    t.mName = "right";
    const PeriodicPairs p = match_periodic_nodes(m, s, t, o);
    EXPECT_EQ(p.mSlave.Size(), 4u);
    EXPECT_EQ(p.mUnmatched.Size(), 0u);
}

TEST(Periodic, RotationFindsTheAxisAsAFixedPoint) {
    // A quarter disk: nodes on the x axis map onto the y axis under +90 deg
    // about z; the origin maps onto itself.
    Mesh m = make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {0, 2, 0}, {1, 1, 0}},
                       "triangle", {{0, 1, 5}, {0, 5, 3}, {1, 2, 5}, {3, 5, 4}});
    m.AddRegion(region("x", RegionKind::Point, {0, 1, 2}));
    m.AddRegion(region("y", RegionKind::Point, {0, 3, 4}));
    PeriodicOptions o;
    o.mTransform = transform_rotation(0.0, 0.0, 1.0, 3.14159265358979323846 / 2.0);
    o.mAtol = 1e-9;
    RegionSelector s, t;
    s.mName = "x";
    t.mName = "y";
    const PeriodicPairs p = match_periodic_nodes(m, s, t, o);
    ASSERT_EQ(p.mSlave.Size(), 2u);
    EXPECT_EQ(p.mSlave.As<std::int64_t>()[0], 1);
    EXPECT_EQ(p.mMaster.As<std::int64_t>()[0], 3);
    EXPECT_EQ(p.mSlave.As<std::int64_t>()[1], 2);
    EXPECT_EQ(p.mMaster.As<std::int64_t>()[1], 4);
    EXPECT_EQ(p.mNumFixed, 1);
    EXPECT_LT(p.mMaxResidual, 1e-15);
}

TEST(Periodic, IncompleteAndDoubleClaimsThrow) {
    Mesh m = quad_grid(3);
    m.AddRegion(region("left", RegionKind::Point, {0, 4, 8, 12}));
    m.AddRegion(region("right", RegionKind::Point, {3, 7, 11}));
    PeriodicOptions o;
    o.mTransform = transform_translation(3.0, 0.0, 0.0);
    RegionSelector s, t;
    s.mName = "left";
    t.mName = "right";
    EXPECT_THROW(match_periodic_nodes(m, s, t, o), std::invalid_argument);
    o.mRequireComplete = false;
    const PeriodicPairs p = match_periodic_nodes(m, s, t, o);
    ASSERT_EQ(p.mUnmatched.Size(), 1u);
    EXPECT_EQ(p.mUnmatched.As<std::int64_t>()[0], 12);
    // A tolerance wide enough to reach two rows claims masters twice.
    o.mAtol = 1.5;
    o.mTransform = transform_translation(3.0, 0.5, 0.0);
    EXPECT_THROW(match_periodic_nodes(m, s, t, o), std::invalid_argument);
    o.mAtol = 0.0;
    EXPECT_THROW(match_periodic_nodes(m, s, t, o), std::invalid_argument);
}
