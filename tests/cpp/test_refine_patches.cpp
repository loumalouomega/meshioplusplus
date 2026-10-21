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
// refine_patches: variational shape approximation of a segmentation.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/refine_patches.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "analytic_surfaces.hpp"
#include "fit_cad_fixtures.hpp"
#include "golden_fixtures.hpp"

using namespace meshioplusplus;

namespace {

Mesh segmented(const Mesh& rSkin, double FeatureAngle = 30.0) {
    SegmentOptions o;
    o.mFeatureAngleDeg = FeatureAngle;
    return segment(rSkin, o).mMesh;
}

// The patch label of every cell, block-major.
std::vector<std::int64_t> labels(const Mesh& rM) {
    std::vector<std::int64_t> out;
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
        const NDArray& a = rM.CellData(kSegmentPatchName, b);
        for (std::size_t i = 0; i < a.Size(); ++i)
            out.push_back(a.As<std::int64_t>()[i]);
    }
    return out;
}

// The rings of every cell of a single-block surface, as point ids.
std::vector<std::vector<std::int64_t>> rings(const Mesh& rM) {
    std::vector<std::vector<std::int64_t>> out;
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
        const NDArray& c = rM.Cells(b).Conn();
        const std::size_t n = c.Shape()[1];
        for (std::size_t i = 0; i < c.Shape()[0]; ++i)
            out.emplace_back(c.As<std::int64_t>() + i * n, c.As<std::int64_t>() + (i + 1) * n);
    }
    return out;
}

// Every patch's cells are connected through shared edges.
bool patches_connected(const Mesh& rM) {
    const std::vector<std::int64_t> lab = labels(rM);
    const std::vector<std::vector<std::int64_t>> r = rings(rM);
    std::map<std::pair<std::int64_t, std::int64_t>, std::vector<std::size_t>> edges;
    for (std::size_t f = 0; f < r.size(); ++f)
        for (std::size_t k = 0; k < r[f].size(); ++k) {
            std::int64_t a = r[f][k], b = r[f][(k + 1) % r[f].size()];
            if (a > b)
                std::swap(a, b);
            edges[{a, b}].push_back(f);
        }
    std::vector<std::vector<std::size_t>> nb(r.size());
    for (const auto& kv : edges)
        if (kv.second.size() == 2) {
            nb[kv.second[0]].push_back(kv.second[1]);
            nb[kv.second[1]].push_back(kv.second[0]);
        }
    std::map<std::int64_t, std::size_t> count;
    for (std::int64_t l : lab)
        if (l >= 0)
            ++count[l];
    std::set<std::int64_t> done;
    for (std::size_t f = 0; f < r.size(); ++f) {
        if (lab[f] < 0 || done.count(lab[f]))
            continue;
        done.insert(lab[f]);
        std::vector<std::uint8_t> seen(r.size(), 0);
        std::vector<std::size_t> stack{f};
        seen[f] = 1;
        std::size_t reached = 0;
        while (!stack.empty()) {
            const std::size_t g = stack.back();
            stack.pop_back();
            ++reached;
            for (std::size_t h : nb[g])
                if (!seen[h] && lab[h] == lab[f]) {
                    seen[h] = 1;
                    stack.push_back(h);
                }
        }
        if (reached != count[lab[f]])
            return false;
    }
    return true;
}

std::int64_t count_proxy(const Mesh& rM, RefineProxy Kind) {
    const NDArray& p = rM.FieldData(kRefinePatchProxyName);
    std::int64_t n = 0;
    for (std::size_t i = 0; i < p.Size(); ++i)
        n += p.As<std::int64_t>()[i] == static_cast<std::int64_t>(Kind) ? 1 : 0;
    return n;
}

}  // namespace

TEST(RefinePatches, UnchangedPartitionIsByteIdentical) {
    // A box segmented by its creases is already six planes: nothing moves.
    const Mesh seg = segmented(analytic::box_skin(3));
    const RefinePatchesResult r = refine_patches(seg);
    EXPECT_EQ(r.mStop, RefineStop::Unchanged);
    EXPECT_EQ(labels(r.mMesh), labels(seg));
    EXPECT_EQ(r.mNumPatches, 6);
    EXPECT_EQ(r.mNumPrimitive, 6);
    EXPECT_EQ(r.mNumSplits + r.mNumMerges + r.mNumAbsorbed, 0);
    EXPECT_NEAR(r.mExplainedAreaFraction, 1.0, 1e-12);
    EXPECT_LE(r.mHausdorffMax, r.mTolerance);
    // The adjacency table and the regions are rewritten the same way.
    const NDArray& adj = r.mMesh.FieldData(kSegmentAdjacencyName);
    const NDArray& adj0 = seg.FieldData(kSegmentAdjacencyName);
    ASSERT_EQ(adj.Size(), adj0.Size());
    for (std::size_t i = 0; i < adj.Size(); ++i)
        EXPECT_EQ(adj.As<std::int64_t>()[i], adj0.As<std::int64_t>()[i]);
    EXPECT_EQ(r.mMesh.NumRegions(), seg.NumRegions());
    // The per-patch report.
    EXPECT_EQ(r.mMesh.FieldData(kRefinePatchIdName).Size(), 6u);
    EXPECT_EQ(count_proxy(r.mMesh, RefineProxy::Plane), 6);
    EXPECT_EQ(r.mMesh.FieldData(kRefineExplainedAreaName).As<double>()[0], 1.0);
}

TEST(RefinePatches, SplitsAPatchHoldingSeveralSurfaces) {
    // With a feature angle above the box's 90 degrees every edge is smooth and
    // the whole box is one patch no primitive explains: RANSAC finds the
    // planes inside it, one or more per iteration.
    const Mesh seg = segmented(analytic::box_skin(3), 100.0);
    ASSERT_EQ(fit_primitives(seg).mNumPatches, 1);
    RefinePatchesOptions o;
    o.mFeatureAngleDeg = 100.0;
    o.mMaxIterations = 20;
    const RefinePatchesResult r = refine_patches(seg, o);
    EXPECT_GT(r.mNumSplits, 0);
    EXPECT_GT(r.mExplainedAreaFraction, r.mInitialExplainedAreaFraction);
    EXPECT_EQ(r.mInitialExplainedAreaFraction, 0.0);
    // Every side comes back: six planes explaining the whole box.
    EXPECT_EQ(r.mNumPrimitive, 6);
    EXPECT_EQ(r.mNumPatches, 6);
    EXPECT_NEAR(r.mExplainedAreaFraction, 1.0, 1e-12);
    EXPECT_EQ(count_proxy(r.mMesh, RefineProxy::Plane), 6);
    EXPECT_TRUE(patches_connected(r.mMesh));
    // Without splitting nothing can be found.
    o.mAllowSplit = false;
    EXPECT_EQ(refine_patches(seg, o).mNumPrimitive, 0);
}

// A CAD skin segmented by angles: tangent junctions leave patches holding
// two surfaces (a wall and its fillet, a plate and the fillet into it), which
// refinement separates until primitives explain all of it, noise included.
TEST(RefinePatches, ExplainsTangentJunctions) {
    struct Case {
        const char* mName;
        Mesh mSkin;
        double mMinExplained;
        double mHausdorffSlack;  // primitives hold 9 points in 10 within the tolerance
    };
    std::vector<Case> cases;
    cases.push_back({"turned_part", analytic::turned_part(), 1.0, 1.0});
    cases.push_back({"filleted_plate", analytic::filleted_plate(), 1.0, 1.0});
    cases.push_back(
        {"noisy_turned_part", fcfix::jittered(analytic::turned_part(), 5e-4), 0.99, 2.0});
    for (const Case& c : cases) {
        SCOPED_TRACE(c.mName);
        const Mesh seg = segmented(c.mSkin);
        const RefinePatchesResult r = refine_patches(seg);
        EXPECT_LT(r.mInitialExplainedAreaFraction, 0.95);
        EXPECT_GE(r.mExplainedAreaFraction, c.mMinExplained - 1e-12);
        EXPECT_GT(r.mNumSplits, 0);
        EXPECT_LE(r.mHausdorffMax, c.mHausdorffSlack * r.mTolerance * (1.0 + 1e-12));
        EXPECT_TRUE(patches_connected(r.mMesh));
    }
}

TEST(RefinePatches, MergesAnOverSegmentedCylinder) {
    // The capped cylinder's wall cut in two along x: one cylinder explains both.
    Mesh seg = segmented(analytic::capped_cylinder());
    const std::int64_t wall = [&] {
        // The wall is the patch with the most cells.
        std::map<std::int64_t, std::int64_t> n;
        for (std::int64_t l : labels(seg))
            ++n[l];
        std::int64_t best = 0, most = -1;
        for (const auto& kv : n)
            if (kv.second > most) {
                best = kv.first;
                most = kv.second;
            }
        return best;
    }();
    std::int64_t next = 0;
    for (std::int64_t l : labels(seg))
        next = std::max(next, l + 1);
    const NDArray& pts = seg.Points();
    for (std::size_t b = 0; b < seg.NumCellBlocks(); ++b) {
        NDArray lab = NDArray::Uninit(DType::Int64, {seg.Cells(b).NumCells()});
        const NDArray& old = seg.CellData(kSegmentPatchName, b);
        const NDArray& c = seg.Cells(b).Conn();
        const std::size_t k = c.Shape()[1];
        for (std::size_t i = 0; i < lab.Size(); ++i) {
            std::int64_t l = old.As<std::int64_t>()[i];
            double y = 0.0;
            for (std::size_t j = 0; j < k; ++j)
                y +=
                    pts.As<double>()[static_cast<std::size_t>(c.As<std::int64_t>()[i * k + j]) * 3 +
                                     1];
            if (l == wall && y > 0.0)
                l = next;
            lab.As<std::int64_t>()[i] = l;
        }
        std::vector<NDArray> blocks;
        for (std::size_t bb = 0; bb < seg.NumCellBlocks(); ++bb)
            blocks.push_back(bb == b ? lab : seg.CellData(kSegmentPatchName, bb));
        seg.AddCellData(kSegmentPatchName, std::move(blocks));
    }
    const RefinePatchesResult r = refine_patches(seg);
    EXPECT_GE(r.mNumMerges, 1);
    EXPECT_EQ(r.mNumPatches, 3);
    EXPECT_EQ(count_proxy(r.mMesh, RefineProxy::Cylinder), 1);
    EXPECT_EQ(count_proxy(r.mMesh, RefineProxy::Plane), 2);
    RefinePatchesOptions off;
    off.mAllowMerge = false;
    EXPECT_EQ(refine_patches(seg, off).mNumMerges, 0);
}

TEST(RefinePatches, NeverCrossesCreasesAndStaysConnected) {
    const Mesh seg = segmented(fcfix::saddle_box(12));
    const RefinePatchesResult r = refine_patches(seg);
    EXPECT_TRUE(patches_connected(r.mMesh));
    // Every patch lies on one side of the box: a crease is never crossed, so
    // its cells' outward normals agree to within the 30 degree feature angle
    // of each other across the whole patch (the saddle top bends less).
    const std::vector<std::int64_t> lab = labels(r.mMesh);
    const std::vector<std::vector<std::int64_t>> rg = rings(r.mMesh);
    const NDArray& p = r.mMesh.Points();
    std::map<std::int64_t, std::vector<std::array<double, 3>>> normals;
    for (std::size_t f = 0; f < rg.size(); ++f) {
        const double* a = p.As<double>() + rg[f][0] * 3;
        const double* b = p.As<double>() + rg[f][1] * 3;
        const double* c = p.As<double>() + rg[f][2] * 3;
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        std::array<double, 3> n = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                   u[0] * v[1] - u[1] * v[0]};
        const double l = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        for (double& x : n)
            x /= l;
        normals[lab[f]].push_back(n);
    }
    for (const auto& kv : normals)
        for (const auto& n : kv.second)
            EXPECT_GT(n[0] * kv.second[0][0] + n[1] * kv.second[0][1] + n[2] * kv.second[0][2], 0.5)
                << "patch " << kv.first;
}

TEST(RefinePatches, EnergyNeverIncreases) {
    // The Lloyd steps and merges never raise the energy; a split trades the
    // freeform rest's cost for the cells a primitive takes, and absorbing a
    // small patch is not an energy move, so both are off here.
    RefinePatchesOptions o;
    o.mMinPatchCells = 1;
    o.mAllowSplit = false;
    for (const Mesh& skin : {analytic::turned_part(), fcfix::saddle_box(16),
                             fcfix::jittered(analytic::turned_part(), 1e-3)}) {
        const RefinePatchesResult r = refine_patches(segmented(skin), o);
        EXPECT_GE(r.mEnergy, 0.0);
        EXPECT_LE(r.mEnergy, r.mInitialEnergy + 1e-9 * (r.mInitialEnergy + 1.0));
        EXPECT_GE(r.mExplainedAreaFraction, 0.0);
        EXPECT_LE(r.mExplainedAreaFraction, 1.0);
        EXPECT_LE(r.mNumIterations, o.mMaxIterations);
    }
}

TEST(RefinePatches, IsDeterministicForASeed) {
    const Mesh seg = segmented(analytic::box_skin(3), 100.0);
    RefinePatchesOptions o;
    o.mFeatureAngleDeg = 100.0;
    const std::uint64_t a = golden::ft_mesh(refine_patches(seg, o).mMesh);
    const std::uint64_t b = golden::ft_mesh(refine_patches(seg, o).mMesh);
    EXPECT_EQ(a, b);
}

TEST(RefinePatches, FitTableIsFitPrimitivesAtTheTolerance) {
    const RefinePatchesResult r = refine_patches(segmented(analytic::turned_part()));
    FitPrimitivesOptions fo;
    fo.mTolerance = r.mTolerance;
    fo.mRelative = false;
    const Mesh again = fit_primitives(r.mMesh, fo).mMesh;
    for (const char* name : {"fit:patch_id", "fit:patch_type", "fit:patch_params"}) {
        const NDArray& x = r.mMesh.FieldData(name);
        const NDArray& y = again.FieldData(name);
        ASSERT_EQ(x.Nbytes(), y.Nbytes()) << name;
        EXPECT_EQ(std::memcmp(x.Data(), y.Data(), x.Nbytes()), 0) << name;
    }
}

TEST(RefinePatches, RefusesBadInput) {
    EXPECT_THROW(refine_patches(analytic::box_skin(2)), std::invalid_argument);
    const Mesh seg = segmented(analytic::box_skin(2));
    const auto bad = [&](auto&& rSet) {
        RefinePatchesOptions o;
        rSet(o);
        EXPECT_THROW(refine_patches(seg, o), std::invalid_argument);
    };
    bad([](RefinePatchesOptions& o) { o.mTolerance = 0.0; });
    bad([](RefinePatchesOptions& o) { o.mNormalWeight = -1.0; });
    bad([](RefinePatchesOptions& o) { o.mEpsilon = -1.0; });
    bad([](RefinePatchesOptions& o) { o.mFeatureAngleDeg = 200.0; });
    bad([](RefinePatchesOptions& o) { o.mMaxIterations = -1; });
    bad([](RefinePatchesOptions& o) { o.mMinPatchCells = 0; });
    bad([](RefinePatchesOptions& o) { o.mTypes = 1u; });
    bad([](RefinePatchesOptions& o) { o.mPrefix.clear(); });
}

TEST(RefinePatches, ZeroIterationsRewritesTheInputPartition) {
    const Mesh seg = segmented(analytic::capped_cylinder());
    RefinePatchesOptions o;
    o.mMaxIterations = 0;
    const RefinePatchesResult r = refine_patches(seg, o);
    EXPECT_EQ(r.mStop, RefineStop::IterationCap);
    EXPECT_EQ(labels(r.mMesh), labels(seg));
    EXPECT_EQ(r.mNumIterations, 0);
    EXPECT_TRUE(r.mMesh.HasFieldData("fit:patch_type"));
}
