// face_support and fit_cad_iterate (roadmap 7.7.1, 7.7.2).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/fit_cad_loop.hpp"
#include "meshioplusplus/operations/repair.hpp"
#include "analytic_surfaces.hpp"
#include "fit_cad_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// A box scanned without its bottom, the hole closed by repair (the fill
// marked repair:hole).
Mesh fcl_open_box() {
    const Mesh box = analytic::box_skin(4);
    std::vector<std::int64_t> conn;
    const auto cb = box.Cells(0);
    const double* p = box.Points().As<double>();
    for (std::size_t c = 0; c < cb.NumCells(); ++c) {
        const std::int64_t* q = cb.Conn().As<std::int64_t>() + 4 * c;
        bool bottom = true;
        for (int k = 0; k < 4; ++k)
            bottom = bottom && p[3 * q[k] + 2] == 0.0;
        if (!bottom)
            conn.insert(conn.end(), q, q + 4);
    }
    Mesh open;
    open.AssignPoints(box.Points());
    NDArray a = NDArray::Uninit(DType::Int64, {conn.size() / 4, 4});
    std::copy(conn.begin(), conn.end(), a.As<std::int64_t>());
    open.AddCellBlock("quad", std::move(a));
    RepairOptions o;
    o.mMaxHoleEdges = 0;
    o.mRecordProvenance = true;
    return repair(open, o).mMesh;
}

// A box whose top bends by `Bend` along its middle: two planes a coarse
// segmentation merges and a finer one separates.
Mesh fcl_bent_box(double Bend) {
    Mesh box = analytic::box_skin(8);
    NDArray q = NDArray::Uninit(DType::Float64, {box.NumPoints(), 3});
    const double* p = box.Points().As<double>();
    for (std::size_t i = 0; i < box.NumPoints(); ++i) {
        q.As<double>()[3 * i] = p[3 * i];
        q.As<double>()[3 * i + 1] = p[3 * i + 1];
        q.As<double>()[3 * i + 2] = p[3 * i + 2] * (1.0 + Bend * std::max(0.0, p[3 * i] - 0.5));
    }
    box.AssignPoints(std::move(q));
    return box;
}

}  // namespace

TEST(FaceSupport, TheBottomAScanMissedIsInferred) {
    const Mesh scan = fcl_open_box();
    ASSERT_TRUE(scan.HasCellData("repair:hole"));
    FitCadOptions fo;
    fo.mBackend = FitCadBackend::Native;
    const FitCadResult fit = fit_cad(fcfix::fitted(scan), "", fo);
    const FaceSupportResult r = face_support(fit.mMesh, scan);
    EXPECT_EQ(r.mNumMeasured, 5);
    EXPECT_EQ(r.mNumInferred, 1);
    EXPECT_EQ(r.mNumPartial, 0);
    EXPECT_TRUE(r.mMesh.HasFieldData(kFitCadFaceSupportName));
    EXPECT_TRUE(r.mMesh.HasCellData(kFitCadSupportName));
    EXPECT_NEAR(r.mShareWithin, 1.0, 1e-12);
    // Without the fill marks, every face is measured.
    Mesh unmarked = detail::clone_mesh(scan);
    FaceSupportOptions so;
    so.mHoleData = "none";
    EXPECT_EQ(face_support(fit.mMesh, scan, so).mNumMeasured, 6);
    EXPECT_THROW(face_support(scan, scan), std::invalid_argument);
}

TEST(FitCadIterate, StopsAtTheTargetOnAPlainBox) {
    FitCadIterateOptions o;
    o.mFitCad.mBackend = FitCadBackend::Native;
    const FitCadIterateResult r = fit_cad_iterate(analytic::box_skin(4), "", o);
    ASSERT_EQ(r.mRounds.size(), 1u);
    EXPECT_EQ(r.mStopReason, "target");
    EXPECT_TRUE(r.mTargetMet);
    EXPECT_EQ(r.mFit.mCounts.mNumFaces, 6);
    o.mPartClass = fit_cad_part_class_from_name("cast");
    EXPECT_NEAR(fit_cad_iterate(analytic::box_skin(4), "", o).mTargetWithin, 0.85, 1e-15);
    EXPECT_THROW(fit_cad_part_class_from_name("plastic"), std::invalid_argument);
    o.mRefineFactor = 1.0;
    EXPECT_THROW(fit_cad_iterate(analytic::box_skin(2), "", o), std::invalid_argument);
}

TEST(FitCadIterate, ARefinedSegmentationExplainsMore) {
    FitCadIterateOptions o;
    o.mFitCad.mBackend = FitCadBackend::Native;
    o.mSegmentTolerance = 0.02;
    o.mFitTolerance = 1e-3;
    o.mTargetExplained = 0.999;
    o.mRefineFactor = 0.1;
    const Mesh bent = fcl_bent_box(0.06);
    const FitCadIterateResult r = fit_cad_iterate(bent, "", o);
    ASSERT_GE(r.mRounds.size(), 2u);
    EXPECT_LT(r.mRounds[0].mExplained, 0.999);
    EXPECT_FALSE(r.mRounds[0].mResegmented.empty());
    // The best round's explained share is never below an earlier one's.
    EXPECT_GE(r.mRounds[static_cast<std::size_t>(r.mBestRound)].mExplained,
              r.mRounds[0].mExplained);
    EXPECT_GT(r.mRounds[1].mExplained, r.mRounds[0].mExplained);
    // Deterministic.
    const FitCadIterateResult again = fit_cad_iterate(bent, "", o);
    ASSERT_EQ(again.mRounds.size(), r.mRounds.size());
    for (std::size_t k = 0; k < r.mRounds.size(); ++k)
        EXPECT_EQ(again.mRounds[k].mExplained, r.mRounds[k].mExplained);
}
