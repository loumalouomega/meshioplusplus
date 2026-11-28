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

// System includes
#include <atomic>
#include <string>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/adapt_volume.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/progress.hpp"

using namespace meshioplusplus;

namespace {

struct ProgRecord {
    std::vector<std::pair<double, std::string>> mEvents;
    int mCancelAfter = -1;
};

bool prog_test_cb(double Fraction, const char* pStage, void* pUserData) {
    auto* pRec = static_cast<ProgRecord*>(pUserData);
    pRec->mEvents.emplace_back(Fraction, pStage ? pStage : "");
    if (pRec->mCancelAfter >= 0 && static_cast<int>(pRec->mEvents.size()) >= pRec->mCancelAfter)
        return false;
    return true;
}

}  // namespace

TEST(Progress, ScopedGuardRestoresPreviousContext) {
    EXPECT_EQ(current_progress_context(), nullptr);
    EXPECT_FALSE(progress_is_cancelled());
    EXPECT_NO_THROW(progress_checkpoint(0.5, "noop"));

    std::atomic<bool> flag1{false};
    ProgressContext ctx1;
    ctx1.mCancel = &flag1;

    {
        ScopedProgressContext g1(&ctx1);
        EXPECT_EQ(current_progress_context(), &ctx1);
        EXPECT_FALSE(progress_is_cancelled());
        EXPECT_NO_THROW(progress_checkpoint(0.25, "outer"));

        std::atomic<bool> flag2{true};
        ProgressContext ctx2;
        ctx2.mCancel = &flag2;
        {
            ScopedProgressContext g2(&ctx2);
            EXPECT_EQ(current_progress_context(), &ctx2);
            EXPECT_TRUE(progress_is_cancelled());
            EXPECT_THROW(progress_checkpoint(0.5, "inner"), CancelledError);
        }
        EXPECT_EQ(current_progress_context(), &ctx1);
        EXPECT_FALSE(progress_is_cancelled());
    }
    EXPECT_EQ(current_progress_context(), nullptr);
}

TEST(Progress, OctreeMeshReportsStagesAndCancels) {
    const Mesh vol = box({2, 2, 2}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");
    const Mesh skin = extract_surface(vol);
    OctreeMeshOptions opts;
    opts.mMinSize = 0.5;
    opts.mMaxSize = 0.5;
    opts.mOutput = OctreeMeshOutput::Tetra;

    ProgRecord rec;
    ProgressContext ctx;
    ctx.mCallback = &prog_test_cb;
    ctx.mUserData = &rec;
    {
        ScopedProgressContext guard(&ctx);
        const OctreeMeshResult res = octree_mesh(skin, opts);
        EXPECT_GT(res.mMesh.NumCellBlocks(), 0u);
    }
    ASSERT_GE(rec.mEvents.size(), 3u);
    EXPECT_EQ(rec.mEvents.front().second, "octree_mesh:prepare");
    EXPECT_EQ(rec.mEvents.back().second, "octree_mesh:done");

    // Pre-cancelled token stops immediately at the first checkpoint.
    std::atomic<bool> cancel{true};
    ProgressContext cancel_ctx;
    cancel_ctx.mCancel = &cancel;
    {
        ScopedProgressContext guard(&cancel_ctx);
        EXPECT_THROW(octree_mesh(skin, opts), CancelledError);
    }

    // Callback returning false after 2 events cancels mid-pipeline.
    ProgRecord cancel_rec;
    cancel_rec.mCancelAfter = 2;
    ProgressContext cb_cancel_ctx;
    cb_cancel_ctx.mCallback = &prog_test_cb;
    cb_cancel_ctx.mUserData = &cancel_rec;
    {
        ScopedProgressContext guard(&cb_cancel_ctx);
        EXPECT_THROW(octree_mesh(skin, opts), CancelledError);
    }
    EXPECT_EQ(cancel_rec.mEvents.size(), 2u);
}

TEST(Progress, OptimizeVolumeAndAdaptVolumeCancel) {
    const Mesh vol = box({2, 2, 2}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");

    std::atomic<bool> cancel{true};
    ProgressContext cancel_ctx;
    cancel_ctx.mCancel = &cancel;

    OptimizeVolumeOptions ov_opts;
    ov_opts.mMaxIterations = 3;
    {
        ScopedProgressContext guard(&cancel_ctx);
        EXPECT_THROW(optimize_volume(vol, ov_opts), CancelledError);
    }

    AdaptVolumeOptions av_opts;
    av_opts.mHmax = 0.4;
    av_opts.mMaxIterations = 3;
    {
        ScopedProgressContext guard(&cancel_ctx);
        EXPECT_THROW(adapt_volume(vol, av_opts), CancelledError);
    }
}
