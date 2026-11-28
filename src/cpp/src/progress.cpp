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

// Project includes
#include "meshioplusplus/progress.hpp"

namespace meshioplusplus {

namespace {

thread_local const ProgressContext* gProgCurrentContext = nullptr;

}  // namespace

const ProgressContext* current_progress_context() noexcept {
    return gProgCurrentContext;
}

const ProgressContext* set_progress_context(const ProgressContext* pContext) noexcept {
    const ProgressContext* prev = gProgCurrentContext;
    gProgCurrentContext = pContext;
    return prev;
}

ScopedProgressContext::ScopedProgressContext(const ProgressContext* pContext) noexcept
    : mPrevious(set_progress_context(pContext)) {}

ScopedProgressContext::~ScopedProgressContext() noexcept {
    set_progress_context(mPrevious);
}

bool progress_is_cancelled() noexcept {
    const ProgressContext* ctx = gProgCurrentContext;
    return ctx != nullptr && ctx->mCancel != nullptr &&
           ctx->mCancel->load(std::memory_order_relaxed);
}

void progress_checkpoint(double Fraction, const char* pStage) {
    const ProgressContext* ctx = gProgCurrentContext;
    if (ctx == nullptr)
        return;
    if (ctx->mCancel != nullptr && ctx->mCancel->load(std::memory_order_relaxed))
        throw CancelledError("operation cancelled");
    if (ctx->mCallback != nullptr) {
        const char* stage = pStage != nullptr ? pStage : "";
        if (!ctx->mCallback(Fraction, stage, ctx->mUserData))
            throw CancelledError("operation cancelled");
    }
}

}  // namespace meshioplusplus
