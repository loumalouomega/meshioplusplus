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
#pragma once

/**
 * @file progress.hpp
 * @brief Thread-local cooperative cancellation and progress reporting.
 *
 * Long-running operations call `progress_checkpoint()` at outer-pass and
 * stage boundaries (never inside per-cell inner loops). When no
 * `ProgressContext` is active on the calling thread, the check reduces to a
 * single null-pointer test with negligible overhead. Passing the context via
 * a scoped thread-local guard (`ScopedProgressContext`) avoids adding fields
 * to every operation's options struct.
 */

// System includes
#include <atomic>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {

/**
 * @brief C-compatible progress callback signature.
 *
 * @param Fraction Progress in `[0.0, 1.0]` (or `< 0.0` when indeterminate).
 * @param pStage Non-null NUL-terminated stage label (may be empty `""`).
 * @param pUserData Caller-owned opaque pointer forwarded unchanged.
 * @return `true` to continue, `false` to request cancellation.
 */
using ProgressCallbackFn = bool (*)(double Fraction, const char* pStage, void* pUserData);

/**
 * @brief Lightweight, non-owning cancellation and progress context.
 *
 * Both members are optional; a context with both null is a no-op.
 */
struct ProgressContext {
    std::atomic<bool>* mCancel = nullptr;
    ProgressCallbackFn mCallback = nullptr;
    void* mUserData = nullptr;
};

/**
 * @brief Returns the active thread-local `ProgressContext`, or `nullptr`.
 */
MESHIOPLUSPLUS_API const ProgressContext* current_progress_context() noexcept;

/**
 * @brief Replaces the active thread-local `ProgressContext` and returns the
 * previous one. Prefer `ScopedProgressContext` in C++ code.
 */
MESHIOPLUSPLUS_API const ProgressContext* set_progress_context(
    const ProgressContext* pContext) noexcept;

/**
 * @brief RAII guard installing a `ProgressContext` for the current thread.
 *
 * Non-copyable and non-movable; restores the previous context on destruction.
 * The referenced `ProgressContext` must outlive the guard.
 */
class MESHIOPLUSPLUS_API ScopedProgressContext {
public:
    explicit ScopedProgressContext(const ProgressContext* pContext) noexcept;
    ~ScopedProgressContext() noexcept;

    ScopedProgressContext(const ScopedProgressContext&) = delete;
    ScopedProgressContext& operator=(const ScopedProgressContext&) = delete;
    ScopedProgressContext(ScopedProgressContext&&) = delete;
    ScopedProgressContext& operator=(ScopedProgressContext&&) = delete;

private:
    const ProgressContext* mPrevious = nullptr;
};

/**
 * @brief Returns `true` if the active context's cancel flag is currently set.
 *
 * Does not invoke the progress callback and never throws.
 */
MESHIOPLUSPLUS_API bool progress_is_cancelled() noexcept;

/**
 * @brief Reports progress and checks for cancellation at an outer-loop boundary.
 *
 * Throws `CancelledError` if the active context's `mCancel` flag is set or if
 * `mCallback` returns `false`. No-op when no context is installed.
 */
MESHIOPLUSPLUS_API void progress_checkpoint(double Fraction, const char* pStage = "");

}  // namespace meshioplusplus
