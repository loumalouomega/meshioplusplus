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

#pragma once

/**
 * @file detail/thread_annotations.hpp
 * @brief Clang `-Wthread-safety` capability macros and thin mutex/lock
 * wrappers (roadmap 3.10.3).
 *
 * libstdc++'s `std::mutex` carries no Clang capability annotations, so raw
 * `std::mutex` cannot be named in `MIO_GUARDED_BY(...)`. `AnnotatedMutex`,
 * `AnnotatedRecursiveMutex`, `AnnotatedLockGuard`, and `AnnotatedUniqueLock`
 * wrap the standard primitives with zero runtime overhead so Clang checks
 * lock ownership on both libstdc++ and libc++ builds; on GCC and MSVC the
 * macros expand to nothing.
 *
 * A **core-private** header (not installed): no API or ABI surface.
 */

// System includes
#include <condition_variable>
#include <mutex>

#if defined(__clang__) && (!defined(SWIG))
#define MIO_THREAD_ANNOTATION(x) __attribute__((x))
#else
#define MIO_THREAD_ANNOTATION(x)
#endif

#define MIO_CAPABILITY(x) MIO_THREAD_ANNOTATION(capability(x))
#define MIO_SCOPED_CAPABILITY MIO_THREAD_ANNOTATION(scoped_lockable)
#define MIO_GUARDED_BY(x) MIO_THREAD_ANNOTATION(guarded_by(x))
#define MIO_PT_GUARDED_BY(x) MIO_THREAD_ANNOTATION(pt_guarded_by(x))
#define MIO_REQUIRES(...) MIO_THREAD_ANNOTATION(requires_capability(__VA_ARGS__))
#define MIO_ACQUIRE(...) MIO_THREAD_ANNOTATION(acquire_capability(__VA_ARGS__))
#define MIO_RELEASE(...) MIO_THREAD_ANNOTATION(release_capability(__VA_ARGS__))
#define MIO_TRY_ACQUIRE(...) MIO_THREAD_ANNOTATION(try_acquire_capability(__VA_ARGS__))
#define MIO_EXCLUDES(...) MIO_THREAD_ANNOTATION(locks_excluded(__VA_ARGS__))
#define MIO_RETURN_CAPABILITY(x) MIO_THREAD_ANNOTATION(lock_returned(x))
#define MIO_NO_THREAD_SAFETY_ANALYSIS MIO_THREAD_ANNOTATION(no_thread_safety_analysis)

namespace meshioplusplus {
namespace detail {

/// Zero-overhead `std::mutex` wrapper annotated as a Clang capability.
class MIO_CAPABILITY("mutex") AnnotatedMutex {
public:
    AnnotatedMutex() = default;
    AnnotatedMutex(const AnnotatedMutex&) = delete;
    AnnotatedMutex& operator=(const AnnotatedMutex&) = delete;

    void lock() MIO_ACQUIRE() { mMutex.lock(); }
    bool try_lock() MIO_TRY_ACQUIRE(true) { return mMutex.try_lock(); }
    void unlock() MIO_RELEASE() { mMutex.unlock(); }

    std::mutex& Native() { return mMutex; }

private:
    std::mutex mMutex;
};

/// Zero-overhead `std::recursive_mutex` wrapper annotated as a Clang capability.
class MIO_CAPABILITY("mutex") AnnotatedRecursiveMutex {
public:
    AnnotatedRecursiveMutex() = default;
    AnnotatedRecursiveMutex(const AnnotatedRecursiveMutex&) = delete;
    AnnotatedRecursiveMutex& operator=(const AnnotatedRecursiveMutex&) = delete;

    void lock() MIO_ACQUIRE() { mMutex.lock(); }
    bool try_lock() MIO_TRY_ACQUIRE(true) { return mMutex.try_lock(); }
    void unlock() MIO_RELEASE() { mMutex.unlock(); }

private:
    std::recursive_mutex mMutex;  // NOSONAR(cpp:S8462): explicit recursive lock
};

/// RAII scoped lock for `AnnotatedMutex` / `AnnotatedRecursiveMutex`.
template <class TMutex>
class MIO_SCOPED_CAPABILITY AnnotatedLockGuard {
public:
    explicit AnnotatedLockGuard(TMutex& rMutex) MIO_ACQUIRE(rMutex) : mrMutex(rMutex) {
        mrMutex.lock();
    }
    ~AnnotatedLockGuard() MIO_RELEASE() { mrMutex.unlock(); }
    AnnotatedLockGuard(const AnnotatedLockGuard&) = delete;
    AnnotatedLockGuard& operator=(const AnnotatedLockGuard&) = delete;

private:
    TMutex& mrMutex;
};

/// Movable-free unique lock over `AnnotatedMutex` supporting `condition_variable::wait`.
class MIO_SCOPED_CAPABILITY AnnotatedUniqueLock {
public:
    explicit AnnotatedUniqueLock(AnnotatedMutex& rMutex) MIO_ACQUIRE(rMutex)
        : mLock(rMutex.Native()) {}
    ~AnnotatedUniqueLock() MIO_RELEASE() = default;

    AnnotatedUniqueLock(const AnnotatedUniqueLock&) = delete;
    AnnotatedUniqueLock& operator=(const AnnotatedUniqueLock&) = delete;

    void lock() MIO_ACQUIRE() { mLock.lock(); }
    void unlock() MIO_RELEASE() { mLock.unlock(); }

    template <class TPred>
    void Wait(std::condition_variable& rCv, TPred&& rPred) MIO_NO_THREAD_SAFETY_ANALYSIS {
        rCv.wait(mLock, std::forward<TPred>(rPred));
    }

private:
    std::unique_lock<std::mutex> mLock;
};

}  // namespace detail
}  // namespace meshioplusplus
