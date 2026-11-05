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
 * @file parallel_threads.cpp
 * @brief The `THREADS` parallel backend: a process-wide pool of `std::thread`
 * workers behind `detail::threads_run` (roadmap 5.2).
 *
 * The backend needs no runtime library beyond the C++ standard one, which is
 * why the Linux wheels can link it statically where libgomp cannot be.
 *
 * - **One pool per process**, created on first use and never destroyed: its
 *   workers are never joined from a static destructor (which deadlocks under
 *   the Windows loader lock) and simply end with the process.
 * - **The caller participates**: a job of `P` participants wakes `P - 1`
 *   workers; every participant claims chunks of the range through one atomic
 *   counter. The iterations are independent, so which thread runs which
 *   chunk changes nothing in the result.
 * - **Nested calls run inline**: a `threads_run` from inside a job (a
 *   `parallel_for` within a `parallel_for` body) runs its range on the calling
 *   participant, so the pool can never wait on itself.
 * - **Concurrent callers run inline**: a second thread that finds the pool
 *   busy (two Python threads with the GIL released) runs its range itself
 *   instead of queueing, so it never waits on a job it is not part of.
 * - **Fork**: the child of a `fork()` inherits the pool's memory but none of
 *   its threads, so an `atfork` handler drops the pool there (leaked, as its
 *   threads do not exist) and the child's first call builds a new one.
 * - **Size**: `MESHIOPLUSPLUS_NUM_THREADS`, else `OMP_NUM_THREADS` (so the
 *   OpenMP sweeps of `tools/bench_ops.sh` cap it too), else
 *   `std::thread::hardware_concurrency()`; `threads_set_concurrency` overrides
 *   it at run time. The pool only grows: lowering the count leaves workers
 *   idle rather than stopping them.
 */

#if defined(MESHIOPLUSPLUS_PARALLEL_THREADS)

// System includes
#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <mutex>
#include <thread>
#include <vector>

#if !defined(_WIN32)
#include <pthread.h>
#endif

// Project includes
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus::detail {
namespace {

/// The most threads a count from the environment may ask for.
constexpr unsigned kPthrMaxThreads = 1024;

/// Nonzero while this thread runs a chunk of a job (a worker always; a caller
/// while it participates): a `threads_run` made there runs inline.
thread_local unsigned tPthrDepth = 0;

/// A positive integer from the start of `pText` (OpenMP's `OMP_NUM_THREADS`
/// may be a comma-separated list, of which the first entry is the outer
/// level's), or 0.
unsigned pthr_parse_count(const char* pText) {
    if (!pText)
        return 0;
    unsigned long v = 0;
    const char* p = pText;
    while (*p == ' ' || *p == '\t')
        ++p;
    if (*p < '0' || *p > '9')
        return 0;
    for (; *p >= '0' && *p <= '9'; ++p) {
        v = v * 10 + static_cast<unsigned long>(*p - '0');
        if (v > kPthrMaxThreads)
            return kPthrMaxThreads;
    }
    return static_cast<unsigned>(v);
}

/// The count the environment asks for, else the hardware's; at least 1.
unsigned pthr_default_concurrency() {
    if (const unsigned n = pthr_parse_count(std::getenv("MESHIOPLUSPLUS_NUM_THREADS")))
        return n;
    if (const unsigned n = pthr_parse_count(std::getenv("OMP_NUM_THREADS")))
        return n;
    return std::max(1u, std::min(kPthrMaxThreads, std::thread::hardware_concurrency()));
}

/// 0 = the default; set by `threads_set_concurrency`.
std::atomic<unsigned> gPthrOverride{0};

/// One `threads_run` call: the range, its chunk size and the body.
struct PthrJob {
    ThreadsChunkFn mFn;
    void* mpCtx;
    std::size_t mN;
    std::size_t mChunk;
    std::atomic<std::size_t> mNext{0};

    /// Claims and runs chunks until none are left.
    void Work() {
        for (;;) {
            const std::size_t b = mNext.fetch_add(mChunk, std::memory_order_relaxed);
            if (b >= mN)
                return;
            mFn(mpCtx, b, std::min(mN, b + mChunk));
        }
    }
};

class PthrPool {
public:
    /// Runs `rJob` on the caller and up to `Helpers` workers; false, running
    /// nothing, when another thread's job holds the pool.
    bool TryRun(PthrJob& rJob, unsigned Helpers) {
        std::unique_lock<std::mutex> submit(mSubmit, std::try_to_lock);
        if (!submit.owns_lock())
            return false;
        Helpers = Grow(Helpers);
        {
            std::lock_guard<std::mutex> lock(mMutex);
            mpJob = &rJob;
            mWanted = Helpers;
            mJoined = 0;
            ++mGeneration;
        }
        if (Helpers)
            mWake.notify_all();
        ++tPthrDepth;
        rJob.Work();
        --tPthrDepth;
        std::unique_lock<std::mutex> lock(mMutex);
        // Every chunk is claimed; wait for the workers still inside one, and
        // close the job to the ones that have not joined yet.
        mDone.wait(lock, [&] { return mRunning == 0; });
        mpJob = nullptr;
        mWanted = 0;
        return true;
    }

private:
    /// Starts workers until there are `Count`; called under `mSubmit`. A
    /// thread the system refuses leaves the job with the workers there are.
    unsigned Grow(unsigned Count) {
        try {
            while (mThreads.size() < Count)
                mThreads.emplace_back([this] { Loop(); });
        } catch (...) {
            return static_cast<unsigned>(mThreads.size());
        }
        return Count;
    }

    void Loop() {
        tPthrDepth = 1;
        std::uint64_t seen = 0;
        std::unique_lock<std::mutex> lock(mMutex);
        for (;;) {
            mWake.wait(lock, [&] { return mGeneration != seen; });
            seen = mGeneration;
            if (!mpJob || mJoined >= mWanted)
                continue;
            ++mJoined;
            ++mRunning;
            PthrJob* p_job = mpJob;
            lock.unlock();
            p_job->Work();
            lock.lock();
            if (--mRunning == 0)
                mDone.notify_one();
        }
    }

    std::mutex mSubmit;  ///< held by the one caller whose job runs
    std::mutex mMutex;   ///< guards everything below
    std::condition_variable mWake;
    std::condition_variable mDone;
    PthrJob* mpJob = nullptr;
    std::uint64_t mGeneration = 0;
    unsigned mWanted = 0;   ///< workers the current job may take
    unsigned mJoined = 0;   ///< workers that took it
    unsigned mRunning = 0;  ///< workers inside it now
    std::vector<std::thread> mThreads;
};

std::mutex gPthrCreate;
std::atomic<PthrPool*> gPthrPool{nullptr};

#if !defined(_WIN32)
void pthr_atfork_prepare() {
    gPthrCreate.lock();
}
void pthr_atfork_parent() {
    gPthrCreate.unlock();
}
void pthr_atfork_child() {
    // The forking thread is the child's only thread and holds gPthrCreate;
    // the pool's workers did not survive, so it is abandoned (not freed: its
    // mutexes may be held by threads that no longer exist).
    gPthrPool.store(nullptr, std::memory_order_relaxed);
    tPthrDepth = 0;
    gPthrCreate.unlock();
}
#endif

PthrPool& pthr_pool() {
    PthrPool* p = gPthrPool.load(std::memory_order_acquire);
    if (p)
        return *p;
    static std::once_flag registered;
    std::call_once(registered, [] {
#if !defined(_WIN32)
        pthread_atfork(pthr_atfork_prepare, pthr_atfork_parent, pthr_atfork_child);
#endif
    });
    std::lock_guard<std::mutex> lock(gPthrCreate);
    p = gPthrPool.load(std::memory_order_acquire);
    if (!p) {
        p = new PthrPool();  // never deleted: see the file comment
        gPthrPool.store(p, std::memory_order_release);
    }
    return *p;
}

}  // namespace

unsigned threads_concurrency() {
    const unsigned set = gPthrOverride.load(std::memory_order_relaxed);
    if (set)
        return set;
    static const unsigned from_env = pthr_default_concurrency();
    return from_env;
}

void threads_set_concurrency(unsigned N) {
    gPthrOverride.store(std::min(N, kPthrMaxThreads), std::memory_order_relaxed);
}

void threads_run(std::size_t N, std::size_t Chunk, unsigned MaxThreads, ThreadsChunkFn Fn,
                 void* pCtx) {
    if (N == 0)
        return;
    Chunk = std::max<std::size_t>(Chunk, 1);
    const std::size_t nchunks = (N + Chunk - 1) / Chunk;
    unsigned participants = threads_concurrency();
    if (MaxThreads)
        participants = std::min(participants, MaxThreads);
    participants = static_cast<unsigned>(std::min<std::size_t>(participants, nchunks));
    if (participants > 1 && tPthrDepth == 0) {
        PthrJob job{Fn, pCtx, N, Chunk};
        if (pthr_pool().TryRun(job, participants - 1))
            return;
    }
    Fn(pCtx, 0, N);
}

}  // namespace meshioplusplus::detail

#endif  // MESHIOPLUSPLUS_PARALLEL_THREADS
