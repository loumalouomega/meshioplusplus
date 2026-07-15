#pragma once
//
// meshioplusplus::parallel_for — a minimal parallel loop with a compile-time-selected
// backend, in the spirit of the C++17 parallel STL:
//
//     meshioplusplus::parallel_for(n, [&](std::size_t i) { out[i] = f(in[i]); });
//
// The backend is chosen at *configure* time via the CMake option
// MESHIOPLUSPLUS_PARALLEL_BACKEND (SEQ | STL | OPENMP | TBB), which defines exactly
// one of:
//
//     MESHIOPLUSPLUS_PARALLEL_SEQ      plain sequential loop (always correct)
//     MESHIOPLUSPLUS_PARALLEL_STL      C++17 std::execution::par (default; on
//                              libstdc++ this requires linking TBB — CMake
//                              probes for it and falls back to SEQ)
//     MESHIOPLUSPLUS_PARALLEL_OPENMP   #pragma omp parallel for
//     MESHIOPLUSPLUS_PARALLEL_TBB      tbb::parallel_for
//
// Adding a new backend (Kokkos, HPX, ...) is a self-contained change:
//   1. add a MESHIOPLUSPLUS_PARALLEL_<NAME> branch to the CMake backend block,
//   2. add a matching #elif block in parallel_for_impl() below, and
//   3. extend parallel_backend_name().
// Nothing else in the codebase refers to the backend.
//
// Contract:
//   * f(i) is invoked exactly once for every i in [0, n), in unspecified
//     order, possibly concurrently. Iterations must be independent.
//   * When n <= grain (default 2048) the loop runs sequentially on the
//     calling thread — small workloads never pay threading overhead.
//   * If f throws, the first exception is captured and rethrown on the
//     calling thread after all iterations complete (the parallel STL would
//     otherwise std::terminate, and OpenMP regions must not leak exceptions).

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <utility>

#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
#include <execution>
#include <thread>
#include <vector>
#elif defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
#include <omp.h>
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#endif

namespace meshioplusplus {

inline constexpr std::size_t parallel_grain_default = 2048;

// Memory-bandwidth-bound loops (byte-swap, transpose, gather) saturate a
// socket's bandwidth with only a few threads and then *regress* as thread
// overhead and cache contention grow — unlike compute-bound loops (zlib,
// base64) which scale to all cores. Cap the bandwidth-bound loops here.
inline constexpr unsigned parallel_bandwidth_threads = 4;

constexpr const char* parallel_backend_name() {
#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
    return "stl";
#elif defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    return "openmp";
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    return "tbb";
#else
    return "seq";
#endif
}

namespace detail {

// Captures the first exception thrown by any iteration; rethrown by the
// caller after the parallel region joins.
class first_exception {
public:
    template <class Body>
    void run(Body&& body) noexcept {
        try {
            body();
        } catch (...) {
            if (!raised_.test_and_set(std::memory_order_acq_rel))
                eptr_ = std::current_exception();
        }
    }
    void rethrow_if_any() {
        if (eptr_) std::rethrow_exception(eptr_);
    }

private:
    std::atomic_flag raised_ = ATOMIC_FLAG_INIT;
    std::exception_ptr eptr_;
};

template <class F>
void parallel_for_impl(std::size_t n, F& f, std::size_t grain, unsigned max_threads) {
#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
    struct Chunk {
        std::size_t begin, end;
    };
    const std::size_t hw = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    std::size_t max_chunks = hw * 4;
    if (max_threads) max_chunks = std::min<std::size_t>(max_chunks, max_threads);
    const std::size_t by_grain = (n + grain - 1) / grain;
    const std::size_t nchunks = std::max<std::size_t>(1, std::min(max_chunks, by_grain));
    const std::size_t per = (n + nchunks - 1) / nchunks;
    // PSTL algorithms require Cpp17ForwardIterators (iota_view iterators do
    // not qualify on all implementations), so iterate a small chunk table.
    std::vector<Chunk> chunks;
    chunks.reserve(nchunks);
    for (std::size_t b = 0; b < n; b += per) chunks.push_back({b, std::min(b + per, n)});
    first_exception exc;
    std::for_each(std::execution::par, chunks.begin(), chunks.end(),
                  [&](const Chunk& c) {
                      exc.run([&] {
                          for (std::size_t i = c.begin; i < c.end; ++i) f(i);
                      });
                  });
    exc.rethrow_if_any();
#elif defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    first_exception exc;
    const long long nn = static_cast<long long>(n);
    const int nt = max_threads ? std::min<int>(static_cast<int>(max_threads),
                                               omp_get_max_threads())
                               : omp_get_max_threads();
    // Dynamic scheduling: on hybrid CPUs (P + E cores) a static split makes the
    // slow cores stragglers while the fast ones idle at the join; moderately
    // sized dynamic chunks self-balance with negligible dispatch overhead.
    const long long chunk =
        static_cast<long long>(std::max<std::size_t>(grain / 4, 256));
#pragma omp parallel for schedule(dynamic, chunk) num_threads(nt)
    for (long long i = 0; i < nn; ++i) {
        exc.run([&] { f(static_cast<std::size_t>(i)); });
    }
    exc.rethrow_if_any();
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    first_exception exc;
    auto body = [&] {
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, n, grain),
                          [&](const tbb::blocked_range<std::size_t>& r) {
                              exc.run([&] {
                                  for (std::size_t i = r.begin(); i != r.end(); ++i) f(i);
                              });
                          });
    };
    if (max_threads) {
        tbb::global_control gc(tbb::global_control::max_allowed_parallelism, max_threads);
        body();
    } else {
        body();
    }
    exc.rethrow_if_any();
#else  // MESHIOPLUSPLUS_PARALLEL_SEQ (and the safe default)
    (void)grain;
    (void)max_threads;
    for (std::size_t i = 0; i < n; ++i) f(i);
#endif
}

}  // namespace detail

// max_threads == 0 means "use all available"; pass parallel_bandwidth_threads
// (or parallel_for_bw below) for memory-bandwidth-bound loops.
template <class F>
void parallel_for(std::size_t n, F&& f, std::size_t grain = parallel_grain_default,
                  unsigned max_threads = 0) {
    if (n == 0) return;
    if (n <= grain) {
        for (std::size_t i = 0; i < n; ++i) f(i);
        return;
    }
    detail::parallel_for_impl(n, f, grain, max_threads);
}

// Convenience for memory-bandwidth-bound loops: caps the thread count to avoid
// the over-subscription regression measured for byte-swap/transpose/gather.
template <class F>
void parallel_for_bw(std::size_t n, F&& f, std::size_t grain = parallel_grain_default) {
    parallel_for(n, std::forward<F>(f), grain, parallel_bandwidth_threads);
}

}  // namespace meshioplusplus
