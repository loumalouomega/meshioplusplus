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

#include <atomic>
#include <cstddef>
#include <exception>
#include <utility>

#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
#include <algorithm>
#include <execution>
#include <thread>
#include <vector>
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#endif

namespace meshioplusplus {

inline constexpr std::size_t parallel_grain_default = 2048;

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
void parallel_for_impl(std::size_t n, F& f, std::size_t grain) {
#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
    struct Chunk {
        std::size_t begin, end;
    };
    const std::size_t hw = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    const std::size_t max_chunks = hw * 4;
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
    (void)grain;
    first_exception exc;
    const long long nn = static_cast<long long>(n);
#pragma omp parallel for schedule(static)
    for (long long i = 0; i < nn; ++i) {
        exc.run([&] { f(static_cast<std::size_t>(i)); });
    }
    exc.rethrow_if_any();
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    first_exception exc;
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, n, grain),
                      [&](const tbb::blocked_range<std::size_t>& r) {
                          exc.run([&] {
                              for (std::size_t i = r.begin(); i != r.end(); ++i) f(i);
                          });
                      });
    exc.rethrow_if_any();
#else  // MESHIOPLUSPLUS_PARALLEL_SEQ (and the safe default)
    (void)grain;
    for (std::size_t i = 0; i < n; ++i) f(i);
#endif
}

}  // namespace detail

template <class F>
void parallel_for(std::size_t n, F&& f,
                  std::size_t grain = parallel_grain_default) {
    if (n == 0) return;
    if (n <= grain) {
        for (std::size_t i = 0; i < n; ++i) f(i);
        return;
    }
    detail::parallel_for_impl(n, f, grain);
}

}  // namespace meshioplusplus
