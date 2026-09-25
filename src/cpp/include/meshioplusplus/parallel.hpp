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
 * @file parallel.hpp
 * @brief `parallel_for`/`parallel_for_bw`: a backend-agnostic parallel loop
 * over a compile-time-selected SEQ/STL/OpenMP/TBB implementation.
 *
 * The active backend is chosen at compile time by the `MESHIOPLUSPLUS_PARALLEL_*`
 * preprocessor definitions (set from CMake's `MESHIOPLUSPLUS_PARALLEL_BACKEND` =
 * `AUTO|SEQ|STL|OPENMP|TBB|KOKKOS`; `AUTO` prefers OpenMP — portable across
 * manylinux/MSVC/macOS without needing TBB — then falls back to STL(+TBB) if
 * detected, else SEQ; KOKKOS is bring-your-own and never picked by AUTO).
 * `parallel_backend_name()`/`_core.__parallel_backend__`
 * report which one is active. Iterations passed to `parallel_for` must be
 * independent (no cross-iteration state) since they may run concurrently in
 * any order; the first exception thrown by any iteration is captured and
 * rethrown once the parallel region has joined (via `detail::FirstException`),
 * so callers see ordinary C++ exception semantics rather than `std::terminate`
 * or a lost exception.
 *
 * There are two flavors, distinguished by how many threads they are allowed
 * to use:
 *  - `parallel_for` — uses all available cores (up to `max_threads` if
 *    non-zero). Appropriate for compute-bound loops where per-element work
 *    is real computation, e.g. zlib/base64 encode-decode in
 *    `detail/vtu_binary.hpp` and ASCII value formatting.
 *  - `parallel_for_bw` — caps the thread count to `parallel_bandwidth_threads`
 *    (4). Appropriate for memory-bandwidth-bound loops — byte-swap,
 *    transpose, index gather — which saturate a socket's memory bandwidth
 *    with only a few threads and then *regress* as thread count grows
 *    further (more cache contention and dispatch overhead without more
 *    usable bandwidth), unlike compute-bound loops which keep scaling to all
 *    cores.
 *
 * Three primitives build on them (v16.16.0): `parallel_sort`, which needs a
 * comparator that is a **total order** on the elements -- then every backend
 * and thread count yields the same sequence, so no stability is needed;
 * `parallel_reduce`, whose chunk size is a constant the caller picks, never
 * the thread count, with partials combined in chunk order; and
 * `parallel_exclusive_scan` over integers, exact under any chunking. Together
 * they are what the sort-based facet and edge tables (roadmap §4) are built
 * from, and each is deterministic by construction.
 *
 * To add a new backend (e.g. HPX): add one CMake branch that defines
 * a new `MESHIOPLUSPLUS_PARALLEL_<NAME>` macro and links the dependency, then
 * add one `#elif defined(MESHIOPLUSPLUS_PARALLEL_<NAME>)` branch in
 * `detail::parallel_for_impl` below (and extend `parallel_backend_name()`
 * to report it, plus the guarded include block).
 *
 * The KOKKOS backend runs on `Kokkos::DefaultHostExecutionSpace` DELIBERATELY,
 * even in a build whose Kokkos has a device (CUDA/HIP/SYCL) enabled: every
 * `parallel_for` body in this codebase captures host pointers
 * (`NDArray::Data()`, `std::vector`, `std::string`) by reference and several
 * call host-only libraries (zlib), so device dispatch is not meaningful here —
 * GPU data movement is served by the DLPack/CuPy handoff (`doc/gpu.md`)
 * instead. Kokkos is lazily initialized on first use ONLY if the embedding
 * application has not already called `Kokkos::initialize()` itself; in that
 * case a matching `Kokkos::finalize()` is registered with `std::atexit`. An
 * application that wants full control of the Kokkos lifecycle should simply
 * initialize Kokkos before its first meshio++ call.
 */

// System includes
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
#include <execution>
#include <thread>
#endif

// External includes
#if defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
#include <omp.h>
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#elif defined(MESHIOPLUSPLUS_PARALLEL_KOKKOS)
#include <Kokkos_Core.hpp>

#include <cstdint>
#include <cstdlib>
#include <mutex>
#endif

namespace meshioplusplus {

/**
 * @brief Default grain size (minimum iterations per dispatched chunk) for
 * `parallel_for`/`parallel_for_bw` when the caller doesn't override it.
 *
 * Below this many total iterations, `parallel_for` runs sequentially rather
 * than paying parallel dispatch overhead (see the `n <= grain` check in
 * `parallel_for` below). Callers with atypically coarse or fine per-iteration
 * work (e.g. one whole zlib block per iteration) pass an explicit smaller
 * `grain` (often `1`) so each iteration dispatches individually.
 */
inline constexpr std::size_t parallel_grain_default = 2048;

/**
 * @brief Thread cap used by `parallel_for_bw` for memory-bandwidth-bound loops.
 *
 * Memory-bandwidth-bound loops (byte-swap, transpose, gather) saturate a
 * socket's bandwidth with only a few threads and then *regress* as thread
 * overhead and cache contention grow — unlike compute-bound loops (zlib,
 * base64) which scale to all cores. Cap the bandwidth-bound loops here.
 */
inline constexpr unsigned parallel_bandwidth_threads = 4;

/**
 * @brief Name of the parallel backend selected at compile time.
 *
 * Reflects whichever of `MESHIOPLUSPLUS_PARALLEL_STL`/`_OPENMP`/`_TBB` was
 * defined (by CMake, based on `MESHIOPLUSPLUS_PARALLEL_BACKEND`); none of
 * them defined means the sequential fallback. Exposed to Python as
 * `_core.__parallel_backend__` so tests/diagnostics can assert which backend
 * actually built.
 * @return One of `"stl"`, `"openmp"`, `"tbb"`, `"kokkos"`, `"seq"`.
 */
constexpr const char* parallel_backend_name() {
#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
    return "stl";
#elif defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    return "openmp";
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    return "tbb";
#elif defined(MESHIOPLUSPLUS_PARALLEL_KOKKOS)
    return "kokkos";
#else
    return "seq";
#endif
}

namespace detail {

#if defined(MESHIOPLUSPLUS_PARALLEL_KOKKOS)
/**
 * @brief Lazily initializes Kokkos on first `parallel_for`, once per process.
 *
 * A library must not fight its host over the Kokkos lifecycle: if the
 * embedding application already called `Kokkos::initialize()`, this is a
 * no-op and meshio++ never finalizes (the host owns the lifecycle). Only when
 * nobody has initialized Kokkos yet does meshio++ initialize it (honouring
 * the `KOKKOS_*` environment variables) and register a guarded
 * `Kokkos::finalize()` with `std::atexit`. The `is_finalized()` check avoids
 * the init-after-finalize abort if a host tears Kokkos down and a straggling
 * meshio++ call arrives after that.
 */
inline void kokkos_ensure_initialized() {
    static std::once_flag once;
    std::call_once(once, [] {
    // Kokkos::is_finalized() only exists from Kokkos 3.7 on.
#if defined(KOKKOS_VERSION) && KOKKOS_VERSION >= 30700
        const bool can_init = !Kokkos::is_initialized() && !Kokkos::is_finalized();
#else
        const bool can_init = !Kokkos::is_initialized();
#endif
        if (can_init) {
            Kokkos::initialize();
            std::atexit([] {
                if (Kokkos::is_initialized())
                    Kokkos::finalize();
            });
        }
    });
}
#endif

/**
 * @brief Captures the first exception thrown by any parallel iteration, to
 * be rethrown by the caller after the parallel region joins.
 *
 * Iterations run on multiple threads cannot let a C++ exception escape
 * across the parallelism boundary (OpenMP/TBB would `std::terminate`), so
 * each backend wraps its per-iteration body in `Run()`, which catches
 * everything and records only the *first* exception (subsequent ones from
 * other threads are discarded — `mRaised` is a one-shot latch via
 * `std::atomic_flag`). After the parallel region has fully joined, the
 * caller calls `RethrowIfAny()` to surface that exception on the calling
 * thread with normal C++ semantics.
 */
class FirstException {
public:
    template <class Body>
    void Run(Body&& body) noexcept {
        try {
            body();
        } catch (...) {
            if (!mRaised.test_and_set(std::memory_order_acq_rel))
                mEptr = std::current_exception();
        }
    }
    void RethrowIfAny() {
        if (mEptr)
            std::rethrow_exception(mEptr);
    }

private:
    std::atomic_flag mRaised = ATOMIC_FLAG_INIT;
    std::exception_ptr mEptr;
};

/**
 * @brief Backend-specific dispatch of `n` independent iterations of `f`.
 *
 * Exactly one `#if`/`#elif` branch compiles, selected by the
 * `MESHIOPLUSPLUS_PARALLEL_*` macro CMake defined:
 *  - **STL**: splits `[0, n)` into up to `hardware_concurrency() * 4` chunks
 *    (fewer if `grain`/`max_threads` constrain it further) and runs them via
 *    `std::for_each(std::execution::par, ...)` over a small chunk table
 *    (iterated explicitly because PSTL algorithms require
 *    `Cpp17ForwardIterator`s, which `iota_view` iterators don't satisfy on
 *    every implementation).
 *  - **OpenMP**: `#pragma omp parallel for schedule(dynamic, chunk)` with
 *    `chunk = max(grain/4, 1)`. Dynamic (not static) scheduling matters on
 *    hybrid P+E-core CPUs, where a static split would leave slow E-cores as
 *    stragglers while fast P-cores idle at the join; `grain/4` keeps
 *    dispatch overhead negligible for fine-grained loops while still
 *    honouring explicitly coarse callers (e.g. VTU zlib blocks pass
 *    `grain=1` because each iteration is already a whole compress, so
 *    per-iteration dispatch is exactly what's wanted — the chunk size must
 *    never be floored above the caller's `grain`).
 *  - **TBB**: `tbb::parallel_for` over a `blocked_range` of grain size
 *    `grain`, optionally under a `tbb::global_control` limiting
 *    `max_allowed_parallelism` to `max_threads`.
 *  - **KOKKOS**: `Kokkos::parallel_for` over a `RangePolicy` pinned to
 *    `Kokkos::DefaultHostExecutionSpace` (host deliberately — see the file
 *    header). Kokkos has no per-call thread cap, so a non-zero `max_threads`
 *    is honoured the way the STL branch does it: the range is partitioned
 *    into at most `max_threads` coarse chunks, so at most that many threads
 *    have work — which is what preserves `parallel_for_bw`'s bandwidth-cap
 *    semantics. Plain `[&]` lambdas are fine (no `KOKKOS_LAMBDA` needed) —
 *    host-space functors run on ordinary host threads, which is also why
 *    `FirstException` works unchanged.
 *  - **(none, SEQ)**: a plain sequential loop; `grain`/`max_threads` are
 *    unused (cast to `void` to silence warnings).
 *
 * Every branch funnels per-iteration exceptions through a `FirstException`
 * so exactly one is rethrown after the region joins.
 *
 * @tparam F Callable invoked as `f(std::size_t i)` for each `i` in `[0, n)`.
 * @param n Number of iterations.
 * @param rF The per-iteration body (iterations must be independent).
 * @param grain Minimum unit of work per dispatched chunk/task.
 * @param max_threads Cap on threads used (0 = no cap, use all available).
 */
template <class F>
void parallel_for_impl(std::size_t n, F& rF, std::size_t grain, unsigned max_threads) {
#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
    struct Chunk {
        std::size_t mBegin, mEnd;
    };
    const std::size_t hw = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    std::size_t max_chunks = hw * 4;
    if (max_threads)
        max_chunks = std::min<std::size_t>(max_chunks, max_threads);
    const std::size_t by_grain = (n + grain - 1) / grain;
    const std::size_t nchunks = std::max<std::size_t>(1, std::min(max_chunks, by_grain));
    const std::size_t per = (n + nchunks - 1) / nchunks;
    // PSTL algorithms require Cpp17ForwardIterators (iota_view iterators do
    // not qualify on all implementations), so iterate a small chunk table.
    std::vector<Chunk> chunks;
    chunks.reserve(nchunks);
    for (std::size_t b = 0; b < n; b += per)
        chunks.push_back({b, std::min(b + per, n)});
    FirstException exc;
    std::for_each(std::execution::par, chunks.begin(), chunks.end(), [&](const Chunk& c) {
        exc.Run([&] {
            for (std::size_t i = c.mBegin; i < c.mEnd; ++i)
                rF(i);
        });
    });
    exc.RethrowIfAny();
#elif defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    FirstException exc;
    const long long nn = static_cast<long long>(n);
    const int nt = max_threads ? std::min<int>(static_cast<int>(max_threads), omp_get_max_threads())
                               : omp_get_max_threads();
    // Dynamic scheduling: on hybrid CPUs (P + E cores) a static split makes the
    // slow cores stragglers while the fast ones idle at the join; moderately
    // sized dynamic chunks self-balance with negligible dispatch overhead.
    // grain/4 keeps dispatch rare for fine-grained loops while honouring
    // explicitly coarse loops (e.g. the VTU zlib blocks pass grain=1: each
    // iteration is a whole compress, so per-iteration dispatch is ideal).
    const long long chunk = static_cast<long long>(std::max<std::size_t>(grain / 4, 1));
#pragma omp parallel for schedule(dynamic, chunk) num_threads(nt)
    for (long long i = 0; i < nn; ++i) {
        exc.Run([&] { rF(static_cast<std::size_t>(i)); });
    }
    exc.RethrowIfAny();
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    FirstException exc;
    auto body = [&] {
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, n, grain),
                          [&](const tbb::blocked_range<std::size_t>& r) {
                              exc.Run([&] {
                                  for (std::size_t i = r.begin(); i != r.end(); ++i)
                                      rF(i);
                              });
                          });
    };
    if (max_threads) {
        tbb::global_control gc(tbb::global_control::max_allowed_parallelism, max_threads);
        body();
    } else {
        body();
    }
    exc.RethrowIfAny();
#elif defined(MESHIOPLUSPLUS_PARALLEL_KOKKOS)
    kokkos_ensure_initialized();
    FirstException exc;
    using Policy =
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::IndexType<std::int64_t> >;
    if (max_threads) {
        // At most max_threads coarse chunks (the STL branch's idiom), so at
        // most that many threads have work — Kokkos has no per-call cap.
        const std::size_t by_grain = (n + grain - 1) / grain;
        const std::size_t nchunks =
            std::max<std::size_t>(1, std::min<std::size_t>(max_threads, by_grain));
        const std::size_t per = (n + nchunks - 1) / nchunks;
        Kokkos::parallel_for("meshioplusplus::parallel_for_bw",
                             Policy(0, static_cast<std::int64_t>(nchunks)), [&](std::int64_t c) {
                                 exc.Run([&] {
                                     const std::size_t b = static_cast<std::size_t>(c) * per;
                                     const std::size_t e = std::min(n, b + per);
                                     for (std::size_t i = b; i < e; ++i)
                                         rF(i);
                                 });
                             });
    } else {
        Policy pol(0, static_cast<std::int64_t>(n));
        pol.set_chunk_size(static_cast<int>(std::max<std::size_t>(grain / 4, 1)));
        Kokkos::parallel_for("meshioplusplus::parallel_for", pol, [&](std::int64_t i) {
            exc.Run([&] { rF(static_cast<std::size_t>(i)); });
        });
    }
    Kokkos::fence();
    exc.RethrowIfAny();
#else  // MESHIOPLUSPLUS_PARALLEL_SEQ (and the safe default)
    (void)grain;
    (void)max_threads;
    for (std::size_t i = 0; i < n; ++i)
        rF(i);
#endif
}

}  // namespace detail

/**
 * @brief Runs `n` independent iterations of `f(i)`, in parallel when it's
 * worthwhile, using the compile-time-selected backend (see
 * `parallel_backend_name()`).
 *
 * If `n <= grain`, runs sequentially in-line — the fixed cost of dispatching
 * a parallel region isn't worth it for small workloads. Otherwise delegates
 * to `detail::parallel_for_impl`. `f` must be safe to invoke concurrently
 * from multiple threads for different `i` (no shared mutable state without
 * external synchronization); the first exception any invocation throws is
 * captured and rethrown on the calling thread after all iterations
 * complete (partial results/side effects from other iterations are not
 * rolled back).
 *
 * @tparam F Callable invoked as `f(std::size_t i)`.
 * @param n Number of iterations; a no-op if `n == 0`.
 * @param f The per-iteration body.
 * @param grain Minimum number of iterations to bother parallelizing, and
 *              (backend-dependent) the target chunk size once it does;
 *              defaults to `parallel_grain_default` (2048). Pass a small
 *              value (e.g. `1`) when each iteration is already coarse work
 *              (a whole zlib block, a whole compress) so dispatch happens
 *              per-iteration rather than being batched further.
 * @param max_threads Cap on threads used; `0` (the default) means "use all
 *                     available". Pass `parallel_bandwidth_threads`
 *                     (or call `parallel_for_bw` instead) for
 *                     memory-bandwidth-bound loops.
 */
template <class F>
void parallel_for(std::size_t n, F&& f, std::size_t grain = parallel_grain_default,
                  unsigned max_threads = 0) {
    if (n == 0)
        return;
    if (n <= grain) {
        for (std::size_t i = 0; i < n; ++i)
            f(i);
        return;
    }
    detail::parallel_for_impl(n, f, grain, max_threads);
}

/**
 * @brief `parallel_for`, thread-capped for memory-bandwidth-bound loops.
 *
 * Convenience wrapper that forwards to `parallel_for` with
 * `max_threads = parallel_bandwidth_threads` (4). Use this for byte-swap,
 * transpose, and index-gather loops: they saturate a socket's memory
 * bandwidth with only a few threads and then *regress* — more threads add
 * cache contention and dispatch overhead without more usable bandwidth —
 * unlike genuinely compute-bound loops (zlib/base64), which should use
 * plain `parallel_for` to scale across all cores.
 *
 * @tparam F Callable invoked as `f(std::size_t i)`.
 * @param n Number of iterations; a no-op if `n == 0`.
 * @param f The per-iteration body.
 * @param grain Minimum iterations per chunk; see `parallel_for`'s `grain`.
 */
template <class F>
void parallel_for_bw(std::size_t n, F&& f, std::size_t grain = parallel_grain_default) {
    parallel_for(n, std::forward<F>(f), grain, parallel_bandwidth_threads);
}

namespace detail {

/// Below this many elements `parallel_sort` is a plain `std::sort`.
inline constexpr std::size_t parallel_sort_serial_below = std::size_t{1} << 14;

/**
 * @brief Chunk-sort-and-merge over `parallel_for`, for the backends with no
 * parallel sort of their own (OpenMP, Kokkos).
 *
 * The range is cut into a power-of-two number of runs, each sorted by one
 * task, then merged pairwise, one pass per level, every merge of a level in
 * parallel. Needs `comp` to be a total order: then the result is the unique
 * sorted sequence, whatever the run count (which follows the thread count).
 */
template <class It, class Comp>
void parallel_merge_sort(It first, It last, Comp comp, std::size_t max_runs) {
    using T = typename std::iterator_traits<It>::value_type;
    const std::size_t n = static_cast<std::size_t>(last - first);
    std::size_t runs = 1;
    while (runs * 2 <= max_runs && n / (runs * 2) >= parallel_sort_serial_below / 4)
        runs *= 2;
    const std::size_t width = (n + runs - 1) / runs;
    parallel_for(
        runs,
        [&](std::size_t r) {
            const std::size_t b = std::min(n, r * width);
            const std::size_t e = std::min(n, b + width);
            std::sort(first + static_cast<std::ptrdiff_t>(b), first + static_cast<std::ptrdiff_t>(e),
                      comp);
        },
        1);
    if (runs == 1)
        return;
    std::vector<T> buffer(n);
    bool in_buffer = false;  // where the current runs live
    for (std::size_t w = width; w < n; w *= 2) {
        const std::size_t pairs = (n + 2 * w - 1) / (2 * w);
        auto merge_level = [&](auto src, auto dst) {
            parallel_for(
                pairs,
                [&](std::size_t p) {
                    const std::size_t b = p * 2 * w;
                    const std::size_t m = std::min(n, b + w);
                    const std::size_t e = std::min(n, b + 2 * w);
                    std::merge(std::make_move_iterator(src + static_cast<std::ptrdiff_t>(b)),
                               std::make_move_iterator(src + static_cast<std::ptrdiff_t>(m)),
                               std::make_move_iterator(src + static_cast<std::ptrdiff_t>(m)),
                               std::make_move_iterator(src + static_cast<std::ptrdiff_t>(e)),
                               dst + static_cast<std::ptrdiff_t>(b), comp);
                },
                1);
        };
        if (in_buffer)
            merge_level(buffer.begin(), first);
        else
            merge_level(first, buffer.begin());
        in_buffer = !in_buffer;
    }
    if (in_buffer)
        parallel_for_bw(n, [&](std::size_t i) {
            first[static_cast<std::ptrdiff_t>(i)] = std::move(buffer[i]);
        });
}

}  // namespace detail

/**
 * @brief Sorts `[first, last)` by `comp`, in parallel where the backend has
 * the means (TBB's `parallel_sort`, the STL's `std::execution::par`, a
 * chunk-sort-and-merge over `parallel_for` for OpenMP and Kokkos), serially
 * otherwise and below `detail::parallel_sort_serial_below` elements.
 *
 * **`comp` must be a total order on the elements**: no two elements that
 * differ may compare equivalent. Then the sorted sequence is unique, and every
 * backend and thread count returns it -- which is how the callers stay
 * byte-identical across builds. Sort (key, unique slot) pairs, never bare
 * keys with ties.
 *
 * @tparam It Random-access iterator.
 * @tparam Comp Strict weak ordering, `bool(const T&, const T&)`.
 */
template <class It, class Comp>
void parallel_sort(It first, It last, Comp comp) {
    const std::size_t n = static_cast<std::size_t>(last - first);
    if (n < detail::parallel_sort_serial_below) {
        std::sort(first, last, comp);
        return;
    }
#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
    std::sort(std::execution::par, first, last, comp);
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    tbb::parallel_sort(first, last, comp);
#elif defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    detail::parallel_merge_sort(first, last, comp,
                                static_cast<std::size_t>(std::max(1, omp_get_max_threads())) * 2);
#elif defined(MESHIOPLUSPLUS_PARALLEL_KOKKOS)
    kokkos_ensure_initialized();
    detail::parallel_merge_sort(
        first, last, comp,
        static_cast<std::size_t>(std::max(1, Kokkos::DefaultHostExecutionSpace().concurrency())) *
            2);
#else
    std::sort(first, last, comp);
#endif
}

/// `parallel_sort` by `operator<`.
template <class It>
void parallel_sort(It first, It last) {
    parallel_sort(first, last, std::less<typename std::iterator_traits<It>::value_type>());
}

/**
 * @brief A reduction in fixed chunks: `chunk_fn(begin, end)` folds each chunk
 * of `Chunk` consecutive indices of `[0, n)` into a partial (in parallel), and
 * the partials are then combined **in chunk order**, left to right from
 * `init`, by `combine(acc, partial)`.
 *
 * The chunk size is the caller's constant, never derived from the thread
 * count, so the result -- floating-point sums included -- is the same on every
 * backend and at every thread count (`sobolev_deform`'s `kSoboChunk` rule).
 * It is not the plain left fold of a serial loop unless `combine` is exact
 * (integer sums, min/max).
 *
 * @tparam T Partial type; default-constructible and movable.
 */
template <class T, class ChunkFn, class Combine>
T parallel_reduce(std::size_t n, std::size_t Chunk, T init, ChunkFn&& chunk_fn,
                  Combine&& combine) {
    if (n == 0)
        return init;
    Chunk = std::max<std::size_t>(Chunk, 1);
    const std::size_t nchunks = (n + Chunk - 1) / Chunk;
    std::vector<T> partial(nchunks);
    parallel_for(
        nchunks,
        [&](std::size_t c) { partial[c] = chunk_fn(c * Chunk, std::min(n, (c + 1) * Chunk)); }, 1);
    T acc = std::move(init);
    for (T& p : partial)
        acc = combine(std::move(acc), std::move(p));
    return acc;
}

/**
 * @brief Exclusive prefix sum of `n` integers: `pOut[i] = init + pIn[0] + ...
 * + pIn[i-1]`; returns the total `init + pIn[0] + ... + pIn[n-1]`. `pOut` may
 * be `pIn` (in place).
 *
 * Two passes over fixed chunks -- per-chunk sums, then per-chunk scans from
 * each chunk's offset -- with a serial prefix over the chunk sums between.
 * Integer addition is exact, so the result is the serial one; for floating
 * point it would depend on the chunk size, which is why this is integer-only.
 */
template <class T, class U>
T parallel_exclusive_scan(const U* pIn, std::size_t n, T* pOut, T init) {
    static_assert(std::is_integral_v<T> && std::is_integral_v<U>,
                  "parallel_exclusive_scan is exact only for integers");
    constexpr std::size_t kChunk = std::size_t{1} << 15;
    const std::size_t nchunks = (n + kChunk - 1) / kChunk;
    std::vector<T> offset(nchunks + 1, T{0});
    parallel_for(
        nchunks,
        [&](std::size_t c) {
            T sum{0};
            for (std::size_t i = c * kChunk, e = std::min(n, (c + 1) * kChunk); i < e; ++i)
                sum += static_cast<T>(pIn[i]);
            offset[c + 1] = sum;
        },
        1);
    offset[0] = init;
    for (std::size_t c = 0; c < nchunks; ++c)
        offset[c + 1] += offset[c];
    parallel_for(
        nchunks,
        [&](std::size_t c) {
            T acc = offset[c];
            for (std::size_t i = c * kChunk, e = std::min(n, (c + 1) * kChunk); i < e; ++i) {
                const T v = static_cast<T>(pIn[i]);
                pOut[i] = acc;
                acc += v;
            }
        },
        1);
    return nchunks ? offset[nchunks] : init;
}

}  // namespace meshioplusplus
