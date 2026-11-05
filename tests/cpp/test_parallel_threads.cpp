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

// The THREADS backend's pool (roadmap 5.2): nesting, exceptions, concurrent
// callers, the bandwidth cap, results independent of the thread count, and a
// fork. Every test is skipped under the other backends.

// System includes
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <mutex>
#include <numeric>
#include <set>
#include <stdexcept>
#include <thread>
#include <vector>

#if !defined(_WIN32)
#include <sys/wait.h>
#include <unistd.h>
#endif

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/parallel.hpp"
#include "thread_sweep.hpp"

namespace {

namespace mio = meshioplusplus;

#if defined(MESHIOPLUSPLUS_PARALLEL_THREADS)
constexpr bool kPtThreads = true;
#else
constexpr bool kPtThreads = false;
#endif

#define PT_REQUIRE_THREADS()                                     \
    do {                                                         \
        if (!kPtThreads)                                         \
            GTEST_SKIP() << "needs the THREADS parallel backend"; \
    } while (0)

/// Restores the thread count a test changed.
struct PtCount {
    explicit PtCount(int Count) : mSaved(mio_test::max_threads()) { mio_test::set_threads(Count); }
    ~PtCount() { mio_test::set_threads(mSaved); }
    int mSaved;
};

TEST(ParallelThreads, NestedLoopsComplete) {
    PT_REQUIRE_THREADS();
    PtCount count(4);
    std::vector<std::int64_t> out(64 * 5000, 0);
    mio::parallel_for(
        64,
        [&](std::size_t i) {
            mio::parallel_for(5000, [&](std::size_t j) {
                out[i * 5000 + j] = static_cast<std::int64_t>(i * 5000 + j);
            }, 16);
        },
        1);
    for (std::size_t k = 0; k < out.size(); ++k)
        ASSERT_EQ(out[k], static_cast<std::int64_t>(k));
}

TEST(ParallelThreads, FirstExceptionReachesTheCaller) {
    PT_REQUIRE_THREADS();
    PtCount count(4);
    std::atomic<int> ran{0};
    EXPECT_THROW(mio::parallel_for(
                     100000,
                     [&](std::size_t i) {
                         ++ran;
                         if (i == 77777)
                             throw std::runtime_error("boom");
                     },
                     64),
                 std::runtime_error);
    EXPECT_GT(ran.load(), 0);
    // The pool still works after a job that threw.
    std::atomic<std::size_t> sum{0};
    mio::parallel_for(100000, [&](std::size_t i) { sum += i; }, 64);
    EXPECT_EQ(sum.load(), std::size_t{100000} * 99999 / 2);
}

TEST(ParallelThreads, ConcurrentCallers) {
    PT_REQUIRE_THREADS();
    PtCount count(4);
    constexpr std::size_t kN = 200000;
    std::vector<std::uint64_t> a(kN), b(kN);
    auto fill = [](std::vector<std::uint64_t>& rV, std::uint64_t Mul) {
        for (int rep = 0; rep < 20; ++rep)
            mio::parallel_for(rV.size(), [&](std::size_t i) { rV[i] = Mul * i + rep; }, 256);
    };
    std::thread t1([&] { fill(a, 3); });
    std::thread t2([&] { fill(b, 5); });
    t1.join();
    t2.join();
    for (std::size_t i = 0; i < kN; ++i) {
        ASSERT_EQ(a[i], 3 * i + 19);
        ASSERT_EQ(b[i], 5 * i + 19);
    }
}

TEST(ParallelThreads, BandwidthLoopsUseAtMostFourThreads) {
    PT_REQUIRE_THREADS();
    PtCount count(8);
    std::mutex m;
    std::set<std::thread::id> ids;
    mio::parallel_for_bw(
        1 << 16,
        [&](std::size_t) {
            std::lock_guard<std::mutex> lock(m);
            ids.insert(std::this_thread::get_id());
        },
        16);
    EXPECT_LE(ids.size(), mio::parallel_bandwidth_threads);
    EXPECT_GE(ids.size(), 1u);
}

TEST(ParallelThreads, ResultsDoNotDependOnTheThreadCount) {
    PT_REQUIRE_THREADS();
    constexpr std::size_t kN = 300000;
    std::vector<std::uint64_t> keys(kN);
    std::uint64_t s = 0x9E3779B97F4A7C15ull;
    for (auto& k : keys) {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        k = s % 100000;
    }
    std::vector<std::uint64_t> ref_sorted;
    std::vector<std::int64_t> ref_scan(kN);
    double ref_sum = 0.0;
    for (int threads : {1, 2, 4, 3}) {
        PtCount count(threads);
        std::vector<std::pair<std::uint64_t, std::size_t>> v(kN);
        for (std::size_t i = 0; i < kN; ++i)
            v[i] = {keys[i], i};
        mio::parallel_sort(v.begin(), v.end());
        std::vector<std::uint64_t> sorted(kN);
        for (std::size_t i = 0; i < kN; ++i)
            sorted[i] = v[i].first * kN + v[i].second;
        std::vector<std::int64_t> scan(kN);
        mio::parallel_exclusive_scan(keys.data(), kN, scan.data(), std::int64_t{0});
        const double sum = mio::parallel_reduce(
            kN, 4096, 0.0,
            [&](std::size_t b, std::size_t e) {
                double acc = 0.0;
                for (std::size_t i = b; i < e; ++i)
                    acc += 1.0 / static_cast<double>(keys[i] + 1);
                return acc;
            },
            [](double x, double y) { return x + y; });
        if (ref_sorted.empty()) {
            ref_sorted = sorted;
            ref_scan = scan;
            ref_sum = sum;
            EXPECT_TRUE(std::is_sorted(sorted.begin(), sorted.end()));
        } else {
            EXPECT_EQ(sorted, ref_sorted) << threads;
            EXPECT_EQ(scan, ref_scan) << threads;
            EXPECT_EQ(sum, ref_sum) << threads;
        }
    }
}

#if !defined(_WIN32)
TEST(ParallelThreads, ForkedChildRebuildsThePool) {
    PT_REQUIRE_THREADS();
    PtCount count(4);
    std::atomic<std::size_t> warm{0};
    mio::parallel_for(100000, [&](std::size_t i) { warm += i & 1; }, 64);
    ASSERT_EQ(warm.load(), 50000u);
    const pid_t pid = fork();
    ASSERT_NE(pid, -1);
    if (pid == 0) {
        std::vector<std::uint64_t> v(100000, 0);
        mio::parallel_for(v.size(), [&](std::size_t i) { v[i] = i * 2; }, 64);
        for (std::size_t i = 0; i < v.size(); ++i)
            if (v[i] != i * 2)
                _exit(2);
        _exit(0);
    }
    int status = 0;
    ASSERT_EQ(waitpid(pid, &status, 0), pid);
    ASSERT_TRUE(WIFEXITED(status));
    EXPECT_EQ(WEXITSTATUS(status), 0);
}
#endif

}  // namespace
