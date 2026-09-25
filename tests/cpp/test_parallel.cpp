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
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/parallel.hpp"

TEST(Parallel, BackendName) {
    std::string name = meshioplusplus::parallel_backend_name();
    EXPECT_TRUE(name == "seq" || name == "stl" || name == "openmp" || name == "tbb" ||
                name == "kokkos")
        << name;
}

TEST(Parallel, EmptyAndSingle) {
    std::atomic<int> count{0};
    meshioplusplus::parallel_for(0, [&](std::size_t) { ++count; });
    EXPECT_EQ(count.load(), 0);
    meshioplusplus::parallel_for(1, [&](std::size_t) { ++count; });
    EXPECT_EQ(count.load(), 1);
}

TEST(Parallel, ScatterSmallAndLarge) {
    // Small (below the grain -> sequential path) and large (parallel path):
    // every index written exactly once, correct values.
    for (std::size_t n : {std::size_t{17}, std::size_t{100000}}) {
        std::vector<std::int64_t> out(n, -1);
        meshioplusplus::parallel_for(
            n, [&](std::size_t i) { out[i] = static_cast<std::int64_t>(i) * 3; });
        for (std::size_t i = 0; i < n; ++i)
            ASSERT_EQ(out[i], static_cast<std::int64_t>(i) * 3) << "i=" << i;
    }
}

TEST(Parallel, EveryIndexExactlyOnce) {
    const std::size_t n = 50000;
    std::vector<std::atomic<int>> hits(n);
    meshioplusplus::parallel_for(n, [&](std::size_t i) { ++hits[i]; });
    for (std::size_t i = 0; i < n; ++i)
        ASSERT_EQ(hits[i].load(), 1) << "i=" << i;
}

TEST(Parallel, CustomGrain) {
    const std::size_t n = 10000;
    std::vector<int> out(n, 0);
    meshioplusplus::parallel_for(n, [&](std::size_t i) { out[i] = 1; }, /*grain=*/64);
    EXPECT_EQ(std::accumulate(out.begin(), out.end(), 0), static_cast<int>(n));
}

TEST(Parallel, ExceptionPropagates) {
    const std::size_t n = 100000;  // above the grain -> parallel path
    EXPECT_THROW(meshioplusplus::parallel_for(n,
                                              [&](std::size_t i) {
                                                  if (i == n / 2)
                                                      throw std::runtime_error("boom");
                                              }),
                 std::runtime_error);
    // Sequential path (n below grain) propagates too.
    EXPECT_THROW(
        meshioplusplus::parallel_for(10, [&](std::size_t) { throw std::runtime_error("boom"); }),
        std::runtime_error);
}

namespace {

std::uint64_t par_xorshift(std::uint64_t& rState) {
    rState ^= rState << 13;
    rState ^= rState >> 7;
    rState ^= rState << 17;
    return rState;
}

}  // namespace

// (key, unique slot) records, heavy key ties: every backend must return the
// one sorted sequence std::sort gives, above and below the serial threshold.
TEST(Parallel, SortOfATotalOrderMatchesStdSort) {
    for (const std::size_t n : {std::size_t{0}, std::size_t{1}, std::size_t{1000},
                                std::size_t{200000}, std::size_t{1000003}}) {
        std::uint64_t state = 0x9E3779B97F4A7C15ull + n;
        std::vector<std::array<std::int64_t, 3>> recs(n);
        for (std::size_t i = 0; i < n; ++i)
            recs[i] = {static_cast<std::int64_t>(par_xorshift(state) % 97),
                       static_cast<std::int64_t>(par_xorshift(state) % 5),
                       static_cast<std::int64_t>(i)};
        auto expect = recs;
        std::sort(expect.begin(), expect.end());
        meshioplusplus::parallel_sort(recs.begin(), recs.end());
        EXPECT_EQ(recs, expect) << n;
    }
}

TEST(Parallel, SortTakesAComparator) {
    std::vector<std::int64_t> v(300000);
    std::iota(v.begin(), v.end(), std::int64_t{0});
    meshioplusplus::parallel_sort(v.begin(), v.end(), std::greater<std::int64_t>());
    for (std::size_t i = 0; i < v.size(); ++i)
        ASSERT_EQ(v[i], static_cast<std::int64_t>(v.size() - 1 - i));
}

TEST(Parallel, ExclusiveScanMatchesTheSerialScan) {
    for (const std::size_t n : {std::size_t{0}, std::size_t{1}, std::size_t{32767},
                                std::size_t{32768}, std::size_t{1000001}}) {
        std::uint64_t state = 12345 + n;
        std::vector<std::uint8_t> flags(n);
        for (auto& f : flags)
            f = static_cast<std::uint8_t>(par_xorshift(state) % 3);
        std::vector<std::int64_t> expect(n);
        std::int64_t acc = 7;
        for (std::size_t i = 0; i < n; ++i) {
            expect[i] = acc;
            acc += flags[i];
        }
        std::vector<std::int64_t> got(n);
        EXPECT_EQ(
            meshioplusplus::parallel_exclusive_scan(flags.data(), n, got.data(), std::int64_t{7}),
            acc);
        EXPECT_EQ(got, expect) << n;
        // In place.
        std::vector<std::int64_t> inplace(flags.begin(), flags.end());
        meshioplusplus::parallel_exclusive_scan(inplace.data(), n, inplace.data(), std::int64_t{7});
        EXPECT_EQ(inplace, expect) << n;
    }
}

// Chunk partials combined in chunk order: the same bits as a hand-written
// chunked loop, whatever the backend and thread count.
TEST(Parallel, ReduceCombinesFixedChunksInOrder) {
    const std::size_t n = 1000003, chunk = 4096;
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i)
        v[i] = 1.0 / static_cast<double>(i + 1);
    double expect = 0.25;
    for (std::size_t b = 0; b < n; b += chunk) {
        double part = 0.0;
        for (std::size_t i = b; i < std::min(n, b + chunk); ++i)
            part += v[i];
        expect += part;
    }
    const double got = meshioplusplus::parallel_reduce(
        n, chunk, 0.25,
        [&](std::size_t b, std::size_t e) {
            double part = 0.0;
            for (std::size_t i = b; i < e; ++i)
                part += v[i];
            return part;
        },
        [](double a, double p) { return a + p; });
    EXPECT_EQ(got, expect);
    EXPECT_EQ(meshioplusplus::parallel_reduce(
                  0, chunk, 3, [](std::size_t, std::size_t) { return 1; },
                  [](int a, int p) { return a + p; }),
              3);
}
