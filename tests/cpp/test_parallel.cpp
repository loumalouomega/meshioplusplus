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
