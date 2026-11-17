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
// Unit tests for detail/hash_utils.hpp -- the hash combiner and the key hashes
// that replaced the hand-written loops (roadmap 6.3). The values must equal the
// loops' bit for bit: they set the iteration order of unordered containers.

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/hash_utils.hpp"

using namespace meshioplusplus::detail;

namespace {

std::size_t loop_hash(const std::vector<std::int64_t>& rKey, std::size_t Seed) {
    std::size_t h = Seed;
    for (std::int64_t v : rKey)
        h ^= std::hash<std::int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

}  // namespace

TEST(HashUtils, CombineIsTheBoostFormula) {
    const std::uint64_t h = 0x123456789abcdefull, v = 42;
    EXPECT_EQ(HashUtils::Combine<std::uint64_t>(h, v),
              h ^ (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2)));
}

TEST(HashUtils, ArrayHashMatchesTheLoopFromSeedZero) {
    for (std::int64_t a = -3; a < 40; a += 7)
        for (std::int64_t b = -2; b < 90; b += 11) {
            const std::array<std::int64_t, 2> k2{a, b};
            const std::array<std::int64_t, 4> k4{a, b, a * b, -a};
            EXPECT_EQ(ArrayHash<2>{}(k2), loop_hash({a, b}, 0));
            EXPECT_EQ(ArrayHash<4>{}(k4), loop_hash({a, b, a * b, -a}, 0));
        }
}

TEST(HashUtils, SizedVectorHashSeedsWithTheSize) {
    const std::vector<std::int64_t> k{5, 9, 1000000007, -4};
    EXPECT_EQ(SizedVectorHash{}(k), loop_hash(k, k.size()));
    EXPECT_NE(SizedVectorHash{}(k), ArrayHash<4>{}({5, 9, 1000000007, -4}));
}
