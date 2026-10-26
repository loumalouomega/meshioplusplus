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
// The maximum-weight matching of roadmap 6.8.1 (Edmonds, Galil 1986): every
// matching it returns is valid, and its weight (and, with MaxCardinality,
// its size first) equals a brute-force maximum on small random graphs.

// System includes
#include <algorithm>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/weighted_matching.hpp"

namespace {

namespace mdt = meshioplusplus::detail;

struct Score {
    std::int64_t mCount = 0;
    std::int64_t mWeight = 0;
};

/// The matching's size and weight, checking it is one: mates are mutual and
/// joined by an edge (the heaviest of parallel ones counted).
Score score(std::int64_t N, const std::vector<mdt::MatchingEdge>& rE,
            const std::vector<std::int64_t>& rMate) {
    Score s;
    EXPECT_EQ(static_cast<std::int64_t>(rMate.size()), N);
    for (std::int64_t v = 0; v < N; ++v) {
        const std::int64_t w = rMate[static_cast<std::size_t>(v)];
        if (w < 0 || w < v)
            continue;
        EXPECT_EQ(rMate[static_cast<std::size_t>(w)], v);
        std::int64_t best = 0;
        bool found = false;
        for (const auto& r_e : rE)
            if ((r_e.mA == v && r_e.mB == w) || (r_e.mA == w && r_e.mB == v)) {
                best = found ? std::max(best, r_e.mWeight) : r_e.mWeight;
                found = true;
            }
        EXPECT_TRUE(found) << v << " " << w;
        ++s.mCount;
        s.mWeight += best;
    }
    return s;
}

/// The best (count, weight) over every matching, by recursion on edges.
void brute(const std::vector<mdt::MatchingEdge>& rE, std::size_t K, std::vector<char>& rUsed,
           Score Now, bool MaxCard, Score& rBest) {
    if (K == rE.size()) {
        const bool better = MaxCard ? (Now.mCount > rBest.mCount ||
                                       (Now.mCount == rBest.mCount && Now.mWeight > rBest.mWeight))
                                    : Now.mWeight > rBest.mWeight;
        if (better)
            rBest = Now;
        return;
    }
    brute(rE, K + 1, rUsed, Now, MaxCard, rBest);
    const auto& r_e = rE[K];
    if (r_e.mA != r_e.mB && !rUsed[static_cast<std::size_t>(r_e.mA)] &&
        !rUsed[static_cast<std::size_t>(r_e.mB)]) {
        rUsed[static_cast<std::size_t>(r_e.mA)] = rUsed[static_cast<std::size_t>(r_e.mB)] = 1;
        brute(rE, K + 1, rUsed, {Now.mCount + 1, Now.mWeight + r_e.mWeight}, MaxCard, rBest);
        rUsed[static_cast<std::size_t>(r_e.mA)] = rUsed[static_cast<std::size_t>(r_e.mB)] = 0;
    }
}

}  // namespace

TEST(WeightedMatching, SmallCases) {
    EXPECT_TRUE(mdt::max_weight_matching(0, {}, false).empty());
    EXPECT_EQ(mdt::max_weight_matching(2, {{0, 1, 1}}, false), (std::vector<std::int64_t>{1, 0}));
    // A path of three edges: the outer two weigh more together.
    EXPECT_EQ(mdt::max_weight_matching(4, {{0, 1, 5}, {1, 2, 11}, {2, 3, 5}}, false),
              (std::vector<std::int64_t>{-1, 2, 1, -1}));
    EXPECT_EQ(mdt::max_weight_matching(4, {{0, 1, 5}, {1, 2, 11}, {2, 3, 5}}, true),
              (std::vector<std::int64_t>{1, 0, 3, 2}));
    // Negative weights are never taken unless cardinality asks.
    EXPECT_EQ(mdt::max_weight_matching(2, {{0, 1, -2}}, false),
              (std::vector<std::int64_t>{-1, -1}));
    EXPECT_EQ(mdt::max_weight_matching(2, {{0, 1, -2}}, true), (std::vector<std::int64_t>{1, 0}));
    // A triangle with a pendant: the blossom {0, 1, 2} is shrunk and expanded.
    EXPECT_EQ(mdt::max_weight_matching(4, {{0, 1, 8}, {0, 2, 9}, {1, 2, 10}, {2, 3, 7}}, false),
              (std::vector<std::int64_t>{1, 0, 3, 2}));
}

// Nested blossoms, S-blossoms relabelled as T, and blossoms expanded in the
// middle of a stage: the classic hard cases of the primal-dual algorithm.
TEST(WeightedMatching, NestedAndExpandedBlossoms) {
    // Create a nested S-blossom, use it for augmentation.
    EXPECT_EQ(
        mdt::max_weight_matching(
            6, {{0, 1, 9}, {0, 2, 9}, {1, 2, 10}, {1, 3, 8}, {2, 4, 8}, {3, 4, 10}, {4, 5, 6}},
            false),
        (std::vector<std::int64_t>{2, 3, 0, 1, 5, 4}));
    // An S-blossom relabelled as a T-blossom.
    EXPECT_EQ(mdt::max_weight_matching(
                  8, {{0, 1, 9}, {0, 2, 8}, {1, 2, 10}, {0, 3, 5}, {3, 4, 4}, {0, 5, 3}}, false),
              (std::vector<std::int64_t>{5, 2, 1, 4, 3, 0, -1, -1}));
    // A nested S-blossom, relabelled as T, expanded.
    EXPECT_EQ(mdt::max_weight_matching(8,
                                       {{0, 1, 19},
                                        {0, 2, 20},
                                        {0, 7, 8},
                                        {1, 2, 25},
                                        {1, 3, 18},
                                        {2, 4, 18},
                                        {3, 4, 13},
                                        {3, 6, 7},
                                        {4, 5, 7}},
                                       false),
              (std::vector<std::int64_t>{7, 2, 1, 6, 5, 4, 3, 0}));
}

TEST(WeightedMatching, AgreesWithBruteForce) {
    std::mt19937_64 rng(20261026);
    for (int trial = 0; trial < 600; ++trial) {
        const std::int64_t n = 2 + static_cast<std::int64_t>(rng() % 9);
        const std::size_t m = 1 + rng() % 14;
        std::vector<mdt::MatchingEdge> e;
        for (std::size_t k = 0; k < m; ++k) {
            const auto a = static_cast<std::int64_t>(rng() % static_cast<std::uint64_t>(n));
            const auto b = static_cast<std::int64_t>(rng() % static_cast<std::uint64_t>(n));
            const auto w = static_cast<std::int64_t>(rng() % 21) - (trial % 3 == 0 ? 5 : 0);
            e.push_back({a, b, w});
        }
        for (const bool max_card : {false, true}) {
            const auto mate = mdt::max_weight_matching(n, e, max_card);
            const Score got = score(n, e, mate);
            Score best;
            std::vector<char> used(static_cast<std::size_t>(n), 0);
            brute(e, 0, used, {}, max_card, best);
            if (max_card)
                EXPECT_EQ(got.mCount, best.mCount) << trial;
            EXPECT_EQ(got.mWeight, best.mWeight) << trial << " " << max_card;
        }
    }
}
