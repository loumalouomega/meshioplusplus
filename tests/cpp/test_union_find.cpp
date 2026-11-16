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
// Unit tests for detail/union_find.hpp (roadmap 6.2): the smaller root wins, so
// the partition's representatives do not depend on the order of the joins.

// System includes
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/union_find.hpp"

using meshioplusplus::detail::UnionFind;

TEST(UnionFind, RootIsTheSmallestIndexWhateverTheJoinOrder) {
    const std::vector<std::pair<std::size_t, std::size_t>> joins{{5, 2}, {7, 5}, {1, 9}, {9, 4},
                                                                  {3, 8}, {8, 6}};
    std::vector<std::size_t> roots_a, roots_b;
    for (int pass = 0; pass < 2; ++pass) {
        auto order = joins;
        if (pass == 1)
            std::shuffle(order.begin(), order.end(), std::mt19937(5));
        UnionFind<std::size_t> uf(10);
        for (const auto& j : order)
            uf.Unite(j.first, j.second);
        auto& out = pass == 0 ? roots_a : roots_b;
        for (std::size_t i = 0; i < 10; ++i)
            out.push_back(uf.Find(i));
    }
    EXPECT_EQ(roots_a, roots_b);
    EXPECT_EQ(roots_a[7], 2u);
    EXPECT_EQ(roots_a[9], 1u);
    EXPECT_EQ(roots_a[0], 0u);
}

TEST(UnionFind, UniteReportsWhetherTheSetsWereApart) {
    UnionFind<std::int64_t> uf(4);
    EXPECT_TRUE(uf.Unite(3, 1));
    EXPECT_FALSE(uf.Unite(1, 3));
    EXPECT_EQ(uf.Find(3), 1);
    EXPECT_LE(uf.Parent(3), 3);
}

TEST(UnionFind, FindInHalvesAForeignParentVector) {
    std::vector<std::size_t> parent{0, 0, 1, 2, 3};
    EXPECT_EQ((UnionFind<std::size_t>::FindIn(parent, 4)), 0u);
    EXPECT_EQ(parent[4], 2u);
}
