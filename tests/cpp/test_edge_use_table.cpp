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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// Unit tests for detail/edge_use_table.hpp (roadmap 6.2): the table answers what
// a std::map keyed by the sorted vertex pair did, in the map's key order.

// System includes
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/edge_use_table.hpp"

using meshioplusplus::detail::EdgeUseTable;

namespace {

// A closed tetrahedron: every edge used twice, once each way.
EdgeUseTable<std::int64_t> tetra() {
    EdgeUseTable<std::int64_t> t;
    const std::int64_t faces[4][3] = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
    for (std::int64_t f = 0; f < 4; ++f)
        for (int k = 0; k < 3; ++k)
            t.Add(faces[f][k], faces[f][(k + 1) % 3], f);
    t.Build();
    return t;
}

}  // namespace

TEST(EdgeUseTable, EmptyHasNoEdges) {
    EdgeUseTable<std::size_t> t;
    t.Build();
    EXPECT_EQ(t.NumEdges(), 0u);
    EXPECT_EQ(t.CountOf(1, 2), 0u);
}

TEST(EdgeUseTable, ClosedSurfaceUsesEveryEdgeTwiceOppositely) {
    const auto t = tetra();
    ASSERT_EQ(t.NumEdges(), 6u);
    for (std::size_t e = 0; e < t.NumEdges(); ++e) {
        ASSERT_EQ(t.Count(e), 2u);
        EXPECT_NE(t.At(e, 0).mForward, t.At(e, 1).mForward);
    }
}

TEST(EdgeUseTable, EdgesComeInKeyOrderLikeAMap) {
    const auto t = tetra();
    std::map<std::pair<std::int64_t, std::int64_t>, int> want;
    const std::int64_t faces[4][3] = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
    for (const auto& f : faces)
        for (int k = 0; k < 3; ++k)
            ++want[{std::min(f[k], f[(k + 1) % 3]), std::max(f[k], f[(k + 1) % 3])}];
    std::size_t e = 0;
    for (const auto& [edge, count] : want) {
        EXPECT_EQ(t.Lo(e), edge.first);
        EXPECT_EQ(t.Hi(e), edge.second);
        EXPECT_EQ(static_cast<int>(t.Count(e)), count);
        ++e;
    }
}

TEST(EdgeUseTable, FindFirstUseAndBoundary) {
    EdgeUseTable<std::int64_t> t;
    t.Add(5, 3, 10);
    t.Add(3, 5, 11);
    t.Add(5, 9, 12);
    t.Build();
    EXPECT_EQ(t.CountOf(3, 5), 2u);
    EXPECT_EQ(t.CountOf(9, 5), 1u);
    EXPECT_EQ(t.CountOf(1, 2), 0u);
    const std::size_t e = t.Find(5, 3);
    ASSERT_NE(e, EdgeUseTable<std::int64_t>::kNone);
    EXPECT_EQ(t.At(e, 0).mFace, 10);
    EXPECT_FALSE(t.At(e, 0).mForward);
    EXPECT_EQ(t.At(e, 1).mFace, 11);
    EXPECT_TRUE(t.At(e, 1).mForward);
}

TEST(EdgeUseTable, ForEachSharedFaceJoinsLaterUsesToTheFirst) {
    EdgeUseTable<std::size_t> t;
    t.Add(0, 1, 7);
    t.Add(1, 0, 8);
    t.Add(0, 1, 9);
    t.Add(2, 3, 4);
    t.Build();
    std::vector<std::pair<std::size_t, std::size_t>> joins;
    t.ForEachSharedFace([&](std::size_t A, std::size_t B) { joins.push_back({A, B}); });
    ASSERT_EQ(joins.size(), 2u);
    EXPECT_EQ(joins[0], (std::pair<std::size_t, std::size_t>{7, 8}));
    EXPECT_EQ(joins[1], (std::pair<std::size_t, std::size_t>{7, 9}));
}
