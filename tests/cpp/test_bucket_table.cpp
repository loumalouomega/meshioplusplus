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
#include <cstdint>
#include <map>
#include <tuple>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/detail/bucket_table.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"

namespace {

using meshioplusplus::detail::BucketTable;
using meshioplusplus::detail::BucketView;
using meshioplusplus::detail::GridKey;
using meshioplusplus::detail::SpatialGrid;

struct Box {
    GridKey mLo;
    GridKey mHi;
};

using Visits = std::vector<std::vector<std::int64_t>>;

// A small deterministic generator: the tests must not depend on the platform's
// `std::uniform_int_distribution`.
struct Lcg {
    std::uint64_t mState;
    std::int64_t Next(std::int64_t Lo, std::int64_t Hi) {  // inclusive
        mState = mState * 6364136223846793005ULL + 1442695040888963407ULL;
        const std::uint64_t span = static_cast<std::uint64_t>(Hi - Lo) + 1;
        return Lo + static_cast<std::int64_t>((mState >> 33) % span);
    }
};

// What `SpatialGrid::InsertBox` leaves for boxes inserted serially, one id per
// box, in the order of `rBoxes`.
SpatialGrid serial_grid(const std::vector<Box>& rBoxes) {
    SpatialGrid grid(1.0);
    for (std::size_t i = 0; i < rBoxes.size(); ++i)
        grid.InsertBox(rBoxes[i].mLo, rBoxes[i].mHi, static_cast<std::int64_t>(i));
    return grid;
}

// The same buckets as a `BucketTable`: keys in the order the serial inserts
// first met them, ids ascending, the occupied box the union of the boxes.
BucketTable frozen_grid(const std::vector<Box>& rBoxes, double Cell = 1.0) {
    std::map<std::tuple<std::int64_t, std::int64_t, std::int64_t>, std::size_t> bucket_of;
    std::vector<GridKey> keys;
    std::vector<std::vector<std::int64_t>> members;
    GridKey lo = rBoxes.front().mLo;
    GridKey hi = rBoxes.front().mHi;
    for (std::size_t i = 0; i < rBoxes.size(); ++i) {
        const Box& b = rBoxes[i];
        lo = GridKey{std::min(lo.x, b.mLo.x), std::min(lo.y, b.mLo.y), std::min(lo.z, b.mLo.z)};
        hi = GridKey{std::max(hi.x, b.mHi.x), std::max(hi.y, b.mHi.y), std::max(hi.z, b.mHi.z)};
        for (std::int64_t z = b.mLo.z; z <= b.mHi.z; ++z)
            for (std::int64_t y = b.mLo.y; y <= b.mHi.y; ++y)
                for (std::int64_t x = b.mLo.x; x <= b.mHi.x; ++x) {
                    const auto key = std::make_tuple(x, y, z);
                    auto it = bucket_of.find(key);
                    if (it == bucket_of.end()) {
                        it = bucket_of.emplace(key, keys.size()).first;
                        keys.push_back(GridKey{x, y, z});
                        members.emplace_back();
                    }
                    members[it->second].push_back(static_cast<std::int64_t>(i));
                }
    }
    std::vector<std::uint64_t> offsets{0};
    std::vector<std::int64_t> ids;
    for (const auto& rMembers : members) {
        ids.insert(ids.end(), rMembers.begin(), rMembers.end());
        offsets.push_back(ids.size());
    }
    return BucketTable(Cell, std::move(keys), std::move(offsets), std::move(ids), lo, hi);
}

std::vector<Box> random_boxes(Lcg& rRng, std::size_t Count, std::int64_t Lo, std::int64_t Hi,
                              std::int64_t MaxSpan, std::int64_t Shift = 0) {
    std::vector<Box> boxes;
    for (std::size_t i = 0; i < Count; ++i) {
        GridKey lo{rRng.Next(Lo, Hi) + Shift, rRng.Next(Lo, Hi) + Shift, rRng.Next(Lo, Hi) + Shift};
        GridKey hi{lo.x + rRng.Next(0, MaxSpan), lo.y + rRng.Next(0, MaxSpan),
                   lo.z + rRng.Next(0, MaxSpan)};
        boxes.push_back({lo, hi});
    }
    return boxes;
}

template <class Grid>
Visits shell_visits(const Grid& rGrid, const GridKey& rCentre, std::int64_t R) {
    Visits out;
    rGrid.ForEachInShell(rCentre, R, [&](const auto& rIds) {
        out.emplace_back(rIds.begin(), rIds.end());
    });
    return out;
}

template <class Grid>
Visits box_visits(const Grid& rGrid, const GridKey& rLo, const GridKey& rHi) {
    Visits out;
    rGrid.ForEachInBox(rLo, rHi, [&](const auto& rIds) { out.emplace_back(rIds.begin(), rIds.end()); });
    return out;
}

// Every shell and box query of a sweep around the boxes returns exactly what
// the serial `SpatialGrid` does, in the same order.
void expect_same_traversals(const std::vector<Box>& rBoxes, std::int64_t Lo, std::int64_t Hi,
                            std::int64_t Shift, Lcg& rRng) {
    const SpatialGrid reference = serial_grid(rBoxes);
    const BucketTable table = frozen_grid(rBoxes);
    ASSERT_FALSE(table.Empty());
    EXPECT_EQ(table.OccupiedLo().x, reference.OccupiedLo().x);
    EXPECT_EQ(table.OccupiedLo().y, reference.OccupiedLo().y);
    EXPECT_EQ(table.OccupiedLo().z, reference.OccupiedLo().z);
    EXPECT_EQ(table.OccupiedHi().x, reference.OccupiedHi().x);
    EXPECT_EQ(table.OccupiedHi().y, reference.OccupiedHi().y);
    EXPECT_EQ(table.OccupiedHi().z, reference.OccupiedHi().z);

    for (int q = 0; q < 200; ++q) {
        // Centres inside and well outside the occupied box.
        const GridKey centre{rRng.Next(Lo - 6, Hi + 6) + Shift, rRng.Next(Lo - 6, Hi + 6) + Shift,
                             rRng.Next(Lo - 6, Hi + 6) + Shift};
        for (std::int64_t r = 0; r <= 9; ++r)
            ASSERT_EQ(shell_visits(table, centre, r), shell_visits(reference, centre, r))
                << "centre " << centre.x << ',' << centre.y << ',' << centre.z << " r " << r;
    }
    for (int q = 0; q < 100; ++q) {
        const GridKey lo{rRng.Next(Lo - 3, Hi) + Shift, rRng.Next(Lo - 3, Hi) + Shift,
                         rRng.Next(Lo - 3, Hi) + Shift};
        const GridKey hi{lo.x + rRng.Next(0, 4), lo.y + rRng.Next(0, 4), lo.z + rRng.Next(0, 4)};
        ASSERT_EQ(box_visits(table, lo, hi), box_visits(reference, lo, hi));
    }
}

}  // namespace

TEST(BucketTable, EmptyTableVisitsNothing) {
    const BucketTable table(2.5);
    EXPECT_TRUE(table.Empty());
    EXPECT_EQ(table.NumBuckets(), 0u);
    EXPECT_EQ(table.CellSize(), 2.5);
    EXPECT_TRUE(shell_visits(table, GridKey{0, 0, 0}, 0).empty());
    EXPECT_TRUE(shell_visits(table, GridKey{0, 0, 0}, 3).empty());
    EXPECT_TRUE(box_visits(table, GridKey{-2, -2, -2}, GridKey{2, 2, 2}).empty());
    BucketView view(nullptr, 0);
    EXPECT_FALSE(table.Find(GridKey{0, 0, 0}, view));
}

TEST(BucketTable, KeyOfQuantizesLikeSpatialGrid) {
    const BucketTable table(0.5);
    const SpatialGrid grid(0.5);
    const double points[][3] = {{0, 0, 0}, {0.49, 0.5, -0.01}, {-3.2, 7.77, 1e3}, {-0.5, -0.5, -0.5}};
    for (const auto& rP : points) {
        const GridKey a = table.KeyOf(rP);
        const GridKey b = grid.KeyOf(rP);
        EXPECT_EQ(a.x, b.x);
        EXPECT_EQ(a.y, b.y);
        EXPECT_EQ(a.z, b.z);
    }
}

TEST(BucketTable, SingleBucketIsFoundAndNothingElse) {
    const std::vector<Box> boxes{{GridKey{3, -4, 5}, GridKey{3, -4, 5}}};
    const BucketTable table = frozen_grid(boxes);
    ASSERT_EQ(table.NumBuckets(), 1u);
    BucketView view(nullptr, 0);
    ASSERT_TRUE(table.Find(GridKey{3, -4, 5}, view));
    ASSERT_EQ(view.size(), 1u);
    EXPECT_EQ(view[0], 0);
    EXPECT_FALSE(table.Find(GridKey{3, -4, 6}, view));
    EXPECT_FALSE(table.Find(GridKey{-3, 4, -5}, view));
    // The cell is at Chebyshev distance 5 from the origin: the shells before
    // it are empty and the one at 5 holds it.
    EXPECT_EQ(shell_visits(table, GridKey{0, 0, 0}, 4), (Visits{}));
    EXPECT_EQ(shell_visits(table, GridKey{0, 0, 0}, 5), (Visits{{0}}));
    EXPECT_EQ(shell_visits(table, GridKey{3, -4, 5}, 0), (Visits{{0}}));
}

TEST(BucketTable, ShellsAndBoxesMatchTheSerialGrid) {
    Lcg rng{1};
    const auto boxes = random_boxes(rng, 60, -6, 6, 3);
    expect_same_traversals(boxes, -6, 9, 0, rng);
}

TEST(BucketTable, ManyBucketsExerciseTheProbeSequence) {
    // Thousands of cells in a tight region: the index is at most half full and
    // collides often, and every present key must still be found.
    Lcg rng{7};
    const auto boxes = random_boxes(rng, 400, 0, 12, 2);
    const SpatialGrid reference = serial_grid(boxes);
    const BucketTable table = frozen_grid(boxes);
    EXPECT_GT(table.NumBuckets(), 1000u);
    for (std::int64_t z = -2; z <= 16; ++z)
        for (std::int64_t y = -2; y <= 16; ++y)
            for (std::int64_t x = -2; x <= 16; ++x) {
                const GridKey key{x, y, z};
                const std::vector<std::int64_t>* expected = reference.Find(key);
                BucketView got(nullptr, 0);
                const bool found = table.Find(key, got);
                ASSERT_EQ(found, expected != nullptr) << x << ',' << y << ',' << z;
                if (found)
                    ASSERT_EQ(std::vector<std::int64_t>(got.begin(), got.end()), *expected);
            }
    expect_same_traversals(boxes, 0, 14, 0, rng);
}

TEST(BucketTable, NegativeAndLargeKeys) {
    Lcg rng{3};
    const std::int64_t far = std::int64_t{1} << 40;
    expect_same_traversals(random_boxes(rng, 40, -4, 4, 2, -far), -4, 6, -far, rng);
    expect_same_traversals(random_boxes(rng, 40, -4, 4, 2, far), -4, 6, far, rng);
}

TEST(BucketTable, BucketsKeepTheirIdsAscendingAndInTheCallersOrder) {
    // The ids of a bucket, and the order buckets were handed over in, are what
    // a nearest-triangle search's tie-break rests on.
    const std::vector<Box> boxes{{GridKey{0, 0, 0}, GridKey{1, 0, 0}},
                                 {GridKey{1, 0, 0}, GridKey{2, 0, 0}},
                                 {GridKey{1, 0, 0}, GridKey{1, 0, 0}}};
    const BucketTable table = frozen_grid(boxes);
    BucketView view(nullptr, 0);
    ASSERT_TRUE(table.Find(GridKey{1, 0, 0}, view));
    EXPECT_EQ(std::vector<std::int64_t>(view.begin(), view.end()), (std::vector<std::int64_t>{0, 1, 2}));
    ASSERT_TRUE(table.Find(GridKey{0, 0, 0}, view));
    EXPECT_EQ(std::vector<std::int64_t>(view.begin(), view.end()), (std::vector<std::int64_t>{0}));
    ASSERT_TRUE(table.Find(GridKey{2, 0, 0}, view));
    EXPECT_EQ(std::vector<std::int64_t>(view.begin(), view.end()), (std::vector<std::int64_t>{1}));
    EXPECT_EQ(table.NumBuckets(), 3u);
}
