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
/**
 * @file test_slot_runs.cpp
 * @brief The private sort-based table (`src/cpp/src/detail/slot_runs.hpp`)
 *        against plain serial references, on inputs large enough to take the
 *        chunked counting sort.
 */

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/slot_runs.hpp"

namespace {

namespace det = meshioplusplus::detail;
using Key = std::array<std::int64_t, 2>;

std::vector<Key> sr_keys(std::size_t n, std::int64_t range) {
    std::vector<Key> keys(n);
    std::uint64_t s = 88172645463325252ull;
    for (Key& k : keys) {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        const auto a = static_cast<std::int64_t>(s % static_cast<std::uint64_t>(range));
        const auto b = static_cast<std::int64_t>((s >> 32) % static_cast<std::uint64_t>(range));
        k = {std::min(a, b), std::max(a, b)};
    }
    return keys;
}

// The reference: a stable sort of the slots by (bucket, key), cut into runs.
det::SlotRuns sr_reference(const std::vector<Key>& rKeys, std::size_t NumBuckets) {
    const auto bucket = [&](std::uint64_t i) {
        return static_cast<std::size_t>(rKeys[i][0]) % NumBuckets;
    };
    std::vector<std::uint64_t> order(rKeys.size());
    std::iota(order.begin(), order.end(), std::uint64_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::uint64_t a, std::uint64_t b) {
        if (bucket(a) != bucket(b))
            return bucket(a) < bucket(b);
        return rKeys[a] < rKeys[b];
    });
    det::SlotRuns out;
    out.mStart.clear();
    for (std::size_t i = 0; i < order.size(); ++i)
        if (i == 0 || rKeys[order[i - 1]] != rKeys[order[i]])
            out.mStart.push_back(i);
    out.mStart.push_back(order.size());
    out.mSlots = std::move(order);
    return out;
}

}  // namespace

TEST(SlotRuns, ChunkedCountingSortMatchesAStableSort) {
    for (const std::size_t buckets : {std::size_t{97}, std::size_t{5000}, std::size_t{400000}}) {
        const std::vector<Key> keys = sr_keys(200000, 60000);
        const det::SlotRuns got = det::group_slots(keys, buckets, [&](const Key& rK) {
            return static_cast<std::size_t>(rK[0]) % buckets;
        });
        const det::SlotRuns want = sr_reference(keys, buckets);
        EXPECT_EQ(got.mSlots, want.mSlots) << buckets;
        EXPECT_EQ(got.mStart, want.mStart) << buckets;
    }
}

TEST(SlotRuns, GroupIntRowsKeepsTheFirstOfEachRow) {
    // Rows of 1..5 values from a small range, so many repeat.
    std::vector<std::int64_t> values;
    std::vector<std::uint64_t> offsets{0};
    std::uint64_t s = 12345;
    for (std::size_t r = 0; r < 60000; ++r) {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        const std::size_t len = 1 + (s >> 60) % 5;
        std::vector<std::int64_t> row;
        for (std::size_t k = 0; k < len; ++k) {
            s = s * 6364136223846793005ull + 1442695040888963407ull;
            row.push_back(static_cast<std::int64_t>((s >> 33) % 12));
        }
        std::sort(row.begin(), row.end());
        values.insert(values.end(), row.begin(), row.end());
        offsets.push_back(values.size());
    }
    const std::size_t n = offsets.size() - 1;
    const det::SlotRuns runs = det::group_int_rows(values, offsets, 12);
    const std::vector<std::uint8_t> dup = det::later_duplicates(runs, n);
    std::map<std::vector<std::int64_t>, int> seen;
    for (std::size_t r = 0; r < n; ++r) {
        const std::vector<std::int64_t> row(
            values.begin() + static_cast<std::ptrdiff_t>(offsets[r]),
            values.begin() + static_cast<std::ptrdiff_t>(offsets[r + 1]));
        const bool first = seen.emplace(row, 0).second;
        ASSERT_EQ(dup[r] == 0, first) << r;
    }
}
