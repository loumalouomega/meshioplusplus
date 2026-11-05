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
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes
#include "bucket_build.hpp"

namespace meshioplusplus {
namespace detail {

// One serial sweep over the (id, z, y, x) pairs -- the serial `InsertBox`
// loop's order -- numbering each cell the first time it is met through a flat
// open-addressing table, and counting its ids; then a second sweep places
// every id at its bucket's next slot. The buckets come out in first-seen
// order with their ids ascending: exactly what the serial inserts left in
// `SpatialGrid`, with one probe per pair and no allocation per cell. A
// parallel grouping (`group_slots`) gave the same table, but on one thread it
// cost more than the hash-map inserts it replaced, which the landing gate
// ("SEQ not slower") does not allow (roadmap 5.4.1).
BucketTable build_bucket_table(double Cell, const std::vector<GridKey>& rLo,
                               const std::vector<GridKey>& rHi) {
    const std::size_t n = rLo.size();
    if (n == 0)
        return BucketTable(Cell);
    std::uint64_t npairs = 0;
    for (std::size_t t = 0; t < n; ++t) {
        const std::int64_t nx = rHi[t].x - rLo[t].x + 1;
        const std::int64_t ny = rHi[t].y - rLo[t].y + 1;
        const std::int64_t nz = rHi[t].z - rLo[t].z + 1;
        if (nx > 0 && ny > 0 && nz > 0)
            npairs += static_cast<std::uint64_t>(nx) * static_cast<std::uint64_t>(ny) *
                      static_cast<std::uint64_t>(nz);
    }
    // A power of two at least twice the most buckets there can be.
    std::size_t cap = 16;
    while (cap < 2 * npairs)
        cap <<= 1;
    const std::size_t mask = cap - 1;
    constexpr std::uint64_t kEmpty = ~std::uint64_t{0};
    std::vector<std::uint64_t> slot(cap, kEmpty);
    std::vector<GridKey> bucket_keys;
    std::vector<std::uint64_t> bucket_count;
    std::vector<std::uint64_t> bucket_of_pair(static_cast<std::size_t>(npairs));
    const GridKeyHash hash;
    std::size_t pair = 0;
    for (std::size_t t = 0; t < n; ++t)
        for (std::int64_t z = rLo[t].z; z <= rHi[t].z; ++z)
            for (std::int64_t y = rLo[t].y; y <= rHi[t].y; ++y)
                for (std::int64_t x = rLo[t].x; x <= rHi[t].x; ++x) {
                    const GridKey key{x, y, z};
                    // GridKeyHash's low bits cluster for neighbouring cells;
                    // the splitmix64 finalizer spreads them over the mask.
                    std::uint64_t h = static_cast<std::uint64_t>(hash(key));
                    h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
                    h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
                    h ^= h >> 31;
                    std::size_t at = static_cast<std::size_t>(h) & mask;
                    for (;; at = (at + 1) & mask) {
                        if (slot[at] == kEmpty) {
                            slot[at] = bucket_keys.size();
                            bucket_keys.push_back(key);
                            bucket_count.push_back(0);
                            break;
                        }
                        if (bucket_keys[static_cast<std::size_t>(slot[at])] == key)
                            break;
                    }
                    ++bucket_count[static_cast<std::size_t>(slot[at])];
                    bucket_of_pair[pair++] = slot[at];
                }
    std::vector<std::uint64_t>().swap(slot);
    const std::size_t nb = bucket_keys.size();
    std::vector<std::uint64_t> bucket_offset(nb + 1);
    bucket_offset[0] = 0;
    for (std::size_t b = 0; b < nb; ++b)
        bucket_offset[b + 1] = bucket_offset[b] + bucket_count[b];
    std::vector<std::uint64_t>& next = bucket_count;  // reused as the fill cursor
    std::copy(bucket_offset.begin(), bucket_offset.end() - 1, next.begin());
    std::vector<std::int64_t> bucket_ids(static_cast<std::size_t>(npairs));
    pair = 0;
    for (std::size_t t = 0; t < n; ++t) {
        const std::int64_t nx = rHi[t].x - rLo[t].x + 1;
        const std::int64_t ny = rHi[t].y - rLo[t].y + 1;
        const std::int64_t nz = rHi[t].z - rLo[t].z + 1;
        if (nx <= 0 || ny <= 0 || nz <= 0)
            continue;
        const std::uint64_t count = static_cast<std::uint64_t>(nx) *
                                    static_cast<std::uint64_t>(ny) * static_cast<std::uint64_t>(nz);
        for (std::uint64_t k = 0; k < count; ++k)
            bucket_ids[static_cast<std::size_t>(next[bucket_of_pair[pair++]]++)] =
                static_cast<std::int64_t>(t);
    }
    // The occupied box: InsertBox covers every box's low and high key.
    GridKey lo = rLo[0];
    GridKey hi = rHi[0];
    for (std::size_t t = 0; t < n; ++t) {
        lo = GridKey{std::min(lo.x, std::min(rLo[t].x, rHi[t].x)),
                     std::min(lo.y, std::min(rLo[t].y, rHi[t].y)),
                     std::min(lo.z, std::min(rLo[t].z, rHi[t].z))};
        hi = GridKey{std::max(hi.x, std::max(rLo[t].x, rHi[t].x)),
                     std::max(hi.y, std::max(rLo[t].y, rHi[t].y)),
                     std::max(hi.z, std::max(rLo[t].z, rHi[t].z))};
    }
    return BucketTable(Cell, std::move(bucket_keys), std::move(bucket_offset),
                       std::move(bucket_ids), lo, hi);
}

namespace {

// One pyramid level from the keys of the level below: each key's parent block,
// kept the first time it is met. `rKeys` holds distinct keys, so the level is
// at most as large.
OccupancyLevel bb_level_above(const std::vector<GridKey>& rKeys) {
    OccupancyLevel level;
    std::size_t cap = 8;
    while (cap < 2 * rKeys.size())
        cap <<= 1;
    level.mMask = cap - 1;
    constexpr std::uint64_t kEmpty = ~std::uint64_t{0};
    level.mSlots.assign(cap, kEmpty);
    for (const GridKey& rKey : rKeys) {
        const GridKey up = grid_key_up(rKey, 1);
        for (std::size_t s = grid_key_home(up, level.mMask);; s = (s + 1) & level.mMask) {
            if (level.mSlots[s] == kEmpty) {
                level.mSlots[s] = level.mKeys.size();
                level.mKeys.push_back(up);
                break;
            }
            if (level.mKeys[static_cast<std::size_t>(level.mSlots[s])] == up)
                break;
        }
    }
    level.mLo = level.mKeys.front();
    level.mHi = level.mKeys.front();
    for (const GridKey& rKey : level.mKeys) {
        level.mLo = GridKey{std::min(level.mLo.x, rKey.x), std::min(level.mLo.y, rKey.y),
                            std::min(level.mLo.z, rKey.z)};
        level.mHi = GridKey{std::max(level.mHi.x, rKey.x), std::max(level.mHi.y, rKey.y),
                            std::max(level.mHi.z, rKey.z)};
    }
    return level;
}

// Whether a level spanning the key box [rLo, rHi] is small enough to be the
// top: at most four blocks per axis. In doubles, since a span of garbage keys
// (a NaN coordinate quantizes to the int64 extremes) overflows int64.
bool bb_small_enough(const GridKey& rLo, const GridKey& rHi) {
    const auto span = [](std::int64_t Lo, std::int64_t Hi) {
        return static_cast<double>(Hi) - static_cast<double>(Lo);
    };
    return span(rLo.x, rHi.x) <= 3.0 && span(rLo.y, rHi.y) <= 3.0 && span(rLo.z, rHi.z) <= 3.0;
}

}  // namespace

// Each level from the one below, until the top spans at most four blocks per
// axis. Serial: a level has at most as many keys as the table has buckets, and
// building one is a probe per key, small beside the table's own build.
OccupancyPyramid build_occupancy_pyramid(const BucketTable& rGrid) {
    constexpr std::size_t kMaxLevels = 20;
    OccupancyPyramid pyramid;
    if (rGrid.Empty() || bb_small_enough(rGrid.OccupiedLo(), rGrid.OccupiedHi()))
        return pyramid;
    while (pyramid.mLevels.size() < kMaxLevels) {
        const std::vector<GridKey>& below =
            pyramid.mLevels.empty() ? rGrid.Keys() : pyramid.mLevels.back().mKeys;
        pyramid.mLevels.push_back(bb_level_above(below));
        const OccupancyLevel& top = pyramid.mLevels.back();
        if (bb_small_enough(top.mLo, top.mHi))
            break;
    }
    return pyramid;
}

}  // namespace detail
}  // namespace meshioplusplus
