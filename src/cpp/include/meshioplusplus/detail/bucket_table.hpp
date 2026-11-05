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
#pragma once

/**
 * @file detail/bucket_table.hpp
 * @brief A read-only bucket grid: the frozen counterpart of `SpatialGrid`.
 *
 * `SpatialGrid` (`detail/spatial_hash.hpp`) is a growable map from a quantized
 * cell key to the ids inserted into it, one `std::vector` per cell inside an
 * `unordered_map`. That is what `merge` needs, since it inserts as it welds.
 * `DistanceQuery` never inserts after the build: its (cell, triangle) pairs are
 * grouped once, in parallel, and then only searched. A `BucketTable` holds
 * exactly that: the buckets in compressed-sparse-row form -- the keys in the
 * order the build met them, the offsets of each bucket's ids, and the ids of
 * every bucket in one flat array -- with a flat open-addressing index from a
 * key to its bucket. There is no node and no vector per cell, and a lookup is
 * a probe over contiguous memory instead of a pointer chase.
 *
 * What is observable is `SpatialGrid`'s contract, unchanged: every bucket lists
 * its ids ascending (the caller hands them over that way), `ForEachInShell` and
 * `ForEachInBox` visit cells in the same ascending dz -> dy -> dx / z -> y -> x
 * order, and the occupied-key box is the one the caller gives. Only a *lookup*
 * touches the index, and a lookup returns the same bucket whatever the probe
 * order was, so the hash cannot affect any result.
 *
 * `detail/` header: exempt from the operations layer's anon-namespace prefix
 * rule.
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/spatial_hash.hpp"

namespace meshioplusplus {
namespace detail {

/// The ids of one bucket: a borrowed, ascending range inside a `BucketTable`.
class BucketView {
public:
    BucketView(const std::int64_t* pFirst, std::size_t Count) : mpFirst(pFirst), mCount(Count) {}

    const std::int64_t* begin() const { return mpFirst; }
    const std::int64_t* end() const { return mpFirst + mCount; }
    std::size_t size() const { return mCount; }
    bool empty() const { return mCount == 0; }
    std::int64_t operator[](std::size_t i) const { return mpFirst[i]; }

private:
    const std::int64_t* mpFirst;
    std::size_t mCount;
};

/// The read-only bucket grid: cell key -> ascending ids, built once.
class BucketTable {
public:
    explicit BucketTable(double CellSize = 1.0) : mCell(CellSize) {}

    /**
     * @brief Take ownership of grouped buckets.
     * @param CellSize the quantization size of `KeyOf`.
     * @param Keys one distinct key per bucket, in the caller's bucket order.
     * @param Offsets `Keys.size() + 1` entries: bucket `b` holds
     *        `Ids[Offsets[b] .. Offsets[b + 1])`, ascending.
     * @param Ids every bucket's ids, back to back.
     * @param rLo the occupied key box's low corner (valid when `Keys` is not empty).
     * @param rHi the occupied key box's high corner.
     */
    BucketTable(double CellSize, std::vector<GridKey> Keys, std::vector<std::uint64_t> Offsets,
                std::vector<std::int64_t> Ids, const GridKey& rLo, const GridKey& rHi)
        : mCell(CellSize),
          mKeys(std::move(Keys)),
          mOffsets(std::move(Offsets)),
          mIds(std::move(Ids)),
          mLo(rLo),
          mHi(rHi) {
        BuildIndex();
    }

    double CellSize() const { return mCell; }

    /// The cell key of a (z-padded, 3-component) coordinate.
    GridKey KeyOf(const double* pCoords) const {
        return GridKey{grid_quantize(pCoords[0], mCell), grid_quantize(pCoords[1], mCell),
                       grid_quantize(pCoords[2], mCell)};
    }

    bool Empty() const { return mKeys.empty(); }

    /// The number of non-empty cells.
    std::size_t NumBuckets() const { return mKeys.size(); }

    /// The occupied-key bounding box (valid only when not `Empty()`).
    const GridKey& OccupiedLo() const { return mLo; }
    const GridKey& OccupiedHi() const { return mHi; }

    /// The ids in the cell at `rKey`, or false (leaving @p rOut alone) if it is empty.
    bool Find(const GridKey& rKey, BucketView& rOut) const {
        const std::size_t b = FindBucket(rKey);
        if (b == kNone)
            return false;
        rOut = Bucket(b);
        return true;
    }

    /// Visits the cells at Chebyshev radius exactly `r` around `rCenter`
    /// (r = 0 is the centre cell alone), calling `fn(BucketView)` for each
    /// non-empty cell, in ascending dz -> dy -> dx order -- `SpatialGrid`'s
    /// shell, clamped to the occupied box the same way.
    template <class F>
    void ForEachInShell(const GridKey& rCenter, std::int64_t r, F&& fn) const {
        if (mKeys.empty())
            return;
        const std::int64_t zl = std::max(-r, mLo.z - rCenter.z);
        const std::int64_t zh = std::min(r, mHi.z - rCenter.z);
        const std::int64_t yl = std::max(-r, mLo.y - rCenter.y);
        const std::int64_t yh = std::min(r, mHi.y - rCenter.y);
        const std::int64_t xl = std::max(-r, mLo.x - rCenter.x);
        const std::int64_t xh = std::min(r, mHi.x - rCenter.x);
        if (zl > zh || yl > yh || xl > xh)
            return;
        auto visit = [&](std::int64_t dx, std::int64_t dy, std::int64_t dz) {
            const std::size_t b =
                FindBucket(GridKey{rCenter.x + dx, rCenter.y + dy, rCenter.z + dz});
            if (b != kNone)
                fn(Bucket(b));
        };
        if (r == 0) {
            visit(0, 0, 0);
            return;
        }
        for (std::int64_t dz = zl; dz <= zh; ++dz) {
            if (dz == -r || dz == r) {
                for (std::int64_t dy = yl; dy <= yh; ++dy)
                    for (std::int64_t dx = xl; dx <= xh; ++dx)
                        visit(dx, dy, dz);
            } else {
                for (std::int64_t dy = yl; dy <= yh; ++dy) {
                    if (dy == -r || dy == r) {
                        for (std::int64_t dx = xl; dx <= xh; ++dx)
                            visit(dx, dy, dz);
                    } else {
                        if (xl == -r)
                            visit(-r, dy, dz);
                        if (xh == r)
                            visit(r, dy, dz);
                    }
                }
            }
        }
    }

    /// Visits every non-empty cell in the inclusive key box `[rLo, rHi]`,
    /// calling `fn(BucketView)` for each, in ascending z -> y -> x order.
    template <class F>
    void ForEachInBox(const GridKey& rLo, const GridKey& rHi, F&& fn) const {
        for (std::int64_t z = rLo.z; z <= rHi.z; ++z)
            for (std::int64_t y = rLo.y; y <= rHi.y; ++y)
                for (std::int64_t x = rLo.x; x <= rHi.x; ++x) {
                    const std::size_t b = FindBucket(GridKey{x, y, z});
                    if (b != kNone)
                        fn(Bucket(b));
                }
    }

private:
    static constexpr std::size_t kNone = static_cast<std::size_t>(-1);
    static constexpr std::uint64_t kEmptySlot = ~std::uint64_t{0};

    BucketView Bucket(std::size_t b) const {
        return BucketView(mIds.data() + mOffsets[b],
                          static_cast<std::size_t>(mOffsets[b + 1] - mOffsets[b]));
    }

    /// The slot a key starts probing at: the key's hash, mixed once more so
    /// the low bits that the mask keeps depend on every axis.
    std::size_t Home(const GridKey& rKey) const {
        std::uint64_t h = static_cast<std::uint64_t>(GridKeyHash{}(rKey));
        h ^= h >> 29;
        h *= 0xbf58476d1ce4e5b9ULL;
        h ^= h >> 32;
        return static_cast<std::size_t>(h) & mMask;
    }

    std::size_t FindBucket(const GridKey& rKey) const {
        if (mSlots.empty())
            return kNone;
        for (std::size_t s = Home(rKey);; s = (s + 1) & mMask) {
            const std::uint64_t b = mSlots[s];
            if (b == kEmptySlot)
                return kNone;
            if (mKeys[static_cast<std::size_t>(b)] == rKey)
                return static_cast<std::size_t>(b);
        }
    }

    /// A power-of-two table at most half full, so a probe is short and a miss
    /// (most of a shell's cells are empty) ends at the first free slot.
    void BuildIndex() {
        if (mKeys.empty())
            return;
        std::size_t capacity = 8;
        while (capacity < 2 * mKeys.size())
            capacity *= 2;
        mMask = capacity - 1;
        mSlots.assign(capacity, kEmptySlot);
        for (std::size_t b = 0; b < mKeys.size(); ++b) {
            std::size_t s = Home(mKeys[b]);
            while (mSlots[s] != kEmptySlot)
                s = (s + 1) & mMask;
            mSlots[s] = b;
        }
    }

    double mCell;
    std::vector<GridKey> mKeys;
    std::vector<std::uint64_t> mOffsets;
    std::vector<std::int64_t> mIds;
    std::vector<std::uint64_t> mSlots;
    std::size_t mMask = 0;
    GridKey mLo{0, 0, 0};
    GridKey mHi{0, 0, 0};
};

}  // namespace detail
}  // namespace meshioplusplus
