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
 * @file detail/exact_weld.hpp
 * @brief First-seen numbering of identical coordinate rows: the vertex weld of
 * the STL, glTF and WKT readers.
 *
 * A **core-private** header (the `detail/slot_runs.hpp` precedent). The readers
 * used to key a hash map on a heap-allocated string or vector per vertex; a
 * row is now a fixed-size key of its coordinate bits, grouped by
 * `group_slots` and numbered by `number_first_seen`, which reproduces the ids
 * a serial sweep handing the next id to each unseen row gives (roadmap 5.7.3).
 * Two equalities are offered:
 *  - `ExactWeldMode::Bits`: rows are equal when their bytes are (STL: two
 *    zeros of opposite sign stay two points, a NaN equals the same NaN);
 *  - `ExactWeldMode::Values`: rows are equal when every coordinate compares
 *    `==` (glTF, WKT): `-0.0` and `+0.0` are one point, and a row holding a
 *    NaN equals no row, itself included.
 */

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

// Project includes
#include "slot_runs.hpp"

namespace meshioplusplus {
namespace detail {

enum class ExactWeldMode { Bits, Values };

/// What `exact_weld` returns: each row's id, and each id's first row.
struct ExactWeld {
    std::vector<std::int64_t> mIdOfRow;
    std::vector<std::size_t> mFirstRow;

    std::size_t NumIds() const { return mFirstRow.size(); }
};

/// One row's key: its width, up to four coordinates' bits (zero-padded), and
/// a tiebreak that is the row's own index for a row equal to no other.
struct ExactWeldKey {
    std::array<std::uint64_t, 4> mWords{};
    std::uint64_t mWidth = 0;
    std::uint64_t mAlone = 0;

    bool operator<(const ExactWeldKey& rOther) const {
        if (mWidth != rOther.mWidth)
            return mWidth < rOther.mWidth;
        if (mWords != rOther.mWords)
            return mWords < rOther.mWords;
        return mAlone < rOther.mAlone;
    }
    std::uint64_t Hash() const {
        std::uint64_t h = 1469598103934665603ull ^ mWidth;
        for (const std::uint64_t w : mWords)
            h = (h ^ w) * 1099511628211ull;
        h = (h ^ mAlone) * 1099511628211ull;
        return h ^ (h >> 29);
    }
};

/// The key of a row of `Width` (at most 4) coordinates of type `T` (`float`
/// or `double`) under `Mode`; `Row` is the row's index.
template <class T>
ExactWeldKey exact_weld_key(const T* pRow, std::size_t Width, ExactWeldMode Mode, std::size_t Row) {
    ExactWeldKey key;
    key.mWidth = Width;
    for (std::size_t j = 0; j < Width; ++j) {
        T v = pRow[j];
        if (Mode == ExactWeldMode::Values) {
            if (std::isnan(v))
                key.mAlone = static_cast<std::uint64_t>(Row) + 1;
            if (v == T(0))
                v = T(0);  // -0 -> +0: the two compare equal
        }
        std::memcpy(&key.mWords[j], &v, sizeof(T));
    }
    return key;
}

/// Weld `NumRows` keys (one per row, in row order).
inline ExactWeld exact_weld_keys(const std::vector<ExactWeldKey>& rKeys) {
    ExactWeld out;
    const std::size_t n = rKeys.size();
    if (n == 0)
        return out;
    const std::size_t buckets = n / 4 + 1;
    const SlotRuns runs = group_slots(
        rKeys, buckets, [&](const ExactWeldKey& rKey) { return rKey.Hash() % buckets; });
    FirstSeen seen = number_first_seen(runs, n);
    out.mIdOfRow = std::move(seen.mIdOfSlot);
    out.mFirstRow.resize(seen.NumIds());
    parallel_for(seen.NumIds(), [&](std::size_t id) {
        out.mFirstRow[id] = static_cast<std::size_t>(runs.Head(seen.mRunOfId[id]));
    });
    return out;
}

/**
 * @brief Weld the `NumRows` rows of `Width` (1 to 4) coordinates at `pRows`
 * (row-major, `T` = `float` or `double`): ids in order of first appearance.
 */
template <class T>
ExactWeld exact_weld(const T* pRows, std::size_t NumRows, std::size_t Width, ExactWeldMode Mode) {
    std::vector<ExactWeldKey> keys(NumRows);
    parallel_for(NumRows, [&](std::size_t i) {
        keys[i] = exact_weld_key(pRows + i * Width, Width, Mode, i);
    });
    return exact_weld_keys(keys);
}

}  // namespace detail
}  // namespace meshioplusplus
