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
 * @file detail/slot_runs.hpp
 * @brief The sort-based facet and edge table: equal keys of a slot buffer,
 * grouped, counted and numbered in parallel with serial-sweep results.
 *
 * A **core-private** header (the `formats/gid_common.hpp` precedent): it lives
 * beside the `.cpp` files, no installed header names it, and it adds nothing
 * to the API or the ABI.
 *
 * Six operations fill one key per *slot* -- (block, cell, local facet or
 * edge), in an order that is a pure function of the mesh -- in parallel, and
 * then used to deduplicate them through a single-threaded hash map: a sweep
 * over the slots that counts each key (the boundary is the keys seen once) or
 * numbers each key the first time it is seen. Both answers are a function of
 * the runs of equal keys once (key, slot) pairs are sorted, and sorting a
 * total order gives the same sequence on every backend and thread count:
 *
 *  - `group_slots` groups the slots by a bucket each key names -- one of its
 *    node ids -- with a stable counting sort, then sorts each small bucket by
 *    (key, slot) in parallel: linear, so a SEQ build is not slower than the
 *    hash map; each run's slots ascend, so a run's first slot is where the
 *    serial sweep first met the key (`group_slots_sorted` is the comparison
 *    sort, for a table searched by key later);
 *  - `number_first_seen` flags each run's first slot, and an exclusive
 *    `parallel_exclusive_scan` of the flags in slot order hands out ids in
 *    exactly the order the serial sweep did;
 *  - a run's size is its key's count.
 *
 * Roadmap §4, "One shared, deterministic facet and edge table".
 */

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

/// Runs of equal keys, each holding its slots in ascending order.
struct SlotRuns {
    std::vector<std::uint64_t> mSlots;     ///< every slot, grouped by run
    std::vector<std::uint64_t> mStart{0};  ///< `NumRuns() + 1` offsets into `mSlots`

    std::size_t NumRuns() const { return mStart.size() - 1; }
    /// The run's first slot: where a serial sweep first met its key.
    std::uint64_t Head(std::size_t Run) const { return mSlots[mStart[Run]]; }
    /// How many slots hold the run's key.
    std::size_t Size(std::size_t Run) const {
        return static_cast<std::size_t>(mStart[Run + 1] - mStart[Run]);
    }
    const std::uint64_t* Begin(std::size_t Run) const { return mSlots.data() + mStart[Run]; }
    const std::uint64_t* End(std::size_t Run) const { return mSlots.data() + mStart[Run + 1]; }
};

namespace slot_runs_impl {

/// Cut an ordering of the slots in which equal keys are adjacent into runs.
template <class Key, class Less, class SameGroup>
SlotRuns cut_runs(const std::vector<Key>& rKeys, std::vector<std::uint64_t> order, Less less,
                  SameGroup same_group) {
    const std::size_t n = order.size();
    std::vector<std::uint8_t> starts(n);
    parallel_for(n, [&](std::size_t i) {
        starts[i] = i == 0 || !same_group(order[i - 1], order[i]) ||
                    less(rKeys[order[i - 1]], rKeys[order[i]]);
    });
    std::vector<std::uint64_t> run_of(n);
    const std::uint64_t runs =
        parallel_exclusive_scan(starts.data(), n, run_of.data(), std::uint64_t{0});
    SlotRuns out;
    out.mStart.assign(runs + 1, n);
    parallel_for(n, [&](std::size_t i) {
        if (starts[i])
            out.mStart[run_of[i]] = i;
    });
    out.mSlots = std::move(order);
    return out;
}

}  // namespace slot_runs_impl

/**
 * @brief Group the slots of `rKeys` (one key per slot) into runs of equal keys,
 * **in ascending key order** -- for a table that is binary-searched later.
 * A comparison sort of (key, slot) pairs; prefer the bucketed overload below
 * wherever only the grouping matters.
 * @tparam Key With `Less` a strict weak order whose equivalence is equality.
 */
template <class Key, class Less = std::less<Key>>
SlotRuns group_slots_sorted(const std::vector<Key>& rKeys, Less less = {}) {
    const std::size_t n = rKeys.size();
    std::vector<std::uint64_t> order(n);
    parallel_for_bw(n, [&](std::size_t i) { order[i] = i; });
    // (key, slot) is a total order: every backend returns this one sequence.
    parallel_sort(order.begin(), order.end(), [&](std::uint64_t a, std::uint64_t b) {
        if (less(rKeys[a], rKeys[b]))
            return true;
        if (less(rKeys[b], rKeys[a]))
            return false;
        return a < b;
    });
    return slot_runs_impl::cut_runs(rKeys, std::move(order), less,
                                    [](std::uint64_t, std::uint64_t) { return true; });
}

/**
 * @brief Group the slots of `rKeys` into runs of equal keys, by bucket.
 *
 * `bucket_of(key)` maps a key to `[0, NumBuckets)` -- in practice one of its
 * node ids, so equal keys share a bucket and buckets hold a handful of slots.
 * A stable counting sort puts the slots in (bucket, slot) order in O(n); each
 * bucket is then sorted by (key, slot), in parallel. Runs come in ascending
 * (bucket, key) order -- the global key order too when the bucket is the
 * key's leading component. Linear where a comparison sort is O(n log n), so
 * a SEQ build is not slower than the hash map this replaces.
 */
template <class Key, class BucketOf, class Less = std::less<Key>>
SlotRuns group_slots(const std::vector<Key>& rKeys, std::size_t NumBuckets, BucketOf bucket_of,
                     Less less = {}) {
    const std::size_t n = rKeys.size();
    std::vector<std::uint32_t> bucket32;
    std::vector<std::uint64_t> bucket64;
    const bool narrow = NumBuckets <= 0xFFFFFFFFull;
    if (narrow)
        bucket32.resize(n);
    else
        bucket64.resize(n);
    parallel_for(n, [&](std::size_t i) {
        const std::size_t b = bucket_of(rKeys[i]);
        if (narrow)
            bucket32[i] = static_cast<std::uint32_t>(b);
        else
            bucket64[i] = b;
    });
    const auto bucket = [&](std::uint64_t slot) -> std::uint64_t {
        return narrow ? bucket32[slot] : bucket64[slot];
    };
    // Stable counting sort by bucket: exact integer counts, so the order is
    // (bucket, slot) on every backend.
    std::vector<std::uint64_t> start(NumBuckets + 1, 0);
    for (std::size_t i = 0; i < n; ++i)
        ++start[bucket(i) + 1];
    for (std::size_t b = 0; b < NumBuckets; ++b)
        start[b + 1] += start[b];
    std::vector<std::uint64_t> order(n);
    {
        std::vector<std::uint64_t> cursor(start.begin(), start.end() - 1);
        for (std::size_t i = 0; i < n; ++i)
            order[cursor[bucket(i)]++] = i;
    }
    // Within a bucket the slots ascend already; sort them by (key, slot).
    parallel_for(
        NumBuckets,
        [&](std::size_t b) {
            std::uint64_t* first = order.data() + start[b];
            std::uint64_t* last = order.data() + start[b + 1];
            if (last - first < 2)
                return;
            std::sort(first, last, [&](std::uint64_t x, std::uint64_t y) {
                if (less(rKeys[x], rKeys[y]))
                    return true;
                if (less(rKeys[y], rKeys[x]))
                    return false;
                return x < y;
            });
        },
        256);
    return slot_runs_impl::cut_runs(
        rKeys, std::move(order), less,
        [&](std::uint64_t x, std::uint64_t y) { return bucket(x) == bucket(y); });
}

/// First-seen numbering of a slot buffer's keys (`number_first_seen`).
struct FirstSeen {
    std::vector<std::int64_t> mIdOfSlot;  ///< per slot: its key's id
    std::vector<std::uint64_t> mRunOfId;  ///< per id: its run in the `SlotRuns`

    std::size_t NumIds() const { return mRunOfId.size(); }
};

/**
 * @brief Number the keys in order of their first slot -- the ids a serial
 * sweep over slots `0 .. NumSlots-1` gives when it hands the next id to each
 * key it has not seen. `NumSlots` is the length of the key buffer `rRuns` was
 * built from.
 */
inline FirstSeen number_first_seen(const SlotRuns& rRuns, std::size_t NumSlots) {
    std::vector<std::uint8_t> is_head(NumSlots, 0);
    const std::size_t runs = rRuns.NumRuns();
    parallel_for(runs, [&](std::size_t r) { is_head[rRuns.Head(r)] = 1; });
    std::vector<std::int64_t> id_at(NumSlots);
    parallel_exclusive_scan(is_head.data(), NumSlots, id_at.data(), std::int64_t{0});
    FirstSeen out;
    out.mIdOfSlot.resize(NumSlots);
    out.mRunOfId.resize(runs);
    parallel_for(runs, [&](std::size_t r) {
        const std::int64_t id = id_at[rRuns.Head(r)];
        out.mRunOfId[static_cast<std::size_t>(id)] = r;
        for (const std::uint64_t* p = rRuns.Begin(r); p != rRuns.End(r); ++p)
            out.mIdOfSlot[*p] = id;
    });
    return out;
}

/// A `FacetKey` of at most `FacetKey::kInline` ids, flattened for sorting:
/// its size, then its sorted ids, zero-padded. Two flattened keys are equal
/// exactly when the `FacetKey`s are.
using FlatFacetKey = std::array<std::int64_t, FacetKey::kInline + 1>;

/// The flattened key of the at most `FacetKey::kInline` ids at @p pIds (any
/// order): what `group_facet_slots` flattens a `FacetKey` to.
inline FlatFacetKey flat_facet_key(const std::int64_t* pIds, std::size_t N) {
    FlatFacetKey f{};
    f[0] = static_cast<std::int64_t>(N);
    for (std::size_t j = 0; j < N; ++j)
        f[j + 1] = pIds[j];
    std::sort(f.begin() + 1, f.begin() + 1 + static_cast<std::ptrdiff_t>(N));
    return f;
}

/// Total order on `FacetKey`s: by size, then lexicographically by sorted id.
struct FacetKeyLess {
    bool operator()(const FacetKey& rA, const FacetKey& rB) const {
        if (rA.Size() != rB.Size())
            return rA.Size() < rB.Size();
        return std::lexicographical_compare(rA.Data(), rA.Data() + rA.Size(), rB.Data(),
                                            rB.Data() + rB.Size());
    }
};

/// The bucket of a facet key for `group_slots`: its smallest id, when that is a
/// node of the mesh; keys naming no valid node share bucket `NumPoints`.
inline std::size_t facet_bucket(const std::int64_t* pSortedIds, std::size_t N,
                                std::size_t NumPoints) {
    if (N == 0 || pSortedIds[0] < 0 || static_cast<std::size_t>(pSortedIds[0]) >= NumPoints)
        return NumPoints;
    return static_cast<std::size_t>(pSortedIds[0]);
}

/**
 * @brief `group_slots` for facet keys, bucketed by their smallest node id:
 * flattened to fixed-width arrays when no key has more than
 * `FacetKey::kInline` ids (every non-polyhedral facet), compared through
 * `FacetKeyLess` otherwise. The runs are the same either way.
 * @tparam Rec A record type; `rKeyOf(rec)` returns its `const FacetKey&`.
 */
template <class Rec, class KeyOf>
SlotRuns group_facet_slots(const std::vector<Rec>& rRecs, KeyOf&& rKeyOf, std::size_t NumPoints) {
    const std::size_t n = rRecs.size();
    const std::size_t widest = parallel_reduce(
        n, 4096, std::size_t{0},
        [&](std::size_t b, std::size_t e) {
            std::size_t w = 0;
            for (std::size_t i = b; i < e; ++i)
                w = std::max(w, rKeyOf(rRecs[i]).Size());
            return w;
        },
        [](std::size_t a, std::size_t p) { return std::max(a, p); });
    if (widest <= FacetKey::kInline) {
        std::vector<FlatFacetKey> flat(n);
        parallel_for_bw(n, [&](std::size_t i) {
            const FacetKey& k = rKeyOf(rRecs[i]);
            FlatFacetKey& f = flat[i];
            f.fill(0);
            f[0] = static_cast<std::int64_t>(k.Size());
            for (std::size_t j = 0; j < k.Size(); ++j)
                f[j + 1] = k.Data()[j];
        });
        return group_slots(flat, NumPoints + 1, [&](const FlatFacetKey& rK) {
            return facet_bucket(rK.data() + 1, static_cast<std::size_t>(rK[0]), NumPoints);
        });
    }
    std::vector<FacetKey> keys(n);
    parallel_for_bw(n, [&](std::size_t i) { keys[i] = rKeyOf(rRecs[i]); });
    return group_slots(
        keys, NumPoints + 1,
        [&](const FacetKey& rK) { return facet_bucket(rK.Data(), rK.Size(), NumPoints); },
        FacetKeyLess{});
}

/// Per slot, how many slots hold its key (the count-only rule).
inline std::vector<std::uint32_t> slot_multiplicity(const SlotRuns& rRuns, std::size_t NumSlots) {
    std::vector<std::uint32_t> count(NumSlots, 0);
    parallel_for(rRuns.NumRuns(), [&](std::size_t r) {
        const auto c = static_cast<std::uint32_t>(rRuns.Size(r));
        for (const std::uint64_t* p = rRuns.Begin(r); p != rRuns.End(r); ++p)
            count[*p] = c;
    });
    return count;
}

}  // namespace detail
}  // namespace meshioplusplus
