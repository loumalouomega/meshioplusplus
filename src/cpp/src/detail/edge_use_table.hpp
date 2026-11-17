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

#pragma once

/**
 * @file detail/edge_use_table.hpp
 * @brief `EdgeUseTable`: the uses of each undirected edge of a face list,
 * grouped by edge in `(lo, hi)` order (roadmap 6.2).
 *
 * Fed the edges of a surface one `Add(a, b, face)` at a time and `Build()` once,
 * it answers what a `std::map` keyed by the sorted vertex pair did: how many
 * faces use an edge (one is a boundary edge, two a manifold one), which faces
 * and in which direction, and every edge in key order, so a message or a
 * tie-break that followed the map's order keeps it. The uses of an edge keep
 * their insertion order, so the first use is the first face that met the edge.
 *
 * A sorted flat table: no node per edge, one sort. A **core-private** header
 * (not installed).
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

namespace meshioplusplus::detail {

template <class Index>
class EdgeUseTable {
public:
    static constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

    /// One use of an edge: the face and whether it walks the edge `lo -> hi`.
    struct Use {
        Index mFace;
        bool mForward;
    };

    void Reserve(std::size_t NumUses) { mRecords.reserve(NumUses); }

    /// Record face `Face` walking the edge from `A` to `B`.
    void Add(Index A, Index B, Index Face) {
        const bool forward = A < B;
        mRecords.push_back({forward ? A : B, forward ? B : A, Face, forward});
    }

    /// Group the recorded uses by edge. Call once, after the last `Add`.
    void Build() {
        std::stable_sort(mRecords.begin(), mRecords.end(), [](const Record& rA, const Record& rB) {
            return rA.mLo != rB.mLo ? rA.mLo < rB.mLo : rA.mHi < rB.mHi;
        });
        mStart.assign(1, 0);
        for (std::size_t i = 1; i < mRecords.size(); ++i)
            if (mRecords[i].mLo != mRecords[i - 1].mLo || mRecords[i].mHi != mRecords[i - 1].mHi)
                mStart.push_back(i);
        if (!mRecords.empty())
            mStart.push_back(mRecords.size());
        else
            mStart.clear();
    }

    /// How many distinct edges there are.
    std::size_t NumEdges() const { return mStart.empty() ? 0 : mStart.size() - 1; }
    Index Lo(std::size_t Edge) const { return mRecords[mStart[Edge]].mLo; }
    Index Hi(std::size_t Edge) const { return mRecords[mStart[Edge]].mHi; }

    /// How many faces use the edge.
    std::size_t Count(std::size_t Edge) const { return mStart[Edge + 1] - mStart[Edge]; }

    /// The `K`-th use of the edge, in insertion order.
    Use At(std::size_t Edge, std::size_t K) const {
        const Record& r = mRecords[mStart[Edge] + K];
        return {r.mFace, r.mForward};
    }

    /// The index of the edge `{A, B}`, or `kNone`.
    std::size_t Find(Index A, Index B) const {
        const Index lo = A < B ? A : B, hi = A < B ? B : A;
        std::size_t first = 0, last = NumEdges();
        while (first < last) {
            const std::size_t mid = first + (last - first) / 2;
            const Record& r = mRecords[mStart[mid]];
            if (r.mLo < lo || (r.mLo == lo && r.mHi < hi))
                first = mid + 1;
            else
                last = mid;
        }
        if (first < NumEdges() && mRecords[mStart[first]].mLo == lo &&
            mRecords[mStart[first]].mHi == hi)
            return first;
        return kNone;
    }

    /// How many faces use the edge `{A, B}`; zero when none does.
    std::size_t CountOf(Index A, Index B) const {
        const std::size_t e = Find(A, B);
        return e == kNone ? 0 : Count(e);
    }

    /// Call `Join(first_face, face)` for every later use of every edge: the
    /// faces an edge shares, each joined to the first that met it.
    template <class Fn>
    void ForEachSharedFace(Fn&& Join) const {
        for (std::size_t e = 0; e < NumEdges(); ++e)
            for (std::size_t k = 1; k < Count(e); ++k)
                Join(At(e, 0).mFace, At(e, k).mFace);
    }

private:
    struct Record {
        Index mLo, mHi, mFace;
        bool mForward;
    };
    std::vector<Record> mRecords;
    std::vector<std::size_t> mStart;
};

}  // namespace meshioplusplus::detail
