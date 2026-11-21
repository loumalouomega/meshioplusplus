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

#pragma once

/**
 * @file detail/union_find.hpp
 * @brief `UnionFind<Index>`: disjoint sets with path halving, whose root is
 * always the smallest index of its set (roadmap 6.2).
 *
 * The smaller root wins a join, so the partition, its representatives and any
 * group numbering derived from them do not depend on the order the joins are
 * applied in. Serial by design (the union sweeps are a non-goal to
 * parallelise). A **core-private** header (not installed).
 */

// System includes
#include <cstddef>
#include <numeric>
#include <vector>

// Project includes
#include "always_inline.hpp"

namespace meshioplusplus::detail {

template <class Index>
class UnionFind {
public:
    explicit UnionFind(std::size_t N) : mParent(N) {
        std::iota(mParent.begin(), mParent.end(), Index{0});
    }

    /// The root of `X`'s set, halving the path on the way.
    MIO_ALWAYS_INLINE Index Find(Index X) { return FindIn(mParent, X); }

    /// Join the sets of `A` and `B`, the smaller root winning; true when they were apart.
    MIO_ALWAYS_INLINE bool Unite(Index A, Index B) {
        A = Find(A);
        B = Find(B);
        if (A == B)
            return false;
        if (A < B)
            mParent[static_cast<std::size_t>(B)] = A;
        else
            mParent[static_cast<std::size_t>(A)] = B;
        return true;
    }

    /// The parent link of `X`: never above `X`, since a root is its set's
    /// smallest index and path halving only moves links down.
    Index Parent(std::size_t X) const { return mParent[X]; }

    /// Path-halving find over a caller-owned parent vector, for a union written
    /// with its own rule beside the same walk.
    static MIO_ALWAYS_INLINE Index FindIn(std::vector<Index>& rParent, Index X) {
        while (rParent[static_cast<std::size_t>(X)] != X) {
            Index& p = rParent[static_cast<std::size_t>(X)];
            p = rParent[static_cast<std::size_t>(p)];
            X = p;
        }
        return X;
    }

private:
    std::vector<Index> mParent;
};

}  // namespace meshioplusplus::detail
