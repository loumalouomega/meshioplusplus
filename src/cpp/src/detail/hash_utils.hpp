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
 * @file detail/hash_utils.hpp
 * @brief `HashUtils`: the boost-style hash combiner the key hashes of the core
 * wrote out by hand, and the key hashes built on it (roadmap 6.3).
 *
 * `Combine` is exactly `h ^= v + 0x9e3779b97f4a7c15 + (h << 6) + (h >> 2)`, so
 * a key hash that switches to it keeps every hash value, and with it the
 * iteration order of the unordered container it serves. The integer keys go
 * through `std::hash<std::int64_t>` as before. A **core-private** header (not
 * installed); the two installed key hashes (`SurfaceEdgeKeyHash`,
 * `GridKeyHash`) keep their own inline bodies until an ABI bump is due.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

// Project includes
#include "always_inline.hpp"

namespace meshioplusplus::detail {

/// Static hash-combining helpers.
struct HashUtils {
    /// `Seed ^ (Value + 0x9e3779b97f4a7c15 + (Seed << 6) + (Seed >> 2))`.
    template <class T>
    static MIO_ALWAYS_INLINE T Combine(T Seed, T Value) {
        return Seed ^ (Value + static_cast<T>(0x9e3779b97f4a7c15ULL) + (Seed << 6) + (Seed >> 2));
    }

    /// Combines every element of a range, each through `std::hash<std::int64_t>`.
    template <class Range>
    static std::size_t Range64(const Range& rKeys, std::size_t Seed) {
        std::size_t h = Seed;
        for (const std::int64_t v : rKeys)
            h = Combine<std::size_t>(h, std::hash<std::int64_t>{}(v));
        return h;
    }
};

/// Hash of a fixed-size integer key, seed 0.
template <std::size_t N>
struct ArrayHash {
    std::size_t operator()(const std::array<std::int64_t, N>& rKey) const {
        return HashUtils::Range64(rKey, 0);
    }
};

/// Hash of a vector key, seed its size.
struct SizedVectorHash {
    std::size_t operator()(const std::vector<std::int64_t>& rKey) const {
        return HashUtils::Range64(rKey, rKey.size());
    }
};

}  // namespace meshioplusplus::detail
