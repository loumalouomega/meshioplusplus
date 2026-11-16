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
 * @file detail/splitmix.hpp
 * @brief splitmix64 (Steele, Lea and Flood 2014): the finalizer and the
 * counter-based generator, written once (roadmap 6.4).
 *
 * The same sequence on every platform, unlike a standard-library distribution,
 * so seeded operations (`box_points`, `refine_patches`' RANSAC) stay
 * reproducible. A **core-private** header (not installed).
 */

// System includes
#include <cstdint>

// Project includes
#include "always_inline.hpp"

namespace meshioplusplus::detail {

/// The golden-ratio increment of splitmix64.
inline constexpr std::uint64_t kSplitMix64Gamma = 0x9E3779B97F4A7C15ULL;

/// The splitmix64 finalizer: a bijective mix of one 64-bit word.
MIO_ALWAYS_INLINE std::uint64_t splitmix64_mix(std::uint64_t Z) {
    Z = (Z ^ (Z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    Z = (Z ^ (Z >> 27)) * 0x94D049BB133111EBULL;
    return Z ^ (Z >> 31);
}

/// splitmix64 over a caller-owned state word.
struct SplitMix64 {
    std::uint64_t mState;

    /// The next 64-bit value.
    MIO_ALWAYS_INLINE std::uint64_t Next() {
        mState += kSplitMix64Gamma;
        return splitmix64_mix(mState);
    }

    /// A uniform draw in [0, 1): the top 53 bits, exact in a double.
    MIO_ALWAYS_INLINE double NextUnit() {
        return static_cast<double>(Next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

}  // namespace meshioplusplus::detail
