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
//
// Unit tests for detail/splitmix.hpp -- the splitmix64 finalizer and generator
// that replaced three hand-written copies (roadmap 6.4).

// System includes
#include <cstdint>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/splitmix.hpp"

using meshioplusplus::detail::SplitMix64;
using meshioplusplus::detail::splitmix64_mix;

// The published reference outputs of splitmix64 for seed 0.
TEST(SplitMix, ReferenceValuesForSeedZero) {
    SplitMix64 rng{0};
    EXPECT_EQ(rng.Next(), 0xE220A8397B1DCDAFull);
    EXPECT_EQ(rng.Next(), 0x6E789E6AA1B965F4ull);
    EXPECT_EQ(rng.Next(), 0x06C45D188009454Full);
}

TEST(SplitMix, MixIsTheGeneratorStep) {
    SplitMix64 rng{12345};
    const std::uint64_t want = splitmix64_mix(12345ull + 0x9E3779B97F4A7C15ull);
    EXPECT_EQ(rng.Next(), want);
}

TEST(SplitMix, NextUnitIsTheTop53BitsInHalfOpenUnit) {
    SplitMix64 a{7}, b{7};
    for (int i = 0; i < 1000; ++i) {
        const double u = a.NextUnit();
        EXPECT_GE(u, 0.0);
        EXPECT_LT(u, 1.0);
        EXPECT_EQ(u, static_cast<double>(b.Next() >> 11) * (1.0 / 9007199254740992.0));
    }
}
