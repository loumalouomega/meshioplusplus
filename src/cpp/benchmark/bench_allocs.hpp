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
//
#pragma once

/**
 * @file bench_allocs.hpp
 * @brief Opt-in global `operator new`/`delete` counter for
 * `meshioplusplus_bench_ops --allocs` (roadmap 5.28.1).
 *
 * Linked into `meshioplusplus_bench_ops` only when configured with
 * `-DMESHIOPLUSPLUS_BENCH_ALLOCS=ON`; never linked into the library or into
 * the default instruction-count benchmark binary.
 */

// System includes
#include <cstdint>

namespace meshioplusplus::bench {

struct AllocStats {
    std::uint64_t mAllocs = 0;
    std::uint64_t mTotalBytes = 0;
    std::uint64_t mPeakLiveBytes = 0;
};

void alloc_counter_start() noexcept;
AllocStats alloc_counter_stop() noexcept;

}  // namespace meshioplusplus::bench
