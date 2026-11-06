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
 * @file detail/parallel_workers.hpp
 * @brief How many threads the build's parallel backend runs `parallel_for` on
 * right now. A **core-private** header: for the chunk counts of algorithms
 * whose result does not depend on them (`slot_runs.hpp`'s counting sort, the
 * chunked number parse), and for deciding that a parallel path cannot pay.
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <thread>

// Project includes
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

/// The workers `parallel_for` runs on: 1 on SEQ.
inline std::size_t parallel_workers() {
#if defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    return static_cast<std::size_t>(std::max(1, omp_get_max_threads()));
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    return std::max<std::size_t>(
        1, tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism));
#elif defined(MESHIOPLUSPLUS_PARALLEL_KOKKOS)
    kokkos_ensure_initialized();
    return static_cast<std::size_t>(std::max(1, Kokkos::DefaultHostExecutionSpace().concurrency()));
#elif defined(MESHIOPLUSPLUS_PARALLEL_STL)
    return std::max<std::size_t>(1, std::thread::hardware_concurrency());
#elif defined(MESHIOPLUSPLUS_PARALLEL_THREADS)
    return std::max<std::size_t>(1, threads_concurrency());
#else
    return 1;
#endif
}

}  // namespace detail
}  // namespace meshioplusplus
