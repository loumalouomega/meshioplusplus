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
 * @file thread_sweep.hpp
 * @brief The thread count of the backends whose count a test can set at run
 * time (OpenMP and THREADS), for the determinism tests that run an operation at
 * 1, 2 and 4 threads and compare the results.
 *
 * `MESHIOPLUSPLUS_TEST_THREAD_SWEEP` is 1 under those backends and 0 otherwise;
 * `set_threads`, `max_threads` and `num_procs` stand in for OpenMP's
 * `omp_set_num_threads`, `omp_get_max_threads` and `omp_get_num_procs`.
 */

// System includes
#include <algorithm>
#include <thread>

#if defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
#include <omp.h>
#endif

// Project includes
#include "meshioplusplus/parallel.hpp"

#if defined(MESHIOPLUSPLUS_PARALLEL_OPENMP) || defined(MESHIOPLUSPLUS_PARALLEL_THREADS)
#define MESHIOPLUSPLUS_TEST_THREAD_SWEEP 1
#else
#define MESHIOPLUSPLUS_TEST_THREAD_SWEEP 0
#endif

namespace mio_test {

/// Sets the thread count later parallel loops use.
inline void set_threads(int Count) {
#if defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    omp_set_num_threads(Count);
#elif defined(MESHIOPLUSPLUS_PARALLEL_THREADS)
    meshioplusplus::detail::threads_set_concurrency(static_cast<unsigned>(std::max(Count, 1)));
#else
    (void)Count;
#endif
}

/// The thread count parallel loops use now.
inline int max_threads() {
#if defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    return omp_get_max_threads();
#elif defined(MESHIOPLUSPLUS_PARALLEL_THREADS)
    return static_cast<int>(meshioplusplus::detail::threads_concurrency());
#else
    return 1;
#endif
}

/// The processors available.
inline int num_procs() {
#if defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    return omp_get_num_procs();
#else
    return static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
#endif
}

}  // namespace mio_test
