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
 * @file detail/always_inline.hpp
 * @brief `MIO_ALWAYS_INLINE`: an inline function the compiler must inline.
 *
 * A **core-private** header. A shared helper that replaced a copy can still
 * change how the compiler inlines the rest of its caller (roadmap 3.4.11: a
 * share measured 2-4% more instructions on its rows). Forcing the helper
 * inline puts its body back into each caller before the size-driven inlining
 * decisions run, so the caller compiles as it did with the copy. Use it only
 * on a helper that replaced a copy, and keep it only where the
 * [instruction-count gate](../../../../doc/benchmarks.md#instruction-count-gate)
 * shows it is needed.
 */

#ifndef MIO_ALWAYS_INLINE
#if defined(_MSC_VER)
#define MIO_ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define MIO_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define MIO_ALWAYS_INLINE inline
#endif
#endif

/**
 * `MIO_ALWAYS_INLINE_LAMBDA`, placed after a lambda's parameter list: the same
 * for a lambda's call operator, for a loop body whose inlining into the
 * `parallel_for` that runs it the gate showed a header change undid.
 */
#ifndef MIO_ALWAYS_INLINE_LAMBDA
#if defined(__GNUC__) || defined(__clang__)
#define MIO_ALWAYS_INLINE_LAMBDA __attribute__((always_inline))
#else
#define MIO_ALWAYS_INLINE_LAMBDA
#endif
#endif
