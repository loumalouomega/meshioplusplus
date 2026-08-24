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
 * @file detail/provenance.hpp
 * @brief The single owner of the one-line provenance credit every format
 * writer emits.
 *
 * Before this header existed, ~25 writers spelled their own version of this
 * line inline, and the two engines drifted three ways at once: the C++ side
 * said `(C++ core)` and carried no version, the Python side said `v{version}`
 * and carried no engine marker, and several Python writers still said
 * `meshio` (a stale-fork artifact) rather than `meshio++`. The Python twin
 * lives at `meshioplusplus/_provenance.py`; the two must keep emitting
 * character-identical text (`tests/python/test_provenance.py` pins this),
 * which is only possible if both sides read the tag from one place rather
 * than composing it locally.
 *
 * This is deliberately the *whole* feature for now — see `doc/roadmap.md`
 * §1's "audit and normalize" bullet. It carries the writer version and
 * nothing else: no source format, no operation chain, no timestamp. Those
 * are separate, larger design questions the roadmap defers on purpose.
 */

// Project includes
#include "meshioplusplus/version.hpp"

namespace meshioplusplus {
namespace detail {

/// The canonical one-line provenance tag every writer emits, wrapped in each
/// format's own comment syntax at each writer's existing header position.
inline constexpr const char* kProvenanceTag = "Written by meshio++ v" MESHIOPLUSPLUS_VERSION_STRING;

}  // namespace detail
}  // namespace meshioplusplus
