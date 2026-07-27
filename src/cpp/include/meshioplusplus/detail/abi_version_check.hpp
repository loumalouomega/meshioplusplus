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
 * @file detail/abi_version_check.hpp
 * @brief Makes a C++ ABI-version mismatch across the install boundary a **link**
 * error instead of silent undefined behaviour.
 *
 * The sibling of `detail/mesh_backend_check.hpp`, for the second axis on which a
 * consumer and the library can silently disagree. That header asks "were these
 * two compiled with the same *backend*?"; this one asks "were these two compiled
 * against the same *headers*?".
 *
 * The failure it catches is the one the version pin in `doc/cpp_api.md` is the
 * only other defence against: a consumer whose translation units were compiled
 * against one meshio++'s headers linking against a different meshio++'s library.
 * `Mesh`, `ModelPart`, every `operations/*` options struct and much else are
 * header-defined, so the two then disagree about layout -- an ODR violation that
 * typically shows up as unrelated memory corruption, never as a diagnosable
 * error. It happens in practice via a relaxed `find_package` version request, a
 * distro that rebuilt the library in place, or a pkg-config /
 * `compile_commands.json` setup that picked up an include path and a `-L` from
 * different prefixes.
 *
 * The mechanism is exactly the backend guard's: the library defines a single
 * `abi_version_is_<N>()` whose *name encodes* `MESHIOPLUSPLUS_ABI_VERSION`
 * (`abi_version.hpp`), and every TU that includes `mesh.hpp` references the one
 * its own headers name. Disagree, and the link fails naming the version the TU
 * expected.
 *
 * ### Deliberately a separate symbol from the backend guard
 *
 * Fusing the two into one `meshioplusplus_kratos_v3()` would halve the
 * relocations and ruin the diagnostic: a consumer with the right backend and
 * stale headers would get a message about the backend, which is the one thing
 * that *was* right. Two axes, two symbols, two messages. The honest cost is two
 * relocations per TU instead of one -- still no static initializer, which is the
 * property that matters.
 *
 * ### What it does NOT catch
 *
 * Only the meshio++-vs-meshio++ axis. A consumer built with a different compiler,
 * standard library, `_GLIBCXX_USE_CXX11_ABI`, or set of `-D`s affecting the
 * standard library is equally broken and equally undiagnosable here; that is
 * outside any single library's reach and is why `doc/abi.md` states the
 * same-toolchain precondition explicitly.
 *
 * Cost is one constant-initialized function pointer and one relocation per TU.
 * Define `MESHIOPLUSPLUS_NO_ABI_VERSION_CHECK` to opt out.
 *
 * @note MSVC + a **shared** meshio++ is a documented gap, identically to the
 *       backend guard and for the identical reason. The MSVC arm below is
 *       `#pragma detect_mismatch`, whose `/FAILIFMISMATCH` records live in the
 *       `.obj` files; they are not reliably carried through a DLL's import
 *       library, so there this guard degrades to no check at all. Static MSVC
 *       builds and every GNU/Clang configuration are covered. On that
 *       configuration the `SOVERSION`-derived DLL name and the version pin are
 *       the remaining defences.
 */

// Project includes
#include "meshioplusplus/abi_version.hpp"
#include "meshioplusplus/export.hpp"

// Two levels, so MESHIOPLUSPLUS_ABI_VERSION expands to its number before ## and
// # see it. Without the indirection the symbol would literally be
// `abi_version_is_MESHIOPLUSPLUS_ABI_VERSION`.
#define MESHIOPLUSPLUS_ABI_SYM_(n) abi_version_is_##n
#define MESHIOPLUSPLUS_ABI_SYM(n) MESHIOPLUSPLUS_ABI_SYM_(n)

#define MESHIOPLUSPLUS_ABI_STR_(n) #n
#define MESHIOPLUSPLUS_ABI_STR(n) MESHIOPLUSPLUS_ABI_STR_(n)

namespace meshioplusplus::detail {

/**
 * @brief Defined once by the library; referenced by every TU that includes
 * `mesh.hpp`. Never called -- only its address is taken, so this costs a
 * relocation rather than a static initializer.
 */
MESHIOPLUSPLUS_API void MESHIOPLUSPLUS_ABI_SYM(MESHIOPLUSPLUS_ABI_VERSION)();

#ifndef MESHIOPLUSPLUS_NO_ABI_VERSION_CHECK
#if defined(__GNUC__) || defined(__clang__)
/**
 * @brief Constant-initialized reference that forces the symbol above to resolve.
 *
 * `gnu::used` is **load-bearing, not decoration** -- do not "tidy" it away. This
 * is the lesson `detail/mesh_backend_check.hpp` paid for over a full release: an
 * `inline` variable has vague linkage and is emitted *lazily*, so with nothing
 * reading the pointer no relocation ever reaches the object file and the entire
 * guard is silently inert. That was the backend guard's state through v9.1.0.
 * Dropping `inline`/`const` does not fix it either -- a plain `static const`
 * initializer is still discarded at `-O2`. `nm -uC` on a consumer object file
 * must show this reference; CI asserts exactly that rather than trusting it.
 */
[[maybe_unused, gnu::used]] inline void (*const abi_version_link_check)() =
    &MESHIOPLUSPLUS_ABI_SYM(MESHIOPLUSPLUS_ABI_VERSION);
#elif defined(_MSC_VER)
// MSVC has no `used`, so the pointer above cannot be forced into the object.
// detect_mismatch reaches the same failure at the same moment by a different
// route: a /FAILIFMISMATCH record in every .obj, which the linker rejects when
// two disagree -- naming BOTH values, which beats an undefined symbol. See the
// DLL caveat in the file comment above.
#pragma detect_mismatch("meshioplusplus_abi_version", \
                        MESHIOPLUSPLUS_ABI_STR(MESHIOPLUSPLUS_ABI_VERSION))
#endif
#endif

}  // namespace meshioplusplus::detail
