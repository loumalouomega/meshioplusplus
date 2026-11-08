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
 * @file detail/link_check.hpp
 * @brief The macros the two link-time sentinels share:
 * `detail/abi_version_check.hpp` (same headers?) and
 * `detail/mesh_backend_check.hpp` (same backend?).
 *
 * Each sentinel names one library symbol whose *name encodes* the value it
 * checks, and makes every translation unit reference it, so a disagreement
 * fails the link naming the value the unit expected. This header holds the
 * pieces both spell the same way: the two-level token paste and stringize, and
 * the GNU/Clang reference that forces the symbol to resolve. The MSVC arm
 * (`#pragma detect_mismatch`) stays in each sentinel, which names its own key.
 */

// Two levels, so a macro argument (`MESHIOPLUSPLUS_ABI_VERSION`,
// `MESHIOPLUSPLUS_ACTIVE_BACKEND`) expands before ## and # see it. Without the
// indirection the symbol would literally be `abi_version_is_MESHIOPLUSPLUS_ABI_VERSION`.
#define MESHIOPLUSPLUS_PP_CAT_(a, b) a##b
#define MESHIOPLUSPLUS_PP_CAT(a, b) MESHIOPLUSPLUS_PP_CAT_(a, b)

#define MESHIOPLUSPLUS_PP_STR_(x) #x
#define MESHIOPLUSPLUS_PP_STR(x) MESHIOPLUSPLUS_PP_STR_(x)

#if defined(__GNUC__) || defined(__clang__)
/**
 * @brief `MESHIOPLUSPLUS_LINK_CHECK_REF(Var, Sym)`: a constant-initialized
 * function pointer `Var` holding `&Sym`, which forces `Sym` to resolve at link
 * time. Never called, so it costs a relocation rather than a static initializer.
 *
 * `[[maybe_unused]]` keeps `-Wunused` quiet. `gnu::used` is **load-bearing, not
 * decoration** -- do not "tidy" it away. An `inline` variable has vague linkage
 * and is emitted *lazily*: nothing reads the pointer, so without `used` no
 * relocation ever reaches the object file and the whole guard is inert. That
 * was the backend guard's state through v9.1.0; `nm -uC` on any consumer TU
 * showed no reference at all, and a mismatched consumer got raw mangled
 * undefined references to whatever it called instead of a message naming the
 * backend. Dropping `inline`/`const` does not fix it either -- a plain
 * `static const` initializer is still discarded at `-O2`. `nm -uC` on a
 * consumer object file must show the reference; CI asserts exactly that.
 */
#define MESHIOPLUSPLUS_LINK_CHECK_REF(Var, Sym) \
    [[maybe_unused, gnu::used]] inline void (*const Var)() = &Sym
#endif
