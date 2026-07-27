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
 * @file detail/mesh_backend_check.hpp
 * @brief Makes a mesh-backend mismatch across the install boundary a **link**
 * error instead of silent undefined behaviour.
 *
 * `meshioplusplus::Mesh` is a different type per backend (`mesh.hpp`), so a
 * translation unit compiled against one backend and linked against a library
 * built with another disagrees about the layout of nearly every argument it
 * passes. Nothing about that is diagnosable at run time -- it is an ODR
 * violation that usually manifests as unrelated memory corruption.
 *
 * Two guards, in order of who they catch:
 *
 *  1. `mesh.hpp` `#error`s when more than one `MESHIOPLUSPLUS_MESH_BACKEND_*`
 *     macro is defined. That catches a consumer who passes a `-D` on top of the
 *     one CMake already propagated.
 *  2. This header catches the harder case: **no** macro defined (so `mesh.hpp`
 *     silently falls through to its MESHIO default) while the library is
 *     NATIVE or KRATOS. Each backend's library defines exactly one
 *     `mesh_backend_is_<backend>()` symbol, and every TU including `mesh.hpp`
 *     references the one its own macros select. Disagree and the link fails
 *     naming the backend it expected.
 *
 * CMake consumers are already safe by construction -- the backend macro rides
 * in the exported target's `INTERFACE_COMPILE_DEFINITIONS`, so
 * `target_link_libraries(... meshioplusplus::core_kratos)` compiles the
 * consumer's own sources with `MESHIOPLUSPLUS_MESH_BACKEND_KRATOS` whether it
 * asked to or not. This guard is for everyone else: pkg-config, hand-written
 * makefiles, and IDE/compile_commands setups that pick up the include path but
 * not the definitions.
 *
 * Cost is one constant-initialized function pointer and one relocation per TU.
 * Define `MESHIOPLUSPLUS_NO_BACKEND_LINK_CHECK` to opt out.
 *
 * @note MSVC + a **shared** meshio++ is a documented gap. The MSVC arm below is
 *       `#pragma detect_mismatch`, whose `/FAILIFMISMATCH` records live in the
 *       `.obj` files; they are not reliably carried through a DLL's import
 *       library, so there the guard degrades to no check at all. Static MSVC
 *       builds and every GNU/Clang configuration are covered. A run-time probe
 *       would close it, at the price of a static initializer in every TU --
 *       exactly what this design exists to avoid.
 */

// Project includes
#include "meshioplusplus/export.hpp"

#if defined(MESHIOPLUSPLUS_MESH_BACKEND_NATIVE)
#define MESHIOPLUSPLUS_ACTIVE_BACKEND native
#elif defined(MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)
#define MESHIOPLUSPLUS_ACTIVE_BACKEND kratos
#else
#define MESHIOPLUSPLUS_ACTIVE_BACKEND meshio
#endif

// Two levels: the outer one expands MESHIOPLUSPLUS_ACTIVE_BACKEND before ## sees it.
#define MESHIOPLUSPLUS_BACKEND_SYM_(name) mesh_backend_is_##name
#define MESHIOPLUSPLUS_BACKEND_SYM(name) MESHIOPLUSPLUS_BACKEND_SYM_(name)

// Same two-level trick for the string form the MSVC arm needs.
#define MESHIOPLUSPLUS_BACKEND_STR_(name) #name
#define MESHIOPLUSPLUS_BACKEND_STR(name) MESHIOPLUSPLUS_BACKEND_STR_(name)

namespace meshioplusplus::detail {

/**
 * @brief Defined once per backend library; referenced by every TU that includes
 * `mesh.hpp`. Never called -- only its address is taken, so this costs a
 * relocation rather than a static initializer.
 */
MESHIOPLUSPLUS_API void MESHIOPLUSPLUS_BACKEND_SYM(MESHIOPLUSPLUS_ACTIVE_BACKEND)();

#ifndef MESHIOPLUSPLUS_NO_BACKEND_LINK_CHECK
#if defined(__GNUC__) || defined(__clang__)
/**
 * @brief Constant-initialized reference that forces the symbol above to resolve.
 *
 * `[[maybe_unused]]` keeps `-Wunused` quiet.
 *
 * `gnu::used` is **load-bearing, not decoration** -- do not "tidy" it away. An
 * `inline` variable has vague linkage and is emitted *lazily*: nothing in this
 * header or in `mesh.hpp` reads the pointer, so without `used` no relocation
 * ever reaches the object file and the entire guard is inert. That was the
 * state through v9.1.0; `nm -uC` on any consumer TU showed no reference at all,
 * with or without a backend macro, and a mismatched consumer got raw mangled
 * undefined references to whatever it actually called instead of a message
 * naming the backend. Dropping `inline`/`const` does not fix it either -- a
 * plain `static const` initializer is still discarded at `-O2`.
 */
[[maybe_unused, gnu::used]] inline void (*const mesh_backend_link_check)() =
    &MESHIOPLUSPLUS_BACKEND_SYM(MESHIOPLUSPLUS_ACTIVE_BACKEND);
#elif defined(_MSC_VER)
// MSVC has no `used`, so the pointer above cannot be forced into the object.
// `detect_mismatch` reaches the same failure at the same moment by a different
// route: it embeds a /FAILIFMISMATCH record in every .obj, and the linker
// rejects a link whose records disagree -- naming BOTH values, which is
// strictly more informative than an undefined symbol. It also needs no symbol
// at all, sidestepping the __declspec(dllimport) function-pointer question.
// See the DLL caveat in the file comment above.
#pragma detect_mismatch("meshioplusplus_mesh_backend", \
                        MESHIOPLUSPLUS_BACKEND_STR(MESHIOPLUSPLUS_ACTIVE_BACKEND))
#endif
#endif

}  // namespace meshioplusplus::detail
