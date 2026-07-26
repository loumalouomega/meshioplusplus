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

namespace meshioplusplus::detail {

/**
 * @brief Defined once per backend library; referenced by every TU that includes
 * `mesh.hpp`. Never called -- only its address is taken, so this costs a
 * relocation rather than a static initializer.
 */
MESHIOPLUSPLUS_API void MESHIOPLUSPLUS_BACKEND_SYM(MESHIOPLUSPLUS_ACTIVE_BACKEND)();

#ifndef MESHIOPLUSPLUS_NO_BACKEND_LINK_CHECK
/**
 * @brief Constant-initialized reference that forces the symbol above to resolve.
 *
 * `[[maybe_unused]]` keeps `-Wunused` quiet; the pointer is still emitted
 * because taking a function's address odr-uses it.
 */
[[maybe_unused]] inline void (*const mesh_backend_link_check)() =
    &MESHIOPLUSPLUS_BACKEND_SYM(MESHIOPLUSPLUS_ACTIVE_BACKEND);
#endif

}  // namespace meshioplusplus::detail
