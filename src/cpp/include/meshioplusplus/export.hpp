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
 * @file export.hpp
 * @brief `MESHIOPLUSPLUS_API`, the visibility/import-export attribute for the
 * installable C++ library.
 *
 * The C++ core targets (`meshioplusplus::core`, `meshioplusplus::core_meshio`,
 * `::core_native`, `::core_kratos`) are built with `CXX_VISIBILITY_PRESET hidden`
 * and `VISIBILITY_INLINES_HIDDEN ON`, so nothing is exported unless it says so.
 * Every entity **declared in a header and defined out of line** in
 * `src/cpp/src/` therefore needs this attribute; templates, `inline`
 * functions and anything else the consumer's own translation unit compiles do
 * not (and must not) carry it.
 *
 * A missing annotation is a *link* error for a consumer of a shared build
 * ("undefined symbol" / "unresolved external symbol"), never a silent runtime
 * fault -- which is what makes the annotation set verifiable rather than merely
 * asserted. Static builds resolve everything inside the archive and are
 * unaffected either way.
 *
 * The macro mirrors `MIO_API` in the C header (`bindings/c/include/`):
 * `MESHIOPLUSPLUS_CORE_SHARED` marks "this library is/was built shared" and is
 * attached INTERFACE by CMake so consumers inherit it, while
 * `MESHIOPLUSPLUS_CORE_BUILDING` is PRIVATE to the library's own compilation and
 * is what flips dllexport to dllimport on Windows.
 */

#if defined(_WIN32) && defined(MESHIOPLUSPLUS_CORE_SHARED)
#ifdef MESHIOPLUSPLUS_CORE_BUILDING
#define MESHIOPLUSPLUS_API __declspec(dllexport)
#else
#define MESHIOPLUSPLUS_API __declspec(dllimport)
#endif
#elif defined(__GNUC__)
#define MESHIOPLUSPLUS_API __attribute__((visibility("default")))
#else
#define MESHIOPLUSPLUS_API
#endif
