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
 * @file detail/dynlib.hpp
 * @brief Portable dynamic loading (dlopen / LoadLibrary) for provider plugins.
 *
 * Core-internal. Compiled out under Emscripten, where the functions report that
 * loading is unavailable.
 */

// System includes
#include <string>

namespace meshioplusplus::detail {

/// Whether this build can load shared libraries at run time.
bool dynlib_available() noexcept;

/// Opens the shared library at @p rPath. Returns null and sets @p rError on failure.
void* dynlib_open(const std::string& rPath, std::string& rError);

/// The address of @p pName in @p pHandle, or null.
void* dynlib_symbol(void* pHandle, const char* pName) noexcept;

/// Closes a library opened by `dynlib_open`.
void dynlib_close(void* pHandle) noexcept;

/// Makes the symbols of the image that contains this code visible to libraries
/// opened afterwards, which resolve the core's functions from the host. An
/// executable linked with -rdynamic already exports them; a shared library the
/// application loaded with RTLD_LOCAL (ctypes, Julia's Libdl) does not, so it is
/// promoted to global scope. Best effort; a no-op where unsupported.
void dynlib_export_host_symbols() noexcept;

/// The file extension of a shared library on this platform (".so", ".dylib", ".dll").
const char* dynlib_extension() noexcept;

}  // namespace meshioplusplus::detail
