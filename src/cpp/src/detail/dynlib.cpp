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

// Portable dynamic loading. See detail/dynlib.hpp.

#include "dynlib.hpp"

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#elif !defined(__EMSCRIPTEN__)
#include <dlfcn.h>
#endif

namespace meshioplusplus::detail {

#if defined(__EMSCRIPTEN__)

bool dynlib_available() noexcept {
    return false;
}
void* dynlib_open(const std::string&, std::string& rError) {
    rError = "loading shared libraries is not available in WebAssembly";
    return nullptr;
}
void* dynlib_symbol(void*, const char*) noexcept {
    return nullptr;
}
void dynlib_close(void*) noexcept {}
void dynlib_export_host_symbols() noexcept {}
const char* dynlib_extension() noexcept {
    return ".wasm";
}

#elif defined(_WIN32)

bool dynlib_available() noexcept {
    return true;
}

void* dynlib_open(const std::string& rPath, std::string& rError) {
    // UTF-8 path to UTF-16.
    const int n = MultiByteToWideChar(CP_UTF8, 0, rPath.c_str(), -1, nullptr, 0);
    std::wstring wide(static_cast<std::size_t>(n > 0 ? n : 1), L'\0');
    if (n > 0)
        MultiByteToWideChar(CP_UTF8, 0, rPath.c_str(), -1, wide.data(), n);
    HMODULE h = LoadLibraryW(wide.c_str());
    if (!h)
        rError = "LoadLibrary failed with error " + std::to_string(GetLastError());
    return reinterpret_cast<void*>(h);
}

void* dynlib_symbol(void* pHandle, const char* pName) noexcept {
    return reinterpret_cast<void*>(GetProcAddress(reinterpret_cast<HMODULE>(pHandle), pName));
}

void dynlib_close(void* pHandle) noexcept {
    FreeLibrary(reinterpret_cast<HMODULE>(pHandle));
}

// A DLL resolves imports at link time, against a named import library; a plugin
// cannot borrow the host's symbols the way an ELF or Mach-O one does.
void dynlib_export_host_symbols() noexcept {}

const char* dynlib_extension() noexcept {
    return ".dll";
}

#else

bool dynlib_available() noexcept {
    return true;
}

void* dynlib_open(const std::string& rPath, std::string& rError) {
    void* h = dlopen(rPath.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (!h) {
        const char* p_msg = dlerror();
        rError = p_msg ? p_msg : "dlopen failed";
    }
    return h;
}

void* dynlib_symbol(void* pHandle, const char* pName) noexcept {
    return dlsym(pHandle, pName);
}

void dynlib_close(void* pHandle) noexcept {
    dlclose(pHandle);
}

void dynlib_export_host_symbols() noexcept {
    // The image holding this function: the executable (whose -rdynamic exports
    // already serve a plugin) or a shared library, which RTLD_NOLOAD | RTLD_GLOBAL
    // promotes without loading a second copy.
    Dl_info info;
    if (dladdr(reinterpret_cast<void*>(&dynlib_export_host_symbols), &info) && info.dli_fname)
        dlopen(info.dli_fname, RTLD_NOW | RTLD_GLOBAL | RTLD_NOLOAD);
}

const char* dynlib_extension() noexcept {
#if defined(__APPLE__)
    return ".dylib";
#else
    return ".so";
#endif
}

#endif

}  // namespace meshioplusplus::detail
