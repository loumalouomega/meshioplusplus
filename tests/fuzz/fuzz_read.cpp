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

/**
 * @file fuzz_read.cpp
 * @brief libFuzzer harness over every native registry reader (doc/fuzzing.md).
 *
 * The registry is path-based, so each input is written to a per-process
 * scratch file named the way the format expects to be found (its default
 * extension, or a fixed basename for `z88`, `d3plot`, `binout`, ...), then
 * read through `registry_readers()` -- the same entry point the native CLI,
 * the C API and WASM call.
 *
 * The contract under test: a reader either returns a mesh or throws
 * `ReadError` (every parser `std::` exception is translated to one by
 * detail/read_guard.hpp). Anything else is a finding and aborts: a crash, a
 * sanitizer report, a stray `std::exception` or a `std::bad_alloc` (a header
 * count trusted before the bytes that should back it). Set
 * `MIO_FUZZ_CRASH_ONLY=1` to tolerate every C++ exception while triaging
 * memory errors alone.
 *
 * The format is taken from, in order: `MIO_FUZZ_FORMAT`, a `-format=<name>`
 * argument, or the binary's name (`meshioplusplus_fuzz_read_<name>`).
 */

// System includes
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <map>
#include <new>
#include <string>
#include <typeinfo>

#include <unistd.h>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

std::string gFormat;
std::filesystem::path gDir;
std::filesystem::path gPath;
bool gCrashOnly = false;

/** @brief The file name a format's reader is found under. */
std::string fuzz_file_name(const std::string& rFormat) {
    // Formats found by basename, not extension (resolve_format / sniff_format).
    static const std::map<std::string, std::string> fixed = {
        {"z88", "z88i1.txt"},        {"lsdyna_d3plot", "d3plot"}, {"lsdyna_binout", "binout"},
        {"radioss_anim", "runA001"}, {"radioss_th", "runT01"},    {"ansys_rst_cyclic", "input.rst"},
    };
    if (auto it = fixed.find(rFormat); it != fixed.end())
        return it->second;
    for (const auto& [ext, fmt] : meshioplusplus::registry_extension_defaults())
        if (fmt == rFormat)
            return "input" + ext;
    return "input.dat";
}

std::string fuzz_format_from_args(int argc, char** argv) {
    if (const char* env = std::getenv("MIO_FUZZ_FORMAT"); env && *env)
        return env;
    for (int i = 1; i < argc; ++i)
        if (std::strncmp(argv[i], "-format=", 8) == 0)
            return argv[i] + 8;
    if (argc > 0) {
        const std::string name = std::filesystem::path(argv[0]).filename().string();
        const std::string prefix = "meshioplusplus_fuzz_read_";
        if (name.rfind(prefix, 0) == 0)
            return name.substr(prefix.size());
    }
    return {};
}

void fuzz_cleanup() {
    std::error_code ec;
    std::filesystem::remove_all(gDir, ec);
}

}  // namespace

/**
 * @brief Picks the format and creates the per-process scratch directory.
 * @return 0; exits with status 2 on an unknown or compiled-out format.
 */
extern "C" int LLVMFuzzerInitialize(int* pArgc, char*** pArgv) {
    gFormat = fuzz_format_from_args(*pArgc, *pArgv);
    const auto& readers = meshioplusplus::registry_readers();
    if (gFormat.empty() || !readers.count(gFormat)) {
        std::fprintf(stderr,
                     "meshioplusplus_fuzz_read: unknown format '%s'; known:", gFormat.c_str());
        for (const auto& kv : readers)
            std::fprintf(stderr, " %s", kv.first.c_str());
        std::fprintf(stderr, "\n");
        std::exit(2);
    }
    const char* crash_only = std::getenv("MIO_FUZZ_CRASH_ONLY");
    gCrashOnly = crash_only && *crash_only && *crash_only != '0';
    gDir = std::filesystem::temp_directory_path() /
           ("mio-fuzz-" + std::to_string(static_cast<long>(::getpid())));
    std::filesystem::create_directories(gDir);
    gPath = gDir / fuzz_file_name(gFormat);
    std::atexit(fuzz_cleanup);
    return 0;
}

/** @brief Reads one input through the registry; aborts on anything but ReadError. */
extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* pData, std::size_t size) {
    {
        std::FILE* f = std::fopen(gPath.c_str(), "wb");
        if (!f)
            std::abort();
        if (size)
            std::fwrite(pData, 1, size, f);
        std::fclose(f);
    }
    try {
        (void)meshioplusplus::registry_readers().at(gFormat)(gPath.string());
    } catch (const meshioplusplus::ReadError&) {
        // The one expected way to refuse an input.
    } catch (const std::exception& e) {
        if (!gCrashOnly) {
            std::fprintf(stderr, "meshioplusplus_fuzz_read: %s reader leaked %s: %s\n",
                         gFormat.c_str(), typeid(e).name(), e.what());
            std::abort();
        }
    }
    return 0;
}
