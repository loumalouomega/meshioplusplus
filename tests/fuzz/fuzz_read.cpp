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
 * Companion-file readers take a deterministic bundle instead of a single file
 * (`fuzz_bundle.hpp`): XDMF XML/heavy-data sets (`.xdmf` plus `.h5`/`.bin`
 * payloads, including DataItem/reference reconstruction) and ADIOS2 `.bp`
 * directories. Findings stay ordinary bundle bytes that replay unchanged.
 *
 * Library-backed readers (HDF5/netCDF, ADIOS2, TecIO, CGNS MLL) run isolated
 * in a forked child: a library abort inside `H5Ovisit`/`H5Dopen`/`nc_open` (or
 * an ADIOS2/TecIO/cgnslib abort) then attributes clearly instead of corrupting
 * the fuzzer's own HDF5 global state. The parent re-aborts to preserve the
 * finding. Set `MIO_FUZZ_NO_ISOLATE=1` to debug without the fork.
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
#include <fstream>
#include <map>
#include <new>
#include <string>
#include <typeinfo>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

// Project includes
#include "fuzz_bundle.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/mdpa.hpp"
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

/** @brief Readers needing companion files: one fuzz input holds a whole bundle. */
bool fuzz_is_bundle_format(const std::string& rFormat) {
    return rFormat == "xdmf" || rFormat == "vtx";
}

/** @brief Library-backed readers run isolated so a library abort cannot corrupt the fuzzer. */
bool fuzz_is_isolated_format(const std::string& rFormat) {
    return rFormat == "cgns" || rFormat == "h5m" || rFormat == "hmf" || rFormat == "med" ||
           rFormat == "nastran_h5" || rFormat == "vtkhdf" || rFormat == "exodus" ||
           rFormat == "xdmf" || rFormat == "vtx" || rFormat == "szplt";
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

/**
 * @brief Reads one input through the registry; aborts on anything but ReadError.
 *
 * The registry drops the side channels (`MdpaInfo`), and the parsers behind
 * them -- `Properties` bodies, verbatim blocks, lenient skipping -- are
 * reachable only through the info overload, so `mdpa` is read through that
 * one too, strictly and leniently.
 *
 * Bundle formats (`xdmf`, `vtx`) accept either a single file (an XML-only
 * `.xdmf`, or a non-directory `.bp` the reader refuses) or a `MIOB` bundle
 * holding the entry file plus its companions. The bundle's main path is the
 * first `.xdmf` entry for XDMF, else the first entry; for VTX it is the first
 * entry's top-level directory (the `.bp` directory the reader opens).
 */
extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* pData, std::size_t size) {
    std::string main_path;
    const bool is_bundle_format = fuzz_is_bundle_format(gFormat);
    const bool is_bundle =
        is_bundle_format && size > 0 && meshioplusplus::fuzz_bundle_is_bundle(pData, size);
    if (is_bundle) {
        std::vector<meshioplusplus::FuzzBundleEntry> entries;
        try {
            entries = meshioplusplus::fuzz_bundle_decode(pData, size);
        } catch (const meshioplusplus::ReadError&) {
            return 0;  // malformed bundle: the mutator exploring structure
        }
        std::error_code ec;
        std::filesystem::remove_all(gDir, ec);
        std::filesystem::create_directories(gDir, ec);
        for (const auto& [name, data] : entries) {
            const std::filesystem::path out = gDir / name;
            std::filesystem::create_directories(out.parent_path(), ec);
            std::ofstream f(out, std::ios::binary);
            if (!f)
                return 0;
            if (!data.empty())
                f.write(reinterpret_cast<const char*>(data.data()),
                        static_cast<std::streamsize>(data.size()));
        }
        if (gFormat == "xdmf") {
            main_path = {};
            for (const auto& [name, data] : entries) {
                (void)data;
                if (name.size() >= 5 && name.compare(name.size() - 5, 5, ".xdmf") == 0) {
                    main_path = (gDir / name).string();
                    break;
                }
            }
            if (main_path.empty())
                main_path = (gDir / entries.front().first).string();
        } else {  // vtx: the first entry's top-level directory is the .bp to open
            const std::string& first = entries.front().first;
            const std::size_t slash = first.find('/');
            if (slash == std::string::npos)
                return 0;  // no directory: the reader refuses it as ReadError below
            main_path = (gDir / first.substr(0, slash)).string();
        }
    } else {
        if (is_bundle_format && gFormat == "vtx") {
            // A lone file is never a valid `.bp` directory: write it where the
            // reader looks so it refuses it as ReadError rather than crashing
            // on an uninitialized path.
            std::error_code ec;
            std::filesystem::remove_all(gDir, ec);
            std::filesystem::create_directories(gDir, ec);
            std::filesystem::create_directories(gPath, ec);
            main_path = gPath.string();
            const std::filesystem::path probe = gPath / "md.idx";
            std::FILE* f = std::fopen(probe.c_str(), "wb");
            if (!f)
                std::abort();
            if (size)
                std::fwrite(pData, 1, size, f);
            std::fclose(f);
        } else {
            std::FILE* f = std::fopen(gPath.c_str(), "wb");
            if (!f)
                std::abort();
            if (size)
                std::fwrite(pData, 1, size, f);
            std::fclose(f);
            main_path = gPath.string();
        }
    }

    const bool isolate = fuzz_is_isolated_format(gFormat);
    const char* no_isolate = std::getenv("MIO_FUZZ_NO_ISOLATE");
    if (isolate && !(no_isolate && *no_isolate && *no_isolate != '0')) {
        const pid_t pid = ::fork();
        if (pid < 0) {
            // Fork unavailable: fall through to a direct read.
        } else if (pid == 0) {
            try {
                (void)meshioplusplus::registry_readers().at(gFormat)(main_path);
            } catch (const meshioplusplus::ReadError&) {
                _exit(0);
            } catch (const std::exception&) {
                _exit(42);  // meshio defect: the parent re-aborts to keep the artifact
            }
            _exit(0);
        } else {
            int status = 0;
            while (::waitpid(pid, &status, 0) < 0) {
            }
            if (WIFEXITED(status)) {
                if (WEXITSTATUS(status) == 0)
                    return 0;
                std::fprintf(stderr,
                             "meshioplusplus_fuzz_read: %s reader leaked an exception on '%s'\n",
                             gFormat.c_str(), main_path.c_str());
                std::abort();
            }
            if (WIFSIGNALED(status)) {
                std::fprintf(stderr,
                             "meshioplusplus_fuzz_read: %s reader died on signal %d in '%s' "
                             "(library abort inside HDF5/netCDF/ADIOS2/TecIO/cgnslib or a "
                             "memory error; replay with meshioplusplus_fuzz_replay to triage)\n",
                             gFormat.c_str(), WTERMSIG(status), main_path.c_str());
                std::abort();
            }
            std::abort();
        }
    }

    try {
        (void)meshioplusplus::registry_readers().at(gFormat)(main_path);
        if (gFormat == "mdpa") {
            meshioplusplus::MdpaInfo info;
            (void)meshioplusplus::read_mdpa(main_path, info);
            meshioplusplus::ReadOptions lenient;
            lenient.mLenient = true;
            (void)meshioplusplus::read_mdpa(main_path, info, lenient);
        }
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
