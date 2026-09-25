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
 * @file replay_main.cpp
 * @brief A main() for the fuzz harness that replays files instead of fuzzing.
 *
 * `meshioplusplus_fuzz_replay -format=<name> <file-or-dir>...` feeds every
 * file to LLVMFuzzerTestOneInput once, so any compiler (GCC included) can run
 * the committed regression inputs under ctest, and a crash found by the
 * fuzzer can be debugged without libFuzzer. A format this build compiled out
 * (an HDF5 or netCDF one) is skipped with exit status 0.
 * `meshioplusplus_fuzz_replay -list-formats` prints this build's readers, one
 * per line (tools/fuzz/seed_corpus.py and the fuzz workflow use it).
 */

// System includes
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/registry.hpp"

extern "C" int LLVMFuzzerInitialize(int* pArgc, char*** pArgv);
extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* pData, std::size_t size);

namespace {

void replay_one(const std::filesystem::path& rPath) {
    std::vector<char> bytes;
    if (std::FILE* f = std::fopen(rPath.string().c_str(), "rb")) {
        char buf[65536];
        for (std::size_t n; (n = std::fread(buf, 1, sizeof buf, f)) > 0;)
            bytes.insert(bytes.end(), buf, buf + n);
        std::fclose(f);
    }
    std::printf("replay %s (%zu bytes)\n", rPath.string().c_str(), bytes.size());
    std::fflush(stdout);
    LLVMFuzzerTestOneInput(reinterpret_cast<const std::uint8_t*>(bytes.data()), bytes.size());
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::strcmp(argv[1], "-list-formats") == 0) {
        for (const auto& kv : meshioplusplus::registry_readers())
            std::printf("%s\n", kv.first.c_str());
        return 0;
    }
    for (int i = 1; i < argc; ++i) {
        if (std::strncmp(argv[i], "-format=", 8) == 0 &&
            meshioplusplus::registry_compiled_out(argv[i] + 8)) {
            std::printf("format '%s' is compiled out of this build; skipped\n", argv[i] + 8);
            return 0;
        }
    }
    LLVMFuzzerInitialize(&argc, &argv);
    std::size_t count = 0;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-')
            continue;
        const std::filesystem::path p(argv[i]);
        if (std::filesystem::is_directory(p)) {
            std::vector<std::filesystem::path> files;
            for (const auto& e : std::filesystem::directory_iterator(p))
                if (e.is_regular_file())
                    files.push_back(e.path());
            std::sort(files.begin(), files.end());
            for (const auto& f : files) {
                replay_one(f);
                ++count;
            }
        } else {
            replay_one(p);
            ++count;
        }
    }
    std::printf("replayed %zu input(s) without a finding\n", count);
    return 0;
}
