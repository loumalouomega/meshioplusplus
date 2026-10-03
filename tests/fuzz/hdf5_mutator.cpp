// SPDX-License-Identifier: MIT
// Structure-aware mutations of valid HDF5 containers. The reader still sees
// ordinary file bytes, so findings replay with the unmodified replay harness.
// netCDF-4/Exodus is HDF5 too; classic netCDF falls back to byte mutation.
#include <algorithm>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include <hdf5.h>
#include <sanitizer/common_interface_defs.h>
#include <sys/wait.h>
#include <unistd.h>

extern "C" std::size_t LLVMFuzzerMutate(std::uint8_t*, std::size_t, std::size_t);

namespace {

// Visit hard-linked objects only (H5Ovisit does not follow soft/external links).
// Generated valid seeds are tiny; cap metadata collected from other inputs.
#if H5_VERSION_GE(1, 12, 0)
using HdfmutObjectInfo = H5O_info2_t;
#else
using HdfmutObjectInfo = H5O_info_t;
#endif

herr_t hdfmut_visit(hid_t, const char* pName, const HdfmutObjectInfo* pInfo, void* pUser) {
    auto& names = *static_cast<std::vector<std::string>*>(pUser);
    if (pInfo->type == H5O_TYPE_DATASET && names.size() < 256 && std::strlen(pName) < 1024)
        names.emplace_back(pName);
    return names.size() >= 256 ? 1 : 0;
}

bool hdfmut_change(hid_t File, std::mt19937& rRandom) {
    std::vector<std::string> names;
#if H5_VERSION_GE(1, 12, 0)
    const herr_t visited =
        H5Ovisit3(File, H5_INDEX_NAME, H5_ITER_INC, hdfmut_visit, &names, H5O_INFO_BASIC);
#elif H5_VERSION_GE(1, 10, 3)
    const herr_t visited =
        H5Ovisit2(File, H5_INDEX_NAME, H5_ITER_INC, hdfmut_visit, &names, H5O_INFO_BASIC);
#else
    const herr_t visited = H5Ovisit(File, H5_INDEX_NAME, H5_ITER_INC, hdfmut_visit, &names);
#endif
    if (visited < 0 || names.empty())
        return false;
    const std::string name = names[rRandom() % names.size()];
    const unsigned op = rRandom() % 6;
    if (op == 0) {
        // Missing objects, dangling links and cycles exercise library-object
        // walking without corrupting the HDF5 superblock/checksums.
        if (H5Ldelete(File, name.c_str(), H5P_DEFAULT) < 0)
            return false;
        return H5Lcreate_soft(rRandom() % 2 ? "/missing" : name.c_str(), File, name.c_str(),
                              H5P_DEFAULT, H5P_DEFAULT) >= 0;
    }
    hid_t dataset = H5Dopen2(File, name.c_str(), H5P_DEFAULT);
    if (dataset < 0)
        return false;
    hid_t space = H5Dget_space(dataset);
    hid_t type = H5Dget_type(dataset);
    const hssize_t count = H5Sget_simple_extent_npoints(space);
    bool changed = false;
    if (op == 1 && count > 0 && count <= 1024 && H5Tget_class(type) == H5T_INTEGER) {
        std::vector<std::int64_t> values(static_cast<std::size_t>(count));
        if (H5Dread(dataset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0) {
            static constexpr std::int64_t edges[] = {-1, 0, 1, 2, 4, 7, 16, 1024};
            values[rRandom() % values.size()] = edges[rRandom() % std::size(edges)];
            changed = H5Dwrite(dataset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               values.data()) >= 0;
        }
    } else if (op == 2) {
        // Attribute type/rank disagreement. Pick an existing attribute rather
        // than inventing one a reader will never inspect.
        const int attrs = H5Aget_num_attrs(dataset);
        if (attrs > 0) {
            hid_t attr =
                H5Aopen_by_idx(dataset, ".", H5_INDEX_NAME, H5_ITER_INC,
                               rRandom() % static_cast<unsigned>(attrs), H5P_DEFAULT, H5P_DEFAULT);
            char attr_name[1024];
            const ssize_t len = H5Aget_name(attr, sizeof attr_name, attr_name);
            H5Aclose(attr);
            if (len > 0 && static_cast<std::size_t>(len) < sizeof attr_name &&
                H5Adelete(dataset, attr_name) >= 0) {
                hid_t scalar = H5Screate(H5S_SCALAR);
                attr = H5Acreate2(dataset, attr_name, H5T_NATIVE_INT64, scalar, H5P_DEFAULT,
                                  H5P_DEFAULT);
                const std::int64_t value = static_cast<std::int64_t>(rRandom() % 17) - 1;
                changed = attr >= 0 && H5Awrite(attr, H5T_NATIVE_INT64, &value) >= 0;
                if (attr >= 0)
                    H5Aclose(attr);
                H5Sclose(scalar);
            }
        }
    } else {
        // Replace a dataset with a small mismatched rank/shape/type. Keep the
        // original type (including compounds) in one mode to reach Nastran's
        // table walkers. No dimensions/allocations are trusted from input.
        H5Dclose(dataset);
        dataset = -1;
        if (H5Ldelete(File, name.c_str(), H5P_DEFAULT) >= 0) {
            const hsize_t dims[] = {rRandom() % 9, rRandom() % 5};
            hid_t replacement = H5Screate_simple(op == 3 ? 2 : 1, dims, nullptr);
            dataset = H5Dcreate2(File, name.c_str(), op == 4 ? H5T_NATIVE_INT64 : type, replacement,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            changed = dataset >= 0;
            H5Sclose(replacement);
        }
    }
    if (dataset >= 0)
        H5Dclose(dataset);
    H5Tclose(type);
    H5Sclose(space);
    return changed;
}

}  // namespace

extern "C" std::size_t LLVMFuzzerCustomMutator(std::uint8_t* pData, std::size_t Size,
                                               std::size_t MaxSize, unsigned Seed) {
    // Some raw mutation is necessary for the libraries' own container parsers.
    if (Seed % 4 == 0 || Size < 8)
        return LLVMFuzzerMutate(pData, Size, MaxSize);
    // Only HDF5-based containers (HDF5 magic) reach the structure-aware path;
    // classic netCDF and other bytes fall back to raw mutation. The magic check
    // also rejects truncated inputs before any library call.
    static constexpr char kHdf5Magic[] = "\x89HDF\r\n\x1a\n";
    if (Size < 8 || std::memcmp(pData, kHdf5Magic, 8) != 0)
        return LLVMFuzzerMutate(pData, Size, MaxSize);
    const auto path = std::filesystem::temp_directory_path() /
                      ("mio-hdf-mutator-" + std::to_string(static_cast<long>(::getpid())));
    std::FILE* stream = std::fopen(path.c_str(), "wb");
    if (!stream)
        return LLVMFuzzerMutate(pData, Size, MaxSize);
    if (std::fwrite(pData, 1, Size, stream) != Size) {
        std::fclose(stream);
        std::filesystem::remove(path);
        return LLVMFuzzerMutate(pData, Size, MaxSize);
    }
    std::fclose(stream);
    // Isolate HDF5's own container parser: a corrupted input can abort inside
    // H5Ovisit/H5Dopen (dense link tables, huge object headers). A crash in the
    // child falls back to byte mutation instead of killing the fuzzer.
    const pid_t pid = ::fork();
    if (pid < 0) {
        std::filesystem::remove(path);
        return LLVMFuzzerMutate(pData, Size, MaxSize);
    }
    if (pid == 0) {
        // libFuzzer's inherited death callback writes the parent's current
        // input, which is empty while mutating. A mutator failure should only
        // trigger the byte-mutation fallback; reader failures still report.
        __sanitizer_set_death_callback(nullptr);
        std::signal(SIGABRT, SIG_DFL);
        std::signal(SIGALRM, SIG_DFL);
        H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
        hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        bool changed = false;
        if (file >= 0) {
            std::mt19937 random(Seed);
            changed = hdfmut_change(file, random);
            H5Fclose(file);
        }
        std::error_code ec;
        const auto length = std::filesystem::file_size(path, ec);
        _exit(changed && !ec && length <= MaxSize ? 0 : 1);
    }
    int status = 1;
    while (::waitpid(pid, &status, 0) < 0) {
    }
    bool changed = WIFEXITED(status) && WEXITSTATUS(status) == 0;
    std::error_code ec;
    const auto length = std::filesystem::file_size(path, ec);
    if (changed && !ec && length <= MaxSize) {
        stream = std::fopen(path.c_str(), "rb");
        if (stream && std::fread(pData, 1, static_cast<std::size_t>(length), stream) == length) {
            std::fclose(stream);
            std::filesystem::remove(path);
            return static_cast<std::size_t>(length);
        }
        if (stream)
            std::fclose(stream);
    }
    std::filesystem::remove(path);
    return LLVMFuzzerMutate(pData, Size, MaxSize);
}
