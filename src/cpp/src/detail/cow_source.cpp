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

// System includes
#include <cstdlib>
#include <filesystem>
#include <system_error>
#include <utility>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/read_options.hpp"

// Project includes (private, not installed)
#include "cow_source.hpp"

#if defined(__EMSCRIPTEN__)
// No mapping under Emscripten: the virtual FS has nothing to map.
#elif defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace meshioplusplus {
namespace detail {

namespace {

/// The size at which a file is mapped: `FileSource`'s `Auto` rule.
std::size_t cow_threshold() {
    static const std::size_t threshold = [] {
        if (const char* env = std::getenv("MESHIOPLUSPLUS_MMAP_THRESHOLD")) {
            char* end = nullptr;
            const unsigned long long v = std::strtoull(env, &end, 10);
            if (end && end != env)
                return static_cast<std::size_t>(v);
        }
        return mmap_auto_threshold_bytes;
    }();
    return threshold;
}

}  // namespace

CowSource::CowSource(const std::string& rPath, const std::string& rMessage,
                     std::size_t MapAtOrAbove) {
    if (TryMap(rPath, MapAtOrAbove ? MapAtOrAbove : cow_threshold()))
        return;
    std::error_code ec;
    const auto status = std::filesystem::status(rPath, ec);
    if (ec || !std::filesystem::is_regular_file(status))
        throw ReadError(rMessage);
    const std::uintmax_t size = std::filesystem::file_size(rPath, ec);
    if (ec)
        throw ReadError(rMessage);
    auto in = make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError(rMessage);
    mBuffer.assign(static_cast<std::size_t>(size) + 1, '\0');
    if (size > 0 && !in.read(mBuffer.data(), static_cast<std::streamsize>(size)))
        throw ReadError(rMessage);
    mSize = static_cast<std::size_t>(size);
}

CowSource::~CowSource() {
    Release();
}

CowSource::CowSource(CowSource&& rOther) noexcept {
    MoveFrom(std::move(rOther));
}

CowSource& CowSource::operator=(CowSource&& rOther) noexcept {
    if (this != &rOther) {
        Release();
        MoveFrom(std::move(rOther));
    }
    return *this;
}

bool CowSource::TryMap(const std::string& rPath, std::size_t MapAtOrAbove) {
#if defined(__EMSCRIPTEN__)
    (void)rPath;
    (void)MapAtOrAbove;
    return false;
#elif defined(_WIN32)
    HANDLE file = CreateFileA(rPath.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
                              FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file == INVALID_HANDLE_VALUE)
        return false;
    LARGE_INTEGER size;
    if (!GetFileSizeEx(file, &size) || size.QuadPart <= 0 ||
        static_cast<unsigned long long>(size.QuadPart) < MapAtOrAbove) {
        CloseHandle(file);
        return false;
    }
    // A copy-on-write view: writes go to private pages, never to the file.
    HANDLE mapping = CreateFileMappingA(file, nullptr, PAGE_WRITECOPY, 0, 0, nullptr);
    CloseHandle(file);  // the mapping keeps its own reference
    if (!mapping)
        return false;
    void* view = MapViewOfFile(mapping, FILE_MAP_COPY, 0, 0, 0);
    if (!view) {
        CloseHandle(mapping);
        return false;
    }
    mpData = static_cast<char*>(view);
    mSize = static_cast<std::size_t>(size.QuadPart);
    mpWinMapping = mapping;
    return true;
#else
    const int fd = ::open(rPath.c_str(), O_RDONLY);
    if (fd < 0)
        return false;
    struct stat st{};
    if (::fstat(fd, &st) != 0 || !S_ISREG(st.st_mode) || st.st_size <= 0 ||
        static_cast<unsigned long long>(st.st_size) < MapAtOrAbove) {
        ::close(fd);
        return false;
    }
    // MAP_PRIVATE with PROT_WRITE: a page is copied the first time it is
    // written, and the file never sees the write (a read-only descriptor
    // suffices for a private mapping).
    void* addr = ::mmap(nullptr, static_cast<std::size_t>(st.st_size), PROT_READ | PROT_WRITE,
                        MAP_PRIVATE, fd, 0);
    ::close(fd);
    if (addr == MAP_FAILED)
        return false;
#if defined(POSIX_MADV_SEQUENTIAL)
    ::posix_madvise(addr, static_cast<std::size_t>(st.st_size), POSIX_MADV_SEQUENTIAL);
#endif
    mpData = static_cast<char*>(addr);
    mSize = static_cast<std::size_t>(st.st_size);
    return true;
#endif
}

void CowSource::Release() noexcept {
#if defined(__EMSCRIPTEN__)
#elif defined(_WIN32)
    if (mpData)
        UnmapViewOfFile(mpData);
    if (mpWinMapping)
        CloseHandle(static_cast<HANDLE>(mpWinMapping));
#else
    if (mpData && mSize)
        ::munmap(mpData, mSize);
#endif
    mpData = nullptr;
    mpWinMapping = nullptr;
    mSize = 0;
    mBuffer.assign(1, '\0');
}

void CowSource::MoveFrom(CowSource&& rOther) noexcept {
    mBuffer = std::move(rOther.mBuffer);
    rOther.mBuffer.assign(1, '\0');
    mpData = rOther.mpData;
    mSize = rOther.mSize;
    mpWinMapping = rOther.mpWinMapping;
    rOther.mpData = nullptr;
    rOther.mpWinMapping = nullptr;
    rOther.mSize = 0;
}

}  // namespace detail
}  // namespace meshioplusplus
