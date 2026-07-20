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
#include <fstream>
#include <utility>

// Project includes
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/exceptions.hpp"

#if defined(__EMSCRIPTEN__)
// No mapping under Emscripten: the virtual FS has nothing to map.
#elif defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
// Only this translation unit (and the IMPLEMENTATION-guarded slice of the
// amalgamated header) ever includes <windows.h> for FileSource, but NOMINMAX
// is still required: without it, <windows.h> defines min/max macros that
// would break any std::min/std::max call later in the same TU.
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

FileSource::FileSource(const std::string& rPath, Mode mode) {
    if (mode != Mode::Buffered && TryMap(rPath, mode == Mode::Mmap))
        return;
    LoadBuffered(rPath);
}

FileSource::FileSource(const std::string& rPath, MmapMode mmap_mode)
    : FileSource(rPath, FromMmapMode(mmap_mode)) {}

FileSource::~FileSource() {
    Release();
}

FileSource::FileSource(FileSource&& rOther) noexcept {
    MoveFrom(std::move(rOther));
}

FileSource& FileSource::operator=(FileSource&& rOther) noexcept {
    if (this != &rOther) {
        Release();
        MoveFrom(std::move(rOther));
    }
    return *this;
}

FileSource::Mode FileSource::FromMmapMode(MmapMode mmap_mode) {
    switch (mmap_mode) {
        case MmapMode::On:
            return Mode::Mmap;
        case MmapMode::Off:
            return Mode::Buffered;
        default:
            return Mode::Auto;
    }
}

void FileSource::LoadBuffered(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    in.seekg(0, std::ios::end);
    const std::streamoff len = in.tellg();
    in.seekg(0, std::ios::beg);
    if (len > 0) {
        mBuffer.resize(static_cast<std::size_t>(len));
        in.read(mBuffer.data(), len);
    }
    mSize = mBuffer.size();
    mpData = nullptr;
    mMapped = false;
}

bool FileSource::TryMap(const std::string& rPath, bool force) {
#if defined(__EMSCRIPTEN__)
    (void)rPath;
    (void)force;
    return false;  // nothing to map on the virtual filesystem
#elif defined(_WIN32)
    HANDLE file = CreateFileA(rPath.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
                              FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file == INVALID_HANDLE_VALUE)
        return false;
    LARGE_INTEGER size;
    if (!GetFileSizeEx(file, &size) || size.QuadPart <= 0 || !WorthMapping(size.QuadPart, force)) {
        CloseHandle(file);
        return false;
    }
    HANDLE mapping = CreateFileMappingA(file, nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (!mapping) {
        CloseHandle(file);
        return false;
    }
    void* view = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
    if (!view) {
        CloseHandle(mapping);
        CloseHandle(file);
        return false;
    }
    mpData = static_cast<const char*>(view);
    mSize = static_cast<std::size_t>(size.QuadPart);
    mWinFile = file;
    mWinMapping = mapping;
    mMapped = true;
    return true;
#else
    const int fd = ::open(rPath.c_str(), O_RDONLY);
    if (fd < 0)
        return false;
    struct stat st{};
    // Only regular files: a pipe or character device has no meaningful size
    // and cannot be mapped.
    if (::fstat(fd, &st) != 0 || !S_ISREG(st.st_mode) || st.st_size <= 0 ||
        !WorthMapping(static_cast<long long>(st.st_size), force)) {
        ::close(fd);
        return false;
    }
    void* addr =
        ::mmap(nullptr, static_cast<std::size_t>(st.st_size), PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);  // the mapping keeps its own reference to the file
    if (addr == MAP_FAILED)
        return false;
#if defined(POSIX_MADV_SEQUENTIAL)
    // Parsers walk the buffer front to back.
    ::posix_madvise(addr, static_cast<std::size_t>(st.st_size), POSIX_MADV_SEQUENTIAL);
#endif
    mpData = static_cast<const char*>(addr);
    mSize = static_cast<std::size_t>(st.st_size);
    mMapped = true;
    return true;
#endif
}

bool FileSource::HasTerminatorSlack(long long size) {
#if defined(__EMSCRIPTEN__)
    (void)size;
    return false;
#elif defined(_WIN32)
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    const long long page = static_cast<long long>(info.dwPageSize);
    return page > 0 && (size % page) != 0;
#else
    const long long page = static_cast<long long>(::sysconf(_SC_PAGESIZE));
    return page > 0 && (size % page) != 0;
#endif
}

bool FileSource::WorthMapping(long long size, bool force) {
    // Applies even to an explicit request: `Mmap` is a preference, and
    // silently trading a guarantee for it would be the wrong bargain.
    if (!HasTerminatorSlack(size))
        return false;
    if (force)
        return true;
    static const std::size_t threshold = [] {
        if (const char* env = std::getenv("MESHIOPLUSPLUS_MMAP_THRESHOLD")) {
            char* end = nullptr;
            const unsigned long long v = std::strtoull(env, &end, 10);
            if (end && end != env)
                return static_cast<std::size_t>(v);
        }
        return mmap_auto_threshold_bytes;
    }();
    return static_cast<std::size_t>(size) >= threshold;
}

void FileSource::Release() {
#if defined(__EMSCRIPTEN__)
#elif defined(_WIN32)
    if (mMapped && mpData)
        UnmapViewOfFile(mpData);
    if (mWinMapping)
        CloseHandle(static_cast<HANDLE>(mWinMapping));
    if (mWinFile && mWinFile != INVALID_HANDLE_VALUE)
        CloseHandle(static_cast<HANDLE>(mWinFile));
    mWinMapping = nullptr;
    mWinFile = nullptr;
#else
    if (mMapped && mpData && mSize)
        ::munmap(const_cast<char*>(mpData), mSize);
#endif
    mpData = nullptr;
    mSize = 0;
    mMapped = false;
}

void FileSource::MoveFrom(FileSource&& rOther) noexcept {
    mBuffer = std::move(rOther.mBuffer);
    mpData = rOther.mpData;
    mSize = rOther.mSize;
    mMapped = rOther.mMapped;
#if defined(_WIN32) && !defined(__EMSCRIPTEN__)
    mWinFile = rOther.mWinFile;
    mWinMapping = rOther.mWinMapping;
    rOther.mWinFile = nullptr;
    rOther.mWinMapping = nullptr;
#endif
    rOther.mpData = nullptr;
    rOther.mSize = 0;
    rOther.mMapped = false;
}

}  // namespace detail
}  // namespace meshioplusplus
