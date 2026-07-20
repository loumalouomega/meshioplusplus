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
 * @file file_source.hpp
 * @brief Whole-file read access, memory-mapped where that helps and buffered
 *        otherwise.
 *
 * Several readers (gmsh, vtk, ensight, ugrid's ASCII branch, openfoam) slurp
 * the entire file into a heap buffer and then parse from it. `FileSource`
 * replaces that idiom with one that can map the file instead, which removes a
 * full-file copy and -- more importantly -- the peak-RSS doubling that copy
 * causes on multi-GB meshes.
 *
 * **Be clear about what this does and does not buy.** Against one bulk `read()`
 * into a pre-sized buffer, mapping saves a `memcpy` and the second copy of the
 * file in RAM. It does *not* make the page faults go away, so this is a
 * memory-footprint feature far more than a throughput one. The clear exception
 * is a reader that was building an `ostringstream` and calling `.str()` --
 * that pays for two extra copies, and mapping removes both.
 *
 * **Mapping is advisory, never required.** A pipe, a zero-length file, a
 * platform without support, WASM, or any failure in the mapping call itself
 * silently falls back to a buffered read. Constructing a `FileSource` on a file
 * that exists always succeeds; only a genuinely unreadable file throws.
 *
 * **Lifetime is deliberately trivial.** Every reader copies what it parses into
 * owning `NDArray`/`std::string` storage via the uniform mesh API, so nothing
 * in a returned `Mesh` ever points into this buffer. `FileSource` is therefore
 * a function-local RAII object destroyed when the reader returns -- there is no
 * keep-alive to plumb through, no capsule, and no ownership shared with the
 * mesh. The one rule to police: **nothing derived from `View()` may be stored
 * in the returned mesh**. A gtest reads a file, then deletes and overwrites it,
 * and re-validates the mesh, so a violation fails loudly rather than becoming a
 * use-after-free in someone else's process.
 */

// System includes
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <string_view>
#include <utility>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/read_options.hpp"

#if defined(__EMSCRIPTEN__)
// No mapping under Emscripten: the virtual FS has nothing to map.
#elif defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
// This header reaches user translation units through the amalgamation, so
// <windows.h> must not be allowed to define min/max macros there.
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

/**
 * @brief Read-only whole-file access: mapped when worthwhile, copied otherwise.
 *
 * Move-only; the mapping (or buffer) is released in the destructor.
 */
class FileSource {
public:
    /** @brief Mapping preference. `Auto` maps only when it is likely to pay. */
    enum class Mode { Auto, Mmap, Buffered };

    /**
     * @brief Open @p rPath for whole-file reading.
     * @param rPath file to read.
     * @param mode mapping preference; `Auto` maps regular files at or above
     *        `mmap_auto_threshold_bytes`.
     * @throws ReadError only if the file cannot be read at all.
     */
    explicit FileSource(const std::string& rPath, Mode mode = Mode::Auto) {
        if (mode != Mode::Buffered && TryMap(rPath, mode == Mode::Mmap))
            return;
        LoadBuffered(rPath);
    }

    /** @brief Open honouring a `ReadOptions`' mmap preference. */
    FileSource(const std::string& rPath, MmapMode mmap_mode)
        : FileSource(rPath, FromMmapMode(mmap_mode)) {}

    ~FileSource() { Release(); }

    FileSource(const FileSource&) = delete;
    FileSource& operator=(const FileSource&) = delete;

    FileSource(FileSource&& rOther) noexcept { MoveFrom(std::move(rOther)); }

    FileSource& operator=(FileSource&& rOther) noexcept {
        if (this != &rOther) {
            Release();
            MoveFrom(std::move(rOther));
        }
        return *this;
    }

    /** @brief First byte of the file contents; never null (may be empty). */
    const char* Data() const { return mpData ? mpData : mBuffer.data(); }

    /** @brief File size in bytes. */
    std::size_t Size() const { return mSize; }

    /** @brief The whole file as a view. Valid until this object is destroyed. */
    std::string_view View() const { return std::string_view(Data(), mSize); }

    /** @brief Whether the contents are mapped rather than copied. */
    bool IsMapped() const { return mMapped; }

private:
    static Mode FromMmapMode(MmapMode mmap_mode) {
        switch (mmap_mode) {
            case MmapMode::On:
                return Mode::Mmap;
            case MmapMode::Off:
                return Mode::Buffered;
            default:
                return Mode::Auto;
        }
    }

    /** @brief Read the whole file into `mBuffer` with one bulk read. */
    void LoadBuffered(const std::string& rPath) {
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

    /**
     * @brief Attempt to map @p rPath.
     * @param force map regardless of size (still only for regular files).
     * @return false when mapping is unavailable or not worth it -- never throws,
     *         so the caller simply falls back.
     */
    bool TryMap(const std::string& rPath, bool force) {
#if defined(__EMSCRIPTEN__)
        (void)rPath;
        (void)force;
        return false;  // nothing to map on the virtual filesystem
#elif defined(_WIN32)
        HANDLE file = CreateFileA(rPath.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                                  OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (file == INVALID_HANDLE_VALUE)
            return false;
        LARGE_INTEGER size;
        if (!GetFileSizeEx(file, &size) || size.QuadPart <= 0 ||
            !WorthMapping(size.QuadPart, force)) {
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

    /**
     * @brief Whether a mapping of @p size bytes would have a readable zero byte
     *        just past the data.
     *
     * The kernel zero-fills the remainder of the final page, so a mapped buffer
     * normally behaves as if NUL-terminated -- which the C string functions some
     * parsers use (`strtod`) quietly rely on. That slack does **not** exist when
     * the file size is an exact multiple of the page size, and a file ending in
     * digits with no trailing delimiter would then read past the mapping and
     * fault. Refusing to map those (a ~1-in-4096 case) buys the guarantee back
     * for the cost of a buffered read.
     */
    static bool HasTerminatorSlack(long long size) {
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

    /**
     * @brief Whether a file of @p size bytes is worth mapping under `Auto`.
     *
     * Below the threshold the copy is cheap and the mapping's setup plus page
     * faults are not obviously better, so small files keep the simple path.
     * `MESHIOPLUSPLUS_MMAP_THRESHOLD` overrides it (bytes) for benchmarking; 0
     * means always map.
     */
    static bool WorthMapping(long long size, bool force) {
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

    void Release() {
#if defined(__EMSCRIPTEN__)
#elif defined(_WIN32)
        if (mMapped && mpData)
            UnmapViewOfFile(mpData);
        if (mWinMapping)
            CloseHandle(mWinMapping);
        if (mWinFile && mWinFile != INVALID_HANDLE_VALUE)
            CloseHandle(mWinFile);
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

    void MoveFrom(FileSource&& rOther) noexcept {
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

    std::string mBuffer;           ///< Backing store when not mapped.
    const char* mpData = nullptr;  ///< Mapped address, or null when buffered.
    std::size_t mSize = 0;
    bool mMapped = false;
#if defined(_WIN32) && !defined(__EMSCRIPTEN__)
    HANDLE mWinFile = nullptr;
    HANDLE mWinMapping = nullptr;
#endif
};

}  // namespace detail
}  // namespace meshioplusplus
