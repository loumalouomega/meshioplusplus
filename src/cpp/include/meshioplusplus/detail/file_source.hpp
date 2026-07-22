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
 *
 * Implementation (`src/cpp/src/file_source.cpp`) is out-of-line: every platform
 * header (`<windows.h>`, `<sys/mman.h>`, ...) stays there rather than leaking
 * into every translation unit -- and every amalgamation consumer -- that merely
 * needs the declarations below.
 */

// System includes
#include <cstddef>
#include <string>
#include <string_view>

// Project includes
#include "meshioplusplus/read_options.hpp"

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
    explicit FileSource(const std::string& rPath, Mode mode = Mode::Auto);

    /** @brief Open honouring a `ReadOptions`' mmap preference. */
    FileSource(const std::string& rPath, MmapMode mmap_mode);

    ~FileSource();

    FileSource(const FileSource&) = delete;
    FileSource& operator=(const FileSource&) = delete;

    FileSource(FileSource&& rOther) noexcept;
    FileSource& operator=(FileSource&& rOther) noexcept;

    /** @brief First byte of the file contents; never null (may be empty). */
    const char* Data() const { return mpData ? mpData : mBuffer.data(); }

    /** @brief File size in bytes. */
    std::size_t Size() const { return mSize; }

    /** @brief The whole file as a view. Valid until this object is destroyed. */
    std::string_view View() const { return std::string_view(Data(), mSize); }

    /** @brief Whether the contents are mapped rather than copied. */
    bool IsMapped() const { return mMapped; }

private:
    static Mode FromMmapMode(MmapMode mmap_mode);

    /** @brief Read the whole file into `mBuffer` with one bulk read. */
    void LoadBuffered(const std::string& rPath);

    /**
     * @brief Attempt to map @p rPath.
     * @param force map regardless of size (still only for regular files).
     * @return false when mapping is unavailable or not worth it -- never throws,
     *         so the caller simply falls back.
     */
    bool TryMap(const std::string& rPath, bool force);

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
    static bool HasTerminatorSlack(long long size);

    /**
     * @brief Whether a file of @p size bytes is worth mapping under `Auto`.
     *
     * Below the threshold the copy is cheap and the mapping's setup plus page
     * faults are not obviously better, so small files keep the simple path.
     * `MESHIOPLUSPLUS_MMAP_THRESHOLD` overrides it (bytes) for benchmarking; 0
     * means always map.
     */
    static bool WorthMapping(long long size, bool force);

    void Release();
    void MoveFrom(FileSource&& rOther) noexcept;

    std::string mBuffer;           ///< Backing store when not mapped.
    const char* mpData = nullptr;  ///< Mapped address, or null when buffered.
    std::size_t mSize = 0;
    bool mMapped = false;
#if defined(_WIN32) && !defined(__EMSCRIPTEN__)
    // void*, not HANDLE: HANDLE is void* under the hood, and using it directly
    // here would drag <windows.h> into every consumer of this header. The
    // .cpp casts back to HANDLE where it actually calls the Win32 API.
    void* mWinFile = nullptr;
    void* mWinMapping = nullptr;
#endif
};

}  // namespace detail
}  // namespace meshioplusplus
