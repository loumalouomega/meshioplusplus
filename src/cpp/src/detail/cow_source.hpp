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
 * @file detail/cow_source.hpp
 * @brief A whole file as **writable** bytes whose writes never reach the
 * file: a copy-on-write mapping when worthwhile, an owned buffer otherwise.
 *
 * A **core-private** header (the `open_source.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI. It exists for
 * parsers that write into their input -- pugixml's `load_buffer_inplace`
 * terminates and unescapes text where it lies -- so that a mapped file is not
 * copied whole first (roadmap 5.7.4): only the pages the parser writes are
 * copied, by the kernel. `FileSource` stays read-only.
 *
 * Maps on POSIX (`MAP_PRIVATE` with `PROT_WRITE`) and Windows
 * (`FILE_MAP_COPY`) at or above `mmap_auto_threshold_bytes` (or
 * `MESHIOPLUSPLUS_MMAP_THRESHOLD`), as `FileSource`'s `Auto` mode does; reads
 * smaller files, and every file under Emscripten, into an owned buffer.
 */

// System includes
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace meshioplusplus {
namespace detail {

class CowSource {
public:
    /// Opens @p rPath; `ReadError(rMessage)` when it cannot be read. Maps it
    /// at or above @p MapAtOrAbove bytes (0: the `Auto` threshold).
    CowSource(const std::string& rPath, const std::string& rMessage, std::size_t MapAtOrAbove = 0);
    ~CowSource();

    CowSource(const CowSource&) = delete;
    CowSource& operator=(const CowSource&) = delete;
    CowSource(CowSource&& rOther) noexcept;
    CowSource& operator=(CowSource&& rOther) noexcept;

    /// The file's bytes, writable; never null (may be empty).
    char* Data() { return mpData ? mpData : mBuffer.data(); }
    const char* Data() const { return mpData ? mpData : mBuffer.data(); }
    std::size_t Size() const { return mSize; }
    std::string_view View() const { return {Data(), mSize}; }
    /// Whether the bytes are a mapping (not an owned copy).
    bool IsMapped() const { return mpData != nullptr; }

private:
    bool TryMap(const std::string& rPath, std::size_t MapAtOrAbove);
    void Release() noexcept;
    void MoveFrom(CowSource&& rOther) noexcept;

    /// The bytes read, plus a NUL, when not mapped. A vector, so a move keeps
    /// the bytes where a parse in place left pointers into them.
    std::vector<char> mBuffer = std::vector<char>(1, '\0');
    char* mpData = nullptr;  ///< the mapping, or null for the buffer
    std::size_t mSize = 0;
    void* mpWinMapping = nullptr;
};

}  // namespace detail
}  // namespace meshioplusplus
