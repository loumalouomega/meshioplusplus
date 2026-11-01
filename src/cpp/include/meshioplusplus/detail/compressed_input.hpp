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
 * @file detail/compressed_input.hpp
 * @brief Reading a compressed mesh file as if it were not: `cube.msh.gz`,
 * `part.vtu.bz2`, `run.inp.xz`, `skin.stl.zst` and a single-member `case.zip`.
 *
 * The registry strips the compression suffix to infer the format
 * (`resolve_format`), and a read of a compressed path decompresses it into a
 * private temporary directory under the name it had inside (`cube.msh`), reads
 * that, and removes the directory (`DecompressedFile`). A compound extension a
 * format registers itself -- Netgen's `.vol.gz`, libMesh's `.xda.gz`, PERMAS's
 * `.dato.gz` -- is that format's own business and is never stripped.
 *
 * The codecs are the ones the core already links: gzip and zip through zlib
 * (`MESHIOPLUSPLUS_WITH_ZLIB`), bzip2 (`MESHIOPLUSPLUS_WITH_BZIP2`), xz
 * (`MESHIOPLUSPLUS_WITH_LZMA`) and zstd (`MESHIOPLUSPLUS_WITH_ZSTD`). A codec
 * the build lacks is a `ReadError` naming its option, never a misread.
 */

// System includes
#include <cstdint>
#include <string>
#include <string_view>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {
namespace detail {

/// A whole-file compression wrapper.
enum class Compression : std::int32_t {
    None = 0,   ///< Not compressed (or a compound extension a format owns).
    Gzip = 1,   ///< `.gz`, gzip (RFC 1952), members concatenated.
    Bzip2 = 2,  ///< `.bz2`, streams concatenated.
    Xz = 3,     ///< `.xz`, streams concatenated.
    Zstd = 4,   ///< `.zst`, frames concatenated.
    Zip = 5,    ///< `.zip` holding exactly one file (stored or deflated).
};

/** @brief The suffix of @p Kind (`".gz"`...), empty for `None`. */
MESHIOPLUSPLUS_API const char* compression_suffix(Compression Kind);

/** @brief The codec named by @p rPath's last suffix (case-insensitive); `None` otherwise. */
MESHIOPLUSPLUS_API Compression compression_from_suffix(const std::string& rPath);

/** @brief Whether this build can decompress @p Kind. */
MESHIOPLUSPLUS_API bool compression_available(Compression Kind);

/** @brief The CMake option that enables @p Kind (`"MESHIOPLUSPLUS_WITH_ZLIB"`...). */
MESHIOPLUSPLUS_API const char* compression_option(Compression Kind);

/**
 * @brief The bytes @p In decompress to under @p Kind; for `Zip`, the one
 * member's bytes.
 * @throws ReadError on corrupt or truncated input, a zip that does not hold
 *         exactly one file (or uses an unsupported method or Zip64), or a codec
 *         this build lacks (naming its option).
 */
MESHIOPLUSPLUS_API std::string decompress_bytes(Compression Kind, std::string_view In);

/** @brief The name of the one file a zip archive holds (its path inside). */
MESHIOPLUSPLUS_API std::string zip_member_name(std::string_view Archive);

/**
 * @brief The name @p rPath's content had before compression: the path without
 * the compression suffix, or for a zip whose name carries no inner extension
 * (`case.zip`) the member's own basename beside it.
 */
MESHIOPLUSPLUS_API std::string decompressed_name(const std::string& rPath, Compression Kind);

/**
 * @brief A compressed file decompressed into a private temporary directory,
 * removed again when this goes out of scope.
 */
class MESHIOPLUSPLUS_API DecompressedFile {
public:
    /**
     * @param rPath The compressed file.
     * @param Kind Its codec (not `None`).
     * @throws ReadError when it cannot be read or decompressed.
     */
    DecompressedFile(const std::string& rPath, Compression Kind);
    ~DecompressedFile();
    DecompressedFile(const DecompressedFile&) = delete;
    DecompressedFile& operator=(const DecompressedFile&) = delete;

    /** @brief The decompressed file, named as it was before compression. */
    const std::string& Path() const { return mPath; }

private:
    std::string mDir;
    std::string mPath;
};

}  // namespace detail
}  // namespace meshioplusplus
