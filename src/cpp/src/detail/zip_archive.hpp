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
 * @file detail/zip_archive.hpp
 * @brief A ZIP container, core-private: list the members, read one by name, and
 * write stored or deflated members with a central directory (roadmap 2.4.5.1).
 *
 * Read: the end-of-central-directory record is found in the last 64 KiB + 22
 * bytes, every central-directory entry is listed (directories skipped), and a
 * member is read through its local header, stored or deflated (raw deflate through
 * `zlib_inflate`). Zip64 and encrypted members are refused by name.
 *
 * Write: each member is deflated at level 9 (the raw stream of `gzip_compress`)
 * when the build has zlib and that is smaller, else stored; DOS time and date are
 * 1980-01-01 00:00, so the bytes depend on the contents alone. The Python twin is
 * `meshioplusplus._zip`, which writes the same bytes.
 */

// System includes
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace meshioplusplus::detail {

/** @brief One file member of an archive. */
struct ZipMember {
    std::string mName;
    std::uint16_t mMethod = 0;  ///< 0 stored, 8 deflated
    std::uint32_t mCompressedSize = 0;
    std::uint32_t mSize = 0;
    std::uint32_t mLocalOffset = 0;
};

/** @brief The members of the archive in @p Data (its bytes outlive this view). */
class ZipReader {
public:
    /// @throws ReadError when @p Data is not a (non-Zip64) ZIP archive.
    explicit ZipReader(std::string_view Data);

    const std::vector<ZipMember>& Members() const { return mMembers; }

    /// The member named @p rName, or null.
    const ZipMember* Find(const std::string& rName) const;

    /// The bytes of member @p rMember. @throws ReadError on a corrupt or encrypted member.
    std::string Read(const ZipMember& rMember) const;

    /// The bytes of the member named @p rName. @throws ReadError when there is none.
    std::string Read(const std::string& rName) const;

private:
    std::string_view mData;
    std::vector<ZipMember> mMembers;
};

/** @brief Builds an archive member by member. */
class ZipWriter {
public:
    /// Add a member; @p Deflate asks for compression (used only when it helps).
    void Add(const std::string& rName, std::string_view Data, bool Deflate = true);

    /// The finished archive. @throws WriteError above 4 GiB (no Zip64).
    std::string Finish();

private:
    std::string mOut;
    std::string mDirectory;
    std::uint16_t mCount = 0;
};

}  // namespace meshioplusplus::detail
