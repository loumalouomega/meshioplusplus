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
// Compressed input files: the codecs, the single-member zip reader and the
// temporary decompressed copy. See compressed_input.hpp for the contract.

// System includes
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <random>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

// External includes
#ifdef MESHIOPLUSPLUS_HAS_BZIP2
#include <bzlib.h>
#endif
#ifdef MESHIOPLUSPLUS_HAS_LZMA
#include <lzma.h>
#endif
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
#include <zstd.h>
#endif

// Project includes
#include "meshioplusplus/detail/compressed_input.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

std::string cmpin_lower(std::string Text) {
    std::transform(Text.begin(), Text.end(), Text.begin(), [](unsigned char c) {
        return static_cast<char>(c >= 'A' && c <= 'Z' ? c + 32 : c);
    });
    return Text;
}

[[noreturn]] void cmpin_missing(Compression Kind) {
    throw ReadError(std::string("meshio++: reading a ") + compression_suffix(Kind) +
                    " file needs a build with " + compression_option(Kind));
}

std::string cmpin_gunzip(std::string_view In) {
    // Every member, one after the other, as `gzip -d` reads them.
    std::string out;
    while (!In.empty()) {
        std::size_t used = 0;
        out += zlib_inflate(In, 31, &used, "gzip");
        if (used == 0)
            throw ReadError("meshio++: gzip: no progress decompressing the stream");
        In.remove_prefix(used);
        // Trailing zero padding (tar-style) after the last member is not a member.
        if (std::all_of(In.begin(), In.end(), [](char c) { return c == '\0'; }))
            break;
    }
    return out;
}

std::string cmpin_bunzip2(std::string_view In) {
#ifdef MESHIOPLUSPLUS_HAS_BZIP2
    std::string out;
    std::size_t pos = 0;
    std::vector<char> buf(1 << 16);
    while (pos < In.size()) {
        bz_stream z{};
        if (BZ2_bzDecompressInit(&z, 0, 0) != BZ_OK)
            throw ReadError("meshio++: bzip2: initialisation failed");
        z.next_in = const_cast<char*>(In.data() + pos);
        z.avail_in = static_cast<unsigned int>(In.size() - pos);
        int rc = BZ_OK;
        while (rc == BZ_OK) {
            z.next_out = buf.data();
            z.avail_out = static_cast<unsigned int>(buf.size());
            rc = BZ2_bzDecompress(&z);
            out.append(buf.data(), buf.size() - z.avail_out);
            if (rc == BZ_OK && z.avail_in == 0 && z.avail_out != 0)
                rc = BZ_UNEXPECTED_EOF;
        }
        pos = In.size() - z.avail_in;
        BZ2_bzDecompressEnd(&z);
        if (rc != BZ_STREAM_END)
            throw ReadError("meshio++: bzip2: corrupt or truncated stream");
    }
    return out;
#else
    (void)In;
    cmpin_missing(Compression::Bzip2);
#endif
}

std::string cmpin_unxz(std::string_view In) {
#ifdef MESHIOPLUSPLUS_HAS_LZMA
    lzma_stream z = LZMA_STREAM_INIT;
    if (lzma_stream_decoder(&z, UINT64_MAX, LZMA_CONCATENATED) != LZMA_OK)
        throw ReadError("meshio++: xz: initialisation failed");
    z.next_in = reinterpret_cast<const std::uint8_t*>(In.data());
    z.avail_in = In.size();
    std::string out;
    std::vector<std::uint8_t> buf(1 << 16);
    lzma_ret rc = LZMA_OK;
    while (rc == LZMA_OK) {
        z.next_out = buf.data();
        z.avail_out = buf.size();
        rc = lzma_code(&z, z.avail_in == 0 ? LZMA_FINISH : LZMA_RUN);
        out.append(reinterpret_cast<const char*>(buf.data()), buf.size() - z.avail_out);
    }
    lzma_end(&z);
    if (rc != LZMA_STREAM_END)
        throw ReadError("meshio++: xz: corrupt or truncated stream");
    return out;
#else
    (void)In;
    cmpin_missing(Compression::Xz);
#endif
}

std::string cmpin_unzstd(std::string_view In) {
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
    ZSTD_DStream* z = ZSTD_createDStream();
    if (!z)
        throw ReadError("meshio++: zstd: initialisation failed");
    ZSTD_initDStream(z);
    ZSTD_inBuffer in{In.data(), In.size(), 0};
    std::string out;
    std::vector<char> buf(ZSTD_DStreamOutSize());
    std::size_t rc = 0;
    while (in.pos < in.size) {
        ZSTD_outBuffer o{buf.data(), buf.size(), 0};
        rc = ZSTD_decompressStream(z, &o, &in);
        if (ZSTD_isError(rc)) {
            ZSTD_freeDStream(z);
            throw ReadError(std::string("meshio++: zstd: ") + ZSTD_getErrorName(rc));
        }
        out.append(buf.data(), o.pos);
    }
    // Flush what the last frame still holds.
    while (rc != 0) {
        ZSTD_outBuffer o{buf.data(), buf.size(), 0};
        rc = ZSTD_decompressStream(z, &o, &in);
        if (ZSTD_isError(rc) || o.pos == 0) {
            ZSTD_freeDStream(z);
            throw ReadError("meshio++: zstd: corrupt or truncated frame");
        }
        out.append(buf.data(), o.pos);
    }
    ZSTD_freeDStream(z);
    return out;
#else
    (void)In;
    cmpin_missing(Compression::Zstd);
#endif
}

// --- zip: the one member of a single-file archive ---------------------------

std::uint32_t cmpin_u16(std::string_view In, std::size_t At) {
    if (At + 2 > In.size())
        throw ReadError("meshio++: zip: truncated archive");
    return static_cast<std::uint32_t>(static_cast<unsigned char>(In[At])) |
           static_cast<std::uint32_t>(static_cast<unsigned char>(In[At + 1])) << 8;
}

std::uint32_t cmpin_u32(std::string_view In, std::size_t At) {
    return cmpin_u16(In, At) | cmpin_u16(In, At + 2) << 16;
}

struct ZipEntry {
    std::string mName;
    std::uint32_t mMethod = 0;
    std::uint32_t mCompressedSize = 0;
    std::uint32_t mLocalOffset = 0;
};

// The archive's one file entry (directory entries are skipped).
ZipEntry cmpin_zip_entry(std::string_view In) {
    // The end-of-central-directory record: within the last 64 KiB + 22 bytes.
    if (In.size() < 22)
        throw ReadError("meshio++: zip: not a zip archive");
    std::size_t eocd = std::string_view::npos;
    const std::size_t lowest = In.size() > 65557 ? In.size() - 65557 : 0;
    for (std::size_t at = In.size() - 22 + 1; at-- > lowest;) {
        if (cmpin_u32(In, at) == 0x06054b50u) {
            eocd = at;
            break;
        }
    }
    if (eocd == std::string_view::npos)
        throw ReadError("meshio++: zip: no end-of-central-directory record");
    const std::uint32_t count = cmpin_u16(In, eocd + 10);
    std::size_t at = cmpin_u32(In, eocd + 16);
    if (count == 0xFFFFu || at == 0xFFFFFFFFu)
        throw ReadError("meshio++: zip: Zip64 archives are not read");
    std::vector<ZipEntry> files;
    for (std::uint32_t k = 0; k < count; ++k) {
        if (cmpin_u32(In, at) != 0x02014b50u)
            throw ReadError("meshio++: zip: corrupt central directory");
        ZipEntry e;
        e.mMethod = cmpin_u16(In, at + 10);
        e.mCompressedSize = cmpin_u32(In, at + 20);
        const std::uint32_t name_len = cmpin_u16(In, at + 28);
        const std::uint32_t extra_len = cmpin_u16(In, at + 30);
        const std::uint32_t comment_len = cmpin_u16(In, at + 32);
        e.mLocalOffset = cmpin_u32(In, at + 42);
        if (at + 46 + name_len > In.size())
            throw ReadError("meshio++: zip: truncated central directory");
        e.mName.assign(In.substr(at + 46, name_len));
        at += 46 + name_len + extra_len + comment_len;
        if (!e.mName.empty() && e.mName.back() == '/')
            continue;  // a directory
        files.push_back(std::move(e));
    }
    if (files.size() != 1)
        throw ReadError("meshio++: zip: the archive holds " + std::to_string(files.size()) +
                        " files; only a single-file archive is read as a mesh");
    if (files[0].mCompressedSize == 0xFFFFFFFFu || files[0].mLocalOffset == 0xFFFFFFFFu)
        throw ReadError("meshio++: zip: Zip64 archives are not read");
    return files[0];
}

std::string cmpin_unzip(std::string_view In) {
    const ZipEntry e = cmpin_zip_entry(In);
    const std::size_t local = e.mLocalOffset;
    if (cmpin_u32(In, local) != 0x04034b50u)
        throw ReadError("meshio++: zip: corrupt local header");
    if (cmpin_u16(In, local + 6) & 1u)
        throw ReadError("meshio++: zip: encrypted members are not read");
    const std::size_t data = local + 30 + cmpin_u16(In, local + 26) + cmpin_u16(In, local + 28);
    if (data + e.mCompressedSize > In.size())
        throw ReadError("meshio++: zip: truncated member data");
    const std::string_view body = In.substr(data, e.mCompressedSize);
    if (e.mMethod == 0)
        return std::string(body);
    if (e.mMethod == 8) {
        if (!zlib_available())
            cmpin_missing(Compression::Zip);
        return zlib_inflate(body, -15, nullptr, "zip");
    }
    throw ReadError("meshio++: zip: compression method " + std::to_string(e.mMethod) +
                    " is not read (only stored and deflated members)");
}

std::string cmpin_slurp(const std::string& rPath) {
    auto in = make_classic_ifstream(std::filesystem::path(rPath), std::ios::binary);
    if (!in)
        throw ReadError("meshio++: cannot open '" + rPath + "'");
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string cmpin_basename(const std::string& rPath) {
    const std::size_t slash = rPath.find_last_of("/\\");
    return slash == std::string::npos ? rPath : rPath.substr(slash + 1);
}

}  // namespace

const char* compression_suffix(Compression Kind) {
    switch (Kind) {
        case Compression::Gzip:
            return ".gz";
        case Compression::Bzip2:
            return ".bz2";
        case Compression::Xz:
            return ".xz";
        case Compression::Zstd:
            return ".zst";
        case Compression::Zip:
            return ".zip";
        case Compression::None:
            break;
    }
    return "";
}

Compression compression_from_suffix(const std::string& rPath) {
    const std::string base = cmpin_lower(cmpin_basename(rPath));
    const std::size_t dot = base.rfind('.');
    if (dot == std::string::npos || dot == 0)
        return Compression::None;
    const std::string ext = base.substr(dot);
    for (const Compression kind : {Compression::Gzip, Compression::Bzip2, Compression::Xz,
                                   Compression::Zstd, Compression::Zip})
        if (ext == compression_suffix(kind))
            return kind;
    return Compression::None;
}

bool compression_available(Compression Kind) {
    switch (Kind) {
        case Compression::Gzip:
        case Compression::Zip:
            return zlib_available();
        case Compression::Bzip2:
#ifdef MESHIOPLUSPLUS_HAS_BZIP2
            return true;
#else
            return false;
#endif
        case Compression::Xz:
#ifdef MESHIOPLUSPLUS_HAS_LZMA
            return true;
#else
            return false;
#endif
        case Compression::Zstd:
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
            return true;
#else
            return false;
#endif
        case Compression::None:
            return true;
    }
    return false;
}

const char* compression_option(Compression Kind) {
    switch (Kind) {
        case Compression::Gzip:
        case Compression::Zip:
            return "MESHIOPLUSPLUS_WITH_ZLIB";
        case Compression::Bzip2:
            return "MESHIOPLUSPLUS_WITH_BZIP2";
        case Compression::Xz:
            return "MESHIOPLUSPLUS_WITH_LZMA";
        case Compression::Zstd:
            return "MESHIOPLUSPLUS_WITH_ZSTD";
        case Compression::None:
            break;
    }
    return "";
}

std::string decompress_bytes(Compression Kind, std::string_view In) {
    switch (Kind) {
        case Compression::Gzip:
            if (!zlib_available())
                cmpin_missing(Kind);
            return cmpin_gunzip(In);
        case Compression::Bzip2:
            return cmpin_bunzip2(In);
        case Compression::Xz:
            return cmpin_unxz(In);
        case Compression::Zstd:
            return cmpin_unzstd(In);
        case Compression::Zip:
            return cmpin_unzip(In);
        case Compression::None:
            break;
    }
    return std::string(In);
}

std::string zip_member_name(std::string_view Archive) {
    return cmpin_zip_entry(Archive).mName;
}

std::string decompressed_name(const std::string& rPath, Compression Kind) {
    const std::string suffix = compression_suffix(Kind);
    std::string stripped = rPath.substr(0, rPath.size() - std::min(rPath.size(), suffix.size()));
    if (Kind == Compression::Zip && cmpin_basename(stripped).find('.') == std::string::npos) {
        // `case.zip`: the member's own name says what it is.
        const std::string member = cmpin_basename(zip_member_name(cmpin_slurp(rPath)));
        const std::size_t slash = stripped.find_last_of("/\\");
        stripped =
            (slash == std::string::npos ? std::string() : stripped.substr(0, slash + 1)) + member;
    }
    return stripped;
}

DecompressedFile::DecompressedFile(const std::string& rPath, Compression Kind) {
    if (Kind == Compression::None)
        throw ReadError("meshio++: '" + rPath + "' is not compressed");
    const std::string bytes = decompress_bytes(Kind, cmpin_slurp(rPath));
    std::string name = cmpin_basename(decompressed_name(rPath, Kind));
    if (name.empty() || name == "." || name == "..")
        name = "decompressed";

    // A fresh private directory: the decompressed file keeps its inner name,
    // which some formats are recognised by (z88i1.txt, d3plot).
    static std::atomic<unsigned> counter{0};
    std::random_device rd;
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    std::error_code ec;
    const std::filesystem::path tmp = std::filesystem::temp_directory_path(ec);
    if (ec)
        throw ReadError("meshio++: no temporary directory to decompress '" + rPath + "' into");
    for (int attempt = 0; attempt < 16; ++attempt) {
        const std::filesystem::path dir =
            tmp / ("meshioplusplus-" + std::to_string(rd()) + "-" + std::to_string(stamp) + "-" +
                   std::to_string(counter++));
        if (std::filesystem::create_directory(dir, ec)) {
            mDir = dir.string();
            break;
        }
    }
    if (mDir.empty())
        throw ReadError("meshio++: could not create a temporary directory to decompress '" + rPath +
                        "' into");
    mPath = (std::filesystem::path(mDir) / name).string();
    auto out = make_classic_ofstream(std::filesystem::path(mPath), std::ios::binary);
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    if (!out) {
        out.close();
        std::filesystem::remove_all(mDir, ec);
        throw ReadError("meshio++: could not write the decompressed copy of '" + rPath + "'");
    }
}

DecompressedFile::~DecompressedFile() {
    std::error_code ec;
    if (!mDir.empty())
        std::filesystem::remove_all(mDir, ec);
}

}  // namespace detail
}  // namespace meshioplusplus
