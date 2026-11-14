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
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <stdexcept>
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
#include "subprocess.hpp"
#include "zip_archive.hpp"
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

// The archive's one file member, through the shared ZIP reader (detail/zip_archive.hpp).
std::string cmpin_unzip(std::string_view In) {
    const ZipReader zip(In);
    const auto& files = zip.Members();
    if (files.size() != 1)
        throw ReadError("meshio++: zip: the archive holds " + std::to_string(files.size()) +
                        " files; only a single-file archive is read as a mesh");
    if (files[0].mMethod == 8 && !zlib_available())
        cmpin_missing(Compression::Zip);
    return zip.Read(files[0]);
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
    const ZipReader zip(Archive);
    if (zip.Members().size() != 1)
        throw ReadError("meshio++: zip: the archive holds " + std::to_string(zip.Members().size()) +
                        " files; only a single-file archive is read as a mesh");
    return zip.Members().front().mName;
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

    // A fresh private directory (make_temp_dir: created exclusively, owner
    // only): the decompressed file keeps its inner name, which some formats
    // are recognised by (z88i1.txt, d3plot).
    try {
        mDir = make_temp_dir("decompress");
    } catch (const std::runtime_error&) {
        throw ReadError("meshio++: could not create a temporary directory to decompress '" + rPath +
                        "' into");
    }
    std::error_code ec;
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
