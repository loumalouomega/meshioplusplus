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
// A ZIP container: list, read and write members. See detail/zip_archive.hpp.

// System includes
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "zip_archive.hpp"
#include "gzip_deflate.hpp"
#include "png_write.hpp"
#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus::detail {

namespace {

std::uint32_t zip_u16(std::string_view In, std::size_t At) {
    if (At + 2 > In.size())
        throw ReadError("meshio++: zip: truncated archive");
    return static_cast<unsigned char>(In[At]) |
           static_cast<std::uint32_t>(static_cast<unsigned char>(In[At + 1])) << 8;
}

std::uint32_t zip_u32(std::string_view In, std::size_t At) {
    return zip_u16(In, At) | zip_u16(In, At + 2) << 16;
}

void zip_put16(std::string& rOut, std::uint32_t V) {
    rOut += static_cast<char>(V & 0xff);
    rOut += static_cast<char>((V >> 8) & 0xff);
}

void zip_put32(std::string& rOut, std::uint32_t V) {
    zip_put16(rOut, V & 0xffff);
    zip_put16(rOut, V >> 16);
}

}  // namespace

ZipReader::ZipReader(std::string_view Data) : mData(Data) {
    if (Data.size() < 22)
        throw ReadError("meshio++: zip: not a zip archive");
    std::size_t eocd = std::string_view::npos;
    const std::size_t lowest = Data.size() > 65557 ? Data.size() - 65557 : 0;
    for (std::size_t at = Data.size() - 22 + 1; at-- > lowest;) {
        if (zip_u32(Data, at) == 0x06054b50u) {
            eocd = at;
            break;
        }
    }
    if (eocd == std::string_view::npos)
        throw ReadError("meshio++: zip: no end-of-central-directory record");
    const std::uint32_t count = zip_u16(Data, eocd + 10);
    std::size_t at = zip_u32(Data, eocd + 16);
    if (count == 0xFFFFu || at == 0xFFFFFFFFu)
        throw ReadError("meshio++: zip: Zip64 archives are not read");
    for (std::uint32_t k = 0; k < count; ++k) {
        if (zip_u32(Data, at) != 0x02014b50u)
            throw ReadError("meshio++: zip: corrupt central directory");
        ZipMember m;
        m.mMethod = static_cast<std::uint16_t>(zip_u16(Data, at + 10));
        m.mCompressedSize = zip_u32(Data, at + 20);
        m.mSize = zip_u32(Data, at + 24);
        const std::uint32_t name_len = zip_u16(Data, at + 28);
        const std::uint32_t extra_len = zip_u16(Data, at + 30);
        const std::uint32_t comment_len = zip_u16(Data, at + 32);
        m.mLocalOffset = zip_u32(Data, at + 42);
        if (at + 46 + name_len > Data.size())
            throw ReadError("meshio++: zip: truncated central directory");
        m.mName.assign(Data.substr(at + 46, name_len));
        at += 46 + name_len + extra_len + comment_len;
        if (!m.mName.empty() && m.mName.back() == '/')
            continue;  // a directory
        if (m.mCompressedSize == 0xFFFFFFFFu || m.mLocalOffset == 0xFFFFFFFFu ||
            m.mSize == 0xFFFFFFFFu)
            throw ReadError("meshio++: zip: Zip64 archives are not read");
        mMembers.push_back(std::move(m));
    }
}

const ZipMember* ZipReader::Find(const std::string& rName) const {
    for (const ZipMember& m : mMembers)
        if (m.mName == rName)
            return &m;
    return nullptr;
}

std::string ZipReader::Read(const ZipMember& rMember) const {
    const std::size_t local = rMember.mLocalOffset;
    if (zip_u32(mData, local) != 0x04034b50u)
        throw ReadError("meshio++: zip: corrupt local header");
    if (zip_u16(mData, local + 6) & 1u)
        throw ReadError("meshio++: zip: encrypted members are not read");
    const std::size_t data = local + 30 + zip_u16(mData, local + 26) + zip_u16(mData, local + 28);
    if (data + rMember.mCompressedSize > mData.size())
        throw ReadError("meshio++: zip: truncated member data");
    const std::string_view body = mData.substr(data, rMember.mCompressedSize);
    if (rMember.mMethod == 0)
        return std::string(body);
    if (rMember.mMethod == 8)
        return zlib_inflate(body, -15, nullptr, "zip");
    throw ReadError("meshio++: zip: compression method " + std::to_string(rMember.mMethod) +
                    " is not read (only stored and deflated members)");
}

std::string ZipReader::Read(const std::string& rName) const {
    const ZipMember* m = Find(rName);
    if (!m)
        throw ReadError("meshio++: zip: the archive has no member '" + rName + "'");
    return Read(*m);
}

void ZipWriter::Add(const std::string& rName, std::string_view Data, bool Deflate) {
    std::string body(Data);
    std::uint16_t method = 0;
    if (Deflate && gzip_compress_available() && !Data.empty()) {
        // gzip_compress is a 10-byte header, the raw deflate stream and 8 bytes of trailer.
        const std::string gz = gzip_compress(Data, "zip");
        const std::string raw = gz.substr(10, gz.size() - 18);
        if (raw.size() < Data.size()) {
            body = raw;
            method = 8;
        }
    }
    const std::uint32_t crc =
        png_crc32(0, reinterpret_cast<const unsigned char*>(Data.data()), Data.size());
    const std::uint64_t offset = mOut.size();
    if (offset + body.size() > 0xFFFFFFFFull || Data.size() > 0xFFFFFFFFull || mCount == 0xFFFF)
        throw WriteError("meshio++: zip: the archive would need Zip64, which is not written");
    const auto header = [&](std::string& rTo, bool Central) {
        zip_put32(rTo, Central ? 0x02014b50u : 0x04034b50u);
        if (Central)
            zip_put16(rTo, 20);  // version made by
        zip_put16(rTo, 20);      // version needed
        zip_put16(rTo, 0);       // flags
        zip_put16(rTo, method);
        zip_put16(rTo, 0);     // time 00:00:00
        zip_put16(rTo, 0x21);  // date 1980-01-01
        zip_put32(rTo, crc);
        zip_put32(rTo, static_cast<std::uint32_t>(body.size()));
        zip_put32(rTo, static_cast<std::uint32_t>(Data.size()));
        zip_put16(rTo, static_cast<std::uint32_t>(rName.size()));
        zip_put16(rTo, 0);  // extra
        if (Central) {
            zip_put16(rTo, 0);  // comment
            zip_put16(rTo, 0);  // disk
            zip_put16(rTo, 0);  // internal attributes
            zip_put32(rTo, 0);  // external attributes
            zip_put32(rTo, static_cast<std::uint32_t>(offset));
        }
        rTo += rName;
    };
    header(mOut, false);
    mOut += body;
    header(mDirectory, true);
    ++mCount;
}

std::string ZipWriter::Finish() {
    std::string out = mOut + mDirectory;
    if (out.size() > 0xFFFFFFFFull)
        throw WriteError("meshio++: zip: the archive would need Zip64, which is not written");
    zip_put32(out, 0x06054b50u);
    zip_put16(out, 0);
    zip_put16(out, 0);
    zip_put16(out, mCount);
    zip_put16(out, mCount);
    zip_put32(out, static_cast<std::uint32_t>(mDirectory.size()));
    zip_put32(out, static_cast<std::uint32_t>(mOut.size()));
    zip_put16(out, 0);
    return out;
}

}  // namespace meshioplusplus::detail
