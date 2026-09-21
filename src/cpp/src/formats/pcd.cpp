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
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <ios>
#include <iterator>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/pcd.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr std::size_t kPcdMaxLiteral = 32;
constexpr std::size_t kPcdMaxRef = 264;  // 2 + 7 + 255
constexpr std::size_t kPcdMaxOffset = std::size_t(1) << 13;
constexpr unsigned kPcdHashBits = 14;
constexpr double kPcdIdentityViewpoint[7] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};

// ---------------------------------------------------------------------------
// LZF: the liblzf stream format (the Python twin is pcd/_lzf.py)

std::string pcd_lzf_decompress(const unsigned char* pSrc, std::size_t n, std::size_t out_size) {
    std::string out(out_size, '\0');
    std::size_t ip = 0, op = 0;
    while (ip < n) {
        const unsigned ctrl = pSrc[ip++];
        if (ctrl < 32) {
            const std::size_t length = ctrl + 1;
            if (ip + length > n)
                throw ReadError("PCD: truncated LZF literal run");
            if (op + length > out_size)
                throw ReadError("PCD: LZF stream overruns the declared size");
            std::memcpy(&out[op], pSrc + ip, length);
            ip += length;
            op += length;
        } else {
            std::size_t length = ctrl >> 5;
            if (length == 7) {
                if (ip >= n)
                    throw ReadError("PCD: truncated LZF back reference");
                length += pSrc[ip++];
            }
            if (ip >= n)
                throw ReadError("PCD: truncated LZF back reference");
            const std::size_t back = (static_cast<std::size_t>(ctrl & 0x1F) << 8) + pSrc[ip++] + 1;
            length += 2;
            if (back > op)
                throw ReadError("PCD: LZF back reference before the start of the data");
            if (op + length > out_size)
                throw ReadError("PCD: LZF stream overruns the declared size");
            const std::size_t ref = op - back;
            for (std::size_t k = 0; k < length; ++k)  // may overlap: the run repeats
                out[op + k] = out[ref + k];
            op += length;
        }
    }
    if (op != out_size)
        throw ReadError("PCD: LZF stream decoded to " + std::to_string(op) + " bytes, expected " +
                        std::to_string(out_size));
    return out;
}

std::string pcd_lzf_compress(const unsigned char* pSrc, std::size_t n) {
    std::string out;
    std::vector<std::int64_t> table(std::size_t(1) << kPcdHashBits, -1);
    std::size_t lit_start = 0, ip = 0;
    auto flush = [&](std::size_t end) {
        std::size_t pos = lit_start;
        while (pos < end) {
            const std::size_t run = std::min(kPcdMaxLiteral, end - pos);
            out.push_back(static_cast<char>(run - 1));
            out.append(reinterpret_cast<const char*>(pSrc) + pos, run);
            pos += run;
        }
    };
    while (ip + 2 < n) {
        std::uint64_t h = (static_cast<std::uint64_t>(pSrc[ip]) << 16) |
                          (static_cast<std::uint64_t>(pSrc[ip + 1]) << 8) | pSrc[ip + 2];
        h = ((h * 2654435761ULL) >> 8) & ((std::uint64_t(1) << kPcdHashBits) - 1);
        const std::int64_t ref = table[h];
        table[h] = static_cast<std::int64_t>(ip);
        if (ref >= 0 && ip - static_cast<std::size_t>(ref) <= kPcdMaxOffset &&
            pSrc[ref] == pSrc[ip] && pSrc[ref + 1] == pSrc[ip + 1] &&
            pSrc[ref + 2] == pSrc[ip + 2]) {
            std::size_t length = 3;
            const std::size_t limit = std::min(kPcdMaxRef, n - ip);
            while (length < limit && pSrc[ref + length] == pSrc[ip + length])
                ++length;
            flush(ip);
            const std::size_t off = ip - static_cast<std::size_t>(ref) - 1;
            const std::size_t coded = length - 2;
            if (coded < 7) {
                out.push_back(static_cast<char>((coded << 5) | (off >> 8)));
            } else {
                out.push_back(static_cast<char>((7 << 5) | (off >> 8)));
                out.push_back(static_cast<char>(coded - 7));
            }
            out.push_back(static_cast<char>(off & 0xFF));
            ip += length;
            lit_start = ip;
        } else {
            ++ip;
        }
    }
    flush(n);
    return out;
}

// ---------------------------------------------------------------------------
// shared helpers

bool pcd_host_big_endian() {
    const std::uint16_t v = 1;
    return *reinterpret_cast<const unsigned char*>(&v) == 0;
}

// Moves one element between host order and PCD's little-endian byte order (the swap is
// its own inverse, so one helper serves both directions).
void pcd_copy_le(void* pDst, const void* pSrc, std::size_t size) {
    if (pcd_host_big_endian())
        detail::bswap_copy(static_cast<char*>(pDst), static_cast<const char*>(pSrc),
                           static_cast<int>(size));
    else
        std::memcpy(pDst, pSrc, size);
}

DType pcd_dtype(char type, std::size_t size) {
    switch (type) {
        case 'F':
            if (size == 4)
                return DType::Float32;
            if (size == 8)
                return DType::Float64;
            break;
        case 'I':
            if (size == 1)
                return DType::Int8;
            if (size == 2)
                return DType::Int16;
            if (size == 4)
                return DType::Int32;
            if (size == 8)
                return DType::Int64;
            break;
        case 'U':
            if (size == 1)
                return DType::UInt8;
            if (size == 2)
                return DType::UInt16;
            if (size == 4)
                return DType::UInt32;
            if (size == 8)
                return DType::UInt64;
            break;
        default:
            break;
    }
    throw ReadError(std::string("PCD: unsupported field type ") + type + std::to_string(size));
}

bool pcd_is_float(DType dt) {
    return dt == DType::Float32 || dt == DType::Float64;
}

bool pcd_is_signed(DType dt) {
    return dt == DType::Int8 || dt == DType::Int16 || dt == DType::Int32 || dt == DType::Int64;
}

char pcd_type_letter(DType dt) {
    return pcd_is_float(dt) ? 'F' : pcd_is_signed(dt) ? 'I' : 'U';
}

// double -> float without the undefined behaviour of an out-of-range conversion.
float pcd_to_float(double v) {
    if (std::isfinite(v) && std::fabs(v) > static_cast<double>(FLT_MAX))
        return v < 0 ? -HUGE_VALF : HUGE_VALF;
    return static_cast<float>(v);
}

void pcd_store(NDArray& rColumn, std::size_t index, double value, long long ivalue) {
    const auto uvalue = static_cast<unsigned long long>(ivalue);
    switch (rColumn.Dtype()) {
        case DType::Float32:
            rColumn.As<float>()[index] = pcd_to_float(value);
            break;
        case DType::Float64:
            rColumn.As<double>()[index] = value;
            break;
        case DType::Int8:
            rColumn.As<std::int8_t>()[index] = static_cast<std::int8_t>(ivalue);
            break;
        case DType::Int16:
            rColumn.As<std::int16_t>()[index] = static_cast<std::int16_t>(ivalue);
            break;
        case DType::Int32:
            rColumn.As<std::int32_t>()[index] = static_cast<std::int32_t>(ivalue);
            break;
        case DType::Int64:
            rColumn.As<std::int64_t>()[index] = static_cast<std::int64_t>(ivalue);
            break;
        case DType::UInt8:
            rColumn.As<std::uint8_t>()[index] = static_cast<std::uint8_t>(uvalue);
            break;
        case DType::UInt16:
            rColumn.As<std::uint16_t>()[index] = static_cast<std::uint16_t>(uvalue);
            break;
        case DType::UInt32:
            rColumn.As<std::uint32_t>()[index] = static_cast<std::uint32_t>(uvalue);
            break;
        case DType::UInt64:
            rColumn.As<std::uint64_t>()[index] = static_cast<std::uint64_t>(uvalue);
            break;
    }
}

std::string pcd_upper(std::string s) {
    for (char& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::string pcd_strip(const std::string& rS) {
    std::size_t b = 0, e = rS.size();
    while (b < e && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(rS[e - 1])))
        --e;
    return rS.substr(b, e - b);
}

std::vector<std::string> pcd_words(const std::string& rS) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < rS.size()) {
        while (i < rS.size() && std::isspace(static_cast<unsigned char>(rS[i])))
            ++i;
        std::size_t j = i;
        while (j < rS.size() && !std::isspace(static_cast<unsigned char>(rS[j])))
            ++j;
        if (j > i)
            out.push_back(rS.substr(i, j - i));
        i = j;
    }
    return out;
}

bool pcd_parse_int(const std::string& rTok, long long& rOut) {
    if (rTok.empty())
        return false;
    errno = 0;
    char* end = nullptr;
    rOut = std::strtoll(rTok.c_str(), &end, 10);
    return errno == 0 && *end == '\0';
}

// ---------------------------------------------------------------------------
// reader

struct PcdHeader {
    std::vector<std::string> mFields;
    std::vector<char> mTypes;
    std::vector<std::size_t> mSizes;
    std::vector<std::size_t> mCounts;
    long long mWidth = 0;
    long long mHeight = 1;
    long long mPoints = 0;
    double mViewpoint[7] = {0, 0, 0, 1, 0, 0, 0};
    std::string mData;
    std::size_t mBody = 0;
};

PcdHeader pcd_parse_header(const std::string& rRaw) {
    std::unordered_map<std::string, std::vector<std::string>> header;
    std::size_t pos = 0;
    for (;;) {
        const std::size_t end = rRaw.find('\n', pos);
        if (end == std::string::npos)
            throw ReadError("PCD: no DATA line found in the header");
        const std::string line = pcd_strip(rRaw.substr(pos, end - pos));
        pos = end + 1;
        if (line.empty() || line[0] == '#')
            continue;
        std::vector<std::string> words = pcd_words(line);
        const std::string key = pcd_upper(words[0]);
        words.erase(words.begin());
        header[key] = std::move(words);
        if (key == "DATA")
            break;
    }

    PcdHeader h;
    h.mBody = pos;
    if (!header.count("FIELDS"))
        throw ReadError("PCD: the header has no FIELDS line");
    h.mFields = header["FIELDS"];
    const std::size_t nfields = h.mFields.size();
    const std::vector<std::string> sizes = header["SIZE"], types = header["TYPE"];
    const std::vector<std::string> counts =
        header.count("COUNT") ? header["COUNT"] : std::vector<std::string>(nfields, "1");
    if (!(sizes.size() == nfields && types.size() == nfields && counts.size() == nfields))
        throw ReadError("PCD: FIELDS, SIZE, TYPE and COUNT disagree on the field count");
    for (std::size_t i = 0; i < nfields; ++i) {
        long long s = 0, c = 0;
        if (!pcd_parse_int(sizes[i], s) || !pcd_parse_int(counts[i], c) || s < 0 || c < 0)
            throw ReadError("PCD: malformed SIZE/COUNT in the header");
        const std::string t = pcd_upper(types[i]);
        h.mTypes.push_back(t.empty() ? '?' : t[0]);
        h.mSizes.push_back(static_cast<std::size_t>(s));
        h.mCounts.push_back(static_cast<std::size_t>(c));
        pcd_dtype(h.mTypes.back(), h.mSizes.back());  // rejects an unsupported type/size pair
    }

    auto first_int = [&](const char* pKey, long long& rOut) {
        auto it = header.find(pKey);
        if (it == header.end())
            return false;
        if (it->second.empty() || !pcd_parse_int(it->second[0], rOut) || rOut < 0)
            throw ReadError("PCD: malformed WIDTH/HEIGHT/POINTS in the header");
        return true;
    };
    long long width = 0, points = 0;
    const bool has_width = first_int("WIDTH", width);
    if (!first_int("HEIGHT", h.mHeight))
        h.mHeight = 1;
    const bool has_points = first_int("POINTS", points);
    if (!has_points) {
        if (!has_width)
            throw ReadError("PCD: the header has neither POINTS nor WIDTH");
        points = width * h.mHeight;
    }
    if (!has_width)
        width = points;
    if (width * h.mHeight != points)
        log::warn("PCD: WIDTH*HEIGHT ({}*{}) != POINTS ({}); using POINTS", width, h.mHeight,
                  points);
    h.mWidth = width;
    h.mPoints = points;

    auto vp = header.find("VIEWPOINT");
    if (vp != header.end()) {
        if (vp->second.size() != 7)
            throw ReadError("PCD: VIEWPOINT needs 7 values (tx ty tz qw qx qy qz)");
        for (std::size_t i = 0; i < 7; ++i) {
            const char* end = nullptr;
            h.mViewpoint[i] = detail::parse_double(vp->second[i].c_str(), end);
            if (end == vp->second[i].c_str() || *end != '\0')
                throw ReadError("PCD: malformed VIEWPOINT");
        }
    }

    std::string mode = header["DATA"].empty() ? "" : header["DATA"][0];
    for (char& c : mode)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (mode != "ascii" && mode != "binary" && mode != "binary_compressed")
        throw ReadError("PCD: unknown DATA mode '" + mode + "'");
    h.mData = mode;
    return h;
}

using PcdColumns = std::vector<NDArray>;  // one (npoints, count) array per field

PcdColumns pcd_make_columns(const PcdHeader& rH, std::size_t npoints) {
    PcdColumns columns;
    for (std::size_t f = 0; f < rH.mFields.size(); ++f)
        columns.emplace_back(pcd_dtype(rH.mTypes[f], rH.mSizes[f]),
                             std::vector<std::size_t>{npoints, rH.mCounts[f]});
    return columns;
}

std::size_t pcd_row_values(const PcdHeader& rH) {
    std::size_t total = 0;
    for (std::size_t c : rH.mCounts)
        total += c;
    return total;
}

std::size_t pcd_row_bytes(const PcdHeader& rH) {
    std::size_t total = 0;
    for (std::size_t f = 0; f < rH.mFields.size(); ++f)
        total += rH.mSizes[f] * rH.mCounts[f];
    return total;
}

PcdColumns pcd_read_ascii(const PcdHeader& rH, const std::string& rRaw, std::size_t npoints) {
    const char* p = rRaw.data() + rH.mBody;
    const char* end = rRaw.data() + rRaw.size();
    const std::size_t expected = npoints * pcd_row_values(rH);

    std::size_t found = 0;
    for (const char* q = p; q < end;) {
        while (q < end && std::isspace(static_cast<unsigned char>(*q)))
            ++q;
        if (q >= end)
            break;
        ++found;
        while (q < end && !std::isspace(static_cast<unsigned char>(*q)))
            ++q;
    }
    if (found != expected)
        throw ReadError("PCD: expected " + std::to_string(expected) + " values, found " +
                        std::to_string(found));

    PcdColumns columns = pcd_make_columns(rH, npoints);
    std::string tok;
    for (std::size_t i = 0; i < npoints; ++i) {
        for (std::size_t f = 0; f < columns.size(); ++f) {
            const DType dt = columns[f].Dtype();
            for (std::size_t e = 0; e < rH.mCounts[f]; ++e) {
                while (std::isspace(static_cast<unsigned char>(*p)))
                    ++p;
                const char* q = p;
                while (q < end && !std::isspace(static_cast<unsigned char>(*q)))
                    ++q;
                tok.assign(p, static_cast<std::size_t>(q - p));
                p = q;
                const std::size_t index = i * rH.mCounts[f] + e;
                const char* stop = nullptr;
                if (pcd_is_float(dt)) {
                    const double v = detail::parse_double(tok.c_str(), stop);
                    if (stop == tok.c_str() || *stop != '\0')
                        throw ReadError("PCD: non-numeric value in the ASCII data");
                    pcd_store(columns[f], index, v, 0);
                } else {
                    errno = 0;
                    char* istop = nullptr;
                    const long long v =
                        pcd_is_signed(dt)
                            ? std::strtoll(tok.c_str(), &istop, 10)
                            : static_cast<long long>(std::strtoull(tok.c_str(), &istop, 10));
                    if (errno != 0 || istop == tok.c_str() || *istop != '\0')
                        throw ReadError("PCD: non-numeric value in the ASCII data");
                    pcd_store(columns[f], index, 0.0, v);
                }
            }
        }
    }
    return columns;
}

PcdColumns pcd_read_binary(const PcdHeader& rH, const std::string& rRaw, std::size_t npoints) {
    const std::size_t stride = pcd_row_bytes(rH);
    const std::size_t have = rRaw.size() - rH.mBody;
    if (stride != 0 && npoints > have / stride)
        throw ReadError("PCD: binary data is shorter than the header declares");
    PcdColumns columns = pcd_make_columns(rH, npoints);
    std::vector<std::size_t> offsets(columns.size());
    std::size_t off = 0;
    for (std::size_t f = 0; f < columns.size(); ++f) {
        offsets[f] = off;
        off += rH.mSizes[f] * rH.mCounts[f];
    }
    const char* base = rRaw.data() + rH.mBody;
    parallel_for_bw(npoints, [&](std::size_t i) {
        for (std::size_t f = 0; f < columns.size(); ++f) {
            const std::size_t size = rH.mSizes[f], count = rH.mCounts[f];
            char* dst = reinterpret_cast<char*>(columns[f].Data()) + i * count * size;
            const char* src = base + i * stride + offsets[f];
            for (std::size_t e = 0; e < count; ++e)
                pcd_copy_le(dst + e * size, src + e * size, size);
        }
    });
    return columns;
}

PcdColumns pcd_read_compressed(const PcdHeader& rH, const std::string& rRaw, std::size_t npoints) {
    const std::size_t have = rRaw.size() - rH.mBody;
    if (have < 8)
        throw ReadError("PCD: binary_compressed data is missing its size prefix");
    const char* payload = rRaw.data() + rH.mBody;
    std::uint32_t comp = 0, uncomp = 0;
    pcd_copy_le(&comp, payload, 4);
    pcd_copy_le(&uncomp, payload + 4, 4);
    if (static_cast<std::size_t>(comp) > have - 8)
        throw ReadError("PCD: binary_compressed data is shorter than its size prefix");
    const std::size_t expected = npoints * pcd_row_bytes(rH);
    if (uncomp != expected)
        throw ReadError("PCD: binary_compressed declares " + std::to_string(uncomp) +
                        " bytes, the header implies " + std::to_string(expected));
    const std::string raw =
        pcd_lzf_decompress(reinterpret_cast<const unsigned char*>(payload + 8), comp, uncomp);
    PcdColumns columns = pcd_make_columns(rH, npoints);
    std::size_t offset = 0;
    for (std::size_t f = 0; f < columns.size(); ++f) {
        const std::size_t size = rH.mSizes[f], values = npoints * rH.mCounts[f];
        char* dst = reinterpret_cast<char*>(columns[f].Data());
        for (std::size_t e = 0; e < values; ++e)
            pcd_copy_le(dst + e * size, raw.data() + offset + e * size, size);
        offset += values * size;
    }
    return columns;
}

NDArray pcd_take_rows(const NDArray& rSrc, const std::vector<std::size_t>& rRows) {
    std::vector<std::size_t> shape = rSrc.Shape();
    const std::size_t old_rows = shape.empty() ? 0 : shape[0];
    shape[0] = rRows.size();
    NDArray out(rSrc.Dtype(), shape);
    if (old_rows == 0 || rRows.empty())
        return out;
    const std::size_t row_bytes = rSrc.Nbytes() / old_rows;
    for (std::size_t k = 0; k < rRows.size(); ++k)
        std::memcpy(out.Data() + k * row_bytes, rSrc.Data() + rRows[k] * row_bytes, row_bytes);
    return out;
}

// PCL's rgb/rgba: a uint32 (0x00RRGGBB / 0xAARRGGBB), stored in a float32 slot when the
// TYPE is F. Unpacked by bit-cast, never by value.
NDArray pcd_unpack_colour(const NDArray& rColumn, char type, bool alpha) {
    const std::size_t n = rColumn.Shape()[0];
    NDArray out(DType::UInt8, {n, alpha ? std::size_t(4) : std::size_t(3)});
    std::uint8_t* dst = out.As<std::uint8_t>();
    for (std::size_t i = 0; i < n; ++i) {
        std::uint32_t packed = 0;
        if (type == 'F')
            std::memcpy(&packed, rColumn.Data() + i * sizeof(float), sizeof(packed));
        else
            packed = static_cast<std::uint32_t>(detail::read_int(rColumn, i));
        std::uint8_t* row = dst + i * (alpha ? 4 : 3);
        row[0] = static_cast<std::uint8_t>((packed >> 16) & 255);
        row[1] = static_cast<std::uint8_t>((packed >> 8) & 255);
        row[2] = static_cast<std::uint8_t>(packed & 255);
        if (alpha)
            row[3] = static_cast<std::uint8_t>((packed >> 24) & 255);
    }
    return out;
}

Mesh pcd_build_mesh(const PcdHeader& rH, PcdColumns& rColumns, bool drop_invalid) {
    const std::size_t npoints = static_cast<std::size_t>(rH.mPoints);
    std::unordered_map<std::string, std::size_t> by_name;  // a repeated name: the last wins
    for (std::size_t f = 0; f < rH.mFields.size(); ++f)
        by_name[rH.mFields[f]] = f;

    const char axes[3] = {'x', 'y', 'z'};
    bool single = true;
    std::size_t axis_field[3];
    for (int a = 0; a < 3; ++a) {
        auto it = by_name.find(std::string(1, axes[a]));
        if (it == by_name.end() || rH.mCounts[it->second] != 1)
            throw ReadError(std::string("PCD: the cloud has no scalar '") + axes[a] + "' field");
        axis_field[a] = it->second;
        single = single && rColumns[it->second].Dtype() == DType::Float32;
    }
    NDArray points(single ? DType::Float32 : DType::Float64, {npoints, std::size_t(3)});
    for (int a = 0; a < 3; ++a)
        for (std::size_t i = 0; i < npoints; ++i)
            pcd_store(points, i * 3 + a, detail::read_double(rColumns[axis_field[a]], i), 0);

    std::unordered_map<std::string, NDArray> point_data;
    std::unordered_set<std::string> consumed = {"x", "y", "z"};
    const char* normal_axes[3] = {"normal_x", "normal_y", "normal_z"};
    bool have_normals = true;
    for (const char* name : normal_axes) {
        auto it = by_name.find(name);
        have_normals = have_normals && it != by_name.end() && rH.mCounts[it->second] == 1;
    }
    if (have_normals) {
        bool single_n = true;
        for (const char* name : normal_axes)
            single_n = single_n && rColumns[by_name[name]].Dtype() == DType::Float32;
        NDArray normals(single_n ? DType::Float32 : DType::Float64, {npoints, std::size_t(3)});
        for (int a = 0; a < 3; ++a)
            for (std::size_t i = 0; i < npoints; ++i)
                pcd_store(normals, i * 3 + a,
                          detail::read_double(rColumns[by_name[normal_axes[a]]], i), 0);
        point_data.emplace("normals", std::move(normals));
        for (const char* name : normal_axes)
            consumed.insert(name);
    }
    for (const std::string& name : rH.mFields) {
        if (consumed.count(name) || name == "_" || point_data.count(name))
            continue;
        const std::size_t f = by_name[name];
        NDArray& column = rColumns[f];
        if ((name == "rgb" || name == "rgba") && rH.mCounts[f] == 1 && rH.mSizes[f] == 4) {
            point_data.emplace(name, pcd_unpack_colour(column, rH.mTypes[f], name == "rgba"));
        } else {
            if (rH.mCounts[f] == 1)
                column.Reshape({npoints});
            point_data.emplace(name, std::move(column));
        }
    }

    const bool organised = rH.mHeight > 1 && rH.mWidth * rH.mHeight == rH.mPoints;
    bool keep_organisation = organised;
    std::size_t n = npoints;
    if (drop_invalid) {
        std::vector<std::size_t> keep;
        keep.reserve(npoints);
        for (std::size_t i = 0; i < npoints; ++i) {
            bool ok = true;
            for (int a = 0; a < 3; ++a)
                ok = ok && std::isfinite(detail::read_double(points, i * 3 + a));
            if (ok)
                keep.push_back(i);
        }
        if (keep.size() != npoints) {
            points = pcd_take_rows(points, keep);
            for (auto& entry : point_data)
                entry.second = pcd_take_rows(entry.second, keep);
            keep_organisation = false;
            n = keep.size();
        }
    }

    Mesh mesh;
    mesh.AssignPoints(std::move(points));
    NDArray cells(DType::Int64, {n, std::size_t(1)});
    for (std::size_t i = 0; i < n; ++i)
        cells.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    mesh.AddCellBlock("vertex", std::move(cells));
    std::vector<std::string> names;
    for (const auto& entry : point_data)
        names.push_back(entry.first);
    std::sort(names.begin(), names.end());
    for (const std::string& name : names)
        mesh.AddPointData(name, std::move(point_data.at(name)));
    if (keep_organisation) {
        NDArray w(DType::Int64, {std::size_t(1)}), h(DType::Int64, {std::size_t(1)});
        w.As<std::int64_t>()[0] = rH.mWidth;
        h.As<std::int64_t>()[0] = rH.mHeight;
        mesh.AddFieldData("pcd:width", std::move(w));
        mesh.AddFieldData("pcd:height", std::move(h));
    }
    if (!std::equal(rH.mViewpoint, rH.mViewpoint + 7, kPcdIdentityViewpoint)) {
        NDArray vp(DType::Float64, {std::size_t(7)});
        std::memcpy(vp.Data(), rH.mViewpoint, sizeof(rH.mViewpoint));
        mesh.AddFieldData("pcd:viewpoint", std::move(vp));
    }
    return mesh;
}

// ---------------------------------------------------------------------------
// writer

struct PcdOutField {
    std::string mName;
    char mType;
    std::size_t mSize;
    std::size_t mCount;
    NDArray mData;  // (n, count), already in the stored dtype
};

std::string pcd_unique(std::string name, std::unordered_set<std::string>& rUsed) {
    for (char& c : name)
        if (std::isspace(static_cast<unsigned char>(c)))
            c = '_';
    if (name.empty())
        name = "field";
    std::string candidate = name;
    for (int k = 2; rUsed.count(candidate); ++k)
        candidate = name + "_" + std::to_string(k);
    rUsed.insert(candidate);
    return candidate;
}

// Columns [first, first + count) of a row-major (n, stride) array, as dtype `target`.
NDArray pcd_extract(const NDArray& rSrc, std::size_t n, std::size_t stride, std::size_t first,
                    std::size_t count, DType target) {
    NDArray out(target, {n, count});
    const bool to_float = pcd_is_float(target);
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t j = 0; j < count; ++j) {
            const std::size_t at = r * stride + first + j;
            if (to_float)
                pcd_store(out, r * count + j, detail::read_double(rSrc, at), 0);
            else
                pcd_store(out, r * count + j, 0.0, detail::read_int(rSrc, at));
        }
    return out;
}

NDArray pcd_pack_colour(const NDArray& rSrc, std::size_t n, std::size_t width, bool as_float) {
    NDArray out(as_float ? DType::Float32 : DType::UInt32, {n, std::size_t(1)});
    for (std::size_t r = 0; r < n; ++r) {
        std::uint32_t channel[4] = {0, 0, 0, 0};
        for (std::size_t j = 0; j < width; ++j) {
            double v = std::nearbyint(detail::read_double(rSrc, r * width + j));
            v = std::isnan(v) ? 0.0 : std::min(255.0, std::max(0.0, v));
            channel[j] = static_cast<std::uint32_t>(v);
        }
        std::uint32_t packed = (channel[0] << 16) | (channel[1] << 8) | channel[2];
        if (width == 4)
            packed |= channel[3] << 24;
        std::memcpy(out.Data() + r * 4, &packed, 4);
    }
    return out;
}

std::string pcd_ascii_value(const NDArray& rData, std::size_t index, char type) {
    char buf[64];
    if (type == 'F') {
        const double v = detail::read_double(rData, index);
        if (std::isnan(v))
            return "nan";
        detail::snprintf_c(buf, sizeof(buf), rData.Dtype() == DType::Float32 ? "%.9g" : "%.17g", v);
        return buf;
    }
    if (rData.Dtype() == DType::UInt64)
        return std::to_string(rData.As<std::uint64_t>()[index]);
    return std::to_string(detail::read_int(rData, index));
}

}  // namespace

Mesh read_pcd(const std::string& rPath, const PcdReadOptions& rOptions) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    const std::string raw((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    const PcdHeader header = pcd_parse_header(raw);
    const std::size_t npoints = static_cast<std::size_t>(header.mPoints);
    PcdColumns columns = header.mData == "ascii"    ? pcd_read_ascii(header, raw, npoints)
                         : header.mData == "binary" ? pcd_read_binary(header, raw, npoints)
                                                    : pcd_read_compressed(header, raw, npoints);
    return pcd_build_mesh(header, columns, rOptions.mDropInvalid);
}

void write_pcd(const std::string& rPath, const Mesh& rMesh, PcdData data, bool float64_points) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();
    if (dim < 3) {
        log::warn("PCD requires 3D points; padding with zeros.");
        detail::provenance_note("point-padding", "points padded with zero coordinates to 3D");
    }
    if (!float64_points && points.Dtype() != DType::Float32) {
        bool changed = false;
        for (std::size_t i = 0; i < n * dim && !changed; ++i) {
            const double v = detail::read_double(points, i);
            if (std::isnan(v))
                continue;
            changed = static_cast<double>(pcd_to_float(v)) != v;
        }
        if (changed)
            detail::provenance_note(
                "dtype", "float64 coordinates written as float32 (point_dtype='float32')");
    }

    std::string skipped;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() == "vertex")
            continue;
        if (!skipped.empty())
            skipped += ", ";
        skipped += cb.Type();
    }
    if (!skipped.empty()) {
        log::warn("PCD holds points only. Skipping {} cells.", skipped);
        detail::provenance_note("cells-dropped",
                                "cell block(s) of type " + skipped + " have no PCD equivalent");
    }
    if (rMesh.NumCellData() > 0)
        detail::provenance_note("data-dropped", "cell data has no PCD equivalent");
    std::string extra;
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.rfind("pcd:", 0) == 0)
            continue;
        if (!extra.empty())
            extra += ", ";
        extra += name;
    }
    if (!extra.empty())
        detail::provenance_note("data-dropped", "field data has no PCD equivalent: " + extra);

    const DType float_type = float64_points ? DType::Float64 : DType::Float32;
    const std::size_t float_size = float64_points ? 8 : 4;
    std::unordered_set<std::string> used = {"x", "y", "z"};
    std::vector<PcdOutField> fields;
    const char* axis_names[3] = {"x", "y", "z"};
    for (std::size_t a = 0; a < 3; ++a) {
        if (a < dim)
            fields.push_back(
                {axis_names[a], 'F', float_size, 1, pcd_extract(points, n, dim, a, 1, float_type)});
        else
            fields.push_back(
                {axis_names[a], 'F', float_size, 1, NDArray(float_type, {n, std::size_t(1)})});
    }
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& array = rMesh.PointData(name);
        const auto& shape = array.Shape();
        const bool two_d = shape.size() == 2 && shape[0] == n;
        if (name == "normals" && two_d && shape[1] == 3) {
            const char* axes[3] = {"normal_x", "normal_y", "normal_z"};
            for (std::size_t a = 0; a < 3; ++a)
                fields.push_back({pcd_unique(axes[a], used), 'F', float_size, 1,
                                  pcd_extract(array, n, 3, a, 1, float_type)});
            continue;
        }
        if ((name == "rgb" || name == "rgba") && two_d && shape[1] == (name == "rgba" ? 4u : 3u)) {
            const bool is_rgb = name == "rgb";
            fields.push_back({pcd_unique(name, used), is_rgb ? 'F' : 'U', 4, 1,
                              pcd_pack_colour(array, n, shape[1], is_rgb)});
            continue;
        }
        std::size_t count = 1;
        for (std::size_t d = 1; d < shape.size(); ++d)
            count *= shape[d];
        DType stored = array.Dtype();
        if ((name == "curvature" || name == "intensity") && pcd_is_float(stored))
            stored = float_type;
        fields.push_back({pcd_unique(name, used), pcd_type_letter(stored), dtype_size(stored),
                          count, pcd_extract(array, n, count, 0, count, stored)});
    }

    long long width = static_cast<long long>(n), height = 1;
    if (rMesh.HasFieldData("pcd:width") && rMesh.HasFieldData("pcd:height") && n) {
        const long long w = detail::read_int(rMesh.FieldData("pcd:width"), 0);
        const long long h = detail::read_int(rMesh.FieldData("pcd:height"), 0);
        if (w * h == static_cast<long long>(n)) {
            width = w;
            height = h;
        }
    }
    double viewpoint[7];
    std::copy(kPcdIdentityViewpoint, kPcdIdentityViewpoint + 7, viewpoint);
    if (rMesh.HasFieldData("pcd:viewpoint")) {
        const NDArray& vp = rMesh.FieldData("pcd:viewpoint");
        if (vp.Size() == 7)
            for (std::size_t i = 0; i < 7; ++i)
                viewpoint[i] = detail::read_double(vp, i);
    }

    const char* mode = data == PcdData::Ascii    ? "ascii"
                       : data == PcdData::Binary ? "binary"
                                                 : "binary_compressed";
    std::string head = "# .PCD v0.7 - Point Cloud Data file format\n";
    head += detail::provenance_render_lines(detail::SlotTier::Block, "# ");
    head += "VERSION 0.7\nFIELDS";
    for (const auto& f : fields)
        head += " " + f.mName;
    head += "\nSIZE";
    for (const auto& f : fields)
        head += " " + std::to_string(f.mSize);
    head += "\nTYPE";
    for (const auto& f : fields)
        head += std::string(" ") + f.mType;
    head += "\nCOUNT";
    for (const auto& f : fields)
        head += " " + std::to_string(f.mCount);
    head +=
        "\nWIDTH " + std::to_string(width) + "\nHEIGHT " + std::to_string(height) + "\nVIEWPOINT";
    char buf[64];
    for (double v : viewpoint) {
        detail::snprintf_c(buf, sizeof(buf), "%.17g", v);
        head += std::string(" ") + buf;
    }
    head += "\nPOINTS " + std::to_string(n) + "\nDATA " + mode + "\n";

    std::string body;
    if (data == PcdData::Ascii) {
        for (std::size_t i = 0; i < n; ++i) {
            bool first = true;
            for (const auto& f : fields)
                for (std::size_t e = 0; e < f.mCount; ++e) {
                    if (!first)
                        body += ' ';
                    first = false;
                    body += pcd_ascii_value(f.mData, i * f.mCount + e, f.mType);
                }
            body += '\n';
        }
    } else if (data == PcdData::Binary) {
        std::size_t stride = 0;
        std::vector<std::size_t> offsets;
        for (const auto& f : fields) {
            offsets.push_back(stride);
            stride += f.mSize * f.mCount;
        }
        body.assign(n * stride, '\0');
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t k = 0; k < fields.size(); ++k) {
                const auto& f = fields[k];
                const char* src =
                    reinterpret_cast<const char*>(f.mData.Data()) + i * f.mCount * f.mSize;
                char* dst = &body[i * stride + offsets[k]];
                for (std::size_t e = 0; e < f.mCount; ++e)
                    pcd_copy_le(dst + e * f.mSize, src + e * f.mSize, f.mSize);
            }
        });
    } else {
        std::string blocks;
        for (const auto& f : fields) {
            const std::size_t values = n * f.mCount;
            const std::size_t at = blocks.size();
            blocks.resize(at + values * f.mSize);
            const char* src = reinterpret_cast<const char*>(f.mData.Data());
            for (std::size_t e = 0; e < values; ++e)
                pcd_copy_le(&blocks[at + e * f.mSize], src + e * f.mSize, f.mSize);
        }
        const std::string packed =
            pcd_lzf_compress(reinterpret_cast<const unsigned char*>(blocks.data()), blocks.size());
        const std::uint32_t sizes[2] = {static_cast<std::uint32_t>(packed.size()),
                                        static_cast<std::uint32_t>(blocks.size())};
        body.resize(8);
        pcd_copy_le(&body[0], &sizes[0], 4);
        pcd_copy_le(&body[4], &sizes[1], 4);
        body += packed;
    }

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);
    os.write(head.data(), static_cast<std::streamsize>(head.size()));
    os.write(body.data(), static_cast<std::streamsize>(body.size()));
}

}  // namespace meshioplusplus
