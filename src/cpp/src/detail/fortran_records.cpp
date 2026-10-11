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
#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/detail/binary_stream.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The marker at `Offset`, or -1 when it does not fit.
std::int64_t fortran_marker(const char* pData, std::size_t Size, std::size_t Offset,
                            const FortranRecordLayout& rLayout) {
    const std::size_t w = static_cast<std::size_t>(rLayout.mMarkerBytes);
    if (Offset > Size || Size - Offset < w)
        return -1;
    ByteCursor c(pData + Offset, w, rLayout.mBigEndian, "Fortran records");
    const std::int64_t v = w == 4 ? static_cast<std::int64_t>(c.I32()) : c.I64();
    return v;
}

}  // namespace

std::optional<FortranRecordLayout> sniff_fortran_records(const char* pData, std::size_t Size) {
    for (int bytes : {4, 8}) {
        for (bool big : {false, true}) {
            const FortranRecordLayout layout{bytes, big};
            const std::int64_t n = fortran_marker(pData, Size, 0, layout);
            if (n <= 0)
                continue;
            const std::size_t w = static_cast<std::size_t>(bytes);
            const std::size_t len = static_cast<std::size_t>(n);
            if (len > Size - w)
                continue;
            if (fortran_marker(pData, Size, w + len, layout) == n)
                return layout;
        }
    }
    return std::nullopt;
}

std::vector<FortranRecord> fortran_records(const char* pData, std::size_t Size,
                                           const FortranRecordLayout& rLayout,
                                           const std::string& rWhat) {
    return fortran_records(pData, Size, rLayout, rWhat, /*Lenient=*/false);
}

std::vector<FortranRecord> fortran_records(const char* pData, std::size_t Size,
                                           const FortranRecordLayout& rLayout,
                                           const std::string& rWhat, bool Lenient) {
    std::vector<FortranRecord> out;
    const std::size_t w = static_cast<std::size_t>(rLayout.mMarkerBytes);
    std::size_t pos = 0;
    while (pos < Size) {
        const std::int64_t n = fortran_marker(pData, Size, pos, rLayout);
        if (n < 0 || static_cast<std::uint64_t>(n) > Size - pos - w) {
            if (Lenient)
                break;
            throw ReadError(rWhat + ": Fortran record at offset " + std::to_string(pos) +
                            " runs past the end of the file");
        }
        const std::size_t len = static_cast<std::size_t>(n);
        if (fortran_marker(pData, Size, pos + w + len, rLayout) != n) {
            if (Lenient)
                break;
            throw ReadError(rWhat + ": Fortran record at offset " + std::to_string(pos) +
                            " has mismatched length markers");
        }
        out.push_back({pos + w, len});
        pos += 2 * w + len;
    }
    return out;
}

namespace {

// `Bytes` low-order bytes of `Value`, in the requested byte order.
void fortran_put_bytes(std::string& rBuf, std::uint64_t Value, int Bytes, bool BigEndian) {
    for (int i = 0; i < Bytes; ++i) {
        const int shift = BigEndian ? 8 * (Bytes - 1 - i) : 8 * i;
        rBuf.push_back(static_cast<char>(static_cast<unsigned char>((Value >> shift) & 0xFFU)));
    }
}

}  // namespace

void fortran_put_int32(std::string& rBuf, std::int32_t Value, bool BigEndian) {
    fortran_put_bytes(rBuf, static_cast<std::uint32_t>(Value), 4, BigEndian);
}

void fortran_put_int64(std::string& rBuf, std::int64_t Value, bool BigEndian) {
    fortran_put_bytes(rBuf, static_cast<std::uint64_t>(Value), 8, BigEndian);
}

void fortran_put_float32(std::string& rBuf, float Value, bool BigEndian) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, &Value, sizeof bits);
    fortran_put_bytes(rBuf, bits, 4, BigEndian);
}

void fortran_put_float64(std::string& rBuf, double Value, bool BigEndian) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &Value, sizeof bits);
    fortran_put_bytes(rBuf, bits, 8, BigEndian);
}

void fortran_append_record(std::string& rOut, const char* pPayload, std::size_t Size,
                           const FortranRecordLayout& rLayout) {
    if (rLayout.mMarkerBytes == 4 && Size > static_cast<std::size_t>(INT32_MAX))
        throw WriteError("Fortran records: a record of " + std::to_string(Size) +
                         " bytes does not fit a 4-byte length marker");
    rOut.reserve(rOut.size() + Size + 2 * static_cast<std::size_t>(rLayout.mMarkerBytes));
    fortran_put_bytes(rOut, Size, rLayout.mMarkerBytes, rLayout.mBigEndian);
    if (Size > 0)
        rOut.append(pPayload, Size);
    fortran_put_bytes(rOut, Size, rLayout.mMarkerBytes, rLayout.mBigEndian);
}

}  // namespace detail
}  // namespace meshioplusplus
