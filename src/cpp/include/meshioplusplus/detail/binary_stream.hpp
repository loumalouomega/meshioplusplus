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
 * @file detail/binary_stream.hpp
 * @brief A bounds-checked cursor over a byte buffer in a chosen byte order.
 *
 * Binary readers that walk a file value by value (libMesh `.xdr`, Abaqus `.fil`)
 * read through `ByteCursor`: every read checks the remaining length and throws a
 * `ReadError` naming the format instead of running past the buffer, and the
 * file's byte order is fixed once at construction. The cursor never owns the
 * buffer.
 */

// System includes
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>

// Project includes
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

/// Sequential reads of fixed-width values from a buffer in one byte order.
class ByteCursor {
public:
    /**
     * @brief A cursor at the start of `[pData, pData + Size)`.
     * @param pData The buffer (not owned; must outlive the cursor).
     * @param Size Its length in bytes.
     * @param BigEndian Whether the values are stored big-endian.
     * @param rWhat The format name used in error messages.
     */
    ByteCursor(const char* pData, std::size_t Size, bool BigEndian, std::string rWhat)
        : mpData(pData),
          mSize(Size),
          mSwap(BigEndian != (std::endian::native == std::endian::big)),
          mWhat(std::move(rWhat)) {}

    /// Bytes consumed so far.
    std::size_t Offset() const { return mPos; }
    /// Bytes left.
    std::size_t Remaining() const { return mSize - mPos; }
    /// Whether every byte has been consumed.
    bool AtEnd() const { return mPos >= mSize; }
    /// Move to absolute offset @p Pos (at most the buffer length).
    void Seek(std::size_t Pos) {
        if (Pos > mSize)
            Fail(Pos - mPos);
        mPos = Pos;
    }
    /// Skip @p N bytes.
    void Skip(std::size_t N) {
        Need(N);
        mPos += N;
    }

    std::uint32_t U32() { return Read<std::uint32_t>(); }
    std::int32_t I32() { return static_cast<std::int32_t>(Read<std::uint32_t>()); }
    std::uint64_t U64() { return Read<std::uint64_t>(); }
    std::int64_t I64() { return static_cast<std::int64_t>(Read<std::uint64_t>()); }
    float F32() { return std::bit_cast<float>(Read<std::uint32_t>()); }
    double F64() { return std::bit_cast<double>(Read<std::uint64_t>()); }

    /// @p N raw bytes as a string.
    std::string Bytes(std::size_t N) {
        Need(N);
        std::string s(mpData + mPos, N);
        mPos += N;
        return s;
    }

    /// Throw unless @p N more bytes are available.
    void Need(std::size_t N) const {
        if (N > mSize - mPos)
            Fail(N);
    }

private:
    template <class T>
    T Read() {
        Need(sizeof(T));
        T v;
        std::memcpy(&v, mpData + mPos, sizeof(T));
        mPos += sizeof(T);
        if (mSwap) {
            if constexpr (sizeof(T) == 4)
                v = bswap32(v);
            else
                v = bswap64(v);
        }
        return v;
    }

    [[noreturn]] void Fail(std::size_t N) const {
        throw ReadError(mWhat + ": file is truncated (needs " + std::to_string(N) +
                        " more bytes at offset " + std::to_string(mPos) + " of " +
                        std::to_string(mSize) + ")");
    }

    const char* mpData;
    std::size_t mSize;
    std::size_t mPos = 0;
    bool mSwap;
    std::string mWhat;
};

}  // namespace detail
}  // namespace meshioplusplus
