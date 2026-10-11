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
 * @file detail/fortran_records.hpp
 * @brief Split a Fortran sequential unformatted file into its records, and
 * frame records to write one.
 *
 * A Fortran `WRITE` to a sequential unformatted unit frames each record with its
 * byte length before and after it. The marker is 4 bytes with gfortran and
 * Intel Fortran (8 with old g77 `-frecord-marker=8` builds) and is stored in the
 * writing machine's byte order. `sniff_fortran_records` finds the layout from
 * the first record -- a marker whose twin sits exactly after the payload -- and
 * `fortran_records` returns every record's payload span, checking each twin.
 * Abaqus `.fil` is the first user; OP2 and EnSight Fortran binary share the
 * framing, and the OpenRadioss time-history file (`radioss_th`, big-endian
 * 4-byte markers) reads and writes through it.
 *
 * The writing half is `fortran_append_record`, which frames one payload in
 * either marker width and byte order, and the `fortran_put_*` helpers that
 * append a scalar in the file's byte order to build a payload. The Python twin
 * is `meshioplusplus._fortran_records`.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {
namespace detail {

/// The framing of a Fortran sequential unformatted file.
struct FortranRecordLayout {
    int mMarkerBytes = 4;     ///< 4 or 8
    bool mBigEndian = false;  ///< the markers' byte order
};

/// One record's payload: `[mOffset, mOffset + mSize)` in the file buffer.
struct FortranRecord {
    std::size_t mOffset = 0;
    std::size_t mSize = 0;
};

/**
 * @brief The framing of the file that starts at @p pData, if its first record is
 * a well-formed Fortran record.
 *
 * Tries 4-byte then 8-byte markers, each in both byte orders; a candidate is
 * accepted when its length is positive, fits the buffer and the same value
 * follows the payload.
 * @param pData The start of the file.
 * @param Size Bytes available.
 * @return the layout, or nothing when no candidate frames the first record.
 */
MESHIOPLUSPLUS_API std::optional<FortranRecordLayout> sniff_fortran_records(const char* pData,
                                                                            std::size_t Size);

/**
 * @brief Every record of a Fortran sequential unformatted file.
 * @param pData The file buffer.
 * @param Size Its length.
 * @param rLayout The framing (see `sniff_fortran_records`).
 * @param rWhat The format name used in error messages.
 * @return the payload spans in file order.
 * @throws ReadError when a record runs past the end or its trailing marker does
 *         not match the leading one.
 */
MESHIOPLUSPLUS_API std::vector<FortranRecord> fortran_records(const char* pData, std::size_t Size,
                                                              const FortranRecordLayout& rLayout,
                                                              const std::string& rWhat);

/**
 * @brief `fortran_records`, optionally tolerant of a truncated tail.
 *
 * With @p Lenient the first record that runs past the end of the buffer, or whose
 * trailing marker does not match its leading one, ends the list instead of
 * throwing: the records before it are returned. A solver killed mid-write leaves
 * such a file (the OpenRadioss time-history reader drops the incomplete last
 * output this way). With `Lenient` false this is the four-argument overload.
 * @param pData The file buffer.
 * @param Size Its length.
 * @param rLayout The framing (see `sniff_fortran_records`).
 * @param rWhat The format name used in error messages.
 * @param Lenient Stop at the first malformed record rather than throw.
 * @return the payload spans in file order.
 * @throws ReadError on a malformed record unless @p Lenient.
 */
MESHIOPLUSPLUS_API std::vector<FortranRecord> fortran_records(const char* pData, std::size_t Size,
                                                              const FortranRecordLayout& rLayout,
                                                              const std::string& rWhat,
                                                              bool Lenient);

/**
 * @brief Append one framed record to @p rOut: the length marker, the payload,
 * the same marker again.
 * @param rOut The file being built.
 * @param pPayload The payload bytes (may be null when @p Size is 0).
 * @param Size The payload length.
 * @param rLayout The marker width and byte order.
 * @throws WriteError when @p Size does not fit the marker.
 */
MESHIOPLUSPLUS_API void fortran_append_record(std::string& rOut, const char* pPayload,
                                              std::size_t Size, const FortranRecordLayout& rLayout);

/// Append @p Value as 4 bytes in the given byte order (payload building).
MESHIOPLUSPLUS_API void fortran_put_int32(std::string& rBuf, std::int32_t Value, bool BigEndian);
/// Append @p Value as 8 bytes in the given byte order.
MESHIOPLUSPLUS_API void fortran_put_int64(std::string& rBuf, std::int64_t Value, bool BigEndian);
/// Append @p Value as an IEEE single in the given byte order.
MESHIOPLUSPLUS_API void fortran_put_float32(std::string& rBuf, float Value, bool BigEndian);
/// Append @p Value as an IEEE double in the given byte order.
MESHIOPLUSPLUS_API void fortran_put_float64(std::string& rBuf, double Value, bool BigEndian);

}  // namespace detail
}  // namespace meshioplusplus
