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
 * @brief Split a Fortran sequential unformatted file into its records.
 *
 * A Fortran `WRITE` to a sequential unformatted unit frames each record with its
 * byte length before and after it. The marker is 4 bytes with gfortran and
 * Intel Fortran (8 with old g77 `-frecord-marker=8` builds) and is stored in the
 * writing machine's byte order. `sniff_fortran_records` finds the layout from
 * the first record -- a marker whose twin sits exactly after the payload -- and
 * `fortran_records` returns every record's payload span, checking each twin.
 * Abaqus `.fil` is the first user; OP2 and EnSight Fortran binary share the
 * framing.
 */

// System includes
#include <cstddef>
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

}  // namespace detail
}  // namespace meshioplusplus
