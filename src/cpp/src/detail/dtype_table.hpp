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

#pragma once

/**
 * @file detail/dtype_table.hpp
 * @brief `DTypeUtils`: one value per `DType`, read from a table in enumerator
 * order instead of a ten-case `switch` (roadmap §8).
 *
 * Every format that names element types (netCDF, HDF5, XDMF, VTK, PLY) mapped
 * a `DType` with the same ten-case switch falling back to the `Float64` entry;
 * each now lists its ten names once, in the order `DType` declares them, and
 * `DTypeUtils::Pick` reads them. A static assertion pins that order.
 *
 * A **core-private** header (not installed): no API or ABI surface.
 */

// System includes
#include <array>
#include <cstddef>

// Project includes
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus::detail {

/// Static per-`DType` table lookups.
struct DTypeUtils {
    /// How many `DType` enumerators there are.
    static constexpr std::size_t kNumDTypes = 10;

    /// A value per `DType`, in enumerator order: Float32, Float64, Int8, Int16,
    /// Int32, Int64, UInt8, UInt16, UInt32, UInt64.
    template <class TValue>
    using Table = std::array<TValue, kNumDTypes>;

    /// `rTable`'s entry for `Dt`; a value outside the enumeration (a corrupted
    /// or future `DType`) reads the `Float64` entry, as the switches did.
    template <class TValue>
    static TValue Pick(DType Dt, const Table<TValue>& rTable) {
        const auto i = static_cast<std::size_t>(Dt);
        return i < kNumDTypes ? rTable[i] : rTable[static_cast<std::size_t>(DType::Float64)];
    }
};

static_assert(static_cast<int>(DType::Float32) == 0 && static_cast<int>(DType::Float64) == 1 &&
                  static_cast<int>(DType::Int8) == 2 && static_cast<int>(DType::Int16) == 3 &&
                  static_cast<int>(DType::Int32) == 4 && static_cast<int>(DType::Int64) == 5 &&
                  static_cast<int>(DType::UInt8) == 6 && static_cast<int>(DType::UInt16) == 7 &&
                  static_cast<int>(DType::UInt32) == 8 && static_cast<int>(DType::UInt64) == 9,
              "DTypeUtils tables list the DType enumerators in declaration order");

}  // namespace meshioplusplus::detail
