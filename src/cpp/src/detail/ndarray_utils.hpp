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
 * @file detail/ndarray_utils.hpp
 * @brief `NDArrayUtils`: owning 1-D `NDArray`s built from `std::vector`s, the
 * helper half a dozen operations and formats each wrote out (roadmap §8).
 *
 * A Kratos-style utility class of static inline members. A core-private
 * header (not installed): no API or ABI surface.
 */

// System includes
#include <cstdint>
#include <cstring>
#include <vector>

// Project includes
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus::detail {

/// Static constructors of owning `NDArray`s.
struct NDArrayUtils {
    /// @p rValues as an owning 1-D `Int64` array of the same length.
    static NDArray Int64Vector(const std::vector<std::int64_t>& rValues) {
        NDArray out = NDArray::Uninit(DType::Int64, {rValues.size()});
        if (!rValues.empty())
            std::memcpy(out.Data(), rValues.data(), rValues.size() * sizeof(std::int64_t));
        return out;
    }

    /// @p rValues as an owning 1-D `Float64` array of the same length.
    static NDArray Float64Vector(const std::vector<double>& rValues) {
        NDArray out = NDArray::Uninit(DType::Float64, {rValues.size()});
        if (!rValues.empty())
            std::memcpy(out.Data(), rValues.data(), rValues.size() * sizeof(double));
        return out;
    }
};

}  // namespace meshioplusplus::detail
