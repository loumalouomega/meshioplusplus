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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

/// @file integer_precision.hpp
/// @brief A warning for integer arrays a format stores as floating point.
///
/// Some formats have no integer data field: EnSight, FRD, Marc T19, glTF, GiD and
/// the Radioss and XPLT animations write every array as single-precision floats or
/// as a few significant digits of text. An integer above 2**24 is no longer exact
/// there. That is the format's own limit, not a defect, but it was silent
/// (roadmap 3.19); `warn_integer_precision` says so once per array. The Python
/// twin is `meshioplusplus._helpers._warn_integer_precision`, and
/// `tests/python/test_integer_precision.py` holds the two format tables equal.

#ifndef MESHIOPLUSPLUS_DETAIL_INTEGER_PRECISION_HPP
#define MESHIOPLUSPLUS_DETAIL_INTEGER_PRECISION_HPP

#include <string>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// The largest magnitude every integer up to which a float32 holds exactly.
inline constexpr long long kExactFloatInteger = 1LL << 24;

/// Warn, once per array, for each integer point or cell array of @p rMesh whose
/// values the format @p rFormat cannot keep exactly. Does nothing for a format
/// that stores integers as integers, for float arrays and for small values.
void warn_integer_precision(const std::string& rFormat, const Mesh& rMesh);

}  // namespace detail
}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_DETAIL_INTEGER_PRECISION_HPP
