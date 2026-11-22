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

#pragma once

/**
 * @file detail/fortran_number.hpp
 * @brief The Fortran double-precision exponent letter (`1.5D+03`), turned into
 * the one `parse_double` reads, written once (roadmap 6.7).
 *
 * Fortran writes `D` (or `d`) where C writes `e`; the text readers of formats
 * that Fortran solvers wrote each looped over their token to replace it. A
 * **core-private** header (not installed). `detail/keyword_card.cpp` keeps its
 * own bounded number, which also reads an exponent written without a letter.
 */

// System includes
#include <cstddef>
#include <string>

// Project includes
#include "always_inline.hpp"

namespace meshioplusplus::detail {

/// Replace every `D` or `d` of `[pFirst, pFirst + N)` by `Exponent` (`'E'` or `'e'`).
MIO_ALWAYS_INLINE void fortran_exponent_to_e(char* pFirst, std::size_t N, char Exponent = 'E') {
    for (std::size_t i = 0; i < N; ++i)
        if (pFirst[i] == 'D' || pFirst[i] == 'd')
            pFirst[i] = Exponent;
}

/// Replace every `D` or `d` of a number's text by `Exponent`.
MIO_ALWAYS_INLINE void fortran_exponent_to_e(std::string& rText, char Exponent = 'E') {
    fortran_exponent_to_e(rText.data(), rText.size(), Exponent);
}

}  // namespace meshioplusplus::detail
