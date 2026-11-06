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
 * @file detail/utf8.hpp
 * @brief Append one code point to a string as UTF-8.
 *
 * A **core-private** header: no installed header names it, and it adds
 * nothing to the API or the ABI.
 */

// System includes
#include <cstdint>
#include <string>

namespace meshioplusplus {
namespace detail {

/**
 * @brief Append the UTF-8 encoding of a code point (1 to 4 bytes).
 * @param rOut String to append to.
 * @param Cp The code point; no surrogate or range check is made.
 */
inline void utf8_append(std::string& rOut, std::uint32_t Cp) {
    if (Cp < 0x80) {
        rOut += static_cast<char>(Cp);
    } else if (Cp < 0x800) {
        rOut += static_cast<char>(0xC0 | (Cp >> 6));
        rOut += static_cast<char>(0x80 | (Cp & 0x3F));
    } else if (Cp < 0x10000) {
        rOut += static_cast<char>(0xE0 | (Cp >> 12));
        rOut += static_cast<char>(0x80 | ((Cp >> 6) & 0x3F));
        rOut += static_cast<char>(0x80 | (Cp & 0x3F));
    } else {
        rOut += static_cast<char>(0xF0 | (Cp >> 18));
        rOut += static_cast<char>(0x80 | ((Cp >> 12) & 0x3F));
        rOut += static_cast<char>(0x80 | ((Cp >> 6) & 0x3F));
        rOut += static_cast<char>(0x80 | (Cp & 0x3F));
    }
}

}  // namespace detail
}  // namespace meshioplusplus
