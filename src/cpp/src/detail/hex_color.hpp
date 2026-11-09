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
 * @file detail/hex_color.hpp
 * @brief `#rrggbb` to an 8-bit RGB triple, the colour syntax every writer that
 * takes a `nan_color` accepts (`gltf`, `colorize`).
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 */

// System includes
#include <cstdint>
#include <stdexcept>
#include <string>

// Project includes
#include "meshioplusplus/detail/colormap.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Parse `#rrggbb` (either case).
 * @param rPrefix prepended to the error message, which names the option.
 * @throws std::invalid_argument on anything else.
 */
inline Rgb parse_hex_color(const std::string& rText, const std::string& rPrefix) {
    auto hex = [](char c) -> int {
        if (c >= '0' && c <= '9')
            return c - '0';
        if (c >= 'a' && c <= 'f')
            return c - 'a' + 10;
        if (c >= 'A' && c <= 'F')
            return c - 'A' + 10;
        return -1;
    };
    bool ok = rText.size() == 7 && rText[0] == '#';
    for (std::size_t i = 1; ok && i < 7; ++i)
        ok = hex(rText[i]) >= 0;
    if (!ok)
        throw std::invalid_argument(rPrefix + "nan_color '" + rText +
                                    "' is not a '#rrggbb' colour");
    auto byte = [&](std::size_t i) {
        return static_cast<std::uint8_t>(hex(rText[i]) * 16 + hex(rText[i + 1]));
    };
    return Rgb{byte(1), byte(3), byte(5)};
}

}  // namespace detail
}  // namespace meshioplusplus
