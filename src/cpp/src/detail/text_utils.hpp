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

#pragma once

/**
 * @file detail/text_utils.hpp
 * @brief `TextUtils`: the trim the text readers each wrote out under their own
 * name, with the whitespace sets they used (roadmap 6.1).
 *
 * A reader strips a different set of characters (a keyword deck blanks and
 * line ends, a column format blanks only), so the set is an argument and the
 * three that occur are named. Each member does exactly what the copies did and
 * never calls `<cctype>`, so the result does not depend on the global locale.
 * A **core-private** header (not installed).
 */

// System includes
#include <cstddef>
#include <string>
#include <string_view>

// Project includes
#include "always_inline.hpp"

namespace meshioplusplus::detail {

struct TextUtils {
    /// Space, tab, carriage return and newline: what the keyword decks strip.
    static constexpr std::string_view kBlanks = " \t\r\n";
    /// Space and tab only: the fixed-column formats.
    static constexpr std::string_view kSpaceTab = " \t";
    /// Space, tab and carriage return, not the newline.
    static constexpr std::string_view kSpaceTabCr = " \t\r";
    /// `std::isspace` in the classic locale: also form feed and vertical tab.
    static constexpr std::string_view kClassicSpace = " \t\n\v\f\r";

    /// `Text` without the characters of `Set` at either end.
    static MIO_ALWAYS_INLINE std::string_view TrimView(std::string_view Text,
                                                       std::string_view Set = kBlanks) {
        const std::size_t b = Text.find_first_not_of(Set);
        if (b == std::string_view::npos)
            return {};
        return Text.substr(b, Text.find_last_not_of(Set) - b + 1);
    }

    /// `TrimView` as an owned string.
    static MIO_ALWAYS_INLINE std::string Trim(std::string_view Text,
                                              std::string_view Set = kBlanks) {
        return std::string(TrimView(Text, Set));
    }
};

}  // namespace meshioplusplus::detail
