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

    /// `Text` without the characters of `Set` at its start.
    static MIO_ALWAYS_INLINE std::string_view LeftTrimView(std::string_view Text,
                                                           std::string_view Set = kBlanks) {
        const std::size_t b = Text.find_first_not_of(Set);
        return b == std::string_view::npos ? std::string_view{} : Text.substr(b);
    }

    /// `Text` without the characters of `Set` at its end.
    static MIO_ALWAYS_INLINE std::string_view RightTrimView(std::string_view Text,
                                                            std::string_view Set = kBlanks) {
        const std::size_t e = Text.find_last_not_of(Set);
        return e == std::string_view::npos ? std::string_view{} : Text.substr(0, e + 1);
    }

    /// `C` with an ASCII capital turned into its small letter; no locale.
    static MIO_ALWAYS_INLINE char LowerChar(char C) {
        return C >= 'A' && C <= 'Z' ? static_cast<char>(C - 'A' + 'a') : C;
    }

    /// `C` with an ASCII small letter turned into its capital; no locale.
    static MIO_ALWAYS_INLINE char UpperChar(char C) {
        return C >= 'a' && C <= 'z' ? static_cast<char>(C - 'a' + 'A') : C;
    }

    /// `Text` with its ASCII capitals turned into small letters.
    static MIO_ALWAYS_INLINE std::string Lower(std::string_view Text) {
        std::string out(Text);
        for (char& c : out)
            c = LowerChar(c);
        return out;
    }

    /// `Text` with its ASCII small letters turned into capitals.
    static MIO_ALWAYS_INLINE std::string Upper(std::string_view Text) {
        std::string out(Text);
        for (char& c : out)
            c = UpperChar(c);
        return out;
    }

    /// Whether the end of `Text`, with its capitals lowered, is the lower-case
    /// `LowerSuffix` (a file extension compared without regard to case).
    static MIO_ALWAYS_INLINE bool EndsWithLowered(std::string_view Text,
                                                  std::string_view LowerSuffix) {
        if (Text.size() < LowerSuffix.size())
            return false;
        const std::string_view tail = Text.substr(Text.size() - LowerSuffix.size());
        for (std::size_t i = 0; i < tail.size(); ++i)
            if (LowerChar(tail[i]) != LowerSuffix[i])
                return false;
        return true;
    }

    /// `TrimView` as an owned string.
    static MIO_ALWAYS_INLINE std::string Trim(std::string_view Text,
                                              std::string_view Set = kBlanks) {
        return std::string(TrimView(Text, Set));
    }
};

}  // namespace meshioplusplus::detail
