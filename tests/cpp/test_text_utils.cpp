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
// Unit tests for detail/text_utils.hpp (roadmap 6.1): the trim the readers each
// wrote out, with the whitespace sets they used.

// System includes
#include <string>
#include <string_view>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/text_utils.hpp"

using meshioplusplus::detail::TextUtils;

TEST(TextUtils, TrimViewStripsTheSetAtBothEnds) {
    EXPECT_EQ(TextUtils::TrimView("  \t a b \r\n"), "a b");
    EXPECT_EQ(TextUtils::TrimView(""), "");
    EXPECT_EQ(TextUtils::TrimView(" \t\r\n"), "");
    EXPECT_EQ(TextUtils::TrimView("abc"), "abc");
}

TEST(TextUtils, EachSetKeepsWhatItDoesNotName) {
    EXPECT_EQ(TextUtils::TrimView(" a \r\n", TextUtils::kSpaceTab), "a \r\n");
    EXPECT_EQ(TextUtils::TrimView(" a \r\n", TextUtils::kSpaceTabCr), "a \r\n");
    EXPECT_EQ(TextUtils::TrimView(" a\r \n", TextUtils::kSpaceTabCr), "a\r \n");
    EXPECT_EQ(TextUtils::TrimView("\v\fa\f\v", TextUtils::kBlanks), "\v\fa\f\v");
    EXPECT_EQ(TextUtils::TrimView("\v\fa\f\v", TextUtils::kClassicSpace), "a");
}

TEST(TextUtils, ClassicSpaceIsStdIsspaceInTheCLocale) {
    for (int c = 0; c < 256; ++c) {
        const bool in_set = TextUtils::kClassicSpace.find(static_cast<char>(c)) !=
                            std::string_view::npos && c != 0;
        const bool is_space = c == ' ' || (c >= '\t' && c <= '\r');
        EXPECT_EQ(in_set, is_space) << c;
    }
}

TEST(TextUtils, TrimReturnsAnOwnedCopy) {
    std::string s = "  x  ";
    const std::string t = TextUtils::Trim(s);
    s = "changed";
    EXPECT_EQ(t, "x");
}
