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

TEST(TextUtils, OneSidedTrims) {
    EXPECT_EQ(TextUtils::LeftTrimView("  a b  "), "a b  ");
    EXPECT_EQ(TextUtils::RightTrimView("  a b  "), "  a b");
    EXPECT_EQ(TextUtils::LeftTrimView("   "), "");
    EXPECT_EQ(TextUtils::RightTrimView(""), "");
    EXPECT_EQ(TextUtils::RightTrimView("a \r", TextUtils::kSpaceTab), "a \r");
}

TEST(TextUtils, CaseFoldIsAsciiOnly) {
    EXPECT_EQ(TextUtils::Upper("Node-1 abz"), "NODE-1 ABZ");
    EXPECT_EQ(TextUtils::Lower("NODE-1 ABZ"), "node-1 abz");
    EXPECT_EQ(TextUtils::Lower(""), "");
    // A byte above 0x7F (Latin-1 "A grave" and "a grave") is left alone.
    EXPECT_EQ(TextUtils::Upper("\xE0"), "\xE0");
    EXPECT_EQ(TextUtils::Lower("\xC0"), "\xC0");
    for (int c = 0; c < 256; ++c) {
        const char ch = static_cast<char>(c);
        const bool upper = c >= 'A' && c <= 'Z', lower = c >= 'a' && c <= 'z';
        EXPECT_EQ(TextUtils::LowerChar(ch), upper ? static_cast<char>(c + 32) : ch);
        EXPECT_EQ(TextUtils::UpperChar(ch), lower ? static_cast<char>(c - 32) : ch);
    }
}

TEST(TextUtils, EndsWithLoweredComparesTheTailWithoutCase) {
    EXPECT_TRUE(TextUtils::EndsWithLowered("model.GLB", ".glb"));
    EXPECT_TRUE(TextUtils::EndsWithLowered("model.glb", ".glb"));
    EXPECT_FALSE(TextUtils::EndsWithLowered("model.glt", ".glb"));
    EXPECT_FALSE(TextUtils::EndsWithLowered("b", ".glb"));
    EXPECT_TRUE(TextUtils::EndsWithLowered(".glb", ".glb"));
}
