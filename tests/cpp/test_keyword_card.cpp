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
// Tests for detail/keyword_card's Fortran format-line parsing (the layout lines
// of ANSYS `.cdb` blocks).

// External includes
#include <gtest/gtest.h>

// System includes
#include <string>
#include <bit>
#include <cstdint>
#include <cmath>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "../../src/cpp/src/detail/keyword_card_view.hpp"

using meshioplusplus::ReadError;
using meshioplusplus::detail::CardField;
using meshioplusplus::detail::parse_fortran_format;
using meshioplusplus::detail::split_fixed;

namespace {

std::string describe(const std::vector<CardField>& rFields) {
    std::string out;
    for (const CardField& f : rFields)
        out += std::string(1, f.mKind) + std::to_string(f.mWidth) + " ";
    return out;
}

}  // namespace

TEST(FortranFormat, ExpandsTheLinesAnsysWrites) {
    EXPECT_EQ(describe(parse_fortran_format("(3i8,6e20.13)")), "i8 i8 i8 r20 r20 r20 r20 r20 r20 ");
    EXPECT_EQ(describe(parse_fortran_format("(3i9,6e21.13e3)")),
              "i9 i9 i9 r21 r21 r21 r21 r21 r21 ");
    EXPECT_EQ(describe(parse_fortran_format("(1i7,2i9,6e21.13)")),
              "i7 i9 i9 r21 r21 r21 r21 r21 r21 ");
    EXPECT_EQ(parse_fortran_format("(19i10)").size(), 19u);
    const auto etblock = parse_fortran_format("(2i9,19a9)");
    EXPECT_EQ(etblock.size(), 21u);
    EXPECT_EQ(etblock[2].mKind, 'a');
    EXPECT_EQ(describe(parse_fortran_format("(2i8,6g16.9)")), "i8 i8 r16 r16 r16 r16 r16 r16 ");
}

TEST(KeywordCardViews, TokenizersMatchTheOwningApi) {
    namespace det = meshioplusplus::detail;
    const std::vector<CardField> layout = {{'i', 8}, {'r', 16}, {'x', 2}, {'a', 8}};
    const std::vector<std::string> lines = {"",
                                            " ",
                                            "1",
                                            " 1, 2.5D+1, ,name,",
                                            "1,2,3,4,5",
                                            "       1       2.500D+1  name    \r\n",
                                            std::string("1,2\0tail,3", 10)};
    for (const auto& line : lines) {
        for (const auto mode : {det::CardMode::Standard, det::CardMode::Long, det::CardMode::I10}) {
            const auto owned = det::split_card(line, layout, mode);
            const auto views = det::split_card_view(line, layout, mode);
            ASSERT_EQ(owned.size(), views.size());
            for (std::size_t i = 0; i < owned.size(); ++i)
                EXPECT_EQ(owned[i], views[i]);
        }
        const auto owned = det::split_fixed(line, layout);
        const auto views = det::split_fixed_view(line, layout);
        ASSERT_EQ(owned.size(), views.size());
        for (std::size_t i = 0; i < owned.size(); ++i)
            EXPECT_EQ(owned[i], views[i]);
    }
}

TEST(KeywordCardViews, NumbersAndErrorsMatchTheOwningApi) {
    namespace det = meshioplusplus::detail;
    std::vector<std::string> tokens = {"",
                                       " ",
                                       "1",
                                       "+1",
                                       "-0",
                                       "-1",
                                       "1suffix",
                                       "1 ",
                                       " 1",
                                       "0x1p2",
                                       ".5-3",
                                       "-1.5+3",
                                       "1D+2",
                                       "1d-2",
                                       "1e999",
                                       "1e-999",
                                       "inf",
                                       "nan",
                                       "1e",
                                       "+",
                                       ".",
                                       "9223372036854775807",
                                       "9223372036854775808",
                                       "-9223372036854775808",
                                       "-9223372036854775809",
                                       std::string("1\0tail", 6),
                                       std::string("\0tail", 5)};
    for (const auto size : {62u, 63u, 64u, 65u, 93u, 94u, 95u, 96u, 97u, 1000u}) {
        tokens.push_back(std::string(size, '0') + "1");
        tokens.push_back("1." + std::string(size, '0') + "-2");
        tokens.push_back(std::string(size, '9'));
    }
    for (const auto& text : tokens) {
        // The view ends before sentinel bytes; no parser may read those bytes.
        const std::string storage = text + "999";
        const std::string_view view(storage.data(), text.size());
        SCOPED_TRACE(text);
        std::string int_error, real_error;
        std::int64_t integer = 0;
        double real = 0.0;
        try {
            integer = det::card_to_int(text, " at test", "test");
        } catch (const ReadError& exc) {
            int_error = exc.what();
        }
        try {
            const auto actual = det::card_to_int_view(view, " at test", "test");
            EXPECT_TRUE(int_error.empty());
            EXPECT_EQ(integer, actual);
        } catch (const ReadError& exc) {
            EXPECT_FALSE(int_error.empty());
            EXPECT_EQ(int_error, exc.what());
        }
        try {
            real = det::card_to_real(text, " at test", "test");
        } catch (const ReadError& exc) {
            real_error = exc.what();
        }
        try {
            const auto actual = det::card_to_real_view(view, " at test", "test");
            EXPECT_TRUE(real_error.empty());
            if (std::isnan(real))
                EXPECT_TRUE(std::isnan(actual));
            else
                EXPECT_EQ(std::bit_cast<std::uint64_t>(real), std::bit_cast<std::uint64_t>(actual));
        } catch (const ReadError& exc) {
            EXPECT_FALSE(real_error.empty());
            EXPECT_EQ(real_error, exc.what());
        }
    }
}

TEST(FortranFormat, GroupsScaleFactorsAndSkips) {
    EXPECT_EQ(describe(parse_fortran_format("(2(i8,e16.9))")), "i8 r16 i8 r16 ");
    EXPECT_EQ(describe(parse_fortran_format("(1P,3E20.12)")), "r20 r20 r20 ");
    EXPECT_EQ(describe(parse_fortran_format("( I5 , 2X , ES12.4 )")), "i5 x2 r12 ");
    EXPECT_EQ(describe(parse_fortran_format("i9")), "i9 ");
    for (const char* bad : {"", "()", "(3q8)", "(i)", "(2i8", "(i8))"})
        EXPECT_THROW(parse_fortran_format(bad), ReadError) << bad;
}

TEST(FortranFormat, SplitsTouchingFixedWidthFields) {
    const auto fields = parse_fortran_format("(3i8,6e20.13)");
    const auto parts = split_fixed(
        "    5609       0       0 3.9797161316330E+00 2.5147820926190E-01"
        "-5.1500799817626E-01\r",
        fields);
    EXPECT_EQ(parts, (std::vector<std::string>{"5609", "0", "0", "3.9797161316330E+00",
                                               "2.5147820926190E-01", "-5.1500799817626E-01"}));
    EXPECT_EQ(split_fixed("   12   3", parse_fortran_format("(i5,2x,i2)")),
              (std::vector<std::string>{"12", "3"}));
    EXPECT_TRUE(split_fixed("", fields).empty());
}
