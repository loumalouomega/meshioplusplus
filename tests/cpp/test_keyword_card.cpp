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
#include <vector>

// Project includes
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/exceptions.hpp"

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
