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
/**
 * @file test_format_infra.cpp
 * @brief The shared reader helpers of v16.7.0: the byte cursor, Fortran record
 *        framing, the degenerate-brick collapse and the Abaqus element types.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <array>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/abaqus_types.hpp"
#include "meshioplusplus/detail/binary_stream.hpp"
#include "meshioplusplus/detail/degenerate_solid.hpp"
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace {

namespace detail = meshioplusplus::detail;
using meshioplusplus::ReadError;

void put_u32(std::string& rOut, std::uint32_t Value, bool Big) {
    for (int b = 0; b < 4; ++b) {
        const int shift = Big ? 8 * (3 - b) : 8 * b;
        rOut += static_cast<char>((Value >> shift) & 0xff);
    }
}

std::string record(const std::string& rPayload, bool Big) {
    std::string out;
    put_u32(out, static_cast<std::uint32_t>(rPayload.size()), Big);
    out += rPayload;
    put_u32(out, static_cast<std::uint32_t>(rPayload.size()), Big);
    return out;
}

}  // namespace

TEST(ByteCursor, ReadsBothByteOrdersAndRefusesToOverrun) {
    const std::string be("\x00\x00\x00\x2a\x3f\xf0\x00\x00\x00\x00\x00\x00", 12);
    detail::ByteCursor big(be.data(), be.size(), true, "test");
    EXPECT_EQ(big.U32(), 42u);
    EXPECT_EQ(big.F64(), 1.0);
    EXPECT_TRUE(big.AtEnd());
    EXPECT_THROW(big.U32(), ReadError);
    const std::string le("\x2a\x00\x00\x00", 4);
    detail::ByteCursor little(le.data(), le.size(), false, "test");
    EXPECT_EQ(little.I32(), 42);
}

TEST(FortranRecords, SniffsTheFramingAndSplitsTheRecords) {
    for (const bool big : {false, true}) {
        const std::string file = record("abcd", big) + record(std::string(12, 'x'), big);
        const auto layout = detail::sniff_fortran_records(file.data(), file.size());
        ASSERT_TRUE(layout.has_value());
        EXPECT_EQ(layout->mMarkerBytes, 4);
        EXPECT_EQ(layout->mBigEndian, big);
        const auto recs = detail::fortran_records(file.data(), file.size(), *layout, "test");
        ASSERT_EQ(recs.size(), 2u);
        EXPECT_EQ(file.substr(recs[0].mOffset, recs[0].mSize), "abcd");
        EXPECT_EQ(recs[1].mSize, 12u);
    }
    const std::string text = "not a fortran file";
    EXPECT_FALSE(detail::sniff_fortran_records(text.data(), text.size()).has_value());
    std::string bad = record("abcd", false);
    bad[bad.size() - 4] = 9;  // trailing marker disagrees
    const detail::FortranRecordLayout layout{4, false};
    EXPECT_THROW(detail::fortran_records(bad.data(), bad.size(), layout, "test"), ReadError);
}

namespace {

// A length marker of `Width` bytes, written independently of the library.
void put_marker(std::string& rOut, std::uint64_t Value, int Width, bool Big) {
    for (int b = 0; b < Width; ++b) {
        const int shift = Big ? 8 * (Width - 1 - b) : 8 * b;
        rOut += static_cast<char>((Value >> shift) & 0xff);
    }
}

}  // namespace

TEST(FortranRecords, AppendRecordFramesInEveryLayoutAndReadsBack) {
    for (const int width : {4, 8}) {
        for (const bool big : {false, true}) {
            const detail::FortranRecordLayout layout{width, big};
            std::string expected, built;
            std::vector<std::string> payloads{"abcd", "", std::string(13, 'x')};
            for (const std::string& p : payloads) {
                put_marker(expected, p.size(), width, big);
                expected += p;
                put_marker(expected, p.size(), width, big);
                detail::fortran_append_record(built, p.data(), p.size(), layout);
            }
            EXPECT_EQ(built, expected) << width << (big ? " big" : " little");
            const auto sniffed = detail::sniff_fortran_records(built.data(), built.size());
            ASSERT_TRUE(sniffed.has_value());
            EXPECT_EQ(sniffed->mMarkerBytes, width);
            EXPECT_EQ(sniffed->mBigEndian, big);
            const auto recs = detail::fortran_records(built.data(), built.size(), layout, "test");
            ASSERT_EQ(recs.size(), payloads.size());
            for (std::size_t i = 0; i < recs.size(); ++i)
                EXPECT_EQ(built.substr(recs[i].mOffset, recs[i].mSize), payloads[i]);

            // A cut-short tail: an error, unless the reader is lenient.
            const std::string cut = built.substr(0, built.size() - 3);
            EXPECT_THROW(detail::fortran_records(cut.data(), cut.size(), layout, "test"),
                         ReadError);
            EXPECT_EQ(detail::fortran_records(cut.data(), cut.size(), layout, "test", true).size(),
                      2u);
            const std::string runs_past = built.substr(0, built.size() - width - 8);
            EXPECT_THROW(
                detail::fortran_records(runs_past.data(), runs_past.size(), layout, "test"),
                ReadError);
            std::string bad = built;
            bad.back() = 0x7f;  // the last trailing marker disagrees
            EXPECT_THROW(detail::fortran_records(bad.data(), bad.size(), layout, "test"),
                         ReadError);
            EXPECT_EQ(detail::fortran_records(bad.data(), bad.size(), layout, "test", true).size(),
                      2u);
        }
    }
}

TEST(FortranRecords, PutHelpersWriteBothByteOrders) {
    std::string be, le;
    detail::fortran_put_int32(be, 0x01020304, true);
    detail::fortran_put_int32(le, 0x01020304, false);
    EXPECT_EQ(be, std::string("\x01\x02\x03\x04", 4));
    EXPECT_EQ(le, std::string("\x04\x03\x02\x01", 4));
    be.clear();
    detail::fortran_put_float32(be, 1.0F, true);
    EXPECT_EQ(be, std::string("\x3f\x80\x00\x00", 4));
    be.clear();
    detail::fortran_put_float64(be, 1.0, true);
    EXPECT_EQ(be, std::string("\x3f\xf0\x00\x00\x00\x00\x00\x00", 8));
    le.clear();
    detail::fortran_put_int64(le, -2, false);
    EXPECT_EQ(le, std::string("\xfe\xff\xff\xff\xff\xff\xff\xff", 8));
    std::string out;
    EXPECT_NO_THROW(detail::fortran_append_record(out, nullptr, 0, {4, true}));
    EXPECT_EQ(out, std::string(8, '\0'));
}

TEST(DegenerateSolid, CollapsesEveryDocumentedPattern) {
    using A = std::array<std::int64_t, 8>;
    EXPECT_STREQ(detail::collapse_brick(A{1, 2, 3, 4, 4, 4, 4, 4}).mType, "tetra");
    EXPECT_STREQ(detail::collapse_brick(A{1, 2, 3, 3, 4, 4, 4, 4}).mType, "tetra");
    EXPECT_STREQ(detail::collapse_brick(A{1, 2, 3, 4, 5, 5, 5, 5}).mType, "pyramid");
    EXPECT_STREQ(detail::collapse_brick(A{1, 2, 3, 4, 5, 6, 7, 8}).mType, "hexahedron");
    const auto lsdyna = detail::collapse_brick(A{1, 2, 3, 3, 5, 6, 7, 7});
    EXPECT_STREQ(lsdyna.mType, "wedge");
    EXPECT_EQ(lsdyna.mNodes, (std::vector<std::int64_t>{1, 2, 3, 5, 6, 7}));
    // Radioss: the side edge 4-1 collapsed in both faces.
    const auto radioss = detail::collapse_brick(A{1, 2, 3, 1, 5, 6, 7, 5});
    EXPECT_STREQ(radioss.mType, "wedge");
    EXPECT_EQ(radioss.mNodes, (std::vector<std::int64_t>{1, 2, 3, 5, 6, 7}));
}

TEST(AbaqusTypes, TableThenFamilyAndNodeCount) {
    EXPECT_EQ(detail::abaqus_cell_type("C3D4H", 4), "tetra");  // was the invalid "tetra4"
    EXPECT_EQ(detail::abaqus_cell_type("C3D20R", 20), "hexahedron20");
    EXPECT_EQ(detail::abaqus_cell_type("CPE8R", 8), "quad8");
    EXPECT_EQ(detail::abaqus_cell_type("DC3D10", 10), "tetra10");
    EXPECT_EQ(detail::abaqus_cell_type("SC8R", 8), "hexahedron");
    EXPECT_EQ(detail::abaqus_cell_type("S4R", 4), "quad");
    EXPECT_EQ(detail::abaqus_cell_type("M3D9", 9), "quad9");
    EXPECT_EQ(detail::abaqus_cell_type("B32", 3), "line3");
    EXPECT_EQ(detail::abaqus_cell_type("SAX1", 2), "line");
    EXPECT_EQ(detail::abaqus_cell_type("U1", 4), "");
    EXPECT_EQ(detail::abaqus_cell_type("C3D8", 7), "");
    // The writer's inverse still names C3D4 for a tetra.
    std::string last;
    for (const auto& [abq, type] : detail::abaqus_type_table())
        if (type == "tetra")
            last = abq;
    EXPECT_EQ(last, "C3D4");
}
