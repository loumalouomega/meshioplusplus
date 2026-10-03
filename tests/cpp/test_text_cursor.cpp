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
 * @file test_text_cursor.cpp
 * @brief The private text tokenizer (`src/cpp/src/detail/text_cursor.hpp`)
 *        against what it replaces: `std::getline`, `istringstream >>`, and
 *        the `strtoll`/`strtoull`/`parse_double` token parses.
 */

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <type_traits>
#include <cmath>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/text_cursor.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"

namespace {

namespace det = meshioplusplus::detail;

static_assert(std::is_constructible_v<det::TextStream, std::string&&>);
static_assert(!std::is_constructible_v<det::TextStream, const std::string&&>);
static_assert(std::is_constructible_v<det::TextStream, std::string&>);
static_assert(std::is_constructible_v<det::TextStream, const std::string&>);
static_assert(std::is_constructible_v<det::TextStream, std::string_view>);
static_assert(std::is_constructible_v<det::TextStream, const char*>);
static_assert(!std::is_constructible_v<det::TextCursor, std::string&&>);
static_assert(std::is_constructible_v<det::TextCursor, std::string_view>);

TEST(TextCursor, StreamOwnsShortAndLongTemporaryStrings) {
    for (const std::size_t padding : {0u, 200u}) {
        det::TextStream stream(std::string(padding, ' ') + "1.25 -7 a_long_owned_token");
        // Both SSO and heap-backed storage survive destruction of the input
        // temporary and unrelated allocation/stack activity before extraction.
        const std::vector<std::string> churn(100, std::string(200, 'x'));
        double real = 0;
        int integer = 0;
        std::string token;
        stream >> real >> integer >> token;
        ASSERT_TRUE(stream);
        EXPECT_DOUBLE_EQ(real, 1.25);
        EXPECT_EQ(integer, -7);
        EXPECT_EQ(token, "a_long_owned_token");
        EXPECT_EQ(churn.front().size(), 200u);
    }
}

// Random lines over an alphabet that makes every kind of edge case likely.
std::vector<std::string> tc_lines(std::size_t count) {
    static const char kAlphabet[] = "0123456789+-.eExXinfaINFA \t,/";
    std::vector<std::string> out = {"",
                                    " ",
                                    "1",
                                    "+1",
                                    "-",
                                    ".",
                                    "1e",
                                    "1e+",
                                    "0x1p3",
                                    "inf",
                                    "nan",
                                    "1e999",
                                    "-1e999",
                                    "1e-999",
                                    "9223372036854775807",
                                    "9223372036854775808",
                                    "-9223372036854775808",
                                    "-9223372036854775809",
                                    "2147483648",
                                    "-2147483649",
                                    "1.5abc 2 3",
                                    "12 34 56",
                                    "  7\t8  ",
                                    "1..2",
                                    "+-5",
                                    "00012",
                                    ".5e-3"};
    std::uint64_t s = 1234567;
    for (std::size_t i = 0; i < count; ++i) {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        const std::size_t len = (s >> 58) % 24;
        std::string line;
        for (std::size_t k = 0; k < len; ++k) {
            s = s * 6364136223846793005ull + 1442695040888963407ull;
            line += kAlphabet[(s >> 33) % (sizeof kAlphabet - 1)];
        }
        out.push_back(line);
    }
    return out;
}

template <class T>
void tc_expect_same_sequence(const std::string& rLine) {
    auto ref = meshioplusplus::detail::make_classic_istringstream(rLine);
    det::TextStream got(rLine);
    for (int k = 0; k < 6; ++k) {
        T a{static_cast<T>(42)}, b{static_cast<T>(42)};
        ref >> a;
        got >> b;
        ASSERT_EQ(static_cast<bool>(ref), static_cast<bool>(got)) << "'" << rLine << "' #" << k;
        if constexpr (std::is_floating_point_v<T>)
            ASSERT_TRUE(a == b || (a != a && b != b))
                << "'" << rLine << "' #" << k << " " << a << " vs " << b;
        else
            ASSERT_EQ(a, b) << "'" << rLine << "' #" << k;
    }
}

}  // namespace

TEST(TextCursor, StreamExtractionMatchesIstringstream) {
    for (const std::string& line : tc_lines(20000)) {
        tc_expect_same_sequence<double>(line);
        tc_expect_same_sequence<long long>(line);
        tc_expect_same_sequence<long>(line);
        tc_expect_same_sequence<int>(line);
        tc_expect_same_sequence<unsigned long>(line);
        tc_expect_same_sequence<unsigned long long>(line);
        tc_expect_same_sequence<unsigned>(line);
        // Mixed: a word, then numbers, as `tag x y z` lines read.
        auto ref = meshioplusplus::detail::make_classic_istringstream(line);
        det::TextStream got(line);
        std::string w1, w2;
        double d1 = 42, d2 = 42;
        ref >> w1 >> d1;
        got >> w2 >> d2;
        ASSERT_EQ(w1, w2) << line;
        ASSERT_TRUE(d1 == d2 || (d1 != d1 && d2 != d2)) << line;
        ASSERT_EQ(static_cast<bool>(ref), static_cast<bool>(got)) << line;
    }
}

TEST(TextCursor, BoundedPrefixCursorMatchesTheCLibrary) {
    auto inputs = tc_lines(2000);
    inputs.insert(inputs.end(), {"2.0 1e3", "1.5-2.5", "0x1p2suffix", "nan(payload)",
                                 std::string("1\0tail", 6), std::string(1000, '0') + "1"});
    for (const auto& text : inputs) {
        const auto storage = text + "999";
        det::TextCursor cursor(std::string_view(storage.data(), text.size()));
        const char* end = nullptr;
        const auto expected = det::parse_double(text.c_str(), end);
        double actual = 0.0;
        const bool converted = cursor.DoublePrefix(actual);
        ASSERT_EQ(converted, end != text.c_str()) << text;
        if (converted) {
            EXPECT_EQ(cursor.Pos(), static_cast<std::size_t>(end - text.c_str())) << text;
            if (std::isnan(expected))
                EXPECT_TRUE(std::isnan(actual));
            else {
                EXPECT_EQ(expected, actual) << text;
                EXPECT_EQ(std::signbit(expected), std::signbit(actual)) << text;
            }
        } else {
            EXPECT_EQ(cursor.Pos(), 0u);
        }
        cursor.Seek(0);
        char* integer_end = nullptr;
        const auto integer = std::strtoll(text.c_str(), &integer_end, 10);
        std::int64_t actual_integer = 0;
        ASSERT_EQ(cursor.IntPrefix(actual_integer), integer_end != text.c_str()) << text;
        if (integer_end != text.c_str()) {
            EXPECT_EQ(integer, actual_integer) << text;
            EXPECT_EQ(cursor.Pos(), static_cast<std::size_t>(integer_end - text.c_str())) << text;
        }
    }
    det::TextCursor adjacent(" 2.0-3.0 1e3");
    for (const auto expected : {2.0, -3.0, 1000.0}) {
        double actual = 0.0;
        ASSERT_TRUE(adjacent.DoublePrefix(actual));
        EXPECT_DOUBLE_EQ(expected, actual);
    }
    EXPECT_TRUE(adjacent.AtEnd());
}

TEST(TextCursor, SharedLineAndRecordPositioning) {
    det::TextCursor cursor("first\r\n\nlast");
    EXPECT_EQ(cursor.Line(true), "first");
    const auto saved = cursor.Pos();
    EXPECT_EQ(cursor.Line(true), "");
    EXPECT_EQ(cursor.Line(true), "last");
    EXPECT_TRUE(cursor.AtEnd());
    EXPECT_EQ(cursor.Line(), "");
    cursor.Seek(saved);
    EXPECT_EQ(cursor.Line(), "");
    EXPECT_EQ(cursor.Line(), "last");
    const std::vector<std::string_view> records = {"one", "", "three"};
    det::RecordCursor<std::string_view> rec(records);
    EXPECT_EQ(rec.Remaining(), 3u);
    EXPECT_EQ(rec.Peek(2), "three");
    EXPECT_EQ(rec.Next(), "one");
    EXPECT_EQ(rec.Pos(), 1u);
    EXPECT_EQ(rec.Next(), "");
    EXPECT_EQ(rec.Next(), "three");
    EXPECT_TRUE(rec.Done());
    EXPECT_EQ(rec.Remaining(), 0u);
}

TEST(TextCursor, GetlineMatchesTheStream) {
    for (const std::string text :
         {std::string(""), std::string("a"), std::string("a\n"), std::string("a\nb"),
          std::string("\n\n"), std::string("a,b,,c"), std::string("x\r\ny\r\n")})
        for (const char delim : {'\n', ','}) {
            std::istringstream ref(text);
            det::TextStream got(text);
            for (int k = 0; k < 6; ++k) {
                std::string a = "?", b = "?";
                std::getline(ref, a, delim);
                getline(got, b, delim);
                ASSERT_EQ(static_cast<bool>(ref), static_cast<bool>(got)) << text << " #" << k;
                if (ref)
                    ASSERT_EQ(a, b) << text << " #" << k;
            }
        }
}

TEST(TextCursor, SplittingMatchesGetlineAndStreams) {
    for (const std::string text :
         {std::string(""), std::string("a"), std::string("a\n"), std::string("a\nb"),
          std::string("\n\n"), std::string("x\r\ny\r\n")}) {
        std::istringstream in(text);
        std::vector<std::string> want;
        for (std::string l; std::getline(in, l);)
            want.push_back(l);
        const std::vector<std::string_view> got = det::split_lines(text);
        ASSERT_EQ(got.size(), want.size()) << text;
        for (std::size_t i = 0; i < got.size(); ++i)
            EXPECT_EQ(std::string(got[i]), want[i]);
    }
    for (const std::string& line : tc_lines(2000)) {
        std::istringstream in(line);
        std::vector<std::string> want;
        for (std::string t; in >> t;)
            want.push_back(t);
        const std::vector<std::string_view> got = det::split_blanks(line);
        ASSERT_EQ(got.size(), want.size()) << line;
        for (std::size_t i = 0; i < got.size(); ++i)
            EXPECT_EQ(std::string(got[i]), want[i]);
    }
}

TEST(TextCursor, TokenParsesMatchTheCLibrary) {
    for (const std::string& tok : tc_lines(20000)) {
        EXPECT_EQ(det::strtoll_token(tok), std::strtoll(tok.c_str(), nullptr, 10)) << tok;
        EXPECT_EQ(det::strtoull_token(tok), std::strtoull(tok.c_str(), nullptr, 10)) << tok;
        const double a = det::parse_double_prefix(tok);
        const double b = det::parse_double(tok);
        EXPECT_TRUE(a == b || (a != a && b != b)) << tok;
        char* end = nullptr;
        const long long ll = tok.empty() ? 0 : std::strtoll(tok.c_str(), &end, 10);
        std::int64_t v = 0;
        const bool whole = !tok.empty() && end == tok.c_str() + tok.size();
        EXPECT_EQ(det::parse_int_token(tok, v), whole) << tok;
        if (whole)
            EXPECT_EQ(v, ll) << tok;
    }
}
