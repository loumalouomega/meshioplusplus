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

// External includes
#include <gtest/gtest.h>

// System includes
#include <cctype>
#include <clocale>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"
#include "comma_locale.hpp"

using meshioplusplus::detail::is_c_decimal_point;
using meshioplusplus::detail::parse_double;
using meshioplusplus::detail::snprintf_c;
using meshioplusplus::detail::stod_c;

namespace {

using mio_test::find_comma_locale;
using mio_test::LcNumericGuard;

}  // namespace

TEST(FastNumber, ParsesBasicLiteralsInTheCLocale) {
    LcNumericGuard guard;
    std::setlocale(LC_NUMERIC, "C");

    struct Case {
        const char* mText;
        double mExpected;
    };
    const std::vector<Case> cases = {
        {"1.5", 1.5},         {"+1.5", 1.5},    {" 1.5", 1.5},
        {"-42.125", -42.125}, {"0", 0.0},       {"-0.0", -0.0},
        {"1e10", 1e10},       {"1e-10", 1e-10}, {"1.7976931348623157e308", 1.7976931348623157e308},
    };
    for (const Case& c : cases) {
        const char* end = nullptr;
        double v = parse_double(c.mText, end);
        EXPECT_DOUBLE_EQ(v, c.mExpected) << c.mText;
        EXPECT_NE(end, c.mText) << c.mText;
    }
}

TEST(FastNumber, RejectsGarbageWithoutAdvancing) {
    LcNumericGuard guard;
    std::setlocale(LC_NUMERIC, "C");

    const char* text = "abc";
    const char* end = nullptr;
    double v = parse_double(text, end);
    EXPECT_EQ(v, 0.0);
    EXPECT_EQ(end, text);

    const char* empty = "";
    const char* end2 = nullptr;
    double v2 = parse_double(empty, end2);
    EXPECT_EQ(v2, 0.0);
    EXPECT_EQ(end2, empty);
}

TEST(FastNumber, HandlesHexFloatsAndInfNan) {
    LcNumericGuard guard;
    std::setlocale(LC_NUMERIC, "C");

    // Every tier must agree with plain strtod on these, since the fallback
    // chain falls through to it for anything from_chars declines.
    const char* cases[] = {"0x1p3", "inf", "-inf", "nan", "1e400", "1e-400"};
    for (const char* text : cases) {
        const char* end = nullptr;
        double got = parse_double(text, end);
        char* strtod_end = nullptr;
        double want = std::strtod(text, &strtod_end);
        if (std::isnan(want)) {
            EXPECT_TRUE(std::isnan(got)) << text;
        } else {
            EXPECT_DOUBLE_EQ(got, want) << text;
        }
        EXPECT_EQ(end - text, strtod_end - text) << text;
    }
}

TEST(FastNumber, EndsWhereStrtodEndsNextToADelimiter) {
    LcNumericGuard guard;
    std::setlocale(LC_NUMERIC, "C");

    // The from_chars tier bounds its input to the span a literal can cover;
    // the value and the end must still be exactly strtod's.
    const char* cases[] = {"1.5,2",      "2.0)", "1e5x",   "1e+",      "1e-7;",
                           "-.5e-3z",    ".5 ",  "5. ",    "-inf ",    "infinity2",
                           "nan(abc_1)", "nan(", "1a2a3a", "-0.0e0\n", "12345678901234567890"};
    for (const char* text : cases) {
        const char* end = nullptr;
        double got = parse_double(text, end);
        char* strtod_end = nullptr;
        double want = std::strtod(text, &strtod_end);
        if (std::isnan(want)) {
            EXPECT_TRUE(std::isnan(got)) << text;
        } else {
            EXPECT_EQ(std::signbit(got), std::signbit(want)) << text;
            EXPECT_EQ(got, want) << text;
        }
        EXPECT_EQ(end - text, strtod_end - text) << text;
    }
}

#ifdef MESHIOPLUSPLUS_HAS_FAST_FROM_CHARS
TEST(FastNumber, ReadsOnlyThePrefixItParsesNotTheRestOfTheBuffer) {
    LcNumericGuard guard;
    std::setlocale(LC_NUMERIC, "C");

    // A buffer with no NUL terminator: parsing its first number must read
    // no further than the delimiter after it. A parser that measured the
    // rest of the buffer (strlen) on every call -- once made every reader
    // parsing a whole file this way quadratic -- overruns this allocation,
    // which the sanitizer leg reports.
    const char text[] = {'1', '.', '5', ' ', '9', '9'};
    std::unique_ptr<char[]> buf(new char[sizeof(text)]);
    std::memcpy(buf.get(), text, sizeof(text));
    const char* end = nullptr;
    EXPECT_EQ(parse_double(buf.get(), end), 1.5);
    EXPECT_EQ(end, buf.get() + 3);
}
#endif

TEST(FastNumber, StringOverloadMatchesCStringOverload) {
    LcNumericGuard guard;
    std::setlocale(LC_NUMERIC, "C");
    EXPECT_DOUBLE_EQ(parse_double(std::string("3.25")), 3.25);
    EXPECT_DOUBLE_EQ(parse_double(std::string("not a number")), 0.0);
}

TEST(FastNumber, StodCMatchesStdStodOnValidInputAndThrowsOnInvalid) {
    // stod_c exists for call sites (the native CLI's option parsing) that
    // want std::stod's throw-on-failure contract rather than parse_double's
    // lenient 0.0, migrated by a plain name substitution.
    LcNumericGuard guard;
    std::setlocale(LC_NUMERIC, "C");
    EXPECT_DOUBLE_EQ(stod_c("1.5"), std::stod("1.5"));
    EXPECT_DOUBLE_EQ(stod_c("-42.125"), std::stod("-42.125"));
    EXPECT_THROW(stod_c("not a number"), std::invalid_argument);
    EXPECT_THROW(stod_c(""), std::invalid_argument);
}

TEST(FastNumber, SnprintfCMatchesPlainSnprintfInTheCLocale) {
    LcNumericGuard guard;
    std::setlocale(LC_NUMERIC, "C");

    char a[64];
    char b[64];
    int na = snprintf_c(a, sizeof(a), "%.17g", 1.5);
    int nb = std::snprintf(b, sizeof(b), "%.17g", 1.5);
    EXPECT_EQ(na, nb);
    EXPECT_STREQ(a, b);

    // Mixed int/float/width format, matching real call-site shapes.
    int nc = snprintf_c(a, sizeof(a), "%6d %25.16E\n", 42, 3.14159265358979);
    int nd = std::snprintf(b, sizeof(b), "%6d %25.16E\n", 42, 3.14159265358979);
    EXPECT_EQ(nc, nd);
    EXPECT_STREQ(a, b);

    // A pure-integer format has nothing to repair.
    int ne = snprintf_c(a, sizeof(a), "%6d\n", 7);
    int nf = std::snprintf(b, sizeof(b), "%6d\n", 7);
    EXPECT_EQ(ne, nf);
    EXPECT_STREQ(a, b);
}

TEST(FastNumber, UnderACommaLocaleParseIgnoresTheLocaleSeparator) {
    LcNumericGuard guard;
    const std::string loc = find_comma_locale();
    if (loc.empty())
        GTEST_SKIP() << "no comma-decimal locale available on this machine";

    ASSERT_EQ(std::string(std::localeconv()->decimal_point), ",");

    // Baseline: prove the locale actually changes plain strtod's behaviour,
    // so a no-op fast_number implementation could not pass this test
    // vacuously.
    EXPECT_EQ(std::strtod("1.5", nullptr), 1.0)
        << "the comma locale did not take effect for plain strtod";

    const char* end = nullptr;
    EXPECT_DOUBLE_EQ(parse_double("1.5", end), 1.5);
    // The locale's own comma separator must NOT be honoured as a decimal
    // point -- "1,5" parses as the integer prefix "1", exactly like strtod
    // under the C locale would.
    EXPECT_DOUBLE_EQ(parse_double("1,5", end), 1.0);
}

TEST(FastNumber, UnderACommaLocaleSnprintfCEmitsADot) {
    LcNumericGuard guard;
    const std::string loc = find_comma_locale();
    if (loc.empty())
        GTEST_SKIP() << "no comma-decimal locale available on this machine";

    // Baseline: plain snprintf really does emit a comma here.
    char baseline[64];
    std::snprintf(baseline, sizeof(baseline), "%.1f", 1.5);
    EXPECT_NE(std::strchr(baseline, ','), nullptr)
        << "the comma locale did not take effect for plain snprintf";

    char buf[64];
    int n = snprintf_c(buf, sizeof(buf), "%.1f", 1.5);
    EXPECT_GT(n, 0);
    EXPECT_STREQ(buf, "1.5");
    EXPECT_EQ(std::strchr(buf, ','), nullptr);
}

// Several numbers in one format, a literal comma beside a number, and a
// truncated buffer: before v19.1 only the first separator of the text was
// repaired, so the second number kept its comma and a literal ',' became '.'.
TEST(FastNumber, UnderACommaLocaleSnprintfCFormatsEveryNumber) {
    LcNumericGuard guard;
    if (find_comma_locale().empty())
        GTEST_SKIP() << "no comma-decimal locale available on this machine";
    char buf[96];
    int n = snprintf_c(buf, sizeof(buf), "%.1f %.1f %.17g\n", 1.5, -2.25, 0.1);
    EXPECT_EQ(n, static_cast<int>(std::strlen(buf)));
    EXPECT_STREQ(buf, "1.5 -2.2 0.10000000000000001\n");
    n = snprintf_c(buf, sizeof(buf), ",%.3e,%d,%.2f", 1.5, 7, 3.0);
    EXPECT_STREQ(buf, ",1.500e+00,7,3.00");
    char small[6];
    n = snprintf_c(small, sizeof(small), "%.4f", 2.5);
    EXPECT_EQ(n, 6);
    EXPECT_STREQ(small, "2.500");
}
