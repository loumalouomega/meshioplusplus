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
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"

using meshioplusplus::detail::is_c_decimal_point;
using meshioplusplus::detail::parse_double;
using meshioplusplus::detail::snprintf_c;
using meshioplusplus::detail::stod_c;

namespace {

// RAII guard: saves/restores LC_NUMERIC around a test, mandatory since every
// gtest in this binary shares one process -- a leaked comma locale would
// corrupt every later test's floating-point parsing/formatting.
class LcNumericGuard {
public:
    LcNumericGuard() {
        const char* cur = std::setlocale(LC_NUMERIC, nullptr);
        mSaved = cur ? cur : "C";
    }
    ~LcNumericGuard() { std::setlocale(LC_NUMERIC, mSaved.c_str()); }
    LcNumericGuard(const LcNumericGuard&) = delete;
    LcNumericGuard& operator=(const LcNumericGuard&) = delete;

private:
    std::string mSaved;
};

// The first candidate whose decimal_point is a comma, or empty if none of
// these locales are installed on this machine.
std::string find_comma_locale() {
    static const char* kCandidates[] = {
        "de_DE.UTF-8", "de_DE.utf8",  "fr_FR.UTF-8", "fr_FR.utf8",          "es_ES.UTF-8",
        "es_ES.utf8",  "nl_NL.UTF-8", "nl_NL.utf8",  "German_Germany.1252", "French_France.1252",
    };
    for (const char* cand : kCandidates) {
        if (std::setlocale(LC_NUMERIC, cand) != nullptr) {
            if (is_c_decimal_point(std::localeconv()->decimal_point))
                continue;  // some libc builds accept the name but keep "C"
            return cand;
        }
    }
    return {};
}

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
