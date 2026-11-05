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
 * @file test_row_writer_to_chars.cpp
 * @brief `CNumber::Append`'s `to_chars` fast path (roadmap 5.7.1) against
 *        `snprintf` itself, byte for byte, over random bit patterns:
 *        subnormals, signed zeros, the largest and smallest exponents and
 *        the formats the writers use, prefixes and suffixes included.
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/row_writer.hpp"

namespace {

std::string rwtc_printf(const char* pFmt, double V) {
    char buf[64];
    const int n = std::snprintf(buf, sizeof buf, pFmt, V);
    return n <= 0 ? std::string() : std::string(buf, std::min<std::size_t>(n, sizeof buf - 1));
}

std::vector<double> rwtc_values() {
    std::vector<double> out = {0.0,
                               -0.0,
                               1.0,
                               -1.0,
                               0.1,
                               1e-300,
                               1e300,
                               123456789.0,
                               1.5e-7,
                               9.999999999999999e22,
                               std::numeric_limits<double>::min(),
                               std::numeric_limits<double>::denorm_min(),
                               -std::numeric_limits<double>::denorm_min(),
                               std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::lowest(),
                               std::numeric_limits<double>::infinity(),
                               -std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN(),
                               -std::numeric_limits<double>::quiet_NaN()};
    std::uint64_t s = 0x9e3779b97f4a7c15ull;
    for (int i = 0; i < 200000; ++i) {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        std::uint64_t bits = s;
        if (i % 4 == 1)
            bits &= 0x800FFFFFFFFFFFFFull;  // a subnormal (or a signed zero)
        else if (i % 4 == 2)
            bits = (bits & 0x800FFFFFFFFFFFFFull) | (std::uint64_t(1023 + (i % 40) - 20) << 52);
        double v;
        std::memcpy(&v, &bits, sizeof v);
        out.push_back(v);
        out.push_back(std::round(v * 1e3) / 1e3);  // short decimal expansions
    }
    return out;
}

}  // namespace

TEST(RowWriterToChars, AppendIsSnprintfByteForByte) {
    const meshioplusplus::detail::CNumber number;
    const char* formats[] = {"%.17g",  "%.16e",  "%.16g",   "%.11e",  "%.15g", " %.17g",
                             ",%.17g", " %.16e", "%.16e ",  "%.6g",   "%.0e",  "%.1g",
                             "%.20g",  "%g",     "%25.16E", "%+.16E", "%e"};
    const std::vector<double> values = rwtc_values();
    for (const char* fmt : formats) {
        std::size_t mismatches = 0;
        for (const double v : values) {
            std::string got;
            number.Append(got, fmt, v);
            const std::string want = rwtc_printf(fmt, v);
            if (got != want && mismatches++ < 5)
                ADD_FAILURE() << "format '" << fmt << "', value bits " << std::hexfloat << v
                              << ": got '" << got << "', want '" << want << "'";
        }
        EXPECT_EQ(mismatches, 0u) << "format '" << fmt << "'";
    }
}
