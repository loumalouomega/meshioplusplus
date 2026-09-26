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
 * @file test_row_writer.cpp
 * @brief The private ASCII row writer (`src/cpp/src/detail/row_writer.hpp`):
 *        `CNumber` against `snprintf_c`, and chunked output against a serial
 *        loop, byte for byte.
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/row_writer.hpp"
#include "meshioplusplus/detail/fast_number.hpp"

namespace {

namespace det = meshioplusplus::detail;

std::vector<double> rw_values() {
    std::vector<double> v = {0.0,
                             -0.0,
                             1.0,
                             -1.5,
                             0.1,
                             1e-300,
                             -1e300,
                             123456789.123456789,
                             std::numeric_limits<double>::min(),
                             std::numeric_limits<double>::denorm_min(),
                             std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()};
    std::uint64_t s = 42;
    for (int i = 0; i < 2000; ++i) {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        const double m = static_cast<double>(s >> 11) / 9007199254740992.0 - 0.5;
        v.push_back(m * std::pow(10.0, static_cast<int>((s >> 3) % 40) - 20));
    }
    return v;
}

}  // namespace

TEST(RowWriter, CNumberMatchesSnprintfC) {
    const det::CNumber num;
    for (const char* fmt : {"%.11e", "%.17g", "%.16e", "%.16g", "%21.13E", ", %.16e", "%.16e "})
        for (double v : rw_values()) {
            char buf[64];
            det::snprintf_c(buf, sizeof buf, fmt, v);
            std::string got;
            num.Append(got, fmt, v);
            ASSERT_EQ(got, buf) << fmt << " " << v;
            char buf2[64];
            num.Print(buf2, sizeof buf2, fmt, v);
            ASSERT_STREQ(buf2, buf) << fmt << " " << v;
        }
}

TEST(RowWriter, AppendIntMatchesTheStream) {
    for (std::int64_t v :
         {std::int64_t{0}, std::int64_t{-1}, std::int64_t{42},
          std::numeric_limits<std::int64_t>::min(), std::numeric_limits<std::int64_t>::max()}) {
        std::ostringstream os;
        os << v;
        std::string got;
        det::append_int(got, v);
        EXPECT_EQ(got, os.str());
    }
}

TEST(RowWriter, ChunkedOutputEqualsTheSerialLoop) {
    const std::vector<double> v = rw_values();
    std::string serial;
    const det::CNumber num;
    for (std::size_t i = 0; i < v.size(); ++i) {
        num.Append(serial, "%.17g", v[i]);
        serial += (i + 1) % 20 == 0 ? '\n' : ' ';
    }
    for (std::size_t chunk :
         {std::size_t{1}, std::size_t{7}, std::size_t{64}, std::size_t{100000}}) {
        std::ostringstream os;
        det::write_row_chunks(
            os, v.size(),
            [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                for (std::size_t i = First; i < Last; ++i) {
                    num.Append(rBuf, "%.17g", v[i]);
                    rBuf += (i + 1) % 20 == 0 ? '\n' : ' ';
                }
            },
            chunk);
        EXPECT_EQ(os.str(), serial) << chunk;
        std::string appended = "head:";
        det::append_row_chunks(
            appended, v.size(),
            [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                for (std::size_t i = First; i < Last; ++i) {
                    num.Append(rBuf, "%.17g", v[i]);
                    rBuf += (i + 1) % 20 == 0 ? '\n' : ' ';
                }
            },
            chunk);
        EXPECT_EQ(appended, "head:" + serial) << chunk;
    }
}
