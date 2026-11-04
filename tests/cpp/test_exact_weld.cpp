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
 * @file test_exact_weld.cpp
 * @brief The private vertex weld of the STL, glTF and WKT readers
 *        (`src/cpp/src/detail/exact_weld.hpp`, roadmap 5.7.3) against the
 *        serial hash-map sweeps it replaced, and its two equalities.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/exact_weld.hpp"

namespace {

namespace det = meshioplusplus::detail;

std::vector<double> ew_rows(std::size_t n, std::uint64_t seed) {
    // Many repeats: values drawn from a small set, so rows collide often.
    std::vector<double> out(n * 3);
    for (double& v : out) {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        v = static_cast<double>(seed % 5) * 0.25;
    }
    return out;
}

}  // namespace

TEST(ExactWeld, BitsMatchesTheSerialByteMapSweep) {
    const std::size_t n = 50000;  // large enough for the chunked counting sort
    const std::vector<double> rows = ew_rows(n, 88172645463325252ull);
    const det::ExactWeld weld = det::exact_weld(rows.data(), n, 3, det::ExactWeldMode::Bits);

    std::unordered_map<std::string, std::int64_t> seen;
    std::vector<std::size_t> first;
    ASSERT_EQ(weld.mIdOfRow.size(), n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::string key(reinterpret_cast<const char*>(rows.data() + 3 * i), 24);
        const auto it = seen.emplace(key, static_cast<std::int64_t>(first.size()));
        if (it.second)
            first.push_back(i);
        ASSERT_EQ(weld.mIdOfRow[i], it.first->second) << "row " << i;
    }
    EXPECT_EQ(weld.mFirstRow, first);
}

TEST(ExactWeld, BitsKeepsSignedZerosApartAndANaNWithItself) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> rows = {0.0, 1.0, 2.0, -0.0, 1.0, 2.0, nan, 0.0,
                                      0.0, nan, 0.0, 0.0,  0.0, 1.0, 2.0};
    const det::ExactWeld weld = det::exact_weld(rows.data(), 5, 3, det::ExactWeldMode::Bits);
    EXPECT_EQ(weld.mIdOfRow, (std::vector<std::int64_t>{0, 1, 2, 2, 0}));
}

TEST(ExactWeld, ValuesMergesSignedZerosAndNeverANaN) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> rows = {0.0, 1.0, 2.0, -0.0, 1.0, 2.0, nan, 0.0,
                                      0.0, nan, 0.0, 0.0,  0.0, 1.0, 2.0};
    const det::ExactWeld weld = det::exact_weld(rows.data(), 5, 3, det::ExactWeldMode::Values);
    EXPECT_EQ(weld.mIdOfRow, (std::vector<std::int64_t>{0, 0, 1, 2, 0}));
    // The first row of an id is the one kept, its bits included.
    EXPECT_EQ(weld.mFirstRow, (std::vector<std::size_t>{0, 2, 3}));
}

TEST(ExactWeld, RowsOfDifferentWidthsNeverMeet) {
    std::vector<det::ExactWeldKey> keys;
    const double a[4] = {1.0, 2.0, 0.0, 0.0};
    keys.push_back(det::exact_weld_key(a, 2, det::ExactWeldMode::Values, 0));
    keys.push_back(det::exact_weld_key(a, 3, det::ExactWeldMode::Values, 1));
    keys.push_back(det::exact_weld_key(a, 2, det::ExactWeldMode::Values, 2));
    const det::ExactWeld weld = det::exact_weld_keys(keys);
    EXPECT_EQ(weld.mIdOfRow, (std::vector<std::int64_t>{0, 1, 0}));
}

TEST(ExactWeld, NoRowsNoIds) {
    const det::ExactWeld weld =
        det::exact_weld(static_cast<const double*>(nullptr), 0, 3, det::ExactWeldMode::Bits);
    EXPECT_TRUE(weld.mIdOfRow.empty());
    EXPECT_EQ(weld.NumIds(), 0u);
}
