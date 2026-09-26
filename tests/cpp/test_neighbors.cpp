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
 * @file test_neighbors.cpp
 * @brief `neighbor_pairs` against brute force: the same pairs, radius
 *        inclusive, k-nearest ties to the lower index, minimum images.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/neighbors.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::NDArray;
using meshioplusplus::NeighborMethod;
using meshioplusplus::NeighborOptions;

NDArray nb_cloud(std::size_t n, std::size_t d, std::uint64_t seed, double scale, bool snap) {
    NDArray p = NDArray::Uninit(DType::Float64, {n, d});
    std::uint64_t s = seed;
    for (std::size_t i = 0; i < n * d; ++i) {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        double v = static_cast<double>(s >> 11) / 9007199254740992.0 * scale;
        if (snap)  // a coarse lattice: many exact ties and exact-cutoff pairs
            v = std::floor(v * 8.0) / 8.0;
        p.As<double>()[i] = v;
    }
    return p;
}

double nb_brute_d2(const NDArray& rP, std::size_t i, std::size_t j,
                   const std::vector<double>& box) {
    const std::size_t d = rP.Shape()[1];
    double s = 0.0;
    for (std::size_t a = 0; a < d; ++a) {
        double x = rP.As<double>()[i * d + a] - rP.As<double>()[j * d + a];
        if (!box.empty())
            x -= box[a] * std::nearbyint(x / box[a]);
        s = a == 0 ? x * x : s + x * x;
    }
    return s;
}

std::set<std::pair<std::int64_t, std::int64_t>> nb_pairs(const meshioplusplus::NeighborPairs& rR) {
    std::set<std::pair<std::int64_t, std::int64_t>> out;
    for (std::size_t e = 0; e < rR.mSource.Size(); ++e)
        EXPECT_TRUE(out.insert({rR.mSource.As<std::int64_t>()[e], rR.mTarget.As<std::int64_t>()[e]})
                        .second);  // never the same pair twice
    return out;
}

TEST(Neighbors, RadiusMatchesBruteForce) {
    for (std::size_t d = 1; d <= 3; ++d)
        for (const bool snap : {false, true})
            for (const bool periodic : {false, true}) {
                const NDArray p = nb_cloud(1500, d, 17 + d, 1.0, snap);
                NeighborOptions o;
                o.mMethod = NeighborMethod::Radius;
                o.mRadius = d == 1 ? 0.004 : (d == 2 ? 0.05 : 0.125);
                if (periodic)
                    o.mBox.assign(d, 1.0);
                std::set<std::pair<std::int64_t, std::int64_t>> expect;
                for (std::size_t i = 0; i < 1500; ++i)
                    for (std::size_t j = i + 1; j < 1500; ++j)
                        if (nb_brute_d2(p, i, j, o.mBox) <= o.mRadius * o.mRadius)
                            expect.insert(
                                {static_cast<std::int64_t>(i), static_cast<std::int64_t>(j)});
                EXPECT_EQ(nb_pairs(meshioplusplus::neighbor_pairs(p, o)), expect)
                    << "d=" << d << " snap=" << snap << " periodic=" << periodic;
            }
}

TEST(Neighbors, KNearestMatchesBruteForceWithTiesToTheLowerIndex) {
    for (std::size_t d = 1; d <= 3; ++d)
        for (const bool snap : {false, true})
            for (const bool periodic : {false, true}) {
                const std::size_t n = 1200;
                const NDArray p = nb_cloud(n, d, 5 + d, 1.0, snap);
                NeighborOptions o;
                o.mMethod = NeighborMethod::KNearest;
                o.mK = 7;
                if (periodic)
                    o.mBox.assign(d, 1.0);
                const auto got = meshioplusplus::neighbor_pairs(p, o);
                ASSERT_EQ(got.mSource.Size(), n * 7);
                for (std::size_t q = 0; q < n; ++q) {
                    std::vector<std::pair<double, std::int64_t>> all;
                    for (std::size_t j = 0; j < n; ++j)
                        if (j != q)
                            all.emplace_back(nb_brute_d2(p, q, j, o.mBox),
                                             static_cast<std::int64_t>(j));
                    std::sort(all.begin(), all.end());
                    for (std::size_t t = 0; t < 7; ++t) {
                        ASSERT_EQ(got.mSource.As<std::int64_t>()[q * 7 + t],
                                  static_cast<std::int64_t>(q));
                        ASSERT_EQ(got.mTarget.As<std::int64_t>()[q * 7 + t], all[t].second)
                            << "d=" << d << " snap=" << snap << " periodic=" << periodic
                            << " q=" << q << " rank " << t;
                    }
                }
            }
}

TEST(Neighbors, TheCellSizeNeverChangesTheAnswer) {
    const NDArray p = nb_cloud(3000, 3, 99, 1.0, false);
    NeighborOptions o;
    o.mMethod = NeighborMethod::KNearest;
    o.mK = 9;
    const auto reference = nb_pairs(meshioplusplus::neighbor_pairs(p, o));
    for (const double cell : {0.01, 0.07, 0.4, 3.0}) {
        o.mCellSize = cell;
        EXPECT_EQ(nb_pairs(meshioplusplus::neighbor_pairs(p, o)), reference) << cell;
    }
}

TEST(Neighbors, BadRequestsAreRefused) {
    const NDArray p = nb_cloud(10, 3, 1, 1.0, false);
    NeighborOptions o;
    o.mRadius = 0.0;
    EXPECT_THROW(meshioplusplus::neighbor_pairs(p, o), std::invalid_argument);
    o.mRadius = 0.6;
    o.mBox = {1.0, 1.0, 1.0};
    EXPECT_THROW(meshioplusplus::neighbor_pairs(p, o), std::invalid_argument);
    o.mBox = {1.0, 1.0};
    o.mRadius = 0.1;
    EXPECT_THROW(meshioplusplus::neighbor_pairs(p, o), std::invalid_argument);
    NDArray bad = NDArray::Uninit(DType::Float64, {4, 4});
    EXPECT_THROW(meshioplusplus::neighbor_pairs(bad, NeighborOptions{}), std::invalid_argument);
}

}  // namespace
