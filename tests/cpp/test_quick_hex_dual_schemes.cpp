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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// Small versions of the hexahedral dual tests (`test_hex_dual_schemes.cpp`). The full suite
// runs deeper trees and twelve random ones, which take minutes in the instrumented coverage
// build, so a pull request skips it there; these build the trees of the five transition
// kinds (flat, concave and convex edge, concave and convex corner) one level shallower, and
// three random trees of three levels, and check what the full tests check: one hexahedron of
// positive volume per interior vertex, a conforming manifold mesh. Their suite name starts
// with `Quick`, which the job's exclusion does not match.

// System includes
#include <cmath>
#include <cstdint>
#include <exception>
#include <random>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "hex_dual_fixtures.hpp"

namespace {

namespace d = meshioplusplus::detail;
using namespace hdsfix;

}  // namespace

// A uniform tree is its own dual: no transition vertex, no scheme.
TEST(QuickHexDualSchemes, AUniformTreeIsItsOwnDual) {
    const d::LinearOctree t =
        hds_tree({{2, 2, 2}}, 3, 0, [](double, double, double) { return false; });
    const d::HexDualMesh dual = hds_check(t);
    EXPECT_EQ(dual.mNumTransitionVertices, 0);
    EXPECT_EQ(hds_schemes(dual), 0);
    EXPECT_EQ(dual.mHexes.size(), 27u);
}

// The four transition kinds one refinement of a corner of the forest makes, each with
// its scheme used.
TEST(QuickHexDualSchemes, EveryTransitionKindIsUsed) {
    struct Case {
        d::HexTransition mKind;
        bool (*mRefine)(double, double, double);
    };
    const Case cases[] = {
        {d::HexTransition::Flat, [](double X, double, double) { return X < 0.5; }},
        {d::HexTransition::Edge, [](double X, double Y, double) { return X < 0.5 || Y < 0.5; }},
        {d::HexTransition::Convex1, [](double X, double Y, double) { return X < 0.5 && Y < 0.5; }},
        {d::HexTransition::VertCenter,
         [](double X, double Y, double Z) { return X < 0.5 || Y < 0.5 || Z < 0.5; }},
    };
    for (const Case& c : cases) {
        const d::LinearOctree t = hds_tree({{2, 2, 2}}, 4, 1, c.mRefine);
        const d::HexDualMesh dual = hds_check(t);
        EXPECT_GT(dual.mNumTransitionVertices, 0);
        EXPECT_GT(dual.mNumSchemes[static_cast<std::size_t>(c.mKind)], 0)
            << static_cast<int>(c.mKind);
    }
    // The convex corner.
    hds_check(hds_tree({{2, 2, 2}}, 4, 1,
                       [](double X, double Y, double Z) { return X < 0.5 && Y < 0.5 && Z < 0.5; }));
}

// A thin shell of refinement two levels deep: conforming, and the transition hexahedra
// stay fair without fitting.
TEST(QuickHexDualSchemes, AShellTwoLevelsDownIsFair) {
    const d::LinearOctree t = hds_tree({{2, 2, 2}}, 4, 2, [](double X, double Y, double Z) {
        const double r =
            std::sqrt((X - 0.5) * (X - 0.5) + (Y - 0.5) * (Y - 0.5) + (Z - 0.5) * (Z - 0.5));
        return std::fabs(r - 0.3) < 0.12;
    });
    const d::HexDualMesh dual = hds_check(t);
    double worst = 1.0;
    for (const auto& h : dual.mHexes)
        worst = std::min(worst, d::hex_min_scaled_jacobian(dual.mPoints, h));
    EXPECT_GT(worst, 0.25);
}

// Random trees: every one conforms.
TEST(QuickHexDualSchemes, RandomTreesConform) {
    int failed = 0;
    std::string first;
    for (std::uint32_t seed = 1; seed <= 3; ++seed) {
        std::mt19937 rng(seed);
        const d::LinearOctree t =
            hds_tree({{2, 2, 2}}, 4, 2, [&](double, double, double) { return rng() % 6 == 0; });
        try {
            hds_check(t);
        } catch (const std::exception& e) {
            if (failed++ == 0)
                first = std::to_string(seed) + ": " + e.what();
        }
    }
    EXPECT_EQ(failed, 0) << first;
}
