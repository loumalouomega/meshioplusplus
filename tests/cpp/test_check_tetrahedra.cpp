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
// check_tetrahedra: hand-built tetrahedra whose defects are known by construction.

// System includes
#include <cstdint>
#include <stdexcept>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/check_tetrahedra.hpp"

using meshioplusplus::check_tetrahedra;
using meshioplusplus::Mesh;
using meshioplusplus::TetrahedraCheckOptions;

namespace {

const std::vector<std::vector<double>> kCorner = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

}  // namespace

TEST(CheckTetrahedra, APositiveTetrahedronIsValid) {
    const auto r = check_tetrahedra(mt::make_mesh(kCorner, "tetra", {{0, 1, 2, 3}}));
    EXPECT_EQ(r.mNumTetrahedra, 1);
    EXPECT_EQ(r.mNumInverted, 0);
    EXPECT_EQ(r.mNumDegenerate, 0);
    EXPECT_EQ(r.mNumOverlaps, 0);
    EXPECT_TRUE(r.mValid);
    EXPECT_NEAR(r.mMinVolume, 1.0 / 6.0, 1e-15);
    EXPECT_NEAR(r.mMaxVolume, 1.0 / 6.0, 1e-15);
    EXPECT_TRUE(r.mWarnings.empty());
}

TEST(CheckTetrahedra, SwappingTwoNodesInvertsIt) {
    const auto r = check_tetrahedra(mt::make_mesh(kCorner, "tetra", {{0, 2, 1, 3}}));
    EXPECT_EQ(r.mNumInverted, 1);
    EXPECT_EQ(r.mNumDegenerate, 0);
    EXPECT_FALSE(r.mValid);
    EXPECT_NEAR(r.mMinVolume, -1.0 / 6.0, 1e-15);
}

TEST(CheckTetrahedra, AFlatTetrahedronIsDegenerate) {
    const auto r = check_tetrahedra(
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}, "tetra", {{0, 1, 2, 3}}));
    EXPECT_EQ(r.mNumDegenerate, 1);
    EXPECT_EQ(r.mNumInverted, 0);
    EXPECT_FALSE(r.mValid);
}

TEST(CheckTetrahedra, TheToleranceIsRelativeToTheEdgesAndTheFactorIsAnOption) {
    // Height 1e-4 over a unit base: a volume of 1.7e-5, valid at the default factor and
    // degenerate once the factor reaches 1e-4.
    const Mesh m =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.2, 0.2, 1e-4}}, "tetra", {{0, 1, 2, 3}});
    EXPECT_TRUE(check_tetrahedra(m).mValid);
    TetrahedraCheckOptions options;
    options.mZeroVolumeFactor = 1e-3;
    EXPECT_EQ(check_tetrahedra(m, options).mNumDegenerate, 1);
}

TEST(CheckTetrahedra, OppositeNodesOnTheSameSideOfASharedFaceOverlap) {
    // Two tetrahedra over the face (0, 1, 2): the fourth nodes on opposite sides are a valid
    // pair, on the same side they overlap.
    const Mesh apart = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}},
                                     "tetra", {{0, 1, 2, 3}, {0, 2, 1, 4}});
    const auto ok = check_tetrahedra(apart);
    EXPECT_EQ(ok.mNumOverlaps, 0);
    EXPECT_TRUE(ok.mValid);

    const Mesh overlapping =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2, 0.2, 0.5}}, "tetra",
                      {{0, 1, 2, 3}, {0, 1, 2, 4}});
    const auto bad = check_tetrahedra(overlapping);
    EXPECT_EQ(bad.mNumTetrahedra, 2);
    EXPECT_EQ(bad.mNumOverlaps, 1);
    EXPECT_FALSE(bad.mValid);
}

TEST(CheckTetrahedra, AFaceThreeTetrahedraShareIsNotAPair) {
    // Three tetrahedra on the face (0, 1, 2): not used by exactly two, so not examined.
    const auto r = check_tetrahedra(mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2, 0.2, 0.5}, {0.3, 0.3, 0.4}}, "tetra",
        {{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 2, 5}}));
    EXPECT_EQ(r.mNumTetrahedra, 3);
    EXPECT_EQ(r.mNumOverlaps, 0);
}

TEST(CheckTetrahedra, QuadraticTetrahedraAreExaminedByTheirCorners) {
    const auto r = check_tetrahedra(mt::tet10_mesh());
    EXPECT_EQ(r.mNumTetrahedra, 1);
}

TEST(CheckTetrahedra, ANonTetrahedralMeshIsReportedNotExamined) {
    const auto r = check_tetrahedra(mt::hex_mesh());
    EXPECT_EQ(r.mNumTetrahedra, 0);
    EXPECT_TRUE(r.mValid);
    EXPECT_EQ(r.mWarnings.size(), 1u);
}

TEST(CheckTetrahedra, ABadFactorIsRefused) {
    TetrahedraCheckOptions options;
    options.mZeroVolumeFactor = -1.0;
    EXPECT_THROW(check_tetrahedra(mt::tet_mesh(), options), std::invalid_argument);
}
