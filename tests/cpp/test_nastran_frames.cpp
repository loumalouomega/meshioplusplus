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
 * @file test_nastran_frames.cpp
 * @brief The Nastran coordinate-system kernel the result readers share:
 *        closed-form positions and bases of rectangular, cylindrical and
 *        spherical systems, CORD2 chains through a cylindrical reference,
 *        CORD1 through GRIDs, and systems that cannot be resolved. The
 *        comparison with pyNastran on real files is in test_nastran_h5.py and
 *        test_nastran_op2.py.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/detail/nastran_model.hpp"
#include "meshioplusplus/mesh.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::detail::NastranCoordCard;

NastranCoordCard cord2(std::int64_t Cid, int Type, std::int64_t Rid, std::vector<double> Abc) {
    NastranCoordCard c;
    c.mCid = Cid;
    c.mType = Type;
    c.mRid = Rid;
    for (std::size_t k = 0; k < 9; ++k)
        c.mAbc[k] = Abc[k];
    return c;
}

Mesh points_mesh(const std::vector<double>& rXyz) {
    Mesh m;
    NDArray p(DType::Float64, {rXyz.size() / 3, 3});
    std::copy(rXyz.begin(), rXyz.end(), p.As<double>());
    m.AssignPoints(std::move(p));
    return m;
}

double at(const Mesh& rMesh, std::size_t Point, std::size_t Axis) {
    return rMesh.Points().As<double>()[3 * Point + Axis];
}

}  // namespace

TEST(NastranFrames, CylindricalAndSphericalPositionsAndBases) {
    // CORD2C 1: origin (1, 0, 0), axes of basic. CORD2S 2: axes of basic.
    const std::vector<NastranCoordCard> cards = {
        cord2(1, 2, 0, {1, 0, 0, 1, 0, 1, 2, 0, 0}),
        cord2(2, 3, 0, {0, 0, 0, 0, 0, 1, 1, 0, 0}),
    };
    // (r=2, theta=90, z=1) in 1; (r=1, theta=90, phi=90) in 2.
    Mesh m = points_mesh({2, 90, 1, 1, 90, 90});
    const auto sys =
        meshioplusplus::detail::nastran_apply_frames(m, cards, {10, 20}, {1, 2}, {1, 2}, "test");
    EXPECT_NEAR(at(m, 0, 0), 1.0, 1e-15);
    EXPECT_NEAR(at(m, 0, 1), 2.0, 1e-15);
    EXPECT_NEAR(at(m, 0, 2), 1.0, 1e-15);
    EXPECT_NEAR(at(m, 1, 0), 0.0, 1e-15);
    EXPECT_NEAR(at(m, 1, 1), 1.0, 1e-15);
    EXPECT_NEAR(at(m, 1, 2), 0.0, 1e-15);
    ASSERT_TRUE(m.HasPointData("nastran:cp"));
    ASSERT_TRUE(m.HasPointData("nastran:cd"));

    // Radial at (x=1, y=2) around the axis through (1, 0): +y. Tangential: -x.
    const double p0[3] = {at(m, 0, 0), at(m, 0, 1), at(m, 0, 2)};
    double v[6] = {1, 0, 0, 0, 1, 0};
    ASSERT_TRUE(meshioplusplus::detail::nastran_rotate_to_basic(sys, 1, p0, v, 2));
    EXPECT_NEAR(v[0], 0.0, 1e-15);
    EXPECT_NEAR(v[1], 1.0, 1e-15);
    EXPECT_NEAR(v[3], -1.0, 1e-15);
    EXPECT_NEAR(v[4], 0.0, 1e-15);
    // Spherical at +y: e_r = +y, e_theta = -z, e_phi = -x.
    const double p1[3] = {0, 1, 0};
    double s[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    ASSERT_TRUE(meshioplusplus::detail::nastran_rotate_to_basic(sys, 2, p1, s, 3));
    EXPECT_NEAR(s[1], 1.0, 1e-15);
    EXPECT_NEAR(s[5], -1.0, 1e-15);
    EXPECT_NEAR(s[6], -1.0, 1e-15);
    // CD 0 and unknown systems leave the values alone.
    double u[3] = {1, 2, 3};
    EXPECT_FALSE(meshioplusplus::detail::nastran_rotate_to_basic(sys, 0, p1, u, 1));
    EXPECT_FALSE(meshioplusplus::detail::nastran_rotate_to_basic(sys, 99, p1, u, 1));
    EXPECT_EQ(u[2], 3.0);
}

TEST(NastranFrames, ChainsResolveInAnyOrder) {
    // CORD1R 7 through GRIDs 1, 2, 3, where GRID 2 is in CORD2R 4, which is
    // defined in CORD2C 3 (listed last): every link waits for the next.
    NastranCoordCard c1;
    c1.mCid = 7;
    c1.mByGrids = true;
    c1.mGrids[0] = 1;
    c1.mGrids[1] = 2;
    c1.mGrids[2] = 3;
    const std::vector<NastranCoordCard> cards = {
        c1,
        // origin at (r=1, theta=0, z=0) of 3, i.e. basic (1, 0, 0); +z; +x
        cord2(4, 1, 3, {1, 0, 0, 1, 0, 1, 2, 0, 0}),
        cord2(3, 2, 0, {0, 0, 0, 0, 0, 1, 1, 0, 0}),
    };
    // GRID 1 at basic origin, GRID 2 at (0, 0, 1) of 4 = basic (1, 0, 1),
    // GRID 3 at (0, 1, 0); GRID 4 at (0, 0, 2) of CORD1R 7.
    Mesh m = points_mesh({0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2});
    meshioplusplus::detail::nastran_apply_frames(m, cards, {1, 2, 3, 4}, {0, 4, 0, 7}, {0, 0, 0, 0},
                                                 "test");
    EXPECT_NEAR(at(m, 1, 0), 1.0, 1e-15);
    EXPECT_NEAR(at(m, 1, 2), 1.0, 1e-15);
    // CORD1R 7: z along (1, 0, 1)/sqrt(2) from the origin.
    EXPECT_NEAR(at(m, 3, 0), 2.0 / std::sqrt(2.0), 1e-15);
    EXPECT_NEAR(at(m, 3, 1), 0.0, 1e-15);
    EXPECT_NEAR(at(m, 3, 2), 2.0 / std::sqrt(2.0), 1e-15);
}

TEST(NastranFrames, UnresolvableSystemsKeepTheCoordinates) {
    const std::vector<NastranCoordCard> cards = {
        cord2(8, 1, 9, {0, 0, 0, 0, 0, 1, 1, 0, 0}),  // circular: 8 in 9, 9 in 8
        cord2(9, 1, 8, {0, 0, 0, 0, 0, 1, 1, 0, 0}),
        cord2(10, 1, 0, {0, 0, 0, 0, 0, 0, 1, 0, 0}),  // degenerate: B == A
    };
    Mesh m = points_mesh({1, 2, 3, 4, 5, 6, 7, 8, 9});
    meshioplusplus::detail::nastran_apply_frames(m, cards, {1, 2, 3}, {8, 10, 55}, {0, 0, 0},
                                                 "test");
    for (std::size_t k = 0; k < 9; ++k)
        EXPECT_EQ(m.Points().As<double>()[k], static_cast<double>(k + 1));
    ASSERT_TRUE(m.HasPointData("nastran:cp"));
    EXPECT_FALSE(m.HasPointData("nastran:cd"));
}
