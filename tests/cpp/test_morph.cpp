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
// Boundary-driven morphing by harmonic extension (roadmap 7.4.1).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/morph.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// An n^3 hexahedral unit box; the displacement lifts its top face by Lift.
Mesh mp_box(int N, double Lift) {
    std::vector<std::vector<double>> pts;
    for (int k = 0; k <= N; ++k)
        for (int j = 0; j <= N; ++j)
            for (int i = 0; i <= N; ++i)
                pts.push_back({double(i) / N, double(j) / N, double(k) / N});
    const auto id = [&](int i, int j, int k) {
        return std::int64_t((k * (N + 1) + j) * (N + 1) + i);
    };
    std::vector<std::vector<std::int64_t>> hexes;
    for (int k = 0; k < N; ++k)
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i)
                hexes.push_back({id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k),
                                 id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1),
                                 id(i, j + 1, k + 1)});
    Mesh m = mt::make_mesh(pts, "hexahedron", hexes);
    NDArray d(DType::Float64, {pts.size(), std::size_t{3}});
    for (std::size_t p = 0; p < pts.size(); ++p) {
        d.As<double>()[3 * p] = 0.0;
        d.As<double>()[3 * p + 1] = 0.0;
        // a bump on the top face: lifted most at its centre
        const double x = pts[p][0], y = pts[p][1];
        d.As<double>()[3 * p + 2] = pts[p][2] == 1.0 ? Lift * 16 * x * (1 - x) * y * (1 - y) : 0.0;
    }
    m.AddPointData("d", std::move(d));
    return m;
}

}  // namespace

TEST(Morph, LiftingAHexBoxTopInvertsNothing) {
    // The roadmap probe: lift the top face by 0.1; the interior follows.
    MorphOptions o;
    o.mArray = "d";
    const MorphResult r = morph(mp_box(6, 0.1), o);
    EXPECT_EQ(r.mNumInverted, 0);
    EXPECT_TRUE(r.mConverged);
    EXPECT_EQ(r.mNumFree, 5 * 5 * 5);
    EXPECT_NEAR(r.mMaxDisplacement, 0.1, 1e-12);
    // The centre column rises monotonically from 0 at the bottom.
    const double* p = r.mMesh.Points().As<double>();
    double last = -1.0;
    for (int k = 0; k <= 6; ++k) {
        const std::size_t c = static_cast<std::size_t>((k * 7 + 3) * 7 + 3);
        const double rise = p[3 * c + 2] - k / 6.0;
        EXPECT_GE(rise, last - 1e-12);
        last = rise;
    }
    EXPECT_NEAR(last, 0.1, 1e-12);
}

TEST(Morph, TooLargeAMotionIsRefusedThenAllowed) {
    MorphOptions o;
    o.mArray = "d";
    // Pushing the top face down by 1.6 times the box's height turns cells over.
    Mesh m = mp_box(4, -1.6);
    EXPECT_THROW(morph(m, o), std::invalid_argument);
    o.mAllowInverted = true;
    const MorphResult r = morph(m, o);
    EXPECT_GT(r.mNumInverted, 0);
}

TEST(Morph, RigidTranslationOfEveryBoundaryPointMovesEverything) {
    Mesh m = mp_box(3, 0.0);
    NDArray d(DType::Float64, {m.NumPoints(), std::size_t{3}});
    for (std::size_t p = 0; p < m.NumPoints(); ++p) {
        d.As<double>()[3 * p] = 0.25;
        d.As<double>()[3 * p + 1] = -0.5;
        d.As<double>()[3 * p + 2] = 0.125;
    }
    m.AddPointData("d", std::move(d));
    MorphOptions o;
    o.mArray = "d";
    o.mSteps = 3;
    const MorphResult r = morph(m, o);
    const double* p0 = m.Points().As<double>();
    const double* p1 = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        EXPECT_NEAR(p1[3 * i] - p0[3 * i], 0.25, 1e-9);
        EXPECT_NEAR(p1[3 * i + 1] - p0[3 * i + 1], -0.5, 1e-9);
        EXPECT_NEAR(p1[3 * i + 2] - p0[3 * i + 2], 0.125, 1e-9);
    }
}

TEST(Morph, Refusals) {
    Mesh m = mp_box(2, 0.1);
    MorphOptions o;
    EXPECT_THROW(morph(m, o), std::invalid_argument);  // no array
    o.mArray = "missing";
    EXPECT_THROW(morph(m, o), std::invalid_argument);
    o.mArray = "d";
    o.mSteps = 0;
    EXPECT_THROW(morph(m, o), std::invalid_argument);
    o.mSteps = 1;
    o.mStiffening = -1.0;
    EXPECT_THROW(morph(m, o), std::invalid_argument);
    o.mStiffening = 1.0;
    EXPECT_THROW(morph(mt::line_mesh(), o), std::invalid_argument);
}
