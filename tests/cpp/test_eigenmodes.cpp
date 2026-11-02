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
// Laplace-Beltrami eigenmodes of a triangle surface (roadmap 7.3.2).

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/eigenmodes.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

TEST(Eigenmodes, UnitSquareNeumann) {
    // A unit square, n x n quads cut into triangles: Neumann eigenvalues
    // pi^2 (a^2 + b^2): 0, pi^2, pi^2, 2 pi^2, 4 pi^2, 4 pi^2.
    const int m = 24;
    std::vector<std::vector<double>> pts;
    for (int j = 0; j <= m; ++j)
        for (int i = 0; i <= m; ++i)
            pts.push_back({double(i) / m, double(j) / m, 0.0});
    std::vector<std::vector<std::int64_t>> tris;
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < m; ++i) {
            const std::int64_t a = j * (m + 1) + i, b = a + 1, c = a + m + 2, d = a + m + 1;
            tris.push_back({a, b, c});
            tris.push_back({a, c, d});
        }
    EigenmodeOptions o;
    o.mNumModes = 6;
    const EigenmodeResult r = laplacian_eigenmodes(mt::make_mesh(pts, "triangle", tris), o);
    ASSERT_TRUE(r.mConverged);
    const double pi2 = M_PI * M_PI;
    const std::array<double, 6> exact{0.0, pi2, pi2, 2 * pi2, 4 * pi2, 4 * pi2};
    EXPECT_EQ(r.mEigenvalues[0], 0.0);
    for (std::size_t j = 1; j < 6; ++j)
        EXPECT_NEAR(r.mEigenvalues[j] / exact[j], 1.0, 0.02) << j;
    // M-orthonormal: the constant mode is 1 / sqrt(area) = 1.
    const double* modes = r.mMesh.PointData(kLaplacianModesName).As<double>();
    EXPECT_NEAR(modes[0], 1.0, 1e-8);
}

TEST(Eigenmodes, DeterministicAndRefusals) {
    const Mesh tri = mt::tri_mesh();
    EigenmodeOptions o;
    o.mNumModes = 2;
    const EigenmodeResult a = laplacian_eigenmodes(tri, o);
    const EigenmodeResult b = laplacian_eigenmodes(tri, o);
    EXPECT_EQ(a.mEigenvalues, b.mEigenvalues);
    o.mNumModes = 0;
    EXPECT_THROW(laplacian_eigenmodes(tri, o), std::invalid_argument);
    o.mNumModes = 1000;
    EXPECT_THROW(laplacian_eigenmodes(tri, o), std::invalid_argument);
    EXPECT_THROW(laplacian_eigenmodes(mt::tet_mesh()), std::invalid_argument);
    EXPECT_THROW(laplacian_eigenmodes(mt::quad_mesh()), std::invalid_argument);
}
