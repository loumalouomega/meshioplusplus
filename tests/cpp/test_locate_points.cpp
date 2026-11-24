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
// Point location (roadmap 9.3): every volume cell type through simplexification,
// planar triangles, the shared-face rule, `interpolate` agreement and the
// refusals.

// System includes
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/locate_points.hpp"

namespace {

using meshioplusplus::locate_points;
using meshioplusplus::LocatePointsOptions;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

// Two unit hexahedra stacked along z (0 <= z <= 2).
Mesh lp_two_hexes() {
    std::vector<std::vector<double>> pts;
    for (int z = 0; z < 3; ++z)
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
                pts.push_back({double(x), double(y), double(z)});
    return mt::make_mesh(pts, "hexahedron",
                         {{0, 1, 3, 2, 4, 5, 7, 6}, {4, 5, 7, 6, 8, 9, 11, 10}});
}

// The weights reproduce the point from the nodes' coordinates.
void expect_reproduces(const Mesh& rMesh, const double* pP, const std::int64_t* pNodes,
                       const double* pW) {
    const NDArray& pts = rMesh.Points();
    const double* c = pts.As<double>();
    const std::size_t dim = rMesh.PointDim();
    double sum_w = 0.0;
    double q[3] = {0, 0, 0};
    for (int k = 0; k < 4; ++k) {
        if (pNodes[k] < 0)
            continue;
        sum_w += pW[k];
        for (std::size_t a = 0; a < dim; ++a)
            q[a] += pW[k] * c[static_cast<std::size_t>(pNodes[k]) * dim + a];
    }
    EXPECT_NEAR(sum_w, 1.0, 1e-12);
    for (std::size_t a = 0; a < dim; ++a)
        EXPECT_NEAR(q[a], pP[a], 1e-12);
}

}  // namespace

TEST(LocatePoints, HexahedraReportTheirInputCell) {
    const Mesh m = lp_two_hexes();
    const std::vector<double> p = {0.5, 0.5, 0.5, 0.25, 0.75, 1.5, 3.0, 3.0, 3.0};
    const auto r = locate_points(m, p);
    ASSERT_EQ(r.mCell.size(), 3u);
    EXPECT_EQ(r.mCell[0], 0);
    EXPECT_EQ(r.mCell[1], 1);
    EXPECT_EQ(r.mCell[2], -1);
    EXPECT_EQ(r.mNumLocated, 2);
    for (int i = 0; i < 2; ++i)
        expect_reproduces(m, &p[3 * i], &r.mNodes[4 * i], &r.mBarycentric[4 * i]);
    // Outside: no nodes, zero weights.
    for (int k = 0; k < 4; ++k) {
        EXPECT_EQ(r.mNodes[8 + k], -1);
        EXPECT_EQ(r.mBarycentric[8 + k], 0.0);
    }
}

TEST(LocatePoints, ASharedFaceBelongsToTheLowestCell) {
    const Mesh m = lp_two_hexes();
    // The plane z = 1 is shared by both hexahedra.
    const auto r = locate_points(m, {0.5, 0.5, 1.0});
    EXPECT_EQ(r.mCell[0], 0);
}

TEST(LocatePoints, TetrahedraAndTheTolerance) {
    const Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra",
                                 {{0, 1, 2, 3}});
    const auto in = locate_points(m, {0.1, 0.1, 0.1});
    EXPECT_EQ(in.mCell[0], 0);
    // Just beyond the slanted face x + y + z = 1, inside the bounding box: the
    // weight of the origin vertex is -1e-9.
    const std::vector<double> near_face = {0.4, 0.4, 0.2 + 1e-9};
    EXPECT_EQ(locate_points(m, near_face).mCell[0], -1);
    LocatePointsOptions loose;
    loose.mTolerance = 1e-6;
    EXPECT_EQ(locate_points(m, near_face, loose).mCell[0], 0);
}

TEST(LocatePoints, PlanarTrianglesAreSearchedInTheirPlane) {
    const Mesh m = mt::make_mesh({{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}}, "triangle",
                                 {{0, 1, 2}, {0, 2, 3}});
    const auto r = locate_points(m, {1.5, 0.5, 0.0, 0.5, 1.5, 0.0, 3.0, 3.0, 0.0});
    EXPECT_EQ(r.mCell[0], 0);
    EXPECT_EQ(r.mCell[1], 1);
    EXPECT_EQ(r.mCell[2], -1);
    EXPECT_EQ(r.mNodes[3], -1);  // a triangle has three nodes
    expect_reproduces(m, std::array<double, 3>{1.5, 0.5, 0.0}.data(), &r.mNodes[0],
                      &r.mBarycentric[0]);
}

TEST(LocatePoints, QuadrilateralsMapToTheirOwnCell) {
    const Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}},
                                 "quad", {{0, 1, 4, 3}, {1, 2, 5, 4}});
    const auto r = locate_points(m, {0.5, 0.5, 0.0, 1.5, 0.5, 0.0});
    EXPECT_EQ(r.mCell[0], 0);
    EXPECT_EQ(r.mCell[1], 1);
}

TEST(LocatePoints, AgreesWithInterpolateBitForBit) {
    Mesh src = lp_two_hexes();
    std::vector<double> f;
    for (std::size_t i = 0; i < src.NumPoints(); ++i)
        f.push_back(0.37 * i * i - 1.9 * i + 0.11);
    src.AddPointData("f", mt::data_array(f));
    const std::vector<double> p = {0.3, 0.6, 0.4, 0.8, 0.1, 1.7, 0.5, 0.5, 1.0};
    const Mesh target = mt::make_mesh({{p[0], p[1], p[2]}, {p[3], p[4], p[5]}, {p[6], p[7], p[8]}},
                                      "vertex", {{0}, {1}, {2}});
    meshioplusplus::InterpolateOptions io;
    io.mMethod = meshioplusplus::InterpolateMethod::Barycentric;
    const Mesh out = meshioplusplus::interpolate(src, target, io);
    const auto r = locate_points(src, p);
    for (int i = 0; i < 3; ++i) {
        ASSERT_GE(r.mCell[i], 0);
        double acc = r.mBarycentric[4 * i] * f[static_cast<std::size_t>(r.mNodes[4 * i])];
        for (int k = 1; k < 4; ++k)
            if (r.mNodes[4 * i + k] >= 0)
                acc += r.mBarycentric[4 * i + k] * f[static_cast<std::size_t>(r.mNodes[4 * i + k])];
        EXPECT_EQ(acc, out.PointData("f").As<double>()[i]) << i;
    }
}

TEST(LocatePoints, Refusals) {
    const Mesh m = lp_two_hexes();
    EXPECT_THROW(locate_points(m, {0, 0}), std::invalid_argument);
    EXPECT_THROW(locate_points(m, {std::nan(""), 0, 0}), std::invalid_argument);
    LocatePointsOptions bad;
    bad.mTolerance = -1.0;
    EXPECT_THROW(locate_points(m, {0, 0, 0}, bad), std::invalid_argument);
    const Mesh lines = mt::make_mesh({{0, 0, 0}, {1, 0, 0}}, "line", {{0, 1}});
    EXPECT_THROW(locate_points(lines, {0.5, 0, 0}), std::invalid_argument);
    EXPECT_TRUE(locate_points(m, {}).mCell.empty());
}
