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

// The second tier of parametric generators (curves, point clouds, surfaces and
// solids) and the full-Lagrange / VTK Lagrange orders every generator builds.

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/operations/primitives.hpp"

namespace mio = meshioplusplus;
using meshioplusplus::Mesh;

namespace {

using PsVec = std::array<double, 3>;

PsVec ps_point(const Mesh& rMesh, std::int64_t I) {
    const double* p = rMesh.Points().As<double>();
    return {p[3 * I], p[3 * I + 1], p[3 * I + 2]};
}

PsVec ps_sub(const PsVec& A, const PsVec& B) {
    return {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
}

PsVec ps_cross(const PsVec& A, const PsVec& B) {
    return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}

double ps_dot(const PsVec& A, const PsVec& B) {
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

double ps_norm(const PsVec& A) {
    return std::sqrt(ps_dot(A, A));
}

std::size_t ps_num_cells(const Mesh& rMesh) {
    std::size_t n = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        n += rMesh.Cells(b).NumCells();
    return n;
}

// Cell C of a block and its node count, rectangular or ragged.
template <typename View>
const std::int64_t* ps_row(const View& rView, std::size_t C, std::size_t& rSize) {
    if (rView.IsRagged()) {
        rSize = rView.RowSize(C);
        return rView.Row(C);
    }
    rSize = rView.NodesPerCell();
    return rView.Conn().template As<std::int64_t>() + C * rSize;
}

// The volume a closed surface of triangles, quads and pentagons encloses
// (positive when it is oriented outward), by fanning each face from node 0.
double ps_enclosed_volume(const Mesh& rMesh) {
    double v = 0.0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto view = rMesh.Cells(b);
        const std::string type = view.Type();
        for (std::size_t c = 0; c < view.NumCells(); ++c) {
            std::size_t width = 0;
            const std::int64_t* row = ps_row(view, c, width);
            const std::size_t corners = type.rfind("triangle", 0) == 0 ? 3
                                        : type.rfind("quad", 0) == 0   ? 4
                                                                       : width;
            for (std::size_t k = 1; k + 1 < corners; ++k)
                v += ps_dot(ps_point(rMesh, row[0]),
                            ps_cross(ps_point(rMesh, row[k]), ps_point(rMesh, row[k + 1])));
        }
    }
    return v / 6.0;
}

// V - E + F over the corners of the surface cells.
long ps_euler(const Mesh& rMesh) {
    std::set<std::int64_t> verts;
    std::set<std::pair<std::int64_t, std::int64_t>> edges;
    long faces = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto view = rMesh.Cells(b);
        const std::string type = view.Type();
        for (std::size_t c = 0; c < view.NumCells(); ++c, ++faces) {
            std::size_t width = 0;
            const std::int64_t* row = ps_row(view, c, width);
            const std::size_t k = type.rfind("triangle", 0) == 0 ? 3
                                  : type.rfind("quad", 0) == 0   ? 4
                                                                 : width;
            for (std::size_t i = 0; i < k; ++i) {
                const std::int64_t a = row[i], e = row[(i + 1) % k];
                verts.insert(a);
                edges.insert({std::min(a, e), std::max(a, e)});
            }
        }
    }
    return static_cast<long>(verts.size()) - static_cast<long>(edges.size()) + faces;
}

// The smallest corner Jacobian of the tetrahedra and hexahedra.
double ps_min_jacobian(const Mesh& rMesh) {
    static const int kHexCorners[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5}, {2, 3, 1, 6}, {3, 0, 2, 7},
                                          {4, 7, 5, 0}, {5, 4, 6, 1}, {6, 5, 7, 2}, {7, 6, 4, 3}};
    static const int kTet[1][4] = {{0, 1, 2, 3}};
    double worst = 1e300;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto view = rMesh.Cells(b);
        const std::string type = view.Type();
        const bool hex = type.rfind("hexahedron", 0) == 0;
        if (!hex && type.rfind("tetra", 0) != 0)
            continue;
        for (std::size_t c = 0; c < view.NumCells(); ++c) {
            std::size_t width = 0;
            const std::int64_t* row = ps_row(view, c, width);
            const auto& corners = hex ? kHexCorners : kTet;
            for (std::size_t q = 0; q < (hex ? 8u : 1u); ++q) {
                const PsVec a = ps_point(rMesh, row[corners[q][0]]);
                const double j = ps_dot(ps_sub(ps_point(rMesh, row[corners[q][1]]), a),
                                        ps_cross(ps_sub(ps_point(rMesh, row[corners[q][2]]), a),
                                                 ps_sub(ps_point(rMesh, row[corners[q][3]]), a)));
                worst = std::min(worst, j);
            }
        }
    }
    return worst;
}

std::size_t ps_distinct_points(const Mesh& rMesh) {
    std::set<std::array<long long, 3>> seen;
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        const PsVec p = ps_point(rMesh, static_cast<std::int64_t>(i));
        seen.insert({std::llround(p[0] * 1e9), std::llround(p[1] * 1e9), std::llround(p[2] * 1e9)});
    }
    return seen.size();
}

std::size_t ps_region_size(const Mesh& rMesh, const std::string& rName,
                           meshioplusplus::RegionKind Kind) {
    const std::size_t i = rMesh.FindRegion(rName, Kind);
    EXPECT_NE(i, Mesh::npos) << rName;
    return i == Mesh::npos ? 0 : rMesh.Region(i).NumEntries();
}

constexpr double kPsPi = 3.14159265358979323846;

}  // namespace

TEST(PrimitiveShapes, ZeroCountsGiveAnEmptyMesh) {
    for (const Mesh& m :
         {mio::line(0), mio::circle(1.0, 0), mio::sphere_points(0), mio::box_points(0),
          mio::cube_sphere(1.0, 0), mio::geo_sphere(1.0, 0), mio::torus(1.0, 0.25, 0),
          mio::moebius(1.0, 0.5, 32, 0), mio::tube(1.0, 1.0, 16, 0), mio::ball(1.0, 0),
          mio::fichera(0), mio::solid_torus(1.0, 0.25, 4, 0)}) {
        EXPECT_EQ(m.NumPoints(), 0u);
        EXPECT_EQ(m.NumCellBlocks(), 0u);
    }
}

TEST(PrimitiveShapes, BadArgumentsRaiseByName) {
    auto raises = [](auto&& rCall, const std::string& rText) {
        try {
            rCall();
        } catch (const std::invalid_argument& rErr) {
            EXPECT_NE(std::string(rErr.what()).find(rText), std::string::npos) << rErr.what();
            return;
        }
        ADD_FAILURE() << "did not raise: " << rText;
    };
    raises([] { mio::circle(1.0, 2); }, "n must be at least 3");
    raises([] { mio::helix(1.0, 1.0, 2.0, 4); }, "n must exceed 2 * turns");
    raises([] { mio::torus(1.0, 2.0); }, "minor_radius must be less than");
    raises([] { mio::moebius(1.0, 3.0); }, "width must be less than twice");
    raises([] { mio::platonic("sphere"); }, "unknown name 'sphere'");
    raises([] { mio::spherical_shell(1.0, 0.5); }, "outer_radius must exceed");
    raises([] { mio::box({{1, 1, 1}}, {{1, 1, 1}}, {{0, 0, 0}}, "hexahedron", 11); },
           "order must be between 1 and 10");
    raises([] { mio::ball(1.0, 30, {{0, 0, 0}}, "hexahedron", 1, 50); }, "above the limit");
}

TEST(PrimitiveShapes, CurvesAndPointClouds) {
    const Mesh l = mio::line(4, {{1.0, 0.0, 0.0}}, {{1.0, 2.0, 4.0}});
    EXPECT_EQ(l.NumPoints(), 5u);
    EXPECT_EQ(ps_point(l, 4)[2], 4.0);
    EXPECT_EQ(ps_region_size(l, "start", meshioplusplus::RegionKind::Point), 1u);

    const Mesh c = mio::circle(2.0, 9, {{0.0, 0.0, 0.0}}, 2);
    EXPECT_EQ(c.Cells(0).Type(), "line3");
    for (std::size_t i = 0; i < c.NumPoints(); ++i)
        EXPECT_NEAR(ps_norm(ps_point(c, static_cast<std::int64_t>(i))), 2.0, 1e-12);

    const Mesh h = mio::helix(1.0, 0.5, 3.0, 30);
    EXPECT_NEAR(ps_point(h, 30)[0], 1.0, 1e-12);
    EXPECT_NEAR(ps_point(h, 30)[2], 1.5, 1e-12);

    const Mesh s = mio::sphere_points(200, 3.0);
    EXPECT_EQ(s.Cells(0).Type(), "vertex");
    for (std::size_t i = 0; i < s.NumPoints(); ++i)
        EXPECT_NEAR(ps_norm(ps_point(s, static_cast<std::int64_t>(i))), 3.0, 1e-12);

    // splitmix64 seeded with 0 first returns 0xE220A8397B1DCDAF.
    const Mesh b = mio::box_points(2);
    EXPECT_EQ(ps_point(b, 0)[0],
              static_cast<double>(0xE220A8397B1DCDAFULL >> 11) / 9007199254740992.0);
}

TEST(PrimitiveShapes, ClosedSurfacesAreOutwardWeldedAndOfTheRightGenus) {
    struct Case {
        Mesh mMesh;
        double mVolume;
        double mRtol;
        long mEuler;
    };
    const double ball = 4.0 / 3.0 * kPsPi;
    std::vector<Case> cases;
    cases.push_back({mio::tetra_sphere(), ball, 0.07, 2});
    cases.push_back({mio::octa_sphere(), ball, 0.1, 2});
    cases.push_back({mio::cube_sphere(1.0, 6), ball, 0.04, 2});
    cases.push_back({mio::geo_sphere(1.0, 5), ball, 0.03, 2});
    cases.push_back({mio::torus(), 2.0 * kPsPi * kPsPi * 0.0625, 0.04, 0});
    cases.push_back({mio::cone(1.0, 1.0, 64), kPsPi / 3.0, 0.01, 2});
    cases.push_back({mio::capsule(0.5, 1.0, 64, 4, 16), kPsPi * 0.25 + kPsPi / 6.0, 0.01, 2});
    cases.push_back({mio::box_surface({{2, 3, 4}}, {{1.0, 2.0, 3.0}}), 6.0, 1e-12, 2});
    cases.push_back({mio::platonic("dodecahedron"), 2.785163863, 1e-6, 2});
    for (const Case& c : cases) {
        EXPECT_EQ(ps_distinct_points(c.mMesh), c.mMesh.NumPoints());
        EXPECT_NEAR(ps_enclosed_volume(c.mMesh), c.mVolume, c.mRtol * c.mVolume);
        EXPECT_EQ(ps_euler(c.mMesh), c.mEuler);
    }
    EXPECT_EQ(ps_euler(mio::tube()), 0);
    EXPECT_EQ(ps_euler(mio::moebius()), 0);
}

TEST(PrimitiveShapes, RegionsNameTheBoundaries) {
    using K = meshioplusplus::RegionKind;
    const Mesh t = mio::tube(1.0, 1.0, 8, 3);
    EXPECT_EQ(ps_region_size(t, "bottom", K::Side), 8u);
    const Mesh c = mio::cone(1.0, 1.0, 8, 3);
    EXPECT_EQ(ps_region_size(c, "side", K::Cell) + ps_region_size(c, "base", K::Cell),
              ps_num_cells(c));
    const Mesh bs = mio::box_surface({{2, 3, 4}});
    EXPECT_EQ(ps_region_size(bs, "xmin", K::Cell), 12u);
    EXPECT_EQ(ps_region_size(bs, "zmax", K::Cell), 6u);
    EXPECT_EQ(ps_region_size(mio::moebius(1.0, 0.5, 10, 2), "boundary", K::Side), 20u);
    EXPECT_EQ(ps_region_size(mio::ball(1.0, 2), "surface", K::Side), 96u);
    const Mesh p = mio::pipe(0.5, 1.0, 1.0, 8, 2, 3);
    EXPECT_EQ(ps_region_size(p, "inner", K::Side), 24u);
    EXPECT_EQ(ps_region_size(p, "top", K::Side), 16u);
    const Mesh f = mio::fichera(2);
    EXPECT_EQ(ps_num_cells(f), 56u);
    EXPECT_EQ(ps_region_size(f, "notch_z", K::Side), 4u);
}

TEST(PrimitiveShapes, SolidsArePositiveAndWelded) {
    for (const Mesh& m :
         {mio::ball(1.0, 3), mio::ball(1.0, 2, {{0, 0, 0}}, "tetra"),
          mio::spherical_shell(0.5, 1.0, 4, 2, 2.0), mio::pipe(0.5, 1.0, 1.0, 16, 2, 4, 0.5),
          mio::fichera(2), mio::fichera(2, 2.0, {{0, 0, 0}}, "tetra"), mio::solid_torus()}) {
        EXPECT_GT(ps_min_jacobian(m), 0.0);
        EXPECT_EQ(ps_distinct_points(m), m.NumPoints());
    }
}

TEST(PrimitiveShapes, CompleteBuildsTheFullLagrangeTypes) {
    EXPECT_EQ(mio::rectangle({{2, 2}}, {{1, 1}}, {{0, 0, 0}}, "quad", "zigzag", 2, 0, true)
                  .Cells(0)
                  .Type(),
              "quad9");
    const Mesh hex = mio::box({{2, 2, 2}}, {{1, 1, 1}}, {{0, 0, 0}}, "hexahedron", 2, 0, true);
    EXPECT_EQ(hex.Cells(0).Type(), "hexahedron27");
    EXPECT_EQ(hex.NumPoints(), 125u);  // the shared 5^3 lattice
    EXPECT_EQ(mio::box({{1, 1, 1}}, {{1, 1, 1}}, {{0, 0, 0}}, "wedge", 2, 0, true).Cells(0).Type(),
              "wedge18");
    const Mesh cs = mio::cube_sphere(1.0, 2, {{0, 0, 0}}, "quad", 2, 0, true);
    for (std::size_t i = 0; i < cs.NumPoints(); ++i)
        EXPECT_NEAR(ps_norm(ps_point(cs, static_cast<std::int64_t>(i))), 1.0, 1e-12);
}

TEST(PrimitiveShapes, HigherOrdersBuildSharedVtkLagrangeCells) {
    for (int order : {3, 5}) {
        const Mesh g = mio::geo_sphere(1.0, 2, {{0, 0, 0}}, order);
        EXPECT_EQ(g.Cells(0).Type(), "VTK_LAGRANGE_TRIANGLE");
        const std::size_t faces = g.Cells(0).NumCells();
        EXPECT_EQ(g.NumPoints(), faces * order * order / 2 + 2);
        for (std::size_t i = 0; i < g.NumPoints(); ++i)
            EXPECT_NEAR(ps_norm(ps_point(g, static_cast<std::int64_t>(i))), 1.0, 1e-12);
        const Mesh b = mio::ball(1.0, 1, {{0, 0, 0}}, "hexahedron", order);
        EXPECT_EQ(b.Cells(0).Type(), "VTK_LAGRANGE_HEXAHEDRON");
        EXPECT_EQ(b.NumPoints(),
                  static_cast<std::size_t>((2 * order + 1) * (2 * order + 1) * (2 * order + 1)));
    }
    // Side regions survive: found on the linear cells, the same facets.
    const Mesh linear = mio::box({{2, 2, 2}});
    const Mesh cubic = mio::box({{2, 2, 2}}, {{1, 1, 1}}, {{0, 0, 0}}, "hexahedron", 3, 0, false);
    using K = meshioplusplus::RegionKind;
    EXPECT_EQ(ps_region_size(linear, "zmax", K::Side), ps_region_size(cubic, "zmax", K::Side));
    EXPECT_EQ(mio::box({{1, 1, 1}}, {{1, 1, 1}}, {{0, 0, 0}}, "tetra", 4).Cells(0).NodesPerCell(),
              35u);
    EXPECT_EQ(mio::line(3, {{0, 0, 0}}, {{1, 0, 0}}, 4).Cells(0).Type(), "VTK_LAGRANGE_CURVE");
}
