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

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/operations/remesh.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

/// A closed icosahedron -- the smallest genuinely 3D closed triangulation, and
/// unlike a cube it has no coplanar faces to make a clustering artificially
/// easy.
Mesh icosahedron() {
    const double t = (1.0 + std::sqrt(5.0)) / 2.0;
    std::vector<std::vector<double>> p = {
        {-1, t, 0}, {1, t, 0},  {-1, -t, 0}, {1, -t, 0}, {0, -1, t},  {0, 1, t},
        {0, -1, -t}, {0, 1, -t}, {t, 0, -1},  {t, 0, 1},  {-t, 0, -1}, {-t, 0, 1}};
    for (auto& v : p) {
        const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        for (auto& x : v)
            x /= n;
    }
    std::vector<std::vector<std::int64_t>> f = {
        {0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11}, {1, 5, 9},  {5, 11, 4},
        {11, 10, 2}, {10, 7, 6}, {7, 1, 8},  {3, 9, 4},  {3, 4, 2},   {3, 2, 6},  {3, 6, 8},
        {3, 8, 9},  {4, 9, 5},  {2, 4, 11},  {6, 2, 10}, {8, 6, 7},   {9, 8, 1}};
    return mt::make_mesh(p, "triangle", f);
}

/// A triangulated unit cube's boundary, welded at the seams (a genuine
/// closed 2-manifold, unlike six independent per-face grids) -- deliberately
/// sharp-cornered, for the feature-preservation tests: quadric mode should
/// keep output vertices near the true cube surface, isotropic mode should
/// round the corners. Vertices are quantized (0..per_edge)^3 lattice points
/// restricted to the boundary shell, so two faces sharing an edge share the
/// exact same point ids -- the adjacency a cluster needs to cross a seam.
Mesh cube_mesh(int per_edge) {
    const int n = per_edge;
    std::vector<std::vector<double>> p;
    std::map<std::array<int, 3>, std::int64_t> id_of;
    auto vid = [&](int i, int j, int k) -> std::int64_t {
        const std::array<int, 3> key = {i, j, k};
        auto it = id_of.find(key);
        if (it != id_of.end())
            return it->second;
        const std::int64_t id = static_cast<std::int64_t>(p.size());
        p.push_back({static_cast<double>(i) / n, static_cast<double>(j) / n,
                    static_cast<double>(k) / n});
        id_of.emplace(key, id);
        return id;
    };
    std::vector<std::vector<std::int64_t>> f;
    auto add_face = [&](auto grid_to_xyz) {
        // grid_to_xyz(a, b) -> (i, j, k) lattice coords on this face, for
        // a, b in [0, n].
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b) {
                const auto [i0, j0, k0] = grid_to_xyz(a, b);
                const auto [i1, j1, k1] = grid_to_xyz(a + 1, b);
                const auto [i2, j2, k2] = grid_to_xyz(a + 1, b + 1);
                const auto [i3, j3, k3] = grid_to_xyz(a, b + 1);
                const std::int64_t v0 = vid(i0, j0, k0), v1 = vid(i1, j1, k1),
                                   v2 = vid(i2, j2, k2), v3 = vid(i3, j3, k3);
                f.push_back({v0, v1, v2});
                f.push_back({v0, v2, v3});
            }
    };
    add_face([&](int a, int b) { return std::array<int, 3>{a, b, 0}; });
    add_face([&](int a, int b) { return std::array<int, 3>{a, b, n}; });
    add_face([&](int a, int b) { return std::array<int, 3>{a, 0, b}; });
    add_face([&](int a, int b) { return std::array<int, 3>{a, n, b}; });
    add_face([&](int a, int b) { return std::array<int, 3>{0, a, b}; });
    add_face([&](int a, int b) { return std::array<int, 3>{n, a, b}; });
    return mt::make_mesh(p, "triangle", f);
}

/// Every undirected edge's use count. A closed 2-manifold triangulation uses
/// each edge exactly twice; this is the watertightness oracle.
std::map<std::pair<std::int64_t, std::int64_t>, int> edge_uses(const Mesh& rM) {
    std::map<std::pair<std::int64_t, std::int64_t>, int> uses;
    for (const auto cb : rM.CellRange()) {
        // Boundary-preserving remesh may attach a second, 2-noded `line`
        // block alongside the triangles; this oracle is triangle-only.
        if (cb.Type() != "triangle")
            continue;
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            for (int k = 0; k < 3; ++k) {
                std::int64_t a = detail::read_int(conn, c * 3 + k);
                std::int64_t b = detail::read_int(conn, c * 3 + (k + 1) % 3);
                if (a > b)
                    std::swap(a, b);
                ++uses[{a, b}];
            }
    }
    return uses;
}

/// Smallest triangle angle in the mesh, in degrees -- the quality figure this
/// whole operation exists to raise.
double min_angle_deg(const Mesh& rM) {
    double worst = 180.0;
    const NDArray& pts = rM.Points();
    const std::size_t dim = rM.PointDim();
    auto pt = [&](std::int64_t i, int d) {
        return detail::read_double(pts, static_cast<std::size_t>(i) * dim + d);
    };
    for (const auto cb : rM.CellRange()) {
        // Same reasoning as edge_uses: skip the optional boundary `line`
        // block, which has 2 nodes per cell rather than 3.
        if (cb.Type() != "triangle")
            continue;
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::array<std::int64_t, 3> v = {detail::read_int(conn, c * 3),
                                             detail::read_int(conn, c * 3 + 1),
                                             detail::read_int(conn, c * 3 + 2)};
            for (int k = 0; k < 3; ++k) {
                double e0[3], e1[3];
                for (int d = 0; d < 3; ++d) {
                    e0[d] = pt(v[(k + 1) % 3], d) - pt(v[k], d);
                    e1[d] = pt(v[(k + 2) % 3], d) - pt(v[k], d);
                }
                const double n0 = std::sqrt(e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2]);
                const double n1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
                if (n0 <= 0.0 || n1 <= 0.0)
                    continue;
                double dot = (e0[0] * e1[0] + e0[1] * e1[1] + e0[2] * e1[2]) / (n0 * n1);
                dot = std::max(-1.0, std::min(1.0, dot));
                worst = std::min(worst, std::acos(dot) * 180.0 / M_PI);
            }
        }
    }
    return worst;
}

/// A unit square triangulated into deliberately awful slivers: every triangle
/// shares the same two base vertices, so the strip is a fan of needles. This
/// is the input QEM decimation cannot help -- it can only delete needles,
/// never replace them.
Mesh sliver_strip(int n) {
    std::vector<std::vector<double>> p;
    for (int i = 0; i <= n; ++i)
        p.push_back({static_cast<double>(i) / n, 0.0, 0.0});
    for (int i = 0; i <= n; ++i)
        p.push_back({static_cast<double>(i) / n, 1.0, 0.0});
    std::vector<std::vector<std::int64_t>> f;
    for (int i = 0; i < n; ++i) {
        f.push_back({i, i + 1, n + 1 + i});
        f.push_back({i + 1, n + 2 + i, n + 1 + i});
    }
    return mt::make_mesh(p, "triangle", f);
}

/// Distance from `(x, y, z)` to the nearest point on the true unit cube's
/// surface (min over the 6 planes' clamped projection) -- the sharp-edge
/// oracle: a rounded corner reports a larger distance than a preserved one.
double dist_to_cube_surface(double x, double y, double z) {
    const double lo[3] = {0.0, 0.0, 0.0}, hi[3] = {1.0, 1.0, 1.0};
    double p[3] = {x, y, z};
    double best = 1e300;
    for (int axis = 0; axis < 3; ++axis)
        for (double target : {lo[axis], hi[axis]}) {
            double q[3] = {p[0], p[1], p[2]};
            q[axis] = target;
            for (int d = 0; d < 3; ++d)
                if (d != axis)
                    q[d] = std::max(lo[d], std::min(hi[d], q[d]));
            const double dx = q[0] - p[0], dy = q[1] - p[1], dz = q[2] - p[2];
            best = std::min(best, std::sqrt(dx * dx + dy * dy + dz * dz));
        }
    return best;
}

double max_deviation_from_cube(const Mesh& rM) {
    double worst = 0.0;
    const NDArray& pts = rM.Points();
    for (std::size_t i = 0; i < rM.NumPoints(); ++i)
        worst = std::max(worst, dist_to_cube_surface(detail::read_double(pts, i * 3),
                                                      detail::read_double(pts, i * 3 + 1),
                                                      detail::read_double(pts, i * 3 + 2)));
    return worst;
}

}  // namespace

TEST(Remesh, ProducesTheRequestedNumberOfVerticesOnASphere) {
    RemeshOptions o;
    o.mNumClusters = 200;
    const RemeshResult r = remesh(icosahedron(), o);

    EXPECT_EQ(r.mNumClusters, 200);
    EXPECT_EQ(static_cast<std::int64_t>(r.mMesh.NumPoints()), r.mNumClusters);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "triangle");
    EXPECT_GT(r.mMesh.Cells(0).NumCells(), 0u);
    // Subdivision must have kicked in on its own: 12 vertices cannot support
    // 200 clusters, so an unsubdivided run would have thrown.
    EXPECT_GT(r.mSubdivideApplied, 0);
}

TEST(Remesh, OutputOfAClosedSurfaceIsItselfClosed) {
    RemeshOptions o;
    o.mNumClusters = 150;
    const RemeshResult r = remesh(icosahedron(), o);

    ASSERT_EQ(r.mNumIsolatedClusters, 0);
    for (const auto& [edge, uses] : edge_uses(r.mMesh))
        EXPECT_EQ(uses, 2) << "edge (" << edge.first << ", " << edge.second << ") used " << uses
                           << " times; the dual of the clustering is not watertight";
}

TEST(Remesh, EulerCharacteristicOfTheSphereIsPreserved) {
    RemeshOptions o;
    o.mNumClusters = 150;
    const RemeshResult r = remesh(icosahedron(), o);

    const std::int64_t v = static_cast<std::int64_t>(r.mMesh.NumPoints());
    const std::int64_t f = static_cast<std::int64_t>(r.mMesh.Cells(0).NumCells());
    const std::int64_t e = static_cast<std::int64_t>(edge_uses(r.mMesh).size());
    // A genus-0 closed surface has V - E + F == 2. This catches a dual that is
    // watertight but has grown a handle -- something edge counting alone
    // cannot see.
    EXPECT_EQ(v - e + f, 2);
}

TEST(Remesh, RaisesQualityOfASliverTriangulation) {
    const Mesh in = sliver_strip(60);
    const double before = min_angle_deg(in);
    ASSERT_LT(before, 5.0) << "the fixture is supposed to be awful";

    RemeshOptions o;
    o.mNumClusters = 60;
    const RemeshResult r = remesh(in, o);
    const double after = min_angle_deg(r.mMesh);

    // The whole justification for this operation over decimate: element shape
    // is a property of the clustering, not of the input.
    EXPECT_GT(after, before * 3.0)
        << "min angle went from " << before << " to " << after << " degrees";
}

TEST(Remesh, IsDeterministic) {
    RemeshOptions o;
    o.mNumClusters = 120;
    const RemeshResult a = remesh(icosahedron(), o);
    const RemeshResult b = remesh(icosahedron(), o);

    ASSERT_EQ(a.mMesh.NumPoints(), b.mMesh.NumPoints());
    ASSERT_EQ(a.mMesh.Cells(0).NumCells(), b.mMesh.Cells(0).NumCells());
    EXPECT_EQ(a.mNumIterations, b.mNumIterations);
    const NDArray& pa = a.mMesh.Points();
    const NDArray& pb = b.mMesh.Points();
    for (std::size_t i = 0; i < a.mMesh.NumPoints() * 3; ++i)
        EXPECT_EQ(detail::read_double(pa, i), detail::read_double(pb, i))
            << "coordinate " << i << " differs between two identical runs";
}

TEST(Remesh, SubdivideZeroSkipsRefinement) {
    RemeshOptions o;
    o.mNumClusters = 10;
    o.mSubdivide = 0;
    const RemeshResult r = remesh(icosahedron(), o);
    EXPECT_EQ(r.mSubdivideApplied, 0);
    EXPECT_LE(r.mNumClusters, 10);
}

TEST(Remesh, TriangulatesQuadInputRatherThanRefusingIt) {
    // A quad grid must go through convert_cells(Simplexify) the way decimate's
    // does, not throw.
    std::vector<std::vector<double>> p;
    for (int j = 0; j <= 6; ++j)
        for (int i = 0; i <= 6; ++i)
            p.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});
    std::vector<std::vector<std::int64_t>> f;
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i) {
            const std::int64_t b = j * 7 + i;
            f.push_back({b, b + 1, b + 8, b + 7});
        }
    RemeshOptions o;
    o.mNumClusters = 20;
    const RemeshResult r = remesh(mt::make_mesh(p, "quad", f), o);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "triangle");
    EXPECT_GT(r.mMesh.Cells(0).NumCells(), 0u);
}

TEST(Remesh, RejectsVolumeAndHigherOrderInputByName) {
    RemeshOptions o;
    o.mNumClusters = 10;
    EXPECT_THROW(
        {
            try {
                remesh(mt::tet_mesh(), o);
            } catch (const std::invalid_argument& e) {
                EXPECT_NE(std::string(e.what()).find("extract_surface"), std::string::npos);
                throw;
            }
        },
        std::invalid_argument);
    EXPECT_THROW(
        {
            try {
                remesh(mt::triangle6_mesh(), o);
            } catch (const std::invalid_argument& e) {
                EXPECT_NE(std::string(e.what()).find("Linearize"), std::string::npos);
                throw;
            }
        },
        std::invalid_argument);
}

TEST(Remesh, RejectsATooSmallClusterCount) {
    RemeshOptions o;
    o.mNumClusters = 3;
    EXPECT_THROW(remesh(icosahedron(), o), std::invalid_argument);
}

TEST(Remesh, CarriesFieldDataAndDropsPointData) {
    Mesh in = icosahedron();
    in.AddFieldData("solver", mt::data_array({1.5}));
    in.AddPointData("temperature", mt::data_array(std::vector<double>(12, 3.0)));

    RemeshOptions o;
    o.mNumClusters = 40;
    const RemeshResult r = remesh(in, o);

    EXPECT_TRUE(r.mMesh.HasFieldData("solver"));
    // Point data has no meaning on new vertices; the documented answer is
    // interpolate, not a silent guess.
    EXPECT_FALSE(r.mMesh.HasPointData("temperature"));
}

TEST(Remesh, MetricNameParsingRoundTrips) {
    EXPECT_EQ(remesh_metric_from_name("isotropic"), RemeshMetric::Isotropic);
    EXPECT_EQ(remesh_metric_from_name("quadric"), RemeshMetric::Quadric);
    EXPECT_EQ(remesh_metric_from_name("anisotropic"), RemeshMetric::Anisotropic);
    EXPECT_THROW(remesh_metric_from_name("bogus"), std::invalid_argument);
}

TEST(Remesh, QuadricModeIsDeterministic) {
    RemeshOptions o;
    o.mNumClusters = 80;
    o.mMetric = RemeshMetric::Quadric;
    const RemeshResult a = remesh(icosahedron(), o);
    const RemeshResult b = remesh(icosahedron(), o);

    ASSERT_EQ(a.mMesh.NumPoints(), b.mMesh.NumPoints());
    const NDArray& pa = a.mMesh.Points();
    const NDArray& pb = b.mMesh.Points();
    for (std::size_t i = 0; i < a.mMesh.NumPoints() * 3; ++i)
        EXPECT_EQ(detail::read_double(pa, i), detail::read_double(pb, i));
}

TEST(Remesh, QuadricModeProducesAValidWatertightMesh) {
    RemeshOptions o;
    o.mNumClusters = 150;
    o.mMetric = RemeshMetric::Quadric;
    const RemeshResult r = remesh(icosahedron(), o);

    for (const auto& [edge, uses] : edge_uses(r.mMesh))
        EXPECT_EQ(uses, 2);
    const std::int64_t v = static_cast<std::int64_t>(r.mMesh.NumPoints());
    const std::int64_t f = static_cast<std::int64_t>(r.mMesh.Cells(0).NumCells());
    const std::int64_t e = static_cast<std::int64_t>(edge_uses(r.mMesh).size());
    EXPECT_EQ(v - e + f, 2);
}

TEST(Remesh, QuadricMetricPreservesSharpCornersBetterThanIsotropic) {
    // per_edge=16 / 150 clusters is a stable, well-resolved point measured by
    // a resolution/cluster-count sweep (per_edge in {6,10,16}, clusters in
    // {60,150,300}): quadric mode's max deviation was consistently ~0.56-0.86x
    // isotropic's at every non-degenerate point in that sweep (a coarser mesh
    // has too few clusters per corner and both metrics land exactly on 0 or
    // near-0 by coincidence, which is not a meaningful comparison).
    const Mesh cube = cube_mesh(16);

    RemeshOptions iso;
    iso.mNumClusters = 150;
    iso.mSubdivide = 0;  // the fixture is already fine enough; keep sizes comparable
    iso.mMetric = RemeshMetric::Isotropic;
    const double iso_dev = max_deviation_from_cube(remesh(cube, iso).mMesh);

    RemeshOptions quad = iso;
    quad.mMetric = RemeshMetric::Quadric;
    const double quad_dev = max_deviation_from_cube(remesh(cube, quad).mMesh);

    // Isotropic clustering rounds corners inward (a corner's centroid is
    // pulled toward the cube's interior); quadric clustering's non-degenerate
    // corner/edge quadrics pin the representative point onto the true
    // surface. This is the operation's entire reason to have two metrics. The
    // margin (0.7, measured 0.556) is deliberately not tight -- it asserts a
    // genuine, reproducible improvement without pinning the sweep's exact
    // ratio, which would be fragile to unrelated tie-break changes.
    EXPECT_LT(quad_dev, iso_dev * 0.7)
        << "quadric max deviation " << quad_dev << " vs isotropic " << iso_dev;
}

/// A Gaussian bump on a flat square: curvature is large near the centre and
/// ~0 near the flat edges -- a spatially *varying* curvature fixture the
/// constant-curvature icosahedron/cube fixtures above cannot discriminate.
Mesh gaussian_bump(int n, double half, double amplitude, double sigma) {
    std::vector<std::vector<double>> p;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            const double x = -half + 2.0 * half * i / n;
            const double y = -half + 2.0 * half * j / n;
            const double r2 = x * x + y * y;
            p.push_back({x, y, amplitude * std::exp(-r2 / (sigma * sigma))});
        }
    std::vector<std::vector<std::int64_t>> f;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const std::int64_t b = j * (n + 1) + i;
            f.push_back({b, b + 1, b + n + 2});
            f.push_back({b, b + n + 2, b + n + 1});
        }
    return mt::make_mesh(p, "triangle", f);
}

/// An open cylindrical tube (no end caps): curvature is 1/radius around the
/// circumference and exactly 0 along the axis -- the canonical anisotropic
/// surface, since a sphere or bump cannot separate "curvature exists" from
/// "curvature is direction-dependent" the way this can. The two rims are
/// left as open boundaries; tests using this disable mPreserveBoundary so
/// boundary pinning cannot be mistaken for the metric's own effect.
Mesh cylinder_mesh(int n_circ, int n_axial, double radius, double height) {
    std::vector<std::vector<double>> p;
    for (int j = 0; j <= n_axial; ++j)
        for (int i = 0; i < n_circ; ++i) {
            const double theta = 2.0 * M_PI * i / n_circ;
            p.push_back(
                {radius * std::cos(theta), radius * std::sin(theta), height * j / n_axial});
        }
    std::vector<std::vector<std::int64_t>> f;
    auto vid = [&](int i, int j) { return j * n_circ + (i % n_circ); };
    for (int j = 0; j < n_axial; ++j)
        for (int i = 0; i < n_circ; ++i) {
            const std::int64_t v0 = vid(i, j), v1 = vid(i + 1, j), v2 = vid(i + 1, j + 1),
                               v3 = vid(i, j + 1);
            f.push_back({v0, v1, v2});
            f.push_back({v0, v2, v3});
        }
    return mt::make_mesh(p, "triangle", f);
}

/// The mean output-edge axial extent divided by the mean output-edge
/// circumferential (arc-length) extent, over a mesh produced from
/// `cylinder_mesh`. Isotropic clustering has no reason to prefer one
/// direction over the other, so this sits near 1; an anisotropic metric
/// that correctly reads "curved around, flat along the axis" should raise
/// it well above 1 -- short edges bridging the curvature, long edges
/// running along the flat direction.
double axial_over_circumferential_ratio(const Mesh& rM, double radius) {
    const NDArray& pts = rM.Points();
    const NDArray& conn = rM.Cells(0).Conn();
    double sum_axial = 0.0, sum_arc = 0.0;
    for (std::size_t c = 0; c < rM.Cells(0).NumCells(); ++c) {
        const std::int64_t tri[3] = {detail::read_int(conn, c * 3),
                                     detail::read_int(conn, c * 3 + 1),
                                     detail::read_int(conn, c * 3 + 2)};
        for (int k = 0; k < 3; ++k) {
            const std::int64_t a = tri[k], b = tri[(k + 1) % 3];
            const double xa = detail::read_double(pts, static_cast<std::size_t>(a) * 3);
            const double ya = detail::read_double(pts, static_cast<std::size_t>(a) * 3 + 1);
            const double za = detail::read_double(pts, static_cast<std::size_t>(a) * 3 + 2);
            const double xb = detail::read_double(pts, static_cast<std::size_t>(b) * 3);
            const double yb = detail::read_double(pts, static_cast<std::size_t>(b) * 3 + 1);
            const double zb = detail::read_double(pts, static_cast<std::size_t>(b) * 3 + 2);
            double dtheta = std::atan2(yb, xb) - std::atan2(ya, xa);
            while (dtheta > M_PI)
                dtheta -= 2.0 * M_PI;
            while (dtheta < -M_PI)
                dtheta += 2.0 * M_PI;
            sum_arc += std::abs(dtheta) * radius;
            sum_axial += std::abs(za - zb);
        }
    }
    return sum_axial / sum_arc;
}

TEST(Remesh, CurvatureGradationConcentratesClustersNearHighCurvature) {
    const Mesh bump = gaussian_bump(40, 1.0, 0.6, 0.25);

    auto count_near_center = [](const Mesh& m, double radius) {
        const NDArray& pts = m.Points();
        std::size_t count = 0;
        for (std::size_t i = 0; i < m.NumPoints(); ++i) {
            const double x = detail::read_double(pts, i * 3);
            const double y = detail::read_double(pts, i * 3 + 1);
            if (x * x + y * y < radius * radius)
                ++count;
        }
        return count;
    };

    RemeshOptions flat;
    flat.mNumClusters = 150;
    flat.mSubdivide = 0;  // the grid is already dense enough; keep positions exact
    flat.mGradation = 0.0;
    const RemeshResult uniform = remesh(bump, flat);

    RemeshOptions graded = flat;
    graded.mGradation = 2.0;
    const RemeshResult curved = remesh(bump, graded);

    // A fixed central disk covering the bulk of the bump's curvature; with
    // gradation weighting each item by area * kappa^gamma, more of the same
    // fixed cluster budget lands inside it than under plain area weighting.
    const double radius = 0.4;
    const std::size_t n_uniform = count_near_center(uniform.mMesh, radius);
    const std::size_t n_curved = count_near_center(curved.mMesh, radius);
    EXPECT_GT(n_curved, n_uniform)
        << "gradation did not concentrate clusters near the curved region: " << n_curved
        << " vs " << n_uniform << " (uniform)";
}

TEST(Remesh, PreservesBoundaryOfAnOpenPatch) {
    // A flat, open square patch: unlike a closed surface, remesh must not
    // silently drop the boundary into whichever interior cluster reaches it
    // first.
    const int n = 24;
    const double half = 1.0;
    std::vector<std::vector<double>> p;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i)
            p.push_back({-half + 2.0 * half * i / n, -half + 2.0 * half * j / n, 0.0});
    std::vector<std::vector<std::int64_t>> f;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const std::int64_t b = j * (n + 1) + i;
            f.push_back({b, b + 1, b + n + 2});
            f.push_back({b, b + n + 2, b + n + 1});
        }
    const Mesh patch = mt::make_mesh(p, "triangle", f);

    RemeshOptions o;
    o.mNumClusters = 80;
    o.mSubdivide = 0;
    const RemeshResult r = remesh(patch, o);

    ASSERT_EQ(r.mMesh.NumCellBlocks(), 2u)
        << "an open input should leave a second, boundary `line` block behind";
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "triangle");
    ASSERT_EQ(r.mMesh.Cells(1).Type(), "line");
    ASSERT_GT(r.mMesh.Cells(1).NumCells(), 0u);

    // The dual boundary polyline should roughly track the true perimeter
    // (4 * 2 * half = 8), not be empty or wildly larger.
    const NDArray& lconn = r.mMesh.Cells(1).Conn();
    const NDArray& pts = r.mMesh.Points();
    double total_length = 0.0;
    for (std::size_t c = 0; c < r.mMesh.Cells(1).NumCells(); ++c) {
        const std::int64_t a = detail::read_int(lconn, c * 2);
        const std::int64_t b = detail::read_int(lconn, c * 2 + 1);
        double d2 = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double diff = detail::read_double(pts, static_cast<std::size_t>(a) * 3 + k) -
                                detail::read_double(pts, static_cast<std::size_t>(b) * 3 + k);
            d2 += diff * diff;
        }
        total_length += std::sqrt(d2);
    }
    const double true_perimeter = 4.0 * 2.0 * half;
    EXPECT_GT(total_length, true_perimeter * 0.5)
        << "boundary polyline length " << total_length << " vs perimeter " << true_perimeter;
    EXPECT_LT(total_length, true_perimeter * 1.5)
        << "boundary polyline length " << total_length << " vs perimeter " << true_perimeter;
}

TEST(Remesh, AnisotropicMetricElongatesClustersAlongTheLowCurvatureAxis) {
    const double radius = 1.0, height = 8.0;
    const Mesh cyl = cylinder_mesh(24, 48, radius, height);

    RemeshOptions iso;
    iso.mNumClusters = 200;
    iso.mSubdivide = 0;
    iso.mPreserveBoundary = false;  // isolate the metric, not boundary pinning
    const double iso_ratio = axial_over_circumferential_ratio(remesh(cyl, iso).mMesh, radius);

    RemeshOptions aniso = iso;
    aniso.mMetric = RemeshMetric::Anisotropic;
    aniso.mMaxAnisotropy = 8.0;
    const double aniso_ratio = axial_over_circumferential_ratio(remesh(cyl, aniso).mMesh, radius);

    // max_anisotropy == 1.0 clamps the shape tensor to the identity
    // everywhere (see remesh.hpp), so it should show no meaningful
    // elongation over plain isotropic clustering -- this isolates "the
    // metric is doing something" from "the clamp actually clamps".
    RemeshOptions aniso_clamped = aniso;
    aniso_clamped.mMaxAnisotropy = 1.0;
    const double clamped_ratio =
        axial_over_circumferential_ratio(remesh(cyl, aniso_clamped).mMesh, radius);

    EXPECT_GT(aniso_ratio, iso_ratio * 1.3)
        << "anisotropic axial/circumferential ratio " << aniso_ratio << " vs isotropic "
        << iso_ratio;
    EXPECT_LT(clamped_ratio, iso_ratio * 1.3)
        << "max_anisotropy=1.0 still elongated: " << clamped_ratio << " vs isotropic "
        << iso_ratio;
}

TEST(Remesh, AnisotropicMetricStaysExactlyFlatOnAFlatSurface) {
    // amplitude 0.0 -> z == 0 identically at every input point, so a
    // correctly-signed, correctly-indexed curvature tensor must reduce to
    // the identity everywhere (hi_kappa == 0 unconditionally) and the
    // z-row of the resulting quadric's linear term must vanish exactly.
    // NOTE this is a regression guard, not a comprehensive oracle for the
    // tensor's normal-pinning term specifically: on a flat, fully symmetric
    // input M == I collapses several distinct bugs (e.g. dropping the
    // tensor's n*n^T term altogether) onto the SAME z == 0 answer, since
    // the quadric solve's ill-conditioning fallback also lands on z == 0
    // here. See AnisotropicMetricKeepsClustersOnACurvedSurface below for the
    // oracle that actually discriminates the normal-pinning term.
    const Mesh flat = gaussian_bump(20, 1.0, 0.0, 0.25);

    RemeshOptions o;
    o.mNumClusters = 60;
    o.mSubdivide = 0;
    o.mMetric = RemeshMetric::Anisotropic;
    const RemeshResult r = remesh(flat, o);

    const NDArray& pts = r.mMesh.Points();
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
        EXPECT_EQ(detail::read_double(pts, i * 3 + 2), 0.0)
            << "vertex " << i << " left the z=0 plane under the anisotropic metric";
}

TEST(Remesh, AnisotropicMetricKeepsClustersOnACurvedSurface) {
    // The tensor's normal eigenvalue is deliberately pinned to the sharper
    // in-plane one (see remesh.hpp) specifically so a cluster elongated
    // in-plane does not also drift off the true surface. A bug dropping
    // that term (e.g. missing the tensor's n*n^T contribution) degenerates
    // the per-cluster quadric solve on any near-flat member and falls back
    // to the plain centroid there -- which is systematically LOWER than the
    // bump's own local height away from its peak, pulling the output
    // visibly under the cap.
    const double amplitude = 0.5, sigma = 0.3;
    const Mesh bump = gaussian_bump(30, 1.0, amplitude, sigma);

    auto max_vertical_deviation = [](const Mesh& m, double amp, double s) {
        const NDArray& pts = m.Points();
        double worst = 0.0;
        for (std::size_t i = 0; i < m.NumPoints(); ++i) {
            const double x = detail::read_double(pts, i * 3);
            const double y = detail::read_double(pts, i * 3 + 1);
            const double z = detail::read_double(pts, i * 3 + 2);
            const double target = amp * std::exp(-(x * x + y * y) / (s * s));
            worst = std::max(worst, std::abs(z - target));
        }
        return worst;
    };

    RemeshOptions o;
    o.mNumClusters = 150;
    o.mSubdivide = 0;
    o.mMetric = RemeshMetric::Anisotropic;
    o.mMaxAnisotropy = 6.0;
    const double dev = max_vertical_deviation(remesh(bump, o).mMesh, amplitude, sigma);

    // Bound has real headroom (measured well under this on a correct
    // implementation) but is tight enough that a broken normal pin, which
    // pulls whole clusters toward the surrounding near-flat height instead
    // of the cap's own local height, blows through it -- verified by
    // sabotage (dropping the tensor's n*n^T term).
    EXPECT_LT(dev, amplitude * 0.5)
        << "max vertical deviation from the true bump surface: " << dev;
}

TEST(Remesh, AnisotropicModeIsDeterministic) {
    RemeshOptions o;
    o.mNumClusters = 120;
    o.mMetric = RemeshMetric::Anisotropic;
    const RemeshResult a = remesh(icosahedron(), o);
    const RemeshResult b = remesh(icosahedron(), o);

    ASSERT_EQ(a.mMesh.NumPoints(), b.mMesh.NumPoints());
    EXPECT_EQ(a.mNumIterations, b.mNumIterations);
    const NDArray& pa = a.mMesh.Points();
    const NDArray& pb = b.mMesh.Points();
    for (std::size_t i = 0; i < a.mMesh.NumPoints() * 3; ++i)
        EXPECT_EQ(detail::read_double(pa, i), detail::read_double(pb, i))
            << "coordinate " << i << " differs between two identical runs";
}

TEST(Remesh, AnisotropicModeProducesAValidWatertightMesh) {
    RemeshOptions o;
    o.mNumClusters = 150;
    o.mMetric = RemeshMetric::Anisotropic;
    const RemeshResult r = remesh(icosahedron(), o);

    for (const auto& [edge, uses] : edge_uses(r.mMesh))
        EXPECT_EQ(uses, 2);
    const std::int64_t v = static_cast<std::int64_t>(r.mMesh.NumPoints());
    const std::int64_t f = static_cast<std::int64_t>(r.mMesh.Cells(0).NumCells());
    const std::int64_t e = static_cast<std::int64_t>(edge_uses(r.mMesh).size());
    EXPECT_EQ(v - e + f, 2);
}

TEST(Remesh, MaxAnisotropyValidation) {
    RemeshOptions below_one;
    below_one.mNumClusters = 60;
    below_one.mMetric = RemeshMetric::Anisotropic;
    below_one.mMaxAnisotropy = 0.5;
    EXPECT_THROW(remesh(icosahedron(), below_one), std::invalid_argument);

    RemeshOptions off_default;
    off_default.mNumClusters = 60;
    off_default.mMetric = RemeshMetric::Isotropic;
    off_default.mMaxAnisotropy = 5.0;
    EXPECT_THROW(remesh(icosahedron(), off_default), std::invalid_argument);

    RemeshOptions off_default_quadric = off_default;
    off_default_quadric.mMetric = RemeshMetric::Quadric;
    EXPECT_THROW(remesh(icosahedron(), off_default_quadric), std::invalid_argument);

    // The default value is always accepted, regardless of metric.
    RemeshOptions ok = off_default;
    ok.mMaxAnisotropy = kRemeshDefaultMaxAnisotropy;
    EXPECT_NO_THROW(remesh(icosahedron(), ok));
}
