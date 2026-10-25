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

// mesh_cad_surface sized by a field (roadmap 6.3: curvature, local feature
// size, grading), with graded poles (6.4.2) and folds undone (6.4.3).

// System includes
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

// Project includes
#include "brep_solids.hpp"
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "../../src/cpp/src/detail/mesh_cad_size.hpp"
#include "../../src/cpp/src/operations/mesh_cad_testing.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::MeshCadSurfaceOptions;
using Vec3 = std::array<double, 3>;
constexpr double kPi = 3.14159265358979323846;

MeshCadSurfaceOptions sized(double H) {
    MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    return o;
}

std::uint64_t digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

Vec3 point(const mio::Mesh& rM, std::int64_t I) {
    const double* x = rM.Points().As<double>();
    return {x[3 * I], x[3 * I + 1], x[3 * I + 2]};
}

double dist(const Vec3& rA, const Vec3& rB) {
    return std::sqrt((rA[0] - rB[0]) * (rA[0] - rB[0]) + (rA[1] - rB[1]) * (rA[1] - rB[1]) +
                     (rA[2] - rB[2]) * (rA[2] - rB[2]));
}

// Every triangle of the first block, as point ids.
std::vector<std::array<std::int64_t, 3>> triangles(const mio::Mesh& rM) {
    std::vector<std::array<std::int64_t, 3>> out;
    const auto cells = rM.Cells(0);
    const mio::NDArray& c = cells.Conn();
    for (std::size_t t = 0; t < cells.NumCells(); ++t)
        out.push_back({mio::detail::read_int(c, 3 * t), mio::detail::read_int(c, 3 * t + 1),
                       mio::detail::read_int(c, 3 * t + 2)});
    return out;
}

// The largest distance from a triangle edge's midpoint to the cylinder of
// radius R about the z axis, over the edges with both ends on it.
double cylinder_chord_error(const mio::Mesh& rM, double R) {
    double worst = 0.0;
    for (const auto& r_t : triangles(rM))
        for (int k = 0; k < 3; ++k) {
            const Vec3 a = point(rM, r_t[k]), b = point(rM, r_t[(k + 1) % 3]);
            if (std::abs(std::hypot(a[0], a[1]) - R) > 1e-9 ||
                std::abs(std::hypot(b[0], b[1]) - R) > 1e-9)
                continue;
            const Vec3 m{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
            worst = std::max(worst, R - std::hypot(m[0], m[1]));
        }
    return worst;
}

// A band of the torus of radii R and r about the z axis, between the tube's
// parallels at -Phi and Phi: a narrow, strongly curved face when r is below h.
mio::BrepInfo torus_band(double R, double r, double Phi) {
    const double rho = R + r * std::cos(Phi), z = r * std::sin(Phi);
    std::vector<mio::BrepCurve> curves = {
        mio::detail::brep_curve_arc({0, 0, -z}, {1, 0, 0}, {0, 1, 0}, rho, 0.0, 2 * kPi),
        mio::detail::brep_curve_arc({0, 0, z}, {1, 0, 0}, {0, 1, 0}, rho, 0.0, 2 * kPi)};
    std::vector<brep_solids::SolidFace> faces(1);
    faces[0].mSurface.mKind = 5;
    faces[0].mSurface.mRadius = R;
    faces[0].mSurface.mMinor = r;
    faces[0].mLoops = {{brep_solids::circle_trim(rho, -z, 2 * kPi, 0, true, false)},
                       {brep_solids::circle_trim(rho, z, -2 * kPi, 1, false, false)}};
    return brep_solids::assemble(faces, curves);
}

}  // namespace

// 6.3.2/6.3.6: a capped cylinder sized by its curvature keeps every chord on
// its side within the deviation asked for, with fewer triangles than the
// uniform mesh that meets the same deviation (the flat caps stay coarse).
TEST(MeshCadSize, CurvatureSizingKeepsTheChordalErrorUnderTheBound) {
    const mio::BrepInfo cylinder = brep_solids::capped_cylinder(1.0, 2.0);
    for (const double d : {0.002, 0.0005}) {
        MeshCadSurfaceOptions o = sized(0.5);
        o.mCurvatureDeviation = d;
        o.mMinAngleDeg = 25.0;
        const auto r = mio::mesh_cad_surface(cylinder, o);
        ASSERT_EQ(r.mNumFailed, 0);
        EXPECT_EQ(r.mNumOpenEdges, 0);
        const double err = cylinder_chord_error(r.mMesh, 1.0);
        std::printf("d %g: %lld triangles, chord error %.6g\n", d,
                    static_cast<long long>(r.mNumTriangles), err);
        EXPECT_LE(err, d) << d;
        // The uniform size whose largest chord just meets d.
        MeshCadSurfaceOptions u = sized(0.5);
        u.mMinAngleDeg = 25.0;
        double hu = 0.5;
        mio::MeshCadSurfaceResult uniform;
        for (int k = 0; k < 60; ++k, hu *= 0.95) {
            u.mTargetSize = hu;
            uniform = mio::mesh_cad_surface(cylinder, u);
            if (cylinder_chord_error(uniform.mMesh, 1.0) <= d)
                break;
        }
        std::printf("   uniform h %.4f: %lld triangles\n", hu,
                    static_cast<long long>(uniform.mNumTriangles));
        EXPECT_LT(r.mNumTriangles, uniform.mNumTriangles) << d;
    }
}

// 6.3.7: a plate t thick gets at least n elements across t on its rims.
TEST(MeshCadSize, ThinWallsGetCellsAcross) {
    const double t = 0.05;
    const mio::BrepInfo plate = brep_solids::box(1.0, 1.0, t);
    for (const double n : {2.0, 3.0}) {
        MeshCadSurfaceOptions o = sized(0.25);
        o.mCellsAcrossGaps = n;
        o.mMinAngleDeg = 25.0;
        const auto r = mio::mesh_cad_surface(plate, o);
        ASSERT_EQ(r.mNumFailed, 0);
        EXPECT_EQ(r.mNumOpenEdges, 0);
        double rim = 0.0;
        for (const auto& r_t : triangles(r.mMesh)) {
            const Vec3 a = point(r.mMesh, r_t[0]), b = point(r.mMesh, r_t[1]),
                       c = point(r.mMesh, r_t[2]);
            const double zmin = std::min({a[2], b[2], c[2]}), zmax = std::max({a[2], b[2], c[2]});
            if (zmax - zmin > 1e-12)  // a rim triangle spans the thickness
                rim = std::max(rim, std::max({dist(a, b), dist(b, c), dist(c, a)}));
        }
        std::printf("n %g: %lld triangles, longest rim edge %.6g (t / n %.6g)\n", n,
                    static_cast<long long>(r.mNumTriangles), rim, t / n);
        EXPECT_LE(rim, t / n * (1.0 + 1e-9)) << n;
    }
}

// 6.3.3: the field itself grows no faster than grading - 1, holds each
// sample's size within its reach, and keeps between its floor and ceiling.
TEST(MeshCadSize, TheFieldIsGradationLimited) {
    std::vector<Vec3> xs;
    std::vector<double> hs;
    for (int i = 0; i < 40; ++i) {
        const double a = 0.37 * i;
        xs.push_back({std::cos(a) * (1.0 + 0.1 * i), std::sin(a), 0.05 * i});
        hs.push_back(0.01 + 0.004 * ((i * 7) % 11));
    }
    const double slope = 0.3, reach = 0.05;
    const mio::detail::MeshCadSizeField field(0.5, 0.005, slope, reach, xs, hs);
    for (std::size_t i = 0; i < xs.size(); ++i)
        EXPECT_LE(field.At(xs[i]), hs[i] + 1e-15);
    for (int k = 0; k < 400; ++k) {
        const Vec3 p{std::sin(1.1 * k) * 3.0, std::cos(0.7 * k) * 2.0, 0.01 * k};
        const Vec3 q{p[0] + 0.03 * std::cos(k), p[1] + 0.02, p[2] - 0.01 * std::sin(k)};
        const double hp = field.At(p), hq = field.At(q);
        EXPECT_LE(std::abs(hp - hq), slope * dist(p, q) + 1e-12) << k;
        EXPECT_GE(hp, 0.005);
        EXPECT_LE(hp, 0.5);
        // The brute-force envelope.
        double best = 0.5;
        for (std::size_t i = 0; i < xs.size(); ++i)
            best = std::min(best, hs[i] + slope * std::max(0.0, dist(p, xs[i]) - reach));
        EXPECT_DOUBLE_EQ(hp, std::max(0.005, best)) << k;
    }
}

// 6.3.3 on a mesh: a gentler grading spreads the fine size farther, so it
// takes more triangles; the size still reaches the target far away.
TEST(MeshCadSize, GradingSpreadsTheSize) {
    const mio::BrepInfo cylinder = brep_solids::capped_cylinder(1.0, 2.0);
    std::int64_t previous = 0;
    for (const double g : {2.0, 1.5, 1.2}) {
        MeshCadSurfaceOptions o = sized(0.5);
        o.mCurvatureDeviation = 0.002;
        o.mGrading = g;
        const auto r = mio::mesh_cad_surface(cylinder, o);
        ASSERT_EQ(r.mNumFailed, 0);
        std::printf("grading %g: %lld triangles\n", g, static_cast<long long>(r.mNumTriangles));
        EXPECT_GT(r.mNumTriangles, previous) << g;
        previous = r.mNumTriangles;
    }
}

// The sized mesh is the same over thread counts and face orders.
TEST(MeshCadSize, SizedMeshesAreDeterministic) {
    MeshCadSurfaceOptions o = sized(0.3);
    o.mCurvatureDeviation = 0.003;
    o.mCellsAcrossGaps = 2.0;
    o.mMinAngleDeg = 25.0;
    const mio::BrepInfo cylinder = brep_solids::capped_cylinder(0.5, 1.0);
    const std::uint64_t ref = digest(mio::mesh_cad_surface(cylinder, o).mMesh);
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (const int n : {1, 3}) {
        omp_set_num_threads(n);
        EXPECT_EQ(digest(mio::mesh_cad_surface(cylinder, o).mMesh), ref) << n;
    }
    omp_set_num_threads(saved);
#endif
    for (const std::uint64_t seed : {1ull, 7ull}) {
        mio::detail::mesh_cad_set_face_order_seed(seed);
        EXPECT_EQ(digest(mio::mesh_cad_surface(cylinder, o).mMesh), ref) << seed;
    }
    mio::detail::mesh_cad_set_face_order_seed(0);
}

// 6.4.2: rings around the sphere's poles make each a star of six, and the
// whole sphere meets min_angle = 25.
TEST(MeshCadSize, PoleRingsGradeThePoles) {
    for (const double h : {0.2, 0.1}) {
        MeshCadSurfaceOptions o = sized(h);
        o.mMinAngleDeg = 25.0;
        o.mPoleRings = 2;
        const auto r = mio::mesh_cad_surface(brep_solids::sphere(1.0), o);
        ASSERT_EQ(r.mNumFailed, 0);
        EXPECT_EQ(r.mNumOpenEdges, 0);
        EXPECT_GE(r.mMinAngleDeg, 25.0) << h;
        for (const double z : {1.0, -1.0}) {
            int valence = 0;
            for (const auto& r_t : triangles(r.mMesh))
                for (const std::int64_t v : r_t)
                    if (dist(point(r.mMesh, v), {0.0, 0.0, z}) < 1e-12)
                        ++valence;
            EXPECT_EQ(valence, 6) << h << " " << z;
        }
        o.mPoleRings = 0;
        const auto fan = mio::mesh_cad_surface(brep_solids::sphere(1.0), o);
        std::printf("h %g: smallest angle %.2f with rings, %.2f without\n", h, r.mMinAngleDeg,
                    fan.mMinAngleDeg);
        EXPECT_LT(fan.mMinAngleDeg, 25.0);
    }
}

// 6.4.3: no fixture folds; a narrow torus band, its tube wrapped nearly all
// the way round and h up to ten times its radius, does not either, and is
// meshed whole.
TEST(MeshCadSize, ANarrowTorusBandDoesNotFold) {
    for (const double phi : {1.0, 2.0, 3.0})
        for (const double h : {0.2, 0.5}) {
            const auto r = mio::mesh_cad_surface(torus_band(1.0, 0.05, phi), sized(h));
            EXPECT_EQ(r.mNumFailed, 0)
                << phi << " " << h << " " << (r.mFaces.empty() ? "" : r.mFaces[0].mReason);
            EXPECT_EQ(r.mNumMisorientedEdges, 0);
        }
}

namespace {

// A cylindrical band of radius R and height H swept by Sweep from the x
// axis, one open face whose four edges are its two arcs and two rulings.
mio::BrepInfo partial_band(double R, double H, double Sweep) {
    using namespace brep_solids;
    const double c = std::cos(Sweep), s = std::sin(Sweep);
    std::vector<mio::BrepCurve> curves = {
        mio::detail::brep_curve_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.0, Sweep),
        mio::detail::brep_curve_line({R * c, R * s, 0}, {R * c, R * s, H}),
        mio::detail::brep_curve_arc({0, 0, H}, {1, 0, 0}, {0, 1, 0}, R, 0.0, Sweep),
        mio::detail::brep_curve_line({R, 0, 0}, {R, 0, H})};
    std::vector<SolidFace> faces(1);
    faces[0].mSurface.mKind = 3;
    faces[0].mSurface.mRadius = R;
    faces[0].mLoops = {{circle_trim(R, 0.0, Sweep, 0, true, false),
                        line_trim({R * c, R * s, 0}, {R * c, R * s, H}, 1, true),
                        circle_trim(R, H, -Sweep, 2, false, false, Sweep),
                        line_trim({R, 0, H}, {R, 0, 0}, 3, false)}};
    return assemble(faces, curves);
}

}  // namespace

// Roadmap 6.4.3: meshed coarser than itself, a band swept 270 degrees or
// more gets one boundary piece per arc, and its two triangles cover the
// short way round, backwards: folded (on the real STEP models of the corpus
// the same happens on small fillets and blade faces). Before v17.18.0 they
// were written; the guard alone failed the face. Now refinement splits a
// fold, asking its arcs to be split across the faces that share them, until
// none is left: every triangle faces out at its corners, the band meshes.
TEST(MeshCadFolds, ABandSweptPastHalfATurnIsUnfolded) {
    for (const double sweep : {270.0, 300.0, 340.0})
        for (const double height : {0.05, 0.3})
            for (const double h : {8.0, 16.0})
                for (const double min_angle : {0.0, 25.0}) {
                    mio::MeshCadSurfaceOptions o;
                    o.mTargetSize = h;
                    o.mMinAngleDeg = min_angle;
                    const auto r =
                        mio::mesh_cad_surface(partial_band(1.0, height, sweep * kPi / 180.0), o);
                    ASSERT_EQ(r.mNumFailed, 0)
                        << "sweep " << sweep << ", h " << h << ": " << r.mFaces[0].mReason;
                    const auto& m = r.mMesh;
                    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
                        if (m.Cells(b).Type() != "triangle")
                            continue;
                        const auto& conn = m.Cells(b).Conn();
                        for (std::size_t c = 0; c < m.Cells(b).NumCells(); ++c) {
                            std::array<Vec3, 3> p{};
                            for (std::size_t k = 0; k < 3; ++k)
                                p[k] = point(m, mio::detail::read_int(conn, 3 * c + k));
                            const std::array<double, 3> e1{p[1][0] - p[0][0], p[1][1] - p[0][1],
                                                           p[1][2] - p[0][2]};
                            const std::array<double, 3> e2{p[2][0] - p[0][0], p[2][1] - p[0][1],
                                                           p[2][2] - p[0][2]};
                            const double nx = e1[1] * e2[2] - e1[2] * e2[1];
                            const double ny = e1[2] * e2[0] - e1[0] * e2[2];
                            // The outward normal at each corner is radial.
                            for (std::size_t k = 0; k < 3; ++k)
                                EXPECT_GT(nx * p[k][0] + ny * p[k][1], 0.0)
                                    << "sweep " << sweep << ", h " << h << ", triangle " << c;
                        }
                    }
                }
}

namespace {

// The median edge length of the triangles `cad:face` puts on face F.
double median_edge_on(const mio::Mesh& rM, std::int64_t F) {
    const mio::NDArray& face = rM.CellData(mio::kMeshCadFaceName, 0);
    const auto tris = triangles(rM);
    std::vector<double> e;
    for (std::size_t t = 0; t < tris.size(); ++t)
        if (mio::detail::read_int(face, t) == F)
            for (int k = 0; k < 3; ++k)
                e.push_back(dist(point(rM, tris[t][k]), point(rM, tris[t][(k + 1) % 3])));
    std::nth_element(e.begin(), e.begin() + e.size() / 2, e.end());
    return e[e.size() / 2];
}

// The faces (indices) each edge of a B-rep bounds, by its topology.
std::vector<std::vector<std::int64_t>> faces_of_edges(const mio::BrepInfo& rB) {
    std::vector<std::vector<std::int64_t>> out(rB.mEdges.size());
    for (std::size_t e = 0; e < rB.mEdges.size(); ++e)
        for (const mio::BrepRef& r_r : rB.mEdges[e].mTopology)
            for (std::size_t f = 0; f < rB.mFaces.size(); ++f)
                if (rB.mFaces[f].mId == r_r.mId)
                    out[e].push_back(static_cast<std::int64_t>(f));
    return out;
}

}  // namespace

// Roadmap 6.3.1 / 6.3.5 / 6.3.8: a face given its own size meshes at it, its
// neighbours grade away within the grading, and the report names exactly the
// edge a finer neighbour limits a face by.
TEST(MeshCadGivenSizes, AFaceMeshesAtItsGivenSize) {
    const mio::BrepInfo box = brep_solids::box();
    mio::MeshCadSizeSources src;
    src.mFaceSizes = {{1, 0.1}};
    const auto r = mio::mesh_cad_surface(box, sized(0.4), src);
    ASSERT_EQ(r.mNumFailed, 0);
    EXPECT_EQ(r.mNumOpenEdges, 0);
    EXPECT_NEAR(median_edge_on(r.mMesh, 1), 0.1, 0.015);
    // The field grows at most grading - 1 = 0.5 per unit length (the faces
    // across the box are 1 away), and never past target_size.
    for (std::int64_t f = 0; f < 6; ++f)
        if (f != 1)
            EXPECT_GT(median_edge_on(r.mMesh, f), 0.1) << f;

    const auto edges = faces_of_edges(box);
    std::int64_t shared = -1, neighbour = -1;
    for (std::size_t e = 0; e < edges.size() && shared < 0; ++e)
        if (edges[e].size() == 2 && (edges[e][0] == 1 || edges[e][1] == 1)) {
            shared = static_cast<std::int64_t>(e);
            neighbour = edges[e][0] == 1 ? edges[e][1] : edges[e][0];
        }
    ASSERT_GE(shared, 0);
    src.mFaceSizes = {{1, 0.1}, {neighbour, 0.4}};
    mio::MeshCadSizeReport report;
    const auto q = mio::mesh_cad_surface(box, sized(0.4), src, &report);
    ASSERT_EQ(q.mNumFailed, 0);
    for (std::size_t f = 0; f < report.mLimitedByEdges.size(); ++f)
        EXPECT_EQ(report.mLimitedByEdges[f], static_cast<std::int64_t>(f) == neighbour
                                                 ? std::vector<std::int64_t>{shared}
                                                 : std::vector<std::int64_t>{})
            << f;
}

TEST(MeshCadGivenSizes, AnEdgeAndABackgroundSizeTheMesh) {
    const mio::BrepInfo box = brep_solids::box();
    mio::MeshCadSizeSources src;
    src.mEdgeSizes = {{4, 0.05}};
    const auto r = mio::mesh_cad_surface(box, sized(0.4), src);
    ASSERT_EQ(r.mNumFailed, 0);
    const mio::NDArray& edge = r.mMesh.CellData(mio::kMeshCadEdgeName, 1);
    const mio::NDArray& conn = r.mMesh.Cells(1).Conn();
    double longest = 0.0;
    for (std::size_t c = 0; c < r.mMesh.Cells(1).NumCells(); ++c)
        if (mio::detail::read_int(edge, c) == 4)
            longest =
                std::max(longest, dist(point(r.mMesh, mio::detail::read_int(conn, 2 * c)),
                                       point(r.mMesh, mio::detail::read_int(conn, 2 * c + 1))));
    EXPECT_LE(longest, 0.05 * 1.000001);

    // A background plane under the box whose size grows along x.
    mio::Mesh bg;
    {
        const int n = 21;
        mio::NDArray pts(mio::DType::Float64, {std::size_t(n * n), 3});
        mio::NDArray size(mio::DType::Float64, {std::size_t(n * n)});
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                const double x = double(i) / (n - 1), y = double(j) / (n - 1);
                double* p = pts.As<double>() + 3 * (j * n + i);
                p[0] = x;
                p[1] = y;
                p[2] = 0.0;
                size.As<double>()[j * n + i] = 0.05 + 0.3 * x;
            }
        mio::NDArray tri(mio::DType::Int64, {std::size_t(2 * (n - 1) * (n - 1)), 3});
        std::int64_t* t = tri.As<std::int64_t>();
        for (int j = 0; j + 1 < n; ++j)
            for (int i = 0; i + 1 < n; ++i) {
                const std::int64_t a = j * n + i, b = a + 1, c = a + n, d = c + 1;
                const std::int64_t row[6] = {a, b, d, a, d, c};
                std::copy(row, row + 6, t);
                t += 6;
            }
        bg.AssignPoints(std::move(pts));
        bg.AddCellBlock("triangle", std::move(tri));
        bg.AddPointData("size", std::move(size));
    }
    mio::MeshCadSizeSources src_bg;
    src_bg.pBackground = &bg;
    src_bg.mBackgroundArray = "size";
    const auto b = mio::mesh_cad_surface(box, sized(0.4), src_bg);
    ASSERT_EQ(b.mNumFailed, 0);
    EXPECT_EQ(b.mNumOpenEdges, 0);
    // On the face at z = 0, the median edge per quarter of x follows the
    // background's size at the quarter's middle within 25%.
    const auto tris = triangles(b.mMesh);
    for (int q = 0; q < 4; ++q) {
        std::vector<double> e;
        for (const auto& r_t : tris) {
            const Vec3 p0 = point(b.mMesh, r_t[0]), p1 = point(b.mMesh, r_t[1]),
                       p2 = point(b.mMesh, r_t[2]);
            const double cx = (p0[0] + p1[0] + p2[0]) / 3.0;
            if (std::abs(p0[2]) + std::abs(p1[2]) + std::abs(p2[2]) > 1e-9 || cx < 0.25 * q ||
                cx >= 0.25 * (q + 1))
                continue;
            e.push_back(dist(p0, p1));
            e.push_back(dist(p1, p2));
            e.push_back(dist(p2, p0));
        }
        std::nth_element(e.begin(), e.begin() + e.size() / 2, e.end());
        const double target = 0.05 + 0.3 * (0.25 * q + 0.125);
        EXPECT_NEAR(e[e.size() / 2], target, 0.25 * target) << q;
    }
    EXPECT_EQ(digest(b.mMesh), digest(mio::mesh_cad_surface(box, sized(0.4), src_bg).mMesh));
}

TEST(MeshCadGivenSizes, EmptySourcesAreTheTwoArgumentOverloadAndBadOnesAreRefused) {
    const mio::BrepInfo box = brep_solids::box();
    EXPECT_EQ(digest(mio::mesh_cad_surface(box, sized(0.3), mio::MeshCadSizeSources{}).mMesh),
              digest(mio::mesh_cad_surface(box, sized(0.3)).mMesh));
    const auto refused = [&](mio::MeshCadSizeSources rS) {
        EXPECT_THROW(mio::mesh_cad_surface(box, sized(0.3), rS), std::invalid_argument);
    };
    mio::MeshCadSizeSources s;
    s.mFaceSizes = {{6, 0.1}};
    refused(s);
    s.mFaceSizes = {{0, 0.0}};
    refused(s);
    s.mFaceSizes = {{0, std::nan("")}};
    refused(s);
    s = {};
    s.mEdgeSizes = {{-1, 0.1}};
    refused(s);
    s = {};
    mio::Mesh no_array = mio::mesh_cad_surface(box, sized(0.5)).mMesh;
    s.pBackground = &no_array;
    s.mBackgroundArray = "size";
    refused(s);
}
