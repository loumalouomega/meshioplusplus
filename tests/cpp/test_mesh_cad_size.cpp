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
