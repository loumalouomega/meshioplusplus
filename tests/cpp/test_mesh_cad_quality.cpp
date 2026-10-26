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
// mesh_cad_surface's quality options (roadmap 6.2): Delaunay refinement to
// a minimum angle measured on the surface, boundary splits shared across
// faces, small input angles and poles reported with the triangles they
// force, guarded smoothing, the order faces are processed in, and the
// options not implemented yet refused. Each probe of roadmap 6.2.7 is a
// named test; 0 keeping the unrefined mesh is pinned in test_op_goldens.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

// Project includes
#include "brep_solids.hpp"
#include "mesh_fixtures.hpp"
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "../../src/cpp/src/detail/cdt2d.hpp"
#include "../../src/cpp/src/operations/mesh_cad_testing.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::MeshCadSurfaceOptions;
using mio::MeshCadSurfaceResult;
using Vec2 = std::array<double, 2>;

constexpr double kPi = 3.14159265358979323846;

MeshCadSurfaceOptions quality(double H, double MinAngle, std::int64_t Smooth = 0) {
    MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    o.mMinAngleDeg = MinAngle;
    o.mSmoothIterations = Smooth;
    return o;
}

std::uint64_t digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

struct Tri {
    std::array<std::int64_t, 3> mIds;
    std::array<const double*, 3> mX;
    double mMinAngle;
};

double angle_at(const double* pO, const double* pP, const double* pQ) {
    const double u[3] = {pP[0] - pO[0], pP[1] - pO[1], pP[2] - pO[2]};
    const double v[3] = {pQ[0] - pO[0], pQ[1] - pO[1], pQ[2] - pO[2]};
    const double c[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                         u[0] * v[1] - u[1] * v[0]};
    return std::atan2(std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]),
                      u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) *
           180.0 / kPi;
}

std::vector<Tri> triangles(const mio::Mesh& rM) {
    std::vector<Tri> out;
    const double* x = rM.Points().As<double>();
    const auto cells = rM.Cells(0);
    const mio::NDArray& c = cells.Conn();
    for (std::size_t t = 0; t < cells.NumCells(); ++t) {
        Tri r;
        for (std::size_t k = 0; k < 3; ++k) {
            r.mIds[k] = mio::detail::read_int(c, 3 * t + k);
            r.mX[k] = x + 3 * r.mIds[k];
        }
        r.mMinAngle =
            std::min({angle_at(r.mX[0], r.mX[1], r.mX[2]), angle_at(r.mX[1], r.mX[2], r.mX[0]),
                      angle_at(r.mX[2], r.mX[0], r.mX[1])});
        out.push_back(r);
    }
    return out;
}

double distance(const double* pA, const double* pB) {
    return std::sqrt((pA[0] - pB[0]) * (pA[0] - pB[0]) + (pA[1] - pB[1]) * (pA[1] - pB[1]) +
                     (pA[2] - pB[2]) * (pA[2] - pB[2]));
}

struct Measures {
    double mArea = 0.0;
    double mVolume = 0.0;
};

Measures measure(const mio::Mesh& rM) {
    Measures m;
    for (const Tri& r_t : triangles(rM)) {
        const double* a = r_t.mX[0];
        const double* b = r_t.mX[1];
        const double* d = r_t.mX[2];
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double v[3] = {d[0] - a[0], d[1] - a[1], d[2] - a[2]};
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                             u[0] * v[1] - u[1] * v[0]};
        m.mArea += 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        m.mVolume += (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
                      a[2] * (b[0] * d[1] - b[1] * d[0])) /
                     6.0;
    }
    return m;
}

void expect_closed(const MeshCadSurfaceResult& rR) {
    EXPECT_EQ(rR.mNumFailed, 0);
    EXPECT_EQ(rR.mNumOpenEdges, 0);
    EXPECT_EQ(rR.mNumNonManifoldEdges, 0);
    EXPECT_EQ(rR.mNumMisorientedEdges, 0);
}

// A planar face on the patch [0, S]^2 at z = 0 (its parameters are its
// coordinates), its loops polygons of straight trims.
mio::BrepInfo planar_face(const std::vector<std::vector<Vec2>>& rLoops, double S) {
    std::ostringstream j;
    j.precision(17);
    j << R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "swapped_surface_normal": false,
        "surface": {"is_trimmed": true, "is_rational": false, "degrees": [1, 1],
        "knot_vectors": [[0, )"
      << S << "], [0, " << S << R"(]], "control_points": [[1, [0, 0, 0, 1]], [2, [)" << S
      << R"(, 0, 0, 1]], [3, [0, )" << S << R"(, 0, 1]], [4, [)" << S << ", " << S
      << R"(, 0, 1]]]}, "boundary_loops": [)";
    int index = 0;
    for (std::size_t l = 0; l < rLoops.size(); ++l) {
        j << (l ? ", " : "") << R"({"loop_type": ")" << (l ? "inner" : "outer")
          << R"(", "trimming_curves": [)";
        for (std::size_t k = 0; k < rLoops[l].size(); ++k) {
            const Vec2& a = rLoops[l][k];
            const Vec2& b = rLoops[l][(k + 1) % rLoops[l].size()];
            j << (k ? ", " : "") << R"({"trim_index": )" << index++
              << R"(, "curve_direction": true, "parameter_curve": {"degree": 1,
                "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[)"
              << a[0] << ", " << a[1] << R"(, 0, 1]], [[)" << b[0] << ", " << b[1]
              << R"(, 0, 1]]]}})";
        }
        j << "]}";
    }
    j << R"(], "embedded_edges": []}], "edges": [], "vertices": []}]})";
    const std::string path = mt::temp_path(".cad.json");
    std::ofstream(path) << j.str();
    return mio::read_cad_json_info(path);
}

// An untrimmed planar bilinear patch on [0, 1]^2 with corners (0, 0),
// (1, 0), (0, 2) and `rCorner`: (1, 2) is a rectangle, its parameters
// affine; (1, 1) gives x = u, y = v (2 - u), which the face's one affine
// map stretches 2.4-fold at that corner.
mio::BrepInfo bilinear_patch(const Vec2& rCorner) {
    std::ostringstream j;
    j.precision(17);
    j << R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "swapped_surface_normal": false,
        "surface": {"is_trimmed": false, "is_rational": false, "degrees": [1, 1],
        "knot_vectors": [[0, 1], [0, 1]], "control_points": [[1, [0, 0, 0, 1]],
        [2, [1, 0, 0, 1]], [3, [0, 2, 0, 1]], [4, [)"
      << rCorner[0] << ", " << rCorner[1] << R"(, 0, 1]]]}, "boundary_loops": [],
        "embedded_edges": []}], "edges": [], "vertices": []}]})";
    const std::string path = mt::temp_path(".cad.json");
    std::ofstream(path) << j.str();
    return mio::read_cad_json_info(path);
}

// The trimmed plate of roadmap 6.2: 4 x 4 with a diamond hole (the
// tests/cpp/test_iga.cpp plate, its trims as straight pieces).
mio::BrepInfo trimmed_plate() {
    return planar_face({{{0, 0}, {4, 0}, {4, 4}, {0, 4}}, {{2, 0.8}, {3.2, 2}, {2, 3.2}, {0.8, 2}}},
                       4.0);
}

// Every triangle under the bound touches a corner the report lists: it has
// a vertex at a small angle's point, or within a pole's reach (at least h).
void expect_bound_or_reported(const MeshCadSurfaceResult& rR, double MinAngle, double H) {
    std::vector<std::pair<std::array<double, 3>, double>> poles;
    std::vector<std::int64_t> corners;
    for (const auto& r_f : rR.mFaces)
        for (const auto& r_c : r_f.mSmallAngleCorners)
            if (r_c.mPole) {
                EXPECT_GE(r_c.mReach, H);
                poles.emplace_back(r_c.mX, r_c.mReach);
            } else {
                corners.push_back(r_c.mPoint);
            }
    std::size_t below = 0;
    for (const Tri& r_t : triangles(rR.mMesh)) {
        if (r_t.mMinAngle >= MinAngle - 1e-9)
            continue;
        ++below;
        bool reported = false;
        for (int k = 0; k < 3 && !reported; ++k) {
            reported = std::find(corners.begin(), corners.end(), r_t.mIds[k]) != corners.end();
            for (const auto& [r_x, reach] : poles)
                reported = reported || distance(r_t.mX[k], r_x.data()) <= reach * (1.0 + 1e-12);
        }
        EXPECT_TRUE(reported) << "a triangle of smallest angle " << r_t.mMinAngle
                              << " touches no reported corner";
    }
    std::int64_t faces_below = 0;
    for (const auto& r_f : rR.mFaces)
        faces_below += r_f.mBelowBound ? 1 : 0;
    EXPECT_EQ(faces_below, rR.mNumFacesBelowBound);
    if (below > 0)
        EXPECT_GT(rR.mNumFacesBelowBound, 0);
}

TEST(MeshCadQuality, PlateCylinderAndSphereMeetTwentyFiveDegrees) {
    // Roadmap 6.2.7: every triangle meets min_angle = 25 except those
    // touching a reported corner. The plate and the cylinder have none; the
    // sphere's two poles are reported, the triangles within h of them left
    // to graded poles (6.4.2).
    for (const double h : {0.5, 0.2, 0.1}) {
        const auto plate = mio::mesh_cad_surface(trimmed_plate(), quality(h, 25.0));
        ASSERT_EQ(plate.mNumFailed, 0) << plate.mFaces[0].mReason;
        EXPECT_NEAR(measure(plate.mMesh).mArea, 16.0 - 2.88, 1e-12) << h;
        EXPECT_GE(plate.mMinAngleDeg, 25.0 - 1e-9) << h;
        EXPECT_EQ(plate.mNumFacesBelowBound, 0) << h;
        EXPECT_TRUE(plate.mFaces[0].mSmallAngleCorners.empty());

        const auto cyl =
            mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), quality(h, 25.0));
        expect_closed(cyl);
        EXPECT_GE(cyl.mMinAngleDeg, 25.0 - 1e-9) << h;
        EXPECT_EQ(cyl.mNumFacesBelowBound, 0) << h;

        const auto sph = mio::mesh_cad_surface(brep_solids::sphere(1.0), quality(h, 25.0));
        expect_closed(sph);
        expect_bound_or_reported(sph, 25.0, h);
        std::size_t poles = 0;
        for (const auto& r_f : sph.mFaces) {
            for (const auto& r_c : r_f.mSmallAngleCorners) {
                EXPECT_TRUE(r_c.mPole);
                EXPECT_NEAR(std::fabs(r_c.mX[2]), 1.0, 1e-12);
                ++poles;
            }
            // Away from its pole every face reaches the bound.
            EXPECT_GE(r_f.mBoundDeg, 25.0 - 1e-9) << h;
        }
        EXPECT_EQ(poles, 2u) << h;
        // Every point still on the exact sphere.
        const double* x = sph.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < sph.mMesh.NumPoints(); ++i)
            EXPECT_NEAR(std::sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] +
                                  x[3 * i + 2] * x[3 * i + 2]),
                        1.0, 1e-12);
    }
}

TEST(MeshCadQuality, AreaAndVolumeStillConvergeAtSecondOrder) {
    const double R = 1.0, H = 2.0;
    const double area = 2.0 * kPi * R * H + 2.0 * kPi * R * R;
    const double volume = kPi * R * R * H;
    double prev_area = 0.0, prev_volume = 0.0;
    for (const double h : {0.2, 0.1, 0.05}) {
        const auto r = mio::mesh_cad_surface(brep_solids::capped_cylinder(R, H), quality(h, 25.0));
        expect_closed(r);
        const Measures m = measure(r.mMesh);
        const double ea = std::fabs(m.mArea - area);
        const double ev = std::fabs(m.mVolume - volume);
        if (prev_area > 0.0) {
            EXPECT_LT(ea, prev_area / 2.5) << h;  // about 4 for second order
            EXPECT_LT(ev, prev_volume / 2.5) << h;
        }
        prev_area = ea;
        prev_volume = ev;
    }
    const double sphere_volume = 4.0 / 3.0 * kPi;
    double prev = 0.0;
    for (const double h : {0.2, 0.1, 0.05}) {
        const auto r = mio::mesh_cad_surface(brep_solids::sphere(1.0), quality(h, 25.0));
        const double ev = std::fabs(measure(r.mMesh).mVolume - sphere_volume);
        if (prev > 0.0)
            EXPECT_LT(ev, prev / 2.5) << h;
        prev = ev;
    }
}

TEST(MeshCadQuality, BoundarySplitsAreSharedByEveryFaceOfTheEdge) {
    // Roadmap 6.2.2. At 30 degrees the cylinder's faces ask for rim splits;
    // each is halved by arc length in the band and in the cap, so the solid
    // stays closed and every rim point is on the exact circle. Splitting by
    // arc length keeps a split piece within Ruppert's straight-segment
    // argument: its sagitta stays a small fraction of its length.
    const double h = 0.2;
    const auto r = mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), quality(h, 30.0));
    expect_closed(r);
    std::int64_t interior = 0;
    for (const auto& r_f : r.mFaces)
        interior += r_f.mNumRefinementPoints;
    EXPECT_GT(r.mNumRefinementPoints, interior);  // edge points too
    const auto plain = mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), quality(h, 0));
    EXPECT_GT(r.mNumEdgeSegments, plain.mNumEdgeSegments);
    // The line cells are the edges' pieces: on a rim of radius 1, a chord c
    // has sagitta 1 - sqrt(1 - c^2 / 4).
    const double* x = r.mMesh.Points().As<double>();
    const auto lines = r.mMesh.Cells(1);
    double worst = 0.0;
    for (std::size_t s = 0; s < lines.NumCells(); ++s) {
        const double* a = x + 3 * mio::detail::read_int(lines.Conn(), 2 * s);
        const double* b = x + 3 * mio::detail::read_int(lines.Conn(), 2 * s + 1);
        if (std::fabs(a[2] - b[2]) > 1e-12)
            continue;  // the seam
        EXPECT_NEAR(std::hypot(a[0], a[1]), 1.0, 1e-12);
        const double c = distance(a, b);
        worst = std::max(worst, (1.0 - std::sqrt(1.0 - 0.25 * c * c)) / c);
    }
    EXPECT_LT(worst, 0.03);
}

TEST(MeshCadQuality, ASmallInputAngleIsReportedWithTheTrianglesItForces) {
    // Roadmap 6.2.3: a planar face whose corner at the origin is 15 degrees.
    // The corner is reported with its angle; every triangle under the bound
    // holds it or spans it (its shortest edge between the corner's two
    // trims), and refinement terminates.
    const double t = std::tan(15.0 * kPi / 180.0);
    const auto info = planar_face({{{0, 0}, {4, 0}, {4, 4 * t}}}, 4.0);
    for (const double h : {0.5, 0.2}) {
        const auto r = mio::mesh_cad_surface(info, quality(h, 25.0));
        ASSERT_EQ(r.mNumFailed, 0) << r.mFaces[0].mReason;
        EXPECT_NEAR(measure(r.mMesh).mArea, 0.5 * 4.0 * 4.0 * t, 1e-12);
        const auto& f = r.mFaces[0];
        ASSERT_EQ(f.mSmallAngleCorners.size(), 1u);
        const auto& c = f.mSmallAngleCorners[0];
        EXPECT_FALSE(c.mPole);
        EXPECT_NEAR(c.mAngleDeg, 15.0, 1e-9);
        EXPECT_NEAR(f.mMinInputAngleDeg, 15.0, 1e-9);
        EXPECT_NEAR(c.mX[0], 0.0, 1e-15);
        EXPECT_NEAR(c.mX[1], 0.0, 1e-15);
        EXPECT_NEAR(c.mU, 0.0, 1e-15);
        EXPECT_TRUE(f.mBelowBound);
        EXPECT_GE(f.mBoundDeg, 25.0 - 1e-9);
        const auto on = [&](const double* pP, double Slope) {
            return pP[0] > 0.0 && std::fabs(pP[1] - Slope * pP[0]) <= 1e-12;
        };
        double forced_min = 180.0;
        for (const Tri& r_t : triangles(r.mMesh)) {
            if (r_t.mMinAngle >= 25.0 - 1e-9)
                continue;
            forced_min = std::min(forced_min, r_t.mMinAngle);
            bool holds = false;
            int sides = 0;
            for (int k = 0; k < 3; ++k) {
                holds = holds || r_t.mIds[k] == c.mPoint;
                sides |= (on(r_t.mX[k], 0.0) ? 1 : 0) | (on(r_t.mX[k], t) ? 2 : 0);
            }
            EXPECT_TRUE(holds || sides == 3) << r_t.mMinAngle;
        }
        EXPECT_NEAR(c.mMinForcedDeg, forced_min, 1e-9);
    }
}

// A planar lens on the patch [0, S]^2 at z = 0: two circular trims (rational
// quadratic arcs) from (X0, Y0) to (X0 + L, Y0), each leaving the chord at
// Theta degrees, so both corners are 2 Theta.
mio::BrepInfo planar_lens(double X0, double Y0, double L, double ThetaDeg, double S) {
    const double th = ThetaDeg * kPi / 180.0;
    const double w = std::cos(th);
    const double dy = 0.5 * L * std::tan(th);
    std::ostringstream j;
    j.precision(17);
    j << R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "swapped_surface_normal": false,
        "surface": {"is_trimmed": true, "is_rational": false, "degrees": [1, 1],
        "knot_vectors": [[0, )"
      << S << "], [0, " << S << R"(]], "control_points": [[1, [0, 0, 0, 1]], [2, [)" << S
      << R"(, 0, 0, 1]], [3, [0, )" << S << R"(, 0, 1]], [4, [)" << S << ", " << S
      << R"(, 0, 1]]]}, "boundary_loops": [{"loop_type": "outer", "trimming_curves": [)";
    const auto arc = [&](int Index, double Ax, double Bx, double My) {
        j << (Index ? ", " : "") << R"({"trim_index": )" << Index
          << R"(, "curve_direction": true, "parameter_curve": {"is_rational": true, "degree": 2,
            "knot_vector": [0, 0, 0, 1, 1, 1], "active_range": [0, 1], "control_points": [[[)"
          << Ax << ", " << Y0 << R"(, 0, 1]], [[)" << X0 + 0.5 * L << ", " << My << ", 0, " << w
          << R"(]], [[)" << Bx << ", " << Y0 << R"(, 0, 1]]]}})";
    };
    arc(0, X0, X0 + L, Y0 - dy);
    arc(1, X0 + L, X0, Y0 + dy);
    j << R"(]}], "embedded_edges": []}], "edges": [], "vertices": []}]})";
    const std::string path = mt::temp_path(".cad.json");
    std::ofstream(path) << j.str();
    return mio::read_cad_json_info(path);
}

TEST(MeshCadQuality, CornerLoppingStatesItsBoundOnALens) {
    // Roadmap 6.2.3, Pav and Walkington (2005): a planar face with two
    // 5-degree corners between circular trims, refined at 20 degrees with
    // its corners lopped. Both are reported with their angle, a ball inside
    // (sqrt 2 - 1) times the corners' distance, and Corollary 3's bound for
    // the trims' turning, which 20 is under; the run ends with no budget
    // warning, and every triangle under 20 degrees lies within twice its
    // corner's radius (Lemma 7).
    const double l = 4.0;
    const auto info = planar_lens(0.5, 1.0, l, 2.5, 5.0);
    for (const double h : {0.05, 0.02}) {
        mio::MeshCadSurfaceExtras extras;
        extras.mLopCorners = true;
        mio::MeshCadSurfaceExtrasReport x;
        const auto r = mio::mesh_cad_surface(info, quality(h, 20.0), mio::MeshCadSizeSources{},
                                             extras, nullptr, &x);
        ASSERT_EQ(r.mNumFailed, 0) << r.mFaces[0].mReason;
        for (const auto& r_w : r.mWarnings)
            EXPECT_EQ(r_w.find("did not finish"), std::string::npos) << r_w;
        ASSERT_EQ(x.mLoppedCorners.size(), 2u) << h;
        const double* pts = r.mMesh.Points().As<double>();
        for (const auto& c : x.mLoppedCorners) {
            EXPECT_EQ(c.mFace, 0);
            // Between the first chords, under the tangents' 5 degrees.
            EXPECT_LE(c.mAngleDeg, 5.0) << h;
            EXPECT_GT(c.mAngleDeg, 4.9) << h;
            EXPECT_NEAR(c.mLfs0, l, 1e-9);
            EXPECT_GT(c.mLopRadius, 0.0);
            EXPECT_LE(c.mLopRadius, (std::sqrt(2.0) - 1.0) * l * (1 + 1e-12));
            EXPECT_GT(c.mStatedBoundDeg, 20.0) << h;
            EXPECT_LE(c.mStatedBoundDeg, mio::detail::cdt2d_lop_bound_deg(0.0));
            EXPECT_LT(c.mTurningDeg, 60.0);
            EXPECT_TRUE(c.mWithinTwoRadii) << h;
            EXPECT_TRUE(c.mGuaranteed) << h << ": " << c.mReason;
            EXPECT_GE(c.mNumForced, 1);
            EXPECT_LT(c.mMinForcedDeg, 20.0);
        }
        for (const Tri& r_t : triangles(r.mMesh)) {
            if (r_t.mMinAngle >= 20.0 - 1e-9)
                continue;
            bool near = false;
            for (const auto& c : x.mLoppedCorners) {
                const double* o = pts + 3 * c.mPoint;
                bool all = true;
                for (int k = 0; k < 3; ++k)
                    all = all && distance(r_t.mX[k], o) <= 2.0 * c.mLopRadius * (1 + 1e-9);
                near = near || all;
            }
            EXPECT_TRUE(near) << "h " << h << ": a triangle of " << r_t.mMinAngle;
        }
    }
}

TEST(MeshCadQuality, CornerLoppingSaysWhyItStatesNoBound) {
    // Above the bound, or on a curved face, the corner is still lopped and
    // reported, with the reason no bound is stated.
    const double t = std::tan(15.0 * kPi / 180.0);
    const auto plane = planar_face({{{0, 0}, {4, 0}, {4, 4 * t}}}, 4.0);
    mio::MeshCadSurfaceExtras extras;
    extras.mLopCorners = true;
    mio::MeshCadSurfaceExtrasReport x;
    const auto r = mio::mesh_cad_surface(plane, quality(0.3, 25.0), mio::MeshCadSizeSources{},
                                         extras, nullptr, &x);
    ASSERT_EQ(r.mNumFailed, 0);
    ASSERT_FALSE(x.mLoppedCorners.empty());
    for (const auto& c : x.mLoppedCorners) {
        EXPECT_NEAR(c.mStatedBoundDeg, mio::detail::cdt2d_lop_bound_deg(0.0), 1e-9);
        EXPECT_FALSE(c.mGuaranteed);
        EXPECT_NE(c.mReason.find("not under the bound"), std::string::npos) << c.mReason;
    }
    const auto r20 = mio::mesh_cad_surface(plane, quality(0.3, 20.0), mio::MeshCadSizeSources{},
                                           extras, nullptr, &x);
    ASSERT_EQ(r20.mNumFailed, 0);
    for (const auto& c : x.mLoppedCorners)
        EXPECT_TRUE(c.mGuaranteed) << c.mReason;
    // Off, the mesh is the shells' byte for byte.
    const auto off = mio::mesh_cad_surface(plane, quality(0.3, 20.0), mio::MeshCadSizeSources{},
                                           mio::MeshCadSurfaceExtras{}, nullptr, &x);
    EXPECT_TRUE(x.mLoppedCorners.empty());
    EXPECT_EQ(digest(off.mMesh), digest(mio::mesh_cad_surface(plane, quality(0.3, 20.0)).mMesh));
}

TEST(MeshCadQuality, TheSpherePinsItsPointCount) {
    // Roadmap 6.2.4: quality is measured on the lifted triangles and the
    // plane only finds the point, so termination is no longer the plane's
    // theorem; the number of points the sphere needs is pinned.
    const auto plain = mio::mesh_cad_surface(brep_solids::sphere(1.0), quality(0.2, 0.0));
    const auto r = mio::mesh_cad_surface(brep_solids::sphere(1.0), quality(0.2, 25.0));
    EXPECT_EQ(plain.mNumPoints, 673);
    EXPECT_EQ(r.mNumPoints, 677);
    // Refinement's points (inside the faces and on their edges); a boundary
    // split also moves the lattice away from it, so the totals differ.
    EXPECT_EQ(r.mNumRefinementPoints, 8);
}

TEST(MeshCadQuality, TheMeshDoesNotDependOnTheFaceOrderOrThreads) {
    // Roadmap 6.2.8: bad triangles are taken worst first, ties broken by
    // vertex ids; faces are meshed independently and their boundary splits
    // merged in (edge, piece) order. Visiting the faces in a shuffled order,
    // or on any number of threads, gives the same mesh.
    const mio::BrepInfo cyl = brep_solids::capped_cylinder(1.0, 2.0);
    const mio::BrepInfo box = brep_solids::box(2, 1, 3);
    const std::uint64_t ref_cyl = digest(mio::mesh_cad_surface(cyl, quality(0.2, 30.0)).mMesh);
    const std::uint64_t ref_box = digest(mio::mesh_cad_surface(box, quality(0.3, 30.0, 2)).mMesh);
    for (const std::uint64_t seed : {1u, 7u, 12345u}) {
        mio::detail::mesh_cad_set_face_order_seed(seed);
        EXPECT_EQ(digest(mio::mesh_cad_surface(cyl, quality(0.2, 30.0)).mMesh), ref_cyl) << seed;
        EXPECT_EQ(digest(mio::mesh_cad_surface(box, quality(0.3, 30.0, 2)).mMesh), ref_box) << seed;
    }
    mio::detail::mesh_cad_set_face_order_seed(0);
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
        EXPECT_EQ(digest(mio::mesh_cad_surface(cyl, quality(0.2, 30.0)).mMesh), ref_cyl) << threads;
    }
    omp_set_num_threads(saved);
#endif
}

TEST(MeshCadQuality, SmoothingKeepsTheBoundaryAndNeverLowersAStar) {
    // Roadmap 6.2.5: interior points move (Laplacian or ODT), boundary
    // points stay where they were, every point stays on the sphere, and the
    // smallest angle away from the poles does not drop.
    const auto plain = mio::mesh_cad_surface(brep_solids::sphere(1.0), quality(0.2, 0.0));
    const double* x0 = plain.mMesh.Points().As<double>();
    for (const auto method : {mio::MeshCadSmoothing::Laplacian, mio::MeshCadSmoothing::Odt}) {
        MeshCadSurfaceOptions o = quality(0.2, 0.0, 3);
        o.mSmoothMethod = method;
        const auto r = mio::mesh_cad_surface(brep_solids::sphere(1.0), o);
        expect_closed(r);
        EXPECT_GT(r.mNumSmoothed, 0);
        ASSERT_EQ(r.mNumPoints, plain.mNumPoints);
        const double* x = r.mMesh.Points().As<double>();
        // The boundary (seam, equator, poles) comes first and is shared.
        const std::size_t boundary =
            static_cast<std::size_t>(plain.mNumPoints - plain.mFaces[0].mNumInteriorPoints -
                                     plain.mFaces[1].mNumInteriorPoints);
        for (std::size_t i = 0; i < 3 * boundary; ++i)
            ASSERT_EQ(x[i], x0[i]);
        for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
            EXPECT_NEAR(std::sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] +
                                  x[3 * i + 2] * x[3 * i + 2]),
                        1.0, 1e-12);
        for (std::size_t f = 0; f < 2; ++f)
            EXPECT_GE(r.mFaces[f].mBoundDeg, plain.mFaces[f].mBoundDeg) << f;
        EXPECT_NEAR(measure(r.mMesh).mVolume, measure(plain.mMesh).mVolume, 0.01);
    }
    // On a plane, refinement then smoothing still meets the bound.
    const auto r = mio::mesh_cad_surface(trimmed_plate(), quality(0.2, 25.0, 3));
    EXPECT_GE(r.mMinAngleDeg, 25.0 - 1e-9);
    EXPECT_NEAR(measure(r.mMesh).mArea, 16.0 - 2.88, 1e-12);
}

TEST(MeshCadQuality, TheHistogramCountsEveryAngle) {
    const auto r =
        mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), quality(0.2, 25.0));
    std::int64_t total = 0;
    for (std::size_t b = 0; b < mio::kMeshCadAngleBins; ++b) {
        total += r.mAngleHistogram[b];
        if (b < 2)
            EXPECT_EQ(r.mAngleHistogram[b], 0) << b;  // nothing under 20 degrees
    }
    EXPECT_EQ(total, 3 * r.mNumTriangles);
    EXPECT_EQ(r.mNumCurvedStraightened, 0);
}

TEST(MeshCadQuality, OptionsNotImplementedYetAreRefused) {
    const auto refused = [](const MeshCadSurfaceOptions& rO, const std::string& rItem) {
        try {
            (void)mio::mesh_cad_surface(brep_solids::box(), rO);
            ADD_FAILURE() << rItem << " accepted";
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find("roadmap " + rItem), std::string::npos)
                << e.what();
        }
    };
    MeshCadSurfaceOptions o = quality(0.5, 0.0);
    o.mInsertion = mio::MeshCadInsertion::Offcentre;
    refused(o, "6.2.9");
    // A varying metric (6.4.1) has landed; a value past the enum is refused.
    o = quality(0.5, 0.0);
    o.mMetric = mio::MeshCadMetric::Varying;
    EXPECT_EQ(mio::mesh_cad_surface(brep_solids::box(), o).mNumFailed, 0);
    o.mMetric = static_cast<mio::MeshCadMetric>(7);
    EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), o), std::invalid_argument);
    // The discrete reparametrisation (6.4.4) has landed: 0 or above 1, order
    // 1 only.
    o = quality(0.5, 0.0);
    o.mReparamConditionNumber = 100.0;
    EXPECT_EQ(mio::mesh_cad_surface(brep_solids::box(), o).mNumFailed, 0);
    for (const double bad : {-1.0, 0.5, 1.0, std::nan("")}) {
        o.mReparamConditionNumber = bad;
        EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), o), std::invalid_argument)
            << bad;
    }
    o.mReparamConditionNumber = 100.0;
    o.mOrder = 2;
    EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), o), std::invalid_argument);
    for (const int bad : {0, 4}) {
        o = quality(0.5, 0.0);
        o.mOrder = bad;
        EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), o), std::invalid_argument)
            << bad;
    }
    for (const double bad : {-1.0, 60.0, std::nan("")})
        EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), quality(0.5, bad)),
                     std::invalid_argument)
            << bad;
    EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), quality(0.5, 0.0, -1)),
                 std::invalid_argument);
    // The sizing options (6.3) and pole rings (6.4.2) have landed; their
    // ranges are checked.
    for (const double bad : {0.5, 1.0, -1.0, std::nan("")}) {
        o = quality(0.5, 0.0);
        o.mGrading = bad;
        EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), o), std::invalid_argument)
            << bad;
    }
    for (const double bad : {-1.0, std::nan("")}) {
        o = quality(0.5, 0.0);
        o.mCurvatureDeviation = bad;
        EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), o), std::invalid_argument);
        o = quality(0.5, 0.0);
        o.mCellsAcrossGaps = bad;
        EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), o), std::invalid_argument);
    }
    for (const int bad : {-1, 9}) {
        o = quality(0.5, 0.0);
        o.mPoleRings = bad;
        EXPECT_THROW((void)mio::mesh_cad_surface(brep_solids::box(), o), std::invalid_argument);
    }
}

TEST(MeshCadQuality, AStretchedCornerDoesNotRunAway) {
    // Found by the fuzz target (tests/fuzz/regressions/op_mesh_cad_surface):
    // where the face's metric departs more than twofold from its plane's,
    // a point placed for the surface can leave triangles as bad as the one
    // it replaced, each smaller, and refinement ran away (198853 points at
    // max_points 200000, h = 0.0588...). There no point lands closer than
    // h / 4 and no triangle under h / 4 is refined; one left under the
    // bound is warned and kept out of bound_deg, as a pole's are. The
    // rectangle, whose parameters are affine, is the other side.
    const auto warned = [](const MeshCadSurfaceResult& rR) {
        for (const std::string& r_w : rR.mWarnings)
            if (r_w.find("roadmap 6.4.1") != std::string::npos)
                return true;
        return false;
    };
    for (const double h : {0.2, 0.1, 0.05884389414473131, 0.05}) {
        for (const double ang : {20.7, 25.0}) {
            const auto flat = mio::mesh_cad_surface(bilinear_patch({1.0, 2.0}), quality(h, ang));
            ASSERT_EQ(flat.mNumFailed, 0);
            EXPECT_GE(flat.mMinAngleDeg, ang - 1e-9) << h << " " << ang;
            EXPECT_FALSE(warned(flat)) << h << " " << ang;

            const auto unrefined =
                mio::mesh_cad_surface(bilinear_patch({1.0, 1.0}), quality(h, 0.0));
            const auto r = mio::mesh_cad_surface(bilinear_patch({1.0, 1.0}), quality(h, ang));
            ASSERT_EQ(r.mNumFailed, 0);
            EXPECT_LT(r.mNumPoints, unrefined.mNumPoints + unrefined.mNumPoints / 10)
                << h << " " << ang;
            for (const std::string& r_w : r.mWarnings)
                EXPECT_EQ(r_w.find("refinement stopped"), std::string::npos) << r_w;
            // Here refinement meets the bound everywhere; a triangle it
            // leaves under the bound would be warned.
            EXPECT_GE(r.mFaces[0].mBoundDeg, ang - 1e-9) << h << " " << ang;
            EXPECT_EQ(r.mMinAngleDeg < ang, warned(r)) << h << " " << ang;
            std::printf(
                "[stretched] h %-8.4g min_angle %4.1f  points %lld (unrefined %lld)  "
                "smallest %.2f  bound_deg %.2f  warned %d\n",
                h, ang, static_cast<long long>(r.mNumPoints),
                static_cast<long long>(unrefined.mNumPoints), r.mMinAngleDeg, r.mFaces[0].mBoundDeg,
                warned(r) ? 1 : 0);
        }
    }
}

TEST(MeshCadQuality, AFaceWhoseRefinementRunsAwayIsMeshedUnrefined) {
    // Found by the fuzz target (tests/fuzz/regressions/op_mesh_cad_surface):
    // the Scordelis-Lo roof with one corner's weight 2 and z lowered by 1,
    // refined to 30 degrees (past the 20.7 Ruppert guarantees), crowds the
    // edge u = 0 without end, and every round of boundary splits spent the
    // whole of max_points again. Refinement now stops at twice the face's
    // unrefined points; the face is meshed unrefined, byte for byte what
    // min_angle = 0 gives, and warned. At 25 degrees it refines normally.
    const char* roof = R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2,
        "swapped_surface_normal": false, "surface": {"is_trimmed": false, "is_rational": true,
        "degrees": [2, 1], "knot_vectors": [[0, 0, 1, 1], [0, 50]],
        "control_points": [[1, [-16.06969024216348, 0, 19.151111077974452, 2]],
          [2, [0, 0, 32.63518223330696, 0.766044443118978]],
          [3, [16.06969024216348, 0, 19.151111077974452, 1]],
          [4, [-16.06969024216348, 50, 18.151111077974452, 1]],
          [5, [0, 50, 32.63518223330696, 0.766044443118978]],
          [6, [16.06969024216348, 50, 19.151111077974452, 1]]]},
        "boundary_loops": [], "embedded_edges": []}], "edges": [], "vertices": []}]})";
    const std::string path = mt::temp_path(".cad.json");
    std::ofstream(path) << roof;
    const mio::BrepInfo info = mio::read_cad_json_info(path);
    const double h = 1.054790343148405;
    auto o = quality(h, 30.0, 2);
    o.mMaxChordalDeviation = 0.01 * h;
    const auto r = mio::mesh_cad_surface(info, o);
    ASSERT_EQ(r.mNumFailed, 0);
    bool warned = false;
    for (const std::string& r_w : r.mWarnings)
        warned = warned || r_w.find("meshed unrefined") != std::string::npos;
    EXPECT_TRUE(warned);
    EXPECT_EQ(r.mFaces[0].mNumRefinementPoints, 0);
    auto u = o;
    u.mMinAngleDeg = 0.0;
    EXPECT_EQ(digest(r.mMesh), digest(mio::mesh_cad_surface(info, u).mMesh));
    o.mMinAngleDeg = 25.0;
    const auto q = mio::mesh_cad_surface(info, o);
    EXPECT_GE(q.mMinAngleDeg, 25.0 - 1e-9);
    EXPECT_GT(q.mNumRefinementPoints, 0);
    for (const std::string& r_w : q.mWarnings)
        EXPECT_EQ(r_w.find("meshed unrefined"), std::string::npos) << r_w;
}

TEST(MeshCadQuality, ATrimLeavingItsSurfaceIsReported) {
    // Found by the fuzz target (tests/fuzz/regressions/op_mesh_cad_surface):
    // a trim that leaves its patch's knot range puts points where the patch
    // is evaluated clamped, off the face, and refinement there ran all its
    // rounds (87 s at 30 degrees). Such a face is unsupported, refined or
    // not; within a thousandth of the range it is still meshed.
    for (const double ang : {0.0, 30.0}) {
        const auto far = mio::mesh_cad_surface(planar_face({{{0, 0}, {4, 0}, {5, 4}, {0, 4}}}, 4.0),
                                               quality(0.25, ang));
        ASSERT_EQ(far.mNumFailed, 1) << ang;
        EXPECT_EQ(far.mFaces[0].mStatus, mio::MeshCadFaceStatus::Unsupported);
        EXPECT_NE(far.mFaces[0].mReason.find("parameter range"), std::string::npos);
        const auto near = mio::mesh_cad_surface(
            planar_face({{{0, 0}, {4.002, 0}, {4.002, 4}, {0, 4}}}, 4.0), quality(0.25, ang));
        EXPECT_EQ(near.mNumFailed, 0) << ang;
    }
}

TEST(MeshCadQuality, TriangulatePolygonTakesTheSameMinAngle) {
    // Roadmap 6.2.6: the planar triangulator refines too; 0 is unchanged.
    const std::vector<std::vector<Vec2>> l = {{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}};
    mio::TriangulatePolygon2dOptions o;
    o.mTargetSize = 0.3;
    const mio::Mesh a = mio::triangulate_polygon_2d(l, o);
    EXPECT_EQ(digest(a), digest(mio::triangulate_polygon_2d(l, 0.3)));
    o.mMinAngleDeg = 30.0;
    const mio::Mesh b = mio::triangulate_polygon_2d(l, o);
    double area = 0.0;
    for (const Tri& r_t : triangles(b)) {
        EXPECT_GE(r_t.mMinAngle, 30.0 - 1e-9);
        area += 0.5 * ((r_t.mX[1][0] - r_t.mX[0][0]) * (r_t.mX[2][1] - r_t.mX[0][1]) -
                       (r_t.mX[1][1] - r_t.mX[0][1]) * (r_t.mX[2][0] - r_t.mX[0][0]));
    }
    EXPECT_NEAR(area, 3.0, 1e-12);
    o.mMinAngleDeg = 61.0;
    EXPECT_THROW((void)mio::triangulate_polygon_2d(l, o), std::invalid_argument);
}

}  // namespace

namespace {

// A cone frustum (half-angle tan 0.5) from z = Z0 to z = 1, its radius
// varying 1 / Z0-fold: one open band.
mio::BrepInfo cone_frustum(double Z0) {
    std::vector<mio::BrepCurve> curves = {
        mio::detail::brep_curve_arc({0, 0, Z0}, {1, 0, 0}, {0, 1, 0}, 0.5 * Z0, 0.0, 2 * kPi),
        mio::detail::brep_curve_arc({0, 0, 1}, {1, 0, 0}, {0, 1, 0}, 0.5, 0.0, 2 * kPi)};
    std::vector<brep_solids::SolidFace> faces(1);
    faces[0].mSurface.mKind = 4;
    faces[0].mSurface.mTanHalfAngle = 0.5;
    faces[0].mLoops = {{brep_solids::circle_trim(0.5 * Z0, Z0, 2 * kPi, 0, true, false)},
                       {brep_solids::circle_trim(0.5, 1.0, -2 * kPi, 1, false, false)}};
    return brep_solids::assemble(faces, curves);
}

// A band of the torus of radii R and r between the tube's parallels at -Phi
// and Phi.
mio::BrepInfo torus_band_q(double R, double r, double Phi) {
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

// Roadmap 6.4.1: on a cone frustum whose radius varies tenfold, a torus
// band whose parallels do eightfold and the bilinear patch whose corner
// stretches its plane, one affine map per face cannot follow the surface.
// Under the mean metric the triangles there stay under the bound and are
// warned; with metric = varying, Lawson's flips and encroachment measured in
// the surface's own metric meet 25 degrees everywhere, unwarned, with no
// more points, the same at every thread count.
TEST(MeshCadVaryingMetric, StretchedFacesMeetTheBound) {
    struct Case {
        const char* mName;
        mio::BrepInfo mB;
        double mH;
        bool mMeanWarned;
    };
    std::vector<Case> cases;
    cases.push_back({"bilinear", bilinear_patch({1.0, 1.0}), 0.05, false});
    cases.push_back({"cone", cone_frustum(0.1), 0.1, true});
    cases.push_back({"cone", cone_frustum(0.05), 0.05, true});
    cases.push_back({"torus", torus_band_q(1.2, 1.0, 2.8), 0.2, true});
    const auto warned = [](const MeshCadSurfaceResult& rR) {
        for (const std::string& r_w : rR.mWarnings)
            if (r_w.find("roadmap 6.4.1") != std::string::npos)
                return true;
        return false;
    };
    for (const auto& r_c : cases) {
        const MeshCadSurfaceOptions mean = quality(r_c.mH, 25.0);
        MeshCadSurfaceOptions varying = mean;
        varying.mMetric = mio::MeshCadMetric::Varying;
        const auto m = mio::mesh_cad_surface(r_c.mB, mean);
        const auto v = mio::mesh_cad_surface(r_c.mB, varying);
        ASSERT_EQ(v.mNumFailed, 0) << r_c.mName;
        EXPECT_EQ(warned(m), r_c.mMeanWarned) << r_c.mName;
        EXPECT_FALSE(warned(v)) << r_c.mName;
        EXPECT_GE(v.mMinAngleDeg, 25.0 - 1e-9) << r_c.mName << " " << r_c.mH;
        EXPECT_LE(v.mNumPoints, m.mNumPoints + m.mNumPoints / 2) << r_c.mName;
        // Unrefined, the metric's flips alone lift the smallest angle.
        MeshCadSurfaceOptions plain = varying;
        plain.mMinAngleDeg = 0.0;
        MeshCadSurfaceOptions plain_mean = mean;
        plain_mean.mMinAngleDeg = 0.0;
        EXPECT_GT(mio::mesh_cad_surface(r_c.mB, plain).mMinAngleDeg,
                  mio::mesh_cad_surface(r_c.mB, plain_mean).mMinAngleDeg + 5.0)
            << r_c.mName;
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
        const std::uint64_t ref = digest(v.mMesh);
        const int saved = omp_get_max_threads();
        for (int threads : {1, 2, 4}) {
            omp_set_num_threads(threads);
            EXPECT_EQ(digest(mio::mesh_cad_surface(r_c.mB, varying).mMesh), ref) << threads;
        }
        omp_set_num_threads(saved);
#endif
    }
    // The solids of the suite mesh under it too, watertight.
    for (const auto& r_b :
         {brep_solids::box(), brep_solids::capped_cylinder(1.0, 2.0), brep_solids::sphere(1.0)}) {
        MeshCadSurfaceOptions o = quality(0.25, 25.0);
        o.mMetric = mio::MeshCadMetric::Varying;
        const auto r = mio::mesh_cad_surface(r_b, o);
        EXPECT_EQ(r.mNumFailed, 0);
        EXPECT_EQ(r.mNumOpenEdges, 0);
        EXPECT_EQ(r.mNumMisorientedEdges, 0);
    }
}
