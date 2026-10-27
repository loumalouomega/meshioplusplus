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
// mesh_cad_surface's discrete reparametrisation (roadmap 6.4.4): a face
// whose parameters no metric rescues is meshed in a flattening of its own
// coarse mesh, every point lifted with the face's surface.

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
#include "../../src/cpp/src/detail/mesh_cad_reparam.hpp"
#include "../../src/cpp/src/detail/param_disc.hpp"
#include "../../src/cpp/src/operations/mesh_cad_testing.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::MeshCadSurfaceOptions;
using mio::MeshCadSurfaceResult;

constexpr double kPi = 3.14159265358979323846;

MeshCadSurfaceOptions quality(double H, double MinAngle) {
    MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    o.mMinAngleDeg = MinAngle;
    return o;
}

std::uint64_t digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

double angle_at(const double* pO, const double* pP, const double* pQ) {
    const double u[3] = {pP[0] - pO[0], pP[1] - pO[1], pP[2] - pO[2]};
    const double v[3] = {pQ[0] - pO[0], pQ[1] - pO[1], pQ[2] - pO[2]};
    const double c[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                         u[0] * v[1] - u[1] * v[0]};
    return std::atan2(std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]),
                      u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) *
           180.0 / kPi;
}

// The smallest angle of the triangles, and the area.
std::array<double, 2> angle_and_area(const mio::Mesh& rM) {
    const double* x = rM.Points().As<double>();
    const auto cells = rM.Cells(0);
    const mio::NDArray& c = cells.Conn();
    double m = 180.0, area = 0.0;
    for (std::size_t t = 0; t < cells.NumCells(); ++t) {
        const double* p[3];
        for (std::size_t k = 0; k < 3; ++k)
            p[k] = x + 3 * mio::detail::read_int(c, 3 * t + k);
        m = std::min({m, angle_at(p[0], p[1], p[2]), angle_at(p[1], p[2], p[0]),
                      angle_at(p[2], p[0], p[1])});
        const double u[3] = {p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]};
        const double v[3] = {p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]};
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                             u[0] * v[1] - u[1] * v[0]};
        area += 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    }
    return {m, area};
}

mio::BrepInfo read_json(const std::string& rJson) {
    const std::string path = mt::temp_path(".cad.json");
    std::ofstream(path) << rJson;
    return mio::read_cad_json_info(path);
}

// A patch of degrees (2, 2) on [0, 1]^2 from its nine control points (u
// fastest), untrimmed or trimmed to the parameter polygon `rTrim`.
mio::BrepInfo quadratic_patch(const std::array<std::array<double, 3>, 9>& rP,
                              const std::vector<std::array<double, 2>>& rTrim = {}) {
    std::ostringstream j;
    j.precision(17);
    j << R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "swapped_surface_normal": false,
        "surface": {"is_trimmed": )"
      << (rTrim.empty() ? "false" : "true") << R"(, "is_rational": false, "degrees": [2, 2],
        "knot_vectors": [[0, 0, 1, 1], [0, 0, 1, 1]], "control_points": [)";
    for (std::size_t i = 0; i < 9; ++i)
        j << (i ? ", " : "") << "[" << i + 1 << ", [" << rP[i][0] << ", " << rP[i][1] << ", "
          << rP[i][2] << ", 1]]";
    j << R"(]}, "boundary_loops": [)";
    if (!rTrim.empty()) {
        j << R"({"loop_type": "outer", "trimming_curves": [)";
        for (std::size_t k = 0; k < rTrim.size(); ++k) {
            const auto& a = rTrim[k];
            const auto& b = rTrim[(k + 1) % rTrim.size()];
            j << (k ? ", " : "") << R"({"trim_index": )" << k
              << R"(, "curve_direction": true, "parameter_curve": {"degree": 1,
                "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[)"
              << a[0] << ", " << a[1] << R"(, 0, 1]], [[)" << b[0] << ", " << b[1]
              << R"(, 0, 1]]]}})";
        }
        j << "]}";
    }
    j << R"(], "embedded_edges": []}], "edges": [], "vertices": []}]})";
    return read_json(j.str());
}

// The planar bilinear patch with corners (0, 0), (1, 0), (0, 2) and (0, 2):
// two corner control points coincide, so its side v = 1 collapses and the
// patch is the triangle (0, 0), (1, 0), (0, 2), its corner there 26.6
// degrees. Written as a quadratic patch: the mid control points on the
// edges.
mio::BrepInfo collapsed_plane() {
    std::array<std::array<double, 3>, 9> p{};
    for (std::size_t jv = 0; jv < 3; ++jv)
        for (std::size_t iu = 0; iu < 3; ++iu) {
            const double u = 0.5 * double(iu), v = 0.5 * double(jv);
            // x = u (1 - v), y = 2 v: bilinear, so its quadratic net is the
            // bilinear one sampled at the halves.
            p[3 * jv + iu] = {u * (1.0 - v), 2.0 * v, 0.0};
        }
    return quadratic_patch(p);
}

// A curved patch whose side v = 1 collapses: a paraboloid cap's sector.
mio::BrepInfo collapsed_dome() {
    std::array<std::array<double, 3>, 9> p{};
    for (std::size_t jv = 0; jv < 3; ++jv)
        for (std::size_t iu = 0; iu < 3; ++iu) {
            const double u = 0.5 * double(iu), v = 0.5 * double(jv);
            p[3 * jv + iu] = {2.0 * u * (1.0 - v), 2.0 * v, 0.6 * u * (1.0 - u) * (1.0 - v)};
        }
    return quadratic_patch(p);
}

void report(const char* pName, const MeshCadSurfaceResult& rR,
            const mio::MeshCadSurfaceExtrasReport* pX = nullptr) {
    const auto aa = angle_and_area(rR.mMesh);
    std::printf("[6.4.4] %-26s failed %lld points %lld min angle %.2f area %.6f warnings %zu",
                pName, static_cast<long long>(rR.mNumFailed), static_cast<long long>(rR.mNumPoints),
                aa[0], aa[1], rR.mWarnings.size());
    if (pX)
        for (const auto& r_f : pX->mReparametrised)
            std::printf("  [face %lld cond %.3g disc %.3g taken %d method %d %s]",
                        static_cast<long long>(r_f.mFace), r_f.mConditionNumber,
                        r_f.mConditionNumberDisc, r_f.mTaken ? 1 : 0, r_f.mMethod,
                        r_f.mReason.c_str());
    std::printf("\n");
    for (const auto& r_w : rR.mWarnings)
        std::printf("      warning: %s\n", r_w.c_str());
}

MeshCadSurfaceResult run(const mio::BrepInfo& rB, const MeshCadSurfaceOptions& rO,
                         mio::MeshCadSurfaceExtrasReport* pX = nullptr) {
    return mio::mesh_cad_surface(rB, rO, mio::MeshCadSizeSources{}, mio::MeshCadSurfaceExtras{},
                                 nullptr, pX);
}

// Three plates in the plane z = 0, each its own patch: [0, 1] x [0, 1], the
// sliver [0, 1] x [1, 1 + W] and [0, 1] x [1 + W, 2 + W], sharing their
// edges (roadmap 6.5.2). With `Fold` the third plate turns up 90 degrees
// about y = 1 + W, so the sliver meets it at a crease.
mio::BrepInfo sliver_strip(double W, bool Fold = false) {
    using brep_solids::line_trim;
    using brep_solids::PatchVec;
    const double y1 = 1.0, y2 = 1.0 + W;
    std::vector<PatchVec> v = {{0, 0, 0},  {1, 0, 0},  {0, y1, 0}, {1, y1, 0},
                               {0, y2, 0}, {1, y2, 0}, {0, y2, 0}, {1, y2, 0}};
    v[6] = Fold ? PatchVec{0, y2, 1} : PatchVec{0, y2 + 1, 0};
    v[7] = Fold ? PatchVec{1, y2, 1} : PatchVec{1, y2 + 1, 0};
    std::vector<std::array<int, 2>> ends;
    std::vector<mio::BrepCurve> curves;
    const auto edge = [&](int A, int B) {
        const int lo = std::min(A, B), hi = std::max(A, B);
        for (std::size_t e = 0; e < ends.size(); ++e)
            if (ends[e][0] == lo && ends[e][1] == hi)
                return static_cast<std::int64_t>(e);
        ends.push_back({lo, hi});
        curves.push_back(mio::detail::brep_curve_line(v[static_cast<std::size_t>(lo)],
                                                      v[static_cast<std::size_t>(hi)]));
        return static_cast<std::int64_t>(ends.size() - 1);
    };
    std::vector<brep_solids::SolidFace> faces;
    const auto quad = [&](std::array<int, 4> rC, const PatchVec& rY, const PatchVec& rZ) {
        brep_solids::SolidFace f;
        f.mSurface = brep_solids::plane(v[static_cast<std::size_t>(rC[0])], {1, 0, 0}, rY, rZ);
        f.mLoops.emplace_back();
        for (std::size_t k = 0; k < 4; ++k) {
            const int a = rC[k], b = rC[(k + 1) % 4];
            f.mLoops[0].push_back(line_trim(v[static_cast<std::size_t>(a)],
                                            v[static_cast<std::size_t>(b)], edge(a, b), a < b));
        }
        faces.push_back(f);
    };
    quad({0, 1, 3, 2}, {0, 1, 0}, {0, 0, 1});
    quad({2, 3, 5, 4}, {0, 1, 0}, {0, 0, 1});
    if (Fold)
        quad({4, 5, 7, 6}, {0, 0, 1}, {0, -1, 0});
    else
        quad({4, 5, 7, 6}, {0, 1, 0}, {0, 0, 1});
    return brep_solids::assemble(faces, curves);
}

}  // namespace

namespace {

// Two bands of radius 1 stacked from z = 0 to 1 and 1 to 2, the upper's seam
// `Turn` radians round from the lower's, so the shared rim's two uses start
// apart and neither can be rotated (each runs between its band's seams). The
// rim is edge 1.
mio::BrepInfo stacked_bands(double Turn) {
    std::vector<mio::BrepCurve> curves = {
        mio::detail::brep_curve_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1.0, 0.0, 2 * kPi),
        mio::detail::brep_curve_arc({0, 0, 1}, {1, 0, 0}, {0, 1, 0}, 1.0, 0.0, 2 * kPi),
        mio::detail::brep_curve_arc({0, 0, 2}, {1, 0, 0}, {0, 1, 0}, 1.0, 0.0, 2 * kPi)};
    std::vector<brep_solids::SolidFace> faces(2);
    faces[0].mSurface.mKind = 3;
    faces[0].mSurface.mRadius = 1.0;
    faces[0].mLoops = {{brep_solids::circle_trim(1.0, 0.0, 2 * kPi, 0, true, false)},
                       {brep_solids::circle_trim(1.0, 1.0, -2 * kPi, 1, false, false)}};
    faces[1].mSurface.mKind = 3;
    faces[1].mSurface.mRadius = 1.0;
    faces[1].mSurface.mX = {std::cos(Turn), std::sin(Turn), 0};
    faces[1].mSurface.mY = {-std::sin(Turn), std::cos(Turn), 0};
    faces[1].mLoops = {{brep_solids::circle_trim(1.0, 1.0, 2 * kPi, 1, true, false, Turn)},
                       {brep_solids::circle_trim(1.0, 2.0, -2 * kPi, 2, false, false, Turn)}};
    return brep_solids::assemble(faces, curves);
}

// Whether every triangle edge used once lies on the strip's outer rectangle
// [0, 1] x [0, Y].
bool open_only_outside(const mio::Mesh& rM, double Y) {
    const double* x = rM.Points().As<double>();
    const auto cells = rM.Cells(0);
    const mio::NDArray& c = cells.Conn();
    std::map<std::pair<std::int64_t, std::int64_t>, int> uses;
    for (std::size_t t = 0; t < cells.NumCells(); ++t)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t a = mio::detail::read_int(c, 3 * t + k);
            const std::int64_t b = mio::detail::read_int(c, 3 * t + (k + 1) % 3);
            ++uses[{std::min(a, b), std::max(a, b)}];
        }
    const auto outside = [&](std::int64_t P) {
        const double* p = x + 3 * P;
        return std::fabs(p[0]) < 1e-12 || std::fabs(p[0] - 1.0) < 1e-12 ||
               std::fabs(p[1]) < 1e-12 || std::fabs(p[1] - Y) < 1e-12;
    };
    for (const auto& [e, n] : uses)
        if (n == 1 && !(outside(e.first) && outside(e.second)))
            return false;
    return true;
}

}  // namespace

// Roadmap 6.4.4's probe: a patch whose corner control points coincide (its
// side v = 1 collapses onto one point), planar and curved. In its own
// parameters the collapsed corner is a pole whose fan neither metric rounds:
// the mean metric leaves triangles under 25 degrees and the varying one, at
// these sizes, too. Meshed in the flattening of its coarse mesh (the
// least-squares conformal map), the corner is a corner, and every triangle
// meets 25 degrees, unwarned, with fewer points.
TEST(MeshCadReparam, ACollapsedPatchMeetsTheBound) {
    struct Case {
        const char* mName;
        mio::BrepInfo mB;
    };
    std::vector<Case> cases;
    cases.push_back({"plane", collapsed_plane()});
    cases.push_back({"dome", collapsed_dome()});
    for (const auto& r_c : cases)
        for (const double h : {0.2, 0.1}) {
            MeshCadSurfaceOptions o = quality(h, 25.0);
            const auto mean = run(r_c.mB, o);
            o.mMetric = mio::MeshCadMetric::Varying;
            const auto varying = run(r_c.mB, o);
            o.mMetric = mio::MeshCadMetric::Mean;
            o.mReparamConditionNumber = 4.0;
            mio::MeshCadSurfaceExtrasReport x;
            const auto r = run(r_c.mB, o, &x);
            report(r_c.mName, r, &x);
            EXPECT_LT(angle_and_area(mean.mMesh)[0], 25.0) << r_c.mName << " " << h;
            EXPECT_LT(angle_and_area(varying.mMesh)[0], 25.0) << r_c.mName << " " << h;
            EXPECT_EQ(r.mNumFailed, 0);
            EXPECT_TRUE(r.mWarnings.empty()) << r_c.mName << " " << h;
            EXPECT_GE(angle_and_area(r.mMesh)[0], 25.0 - 1e-9) << r_c.mName << " " << h;
            ASSERT_EQ(x.mReparametrised.size(), 1u);
            EXPECT_TRUE(x.mReparametrised[0].mTaken);
            EXPECT_EQ(x.mReparametrised[0].mMethod, 1);
            EXPECT_GT(x.mReparametrised[0].mConditionNumber, 4.0);
            EXPECT_LT(x.mReparametrised[0].mConditionNumberDisc,
                      x.mReparametrised[0].mConditionNumber);
            // The area is the surface's: the plane's triangle is exactly 1.
            if (std::string(r_c.mName) == "plane")
                EXPECT_NEAR(angle_and_area(r.mMesh)[1], 1.0, 1e-12);
        }
}

// The same patch trimmed short of its collapsed side: no pole, so no
// exemption, only a corner the parameters stretch thirtyfold. The mean
// metric leaves a 1-degree triangle there and warns; the flattening meets
// 25 degrees with fewer points than the varying metric.
TEST(MeshCadReparam, ACollapsedCornerOffTheTrimMeetsTheBound) {
    std::array<std::array<double, 3>, 9> p{};
    for (std::size_t jv = 0; jv < 3; ++jv)
        for (std::size_t iu = 0; iu < 3; ++iu) {
            const double u = 0.5 * double(iu), v = 0.5 * double(jv);
            p[3 * jv + iu] = {u * (1.0 - v), 2.0 * v, 0.0};
        }
    const mio::BrepInfo b = quadratic_patch(p, {{0, 0}, {1, 0}, {0.03, 0.97}, {0, 0.97}});
    for (const double h : {0.2, 0.1}) {
        MeshCadSurfaceOptions o = quality(h, 25.0);
        const auto mean = run(b, o);
        EXPECT_LT(angle_and_area(mean.mMesh)[0], 5.0);
        EXPECT_FALSE(mean.mWarnings.empty());
        o.mMetric = mio::MeshCadMetric::Varying;
        const auto varying = run(b, o);
        o.mMetric = mio::MeshCadMetric::Mean;
        o.mReparamConditionNumber = 4.0;
        mio::MeshCadSurfaceExtrasReport x;
        const auto r = run(b, o, &x);
        EXPECT_TRUE(r.mWarnings.empty());
        EXPECT_GE(angle_and_area(r.mMesh)[0], 25.0 - 1e-9) << h;
        EXPECT_LT(r.mNumPoints, varying.mNumPoints) << h;
        ASSERT_EQ(x.mReparametrised.size(), 1u);
        EXPECT_TRUE(x.mReparametrised[0].mTaken);
        // Every point on the patch: the plane z = 0 inside its trim.
        const double* px = r.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
            EXPECT_EQ(px[3 * i + 2], 0.0);
    }
}

// The sphere's two hemispheres, each closed through a seam onto a pole, are
// flattened like any other face: the pole is then a point of the seam's cut,
// not a fan's apex, and every triangle meets 25 degrees where the fans and
// even two pole rings do not at h = 0.05. Every point is on the sphere, the
// mesh is watertight, and it is the same over thread counts and a shuffled
// face order.
TEST(MeshCadReparam, TheSphereIsFlattenedThroughItsPoles) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    for (const double h : {0.2, 0.1, 0.05}) {
        MeshCadSurfaceOptions o = quality(h, 25.0);
        o.mReparamConditionNumber = 4.0;
        mio::MeshCadSurfaceExtrasReport x;
        const auto r = run(sphere, o, &x);
        EXPECT_EQ(r.mNumFailed, 0);
        EXPECT_EQ(r.mNumOpenEdges, 0);
        EXPECT_EQ(r.mNumNonManifoldEdges, 0);
        EXPECT_EQ(r.mNumMisorientedEdges, 0);
        EXPECT_GE(angle_and_area(r.mMesh)[0], 25.0 - 1e-9) << h;
        ASSERT_EQ(x.mReparametrised.size(), 2u);
        for (const auto& r_f : x.mReparametrised)
            EXPECT_TRUE(r_f.mTaken) << r_f.mReason;
        const double* px = r.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
            EXPECT_NEAR(std::sqrt(px[3 * i] * px[3 * i] + px[3 * i + 1] * px[3 * i + 1] +
                                  px[3 * i + 2] * px[3 * i + 2]),
                        1.0, 1e-12);
        // Area to second order: 4 pi within the chordal error.
        EXPECT_NEAR(angle_and_area(r.mMesh)[1], 4.0 * kPi, 2.0 * h * h);
    }
    MeshCadSurfaceOptions o = quality(0.2, 25.0);
    o.mReparamConditionNumber = 4.0;
    const std::uint64_t ref = digest(run(sphere, o).mMesh);
    for (const std::uint64_t seed : {1u, 7u}) {
        mio::detail::mesh_cad_set_face_order_seed(seed);
        EXPECT_EQ(digest(run(sphere, o).mMesh), ref) << seed;
    }
    mio::detail::mesh_cad_set_face_order_seed(0);
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
        EXPECT_EQ(digest(run(sphere, o).mMesh), ref) << threads;
    }
    omp_set_num_threads(saved);
#endif
}

// A threshold no face passes changes nothing: the mesh is the option off's
// byte for byte (the coarse snapshot leaves no state behind), and faces
// whose metric is the plane's (a box, a cylinder) are never flattened.
TEST(MeshCadReparam, BelowTheThresholdNothingChanges) {
    for (const auto& r_b : {collapsed_plane(), brep_solids::sphere(1.0), brep_solids::box(),
                            brep_solids::capped_cylinder(1.0, 2.0)}) {
        MeshCadSurfaceOptions o = quality(0.2, 25.0);
        const std::uint64_t off = digest(run(r_b, o).mMesh);
        o.mReparamConditionNumber = 1e300;
        mio::MeshCadSurfaceExtrasReport x;
        EXPECT_EQ(digest(run(r_b, o, &x).mMesh), off);
        EXPECT_TRUE(x.mReparametrised.empty());
    }
    // A plane's metric is its plane's exactly; a rational circle's speed
    // varies, so a cylinder's condition number is about 1.4.
    {
        MeshCadSurfaceOptions o = quality(0.2, 25.0);
        const std::uint64_t off = digest(run(brep_solids::box(), o).mMesh);
        o.mReparamConditionNumber = 1.0 + 1e-9;
        EXPECT_EQ(digest(run(brep_solids::box(), o).mMesh), off);
        o.mReparamConditionNumber = 2.0;
        const auto cyl = brep_solids::capped_cylinder(1.0, 2.0);
        EXPECT_EQ(digest(run(cyl, o).mMesh), digest(run(cyl, quality(0.2, 25.0)).mMesh));
    }
}

// Collapsed coordinates (Duffy 1982): in the coarse triangles that keep a
// corner at a pole, the parameters along an edge out of the pole are the
// same from both sides, and the strip the dropped triangle spanned is
// covered: here the unit square's side v = 1 collapses, the pole's segment
// is [(0, 1), (1, 1)], and the two triangles (P, Q, R) and (P, S, Q) share
// the edge P-Q.
TEST(MeshCadReparam, CollapsedCoordinatesAreContinuousAtAPole) {
    mio::detail::McadCoarseMesh c;
    c.mX = {{0.5, 1.0, 0.0}, {0.5, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    c.mId = {0, 1, 2, 3};
    c.mTris = {{0, 1, 2}, {0, 3, 1}};
    c.mCornerUv = {{{{0.6, 1.0}, {0.5, 0.0}, {1.0, 0.0}}}, {{{0.4, 1.0}, {0.0, 0.0}, {0.5, 0.0}}}};
    c.mOwner = {0, 0};
    c.mPoleCorner = {0, 0};
    c.mReversed = {0};
    c.mPoleSegment[{0, 0}] = {{{0.0, 1.0}, {1.0, 1.0}}};
    mio::detail::McadDisc d;
    d.pCoarse = &c;
    d.mP = {{0.5, 1.0}, {0.5, 0.0}, {1.0, 0.0}, {0.0, 0.0}};
    d.Index();
    // Along P-Q (corners 0 and 1 of the first, 0 and 2 of the second).
    for (const double t : {0.1, 0.5, 0.9}) {
        const auto a = d.UvAt(0, {1.0 - t, t, 0.0});
        const auto b = d.UvAt(1, {1.0 - t, 0.0, t});
        EXPECT_NEAR(a[0], b[0], 1e-15) << t;
        EXPECT_NEAR(a[1], b[1], 1e-15) << t;
        // On the line u = 0.5 from the pole's point (0.5, 1) to Q.
        EXPECT_NEAR(a[0], 0.5, 1e-15);
        EXPECT_NEAR(a[1], t * 0.0 + (1.0 - t) * 1.0, 1e-15);
    }
    // At the pole every side reaches its own point of the segment: the edge
    // P-R ends at u = 1, P-S at u = 0.
    EXPECT_NEAR(d.UvAt(0, {1.0 - 1e-9, 0.0, 1e-9})[0], 1.0, 1e-6);
    EXPECT_NEAR(d.UvAt(1, {1.0 - 1e-9, 1e-9, 0.0})[0], 0.0, 1e-6);
    // Located through the disc.
    const auto hit = d.Locate({0.75, 0.25});
    EXPECT_TRUE(hit.mInside);
    EXPECT_EQ(hit.mTri, 0u);
}

// A point in a bucket no triangle's box meets is outside the flattening: it
// is located against the nearest ring of buckets that holds any, with its
// clamped fallback, not against every triangle (found by the fuzz target:
// each lattice point outside a curved disc cost the whole coarse mesh).
TEST(MeshCadReparam, APointOutsideTheDiscIsLocatedNearby) {
    mio::detail::McadCoarseMesh c;
    // A strip of 20 triangles along y in [0, 1] and one triangle far above.
    for (int i = 0; i <= 10; ++i) {
        c.mX.push_back({double(i), 0.0, 0.0});
        c.mX.push_back({double(i), 1.0, 0.0});
    }
    c.mX.push_back({0.0, 9.0, 0.0});
    c.mX.push_back({1.0, 9.0, 0.0});
    c.mX.push_back({0.0, 10.0, 0.0});
    for (std::size_t i = 0; i < 10; ++i) {
        c.mTris.push_back({2 * i, 2 * i + 2, 2 * i + 3});
        c.mTris.push_back({2 * i, 2 * i + 3, 2 * i + 1});
    }
    c.mTris.push_back({22, 23, 24});
    mio::detail::McadDisc d;
    d.pCoarse = &c;
    for (const auto& r_x : c.mX) {
        c.mId.push_back(static_cast<std::int64_t>(c.mId.size()));
        d.mP.push_back({r_x[0], r_x[1]});
    }
    for (const auto& r_t : c.mTris) {
        std::array<std::array<double, 2>, 3> uv{};
        for (std::size_t k = 0; k < 3; ++k)
            uv[k] = {c.mX[r_t[k]][0], c.mX[r_t[k]][1]};
        c.mCornerUv.push_back(uv);
        c.mOwner.push_back(0);
        c.mPoleCorner.push_back(-1);
    }
    c.mReversed = {0};
    d.Index();  // about two triangles a bucket: a 4 by 4 grid, empty in the middle
    // Inside: found, its parameters the point's own.
    const auto in = d.Locate({4.25, 0.5});
    EXPECT_TRUE(in.mInside);
    EXPECT_NEAR(in.mUv[0], 4.25, 1e-12);
    EXPECT_NEAR(in.mUv[1], 0.5, 1e-12);
    // Far from every triangle: outside, clamped onto one nearby.
    const auto out = d.Locate({6.0, 5.0});
    EXPECT_FALSE(out.mInside);
    ASSERT_LT(out.mTri, c.mTris.size());
    EXPECT_GE(out.mUv[0], 0.0);
    EXPECT_LE(out.mUv[0], 10.0);
    EXPECT_LE(out.mUv[1], 10.0);
}

// Roadmap 6.5.2: three plates in a row, the middle one h / 20 wide. Alone,
// the sliver is a strip of needles and its neighbours are refined against
// it; merged into the neighbour it meets along the most edge, the pair is
// one triangulation that meets 25 degrees with under a third of the
// triangles. Each triangle keeps the face its centroid lies in, the strip is
// open only along its outside, and the merge is reported.
TEST(MeshCadVirtualFaces, ASliverStripIsMeshedWithItsNeighbour) {
    for (const double h : {0.2, 0.1}) {
        const double w = h / 20.0;
        const mio::BrepInfo b = sliver_strip(w);
        const auto alone = run(b, quality(h, 25.0));
        mio::MeshCadSurfaceExtras extras;
        extras.mVirtualFaces = true;
        mio::MeshCadSurfaceExtrasReport x;
        const auto r = mio::mesh_cad_surface(b, quality(h, 25.0), mio::MeshCadSizeSources{}, extras,
                                             nullptr, &x);
        EXPECT_EQ(r.mNumFailed, 0);
        EXPECT_GE(angle_and_area(r.mMesh)[0], 25.0 - 1e-9) << h;
        EXPECT_LT(3 * r.mNumTriangles, alone.mNumTriangles) << h;
        EXPECT_NEAR(angle_and_area(r.mMesh)[1], 2.0 + w, 1e-12);
        EXPECT_EQ(r.mNumMisorientedEdges, 0);
        EXPECT_EQ(r.mNumNonManifoldEdges, 0);
        EXPECT_TRUE(open_only_outside(r.mMesh, 2.0 + w));
        ASSERT_EQ(x.mMergedFaces.size(), 1u);
        const auto& r_m = x.mMergedFaces[0];
        EXPECT_EQ(r_m.mSliver, 1);
        EXPECT_TRUE(r_m.mMerged) << r_m.mReason;
        EXPECT_NEAR(r_m.mWidth, 2.0 * w / (2.0 + 2.0 * w), 1e-12);
        ASSERT_EQ(r_m.mEdges.size(), 1u);
        // cad:face: every face has triangles, each inside its own face.
        const double* px = r.mMesh.Points().As<double>();
        const auto cells = r.mMesh.Cells(0);
        const mio::NDArray& conn = cells.Conn();
        const mio::NDArray& face = r.mMesh.CellData(mio::kMeshCadFaceName, 0);
        std::array<std::int64_t, 3> count{0, 0, 0};
        const std::array<double, 4> ys{0.0, 1.0, 1.0 + w, 2.0 + w};
        for (std::size_t t = 0; t < cells.NumCells(); ++t) {
            const std::int64_t f = mio::detail::read_int(face, t);
            ASSERT_GE(f, 0);
            ASSERT_LT(f, 3);
            ++count[static_cast<std::size_t>(f)];
            double y = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
                y += px[3 * mio::detail::read_int(conn, 3 * t + k) + 1] / 3.0;
            EXPECT_GE(y, ys[static_cast<std::size_t>(f)] - 1e-12);
            EXPECT_LE(y, ys[static_cast<std::size_t>(f) + 1] + 1e-12);
        }
        for (const std::int64_t n : count)
            EXPECT_GT(n, 0);
        EXPECT_EQ(r.mFaces[1].mNumTriangles, count[1]);
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
        const std::uint64_t ref = digest(r.mMesh);
        const int saved = omp_get_max_threads();
        for (int threads : {1, 2, 4}) {
            omp_set_num_threads(threads);
            EXPECT_EQ(
                digest(mio::mesh_cad_surface(b, quality(h, 25.0), mio::MeshCadSizeSources{}, extras)
                           .mMesh),
                ref)
                << threads;
        }
        omp_set_num_threads(saved);
#endif
    }
}

// A sliver at a crease (its far neighbour turned up 90 degrees) still
// merges with the flat one; folded on both sides, it is declined, and the
// mesh is the option off's. Order 2 merges too (roadmap 6.5.2): its nodes
// lie on the strip's planes.
TEST(MeshCadVirtualFaces, ACreaseIsDeclinedAndOptionsAreChecked) {
    const double h = 0.2;
    const mio::BrepInfo b = sliver_strip(h / 20.0, true);
    mio::MeshCadSurfaceExtras extras;
    extras.mVirtualFaces = true;
    mio::MeshCadSurfaceExtrasReport x;
    const auto r =
        mio::mesh_cad_surface(b, quality(h, 25.0), mio::MeshCadSizeSources{}, extras, nullptr, &x);
    ASSERT_EQ(x.mMergedFaces.size(), 1u);
    if (x.mMergedFaces[0].mMerged) {
        EXPECT_EQ(x.mMergedFaces[0].mHost, 0);
    } else {
        EXPECT_NE(x.mMergedFaces[0].mReason.find("tangentially"), std::string::npos);
        EXPECT_EQ(digest(r.mMesh), digest(run(b, quality(h, 25.0)).mMesh));
    }
    MeshCadSurfaceOptions o = quality(h, 0.0);
    o.mOrder = 2;
    const auto r2 = mio::mesh_cad_surface(b, o, mio::MeshCadSizeSources{}, extras);
    EXPECT_EQ(r2.mNumFailed, 0);
    EXPECT_EQ(r2.mMesh.Cells(0).Type(), "triangle6");
}

// Roadmap 6.5.1: two stacked bands whose shared rim's uses start 1e-5 apart
// and neither of which can be rotated: the global tolerance (1e-7 of the
// diagonal) fails the face whose use is not the anchor. The rim declaring a tolerance of 1e-4, as
// an OpenCASCADE read does, the two starts are one point and the faces mesh; the use is reported.
// Without the option, the declaration is ignored.
TEST(MeshCadTolerances, ADeclaredToleranceJoinsAClosedEdgesStarts) {
    mio::BrepInfo b = stacked_bands(1e-5);
    b.mEdges[1].mExtra = R"({"tolerance":0.0001,"end_tolerances":[0.0001,0.0001]})";
    const MeshCadSurfaceOptions o = quality(0.2, 0.0);
    const auto off = run(b, o);
    EXPECT_EQ(off.mNumFailed, 1);
    mio::MeshCadSurfaceExtras extras;
    extras.mModelTolerances = true;
    mio::MeshCadSurfaceExtrasReport x;
    const auto r = mio::mesh_cad_surface(b, o, mio::MeshCadSizeSources{}, extras, nullptr, &x);
    EXPECT_EQ(r.mNumFailed, 0);
    EXPECT_EQ(r.mNumMisorientedEdges, 0);
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    ASSERT_FALSE(x.mToleranceUses.empty());
    EXPECT_EQ(x.mToleranceUses[0].mEdge, 1);
    EXPECT_EQ(x.mToleranceUses[0].mWhere, "closed_start");
    EXPECT_NEAR(x.mToleranceUses[0].mGap, 1e-5, 1e-9);
    EXPECT_EQ(x.mToleranceUses[0].mDeclared, 1e-4);
    // A gap past the declaration still fails; a bad declaration is warned.
    mio::BrepInfo far = stacked_bands(1e-3);
    far.mEdges[1].mExtra = R"({"tolerance":0.0001})";
    EXPECT_EQ(mio::mesh_cad_surface(far, o, mio::MeshCadSizeSources{}, extras).mNumFailed, 1);
    far.mEdges[1].mExtra = R"({"tolerance":-1})";
    const auto bad = mio::mesh_cad_surface(far, o, mio::MeshCadSizeSources{}, extras);
    EXPECT_TRUE(std::any_of(bad.mWarnings.begin(), bad.mWarnings.end(), [](const std::string& rW) {
        return rW.find("declared tolerance") != std::string::npos;
    }));
    // A use that can be rotated is rotated as without the option, though its
    // start lies within the declaration: a declared tolerance accepts only
    // what the global one refuses (on the corpus, keeping such a cap's start
    // folded a triangle).
    {
        mio::BrepInfo cyl = brep_solids::capped_cylinder(1.0, 2.0, 1e-5);
        const std::uint64_t plain = digest(run(cyl, o).mMesh);
        for (auto& r_e : cyl.mEdges)
            r_e.mExtra = R"({"tolerance":0.0001,"end_tolerances":[0.0001,0.0001]})";
        mio::MeshCadSurfaceExtrasReport xc;
        const auto rc = mio::mesh_cad_surface(cyl, o, mio::MeshCadSizeSources{}, extras, nullptr, &xc);
        EXPECT_EQ(digest(rc.mMesh), plain);
        EXPECT_TRUE(xc.mToleranceUses.empty());
    }
    // With nothing declared, the option changes nothing.
    const mio::BrepInfo box = brep_solids::box();
    EXPECT_EQ(
        digest(mio::mesh_cad_surface(box, quality(0.3, 25.0), mio::MeshCadSizeSources{}, extras)
                   .mMesh),
        digest(run(box, quality(0.3, 25.0)).mMesh));
}

// Roadmap 6.5.1: the box without edge topology, one face translated by
// 1e-5 (well past the global tolerance, 1e-7 of the diagonal): its four
// trims do not pair with their neighbours'. The file's own
// `tolerances.model_tolerance` (Kratos CAD JSON) pairs them, and each pair
// is reported with no `BrepEdge`. Without the option it is ignored.
TEST(MeshCadTolerances, TheModelsTolerancePairsTheTrimsNoEdgeNames) {
    mio::BrepInfo b = brep_solids::box();
    b.mEdges.clear();
    mio::NDArray& r_cp = b.mPatches[static_cast<std::size_t>(b.mFaces[0].mPatch)].mControlPoints;
    for (std::size_t i = 0; i < r_cp.Size(); ++i)
        r_cp.As<double>()[i] += 1e-5;
    b.mExtra = R"({"tolerances":{"model_tolerance":0.0001},"version_number":1})";
    const MeshCadSurfaceOptions o = quality(0.2, 0.0);
    const auto off = run(b, o);
    EXPECT_EQ(off.mNumEdges, 16);
    mio::MeshCadSurfaceExtras extras;
    extras.mModelTolerances = true;
    mio::MeshCadSurfaceExtrasReport x;
    const auto r = mio::mesh_cad_surface(b, o, mio::MeshCadSizeSources{}, extras, nullptr, &x);
    EXPECT_EQ(r.mNumFailed, 0);
    EXPECT_EQ(r.mNumEdges, 12);
    EXPECT_EQ(r.mNumOpenEdges, 0);
    ASSERT_EQ(x.mToleranceUses.size(), 4u);
    for (const auto& r_u : x.mToleranceUses) {
        EXPECT_EQ(r_u.mEdge, -1);
        EXPECT_EQ(r_u.mWhere, "pairing");
        EXPECT_NEAR(r_u.mGap, std::sqrt(3.0) * 1e-5, 1e-9);
        EXPECT_EQ(r_u.mDeclared, 1e-4);
    }
    // A model tolerance under the gap pairs nothing; a malformed one is warned.
    b.mExtra = R"({"tolerances":{"model_tolerance":1e-6}})";
    EXPECT_EQ(mio::mesh_cad_surface(b, o, mio::MeshCadSizeSources{}, extras).mNumEdges, 16);
    b.mExtra = R"({"tolerances":{"model_tolerance":"big"}})";
    const auto bad = mio::mesh_cad_surface(b, o, mio::MeshCadSizeSources{}, extras);
    EXPECT_EQ(bad.mNumEdges, 16);
    EXPECT_TRUE(std::any_of(bad.mWarnings.begin(), bad.mWarnings.end(), [](const std::string& rW) {
        return rW.find("model's declared tolerance") != std::string::npos;
    }));
}

// Roadmap 6.5.1: the box with its top face translated by 0.07, a gap past
// what the global tolerance accepts at h = 0.2 (a quarter of h). The top's
// edges declare 0.01, which accepts their projections, but the vertical
// edges declare 1e-6, so their ends fail where they meet the top's corners.
// A vertex declaring 0.01 there accepts the corner for every edge in it.
TEST(MeshCadTolerances, AVertexsToleranceAcceptsItsCorner) {
    mio::BrepInfo b = brep_solids::box();
    mio::NDArray& r_cp = b.mPatches[static_cast<std::size_t>(b.mFaces[0].mPatch)].mControlPoints;
    for (std::size_t i = 0; i < r_cp.Size(); ++i)
        r_cp.As<double>()[i] += 0.07 / std::sqrt(3.0);
    std::vector<bool> top(b.mEdges.size(), false);
    for (std::size_t e = 0; e < b.mEdges.size(); ++e)
        for (const auto& r_r : b.mEdges[e].mTopology)
            if (r_r.mId == b.mFaces[0].mId)
                top[e] = true;
    for (std::size_t e = 0; e < b.mEdges.size(); ++e)
        b.mEdges[e].mExtra = top[e] ? R"({"tolerance":0.01})" : R"({"tolerance":1e-6})";
    const MeshCadSurfaceOptions o = quality(0.2, 0.0);
    mio::MeshCadSurfaceExtras extras;
    extras.mModelTolerances = true;
    const auto edges_only = mio::mesh_cad_surface(b, o, mio::MeshCadSizeSources{}, extras);
    EXPECT_GT(edges_only.mNumFailed, 0);
    // A vertex at every trim's start.
    for (const auto& r_f : b.mFaces)
        for (const auto& r_l : r_f.mLoops)
            for (const auto& r_t : r_l.mTrims) {
                mio::BrepVertex v;
                mio::BrepRef ref;
                ref.mId = r_f.mId;
                ref.mTrimIndex = r_t.mTrimIndex;
                v.mTopology.push_back(ref);
                v.mExtra = R"({"tolerance":0.01})";
                b.mVertices.push_back(std::move(v));
            }
    mio::MeshCadSurfaceExtrasReport x;
    const auto r = mio::mesh_cad_surface(b, o, mio::MeshCadSizeSources{}, extras, nullptr, &x);
    EXPECT_EQ(r.mNumFailed, 0);
    EXPECT_TRUE(std::any_of(x.mToleranceUses.begin(), x.mToleranceUses.end(), [&](const auto& rU) {
        return rU.mWhere == "corner" && rU.mEdge >= 0 && !top[static_cast<std::size_t>(rU.mEdge)] &&
               rU.mDeclared == 0.01;
    }));
    // Without the option the vertices change nothing.
    EXPECT_GT(run(b, o).mNumFailed, 0);
}

namespace {

// A frustum of the cone of half-angle 45 degrees about the z axis, its apex at
// the origin, from z = Z0 to z = Z1: one band closed without a seam, its two
// rims circles of radius Z0 and Z1. In the cone's own parameters (angle,
// height) the metric stretches Z1 / Z0 times from rim to rim.
mio::BrepInfo cone_band(double Z0, double Z1) {
    std::vector<mio::BrepCurve> curves = {
        mio::detail::brep_curve_arc({0, 0, Z0}, {1, 0, 0}, {0, 1, 0}, Z0, 0.0, 2 * kPi),
        mio::detail::brep_curve_arc({0, 0, Z1}, {1, 0, 0}, {0, 1, 0}, Z1, 0.0, 2 * kPi)};
    std::vector<brep_solids::SolidFace> faces(1);
    faces[0].mSurface.mKind = 4;
    faces[0].mSurface.mTanHalfAngle = 1.0;
    faces[0].mLoops = {{brep_solids::circle_trim(Z0, Z0, 2 * kPi, 0, true, false)},
                       {brep_solids::circle_trim(Z1, Z1, -2 * kPi, 1, false, false)}};
    return brep_solids::assemble(faces, curves);
}

}  // namespace

// Roadmap 6.4.4: a band closed without a seam (a cone's frustum, its
// rims full circles) is not declined: its coarse mesh keeps the seam's two
// copies of each point apart, as in the parameters, so it is a disc, and the
// conformal map flattens it. In its own parameters, which stretch tenfold
// from rim to rim, triangles stay at 4 degrees; flattened, it meets 25.
TEST(MeshCadReparam, ASeamlessBandIsFlattenedAsADisc) {
    const mio::BrepInfo band = cone_band(0.1, 1.0);
    for (const double h : {0.1, 0.05}) {
        MeshCadSurfaceOptions o = quality(h, 25.0);
        const auto plain = run(band, o);
        o.mReparamConditionNumber = 4.0;
        mio::MeshCadSurfaceExtrasReport x;
        const auto r = run(band, o, &x);
        report("cone band", r, &x);
        EXPECT_EQ(r.mNumFailed, 0);
        EXPECT_EQ(r.mNumNonManifoldEdges, 0);
        EXPECT_EQ(r.mNumMisorientedEdges, 0);
        ASSERT_EQ(x.mReparametrised.size(), 1u);
        EXPECT_TRUE(x.mReparametrised[0].mTaken) << x.mReparametrised[0].mReason;
        EXPECT_EQ(x.mReparametrised[0].mMethod, 1);
        EXPECT_GE(angle_and_area(r.mMesh)[0], 25.0 - 1e-9) << h;
        EXPECT_LT(angle_and_area(plain.mMesh)[0], 5.0) << h;
        EXPECT_LT(r.mNumPoints, plain.mNumPoints) << h;
        // The frustum's area, pi sqrt 2 (Z1^2 - Z0^2), within the chordal error.
        EXPECT_NEAR(angle_and_area(r.mMesh)[1], kPi * std::sqrt(2.0) * (1.0 - 0.01), 4.0 * h * h);
    }
    // Curved elements on the flattened band (roadmap 6.4.4): the nodes sit at
    // the flattening's coordinates, lifted through the coarse mesh, so every
    // node lies on the cone (x^2 + y^2 = z^2), and every element is proved
    // valid: none is straightened.
    for (const int order : {2, 3}) {
        MeshCadSurfaceOptions o = quality(0.1, 25.0);
        o.mReparamConditionNumber = 4.0;
        o.mOrder = order;
        const auto r = run(band, o);
        EXPECT_EQ(r.mNumFailed, 0);
        EXPECT_EQ(r.mMesh.Cells(0).Type(), order == 2 ? "triangle6" : "triangle10");
        EXPECT_EQ(r.mNumCurvedStraightened, 0) << order;
        const double* x = r.mMesh.Points().As<double>();
        for (std::size_t p = 0; p < r.mMesh.NumPoints(); ++p)
            EXPECT_NEAR(std::hypot(x[3 * p], x[3 * p + 1]), x[3 * p + 2], 1e-9) << p;
    }
}

// Roadmap 6.4.4: a coarse mesh that is an annulus (two loops, Euler
// characteristic 0, as a face with a hole has) laid between two concentric
// circles: the frustum of the cone above welded at its seam, 32 sectors by
// 6 rings. Nothing flips, each rim lies on its circle by arc length, the
// inner at exp(-2 pi w / L) of the outer, and the map does not twist. A disc
// is refused.
TEST(MeshCadReparam, AnAnnulusIsLaidBetweenTwoCircles) {
    const std::size_t ns = 32, nr = 6;
    std::vector<mio::detail::PatchVec> pts;
    for (std::size_t j = 0; j <= nr; ++j) {
        const double z = 0.1 + 0.9 * double(j) / double(nr);
        for (std::size_t i = 0; i < ns; ++i) {
            const double a = 2.0 * kPi * double(i) / double(ns);
            pts.push_back({z * std::cos(a), z * std::sin(a), z});
        }
    }
    std::vector<mio::detail::PatchTri> tris;
    for (std::size_t j = 0; j < nr; ++j)
        for (std::size_t i = 0; i < ns; ++i) {
            const std::size_t a = j * ns + i, b = j * ns + (i + 1) % ns;
            const std::size_t c = a + ns, d = b + ns;
            tris.push_back({a, b, d});
            tris.push_back({a, d, c});
        }
    const mio::detail::DiscParam p = mio::detail::param_annulus(pts, tris);
    ASSERT_TRUE(p.mOk);
    EXPECT_EQ(p.mMethod, 6);
    EXPECT_EQ(p.mNumFlips, 0u);
    EXPECT_EQ(mio::detail::param_flips(tris, p.mU, p.mV), 0u);
    const auto radius = [&](std::size_t V) { return std::hypot(p.mU[V], p.mV[V]); };
    // The longer rim (z = 1) outside, the shorter inside, each on a circle.
    const double ro = radius(nr * ns), ri = radius(0);
    for (std::size_t i = 0; i < ns; ++i) {
        EXPECT_NEAR(radius(nr * ns + i), ro, 1e-12 * ro);
        EXPECT_NEAR(radius(i), ri, 1e-12 * ro);
    }
    // The modulus: the mean rim length L, the area over L for the height.
    double lo = 0.0, li = 0.0, area = 0.0;
    for (std::size_t i = 0; i < ns; ++i) {
        const auto d = [&](std::size_t A, std::size_t B) {
            return std::hypot(pts[A][0] - pts[B][0], pts[A][1] - pts[B][1], pts[A][2] - pts[B][2]);
        };
        lo += d(nr * ns + i, nr * ns + (i + 1) % ns);
        li += d(i, (i + 1) % ns);
    }
    for (const auto& r_t : tris) {
        const auto& a = pts[r_t[0]];
        const auto& b = pts[r_t[1]];
        const auto& c = pts[r_t[2]];
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        area += 0.5 * std::hypot(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                 u[0] * v[1] - u[1] * v[0]);
    }
    const double len = 0.5 * (lo + li);
    EXPECT_NEAR(ri / ro, std::exp(-2.0 * kPi * (area / len) / len), 1e-12);
    // No twist: the inner rim's first vertex (below the outer's first) at the
    // same angle.
    EXPECT_NEAR(std::atan2(p.mV[0], p.mU[0]), std::atan2(p.mV[nr * ns], p.mU[nr * ns]), 1e-12);
    // A disc (one ring of the band cut open) is not an annulus.
    std::vector<mio::detail::PatchTri> open(tris.begin(), tris.begin() + 2 * (ns - 1));
    EXPECT_FALSE(mio::detail::param_annulus(pts, open).mOk);
}

namespace {

// Roadmap 6.5.2's chains: the plate [0, 1] x [0, 1], then `K` slivers of
// width W stacked in y (a fillet split into strips), then, with `Plate`, the
// plate [0, 1] x [1 + K W, 2 + K W]; each its own patch, sharing edges.
// Bands [0, 1] x [ys[i], ys[i + 1]] stacked in y, each its own patch,
// sharing edges.
mio::BrepInfo bands(const std::vector<double>& ys) {
    using brep_solids::line_trim;
    using brep_solids::PatchVec;
    std::vector<PatchVec> v;
    for (const double y : ys) {
        v.push_back({0, y, 0});
        v.push_back({1, y, 0});
    }
    std::vector<std::array<int, 2>> ends;
    std::vector<mio::BrepCurve> curves;
    const auto edge = [&](int A, int B) {
        const int lo = std::min(A, B), hi = std::max(A, B);
        for (std::size_t e = 0; e < ends.size(); ++e)
            if (ends[e][0] == lo && ends[e][1] == hi)
                return static_cast<std::int64_t>(e);
        ends.push_back({lo, hi});
        curves.push_back(mio::detail::brep_curve_line(v[static_cast<std::size_t>(lo)],
                                                      v[static_cast<std::size_t>(hi)]));
        return static_cast<std::int64_t>(ends.size() - 1);
    };
    std::vector<brep_solids::SolidFace> faces;
    for (std::size_t r = 0; r + 1 < ys.size(); ++r) {
        const int a = static_cast<int>(2 * r);
        const std::array<int, 4> c{a, a + 1, a + 3, a + 2};
        brep_solids::SolidFace f;
        f.mSurface =
            brep_solids::plane(v[static_cast<std::size_t>(a)], {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
        f.mLoops.emplace_back();
        for (std::size_t k = 0; k < 4; ++k) {
            const int p = c[k], q = c[(k + 1) % 4];
            f.mLoops[0].push_back(line_trim(v[static_cast<std::size_t>(p)],
                                            v[static_cast<std::size_t>(q)], edge(p, q), p < q));
        }
        faces.push_back(f);
    }
    return brep_solids::assemble(faces, curves);
}

mio::BrepInfo sliver_chain(double W, int K, bool Plate) {
    std::vector<double> ys = {0.0, 1.0};
    for (int k = 1; k <= K; ++k)
        ys.push_back(1.0 + k * W);
    if (Plate)
        ys.push_back(2.0 + K * W);
    return bands(ys);
}

// Every triangle's cad:face, and that its centroid lies in that face's band
// of y (`rYs`, the faces' bounds in order).
std::vector<std::int64_t> count_by_face(const mio::Mesh& rM, const std::vector<double>& rYs) {
    std::vector<std::int64_t> count(rYs.size() - 1, 0);
    const double* px = rM.Points().As<double>();
    const auto cells = rM.Cells(0);
    const mio::NDArray& conn = cells.Conn();
    const mio::NDArray& face = rM.CellData(mio::kMeshCadFaceName, 0);
    for (std::size_t t = 0; t < cells.NumCells(); ++t) {
        const std::int64_t f = mio::detail::read_int(face, t);
        EXPECT_GE(f, 0);
        EXPECT_LT(f, static_cast<std::int64_t>(count.size()));
        if (f < 0 || f >= static_cast<std::int64_t>(count.size()))
            continue;
        ++count[static_cast<std::size_t>(f)];
        double y = 0.0;
        for (std::size_t k = 0; k < 3; ++k)
            y += px[3 * mio::detail::read_int(conn, 3 * t + k) + 1] / 3.0;
        EXPECT_GE(y, rYs[static_cast<std::size_t>(f)] - 1e-12);
        EXPECT_LE(y, rYs[static_cast<std::size_t>(f) + 1] + 1e-12);
    }
    return count;
}

}  // namespace

// Roadmap 6.5.2: three strips of width h / 20 beside a plate. Only the first
// borders a face that is not a sliver; the second and third join the group
// it forms in the rounds that follow, so the plate's triangulation spans all
// four faces, each triangle keeping the face it lies in.
TEST(MeshCadVirtualFaces, AChainOfSliverStripsJoinsOneHost) {
    const double h = 0.2, w = h / 20.0;
    const mio::BrepInfo b = sliver_chain(w, 3, false);
    const auto alone = run(b, quality(h, 25.0));
    mio::MeshCadSurfaceExtras extras;
    extras.mVirtualFaces = true;
    mio::MeshCadSurfaceExtrasReport x;
    const auto r =
        mio::mesh_cad_surface(b, quality(h, 25.0), mio::MeshCadSizeSources{}, extras, nullptr, &x);
    EXPECT_EQ(r.mNumFailed, 0);
    EXPECT_GE(angle_and_area(r.mMesh)[0], 25.0 - 1e-9);
    EXPECT_LT(r.mNumTriangles, alone.mNumTriangles);
    EXPECT_NEAR(angle_and_area(r.mMesh)[1], 1.0 + 3.0 * w, 1e-12);
    EXPECT_EQ(r.mNumMisorientedEdges, 0);
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    EXPECT_TRUE(open_only_outside(r.mMesh, 1.0 + 3.0 * w));
    ASSERT_EQ(x.mMergedFaces.size(), 3u);
    for (std::size_t i = 0; i < 3; ++i) {
        const auto& r_m = x.mMergedFaces[i];
        EXPECT_EQ(r_m.mSliver, static_cast<std::int64_t>(i + 1));
        EXPECT_TRUE(r_m.mMerged) << r_m.mReason;
        EXPECT_EQ(r_m.mHost, 0);
        EXPECT_EQ(r_m.mEdges.size(), 1u);
    }
    const auto count = count_by_face(r.mMesh, {0.0, 1.0, 1.0 + w, 1.0 + 2.0 * w, 1.0 + 3.0 * w});
    for (std::size_t f = 0; f < count.size(); ++f) {
        EXPECT_GT(count[f], 0) << f;
        EXPECT_EQ(r.mFaces[f].mNumTriangles, count[f]) << f;
    }
    // Deterministic across thread counts.
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const std::uint64_t ref = digest(r.mMesh);
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
        EXPECT_EQ(
            digest(mio::mesh_cad_surface(b, quality(h, 25.0), mio::MeshCadSizeSources{}, extras)
                       .mMesh),
            ref)
            << threads;
    }
    omp_set_num_threads(saved);
#endif
}

// Roadmap 6.5.2: three strips with no plate, all narrower than h (the
// middle one the widest). No face but a sliver borders the chain, so its
// widest face hosts it and the other two join: one triangulation spans the
// three, each triangle keeping its own face.
TEST(MeshCadVirtualFaces, AChainOfSliversAloneIsHostedByItsWidest) {
    const double h = 0.2, w = h / 20.0;
    const std::vector<double> ys{0.0, w, 3.0 * w, 4.0 * w};
    const mio::BrepInfo b = bands(ys);
    const auto alone = run(b, quality(h, 25.0));
    mio::MeshCadSurfaceExtras extras;
    extras.mVirtualFaces = true;
    mio::MeshCadSurfaceExtrasReport x;
    const auto r =
        mio::mesh_cad_surface(b, quality(h, 25.0), mio::MeshCadSizeSources{}, extras, nullptr, &x);
    EXPECT_EQ(r.mNumFailed, 0);
    EXPECT_NEAR(angle_and_area(r.mMesh)[1], 4.0 * w, 1e-12);
    EXPECT_EQ(r.mNumMisorientedEdges, 0);
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    EXPECT_TRUE(open_only_outside(r.mMesh, 4.0 * w));
    std::printf(
        "[6.5.2] all-sliver chain: %zu triangles face by face, %zu merged; min angle "
        "%.2f against %.2f\n",
        static_cast<std::size_t>(alone.mNumTriangles), static_cast<std::size_t>(r.mNumTriangles),
        angle_and_area(alone.mMesh)[0], angle_and_area(r.mMesh)[0]);
    ASSERT_EQ(x.mMergedFaces.size(), 3u);
    for (const auto& r_m : x.mMergedFaces) {
        EXPECT_TRUE(r_m.mMerged) << r_m.mSliver << ": " << r_m.mReason;
        EXPECT_EQ(r_m.mHost, 1) << r_m.mSliver;
    }
    const auto count = count_by_face(r.mMesh, ys);
    for (std::size_t f = 0; f < count.size(); ++f)
        EXPECT_GT(count[f], 0) << f;
    // Strips of one width: the tie goes to the lower index.
    const auto uniform = mio::mesh_cad_surface(bands({0.0, w, 2.0 * w, 3.0 * w}), quality(h, 25.0),
                                               mio::MeshCadSizeSources{}, extras, nullptr, &x);
    EXPECT_EQ(uniform.mNumFailed, 0);
    for (const auto& r_m : x.mMergedFaces)
        EXPECT_EQ(r_m.mHost, 0) << r_m.mSliver;
}

// Roadmap 6.5.2: plate, three strips, plate. The outer strips merge with
// their plates in the first round; the middle one borders only those two
// groups and joins the first (the same shared length, the lower host), a
// merge across more than two faces.
TEST(MeshCadVirtualFaces, ASliverBetweenTwoMergedGroupsJoinsOne) {
    const double h = 0.2, w = h / 20.0;
    const mio::BrepInfo b = sliver_chain(w, 3, true);
    mio::MeshCadSurfaceExtras extras;
    extras.mVirtualFaces = true;
    mio::MeshCadSurfaceExtrasReport x;
    const auto r =
        mio::mesh_cad_surface(b, quality(h, 25.0), mio::MeshCadSizeSources{}, extras, nullptr, &x);
    EXPECT_EQ(r.mNumFailed, 0);
    EXPECT_GE(angle_and_area(r.mMesh)[0], 25.0 - 1e-9);
    EXPECT_NEAR(angle_and_area(r.mMesh)[1], 2.0 + 3.0 * w, 1e-12);
    EXPECT_EQ(r.mNumMisorientedEdges, 0);
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    EXPECT_TRUE(open_only_outside(r.mMesh, 2.0 + 3.0 * w));
    ASSERT_EQ(x.mMergedFaces.size(), 3u);
    const std::array<std::int64_t, 3> host{0, 0, 4};
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_TRUE(x.mMergedFaces[i].mMerged) << x.mMergedFaces[i].mReason;
        EXPECT_EQ(x.mMergedFaces[i].mHost, host[i]) << i;
    }
    const auto count =
        count_by_face(r.mMesh, {0.0, 1.0, 1.0 + w, 1.0 + 2.0 * w, 1.0 + 3.0 * w, 2.0 + 3.0 * w});
    for (std::size_t f = 0; f < count.size(); ++f)
        EXPECT_GT(count[f], 0) << f;
}
