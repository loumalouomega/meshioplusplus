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
// mesh is the option off's. Order 2 is refused.
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
    EXPECT_THROW((void)mio::mesh_cad_surface(b, o, mio::MeshCadSizeSources{}, extras),
                 std::invalid_argument);
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
