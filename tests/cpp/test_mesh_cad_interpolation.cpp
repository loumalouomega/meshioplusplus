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
// Roadmap 6.11: a generated mesh measured by what it is for. Analytic fields
// are interpolated linearly on the meshes mesh_cad_surface makes of the
// sphere and the Scordelis-Lo roof, and the L2 and maximum interpolation
// errors are measured under halving of h, so a change to refinement shows
// its effect on accuracy, not only on angles (J. R. Shewchuk, "What is a good
// linear element?", 2002). The error operation (doc/error.md) recovers a
// gradient and integrates nothing, so the quadrature is this file's own: the
// symmetric 7-point rule of degree 5 on each flat triangle (Strang and Fix;
// D. A. Dunavant, IJNME 21, 1985), the field evaluated in space at its
// points. A smooth field converges at order 2 in both norms, pinned for the
// unrefined meshes and for min_angle = 25; a field with a sharp layer is
// measured beside it, for the size fields of 6.3 to improve on. The curved
// elements of 6.6 are measured the same way, the interpolant their own
// Lagrange one on the curved element and the quadrature a collapsed 8 x 8
// Gauss-Legendre rule on it: quadratic and cubic elements converge faster.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/curved_triangle.hpp"
#include "brep_solids.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/detail/sym3_eigen.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/hessian.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"

namespace {

namespace mio = meshioplusplus;
using Field = std::function<double(const double*)>;

/// The 7-point rule of degree 5: barycentric points and weights (summing to
/// one, so the triangle's area multiplies them).
struct Rule {
    std::array<std::array<double, 3>, 7> mPoints;
    std::array<double, 7> mWeights;
};

Rule seven_point_rule() {
    const double a1 = 0.059715871789769820, b1 = 0.470142064105115090;
    const double a2 = 0.797426985353087322, b2 = 0.101286507323456339;
    const double w0 = 0.225;
    const double w1 = 0.132394152788506181;
    const double w2 = 0.125939180544827153;
    Rule r;
    r.mPoints = {{{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
                  {a1, b1, b1},
                  {b1, a1, b1},
                  {b1, b1, a1},
                  {a2, b2, b2},
                  {b2, a2, b2},
                  {b2, b2, a2}}};
    r.mWeights = {w0, w1, w1, w1, w2, w2, w2};
    return r;
}

struct Errors {
    double mL2 = 0.0;
    double mMax = 0.0;
    double mMaxAway = 0.0;  ///< the maximum over triangles out of every pole's reach
    std::size_t mPoints = 0;
    double mMinAngleDeg = 0.0;  ///< the mesh's quality, reported beside its errors
    std::array<std::int64_t, mio::kMeshCadAngleBins> mHistogram{};
};

using Poles = std::vector<std::pair<std::array<double, 3>, double>>;

/// The L2 and maximum errors of the linear interpolant of `rF` on the mesh's
/// triangles, against `rF` at the quadrature points in space; the maximum
/// also away from the poles (their fans are left to graded poles, 6.4.2).
Errors interpolation_errors(const mio::Mesh& rM, const Field& rF, const Poles& rPoles) {
    const Rule rule = seven_point_rule();
    const double* x = rM.Points().As<double>();
    const auto cells = rM.Cells(0);
    const mio::NDArray& c = cells.Conn();
    double sum = 0.0;
    Errors e;
    e.mPoints = rM.NumPoints();
    for (std::size_t t = 0; t < cells.NumCells(); ++t) {
        const double* p[3];
        double f[3];
        for (std::size_t k = 0; k < 3; ++k) {
            p[k] = x + 3 * mio::detail::read_int(c, 3 * t + k);
            f[k] = rF(p[k]);
        }
        const double u[3] = {p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]};
        const double v[3] = {p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]};
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                             u[0] * v[1] - u[1] * v[0]};
        const double area = 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        bool near_pole = false;
        for (const auto& [r_x, reach] : rPoles)
            for (std::size_t k = 0; k < 3; ++k)
                near_pole = near_pole || std::hypot(p[k][0] - r_x[0], p[k][1] - r_x[1],
                                                    p[k][2] - r_x[2]) <= reach;
        for (std::size_t q = 0; q < 7; ++q) {
            const auto& l = rule.mPoints[q];
            double xq[3];
            for (int k = 0; k < 3; ++k)
                xq[k] = l[0] * p[0][k] + l[1] * p[1][k] + l[2] * p[2][k];
            const double err = rF(xq) - (l[0] * f[0] + l[1] * f[1] + l[2] * f[2]);
            sum += rule.mWeights[q] * area * err * err;
            e.mMax = std::max(e.mMax, std::fabs(err));
            if (!near_pole)
                e.mMaxAway = std::max(e.mMaxAway, std::fabs(err));
        }
    }
    e.mL2 = std::sqrt(sum);
    return e;
}

struct Study {
    std::vector<Errors> mErrors;
    std::vector<double> mL2Orders;
    std::vector<double> mMaxOrders;
    std::vector<double> mMaxAwayOrders;
};

Study study(const mio::BrepInfo& rBrep, const std::vector<double>& rSizes, double MinAngle,
            const Field& rF) {
    Study s;
    for (const double h : rSizes) {
        mio::MeshCadSurfaceOptions o;
        o.mTargetSize = h;
        o.mMinAngleDeg = MinAngle;
        o.mEdges = false;
        const auto r = mio::mesh_cad_surface(rBrep, o);
        EXPECT_EQ(r.mNumFailed, 0);
        // The poles' reach as refinement reports it (h on the unrefined
        // mesh, which reports none: the same rule is applied here).
        Poles poles;
        for (const auto& r_f : r.mFaces)
            for (const auto& r_c : r_f.mSmallAngleCorners)
                if (r_c.mPole)
                    poles.emplace_back(r_c.mX, r_c.mReach);
        if (MinAngle == 0.0) {
            const auto refined = [&]() {
                mio::MeshCadSurfaceOptions q = o;
                q.mMinAngleDeg = 20.0;
                return mio::mesh_cad_surface(rBrep, q);
            }();
            for (const auto& r_f : refined.mFaces)
                for (const auto& r_c : r_f.mSmallAngleCorners)
                    if (r_c.mPole)
                        poles.emplace_back(r_c.mX, r_c.mReach);
        }
        s.mErrors.push_back(interpolation_errors(r.mMesh, rF, poles));
        s.mErrors.back().mMinAngleDeg = r.mMinAngleDeg;
        s.mErrors.back().mHistogram = r.mAngleHistogram;
    }
    for (std::size_t k = 1; k < s.mErrors.size(); ++k) {
        s.mL2Orders.push_back(std::log2(s.mErrors[k - 1].mL2 / s.mErrors[k].mL2));
        s.mMaxOrders.push_back(std::log2(s.mErrors[k - 1].mMax / s.mErrors[k].mMax));
        s.mMaxAwayOrders.push_back(std::log2(s.mErrors[k - 1].mMaxAway / s.mErrors[k].mMaxAway));
    }
    return s;
}

void print(const char* pName, double MinAngle, const std::vector<double>& rSizes, const Study& rS) {
    for (std::size_t k = 0; k < rS.mErrors.size(); ++k) {
        // The angle histogram beside the errors: its 10-degree bins from 0
        // to 90, every angle of every triangle on the surface.
        std::string bins;
        for (std::size_t b = 0; b < 9; ++b)
            bins += (b ? " " : "") + std::to_string(rS.mErrors[k].mHistogram[b]);
        std::printf(
            "[6.11] %-28s min_angle %4.1f  h %-6g smallest angle %.2f  angles per 10 "
            "degrees to 90: %s\n",
            pName, MinAngle, rSizes[k], rS.mErrors[k].mMinAngleDeg, bins.c_str());
        std::printf(
            "[6.11] %-28s min_angle %4.1f  h %-6g points %-7zu L2 %.4e  max %.4e  "
            "max away from poles %.4e%s\n",
            pName, MinAngle, rSizes[k], rS.mErrors[k].mPoints, rS.mErrors[k].mL2,
            rS.mErrors[k].mMax, rS.mErrors[k].mMaxAway,
            k == 0 ? ""
                   : ("  orders " + std::to_string(rS.mL2Orders[k - 1]) + " " +
                      std::to_string(rS.mMaxOrders[k - 1]) + " " +
                      std::to_string(rS.mMaxAwayOrders[k - 1]))
                         .c_str());
    }
}

// Smooth on the unit sphere: a mix of every coordinate, nothing aligned with
// the parameters.
double sphere_smooth(const double* pX) {
    return std::sin(2.0 * pX[0]) * std::cos(1.5 * pX[1]) + std::exp(0.5 * pX[2]);
}

// A layer of width 0.05 across the sphere at z = 0.3.
double sphere_layer(const double* pX) {
    return std::tanh((pX[2] - 0.3) / 0.05);
}

TEST(MeshCadInterpolation, TheSphereConvergesAtSecondOrder) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    const std::vector<double> sizes = {0.2, 0.1, 0.05};
    for (const double ang : {0.0, 25.0}) {
        const Study s = study(sphere, sizes, ang, sphere_smooth);
        print("sphere, smooth field", ang, sizes, s);
        for (std::size_t k = 0; k < s.mL2Orders.size(); ++k) {
            EXPECT_NEAR(s.mL2Orders[k], 2.0, 0.2) << "L2, min_angle " << ang;
            // The maximum is the poles' fans' until graded poles land: it is
            // pinned away from them.
            EXPECT_NEAR(s.mMaxAwayOrders[k], 2.0, 0.35) << "max, min_angle " << ang;
        }
    }
    // The layer: measured, not pinned to an order (a size field, 6.3, should
    // reach the same error with fewer points).
    const Study layer = study(sphere, sizes, 25.0, sphere_layer);
    print("sphere, layer at z = 0.3", 25.0, sizes, layer);
    EXPECT_LT(layer.mErrors.back().mL2, layer.mErrors.front().mL2);
}

// Roadmap 6.11's layer: a size field from the field's own Hessian, carried
// as a background (6.3.1), reaches the uniform mesh's L2 error with fewer
// points. The background is a coarse mesh (h = 0.1) of the sphere with the
// layer at its points; `hessian` (composed Green-Gauss, point values) gives
// its second derivatives, projected on the tangent plane (P H P, P = I - n
// n^T), and the size at a point is sqrt(8 eps / |lambda|max), the length at
// which linear interpolation's error is about eps there (Loseille and
// Alauzet 2011), between h / 4 and the target 0.2, with eps = 8e-3. A size
// only refines, so the mesh starts at 0.2 and is refined to the layer. Both
// meshes grade their poles (two rings), so the angles meet 25 degrees.
// Measured: the uniform mesh at h = 0.05 has 11,281 points and an L2 error
// of 1.58e-2; the sized one 3,848 points and 1.39e-2.
TEST(MeshCadInterpolation, AHessianSizedMeshReachesTheLayerWithFewerPoints) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    const auto layer = [](const double* pX) { return sphere_layer(pX); };
    mio::MeshCadSurfaceOptions o;
    o.mTargetSize = 0.05;
    o.mMinAngleDeg = 25.0;
    o.mEdges = false;
    o.mPoleRings = 2;
    const auto uniform = mio::mesh_cad_surface(sphere, o);
    const Errors eu = interpolation_errors(uniform.mMesh, layer, {});
    // The background and its sizes.
    o.mTargetSize = 0.1;
    mio::Mesh bg = mio::mesh_cad_surface(sphere, o).mMesh;
    const std::size_t n = bg.NumPoints();
    {
        mio::NDArray f(mio::DType::Float64, {n});
        const double* x = bg.Points().As<double>();
        for (std::size_t i = 0; i < n; ++i)
            f.As<double>()[i] = sphere_layer(x + 3 * i);
        bg.AddPointData("f", std::move(f));
    }
    mio::HessianOptions ho;
    ho.mArrayName = "f";
    ho.mLocation = mio::DataLocation::Point;
    const mio::Mesh hm = mio::hessian(bg, ho).mMesh;
    const double* hess = hm.PointData("f:hessian").As<double>();
    const double* x = bg.Points().As<double>();
    const double eps = 8.0e-3;
    mio::NDArray size(mio::DType::Float64, {n});
    for (std::size_t i = 0; i < n; ++i) {
        const double* p = x + 3 * i;
        const double len = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        const double nn[3] = {p[0] / len, p[1] / len, p[2] / len};
        double pm[3][3], h[3][3], t[3][3], ph[3][3];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                pm[a][b] = (a == b ? 1.0 : 0.0) - nn[a] * nn[b];
                h[a][b] = 0.5 * (hess[9 * i + 3 * a + b] + hess[9 * i + 3 * b + a]);
            }
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                t[a][b] = 0.0;
                for (int k = 0; k < 3; ++k)
                    t[a][b] += pm[a][k] * h[k][b];
            }
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                ph[a][b] = 0.0;
                for (int k = 0; k < 3; ++k)
                    ph[a][b] += t[a][k] * pm[k][b];
            }
        const double six[6] = {ph[0][0], ph[1][1], ph[2][2], ph[0][1], ph[1][2], ph[2][0]};
        double eig[3];
        mio::detail::sym3_principal(six, eig);
        const double lmax = std::max(std::fabs(eig[0]), std::fabs(eig[2]));
        const double hs = lmax > 0.0 ? std::sqrt(8.0 * eps / lmax) : 0.2;
        size.As<double>()[i] = std::clamp(hs, 0.05 / 4.0, 0.2);
    }
    {
        double lo = 1e300, hi = 0;
        for (std::size_t i = 0; i < n; ++i) {
            lo = std::min(lo, size.As<double>()[i]);
            hi = std::max(hi, size.As<double>()[i]);
        }
        std::printf("[6.11] background sizes %.4g .. %.4g over %zu points\n", lo, hi, n);
    }
    bg.AddPointData("size", std::move(size));
    o.mTargetSize = 0.2;
    mio::MeshCadSizeSources src;
    src.pBackground = &bg;
    const auto sized = mio::mesh_cad_surface(sphere, o, src);
    ASSERT_EQ(sized.mNumFailed, 0);
    const Errors es = interpolation_errors(sized.mMesh, layer, {});
    std::printf(
        "[6.11] layer: uniform h 0.05 points %zu L2 %.4e max %.4e; Hessian-sized "
        "points %zu L2 %.4e max %.4e\n",
        eu.mPoints, eu.mL2, eu.mMax, es.mPoints, es.mL2, es.mMax);
    EXPECT_LE(es.mL2, eu.mL2);
    EXPECT_LT(double(es.mPoints), 0.4 * double(eu.mPoints));
    EXPECT_GE(sized.mMinAngleDeg, 25.0 - 1e-9);
}

// The roof (radius 25, length 50, 80 degrees): an untrimmed rational patch
// whose parameter scales differ fortyfold, written in-test (the C++ suite has
// no fixture path) with the control net of
// tests/python/meshes/cad_json/cylindrical_roof.cad.json.
mio::BrepInfo scordelis_lo_roof() {
    const std::string path = mt::temp_path(".cad.json");
    {
        const double c = std::cos(40.0 * 3.14159265358979323846 / 180.0);
        const double s = std::sin(40.0 * 3.14159265358979323846 / 180.0);
        const double r = 25.0;
        // A rational quadratic arc from -40 to 40 degrees in the x-z plane,
        // extruded 50 along y.
        char json[2048];
        std::snprintf(
            json, sizeof(json),
            R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "swapped_surface_normal": false,
            "surface": {"is_trimmed": false, "is_rational": true, "degrees": [2, 1],
            "knot_vectors": [[0, 0, 1, 1], [0, 50]],
            "control_points": [[1, [%.17g, 0, %.17g, 1]], [2, [0, 0, %.17g, %.17g]],
                               [3, [%.17g, 0, %.17g, 1]], [4, [%.17g, 50, %.17g, 1]],
                               [5, [0, 50, %.17g, %.17g]], [6, [%.17g, 50, %.17g, 1]]]},
            "boundary_loops": [], "embedded_edges": []}], "edges": [], "vertices": []}]})",
            -r * s, r * c, r / c, c, r * s, r * c, -r * s, r * c, r / c, c, r * s, r * c);
        std::ofstream(path) << json;
    }
    return mio::read_cad_json_info(path);
}

// Smooth on the roof, the cylinder of radius 25 about the y axis.
double roof_smooth(const double* pX) {
    return std::sin(pX[0] / 8.0) * std::cos(pX[1] / 10.0) + std::exp(pX[2] / 40.0);
}

TEST(MeshCadInterpolation, TheScordelisLoRoofConvergesAtSecondOrder) {
    const mio::BrepInfo roof = scordelis_lo_roof();
    const Field smooth = roof_smooth;
    const std::vector<double> sizes = {4.0, 2.0, 1.0};
    for (const double ang : {0.0, 25.0}) {
        const Study s = study(roof, sizes, ang, smooth);
        print("Scordelis-Lo roof, smooth", ang, sizes, s);
        for (std::size_t k = 0; k < s.mL2Orders.size(); ++k) {
            EXPECT_NEAR(s.mL2Orders[k], 2.0, 0.2) << "L2, min_angle " << ang;
            EXPECT_NEAR(s.mMaxOrders[k], 2.0, 0.35) << "max, min_angle " << ang;
        }
    }
}

/// The L2 and maximum errors of the Lagrange interpolant of `rF` on the
/// curved elements of order `Order` of the mesh's first block, against `rF`
/// at the points of a collapsed 8 x 8 Gauss-Legendre rule on each element
/// (the area element's components are J against the three axes); then both
/// again over the elements out of every pole's reach.
std::array<double, 4> curved_interpolation_errors(const mio::Mesh& rM, int Order, const Field& rF,
                                                  const Poles& rPoles) {
    namespace mdt = mio::detail;
    static const double kX[8] = {0.0198550717512319, 0.1016667612931866, 0.2372337950418355,
                                 0.4082826787521751, 0.5917173212478249, 0.7627662049581645,
                                 0.8983332387068134, 0.9801449282487681};
    static const double kW[8] = {0.0506142681451881, 0.1111905172266872, 0.1568533229389436,
                                 0.1813418916891810, 0.1813418916891810, 0.1568533229389436,
                                 0.1111905172266872, 0.0506142681451881};
    const auto cells = rM.Cells(0);
    const std::size_t nn = static_cast<std::size_t>(mdt::curved_triangle_num_nodes(Order));
    const double* x = rM.Points().As<double>();
    double sum = 0.0, max = 0.0, sum_away = 0.0, max_away = 0.0;
    std::vector<mdt::CurvedVec3> nodes(nn), values(nn);
    for (std::size_t c = 0; c < cells.NumCells(); ++c) {
        bool near_pole = false;
        for (std::size_t i = 0; i < nn; ++i) {
            const auto id = static_cast<std::size_t>(mdt::read_int(cells.Conn(), nn * c + i));
            nodes[i] = {x[3 * id], x[3 * id + 1], x[3 * id + 2]};
            // The interpolant is linear in its nodal values: the element's
            // own map with the values as first coordinates.
            values[i] = {rF(nodes[i].data()), 0.0, 0.0};
            for (const auto& [r_x, reach] : rPoles)
                near_pole = near_pole || std::hypot(nodes[i][0] - r_x[0], nodes[i][1] - r_x[1],
                                                    nodes[i][2] - r_x[2]) <= reach;
        }
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j) {
                const double xi = kX[i], eta = kX[j] * (1.0 - kX[i]);
                const double w = kW[i] * kW[j] * (1.0 - kX[i]);
                mdt::CurvedVec3 n{};
                for (std::size_t k = 0; k < 3; ++k) {
                    mdt::CurvedVec3 axis{0.0, 0.0, 0.0};
                    axis[k] = 1.0;
                    n[k] = mdt::curved_triangle_jacobian(Order, nodes.data(), axis, xi, eta);
                }
                const auto xq = mdt::curved_triangle_point(Order, nodes.data(), xi, eta);
                const double err =
                    rF(xq.data()) - mdt::curved_triangle_point(Order, values.data(), xi, eta)[0];
                const double e2 =
                    w * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * err * err;
                sum += e2;
                max = std::max(max, std::fabs(err));
                if (!near_pole) {
                    sum_away += e2;
                    max_away = std::max(max_away, std::fabs(err));
                }
            }
    }
    return {std::sqrt(sum), max, std::sqrt(sum_away), max_away};
}

TEST(MeshCadInterpolation, CurvedElementsConvergeFaster) {
    // Roadmap 6.6 against 6.11: on the sphere and the roof, the L2 error of
    // quadratic and cubic interpolation under halving of h, unrefined.
    // Measured: on the roof 2.0, 3.1 and 4.1 to 3.9 for orders 1 to 3; on
    // the sphere away from the poles 2.0, 3.0 and 3.7 to 3.8, and over all of
    // it 2.0, 2.8 to 2.5 and 3.7 to 3.4, the poles' fans (left to graded
    // poles, 6.4.2) holding the curved elements back. Pinned with margin.
    struct Case {
        const char* mName;
        mio::BrepInfo mBrep;
        Field mField;
        std::vector<double> mSizes;
        std::array<double, 3> mMinOrderAll;   ///< per order, the L2 over the mesh
        std::array<double, 3> mMinOrderAway;  ///< per order, away from the poles
    };
    const Case cases[] = {{"sphere, smooth field",
                           brep_solids::sphere(1.0),
                           sphere_smooth,
                           {0.2, 0.1, 0.05},
                           {1.8, 2.3, 3.2},
                           {1.8, 2.8, 3.5}},
                          {"Scordelis-Lo roof, smooth",
                           scordelis_lo_roof(),
                           roof_smooth,
                           {4.0, 2.0, 1.0},
                           {1.8, 2.8, 3.6},
                           {1.8, 2.8, 3.6}}};
    for (const Case& r_c : cases) {
        std::vector<std::array<double, 4>> previous;
        for (const int order : {1, 2, 3}) {
            std::vector<std::array<double, 4>> errors;
            for (const double h : r_c.mSizes) {
                mio::MeshCadSurfaceOptions o;
                o.mTargetSize = h;
                o.mOrder = order;
                o.mEdges = false;
                const auto r = mio::mesh_cad_surface(r_c.mBrep, o);
                ASSERT_EQ(r.mNumFailed, 0);
                // The poles' reach as refinement reports it, as in study().
                mio::MeshCadSurfaceOptions q = o;
                q.mMinAngleDeg = 20.0;
                q.mOrder = 1;
                Poles poles;
                for (const auto& r_f : mio::mesh_cad_surface(r_c.mBrep, q).mFaces)
                    for (const auto& r_p : r_f.mSmallAngleCorners)
                        if (r_p.mPole)
                            poles.emplace_back(r_p.mX, r_p.mReach);
                errors.push_back(curved_interpolation_errors(r.mMesh, order, r_c.mField, poles));
                const std::size_t k = errors.size() - 1;
                std::printf(
                    "[6.11] %-28s order %d  h %-6g points %-7zu L2 %.4e  max %.4e  "
                    "away from poles L2 %.4e max %.4e",
                    r_c.mName, order, h, r.mMesh.NumPoints(), errors[k][0], errors[k][1],
                    errors[k][2], errors[k][3]);
                if (k > 0)
                    std::printf("  orders %.2f %.2f %.2f %.2f",
                                std::log2(errors[k - 1][0] / errors[k][0]),
                                std::log2(errors[k - 1][1] / errors[k][1]),
                                std::log2(errors[k - 1][2] / errors[k][2]),
                                std::log2(errors[k - 1][3] / errors[k][3]));
                std::printf("\n");
                if (k > 0) {
                    const auto index = static_cast<std::size_t>(order - 1);
                    EXPECT_GT(std::log2(errors[k - 1][0] / errors[k][0]), r_c.mMinOrderAll[index])
                        << r_c.mName << ", order " << order << ", h " << h;
                    EXPECT_GT(std::log2(errors[k - 1][2] / errors[k][2]), r_c.mMinOrderAway[index])
                        << r_c.mName << ", order " << order << ", h " << h;
                }
                // A higher order is more accurate at every size.
                if (!previous.empty())
                    EXPECT_LT(errors[k][0], previous[k][0])
                        << r_c.mName << ", order " << order << ", h " << h;
            }
            previous = errors;
        }
    }
}

TEST(MeshCadInterpolation, PoleRingsLiftTheWholeSphere) {
    // Roadmap 6.4.2 against 6.11: with two rings at each pole, measured over
    // every element. Measured: the L2 orders 2.08 and 2.12, 2.74 and 3.11,
    // 3.05 and 4.41 for orders 1 to 3 (without rings 2.01 and 2.03, 2.76 and
    // 2.53, 3.73 and 3.43), and the linear maximum's 1.93 and 1.87 (without,
    // 2.07 and 1.34). The curved maximum stays at the pole's own fan (1.15 at
    // order 2 at the last halving), whose nodes keep the degenerate
    // parameters' placement: measured, not pinned.
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    const std::vector<double> sizes = {0.2, 0.1, 0.05};
    const double min_l2_order[3] = {1.9, 2.5, 2.8};
    for (const int order : {1, 2, 3}) {
        std::vector<std::array<double, 4>> errors;
        for (const double h : sizes) {
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = h;
            o.mOrder = order;
            o.mPoleRings = 2;
            o.mEdges = false;
            const auto r = mio::mesh_cad_surface(sphere, o);
            ASSERT_EQ(r.mNumFailed, 0);
            errors.push_back(curved_interpolation_errors(r.mMesh, order, sphere_smooth, {}));
            const std::size_t k = errors.size() - 1;
            std::printf("[6.11] %-28s order %d  h %-6g points %-7zu L2 %.4e  max %.4e",
                        "sphere, pole rings 2", order, h, r.mMesh.NumPoints(), errors[k][0],
                        errors[k][1]);
            if (k > 0)
                std::printf("  orders %.2f %.2f", std::log2(errors[k - 1][0] / errors[k][0]),
                            std::log2(errors[k - 1][1] / errors[k][1]));
            std::printf("\n");
            if (k > 0) {
                EXPECT_GT(std::log2(errors[k - 1][0] / errors[k][0]),
                          min_l2_order[static_cast<std::size_t>(order - 1)])
                    << "order " << order << ", h " << h;
                if (order == 1)
                    EXPECT_GT(std::log2(errors[k - 1][1] / errors[k][1]), 1.6) << "h " << h;
            }
        }
    }
}

}  // namespace

namespace {

// The aspect ratio of a triangle: its longest edge over the height on it.
double aspect(const double* pA, const double* pB, const double* pC) {
    const double* p[3] = {pA, pB, pC};
    double longest = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double* a = p[k];
        const double* b = p[(k + 1) % 3];
        longest = std::max(longest,
                           std::sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]) +
                                     (b[2] - a[2]) * (b[2] - a[2])));
    }
    const double u[3] = {pB[0] - pA[0], pB[1] - pA[1], pB[2] - pA[2]};
    const double v[3] = {pC[0] - pA[0], pC[1] - pA[1], pC[2] - pA[2]};
    const double c[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                         u[0] * v[1] - u[1] * v[0]};
    const double area2 = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    return area2 > 0.0 ? longest * longest / area2 : 1e300;
}

}  // namespace

// Roadmap 6.10.2: the layer of 6.11, sized by the metric of the same
// Hessian (R |Λ| R^T / (8 eps), Loseille and Alauzet 2011) rather than by
// its largest eigenvalue alone. The triangles stretch along the layer
// (where the field hardly varies), so the mesh reaches the uniform mesh's
// error with a fraction of its points. Measured: 1,008 points, L2 1.53e-2
// and max 5.5e-2, a mean aspect ratio of 5.4 in the layer, against the
// uniform mesh's 11,281 points (L2 1.58e-2, max 9.7e-2) and the isotropic
// Hessian size's 3,848 (L2 1.39e-2, max 1.4e-1). Bad arguments are refused
// by name.
TEST(MeshCadInterpolation, AMetricFromTheHessianStretchesAlongTheLayer) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    const auto layer = [](const double* pX) { return sphere_layer(pX); };
    mio::MeshCadSurfaceOptions o;
    o.mTargetSize = 0.05;
    o.mMinAngleDeg = 25.0;
    o.mEdges = false;
    o.mPoleRings = 2;
    const auto uniform = mio::mesh_cad_surface(sphere, o);
    const Errors eu = interpolation_errors(uniform.mMesh, layer, {});
    o.mTargetSize = 0.1;
    mio::Mesh bg = mio::mesh_cad_surface(sphere, o).mMesh;
    {
        const std::size_t n = bg.NumPoints();
        mio::NDArray f(mio::DType::Float64, {n});
        const double* x = bg.Points().As<double>();
        for (std::size_t i = 0; i < n; ++i)
            f.As<double>()[i] = sphere_layer(x + 3 * i);
        bg.AddPointData("f", std::move(f));
    }
    mio::HessianOptions ho;
    ho.mArrayName = "f";
    ho.mLocation = mio::DataLocation::Point;
    const mio::Mesh hm = mio::hessian(bg, ho).mMesh;
    o.mTargetSize = 0.2;
    mio::MeshCadSizeSources src;
    src.pMetricBackground = &hm;
    src.mMetricArray = "f:hessian";
    src.mMetricError = 8.0e-3;
    const auto r = mio::mesh_cad_surface(sphere, o, src);
    ASSERT_EQ(r.mNumFailed, 0);
    EXPECT_EQ(r.mNumOpenEdges, 0);
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    EXPECT_EQ(r.mNumMisorientedEdges, 0);
    const Errors em = interpolation_errors(r.mMesh, layer, {});
    // The triangles within the layer (|z - 0.3| < 0.05): stretched.
    double sum = 0.0;
    std::size_t in_layer = 0;
    {
        const double* px = r.mMesh.Points().As<double>();
        const auto cells = r.mMesh.Cells(0);
        const mio::NDArray& c = cells.Conn();
        for (std::size_t t = 0; t < cells.NumCells(); ++t) {
            const double* p[3];
            for (std::size_t k = 0; k < 3; ++k)
                p[k] = px + 3 * mio::detail::read_int(c, 3 * t + k);
            const double z = (p[0][2] + p[1][2] + p[2][2]) / 3.0;
            if (std::fabs(z - 0.3) < 0.05) {
                sum += aspect(p[0], p[1], p[2]);
                ++in_layer;
            }
        }
    }
    const double mean_aspect = in_layer ? sum / double(in_layer) : 0.0;
    std::printf(
        "[6.10.2] layer: uniform h 0.05 points %zu L2 %.4e max %.4e; metric points %zu "
        "L2 %.4e max %.4e, mean aspect in the layer %.3g over %zu triangles\n",
        eu.mPoints, eu.mL2, eu.mMax, em.mPoints, em.mL2, em.mMax, mean_aspect, in_layer);
    EXPECT_LE(em.mL2, eu.mL2);
    EXPECT_LE(em.mMax, eu.mMax);
    EXPECT_LT(double(em.mPoints), 0.15 * double(eu.mPoints));
    EXPECT_GT(mean_aspect, 3.0);
    // Deterministic.
    EXPECT_EQ(mio::mesh_cad_surface(sphere, o, src).mNumPoints, r.mNumPoints);
    // Refused by name.
    mio::MeshCadSizeSources bad = src;
    bad.mMetricError = 0.0;
    EXPECT_THROW((void)mio::mesh_cad_surface(sphere, o, bad), std::invalid_argument);
    bad = src;
    bad.mMetricArray = "missing";
    EXPECT_THROW((void)mio::mesh_cad_surface(sphere, o, bad), std::invalid_argument);
    bad = mio::MeshCadSizeSources{};
    bad.mMetricError = 1e-3;
    EXPECT_THROW((void)mio::mesh_cad_surface(sphere, o, bad), std::invalid_argument);
    // The field itself as the source: its Hessian at the points is taken
    // as hessian takes it, so the mesh is the same, byte for byte.
    mio::MeshCadSizeSources from_field = src;
    from_field.pMetricBackground = &bg;
    from_field.mMetricArray = "f";
    mio::MeshCadSurfaceExtrasReport xf, xh;
    const auto rf =
        mio::mesh_cad_surface(sphere, o, from_field, mio::MeshCadSurfaceExtras{}, nullptr, &xf);
    (void)mio::mesh_cad_surface(sphere, o, src, mio::MeshCadSurfaceExtras{}, nullptr, &xh);
    EXPECT_EQ(rf.mNumPoints, r.mNumPoints);
    EXPECT_EQ(rf.mNumTriangles, r.mNumTriangles);
    // The angles in the metric, which refinement judges, meet the bound the
    // surface's own angles do not.
    std::printf("[6.10.2] least angle in the metric %.2f, on the surface %.2f\n",
                xh.mMinMetricAngleDeg, r.mMinAngleDeg);
    EXPECT_NEAR(xf.mMinMetricAngleDeg, xh.mMinMetricAngleDeg, 1e-12);
    EXPECT_GE(xh.mMinMetricAngleDeg, 25.0 - 1e-6);
    EXPECT_LT(r.mMinAngleDeg, 25.0);
    // Curved elements and quads take the metric (roadmap 6.10.2, 6.8.3): the
    // all-quad split asks it for twice its size.
    for (const int order : {2, 3}) {
        mio::MeshCadSurfaceOptions o2 = o;
        o2.mOrder = order;
        const auto r2 = mio::mesh_cad_surface(sphere, o2, src);
        EXPECT_EQ(r2.mNumFailed, 0) << order;
        EXPECT_EQ(r2.mNumOpenEdges, 0) << order;
    }
    mio::MeshCadSurfaceExtras all;
    all.mQuads = mio::MeshCadQuads::All;
    mio::MeshCadSurfaceExtrasReport xq;
    const auto rq = mio::mesh_cad_surface(sphere, o, src, all, nullptr, &xq);
    EXPECT_EQ(rq.mNumFailed, 0);
    EXPECT_EQ(rq.mNumOpenEdges, 0);
    EXPECT_GT(xq.mNumQuads, 0);
    EXPECT_EQ(xq.mNumLeftoverTriangles, 0);
}
