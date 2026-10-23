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
// measured beside it, for the size fields of 6.3 to improve on.

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
#include "brep_solids.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/mesh.hpp"
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

TEST(MeshCadInterpolation, TheScordelisLoRoofConvergesAtSecondOrder) {
    // The roof (radius 25, length 50, 80 degrees): an untrimmed rational
    // patch whose parameter scales differ fortyfold, written in-test (the C++
    // suite has no fixture path) with the control net of
    // tests/python/meshes/cad_json/cylindrical_roof.cad.json.
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
    const mio::BrepInfo roof = mio::read_cad_json_info(path);
    // On the cylinder of radius 25 about the y axis.
    const Field smooth = [](const double* pX) {
        return std::sin(pX[0] / 8.0) * std::cos(pX[1] / 10.0) + std::exp(pX[2] / 40.0);
    };
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

}  // namespace
