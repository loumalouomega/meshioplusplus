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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
/**
 * @file test_iga.cpp
 * @brief The isogeometric layer: knot insertion and degree elevation that
 *        leave the geometry, the fields and a trimmed face's trimming loops
 *        exactly where they were, evaluation with derivatives and normals,
 *        the generic NURBS sampler and its MFEM twin, and the error paths,
 *        on the native path directly (the Python tests cover parity with the
 *        numpy twin on every NURBS fixture).
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/mfem.hpp"
#include "../../src/cpp/src/detail/edge_sampler.hpp"
#include "../../src/cpp/src/detail/surface_project.hpp"
#include "../../src/cpp/src/formats/lagrange_common.hpp"
#include "meshioplusplus/operations/iga.hpp"

namespace {

using meshioplusplus::BrepInfo;
using meshioplusplus::IgaRefineOptions;
using meshioplusplus::IgaSampleOptions;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::NurbsInfo;
using meshioplusplus::NurbsPatch;

NDArray doubles(const std::vector<std::size_t>& rShape, const std::vector<double>& rValues) {
    NDArray a(meshioplusplus::DType::Float64, rShape);
    std::copy(rValues.begin(), rValues.end(), a.As<double>());
    return a;
}

// A patch from its parts: control points are rows of `sdim`, weights optional.
NurbsPatch make_patch(std::vector<std::int32_t> Degrees, std::vector<std::vector<double>> Knots,
                      std::size_t Sdim, const std::vector<double>& rPoints,
                      const std::vector<double>& rWeights = {}) {
    NurbsPatch p;
    p.mDegrees = std::move(Degrees);
    p.mKnots = std::move(Knots);
    const std::size_t n = rPoints.size() / Sdim;
    p.mControlPoints = doubles({n, Sdim}, rPoints);
    p.mWeights = doubles({n}, rWeights.empty() ? std::vector<double>(n, 1.0) : rWeights);
    return p;
}

// The rational quadratic quarter circle of radius 1.
NurbsPatch arc() {
    const double s = std::sqrt(0.5);
    return make_patch({2}, {{0, 0, 0, 1, 1, 1}}, 2, {1, 0, 1, 1, 0, 1}, {1.0, s, 1.0});
}

// The quarter annulus, radii 1 and 2, order 2 x 1 and rational: a surface in 3-D.
NurbsPatch annulus() {
    const double s = std::sqrt(0.5);
    // rows along v (radius 1 then 2), three control points along u each
    std::vector<double> pts, w;
    for (const double r : {1.0, 2.0}) {
        pts.insert(pts.end(), {r, 0, 0, r, r, 0, 0, r, 0});
        w.insert(w.end(), {1.0, s, 1.0});
    }
    return make_patch({2, 1}, {{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}, 3, pts, w);
}

// A trilinear hexahedron with its corners pushed about, polynomial.
NurbsPatch hexahedron() {
    return make_patch(
        {1, 1, 1}, {{0, 0, 1, 1}, {0, 0, 1, 1}, {0, 0, 1, 1}}, 3,
        {0, 0, 0, 1, 0, 0, 0, 1, 0, 1.2, 1.1, 0, 0, 0, 1, 1, 0, 1.3, 0, 1, 1, 1.1, 1.2, 1.4});
}

NDArray grid(const NurbsPatch& rP, std::size_t PerDirection) {
    const std::size_t dim = rP.mDegrees.size();
    std::size_t total = 1;
    for (std::size_t d = 0; d < dim; ++d)
        total *= PerDirection;
    NDArray u(meshioplusplus::DType::Float64, {total, dim});
    double* pu = u.As<double>();
    for (std::size_t i = 0; i < total; ++i) {
        std::size_t rest = i;
        for (std::size_t d = 0; d < dim; ++d) {
            const std::size_t k = rest % PerDirection;
            rest /= PerDirection;
            const double lo = rP.mKnots[d][static_cast<std::size_t>(rP.mDegrees[d])];
            const double hi =
                rP.mKnots[d][rP.mKnots[d].size() - static_cast<std::size_t>(rP.mDegrees[d]) - 1];
            pu[i * dim + d] =
                lo + (hi - lo) * static_cast<double>(k) / static_cast<double>(PerDirection - 1);
        }
    }
    return u;
}

NurbsInfo one_patch(NurbsPatch P) {
    NurbsInfo info;
    info.mPatches.push_back(std::move(P));
    return info;
}

double max_diff(const NDArray& rA, const NDArray& rB) {
    EXPECT_EQ(rA.Size(), rB.Size());
    double d = 0.0;
    for (std::size_t i = 0; i < std::min(rA.Size(), rB.Size()); ++i)
        d = std::max(d, std::abs(rA.As<double>()[i] - rB.As<double>()[i]));
    return d;
}

TEST(Iga, KnotInsertionKeepsTheArcExact) {
    const NurbsPatch p = arc();
    const NurbsPatch q = meshioplusplus::iga_insert_knots(p, 0, {0.25, 0.5, 0.5, 0.75});
    EXPECT_EQ(q.mKnots[0], (std::vector<double>{0, 0, 0, 0.25, 0.5, 0.5, 0.75, 1, 1, 1}));
    EXPECT_EQ(q.mControlPoints.Shape()[0], 7u);
    EXPECT_EQ(q.mWeights.Size(), 7u);
    const NDArray u = grid(p, 41);
    const auto a = meshioplusplus::iga_evaluate(one_patch(p), 0, u);
    const auto b = meshioplusplus::iga_evaluate(one_patch(q), 0, u);
    EXPECT_LT(max_diff(a.mPoints, b.mPoints), 1e-14);
    for (std::size_t i = 0; i < b.mPoints.Shape()[0]; ++i)
        EXPECT_NEAR(std::hypot(b.mPoints.As<double>()[i * 2], b.mPoints.As<double>()[i * 2 + 1]),
                    1.0, 1e-14);
}

TEST(Iga, DegreeElevationKeepsTheArcExactAndRaisesContinuity) {
    const NurbsPatch p = arc();
    const NurbsPatch q = meshioplusplus::iga_elevate_degree(p, 0, 2);
    EXPECT_EQ(q.mDegrees[0], 4);
    EXPECT_EQ(q.mKnots[0], (std::vector<double>{0, 0, 0, 0, 0, 1, 1, 1, 1, 1}));
    EXPECT_EQ(q.mControlPoints.Shape()[0], 5u);  // a quartic Bezier segment
    const NDArray u = grid(p, 31);
    EXPECT_LT(max_diff(meshioplusplus::iga_evaluate(one_patch(p), 0, u).mPoints,
                       meshioplusplus::iga_evaluate(one_patch(q), 0, u).mPoints),
              1e-14);
}

TEST(Iga, ElevationOfAMultiSpanCurveKeepsEveryKnotsContinuity) {
    // C1 at 1: quadratic, one interior knot of multiplicity 1
    const NurbsPatch p = make_patch({2}, {{0, 0, 0, 1, 2, 2, 2}}, 2, {0, 0, 1, 2, 2, 0, 3, 2});
    const NurbsPatch q = meshioplusplus::iga_elevate_degree(p, 0, 2);
    EXPECT_EQ(q.mKnots[0], (std::vector<double>{0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2}));
    const NDArray u = grid(p, 51);
    EXPECT_LT(max_diff(meshioplusplus::iga_evaluate(one_patch(p), 0, u).mPoints,
                       meshioplusplus::iga_evaluate(one_patch(q), 0, u).mPoints),
              1e-14);
}

TEST(Iga, RefinesASurfaceAndItsDerivativesStay) {
    IgaRefineOptions o;
    o.mSubdivisions = {3, 2, 1};
    o.mElevate = {1, 2, 0};
    const NurbsInfo in = one_patch(annulus());
    const NurbsInfo out = meshioplusplus::iga_refine(in, o);
    const NurbsPatch& q = out.mPatches[0];
    EXPECT_EQ(q.mDegrees, (std::vector<std::int32_t>{3, 3}));
    const NDArray u = grid(in.mPatches[0], 9);
    const auto a = meshioplusplus::iga_evaluate(in, 0, u, 1);
    const auto b = meshioplusplus::iga_evaluate(out, 0, u, 1);
    EXPECT_LT(max_diff(a.mPoints, b.mPoints), 1e-13);
    EXPECT_LT(max_diff(a.mDerivatives, b.mDerivatives), 1e-12);
    EXPECT_LT(max_diff(a.mNormals, b.mNormals), 1e-12);
    // the annulus is a 2-D surface in the plane: unit normals along +-z, radius 1..2
    for (std::size_t i = 0; i < a.mNormals.Shape()[0]; ++i) {
        const double* n = a.mNormals.As<double>() + i * 3;
        EXPECT_NEAR(std::abs(n[2]), 1.0, 1e-13);
        EXPECT_NEAR(n[0], 0.0, 1e-13);
    }
    // cleared, because the sampling they described is gone
    EXPECT_EQ(out.mCellPatch.Size(), 0u);
    EXPECT_EQ(out.mBoundaryCellPatch.Size(), 0u);
}

TEST(Iga, DerivativesMatchCentralDifferencesInAVolume) {
    const NurbsInfo in = one_patch(hexahedron());
    IgaRefineOptions o;
    o.mSubdivisions = {2, 2, 2};
    o.mElevate = {1, 1, 1};
    const NurbsInfo out = meshioplusplus::iga_refine(in, o);
    const NDArray u = grid(in.mPatches[0], 5);
    const auto a = meshioplusplus::iga_evaluate(in, 0, u);
    const auto b = meshioplusplus::iga_evaluate(out, 0, u, 1);
    EXPECT_LT(max_diff(a.mPoints, b.mPoints), 1e-13);
    EXPECT_EQ(b.mNormals.Size(), 0u);  // a volume has no normal
    const std::size_t n = u.Shape()[0];
    const double h = 1e-6;
    for (std::size_t d = 0; d < 3; ++d) {
        NDArray fwd = u, bwd = u;
        for (std::size_t i = 0; i < n; ++i) {
            // stay inside the domain: central differences from inner points only
            double* f = fwd.As<double>() + i * 3 + d;
            double* r = bwd.As<double>() + i * 3 + d;
            const double base = *f;
            *f = std::min(base + h, 1.0);
            *r = std::max(base - h, 0.0);
        }
        const auto pf = meshioplusplus::iga_evaluate(out, 0, fwd);
        const auto pb = meshioplusplus::iga_evaluate(out, 0, bwd);
        for (std::size_t i = 0; i < n; ++i) {
            const double span = fwd.As<double>()[i * 3 + d] - bwd.As<double>()[i * 3 + d];
            for (std::size_t c = 0; c < 3; ++c) {
                const double fd =
                    (pf.mPoints.As<double>()[i * 3 + c] - pb.mPoints.As<double>()[i * 3 + c]) /
                    span;
                EXPECT_NEAR(b.mDerivatives.As<double>()[(i * 3 + d) * 3 + c], fd, 1e-6);
            }
        }
    }
}

// The second derivatives against central differences of the first, in every
// parametric direction, on a rational curve, a rational surface and a volume.
namespace {

void expect_second_derivatives_match_differences(const NurbsInfo& rInfo, std::size_t PerDirection,
                                                 double Tolerance) {
    const NurbsPatch& patch = rInfo.mPatches[0];
    const std::size_t dim = patch.mDegrees.size();
    const std::size_t sdim = patch.mControlPoints.Shape()[1];
    NDArray u = grid(patch, PerDirection);
    const std::size_t n = u.Shape()[0];
    // Keep every point two steps inside the domain, so the differences below
    // are central (the ends would make them one-sided).
    const double h = 1e-5;
    for (std::size_t d = 0; d < dim; ++d) {
        const double lo = patch.mKnots[d][static_cast<std::size_t>(patch.mDegrees[d])];
        const double hi =
            patch.mKnots[d][patch.mKnots[d].size() - static_cast<std::size_t>(patch.mDegrees[d]) - 1];
        for (std::size_t i = 0; i < n; ++i) {
            double& v = u.As<double>()[i * dim + d];
            v = std::min(std::max(v, lo + 2 * h), hi - 2 * h);
        }
    }
    const auto ev = meshioplusplus::iga_evaluate(rInfo, 0, u, 2);
    ASSERT_EQ(ev.mSecondDerivatives.Shape(), (std::vector<std::size_t>{n, dim, dim, sdim}));
    // Order 2 carries the points and first derivatives too.
    const auto ev1 = meshioplusplus::iga_evaluate(rInfo, 0, u, 1);
    EXPECT_LT(max_diff(ev.mPoints, ev1.mPoints), 1e-13);
    EXPECT_LT(max_diff(ev.mDerivatives, ev1.mDerivatives), 1e-12);
    for (std::size_t l = 0; l < dim; ++l) {
        const double lo = patch.mKnots[l][static_cast<std::size_t>(patch.mDegrees[l])];
        const double hi =
            patch.mKnots[l][patch.mKnots[l].size() - static_cast<std::size_t>(patch.mDegrees[l]) - 1];
        NDArray fwd = u, bwd = u;
        for (std::size_t i = 0; i < n; ++i) {
            double* f = fwd.As<double>() + i * dim + l;
            double* r = bwd.As<double>() + i * dim + l;
            const double base = *f;
            *f = std::min(base + h, hi);
            *r = std::max(base - h, lo);
        }
        const auto pf = meshioplusplus::iga_evaluate(rInfo, 0, fwd, 1);
        const auto pb = meshioplusplus::iga_evaluate(rInfo, 0, bwd, 1);
        for (std::size_t i = 0; i < n; ++i) {
            const double span = fwd.As<double>()[i * dim + l] - bwd.As<double>()[i * dim + l];
            for (std::size_t k = 0; k < dim; ++k)
                for (std::size_t c = 0; c < sdim; ++c) {
                    const double fd = (pf.mDerivatives.As<double>()[(i * dim + k) * sdim + c] -
                                       pb.mDerivatives.As<double>()[(i * dim + k) * sdim + c]) /
                                      span;
                    EXPECT_NEAR(ev.mSecondDerivatives.As<double>()[((i * dim + k) * dim + l) * sdim + c],
                                fd, Tolerance)
                        << "point " << i << " d" << k << " d" << l << " component " << c;
                }
        }
    }
    // The Hessian is symmetric.
    const double* h2 = ev.mSecondDerivatives.As<double>();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim; ++k)
            for (std::size_t l = 0; l < dim; ++l)
                for (std::size_t c = 0; c < sdim; ++c)
                    EXPECT_EQ(h2[((i * dim + k) * dim + l) * sdim + c],
                              h2[((i * dim + l) * dim + k) * sdim + c]);
}

}  // namespace

TEST(Iga, SecondDerivativesMatchDifferencesOfTheFirst) {
    expect_second_derivatives_match_differences(one_patch(arc()), 7, 1e-6);
    expect_second_derivatives_match_differences(one_patch(annulus()), 5, 1e-5);
    expect_second_derivatives_match_differences(one_patch(hexahedron()), 4, 1e-5);
}

TEST(Iga, SecondDerivativesAgreeWithTheSurfaceKernel) {
    // The general kernel against the surface-only one `BrepSurfaceEval::Eval2` runs.
    const NurbsInfo in = one_patch(annulus());
    const NurbsPatch& patch = in.mPatches[0];
    const meshioplusplus::detail::BrepSurfaceEval surface(patch, true);
    const NDArray u = grid(patch, 6);
    const auto ev = meshioplusplus::iga_evaluate(in, 0, u, 2);
    const std::size_t n = u.Shape()[0];
    for (std::size_t i = 0; i < n; ++i) {
        double x[3], su[3], sv[3], suu[3], suv[3], svv[3];
        surface.Eval2(u.As<double>()[i * 2], u.As<double>()[i * 2 + 1], x, su, sv, suu, suv, svv);
        const double* h = ev.mSecondDerivatives.As<double>() + i * 4 * 3;
        for (std::size_t c = 0; c < 3; ++c) {
            EXPECT_NEAR(h[(0 * 2 + 0) * 3 + c], suu[c], 1e-10);
            EXPECT_NEAR(h[(0 * 2 + 1) * 3 + c], suv[c], 1e-10);
            EXPECT_NEAR(h[(1 * 2 + 1) * 3 + c], svv[c], 1e-10);
        }
    }
}

TEST(Iga, SecondDerivativesOfAQuarterCircleAreCentripetal) {
    // arc() is a rational quarter circle of radius 1 whose speed is not
    // constant, so only the closed identity x . x' = 0 differentiated twice is
    // checked: x . x'' + x' . x' = 0.
    const NurbsInfo in = one_patch(arc());
    const NDArray u = grid(in.mPatches[0], 9);
    const auto ev = meshioplusplus::iga_evaluate(in, 0, u, 2);
    const std::size_t n = u.Shape()[0];
    const std::size_t sdim = ev.mPoints.Shape()[1];
    for (std::size_t i = 0; i < n; ++i) {
        double xx2 = 0.0, d1d1 = 0.0;
        for (std::size_t c = 0; c < sdim; ++c) {
            xx2 += ev.mPoints.As<double>()[i * sdim + c] * ev.mSecondDerivatives.As<double>()[i * sdim + c];
            d1d1 += ev.mDerivatives.As<double>()[i * sdim + c] * ev.mDerivatives.As<double>()[i * sdim + c];
        }
        EXPECT_NEAR(xx2 + d1d1, 0.0, 1e-10) << i;
    }
}

TEST(Iga, FieldsRefineWithTheGeometry) {
    NurbsInfo in = one_patch(arc());
    in.mFieldNames = {"u"};
    // a field with two components, one control value per control point
    in.mPatches[0].mFields.push_back(doubles({3, 2}, {1, 10, 2, 20, 3, 30}));
    IgaRefineOptions o;
    o.mSubdivisions = {4, 1, 1};
    o.mElevate = {1, 0, 0};
    const NurbsInfo out = meshioplusplus::iga_refine(in, o);
    ASSERT_EQ(out.mPatches[0].mFields.size(), 1u);
    EXPECT_EQ(out.mPatches[0].mFields[0].Shape()[0], out.mPatches[0].mWeights.Size());
    EXPECT_EQ(out.mPatches[0].mFields[0].Shape()[1], 2u);
    const NDArray u = grid(in.mPatches[0], 17);
    const auto a = meshioplusplus::iga_evaluate(in, 0, u);
    const auto b = meshioplusplus::iga_evaluate(out, 0, u);
    ASSERT_EQ(b.mFields.size(), 1u);
    EXPECT_LT(max_diff(a.mFields[0], b.mFields[0]), 1e-13);
}

TEST(Iga, APolynomialPatchKeepsUnitWeights) {
    const NurbsPatch q = meshioplusplus::iga_elevate_degree(hexahedron(), 1, 1);
    for (std::size_t i = 0; i < q.mWeights.Size(); ++i)
        EXPECT_EQ(q.mWeights.As<double>()[i], 1.0);
}

TEST(Iga, SamplesAHandBuiltModelAtItsLattice) {
    NurbsInfo in = one_patch(annulus());
    in.mFieldNames = {"f"};
    std::vector<double> f;
    for (std::size_t i = 0; i < 6; ++i)
        f.push_back(static_cast<double>(i));
    in.mPatches[0].mFields.push_back(doubles({6, 1}, f));
    IgaSampleOptions so;
    so.mSubdivisions = 2;
    NurbsInfo mapped;
    const Mesh mesh = meshioplusplus::iga_sample(in, so, &mapped);
    // one span in u, one in v, split in 2 each: 4 quad9 cells, 25 nodes
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "quad9");
    EXPECT_EQ(mesh.Cells(0).NumCells(), 4u);
    EXPECT_EQ(mesh.Points().Shape()[0], 25u);
    EXPECT_TRUE(mesh.HasPointData("f"));
    EXPECT_EQ(mapped.mCellPatch.Size(), 4u);
    EXPECT_EQ(mapped.mCellSpan.Shape()[1], 3u);
    // every node of the annulus lies between radius 1 and 2 on z = 0
    for (std::size_t i = 0; i < 25; ++i) {
        const double* x = mesh.Points().As<double>() + i * 3;
        const double r = std::hypot(x[0], x[1]);
        EXPECT_GE(r, 1.0 - 1e-13);
        EXPECT_LE(r, 2.0 + 1e-13);
        EXPECT_EQ(x[2], 0.0);
    }
    so.mSubdivisions = 65;
    EXPECT_THROW(meshioplusplus::iga_sample(in, so), std::invalid_argument);
}

// --- Bezier extraction (roadmap 9.8.2) ---------------------------------------

namespace {

double bernstein(int N, int I, double T) {
    double c = 1.0;
    for (int k = 1; k <= I; ++k)
        c = c * (N - I + k) / k;
    return c * std::pow(T, I) * std::pow(1.0 - T, N - I);
}

// The point of Bezier cell `Cell` of `rMesh` at local coordinates `pT`
// (rational: the weighted Bernstein quotient).
std::array<double, 3> bezier_point(const Mesh& rMesh, std::size_t Cell, std::size_t Dim, int P,
                                   const double* pT) {
    using meshioplusplus::lagrange::Shape;
    const auto lat = meshioplusplus::lagrange::vtk_lattice(Dim == 1   ? Shape::Line
                                                           : Dim == 2 ? Shape::Quad
                                                                      : Shape::Hexahedron,
                                                           P);
    const NDArray& pts = rMesh.Points();
    const std::size_t cols = pts.Shape()[1];
    const auto cb = *rMesh.CellRange().begin();
    const std::int64_t* conn = cb.Conn().As<std::int64_t>() + Cell * lat.size();
    const bool rational = rMesh.HasPointData("RationalWeights");
    std::array<double, 3> num{0, 0, 0};
    double den = 0.0;
    for (std::size_t n = 0; n < lat.size(); ++n) {
        double b = 1.0;
        for (std::size_t d = 0; d < Dim; ++d)
            b *= bernstein(P, lat[n][d], pT[d]);
        const double w = rational ? rMesh.PointData("RationalWeights").As<double>()[conn[n]] : 1.0;
        for (std::size_t d = 0; d < cols; ++d)
            num[d] += b * w * pts.As<double>()[static_cast<std::size_t>(conn[n]) * cols + d];
        den += b * w;
    }
    for (std::size_t d = 0; d < 3; ++d)
        num[d] /= den;
    return num;
}

// Every cell of the extraction reproduces `iga_evaluate` over its own span box.
void expect_extraction_matches(const NurbsInfo& rInfo, std::size_t ExpectedCells) {
    const NurbsPatch& patch = rInfo.mPatches[0];
    const std::size_t dim = patch.mDegrees.size();
    const int p = patch.mDegrees[0];
    const Mesh mesh = meshioplusplus::iga_bezier_extract(rInfo, 0);
    const auto cb = *mesh.CellRange().begin();
    ASSERT_EQ(cb.NumCells(), ExpectedCells);
    const NDArray& span = mesh.CellData("iga:span", 0);
    const NDArray& deg = mesh.CellData("HigherOrderDegrees", 0);
    std::array<std::vector<double>, 3> br;
    for (std::size_t d = 0; d < dim; ++d)
        for (std::size_t i = 0; i < patch.mKnots[d].size(); ++i)
            if (i == 0 || patch.mKnots[d][i] != patch.mKnots[d][i - 1])
                br[d].push_back(patch.mKnots[d][i]);
    for (std::size_t c = 0; c < ExpectedCells; ++c) {
        EXPECT_EQ(deg.As<double>()[c * 3], p);
        EXPECT_EQ(deg.As<double>()[c * 3 + 2], dim > 2 ? p : 0.0);
        double lo[3] = {0, 0, 0}, hi[3] = {1, 1, 1};
        for (std::size_t d = 0; d < dim; ++d) {
            const std::size_t k = static_cast<std::size_t>(span.As<std::int64_t>()[c * 3 + d]);
            lo[d] = br[d][k];
            hi[d] = br[d][k + 1];
        }
        for (const double t : {0.0, 0.25, 0.6, 1.0}) {
            double tt[3] = {t, 1.0 - t, 0.5 * t};
            double u[3];
            for (std::size_t d = 0; d < dim; ++d)
                u[d] = lo[d] + (hi[d] - lo[d]) * tt[d];
            const auto bz = bezier_point(mesh, c, dim, p, tt);
            NDArray up(meshioplusplus::DType::Float64, {1, dim});
            std::copy(u, u + dim, up.As<double>());
            const auto ev = meshioplusplus::iga_evaluate(rInfo, 0, up, 0);
            for (std::size_t d = 0; d < patch.mControlPoints.Shape()[1]; ++d)
                EXPECT_NEAR(bz[d], ev.mPoints.As<double>()[d], 1e-13);
        }
    }
}

}  // namespace

TEST(Iga, BezierExtractionOfAnArcIsExactAndRational) {
    expect_extraction_matches(one_patch(arc()), 1);
    IgaRefineOptions o;
    o.mSubdivisions = {3, 1, 1};
    expect_extraction_matches(meshioplusplus::iga_refine(one_patch(arc()), o), 3);
    const Mesh m = meshioplusplus::iga_bezier_extract(one_patch(arc()), 0);
    EXPECT_EQ((*m.CellRange().begin()).Type(), "VTK_BEZIER_CURVE");
    EXPECT_EQ(m.NumPoints(), 3u);
    ASSERT_TRUE(m.HasPointData("RationalWeights"));
    EXPECT_NEAR(m.PointData("RationalWeights").As<double>()[2], std::sqrt(0.5), 1e-15);
    EXPECT_EQ(m.CellData("iga:patch", 0).As<std::int64_t>()[0], 0);
}

TEST(Iga, BezierExtractionOfASurfaceAndAVolume) {
    NurbsPatch a = meshioplusplus::iga_elevate_degree(annulus(), 1, 1, true);
    NurbsInfo ai = one_patch(a);
    IgaRefineOptions o;
    o.mSubdivisions = {2, 3, 1};
    expect_extraction_matches(meshioplusplus::iga_refine(ai, o), 6);
    expect_extraction_matches(ai, 1);
    IgaRefineOptions oh;
    oh.mSubdivisions = {2, 2, 2};
    const NurbsInfo hi = meshioplusplus::iga_refine(one_patch(hexahedron()), oh);
    expect_extraction_matches(hi, 8);
    const Mesh mh = meshioplusplus::iga_bezier_extract(hi, 0);
    EXPECT_EQ((*mh.CellRange().begin()).Type(), "VTK_BEZIER_HEXAHEDRON");
    EXPECT_FALSE(mh.HasPointData("RationalWeights"));
}

TEST(Iga, BezierExtractionRefusesUnequalDegreesAndBadPatches) {
    EXPECT_THROW(meshioplusplus::iga_bezier_extract(one_patch(annulus()), 0),
                 std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_bezier_extract(one_patch(arc()), 1), std::invalid_argument);
}

// --- Knot removal (roadmap 9.8.3) ---------------------------------------------

namespace {

NurbsPatch cubic_curve() {
    return make_patch({3}, {{0, 0, 0, 0, 0.5, 1, 1, 1, 1}}, 2,
                      {0, 0, 1, 2, 2, 3.5, 3.5, 3, 4, 0.5});
}

double max_eval_diff(const NurbsInfo& rA, const NurbsInfo& rB) {
    const NDArray u = grid(rA.mPatches[0], 9);
    const auto a = meshioplusplus::iga_evaluate(rA, 0, u);
    const auto b = meshioplusplus::iga_evaluate(rB, 0, u);
    return max_diff(a.mPoints, b.mPoints);
}

}  // namespace

TEST(Iga, RemovingInsertedKnotsRestoresTheArc) {
    const NurbsPatch original = arc();
    const NurbsPatch refined = meshioplusplus::iga_insert_knots(original, 0, {0.3, 0.7}, true);
    ASSERT_EQ(refined.mControlPoints.Shape()[0], 5u);
    meshioplusplus::IgaRemoveKnotsReport rep;
    const NurbsInfo out = meshioplusplus::iga_remove_knots(one_patch(refined), 0, 0, {}, 0.0, &rep);
    const NurbsPatch& back = out.mPatches[0];
    EXPECT_EQ(back.mKnots[0], original.mKnots[0]);
    EXPECT_EQ(rep.mNumRemoved, 2);
    ASSERT_EQ(rep.mKnots.size(), 2u);
    EXPECT_EQ(rep.mRemoved[0], 1);
    EXPECT_EQ(rep.mRemoved[1], 1);
    EXPECT_LT(max_diff(back.mControlPoints, original.mControlPoints), 1e-13);
    EXPECT_LT(max_diff(back.mWeights, original.mWeights), 1e-13);
    // One listed knot removes one; an unknown or end knot is refused.
    const NurbsInfo one =
        meshioplusplus::iga_remove_knots(one_patch(refined), 0, 0, {0.7}, 0.0, &rep);
    EXPECT_EQ(one.mPatches[0].mControlPoints.Shape()[0], 4u);
    EXPECT_EQ(rep.mNumRemoved, 1);
    EXPECT_LT(max_eval_diff(one_patch(refined), one), 1e-13);
}

TEST(Iga, AKnotTheCurveNeedsStaysUnlessTheToleranceAllowsIt) {
    const NurbsInfo in = one_patch(cubic_curve());
    meshioplusplus::IgaRemoveKnotsReport rep;
    const NurbsInfo exact = meshioplusplus::iga_remove_knots(in, 0, 0, {}, 0.0, &rep);
    EXPECT_EQ(rep.mNumRemoved, 0);
    EXPECT_EQ(exact.mPatches[0].mControlPoints.Shape()[0], 5u);
    const NurbsInfo loose = meshioplusplus::iga_remove_knots(in, 0, 0, {0.5}, 100.0, &rep);
    EXPECT_EQ(rep.mNumRemoved, 1);
    EXPECT_EQ(loose.mPatches[0].mControlPoints.Shape()[0], 4u);
    // A tolerance is a bound on the deviation of the curve.
    for (const double tol : {0.01, 0.1, 0.5, 2.0}) {
        const NurbsInfo cut = meshioplusplus::iga_remove_knots(in, 0, 0, {}, tol, &rep);
        EXPECT_LE(max_eval_diff(in, cut), tol * 1.0000001) << "tolerance " << tol;
    }
}

TEST(Iga, KnotRemovalInASurfaceAndAVolumeNeedsEveryLine) {
    NurbsPatch a = meshioplusplus::iga_insert_knots(annulus(), 1, {0.5}, true);
    a = meshioplusplus::iga_insert_knots(a, 0, {0.25, 0.25}, true);
    const NurbsInfo in = one_patch(a);
    meshioplusplus::IgaRemoveKnotsReport rep;
    const NurbsInfo back = meshioplusplus::iga_remove_knots(in, 0, 0, {}, 0.0, &rep);
    EXPECT_EQ(rep.mNumRemoved, 2);  // the double knot 0.25, both copies
    EXPECT_LT(max_eval_diff(in, back), 1e-13);
    const NurbsInfo flat = meshioplusplus::iga_remove_knots(back, 0, 1, {}, 0.0, &rep);
    EXPECT_EQ(rep.mNumRemoved, 1);
    EXPECT_LT(max_eval_diff(in, flat), 1e-13);
    EXPECT_LT(max_diff(flat.mPatches[0].mControlPoints, annulus().mControlPoints), 1e-13);

    NurbsPatch h = meshioplusplus::iga_insert_knots(hexahedron(), 2, {0.5}, false);
    const NurbsInfo hb = meshioplusplus::iga_remove_knots(one_patch(h), 0, 2, {}, 0.0, &rep);
    EXPECT_EQ(rep.mNumRemoved, 1);
    EXPECT_LT(max_diff(hb.mPatches[0].mControlPoints, hexahedron().mControlPoints), 1e-13);
}

TEST(Iga, KnotRemovalRefusesWhatItCannotDo) {
    const NurbsInfo in = one_patch(meshioplusplus::iga_insert_knots(arc(), 0, {0.4, 0.4}, true));
    EXPECT_THROW(meshioplusplus::iga_remove_knots(in, 1, 0, {}, 0.0), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_remove_knots(in, 0, 1, {}, 0.0), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_remove_knots(in, 0, 0, {}, -1.0), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_remove_knots(in, 0, 0, {}, std::nan("")),
                 std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_remove_knots(in, 0, 0, {0.6}, 0.0), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_remove_knots(in, 0, 0, {0.0}, 0.0), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_remove_knots(in, 0, 0, {0.4, 0.4, 0.4}, 0.0),
                 std::invalid_argument);
}

// --- Closest point (roadmap 9.8.4) --------------------------------------------

TEST(Iga, ClosestPointOnTheArcIsTheRadialProjection) {
    const NurbsInfo in = one_patch(arc());
    std::vector<double> pts;
    for (const double th : {0.1, 0.5, 0.9, 1.4})
        for (const double r : {0.3, 2.5})
            pts.insert(pts.end(), {r * std::cos(th), r * std::sin(th)});
    pts.insert(pts.end(), {1.0, -1.0, -1.0, 1.0});  // beyond either end
    const auto res = meshioplusplus::iga_closest_point(in, 0, doubles({pts.size() / 2, 2}, pts));
    const double* x = res.mPoints.As<double>();
    const double* d = res.mDistances.As<double>();
    for (std::size_t i = 0; i < 8; ++i) {
        const double r = std::hypot(pts[i * 2], pts[i * 2 + 1]);
        EXPECT_NEAR(x[i * 2], pts[i * 2] / r, 1e-12);
        EXPECT_NEAR(x[i * 2 + 1], pts[i * 2 + 1] / r, 1e-12);
        EXPECT_NEAR(d[i], std::abs(r - 1.0), 1e-12);
    }
    EXPECT_NEAR(x[16], 1.0, 1e-12);  // clamped to the first end
    EXPECT_NEAR(x[17], 0.0, 1e-12);
    EXPECT_NEAR(x[18], 0.0, 1e-12);  // and to the last
    EXPECT_NEAR(x[19], 1.0, 1e-12);
    EXPECT_NEAR(res.mParams.As<double>()[8], 0.0, 1e-12);
}

TEST(Iga, ClosestPointOnACubicCurveBeatsADenseSampling) {
    const NurbsInfo in = one_patch(cubic_curve());
    std::vector<double> pts;
    for (int i = 0; i < 12; ++i)
        pts.insert(pts.end(), {0.4 * i - 0.5, 4.0 - 0.5 * i});
    const NDArray q = doubles({12, 2}, pts);
    const auto res = meshioplusplus::iga_closest_point(in, 0, q);
    const NDArray u = grid(in.mPatches[0], 20001);
    const auto dense = meshioplusplus::iga_evaluate(in, 0, u);
    for (std::size_t i = 0; i < 12; ++i) {
        double best = 1e300;
        for (std::size_t k = 0; k < 20001; ++k)
            best =
                std::min(best, std::hypot(dense.mPoints.As<double>()[k * 2] - pts[i * 2],
                                          dense.mPoints.As<double>()[k * 2 + 1] - pts[i * 2 + 1]));
        EXPECT_LE(res.mDistances.As<double>()[i], best + 1e-12) << "point " << i;
        EXPECT_NEAR(res.mDistances.As<double>()[i], best, 1e-6) << "point " << i;
    }
}

TEST(Iga, ClosestPointOnASurfaceAgreesWithTheSurfaceProjector) {
    const NurbsPatch patch = annulus();
    const NurbsInfo in = one_patch(patch);
    std::vector<double> pts;
    for (const double th : {0.2, 0.7, 1.2})
        for (const double r : {1.1, 1.5, 1.9})
            for (const double z : {-0.4, 0.3})
                pts.insert(pts.end(), {r * std::cos(th), r * std::sin(th), z});
    pts.insert(pts.end(), {3.0, 3.0, 0.5, -0.5, 0.2, 0.1, 0.3, 0.1, 1.0});  // outside the patch
    const std::size_t n = pts.size() / 3;
    const auto res = meshioplusplus::iga_closest_point(in, 0, doubles({n, 3}, pts));
    const meshioplusplus::detail::BrepSurfaceEval eval(patch, true);
    const meshioplusplus::detail::SurfaceProjector projector(eval, patch);
    for (std::size_t i = 0; i < n; ++i) {
        meshioplusplus::detail::EdgeVec3 px{};
        projector.Project({pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2]}, px);
        for (std::size_t d = 0; d < 3; ++d)
            EXPECT_NEAR(res.mPoints.As<double>()[i * 3 + d], px[d], 1e-8) << "point " << i;
    }
    // A point beyond the inner edge lands on it, exactly (the projector's own
    // polish stops about 1e-9 short there).
    const double rr = std::hypot(0.3, 0.1);
    EXPECT_NEAR(res.mPoints.As<double>()[20 * 3], 0.3 / rr, 1e-12);
    EXPECT_NEAR(res.mPoints.As<double>()[20 * 3 + 1], 0.1 / rr, 1e-12);
    // On the quarter annulus in z = 0 the closest point drops the z.
    for (std::size_t i = 0; i < 18; ++i) {
        EXPECT_NEAR(res.mPoints.As<double>()[i * 3 + 2], 0.0, 1e-12);
        EXPECT_NEAR(res.mDistances.As<double>()[i], std::abs(pts[i * 3 + 2]), 1e-10);
    }
}

TEST(Iga, ClosestPointRefusesWhatItCannotDo) {
    const NDArray p3 = doubles({2, 3}, {0.5, 0.5, 0.5, 2.0, 0.5, 0.5});
    const meshioplusplus::IgaClosestPoint vr =
        meshioplusplus::iga_closest_point(one_patch(hexahedron()), 0, p3);
    EXPECT_NEAR(vr.mDistances.As<double>()[0], 0.0, 1e-12);
    // hexahedron() has a skewed +x face (x in [1.0, 1.2]), so the exterior query clamps to u = 1.
    EXPECT_NEAR(vr.mParams.As<double>()[3], 1.0, 1e-12);
    EXPECT_GT(vr.mDistances.As<double>()[1], 0.8);
    EXPECT_LT(vr.mDistances.As<double>()[1], 1.0);
    // On an exact unit cube, (2, 0.5, 0.5) projects to (1, 0.5, 0.5) with distance 1.0.
    const NurbsPatch unit_cube = make_patch(
        {1, 1, 1}, {{0, 0, 1, 1}, {0, 0, 1, 1}, {0, 0, 1, 1}}, 3,
        {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1});
    const meshioplusplus::IgaClosestPoint vc =
        meshioplusplus::iga_closest_point(one_patch(unit_cube), 0, p3);
    EXPECT_NEAR(vc.mDistances.As<double>()[0], 0.0, 1e-12);
    EXPECT_NEAR(vc.mDistances.As<double>()[1], 1.0, 1e-12);
    EXPECT_NEAR(vc.mPoints.As<double>()[3], 1.0, 1e-12);
    const NurbsInfo in = one_patch(arc());
    EXPECT_THROW(meshioplusplus::iga_closest_point(in, 1, doubles({1, 2}, {0, 0})),
                 std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_closest_point(in, 0, doubles({1, 3}, {0.5, 0.5, 0.5})),
                 std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_closest_point(in, 0, doubles({1, 2}, {std::nan(""), 0.0})),
                 std::invalid_argument);
    NurbsPatch degen;
    degen.mDegrees = {1};
    degen.mKnots = {{0.0, 0.0, 0.0, 0.0}};
    degen.mControlPoints = doubles({2, 2}, {0.0, 0.0, 1.0, 0.0});
    degen.mWeights = doubles({2}, {1.0, 1.0});
    EXPECT_THROW(
        meshioplusplus::iga_closest_point(one_patch(degen), 0, doubles({1, 2}, {0.5, 0.0})),
        std::invalid_argument);
}

TEST(Iga, ErrorsAreInvalidArguments) {
    const NurbsInfo in = one_patch(arc());
    IgaRefineOptions o;
    o.mSubdivisions = {0, 1, 1};
    EXPECT_THROW(meshioplusplus::iga_refine(in, o), std::invalid_argument);
    o.mSubdivisions = {1, 1, 1};
    o.mElevate = {-1, 0, 0};
    EXPECT_THROW(meshioplusplus::iga_refine(in, o), std::invalid_argument);
    o.mElevate = {31, 0, 0};
    EXPECT_THROW(meshioplusplus::iga_refine(in, o), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_insert_knots(arc(), 1, {0.5}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_insert_knots(arc(), 0, {0.0}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_insert_knots(arc(), 0, {0.6, 0.4}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_elevate_degree(arc(), 0, -1), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 1, grid(arc(), 3)), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 0, grid(arc(), 3), 3), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 0, grid(arc(), 3), -1), std::invalid_argument);
    NDArray out_of_range = doubles({1, 1}, {1.5});
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 0, out_of_range), std::invalid_argument);
    NDArray wrong_width = doubles({1, 2}, {0.5, 0.5});
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 0, wrong_width), std::invalid_argument);
    NurbsInfo bad = in;
    bad.mPatches[0].mWeights.As<double>()[1] = 0.0;
    EXPECT_THROW(meshioplusplus::iga_refine(bad, IgaRefineOptions{}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_refine(NurbsInfo{}, IgaRefineOptions{}),
                 std::invalid_argument);
}

// --- MFEM's patch topology ---------------------------------------------------

// The quarter annulus in MFEM's global form (see test_mfem.cpp).
std::string annulus_mesh() {
    const std::string w = "0.70710678118654757";
    const std::string path = mt::temp_path(".mesh");
    std::ofstream(path)
        << "MFEM NURBS mesh v1.0\n"
           "dimension\n2\n"
           "elements\n1\n5 3 0 1 2 3\n"
           "boundary\n0\n"
           "edges\n4\n0 0 1\n0 3 2\n1 0 3\n1 1 2\n"
           "vertices\n4\n"
           "knotvectors\n2\n2 3 0 0 0 1 1 1\n1 2 0 0 1 1\n"
           "weights\n1\n1\n1\n1\n"
        << w << "\n"
        << w
        << "\n"
           "FiniteElementSpace\nFiniteElementCollection: NURBS2\nVDim: 2\nOrdering: 1\n"
           "1 0\n0 1\n0 2\n2 0\n1 1\n2 2\n";
    return path;
}

TEST(Iga, ADefaultSampleOfAReadModelIsTheReadersMesh) {
    NurbsInfo info;
    const Mesh read =
        meshioplusplus::read_mfem(annulus_mesh(), {}, meshioplusplus::ReadOptions{}, info);
    ASSERT_EQ(info.mPatches.size(), 1u);
    NurbsInfo mapped;
    const Mesh sampled = meshioplusplus::iga_sample(info, {}, &mapped);
    EXPECT_EQ(max_diff(read.Points(), sampled.Points()), 0.0);
    EXPECT_EQ(read.NumCellBlocks(), sampled.NumCellBlocks());
    EXPECT_EQ(mapped.mCellPatch.Size(), info.mCellPatch.Size());
    EXPECT_EQ(mapped.mTopology.mNumVertices, info.mTopology.mNumVertices);
    // an explicit subdivision count is the generic sampler: no boundary cells
    IgaSampleOptions so;
    so.mSubdivisions = 1;
    EXPECT_EQ(meshioplusplus::iga_sample(info, so).NumCellBlocks(), 1u);
}

TEST(Iga, ARefinedReadModelWritesBackAsAnMfemNurbsMesh) {
    NurbsInfo info;
    meshioplusplus::read_mfem(annulus_mesh(), {}, meshioplusplus::ReadOptions{}, info);
    IgaRefineOptions o;
    o.mSubdivisions = {3, 3, 1};
    o.mElevate = {1, 1, 0};
    const NurbsInfo refined = meshioplusplus::iga_refine(info, o);
    NurbsInfo mapped;
    const Mesh mesh = meshioplusplus::iga_sample(refined, {}, &mapped);
    const std::string out = mt::temp_path(".mesh");
    meshioplusplus::write_mfem(out, mesh, false, mapped);
    NurbsInfo back;
    const Mesh again = meshioplusplus::read_mfem(out, {}, meshioplusplus::ReadOptions{}, back);
    EXPECT_EQ(max_diff(again.Points(), mesh.Points()), 0.0);
    ASSERT_EQ(back.mPatches.size(), 1u);
    EXPECT_EQ(back.mPatches[0].mDegrees, (std::vector<std::int32_t>{3, 2}));
    EXPECT_EQ(max_diff(back.mPatches[0].mControlPoints, mapped.mPatches[0].mControlPoints), 0.0);
    // refining directions differently would break the shared knot vectors
    o.mSubdivisions = {2, 3, 1};
    EXPECT_THROW(meshioplusplus::iga_refine(info, o), std::invalid_argument);
}

// --- trimmed B-reps ----------------------------------------------------------

// A unit-spaced 4 x 4 degree-1 plate, trimmed by a square loop and a diamond hole.
const char* kPlate = R"({
    "breps": [{
        "brep_id": 1,
        "faces": [{
            "brep_id": 2,
            "swapped_surface_normal": false,
            "surface": {
                "is_trimmed": true,
                "is_rational": false,
                "degrees": [1, 1],
                "knot_vectors": [[0, 4], [0, 4]],
                "control_points": [[1, [0, 0, 0, 1]], [2, [4, 0, 0, 1]],
                                   [3, [0, 4, 0, 1]], [4, [4, 4, 0, 1]]]
            },
            "boundary_loops": [
                {"loop_type": "outer", "trimming_curves": [
                    {"trim_index": 0, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[0, 0, 0, 1]], [[4, 0, 0, 1]]]}},
                    {"trim_index": 1, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[4, 0, 0, 1]], [[4, 4, 0, 1]]]}},
                    {"trim_index": 2, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[4, 4, 0, 1]], [[0, 4, 0, 1]]]}},
                    {"trim_index": 3, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[0, 4, 0, 1]], [[0, 0, 0, 1]]]}}]},
                {"loop_type": "inner", "trimming_curves": [
                    {"trim_index": 4, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 0.25, 0.5, 0.75, 1, 1], "active_range": [0, 1],
                        "control_points": [[[2, 0.8, 0, 1]], [[3.2, 2, 0, 1]], [[2, 3.2, 0, 1]],
                                           [[0.8, 2, 0, 1]], [[2, 0.8, 0, 1]]]}}]}],
            "embedded_edges": []
        }],
        "edges": [],
        "vertices": []
    }]
})";

BrepInfo plate() {
    const std::string path = mt::temp_path(".cad.json");
    std::ofstream(path) << kPlate;
    return meshioplusplus::read_cad_json_info(path);
}

TEST(Iga, ARefinedTrimmedFaceKeepsItsLoopsAndItsInsideTest) {
    const BrepInfo in = plate();
    IgaRefineOptions o;
    o.mSubdivisions = {3, 2, 1};
    o.mElevate = {1, 2, 0};
    const BrepInfo out = meshioplusplus::iga_refine(in, o);
    ASSERT_EQ(out.mFaces.size(), 1u);
    EXPECT_EQ(out.mFaces[0].mControlPointIds.Size(), 0u);
    EXPECT_EQ(out.mFaces[0].mLoops.size(), 2u);
    EXPECT_EQ(out.mPatches[0].mDegrees, (std::vector<std::int32_t>{2, 3}));
    // a grid of parameters: the plate is planar and unrefined geometry is the identity map
    NDArray u(meshioplusplus::DType::Float64, {25, 2});
    for (std::size_t j = 0; j < 5; ++j)
        for (std::size_t i = 0; i < 5; ++i) {
            u.As<double>()[(j * 5 + i) * 2] = 0.1 + 0.8 * static_cast<double>(i);
            u.As<double>()[(j * 5 + i) * 2 + 1] = 0.1 + 0.8 * static_cast<double>(j);
        }
    const auto a = meshioplusplus::iga_evaluate(in, 0, u, 1);
    const auto b = meshioplusplus::iga_evaluate(out, 0, u, 1);
    EXPECT_LT(max_diff(a.mPoints, b.mPoints), 1e-13);
    EXPECT_LT(max_diff(a.mNormals, b.mNormals), 1e-13);
    ASSERT_EQ(a.mInside.Size(), 25u);
    std::size_t inside = 0;
    for (std::size_t i = 0; i < 25; ++i) {
        EXPECT_EQ(a.mInside.As<std::uint8_t>()[i], b.mInside.As<std::uint8_t>()[i]);
        inside += a.mInside.As<std::uint8_t>()[i];
    }
    // the diamond |x-2| + |y-2| < 1.2 holds the centre of the grid and no corner
    EXPECT_EQ(a.mInside.As<std::uint8_t>()[12], 0);
    EXPECT_EQ(a.mInside.As<std::uint8_t>()[0], 1);
    EXPECT_LT(inside, 25u);
    // sampling the refined face is sampling CAD JSON as it does
    IgaSampleOptions so;
    so.mSubdivisions = 8;
    const Mesh mesh = meshioplusplus::iga_sample(out, so);
    EXPECT_GT(mesh.NumCellBlocks(), 0u);
}

TEST(Iga, FacesSharingAPatchMustAgreeOnRational) {
    BrepInfo in = plate();
    in.mFaces.push_back(in.mFaces[0]);
    in.mFaces[1].mRational = true;
    EXPECT_THROW(meshioplusplus::iga_refine(in, IgaRefineOptions{}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 5, grid(in.mPatches[0], 3)),
                 std::invalid_argument);
}

}  // namespace
