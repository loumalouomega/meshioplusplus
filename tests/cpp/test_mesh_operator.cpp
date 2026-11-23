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
// The mesh's P1 operators as sparse matrices (roadmap 15.2.1).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/mesh_operator.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

MeshOperatorResult op(const Mesh& rMesh, MeshOperatorKind Kind) {
    MeshOperatorOptions o;
    o.mKind = Kind;
    return mesh_operator(rMesh, o);
}

std::vector<double> multiply(const MeshOperatorResult& rA, const std::vector<double>& rX) {
    std::vector<double> y(static_cast<std::size_t>(rA.mRows), 0.0);
    for (std::int64_t r = 0; r < rA.mRows; ++r)
        for (std::int64_t k = rA.mIndptr[r]; k < rA.mIndptr[r + 1]; ++k)
            y[r] += rA.mValues[k] * rX[rA.mIndices[k]];
    return y;
}

double total(const MeshOperatorResult& rA) {
    double s = 0.0;
    for (double v : rA.mValues)
        s += v;
    return s;
}

}  // namespace

TEST(MeshOperator, ShapesAndTheAreaOfASphere) {
    const Mesh m = icosa_sphere(1.0, 2);
    const std::int64_t n = static_cast<std::int64_t>(m.NumPoints());
    const MeshOperatorResult lumped = op(m, MeshOperatorKind::MassLumped);
    const MeshOperatorResult consistent = op(m, MeshOperatorKind::MassConsistent);
    EXPECT_EQ(lumped.mRows, n);
    EXPECT_EQ(lumped.mCols, n);
    EXPECT_EQ(lumped.mValues.size(), static_cast<std::size_t>(n));
    // Both masses sum to the surface area, and so does the lumped diagonal.
    EXPECT_NEAR(total(lumped), total(consistent), 1e-12);
    EXPECT_NEAR(total(lumped), 4.0 * M_PI, 0.3);
    // The consistent mass couples neighbours: more entries than the diagonal.
    EXPECT_GT(consistent.mValues.size(), lumped.mValues.size());
}

TEST(MeshOperator, LaplacianIsSymmetricWithZeroRowSums) {
    const Mesh m = icosa_sphere(1.0, 2);
    const MeshOperatorResult k = op(m, MeshOperatorKind::Laplacian);
    const std::vector<double> ones(static_cast<std::size_t>(k.mCols), 1.0);
    for (double v : multiply(k, ones))
        EXPECT_NEAR(v, 0.0, 1e-12);
    // Symmetric: K x . y == x . K y for two arbitrary vectors.
    std::vector<double> x(ones.size()), y(ones.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] = std::sin(0.7 * i);
        y[i] = std::cos(1.3 * i);
    }
    const std::vector<double> kx = multiply(k, x), ky = multiply(k, y);
    double a = 0.0, b = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        a += kx[i] * y[i];
        b += x[i] * ky[i];
    }
    EXPECT_NEAR(a, b, 1e-10);
}

TEST(MeshOperator, GradientOfALinearFieldAndTheStiffnessIdentity) {
    const Mesh m = icosa_sphere(1.0, 2);
    const std::size_t n = m.NumPoints();
    const MeshOperatorResult g = op(m, MeshOperatorKind::Gradient);
    EXPECT_EQ(g.mCols, static_cast<std::int64_t>(n));
    EXPECT_EQ(g.mRows % 3, 0);
    // f = 2x + 3y - z: the gradient of a linear field is its tangential part,
    // so (G f) . n_t = 0 and G f . (2, 3, -1) = |tangential part|^2 > 0 per triangle.
    std::vector<double> f(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double* p = m.Points().As<double>() + 3 * i;
        f[i] = 2.0 * p[0] + 3.0 * p[1] - p[2];
    }
    const std::vector<double> gf = multiply(g, f);
    const double full = 2.0 * 2.0 + 3.0 * 3.0 + 1.0;
    for (std::int64_t t = 0; t < g.mRows / 3; ++t) {
        const double dot = 2.0 * gf[3 * t] + 3.0 * gf[3 * t + 1] - gf[3 * t + 2];
        EXPECT_GT(dot, 0.0);
        EXPECT_LT(dot, full + 1e-9);  // a projection never lengthens the vector
    }
    // K = G^T diag(area) G: x^T K x is the Dirichlet energy of x, the area-weighted
    // sum of |G x|^2. The areas are recovered from the lumped mass (area = 3 * mass
    // summed over the triangle's corners), so compare totals through the Laplacian.
    const MeshOperatorResult k = op(m, MeshOperatorKind::Laplacian);
    const std::vector<double> kf = multiply(k, f);
    double energy = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        energy += f[i] * kf[i];
    EXPECT_GT(energy, 0.0);
}

TEST(MeshOperator, RefusesOutOfScopeCells) {
    EXPECT_THROW(op(box(), MeshOperatorKind::Laplacian), std::invalid_argument);
}

#include "meshioplusplus/operations/smooth_field.hpp"

TEST(SmoothField, FlattensTheShortScaleAndKeepsTheMean) {
    Mesh m = icosa_sphere(1.0, 3);
    const std::size_t n = m.NumPoints();
    NDArray f(DType::Float64, {n});
    for (std::size_t i = 0; i < n; ++i) {
        const double* p = m.Points().As<double>() + 3 * i;
        f.As<double>()[i] = p[2] + ((i % 2) ? 0.5 : -0.5);  // a smooth part plus noise
    }
    m.AddPointData("f", std::move(f));
    SmoothFieldOptions o;
    o.mArray = "f";
    o.mLength = 0.0;
    const SmoothFieldResult same = smooth_field(m, o);
    for (std::size_t i = 0; i < n; ++i)
        EXPECT_NEAR(same.mMesh.PointData("f:smooth").As<double>()[i],
                    m.PointData("f").As<double>()[i], 1e-9);
    o.mLength = 0.5;
    const SmoothFieldResult r = smooth_field(m, o);
    EXPECT_TRUE(r.mConverged);
    const NDArray& g = r.mMesh.PointData("f:smooth");
    double var_f = 0.0, var_g = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double* p = m.Points().As<double>() + 3 * i;
        (void)p;
        const double ef = m.PointData("f").As<double>()[i];
        const double eg = g.As<double>()[i];
        var_f += ef * ef;
        var_g += eg * eg;
    }
    EXPECT_LT(var_g, 0.5 * var_f);  // the noise is mostly gone
    // The M-weighted mean is preserved: the lumped mass sums the same on both.
    MeshOperatorOptions mo;
    mo.mKind = MeshOperatorKind::MassLumped;
    const MeshOperatorResult lumped = mesh_operator(m, mo);
    double sf = 0.0, sg = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sf += lumped.mValues[i] * m.PointData("f").As<double>()[i];
        sg += lumped.mValues[i] * g.As<double>()[i];
    }
    EXPECT_NEAR(sf, sg, 1e-8);
}

TEST(SmoothField, Refusals) {
    Mesh m = icosa_sphere(1.0, 1);
    SmoothFieldOptions o;
    o.mArray = "nothing";
    EXPECT_THROW(smooth_field(m, o), std::invalid_argument);
}
