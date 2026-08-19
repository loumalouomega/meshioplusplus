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
// Tests for the Hessian composition (`operations/hessian.hpp`).
//
// The one MESH-SHAPE-INDEPENDENT guarantee is "a linear field has an
// exactly-zero Hessian" (see hessian.hpp's own header comment): a linear
// field's gradient is EXACTLY constant (not merely linear) after the first
// Green-Gauss + Point-averaging pass, and Green-Gauss of a spatially constant
// field is trivially exact, so the second pass gives exactly zero regardless
// of mesh regularity.
//
// For a genuinely quadratic field, the tests below measure (not assume) the
// honest picture: on a regular axis-aligned hex grid, every INTERIOR node's
// 8-cell neighbourhood is symmetric enough for the intermediate
// Uniform-weighted Point-averaging step to be exact too, so the composed
// Hessian comes back exact to machine precision there
// (QuadraticFieldOnAStructuredGridIsExactAwayFromTheBoundary). The SAME
// mesh's own boundary cells, and any genuinely irregular cell, have a
// one-sided/asymmetric neighbourhood where that averaging step is only
// approximate -- real, bounded, non-zero error, pinned with an
// empirically-calibrated (not derived) tolerance so a regression that makes
// it silently worse is caught.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/hessian.hpp"

namespace {

using meshioplusplus::DataLocation;
using meshioplusplus::DType;
using meshioplusplus::GradientMethod;
using meshioplusplus::hessian;
using meshioplusplus::HessianOptions;
using meshioplusplus::HessianResult;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

constexpr double kTight = 1e-9;

HessianOptions opts(const std::string& rArray, GradientMethod method = GradientMethod::GreenGauss) {
    HessianOptions o;
    o.mArrayName = rArray;
    o.mMethod = method;
    return o;
}

// A structured nx x ny x nz hexahedron grid on [0,nx]x[0,ny]x[0,nz].
Mesh hex_grid(std::size_t nx, std::size_t ny, std::size_t nz) {
    std::vector<std::vector<double>> pts;
    for (std::size_t k = 0; k <= nz; ++k)
        for (std::size_t j = 0; j <= ny; ++j)
            for (std::size_t i = 0; i <= nx; ++i)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    auto vid = [&](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<std::int64_t>((k * (ny + 1) + j) * (nx + 1) + i);
    };
    std::vector<std::vector<std::int64_t>> cells;
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i)
                cells.push_back({vid(i, j, k), vid(i + 1, j, k), vid(i + 1, j + 1, k),
                                 vid(i, j + 1, k), vid(i, j, k + 1), vid(i + 1, j, k + 1),
                                 vid(i + 1, j + 1, k + 1), vid(i, j + 1, k + 1)});
    return mt::make_mesh(std::move(pts), "hexahedron", std::move(cells));
}

// A frustum: a 2x2 base at z=0 tapering to a 1x1 top at z=1 -- an irregular
// cell whose faces are planar trapezoids, used to show the quadratic-field
// composition is approximate (not exact) off a structured grid.
Mesh tapered_hex() {
    return mt::make_mesh({{0.0, 0.0, 0.0},
                          {2.0, 0.0, 0.0},
                          {2.0, 2.0, 0.0},
                          {0.0, 2.0, 0.0},
                          {0.5, 0.5, 1.0},
                          {1.5, 0.5, 1.0},
                          {1.5, 1.5, 1.0},
                          {0.5, 1.5, 1.0}},
                         "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
}

// f = a.x + d, linear, sampled at every point.
NDArray linear_scalar(const std::vector<std::vector<double>>& rPts, double a, double b, double c,
                      double d) {
    std::vector<double> v;
    v.reserve(rPts.size());
    for (const auto& r_p : rPts)
        v.push_back(a * r_p[0] + b * r_p[1] + c * r_p[2] + d);
    return mt::data_array(v);
}

// f = 1/2 x^T H x, a general quadratic whose Hessian is the constant matrix H
// (row-major, symmetric by construction). No linear/constant term, so the
// analytic Hessian is exactly H everywhere.
NDArray quadratic_scalar(const std::vector<std::vector<double>>& rPts, const double H[9]) {
    std::vector<double> v;
    v.reserve(rPts.size());
    for (const auto& r_p : rPts) {
        double acc = 0.0;
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                acc += r_p[i] * H[i * 3 + j] * r_p[j];
        v.push_back(0.5 * acc);
    }
    return mt::data_array(v);
}

std::vector<std::vector<double>> mesh_points(const Mesh& rMesh) {
    std::vector<std::vector<double>> out;
    const NDArray& p = rMesh.Points();
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<double> row(dim);
        for (std::size_t d = 0; d < dim; ++d)
            row[d] = meshioplusplus::detail::read_double(p, i * dim + d);
        out.push_back(std::move(row));
    }
    return out;
}

}  // namespace

// --- linear field: exactly zero Hessian, any mesh -----------------------

TEST(Hessian, LinearFieldOnATaperedHexIsExactlyZero) {
    Mesh m = tapered_hex();
    m.AddPointData("f", linear_scalar(mesh_points(m), 3.0, -2.0, 5.0, 7.0));
    const HessianResult r = hessian(m, opts("f"));
    const NDArray& h = r.mMesh.CellData("f:hessian", 0);
    for (std::size_t i = 0; i < 9; ++i)
        EXPECT_NEAR(meshioplusplus::detail::read_double(h, i), 0.0, kTight) << "entry " << i;
    EXPECT_EQ(r.mNumSkipped, 0);
}

TEST(Hessian, LinearFieldOnAHexGridIsExactlyZero) {
    Mesh m = hex_grid(3, 3, 3);
    m.AddPointData("f", linear_scalar(mesh_points(m), 0.5, -1.25, 2.0, 4.0));
    const HessianResult r = hessian(m, opts("f"));
    const NDArray& h = r.mMesh.CellData("f:hessian", 0);
    const std::size_t ncells = r.mMesh.Cells(0).NumCells();
    for (std::size_t c = 0; c < ncells; ++c)
        for (std::size_t i = 0; i < 9; ++i)
            EXPECT_NEAR(meshioplusplus::detail::read_double(h, c * 9 + i), 0.0, kTight);
}

// --- quadratic field: a good but honestly-approximate estimate -----------

TEST(Hessian, QuadraticFieldOnAStructuredGridIsExactAwayFromTheBoundary) {
    // H = [[2,1,0],[1,3,1],[0,1,4]] (symmetric). On a regular axis-aligned
    // grid, every INTERIOR node has a perfectly symmetric 8-cell
    // neighbourhood, so the Uniform-weighted Point-averaging step between the
    // two gradient() passes reproduces the true (genuinely varying) nodal
    // gradient exactly rather than merely approximately -- measured, not
    // assumed: interior cells come back with zero error to machine precision,
    // strictly stronger than the general "approximate on an irregular mesh"
    // caveat this operation otherwise carries (see the boundary/irregular
    // tests below, where the symmetry breaks down and real error appears).
    const double H[9] = {2, 1, 0, 1, 3, 1, 0, 1, 4};
    const std::size_t n = 8;
    Mesh m = hex_grid(n, n, n);
    m.AddPointData("f", quadratic_scalar(mesh_points(m), H));
    const HessianResult r = hessian(m, opts("f"));
    const NDArray& h = r.mMesh.CellData("f:hessian", 0);
    // At least two rings in from every face, so every incident node's own
    // 1-ring is itself fully interior too.
    for (std::size_t k = 2; k < n - 2; ++k)
        for (std::size_t j = 2; j < n - 2; ++j)
            for (std::size_t i = 2; i < n - 2; ++i) {
                const std::size_t c = (k * n + j) * n + i;
                for (std::size_t comp = 0; comp < 9; ++comp) {
                    const double v = meshioplusplus::detail::read_double(h, c * 9 + comp);
                    EXPECT_NEAR(v, H[comp], kTight) << "cell " << c << " entry " << comp;
                }
            }
}

TEST(Hessian, QuadraticFieldOnAStructuredGridHasBoundedBoundaryError) {
    // The boundary cells of the SAME structured grid, by contrast, sit on
    // asymmetric (one-sided) node neighbourhoods -- this is where the
    // Uniform-averaging approximation actually shows up, even on an
    // otherwise-regular mesh. Bounded and non-zero, not exact: proves the
    // interior test above is measuring a real distinction, not a universal
    // exactness the composition doesn't actually have.
    const double H[9] = {2, 1, 0, 1, 3, 1, 0, 1, 4};
    const std::size_t n = 8;
    Mesh m = hex_grid(n, n, n);
    m.AddPointData("f", quadratic_scalar(mesh_points(m), H));
    const HessianResult r = hessian(m, opts("f"));
    const NDArray& h = r.mMesh.CellData("f:hessian", 0);
    const std::size_t ncells = r.mMesh.Cells(0).NumCells();
    double max_err = 0.0;
    for (std::size_t c = 0; c < ncells; ++c)
        for (std::size_t comp = 0; comp < 9; ++comp) {
            const double v = meshioplusplus::detail::read_double(h, c * 9 + comp);
            max_err = std::max(max_err, std::fabs(v - H[comp]));
        }
    // Empirically calibrated, not a claimed exactness bound.
    EXPECT_GT(max_err, 1e-6);  // proves the boundary error is real
    EXPECT_LT(max_err, 3.0);   // ... but bounded; a regression trip-wire
}

TEST(Hessian, QuadraticFieldOnAnIrregularMeshIsApproximateButBounded) {
    const double H[9] = {2, 0, 0, 0, 2, 0, 0, 0, 2};
    Mesh m = tapered_hex();
    m.AddPointData("f", quadratic_scalar(mesh_points(m), H));
    const HessianResult r = hessian(m, opts("f"));
    const NDArray& h = r.mMesh.CellData("f:hessian", 0);
    double max_err = 0.0;
    for (std::size_t i = 0; i < 9; ++i)
        max_err = std::max(max_err, std::fabs(meshioplusplus::detail::read_double(h, i) - H[i]));
    // A single irregular cell: bounded, but not claimed tight.
    EXPECT_LT(max_err, 3.0);
    EXPECT_GT(max_err, 1e-6);  // proves the approximation is REAL, not accidentally exact
}

// --- shape / contract ------------------------------------------------------

TEST(Hessian, ScalarInputGivesNineComponents) {
    Mesh m = hex_grid(2, 2, 2);
    m.AddPointData("f", linear_scalar(mesh_points(m), 1.0, 0.0, 0.0, 0.0));
    const HessianResult r = hessian(m, opts("f"));
    ASSERT_TRUE(r.mMesh.HasCellData("f:hessian"));
    const NDArray& h = r.mMesh.CellData("f:hessian", 0);
    EXPECT_EQ(h.Shape().size(), 2u);
    EXPECT_EQ(h.Shape()[1], 9u);
}

TEST(Hessian, RejectsAMultiComponentInputByName) {
    Mesh m = hex_grid(2, 2, 2);
    const std::size_t n = m.NumPoints();
    m.AddPointData("u", mt::data_array(std::vector<double>(n * 3, 0.0), 3));
    HessianOptions o = opts("u");
    try {
        hessian(m, o);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("scalar fields only"), std::string::npos);
    }
}

TEST(Hessian, RejectsACellDataInputByName) {
    Mesh m = hex_grid(2, 2, 2);
    const std::size_t nc = m.Cells(0).NumCells();
    m.AddCellData("c", {mt::data_array(std::vector<double>(nc, 0.0))});
    try {
        hessian(m, opts("c"));
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("cell_data_to_point_data"), std::string::npos);
    }
}

TEST(Hessian, RejectsAnUnknownArrayName) {
    Mesh m = hex_grid(2, 2, 2);
    m.AddPointData("f", linear_scalar(mesh_points(m), 1.0, 0.0, 0.0, 0.0));
    EXPECT_THROW(hessian(m, opts("nope")), std::invalid_argument);
}

TEST(Hessian, PointLocationAttachesPointData) {
    Mesh m = hex_grid(3, 3, 3);
    m.AddPointData("f", linear_scalar(mesh_points(m), 1.0, -1.0, 2.0, 0.0));
    HessianOptions o = opts("f");
    o.mLocation = DataLocation::Point;
    const HessianResult r = hessian(m, o);
    ASSERT_TRUE(r.mMesh.HasPointData("f:hessian"));
    EXPECT_FALSE(r.mMesh.HasCellData("f:hessian"));
    const NDArray& h = r.mMesh.PointData("f:hessian");
    EXPECT_EQ(h.Shape()[0], r.mMesh.NumPoints());
    EXPECT_EQ(h.Shape()[1], 9u);
}

TEST(Hessian, OutputNameAndOverwrite) {
    Mesh m = hex_grid(2, 2, 2);
    m.AddPointData("f", linear_scalar(mesh_points(m), 1.0, 0.0, 0.0, 0.0));
    HessianOptions o = opts("f");
    o.mOutputName = "curv";
    const HessianResult r1 = hessian(m, o);
    ASSERT_TRUE(r1.mMesh.HasCellData("curv"));

    // Without overwrite, a second call onto an existing name fails.
    EXPECT_THROW(hessian(r1.mMesh, o), std::invalid_argument);
    o.mOverwrite = true;
    EXPECT_NO_THROW(hessian(r1.mMesh, o));
}

TEST(Hessian, GeometryAndExistingDataPassThroughUnchanged) {
    Mesh m = hex_grid(2, 2, 2);
    m.AddPointData("f", linear_scalar(mesh_points(m), 1.0, 2.0, 3.0, 0.0));
    m.AddCellData("mat", {mt::data_array({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0})});
    const HessianResult r = hessian(m, opts("f"));
    EXPECT_EQ(r.mMesh.NumPoints(), m.NumPoints());
    EXPECT_TRUE(r.mMesh.HasCellData("mat"));
    const NDArray& before = m.CellData("mat", 0);
    const NDArray& after = r.mMesh.CellData("mat", 0);
    for (std::size_t i = 0; i < 8; ++i)
        EXPECT_EQ(meshioplusplus::detail::read_double(before, i),
                 meshioplusplus::detail::read_double(after, i));
}

// --- determinism -----------------------------------------------------------

TEST(Hessian, RepeatRunsByteIdentical) {
    Mesh m = hex_grid(4, 4, 4);
    const double H[9] = {1, 0.3, 0, 0.3, 2, 0.1, 0, 0.1, 1.5};
    m.AddPointData("f", quadratic_scalar(mesh_points(m), H));
    const HessianResult r1 = hessian(m, opts("f"));
    const HessianResult r2 = hessian(m, opts("f"));
    const NDArray& h1 = r1.mMesh.CellData("f:hessian", 0);
    const NDArray& h2 = r2.mMesh.CellData("f:hessian", 0);
    const std::size_t n = h1.Size();
    for (std::size_t i = 0; i < n; ++i)
        EXPECT_EQ(meshioplusplus::detail::read_double(h1, i),
                 meshioplusplus::detail::read_double(h2, i));
}

TEST(Hessian, LeastSquaresMethodRuns) {
    Mesh m = hex_grid(3, 3, 3);
    m.AddPointData("f", linear_scalar(mesh_points(m), 1.0, -2.0, 0.5, 3.0));
    const HessianResult r = hessian(m, opts("f", GradientMethod::LeastSquares));
    ASSERT_TRUE(r.mMesh.HasCellData("f:hessian"));
    const NDArray& h = r.mMesh.CellData("f:hessian", 0);
    const std::size_t ncells = r.mMesh.Cells(0).NumCells();
    for (std::size_t c = 0; c < ncells; ++c)
        for (std::size_t i = 0; i < 9; ++i)
            EXPECT_NEAR(meshioplusplus::detail::read_double(h, c * 9 + i), 0.0, 1e-6);
}
