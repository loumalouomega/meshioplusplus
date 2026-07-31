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
// Tests for the field differential operators (`operations/gradient.hpp`).
//
// Several "obvious" tests are inert here and are deliberately absent or
// strengthened; the fixtures below exist to defeat exactly those traps:
//  - "the gradient of a constant field is zero" proves nothing: the numerator is
//    f_bar * sum(A_j) and sum(A_j) == 0 over ANY closed surface with ANY
//    quadrature weights and ANY volume.
//  - a cube or an axis-aligned hex grid cannot distinguish the face fan from a
//    naive corner-average, because a parallelogram's corner average IS its area
//    centroid. Hence `tapered_hex` (planar trapezoidal side faces) and
//    `warped_hex` (a genuinely non-planar face).
//  - an all-tetra suite cannot see any quad-face bug at all, since a triangle's
//    fan is redundant.
//  - a solenoidal field cannot see a swapped divergence/curl; a curl of the form
//    (0,0,c) cannot see an index permutation. The vector fixtures use three
//    distinct nonzero answers so any permutation or sign slip fails.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/gradient.hpp"

namespace {

using meshioplusplus::DataLocation;
using meshioplusplus::DType;
using meshioplusplus::gradient;
using meshioplusplus::gradient_method_from_name;
using meshioplusplus::gradient_operator_from_name;
using meshioplusplus::GradientMethod;
using meshioplusplus::GradientOperator;
using meshioplusplus::GradientOptions;
using meshioplusplus::GradientResult;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

// A tolerance appropriate for "exact up to rounding": the gradient is one
// division of two quantities each accumulated over a handful of terms.
constexpr double kTight = 1e-12;

GradientOptions opts(const std::string& rArray,
                     GradientMethod method = GradientMethod::GreenGauss,
                     GradientOperator op = GradientOperator::Gradient) {
    GradientOptions o;
    o.mArrayName = rArray;
    o.mMethod = method;
    o.mOperator = op;
    return o;
}

// --- fixtures ---------------------------------------------------------------

// A frustum: a 2x2 base at z = 0 tapering to a 1x1 top at z = 1. Every side face
// is a PLANAR TRAPEZOID, whose corner average is NOT its area centroid -- this
// is the fixture that a naive corner-average face value fails and the face fan
// passes.
std::vector<std::vector<double>> tapered_hex_points() {
    return {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 2.0, 0.0}, {0.0, 2.0, 0.0},
            {0.5, 0.5, 1.0}, {1.5, 0.5, 1.0}, {1.5, 1.5, 1.0}, {0.5, 1.5, 1.0}};
}

// A cube with one top corner lifted, so the top face (4,5,6,7) is genuinely
// non-planar. Green-Gauss stays exact here because the fan surface is closed
// regardless of face planarity -- that is the claim this pins.
std::vector<std::vector<double>> warped_hex_points() {
    return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.4}, {0.0, 1.0, 1.0}};
}

// A planar trapezoid in the xy-plane: the 2D counterpart of the frustum trap.
std::vector<std::vector<double>> trapezoid_points() {
    return {{0.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};
}

// f = a.x + d sampled at every point of `rPts`.
NDArray linear_scalar(const std::vector<std::vector<double>>& rPts, double a, double b, double c,
                      double d) {
    std::vector<double> v;
    v.reserve(rPts.size());
    for (const auto& r_p : rPts)
        v.push_back(a * r_p[0] + b * r_p[1] + c * r_p[2] + d);
    return mt::data_array(v);
}

// u_i = sum_j C[i][j] x_j sampled at every point, as an (n, 3) array.
NDArray linear_vector(const std::vector<std::vector<double>>& rPts, const double C[9]) {
    std::vector<double> v;
    v.reserve(rPts.size() * 3);
    for (const auto& r_p : rPts)
        for (std::size_t i = 0; i < 3; ++i)
            v.push_back(C[i * 3 + 0] * r_p[0] + C[i * 3 + 1] * r_p[1] + C[i * 3 + 2] * r_p[2]);
    return mt::data_array(v, 3);
}

// A structured nx x ny x nz hexahedron grid on [0,nx]x[0,ny]x[0,nz], shifted by
// `off` so the same builder serves the translation-invariance test.
Mesh hex_grid(std::size_t nx, std::size_t ny, std::size_t nz, double off = 0.0) {
    std::vector<std::vector<double>> pts;
    for (std::size_t k = 0; k <= nz; ++k)
        for (std::size_t j = 0; j <= ny; ++j)
            for (std::size_t i = 0; i <= nx; ++i)
                pts.push_back({static_cast<double>(i) + off, static_cast<double>(j) + off,
                               static_cast<double>(k) + off});
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
    Mesh m = mt::make_mesh(std::move(pts), "hexahedron", std::move(cells));
    return m;
}

// Asserts a cell's 3-component gradient equals (a, b, c).
void expect_grad(const Mesh& rMesh, const std::string& rName, std::size_t Block, std::size_t Cell,
                 double a, double b, double c, double Tol = kTight) {
    const NDArray& arr = rMesh.CellData(rName, Block);
    EXPECT_NEAR(meshioplusplus::detail::read_double(arr, Cell * 3 + 0), a, Tol);
    EXPECT_NEAR(meshioplusplus::detail::read_double(arr, Cell * 3 + 1), b, Tol);
    EXPECT_NEAR(meshioplusplus::detail::read_double(arr, Cell * 3 + 2), c, Tol);
}

}  // namespace

// --- the headline invariant: a linear field is differentiated exactly --------

TEST(Gradient, LinearFieldOnATaperedHexIsExact) {
    // The trap this exists for: on a cube every side face is a parallelogram,
    // whose corner average IS its area centroid, so a naive (no-fan) face value
    // would pass. A frustum's trapezoidal faces separate the two.
    const auto pts = tapered_hex_points();
    Mesh m = mt::make_mesh(pts, "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    m.AddPointData("f", linear_scalar(pts, 3.0, -2.0, 5.0, 7.0));

    const GradientResult r = gradient(m, opts("f"));
    EXPECT_EQ(r.mNumSkipped, 0);
    EXPECT_EQ(r.mNumFallback, 0);
    expect_grad(r.mMesh, "f:gradient", 0, 0, 3.0, -2.0, 5.0);
}

TEST(Gradient, LinearFieldOnAWarpedHexIsExact) {
    // Exactness does NOT depend on face planarity: two faces sharing an edge
    // contribute oppositely-wound fan triangles there, so the fan surface is
    // closed and the divergence theorem applies verbatim.
    const auto pts = warped_hex_points();
    Mesh m = mt::make_mesh(pts, "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    m.AddPointData("f", linear_scalar(pts, -1.5, 4.0, 0.25, -3.0));

    const GradientResult r = gradient(m, opts("f"));
    EXPECT_EQ(r.mNumSkipped, 0);
    expect_grad(r.mMesh, "f:gradient", 0, 0, -1.5, 4.0, 0.25);
}

TEST(Gradient, LinearFieldOnATetAndAWedgeIsExact) {
    const std::vector<std::vector<double>> tp = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    Mesh t = mt::make_mesh(tp, "tetra", {{0, 1, 2, 3}});
    t.AddPointData("f", linear_scalar(tp, 2.0, 3.0, 5.0, 1.0));
    expect_grad(gradient(t, opts("f")).mMesh, "f:gradient", 0, 0, 2.0, 3.0, 5.0);

    const std::vector<std::vector<double>> wp = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                                {0.0, 1.0, 0.0}, {0.0, 0.0, 2.0},
                                                {1.0, 0.0, 2.0}, {0.0, 1.0, 2.0}};
    Mesh w = mt::make_mesh(wp, "wedge", {{0, 1, 2, 3, 4, 5}});
    w.AddPointData("f", linear_scalar(wp, -4.0, 6.0, 0.5, 0.0));
    expect_grad(gradient(w, opts("f")).mMesh, "f:gradient", 0, 0, -4.0, 6.0, 0.5);
}

TEST(Gradient, LinearFieldOnAPlanarTrapezoidIsExact) {
    // 2D: the in-plane gradient of a planar cell. The z-component must be
    // exactly zero, not merely small, because the formula projects onto the
    // cell's own plane.
    const auto pts = trapezoid_points();
    Mesh m = mt::make_mesh(pts, "quad", {{0, 1, 2, 3}});
    m.AddPointData("f", linear_scalar(pts, 1.25, -0.75, 0.0, 9.0));

    const GradientResult r = gradient(m, opts("f"));
    EXPECT_EQ(r.mNumSkipped, 0);
    expect_grad(r.mMesh, "f:gradient", 0, 0, 1.25, -0.75, 0.0);
}

TEST(Gradient, LeastSquaresIsExactForALinearField) {
    Mesh m = hex_grid(3, 3, 3);
    std::vector<std::vector<double>> pts;
    for (std::size_t k = 0; k <= 3; ++k)
        for (std::size_t j = 0; j <= 3; ++j)
            for (std::size_t i = 0; i <= 3; ++i)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    m.AddPointData("f", linear_scalar(pts, 0.5, -1.25, 2.0, 4.0));

    const GradientResult r = gradient(m, opts("f", GradientMethod::LeastSquares));
    EXPECT_EQ(r.mNumSkipped, 0);
    EXPECT_EQ(r.mNumFallback, 0) << "a 3x3x3 grid has a full-rank neighbourhood everywhere";
    for (std::size_t c = 0; c < 27; ++c)
        expect_grad(r.mMesh, "f:gradient", 0, c, 0.5, -1.25, 2.0, 1e-10);
}

TEST(Gradient, TranslatingTheMeshDoesNotMoveTheGradient) {
    // The gate for recentring. V = (1/3) sum x_j . A_j only telescopes because
    // sum A_j == 0; without recentring the terms are ~off * area and cancel down
    // to ~area, costing log10(off) digits of a quantity we then divide by -- and
    // the same for the numerator, whose field values are ~off too.
    //
    // The offset and the tolerance are both load-bearing. At off = 1e6 the
    // un-recentred error is ~1e-10 absolute, which a 1e-8 tolerance does not
    // catch: this test was verified to be INERT in that form and to fire at
    // 1e8 / 1e-9, where the un-recentred error is ~1e-8. With recentring the
    // answer is exact to rounding at any offset.
    const double off = 1.0e8;
    for (GradientMethod method : {GradientMethod::GreenGauss, GradientMethod::LeastSquares}) {
        Mesh far = hex_grid(2, 2, 2, off);
        std::vector<std::vector<double>> pts;
        for (std::size_t k = 0; k <= 2; ++k)
            for (std::size_t j = 0; j <= 2; ++j)
                for (std::size_t i = 0; i <= 2; ++i)
                    pts.push_back({static_cast<double>(i) + off, static_cast<double>(j) + off,
                                   static_cast<double>(k) + off});
        far.AddPointData("f", linear_scalar(pts, 3.0, -2.0, 5.0, 0.0));
        const GradientResult r = gradient(far, opts("f", method));
        for (std::size_t c = 0; c < 8; ++c)
            expect_grad(r.mMesh, "f:gradient", 0, c, 3.0, -2.0, 5.0, 1e-9);
    }
}

// --- vector operators -------------------------------------------------------

TEST(Gradient, DivergenceOfAnAnisotropicFieldIsTheExactConstant) {
    // u = (2x, 3y, 5z) -> div = 10 from three DISTINCT nonzero diagonal terms,
    // so a wrong diagonal or a swapped div/curl cannot pass. A solenoidal field
    // (div = 0) would be inert.
    const auto pts = tapered_hex_points();
    Mesh m = mt::make_mesh(pts, "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    const double C[9] = {2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 5.0};
    m.AddPointData("u", linear_vector(pts, C));

    const GradientResult r = gradient(m, opts("u", GradientMethod::GreenGauss,
                                              GradientOperator::Divergence));
    ASSERT_TRUE(r.mMesh.HasCellData("u:divergence"));
    EXPECT_NEAR(meshioplusplus::detail::read_double(r.mMesh.CellData("u:divergence", 0), 0), 10.0,
                kTight);
}

TEST(Gradient, CurlHasThreeDistinctNonzeroComponents) {
    // u = (7z, 11x, 13y) -> curl = (13, 7, 11). Distinct and nonzero in every
    // slot, so ANY index permutation or sign flip in the curl formula fails.
    const auto pts = warped_hex_points();
    Mesh m = mt::make_mesh(pts, "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    const double C[9] = {0.0, 0.0, 7.0, 11.0, 0.0, 0.0, 0.0, 13.0, 0.0};
    m.AddPointData("u", linear_vector(pts, C));

    const GradientResult r =
        gradient(m, opts("u", GradientMethod::GreenGauss, GradientOperator::Curl));
    expect_grad(r.mMesh, "u:curl", 0, 0, 13.0, 7.0, 11.0);
}

TEST(Gradient, TensorLayoutIsComponentMajorThenDerivative) {
    // All nine entries of C are distinct, so the i*3+j <-> j*3+i transpose fails.
    const auto pts = tapered_hex_points();
    Mesh m = mt::make_mesh(pts, "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    const double C[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    m.AddPointData("u", linear_vector(pts, C));

    const GradientResult r = gradient(m, opts("u"));
    const NDArray& g = r.mMesh.CellData("u:gradient", 0);
    ASSERT_EQ(g.Shape().size(), 2u);
    EXPECT_EQ(g.Shape()[1], 9u);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            EXPECT_NEAR(meshioplusplus::detail::read_double(g, i * 3 + j), C[i * 3 + j], kTight)
                << "d(u_" << i << ")/d(x_" << j << ")";
}

TEST(Gradient, ComponentSelectionYieldsThreeComponents) {
    const auto pts = tapered_hex_points();
    Mesh m = mt::make_mesh(pts, "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    const double C[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    m.AddPointData("u", linear_vector(pts, C));

    GradientOptions o = opts("u");
    o.mComponent = 1;
    const GradientResult r = gradient(m, o);
    expect_grad(r.mMesh, "u:gradient", 0, 0, 4.0, 5.0, 6.0);
}

// --- structural invariants --------------------------------------------------

TEST(Gradient, AnInvertedCellMatchesItsPositivelyWoundTwin) {
    // Numerator and denominator flip together, so the answer must be numerically
    // equal -- not merely finite, which would be inert.
    const auto pts = warped_hex_points();
    Mesh good = mt::make_mesh(pts, "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    good.AddPointData("f", linear_scalar(pts, 3.0, -2.0, 5.0, 7.0));
    Mesh flipped = mt::make_mesh(pts, "hexahedron", {{4, 5, 6, 7, 0, 1, 2, 3}});
    flipped.AddPointData("f", linear_scalar(pts, 3.0, -2.0, 5.0, 7.0));

    const GradientResult a = gradient(good, opts("f"));
    const GradientResult b = gradient(flipped, opts("f"));
    for (std::size_t k = 0; k < 3; ++k)
        EXPECT_NEAR(meshioplusplus::detail::read_double(a.mMesh.CellData("f:gradient", 0), k),
                    meshioplusplus::detail::read_double(b.mMesh.CellData("f:gradient", 0), k),
                    kTight);
}

TEST(Gradient, ReversingAQuadRingLeavesTheGradientUnchanged) {
    // A triangle's reversal is nearly symmetric and would be a weak test; a quad
    // exercises the ring order properly.
    const auto pts = trapezoid_points();
    Mesh fwd = mt::make_mesh(pts, "quad", {{0, 1, 2, 3}});
    fwd.AddPointData("f", linear_scalar(pts, 1.25, -0.75, 0.0, 9.0));
    Mesh rev = mt::make_mesh(pts, "quad", {{0, 3, 2, 1}});
    rev.AddPointData("f", linear_scalar(pts, 1.25, -0.75, 0.0, 9.0));

    const GradientResult a = gradient(fwd, opts("f"));
    const GradientResult b = gradient(rev, opts("f"));
    for (std::size_t k = 0; k < 3; ++k)
        EXPECT_NEAR(meshioplusplus::detail::read_double(a.mMesh.CellData("f:gradient", 0), k),
                    meshioplusplus::detail::read_double(b.mMesh.CellData("f:gradient", 0), k),
                    kTight);
}

TEST(Gradient, GeometryAndExistingDataAreUntouched) {
    Mesh m = mt::data_mesh();
    const GradientResult r = gradient(m, opts("T"));
    mt::expect_same_geometry(m, r.mMesh);
    EXPECT_TRUE(r.mMesh.HasPointData("T"));
    EXPECT_TRUE(r.mMesh.HasPointData("v"));
    EXPECT_TRUE(r.mMesh.HasCellData("mat"));
    EXPECT_TRUE(r.mMesh.HasFieldData("meta"));
    ASSERT_EQ(r.mMesh.CellDataNumBlocks("T:gradient"), r.mMesh.NumCellBlocks());
}

TEST(Gradient, DeterministicAcrossRuns) {
    Mesh m = hex_grid(3, 3, 3);
    std::vector<std::vector<double>> pts;
    for (std::size_t k = 0; k <= 3; ++k)
        for (std::size_t j = 0; j <= 3; ++j)
            for (std::size_t i = 0; i <= 3; ++i)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    m.AddPointData("f", linear_scalar(pts, 0.5, -1.25, 2.0, 4.0));

    for (GradientMethod method : {GradientMethod::GreenGauss, GradientMethod::LeastSquares}) {
        const GradientResult a = gradient(m, opts("f", method));
        const GradientResult b = gradient(m, opts("f", method));
        const NDArray& ga = a.mMesh.CellData("f:gradient", 0);
        const NDArray& gb = b.mMesh.CellData("f:gradient", 0);
        ASSERT_EQ(ga.Size(), gb.Size());
        for (std::size_t i = 0; i < ga.Size(); ++i)
            EXPECT_EQ(ga.As<double>()[i], gb.As<double>()[i]) << "run-to-run drift at " << i;
    }
}

// --- the gates must actually fire -------------------------------------------

TEST(Gradient, LeastSquaresFallsBackOnADegenerateNeighbourhood) {
    // Asserting mNumFallback == 0 on a nice mesh would be inert -- it passes if
    // the conditioning test never fires. Two genuinely rank-deficient cases:
    // a lone cell (no neighbours at all, M == 0) and a collinear strip (M has
    // rank 1 in 3D). Both must fall back AND agree with Green-Gauss.
    const auto pts = tapered_hex_points();
    Mesh lone = mt::make_mesh(pts, "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    lone.AddPointData("f", linear_scalar(pts, 3.0, -2.0, 5.0, 7.0));

    const GradientResult lsq = gradient(lone, opts("f", GradientMethod::LeastSquares));
    EXPECT_EQ(lsq.mNumFallback, 1) << "a cell with no neighbours cannot be fitted";
    EXPECT_EQ(lsq.mNumSkipped, 0) << "the fallback produces a real answer, not NaN";
    const GradientResult gg = gradient(lone, opts("f", GradientMethod::GreenGauss));
    for (std::size_t k = 0; k < 3; ++k)
        EXPECT_EQ(meshioplusplus::detail::read_double(lsq.mMesh.CellData("f:gradient", 0), k),
                  meshioplusplus::detail::read_double(gg.mMesh.CellData("f:gradient", 0), k))
            << "the fallback must be exactly the Green-Gauss answer";

    Mesh strip = hex_grid(4, 1, 1);
    std::vector<std::vector<double>> sp;
    for (std::size_t k = 0; k <= 1; ++k)
        for (std::size_t j = 0; j <= 1; ++j)
            for (std::size_t i = 0; i <= 4; ++i)
                sp.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    strip.AddPointData("f", linear_scalar(sp, 1.0, 0.0, 0.0, 0.0));
    const GradientResult sr = gradient(strip, opts("f", GradientMethod::LeastSquares));
    EXPECT_GT(sr.mNumFallback, 0) << "a collinear strip has a rank-1 normal matrix";
}

TEST(Gradient, UnsupportedCellsAreNanAndCounted) {
    // A boundary triangle block on a tet mesh is below the mesh's dimension and
    // must be NaN-and-counted, never approximated.
    Mesh m;
    m.AssignPoints(mt::points_from(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}));
    m.AddCellBlock("tetra", mt::conn_from({{0, 1, 2, 3}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}, {0, 1, 3}}));
    m.AddPointData("f", mt::data_array({0.0, 1.0, 2.0, 3.0}));

    const GradientResult r = gradient(m, opts("f"));
    EXPECT_EQ(r.mNumSkipped, 2) << "both triangles are below the mesh's dimension";
    const NDArray& tri = r.mMesh.CellData("f:gradient", 1);
    for (std::size_t i = 0; i < tri.Size(); ++i)
        EXPECT_TRUE(std::isnan(meshioplusplus::detail::read_double(tri, i)));
    const NDArray& tet = r.mMesh.CellData("f:gradient", 0);
    for (std::size_t i = 0; i < tet.Size(); ++i)
        EXPECT_FALSE(std::isnan(meshioplusplus::detail::read_double(tet, i)));
}

TEST(Gradient, RaggedBlocksAreSkippedNotGuessed) {
    Mesh m;
    m.AssignPoints(mt::points_from(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {2.0, 0.5, 0.0}}));
    m.AddPolygonBlock("polygon", {{0, 1, 2, 3}, {1, 4, 2}});
    m.AddPointData("f", mt::data_array({0.0, 1.0, 2.0, 3.0, 4.0}));

    const GradientResult r = gradient(m, opts("f"));
    EXPECT_EQ(r.mNumSkipped, 2);
    const NDArray& a = r.mMesh.CellData("f:gradient", 0);
    for (std::size_t i = 0; i < a.Size(); ++i)
        EXPECT_TRUE(std::isnan(meshioplusplus::detail::read_double(a, i)));
}

// --- point location ---------------------------------------------------------

TEST(Gradient, PointLocationIsExactUpToRounding) {
    // Deliberately a tolerance, not equality: every contributing cell value is
    // the exact constant, but summing n copies and dividing by n is not
    // bit-exact in IEEE.
    Mesh m = hex_grid(2, 2, 2);
    std::vector<std::vector<double>> pts;
    for (std::size_t k = 0; k <= 2; ++k)
        for (std::size_t j = 0; j <= 2; ++j)
            for (std::size_t i = 0; i <= 2; ++i)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    m.AddPointData("f", linear_scalar(pts, 3.0, -2.0, 5.0, 7.0));

    GradientOptions o = opts("f");
    o.mLocation = DataLocation::Point;
    const GradientResult r = gradient(m, o);
    ASSERT_TRUE(r.mMesh.HasPointData("f:gradient"));
    EXPECT_FALSE(r.mMesh.HasCellData("f:gradient")) << "the intermediate cell array is dropped";
    const NDArray& g = r.mMesh.PointData("f:gradient");
    ASSERT_EQ(g.Shape()[0], m.NumPoints());
    for (std::size_t p = 0; p < m.NumPoints(); ++p) {
        EXPECT_NEAR(meshioplusplus::detail::read_double(g, p * 3 + 0), 3.0, 1e-12);
        EXPECT_NEAR(meshioplusplus::detail::read_double(g, p * 3 + 1), -2.0, 1e-12);
        EXPECT_NEAR(meshioplusplus::detail::read_double(g, p * 3 + 2), 5.0, 1e-12);
    }
}

// --- errors -----------------------------------------------------------------

TEST(Gradient, CellDataFieldThrowsByName) {
    Mesh m = mt::data_mesh();
    try {
        gradient(m, opts("mat"));
        FAIL() << "a cell_data field has no derivative and must be rejected";
    } catch (const std::invalid_argument& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("to-point"), std::string::npos) << msg;
        EXPECT_NE(msg.find("cell_data_to_point_data"), std::string::npos) << msg;
    }
}

TEST(Gradient, BadArgumentsThrow) {
    Mesh m = mt::data_mesh();
    EXPECT_THROW(gradient(m, opts("")), std::invalid_argument);
    EXPECT_THROW(gradient(m, opts("nope")), std::invalid_argument);

    GradientOptions o = opts("v");
    o.mComponent = 7;
    EXPECT_THROW(gradient(m, o), std::invalid_argument);

    GradientOptions div = opts("v", GradientMethod::GreenGauss, GradientOperator::Divergence);
    div.mComponent = 0;
    EXPECT_THROW(gradient(m, div), std::invalid_argument)
        << "a component selection is meaningless for divergence";

    // A scalar cannot have a divergence.
    EXPECT_THROW(
        gradient(m, opts("T", GradientMethod::GreenGauss, GradientOperator::Divergence)),
        std::invalid_argument);

    GradientOptions field = opts("T");
    field.mLocation = DataLocation::Field;
    EXPECT_THROW(gradient(m, field), std::invalid_argument);
}

TEST(Gradient, OutputNameCollisionIsRejectedUnlessOverwriting) {
    Mesh m = mt::data_mesh();
    GradientOptions o = opts("T");
    o.mOutputName = "mat";
    EXPECT_THROW(gradient(m, o), std::invalid_argument);
    o.mOverwrite = true;
    EXPECT_NO_THROW(gradient(m, o));
}

TEST(Gradient, UnknownArrayMessageListsWhatExists) {
    Mesh m = mt::data_mesh();
    try {
        gradient(m, opts("nope"));
        FAIL() << "expected a throw";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("T"), std::string::npos) << e.what();
    }
}

TEST(Gradient, OperatorAndMethodParsers) {
    EXPECT_EQ(gradient_operator_from_name("gradient"), GradientOperator::Gradient);
    EXPECT_EQ(gradient_operator_from_name("grad"), GradientOperator::Gradient);
    EXPECT_EQ(gradient_operator_from_name(""), GradientOperator::Gradient);
    EXPECT_EQ(gradient_operator_from_name("divergence"), GradientOperator::Divergence);
    EXPECT_EQ(gradient_operator_from_name("div"), GradientOperator::Divergence);
    EXPECT_EQ(gradient_operator_from_name("curl"), GradientOperator::Curl);
    EXPECT_THROW(gradient_operator_from_name("laplacian"), std::invalid_argument);

    EXPECT_EQ(gradient_method_from_name("green-gauss"), GradientMethod::GreenGauss);
    EXPECT_EQ(gradient_method_from_name(""), GradientMethod::GreenGauss);
    EXPECT_EQ(gradient_method_from_name("least-squares"), GradientMethod::LeastSquares);
    EXPECT_EQ(gradient_method_from_name("lsq"), GradientMethod::LeastSquares);
    EXPECT_THROW(gradient_method_from_name("magic"), std::invalid_argument);
}

// --- vector-calculus identities ---------------------------------------------

TEST(Gradient, CurlOfAGradientAndDivergenceOfACurlVanish) {
    // Discretization-level identities, so a loose tolerance is honest here.
    Mesh m = hex_grid(3, 3, 3);
    std::vector<std::vector<double>> pts;
    for (std::size_t k = 0; k <= 3; ++k)
        for (std::size_t j = 0; j <= 3; ++j)
            for (std::size_t i = 0; i <= 3; ++i)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    m.AddPointData("f", linear_scalar(pts, 2.0, -3.0, 4.0, 0.0));
    const double C[9] = {0.0, 0.0, 7.0, 11.0, 0.0, 0.0, 0.0, 13.0, 0.0};
    m.AddPointData("u", linear_vector(pts, C));

    // curl(grad f): grad f is constant here, so its curl is identically zero.
    GradientOptions g = opts("f");
    g.mLocation = DataLocation::Point;
    Mesh with_grad = gradient(m, g).mMesh;
    const GradientResult cg =
        gradient(with_grad, opts("f:gradient", GradientMethod::GreenGauss, GradientOperator::Curl));
    const NDArray& cga = cg.mMesh.CellData("f:gradient:curl", 0);
    for (std::size_t i = 0; i < cga.Size(); ++i)
        EXPECT_NEAR(meshioplusplus::detail::read_double(cga, i), 0.0, 1e-9);

    // div(curl u).
    GradientOptions c = opts("u", GradientMethod::GreenGauss, GradientOperator::Curl);
    c.mLocation = DataLocation::Point;
    Mesh with_curl = gradient(m, c).mMesh;
    const GradientResult dc = gradient(
        with_curl, opts("u:curl", GradientMethod::GreenGauss, GradientOperator::Divergence));
    const NDArray& dca = dc.mMesh.CellData("u:curl:divergence", 0);
    for (std::size_t i = 0; i < dca.Size(); ++i)
        EXPECT_NEAR(meshioplusplus::detail::read_double(dca, i), 0.0, 1e-9);
}
