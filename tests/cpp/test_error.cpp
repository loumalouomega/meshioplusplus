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
// Tests for the ZZ recovery-based error estimator (`operations/error.hpp`).
//
// The central oracle: `gradient`'s Green-Gauss is EXACT for a linear field, so
// the raw and recovered gradients agree everywhere and `error:zz` must be
// (up to rounding) zero mesh-wide -- "the estimator is zero when the field it
// is estimating for is exactly resolvable" is the one property every
// recovery-based estimator must have. A quadratic field breaks that exactness
// on purpose, giving a nonzero, non-uniform indicator to mark against.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/error.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::error_marking_from_name;
using meshioplusplus::error_method_from_name;
using meshioplusplus::ErrorMarking;
using meshioplusplus::ErrorMethod;
using meshioplusplus::ErrorOptions;
using meshioplusplus::ErrorResult;
using meshioplusplus::estimate_error;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::detail::lattice_build_mesh;
using meshioplusplus::detail::LatticeSpec;

constexpr double kTight = 1e-9;

// A 3x3x3 hexahedron grid over [0, 3]^3 -- enough cells (27) that a marking
// policy has a real ranking to work with, small enough to stay fast.
Mesh grid_mesh() {
    LatticeSpec spec;
    spec.mOrigin = {{0.0, 0.0, 0.0}};
    spec.mSpacing = {{1.0, 1.0, 1.0}};
    spec.mDims = {{3, 3, 3}};
    return lattice_build_mesh(spec);
}

/// Attaches `point_data["f"] = fn(x, y, z)` to a copy of `rMesh`.
Mesh with_field(Mesh mesh, const std::string& name,
                const std::function<double(double, double, double)>& fn) {
    const NDArray& pts = mesh.Points();
    const std::size_t n = mesh.NumPoints();
    const std::size_t dim = mesh.PointDim();
    NDArray field(DType::Float64, {n});
    double* pf = field.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        double x = meshioplusplus::detail::read_double(pts, i * dim + 0);
        double y = dim > 1 ? meshioplusplus::detail::read_double(pts, i * dim + 1) : 0.0;
        double z = dim > 2 ? meshioplusplus::detail::read_double(pts, i * dim + 2) : 0.0;
        pf[i] = fn(x, y, z);
    }
    mesh.AddPointData(name, std::move(field));
    return mesh;
}

ErrorOptions opts(const std::string& array, ErrorMarking marking = ErrorMarking::None,
                  double marking_value = 0.0) {
    ErrorOptions o;
    o.mArrayName = array;
    o.mMarking = marking;
    o.mMarkingValue = marking_value;
    return o;
}

/// Every entry of a per-block cell_data array, in block-major order.
std::vector<double> flat_double(const Mesh& rMesh, const std::string& name) {
    std::vector<double> out;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const NDArray& a = rMesh.CellData(name, b);
        const double* p = a.As<double>();
        for (std::size_t i = 0; i < meshioplusplus::detail::rows(a); ++i)
            out.push_back(p[i]);
    }
    return out;
}

std::vector<std::int64_t> flat_int(const Mesh& rMesh, const std::string& name) {
    std::vector<std::int64_t> out;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const NDArray& a = rMesh.CellData(name, b);
        const std::int64_t* p = a.As<std::int64_t>();
        for (std::size_t i = 0; i < meshioplusplus::detail::rows(a); ++i)
            out.push_back(p[i]);
    }
    return out;
}

// --- name parsing ------------------------------------------------------------

TEST(ErrorNames, MethodRoundTrips) {
    EXPECT_EQ(error_method_from_name("zz"), ErrorMethod::Zz);
    EXPECT_EQ(error_method_from_name(""), ErrorMethod::Zz);
    EXPECT_THROW(error_method_from_name("kelly"), std::invalid_argument);
}

TEST(ErrorNames, MarkingRoundTrips) {
    EXPECT_EQ(error_marking_from_name(""), ErrorMarking::None);
    EXPECT_EQ(error_marking_from_name("none"), ErrorMarking::None);
    EXPECT_EQ(error_marking_from_name("absolute"), ErrorMarking::Absolute);
    EXPECT_EQ(error_marking_from_name("abs"), ErrorMarking::Absolute);
    EXPECT_EQ(error_marking_from_name("fraction"), ErrorMarking::Fraction);
    EXPECT_EQ(error_marking_from_name("frac"), ErrorMarking::Fraction);
    EXPECT_EQ(error_marking_from_name("dorfler"), ErrorMarking::Dorfler);
    EXPECT_EQ(error_marking_from_name("bulk"), ErrorMarking::Dorfler);
    EXPECT_THROW(error_marking_from_name("median"), std::invalid_argument);
}

// --- validation ---------------------------------------------------------------

TEST(EstimateError, RejectsEmptyArrayName) {
    Mesh m = with_field(grid_mesh(), "f", [](double x, double, double) { return x; });
    EXPECT_THROW(estimate_error(m, opts("")), std::invalid_argument);
}

TEST(EstimateError, RejectsUnknownArray) {
    Mesh m = with_field(grid_mesh(), "f", [](double x, double, double) { return x; });
    EXPECT_THROW(estimate_error(m, opts("nope")), std::invalid_argument);
}

TEST(EstimateError, RejectsCellDataArrayNamingTheFix) {
    Mesh m = with_field(grid_mesh(), "f", [](double x, double, double) { return x; });
    std::vector<NDArray> blocks;
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b)
        blocks.push_back(NDArray(DType::Float64, {m.Cells(b).NumCells()}));
    m.AddCellData("f_cell", std::move(blocks));
    try {
        estimate_error(m, opts("f_cell"));
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("cell_data_to_point_data"), std::string::npos);
    }
}

TEST(EstimateError, RejectsMarkingValueOutOfRangeForFraction) {
    Mesh m = with_field(grid_mesh(), "f", [](double x, double, double) { return x; });
    EXPECT_THROW(estimate_error(m, opts("f", ErrorMarking::Fraction, 0.0)), std::invalid_argument);
    EXPECT_THROW(estimate_error(m, opts("f", ErrorMarking::Fraction, 1.5)), std::invalid_argument);
    EXPECT_THROW(estimate_error(m, opts("f", ErrorMarking::Dorfler, -0.1)), std::invalid_argument);
}

TEST(EstimateError, OutputNameCollisionErrorsWithoutOverwrite) {
    Mesh m = with_field(grid_mesh(), "f", [](double x, double, double) { return x; });
    ErrorResult first = estimate_error(m, opts("f"));
    EXPECT_THROW(estimate_error(first.mMesh, opts("f")), std::invalid_argument);
    ErrorOptions o = opts("f");
    o.mOverwrite = true;
    EXPECT_NO_THROW(estimate_error(first.mMesh, o));
}

// --- the exactness oracle -----------------------------------------------------

TEST(EstimateError, LinearFieldGivesZeroIndicatorEverywhere) {
    Mesh m = with_field(grid_mesh(), "f", [](double x, double y, double z) {
        return 2.0 * x - 3.0 * y + 0.5 * z + 7.0;
    });
    ErrorResult r = estimate_error(m, opts("f"));
    EXPECT_EQ(r.mNumSkipped, 0);
    EXPECT_NEAR(r.mGlobalError, 0.0, kTight);
    for (double v : flat_double(r.mMesh, "error:zz"))
        EXPECT_NEAR(v, 0.0, kTight);
}

TEST(EstimateError, LinearFieldMarksNothingUnderAnyPolicy) {
    Mesh m = with_field(grid_mesh(), "f", [](double x, double y, double z) { return x + y + z; });
    ErrorResult abs_r = estimate_error(m, opts("f", ErrorMarking::Absolute, 1e-9));
    EXPECT_EQ(abs_r.mNumMarked, 0);
    ErrorResult frac_r = estimate_error(m, opts("f", ErrorMarking::Fraction, 0.5));
    // Every indicator is (numerically) zero, so Fraction has no ranking signal
    // to prefer one zero cell over another -- it may still fill its quota. What
    // must hold is the quota itself.
    EXPECT_EQ(static_cast<std::size_t>(frac_r.mNumMarked),
              static_cast<std::size_t>(0.5 * 27 + 0.5));
    ErrorResult dorfler_r = estimate_error(m, opts("f", ErrorMarking::Dorfler, 0.5));
    EXPECT_EQ(dorfler_r.mNumMarked, 0);  // total squared error is 0, target is 0
}

// --- a genuine, non-uniform indicator ------------------------------------------

TEST(EstimateError, QuadraticFieldGivesNonzeroNonuniformIndicator) {
    Mesh m = with_field(grid_mesh(), "f",
                        [](double x, double y, double z) { return x * x + y * y + z * z; });
    ErrorResult r = estimate_error(m, opts("f"));
    EXPECT_EQ(r.mNumSkipped, 0);
    EXPECT_GT(r.mGlobalError, 0.0);
    const std::vector<double> zz = flat_double(r.mMesh, "error:zz");
    EXPECT_EQ(zz.size(), 27u);
    bool any_positive = false;
    double minv = zz[0], maxv = zz[0];
    for (double v : zz) {
        EXPECT_GE(v, 0.0);
        any_positive = any_positive || v > kTight;
        minv = std::min(minv, v);
        maxv = std::max(maxv, v);
    }
    EXPECT_TRUE(any_positive);
    // A quadratic field on a uniform-but-not-centred grid does not give every
    // cell the same curvature-relative indicator -- corner cells differ from
    // face-centre cells. If this ever failed it would mean the recovery step
    // degenerated into a no-op.
    EXPECT_GT(maxv - minv, kTight);

    // mGlobalError is exactly sqrt(sum zz_i^2) -- cross-check independently.
    double sumsq = 0.0;
    for (double v : zz)
        sumsq += v * v;
    EXPECT_NEAR(r.mGlobalError, std::sqrt(sumsq), 1e-9 * std::sqrt(sumsq));
}

// --- marking policies -----------------------------------------------------------

TEST(EstimateError, AbsoluteMarksExactlyCellsAboveThreshold) {
    Mesh m = with_field(grid_mesh(), "f",
                        [](double x, double y, double z) { return x * x + y * y + z * z; });
    ErrorResult plain = estimate_error(m, opts("f"));
    std::vector<double> zz = flat_double(plain.mMesh, "error:zz");
    std::sort(zz.begin(), zz.end());
    const double threshold = zz[zz.size() / 2];  // a value with cells on both sides

    ErrorResult marked = estimate_error(m, opts("f", ErrorMarking::Absolute, threshold));
    const std::vector<double> zz2 = flat_double(marked.mMesh, "error:zz");
    const std::vector<std::int64_t> mk = flat_int(marked.mMesh, "error:marked");
    ASSERT_EQ(zz2.size(), mk.size());
    std::int64_t expected = 0;
    for (std::size_t i = 0; i < zz2.size(); ++i) {
        const bool should = zz2[i] > threshold;
        EXPECT_EQ(mk[i], should ? 1 : 0) << "cell " << i;
        expected += should ? 1 : 0;
    }
    EXPECT_EQ(marked.mNumMarked, expected);
}

TEST(EstimateError, FractionMarksTheLargestIndicatorsAndTheRightCount) {
    Mesh m = with_field(grid_mesh(), "f",
                        [](double x, double y, double z) { return x * x + y * y + z * z; });
    const double fraction = 0.25;
    ErrorResult r = estimate_error(m, opts("f", ErrorMarking::Fraction, fraction));
    const std::vector<double> zz = flat_double(r.mMesh, "error:zz");
    const std::vector<std::int64_t> mk = flat_int(r.mMesh, "error:marked");

    const auto expected_count =
        static_cast<std::size_t>(fraction * static_cast<double>(zz.size()) + 0.5);
    EXPECT_EQ(static_cast<std::size_t>(r.mNumMarked), expected_count);

    // Every marked indicator must be >= every unmarked one (a valid top-k, not
    // merely the right count).
    double min_marked = 1e300;
    double max_unmarked = -1e300;
    for (std::size_t i = 0; i < zz.size(); ++i) {
        if (mk[i] != 0)
            min_marked = std::min(min_marked, zz[i]);
        else
            max_unmarked = std::max(max_unmarked, zz[i]);
    }
    if (expected_count > 0 && expected_count < zz.size())
        EXPECT_GE(min_marked, max_unmarked);
}

TEST(EstimateError, DorflerMarksTheMinimalBulkPrefix) {
    Mesh m = with_field(grid_mesh(), "f",
                        [](double x, double y, double z) { return x * x + y * y + z * z; });
    const double theta = 0.6;
    ErrorResult r = estimate_error(m, opts("f", ErrorMarking::Dorfler, theta));
    const std::vector<double> zz = flat_double(r.mMesh, "error:zz");
    const std::vector<std::int64_t> mk = flat_int(r.mMesh, "error:marked");

    double total_sq = 0.0;
    for (double v : zz)
        total_sq += v * v;
    double marked_sq = 0.0;
    std::int64_t marked_count = 0;
    for (std::size_t i = 0; i < zz.size(); ++i) {
        if (mk[i] != 0) {
            marked_sq += zz[i] * zz[i];
            ++marked_count;
        }
    }
    EXPECT_EQ(r.mNumMarked, marked_count);
    // Sufficiency: the marked set covers at least theta of the total.
    EXPECT_GE(marked_sq, theta * total_sq - 1e-9);
    // Minimality: dropping the single smallest marked indicator would no
    // longer cover the target.
    double min_marked = 1e300;
    for (std::size_t i = 0; i < zz.size(); ++i)
        if (mk[i] != 0)
            min_marked = std::min(min_marked, zz[i]);
    if (marked_count > 0)
        EXPECT_LT(marked_sq - min_marked * min_marked, theta * total_sq - 1e-9);
}

TEST(EstimateError, MarkingNoneAttachesNoMarkedArray) {
    Mesh m = with_field(grid_mesh(), "f",
                        [](double x, double y, double z) { return x * x + y * y + z * z; });
    ErrorResult r = estimate_error(m, opts("f"));
    EXPECT_FALSE(r.mMesh.HasCellData("error:marked"));
    EXPECT_EQ(r.mNumMarked, 0);
}

// --- output naming --------------------------------------------------------------

TEST(EstimateError, CustomOutputNames) {
    Mesh m = with_field(grid_mesh(), "f",
                        [](double x, double y, double z) { return x * x + y * y + z * z; });
    ErrorOptions o = opts("f", ErrorMarking::Absolute, 0.5);
    o.mOutputName = "my_indicator";
    o.mMarkedName = "my_marks";
    ErrorResult r = estimate_error(m, o);
    EXPECT_TRUE(r.mMesh.HasCellData("my_indicator"));
    EXPECT_TRUE(r.mMesh.HasCellData("my_marks"));
    EXPECT_FALSE(r.mMesh.HasCellData("error:zz"));
    EXPECT_FALSE(r.mMesh.HasCellData("error:marked"));
}

// --- geometry / other data pass through untouched --------------------------------

TEST(EstimateError, GeometryAndUnrelatedDataPassThroughUnchanged) {
    Mesh m = with_field(grid_mesh(), "f",
                        [](double x, double y, double z) { return x * x + y * y + z * z; });
    NDArray tag(DType::Int64, {m.Cells(0).NumCells()});
    for (std::size_t i = 0; i < m.Cells(0).NumCells(); ++i)
        tag.As<std::int64_t>()[i] = 42;
    m.AddCellData("tag", {tag});

    ErrorResult r = estimate_error(m, opts("f"));
    EXPECT_EQ(r.mMesh.NumPoints(), m.NumPoints());
    EXPECT_EQ(r.mMesh.NumCellBlocks(), m.NumCellBlocks());
    ASSERT_TRUE(r.mMesh.HasCellData("tag"));
    const NDArray& out_tag = r.mMesh.CellData("tag", 0);
    for (std::size_t i = 0; i < meshioplusplus::detail::rows(out_tag); ++i)
        EXPECT_EQ(out_tag.As<std::int64_t>()[i], 42);
    ASSERT_TRUE(r.mMesh.HasPointData("f"));
}

}  // namespace
