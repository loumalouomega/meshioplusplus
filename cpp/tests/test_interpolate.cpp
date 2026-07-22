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

// System includes
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/interpolate.hpp"

namespace {

using meshioplusplus::interpolate;
using meshioplusplus::interpolate_conflict_from_name;
using meshioplusplus::interpolate_method_from_name;
using meshioplusplus::InterpolateConflict;
using meshioplusplus::InterpolateMethod;
using meshioplusplus::InterpolateOptions;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
namespace md = meshioplusplus::detail;

// --- fixtures ---------------------------------------------------------------

// An n x n x n block of unit-cube hexahedra over [shift, n + shift]^3.
Mesh hex_grid(int n, double shift = 0.0) {
    const int m = n + 1;
    std::vector<std::vector<double>> pts;
    for (int k = 0; k < m; ++k)
        for (int j = 0; j < m; ++j)
            for (int i = 0; i < m; ++i)
                pts.push_back({static_cast<double>(i) + shift, static_cast<double>(j) + shift,
                               static_cast<double>(k) + shift});
    auto pid = [m](int i, int j, int k) { return static_cast<std::int64_t>((k * m + j) * m + i); };
    std::vector<std::vector<std::int64_t>> cells;
    for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                cells.push_back({pid(i, j, k), pid(i + 1, j, k), pid(i + 1, j + 1, k),
                                 pid(i, j + 1, k), pid(i, j, k + 1), pid(i + 1, j, k + 1),
                                 pid(i + 1, j + 1, k + 1), pid(i, j + 1, k + 1)});
    return mt::make_mesh(std::move(pts), "hexahedron", std::move(cells));
}

// An n x n grid of unit quads in the z = 0 plane over [0, n]^2.
Mesh quad_grid_2d(int n) {
    const int m = n + 1;
    std::vector<std::vector<double>> pts;
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < m; ++i)
            pts.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});
    auto pid = [m](int i, int j) { return static_cast<std::int64_t>(j * m + i); };
    std::vector<std::vector<std::int64_t>> cells;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            cells.push_back({pid(i, j), pid(i + 1, j), pid(i + 1, j + 1), pid(i, j + 1)});
    return mt::make_mesh(std::move(pts), "quad", std::move(cells));
}

// A point-cloud target (points only, no cell blocks).
Mesh point_cloud(std::vector<std::vector<double>> pts) {
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    return m;
}

// The linear field 2x + 3y + 5z + 7 evaluated at every point of `rMesh`.
NDArray linear_field(const Mesh& rMesh) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    std::vector<double> vals(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = md::read_double(rMesh.Points(), i * dim + 0);
        const double y = dim > 1 ? md::read_double(rMesh.Points(), i * dim + 1) : 0.0;
        const double z = dim > 2 ? md::read_double(rMesh.Points(), i * dim + 2) : 0.0;
        vals[i] = 2.0 * x + 3.0 * y + 5.0 * z + 7.0;
    }
    return mt::data_array(vals);
}

double linear_at(const std::vector<double>& rP) {
    const double y = rP.size() > 1 ? rP[1] : 0.0;
    const double z = rP.size() > 2 ? rP[2] : 0.0;
    return 2.0 * rP[0] + 3.0 * y + 5.0 * z + 7.0;
}

// --- nearest ----------------------------------------------------------------

TEST(Interpolate, NearestPiecewiseConstant) {
    Mesh src = point_cloud({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}});
    src.AddPointData("f", mt::data_array({10.0, 20.0, 30.0, 40.0}));
    Mesh tgt = point_cloud({{0.1, 0.1, 0}, {0.9, 0.2, 0}, {0.2, 0.8, 0}, {1.4, 1.3, 0}});

    const Mesh out = interpolate(src, tgt);
    const NDArray& f = out.PointData("f");
    EXPECT_DOUBLE_EQ(md::read_double(f, 0), 10.0);
    EXPECT_DOUBLE_EQ(md::read_double(f, 1), 20.0);
    EXPECT_DOUBLE_EQ(md::read_double(f, 2), 30.0);
    EXPECT_DOUBLE_EQ(md::read_double(f, 3), 40.0);
}

TEST(Interpolate, NearestTieBreaksToLowestIndex) {
    Mesh src = point_cloud({{0, 0, 0}, {1, 0, 0}});
    src.AddPointData("f", mt::data_array({10.0, 20.0}));
    Mesh tgt = point_cloud({{0.5, 0, 0}});

    const Mesh out = interpolate(src, tgt);
    EXPECT_DOUBLE_EQ(md::read_double(out.PointData("f"), 0), 10.0);
}

TEST(Interpolate, NearestPreservesIntegerDtype) {
    Mesh src = point_cloud({{0, 0, 0}, {1, 0, 0}});
    src.AddPointData("tag", mt::int_data_array({3, 9}));
    Mesh tgt = point_cloud({{0.1, 0, 0}, {0.9, 0, 0}});

    const Mesh out = interpolate(src, tgt);
    const NDArray& tag = out.PointData("tag");
    EXPECT_FALSE(md::is_float_dtype(tag.Dtype()));
    EXPECT_EQ(md::read_int(tag, 0), 3);
    EXPECT_EQ(md::read_int(tag, 1), 9);
}

// --- barycentric ------------------------------------------------------------

TEST(Interpolate, BarycentricLinearFieldExact) {
    Mesh src = hex_grid(2);
    src.AddPointData("f", linear_field(src));
    const std::vector<std::vector<double>> probes = {
        {0.25, 0.5, 0.75}, {1.1, 0.2, 1.9}, {0.0, 0.0, 0.0}, {2.0, 2.0, 2.0},
        {1.0, 1.0, 1.0},   {1.5, 0.5, 0.5}, {0.7, 1.3, 0.2}};
    Mesh tgt = point_cloud(probes);

    InterpolateOptions opts;
    opts.mMethod = InterpolateMethod::Barycentric;
    const Mesh out = interpolate(src, tgt, opts);
    const NDArray& f = out.PointData("f");
    EXPECT_EQ(f.Dtype(), meshioplusplus::DType::Float64);
    for (std::size_t i = 0; i < probes.size(); ++i)
        EXPECT_NEAR(md::read_double(f, i), linear_at(probes[i]), 1e-12);
}

TEST(Interpolate, BarycentricTriangleSource2D) {
    Mesh src = quad_grid_2d(3);
    src.AddPointData("f", linear_field(src));
    const std::vector<std::vector<double>> probes = {
        {0.5, 0.5, 0}, {2.7, 1.2, 0}, {0.0, 3.0, 0}, {1.0, 1.0, 0}};
    Mesh tgt = point_cloud(probes);

    InterpolateOptions opts;
    opts.mMethod = InterpolateMethod::Barycentric;
    const Mesh out = interpolate(src, tgt, opts);
    for (std::size_t i = 0; i < probes.size(); ++i)
        EXPECT_NEAR(md::read_double(out.PointData("f"), i), linear_at(probes[i]), 1e-12);
}

TEST(Interpolate, BarycentricVectorComponentsExact) {
    Mesh src = hex_grid(2);
    const std::size_t n = src.NumPoints();
    std::vector<double> v(2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = md::read_double(src.Points(), i * 3 + 0);
        const double y = md::read_double(src.Points(), i * 3 + 1);
        v[2 * i + 0] = x + 2.0 * y;
        v[2 * i + 1] = 3.0 * x - y;
    }
    src.AddPointData("v", mt::data_array(v, 2));
    Mesh tgt = point_cloud({{0.3, 1.4, 0.6}});

    InterpolateOptions opts;
    opts.mMethod = InterpolateMethod::Barycentric;
    const Mesh out = interpolate(src, tgt, opts);
    const NDArray& ov = out.PointData("v");
    ASSERT_EQ(ov.Shape().size(), 2u);
    EXPECT_NEAR(md::read_double(ov, 0), 0.3 + 2.0 * 1.4, 1e-12);
    EXPECT_NEAR(md::read_double(ov, 1), 3.0 * 0.3 - 1.4, 1e-12);
}

TEST(Interpolate, OutsideDefaultValue) {
    Mesh src = hex_grid(1);
    src.AddPointData("f", linear_field(src));
    Mesh tgt = point_cloud({{0.5, 0.5, 0.5}, {5.0, 5.0, 5.0}});

    InterpolateOptions opts;
    opts.mMethod = InterpolateMethod::Barycentric;
    opts.mDefaultValue = -123.0;
    const Mesh out = interpolate(src, tgt, opts);
    EXPECT_NEAR(md::read_double(out.PointData("f"), 0), linear_at({0.5, 0.5, 0.5}), 1e-12);
    EXPECT_DOUBLE_EQ(md::read_double(out.PointData("f"), 1), -123.0);
}

TEST(Interpolate, OutsideExtrapolateFallsBackToNearest) {
    Mesh src = hex_grid(1);
    src.AddPointData("f", linear_field(src));
    Mesh tgt = point_cloud({{5.0, 5.0, 5.0}});

    InterpolateOptions opts;
    opts.mMethod = InterpolateMethod::Barycentric;
    opts.mExtrapolate = true;
    const Mesh out = interpolate(src, tgt, opts);
    // Nearest source point to (5,5,5) is the cube corner (1,1,1).
    EXPECT_DOUBLE_EQ(md::read_double(out.PointData("f"), 0), linear_at({1.0, 1.0, 1.0}));
}

TEST(Interpolate, BarycentricIntegerDataBecomesFloat64) {
    Mesh src = hex_grid(1);
    src.AddPointData("tag", mt::int_data_array({1, 2, 3, 4, 5, 6, 7, 8}));
    Mesh tgt = point_cloud({{0.5, 0.5, 0.5}});

    InterpolateOptions opts;
    opts.mMethod = InterpolateMethod::Barycentric;
    const Mesh out = interpolate(src, tgt, opts);
    EXPECT_EQ(out.PointData("tag").Dtype(), meshioplusplus::DType::Float64);
}

TEST(Interpolate, BarycentricNoSimplicesThrows) {
    Mesh src = point_cloud({{0, 0, 0}, {1, 0, 0}});  // no cells at all
    src.AddPointData("f", mt::data_array({1.0, 2.0}));
    Mesh tgt = point_cloud({{0.5, 0, 0}});

    InterpolateOptions opts;
    opts.mMethod = InterpolateMethod::Barycentric;
    EXPECT_THROW(interpolate(src, tgt, opts), std::invalid_argument);
}

// --- cell_data --------------------------------------------------------------

TEST(Interpolate, CellDataAlwaysNearestCell) {
    // Two unit quads side by side, tagged 1 and 2.
    Mesh src = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}},
                             "quad", {{0, 1, 4, 3}, {1, 2, 5, 4}});
    src.AddCellData("mat", {mt::data_array({1.0, 2.0})});
    // One quad centred over the second source quad.
    Mesh tgt = mt::make_mesh({{1.2, 0.2, 0}, {1.8, 0.2, 0}, {1.8, 0.8, 0}, {1.2, 0.8, 0}}, "quad",
                             {{0, 1, 2, 3}});

    for (const InterpolateMethod method :
         {InterpolateMethod::Nearest, InterpolateMethod::Barycentric}) {
        InterpolateOptions opts;
        opts.mMethod = method;
        opts.mArrays = {"mat"};
        const Mesh out = interpolate(src, tgt, opts);
        ASSERT_EQ(out.CellDataNumBlocks("mat"), 1u);
        EXPECT_DOUBLE_EQ(md::read_double(out.CellData("mat", 0), 0), 2.0);
    }
}

TEST(Interpolate, DefaultArraysIsPointDataOnly) {
    Mesh src = quad_grid_2d(1);
    src.AddPointData("f", linear_field(src));
    src.AddCellData("mat", {mt::data_array({4.0})});
    Mesh tgt = quad_grid_2d(1);

    const Mesh out = interpolate(src, tgt);
    EXPECT_TRUE(out.HasPointData("f"));
    EXPECT_FALSE(out.HasCellData("mat"));
}

TEST(Interpolate, ExplicitArrayInBothLocationsTransfersBoth) {
    Mesh src = quad_grid_2d(1);
    src.AddPointData("m", linear_field(src));
    src.AddCellData("m", {mt::data_array({4.0})});
    Mesh tgt = quad_grid_2d(1);

    InterpolateOptions opts;
    opts.mArrays = {"m"};
    const Mesh out = interpolate(src, tgt, opts);
    EXPECT_TRUE(out.HasPointData("m"));
    EXPECT_TRUE(out.HasCellData("m"));
}

// --- conflicts / validation -------------------------------------------------

TEST(Interpolate, ConflictErrorOverwriteSuffix) {
    Mesh src = point_cloud({{0, 0, 0}});
    src.AddPointData("f", mt::data_array({99.0}));
    Mesh tgt = point_cloud({{0.25, 0, 0}});
    tgt.AddPointData("f", mt::data_array({1.0}));

    EXPECT_THROW(interpolate(src, tgt), std::invalid_argument);

    InterpolateOptions ow;
    ow.mOnConflict = InterpolateConflict::Overwrite;
    EXPECT_DOUBLE_EQ(md::read_double(interpolate(src, tgt, ow).PointData("f"), 0), 99.0);

    InterpolateOptions sf;
    sf.mOnConflict = InterpolateConflict::Suffix;
    const Mesh out = interpolate(src, tgt, sf);
    EXPECT_DOUBLE_EQ(md::read_double(out.PointData("f"), 0), 1.0);  // target's own survives
    EXPECT_DOUBLE_EQ(md::read_double(out.PointData("f_interp"), 0), 99.0);

    tgt.AddPointData("f_interp", mt::data_array({2.0}));
    EXPECT_THROW(interpolate(src, tgt, sf), std::invalid_argument);
}

TEST(Interpolate, TargetGeometryAndOwnDataUntouched) {
    Mesh src = quad_grid_2d(2);
    src.AddPointData("srcfield", linear_field(src));
    Mesh tgt = mt::data_mesh();

    const Mesh out = interpolate(src, tgt);
    mt::expect_same_geometry(tgt, out);
    EXPECT_TRUE(out.HasPointData("srcfield"));
    for (const std::string& name : tgt.PointDataNames()) {
        const NDArray& a = tgt.PointData(name);
        const NDArray& b = out.PointData(name);
        ASSERT_EQ(a.Size(), b.Size());
        for (std::size_t i = 0; i < a.Size(); ++i)
            EXPECT_DOUBLE_EQ(md::read_double(a, i), md::read_double(b, i));
    }
    for (const std::string& name : tgt.CellDataNames()) {
        ASSERT_EQ(tgt.CellDataNumBlocks(name), out.CellDataNumBlocks(name));
        for (std::size_t b = 0; b < tgt.CellDataNumBlocks(name); ++b)
            for (std::size_t i = 0; i < tgt.CellData(name, b).Size(); ++i)
                EXPECT_DOUBLE_EQ(md::read_double(tgt.CellData(name, b), i),
                                 md::read_double(out.CellData(name, b), i));
    }
    for (const std::string& name : tgt.FieldDataNames())
        EXPECT_TRUE(out.HasFieldData(name));
}

TEST(Interpolate, UnknownArrayThrows) {
    Mesh src = point_cloud({{0, 0, 0}});
    src.AddPointData("f", mt::data_array({1.0}));
    Mesh tgt = point_cloud({{0, 0, 0}});
    InterpolateOptions opts;
    opts.mArrays = {"nope"};
    EXPECT_THROW(interpolate(src, tgt, opts), std::invalid_argument);
}

TEST(Interpolate, EmptySourceThrows) {
    Mesh src;
    Mesh tgt = point_cloud({{0, 0, 0}});
    EXPECT_THROW(interpolate(src, tgt), std::invalid_argument);
}

TEST(Interpolate, MethodAndConflictFromName) {
    EXPECT_EQ(interpolate_method_from_name("nearest"), InterpolateMethod::Nearest);
    EXPECT_EQ(interpolate_method_from_name("barycentric"), InterpolateMethod::Barycentric);
    EXPECT_THROW(interpolate_method_from_name("bogus"), std::invalid_argument);
    EXPECT_EQ(interpolate_conflict_from_name("error"), InterpolateConflict::Error);
    EXPECT_EQ(interpolate_conflict_from_name("overwrite"), InterpolateConflict::Overwrite);
    EXPECT_EQ(interpolate_conflict_from_name("suffix"), InterpolateConflict::Suffix);
    EXPECT_THROW(interpolate_conflict_from_name("bogus"), std::invalid_argument);
}

// --- determinism ------------------------------------------------------------

TEST(Interpolate, RepeatRunsByteIdentical) {
    Mesh src = hex_grid(3);
    src.AddPointData("f", linear_field(src));
    src.AddCellData("mat", {mt::data_array(std::vector<double>(27, 5.0))});
    // Shifted off-grid so samples land strictly inside cells.
    Mesh tgt = hex_grid(2, 0.37);

    for (const InterpolateMethod method :
         {InterpolateMethod::Nearest, InterpolateMethod::Barycentric}) {
        InterpolateOptions opts;
        opts.mMethod = method;
        opts.mArrays = {"f", "mat"};
        const Mesh a = interpolate(src, tgt, opts);
        const Mesh b = interpolate(src, tgt, opts);
        const NDArray& fa = a.PointData("f");
        const NDArray& fb = b.PointData("f");
        ASSERT_EQ(fa.Nbytes(), fb.Nbytes());
        EXPECT_EQ(0, std::memcmp(fa.Data(), fb.Data(), fa.Nbytes()));
        const NDArray& ma = a.CellData("mat", 0);
        const NDArray& mb = b.CellData("mat", 0);
        ASSERT_EQ(ma.Nbytes(), mb.Nbytes());
        EXPECT_EQ(0, std::memcmp(ma.Data(), mb.Data(), ma.Nbytes()));
    }
}

}  // namespace
