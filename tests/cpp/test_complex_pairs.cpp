//  License: MIT License (SPDX-License-Identifier: MIT)
//  Tests for roadmap §8.1 complex-field convention (`<base>_real` / `<base>_imag`).

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/formats/nastran_op2.hpp"
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_calc.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/transform.hpp"

using namespace meshioplusplus;

namespace {

Mesh make_complex_mesh() {
    Mesh m;
    NDArray pts(DType::Float64, {4, 3});
    double* p = pts.As<double>();
    const double coords[12] = {
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    };
    for (int i = 0; i < 12; ++i)
        p[i] = coords[i];
    m.AssignPoints(std::move(pts));

    NDArray conn(DType::Int64, {2, 3});
    auto* c = conn.As<std::int64_t>();
    c[0] = 0;
    c[1] = 1;
    c[2] = 2;
    c[3] = 0;
    c[4] = 2;
    c[5] = 3;
    m.AddCellBlock("triangle", std::move(conn));

    NDArray p_re(DType::Float64, {4});
    NDArray p_im(DType::Float64, {4});
    for (int i = 0; i < 4; ++i) {
        p_re.As<double>()[i] = 3.0 * (i + 1);
        p_im.As<double>()[i] = 4.0 * (i + 1);
    }
    m.AddPointData("p_real", std::move(p_re));
    m.AddPointData("p_imag", std::move(p_im));

    NDArray u_re(DType::Float64, {4, 3});
    NDArray u_im(DType::Float64, {4, 3});
    for (int i = 0; i < 4; ++i) {
        u_re.As<double>()[3 * i + 0] = 1.0;
        u_re.As<double>()[3 * i + 1] = 0.0;
        u_re.As<double>()[3 * i + 2] = 0.0;
        u_im.As<double>()[3 * i + 0] = 0.0;
        u_im.As<double>()[3 * i + 1] = 2.0;
        u_im.As<double>()[3 * i + 2] = 0.0;
    }
    m.AddPointData("U_real", std::move(u_re));
    m.AddPointData("U_imag", std::move(u_im));

    NDArray c_re(DType::Float64, {2});
    NDArray c_im(DType::Float64, {2});
    c_re.As<double>()[0] = 1.5;
    c_re.As<double>()[1] = -2.5;
    c_im.As<double>()[0] = 0.5;
    c_im.As<double>()[1] = 3.5;
    m.AddCellData("flux_real", {std::move(c_re)});
    m.AddCellData("flux_imag", {std::move(c_im)});

    return m;
}

std::string read_file_text(const std::filesystem::path& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

}  // namespace

TEST(ComplexPairs, UnvRoundTripComplexDataset) {
    const auto dir = std::filesystem::temp_directory_path() / "mio_test_complex_unv";
    std::filesystem::create_directories(dir);
    const auto path = dir / "complex.unv";

    Mesh m = make_complex_mesh();
    write_unv(path.string(), m);

    const std::string text = read_file_text(path);
    // Dataset 2414 with data_type=6 (double precision complex) should be written.
    EXPECT_NE(text.find("         6"), std::string::npos);

    Mesh back = read_unv(path.string());
    ASSERT_TRUE(back.HasPointData("p_real"));
    ASSERT_TRUE(back.HasPointData("p_imag"));
    ASSERT_TRUE(back.HasPointData("U_real"));
    ASSERT_TRUE(back.HasPointData("U_imag"));
    ASSERT_TRUE(back.HasCellData("flux_real"));
    ASSERT_TRUE(back.HasCellData("flux_imag"));

    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_NEAR(back.PointData("p_real").As<double>()[i], 3.0 * (i + 1), 1e-9);
        EXPECT_NEAR(back.PointData("p_imag").As<double>()[i], 4.0 * (i + 1), 1e-9);
    }
    EXPECT_NEAR(back.CellData("flux_real", 0).As<double>()[1], -2.5, 1e-9);
    EXPECT_NEAR(back.CellData("flux_imag", 0).As<double>()[1], 3.5, 1e-9);

    std::filesystem::remove_all(dir);
}

TEST(ComplexPairs, EnsightRoundTripComplexVariables) {
    const auto dir = std::filesystem::temp_directory_path() / "mio_test_complex_ensight";
    std::filesystem::create_directories(dir);
    const auto path = dir / "complex.case";

    Mesh m = make_complex_mesh();
    write_ensight(path.string(), m, /*binary=*/false);

    const std::string case_text = read_file_text(path);
    EXPECT_NE(case_text.find("complex scalar per node: p p_real.scl p_imag.scl"),
              std::string::npos);
    EXPECT_NE(case_text.find("complex vector per node: U U_real.vec U_imag.vec"),
              std::string::npos);
    EXPECT_NE(case_text.find("complex scalar per element: flux flux_real.escl flux_imag.escl"),
              std::string::npos);

    Mesh back = read_ensight(path.string());
    ASSERT_TRUE(back.HasPointData("p_real"));
    ASSERT_TRUE(back.HasPointData("p_imag"));
    ASSERT_TRUE(back.HasPointData("U_real"));
    ASSERT_TRUE(back.HasPointData("U_imag"));
    ASSERT_TRUE(back.HasCellData("flux_real"));
    ASSERT_TRUE(back.HasCellData("flux_imag"));
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_NEAR(back.PointData("p_real").As<double>()[i], 3.0 * (i + 1), 1e-5);
        EXPECT_NEAR(back.PointData("p_imag").As<double>()[i], 4.0 * (i + 1), 1e-5);
    }
    EXPECT_NEAR(back.CellData("flux_imag", 0).As<double>()[1], 3.5, 1e-5);

    std::filesystem::remove_all(dir);
}

TEST(ComplexPairs, DataCalcAbsAndPhase) {
    Mesh m = make_complex_mesh();

    DataCalcOptions mag_opts;
    mag_opts.output = "p_mag";
    Mesh m1 = data_calc(m, "abs(p)", mag_opts);
    ASSERT_TRUE(m1.HasPointData("p_mag"));
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_NEAR(m1.PointData("p_mag").As<double>()[i], 5.0 * (i + 1), 1e-12);

    DataCalcOptions phase_opts;
    phase_opts.output = "p_phase";
    Mesh m2 = data_calc(m, "phase(p)", phase_opts);
    ASSERT_TRUE(m2.HasPointData("p_phase"));
    const double expected_phase = std::atan2(4.0, 3.0);
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_NEAR(m2.PointData("p_phase").As<double>()[i], expected_phase, 1e-12);
}

TEST(ComplexPairs, InterpolateBaseNameExpandsPair) {
    Mesh src = make_complex_mesh();
    Mesh dst;
    dst.AssignPoints(detail::data_owned_copy(src.Points()));
    dst.AddCellBlock("triangle", detail::data_owned_copy(src.Cells(0).Conn()));

    InterpolateOptions opts;
    opts.mArrays = {"p"};
    Mesh out = interpolate(src, dst, opts);
    ASSERT_TRUE(out.HasPointData("p_real"));
    ASSERT_TRUE(out.HasPointData("p_imag"));
    EXPECT_NEAR(out.PointData("p_real").As<double>()[2], 9.0, 1e-12);
    EXPECT_NEAR(out.PointData("p_imag").As<double>()[2], 12.0, 1e-12);
}

TEST(ComplexPairs, TransformRotatesBothVectorHalves) {
    Mesh m = make_complex_mesh();
    const double half_pi = std::atan2(1.0, 0.0);
    const AffineTransform rot = transform_rotation(0.0, 0.0, 1.0, half_pi);
    Mesh out = transform(m, rot, /*rotate_vector_data=*/true);
    // U_real = (1, 0, 0) rotated +90 deg around Z -> (0, 1, 0)
    // U_imag = (0, 2, 0) rotated +90 deg around Z -> (-2, 0, 0)
    EXPECT_NEAR(out.PointData("U_real").As<double>()[0], 0.0, 1e-12);
    EXPECT_NEAR(out.PointData("U_real").As<double>()[1], 1.0, 1e-12);
    EXPECT_NEAR(out.PointData("U_imag").As<double>()[0], -2.0, 1e-12);
    EXPECT_NEAR(out.PointData("U_imag").As<double>()[1], 0.0, 1e-12);
}
