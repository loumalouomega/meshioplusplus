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
// Tests for the mesh quality operation: regular elements give ideal metric
// values, deliberately-inverted/degenerate elements are flagged.

// System includes
#include <cmath>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/quality.hpp"

namespace {

using meshioplusplus::compute_quality;
using meshioplusplus::Mesh;
using meshioplusplus::QualityReport;

double q_val(const QualityReport& rRep, const std::string& rName, std::size_t block,
             std::size_t cell) {
    for (const auto& e : rRep.mCellArrays)
        if (e.first == rName)
            return e.second[block].template As<double>()[cell];
    return std::nan("");
}

const double kSqrt2 = std::sqrt(2.0);

}  // namespace

TEST(Quality, UnitCubeHexIsIdeal) {
    Mesh m = mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    QualityReport r = compute_quality(m);
    EXPECT_NEAR(q_val(r, "quality:volume", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:scaled_jacobian", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:skewness", 0, 0), 0.0, 1e-12);
    EXPECT_EQ(r.mNumInverted, 0);
    EXPECT_EQ(r.mNumDegenerate, 0);
}

TEST(Quality, RegularTetraIsIdeal) {
    const double s = kSqrt2;
    Mesh m = mt::make_mesh({{0, 0, 0},
                            {s, 0, 0},
                            {s / 2, s * std::sqrt(3.0) / 2, 0},
                            {s / 2, s * std::sqrt(3.0) / 6, s * std::sqrt(2.0 / 3.0)}},
                           "tetra", {{0, 1, 2, 3}});
    QualityReport r = compute_quality(m);
    EXPECT_GT(q_val(r, "quality:volume", 0, 0), 0.0);
    EXPECT_NEAR(q_val(r, "quality:scaled_jacobian", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:aspect_ratio", 0, 0), 1.0, 1e-12);
    const double ideal_dih = std::acos(1.0 / 3.0) * 180.0 / M_PI;  // ~70.5288
    EXPECT_NEAR(q_val(r, "quality:min_dihedral", 0, 0), ideal_dih, 1e-9);
    EXPECT_NEAR(q_val(r, "quality:max_dihedral", 0, 0), ideal_dih, 1e-9);
    EXPECT_EQ(r.mNumInverted, 0);
}

TEST(Quality, EquilateralTriangle2D) {
    Mesh m = mt::make_mesh({{0, 0}, {1, 0}, {0.5, std::sqrt(3.0) / 2}}, "triangle", {{0, 1, 2}});
    QualityReport r = compute_quality(m);
    EXPECT_NEAR(q_val(r, "quality:aspect_ratio", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:min_angle", 0, 0), 60.0, 1e-9);
    EXPECT_NEAR(q_val(r, "quality:max_angle", 0, 0), 60.0, 1e-9);
}

TEST(Quality, UnitSquare2D) {
    Mesh m = mt::make_mesh({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, "quad", {{0, 1, 2, 3}});
    QualityReport r = compute_quality(m);
    EXPECT_NEAR(q_val(r, "quality:volume", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:skewness", 0, 0), 0.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:aspect_ratio", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:warpage", 0, 0), 0.0, 1e-12);
}

TEST(Quality, ReferenceWedgeAndPyramidScaledJacobian) {
    const double h = std::sqrt(3.0) / 2;
    Mesh w = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0.5, h, 0}, {0, 0, 1}, {1, 0, 1}, {0.5, h, 1}},
                           "wedge", {{0, 1, 2, 3, 4, 5}});
    EXPECT_NEAR(q_val(compute_quality(w), "quality:scaled_jacobian", 0, 0), 1.0, 1e-12);

    const double ph = 1.0 / kSqrt2;
    Mesh p = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, ph}}, "pyramid",
                           {{0, 1, 2, 3, 4}});
    EXPECT_NEAR(q_val(compute_quality(p), "quality:scaled_jacobian", 0, 0), 1.0, 1e-12);
}

TEST(Quality, InvertedTetraIsFlagged) {
    // A regular tetra with two nodes swapped -> negative signed volume.
    const double s = kSqrt2;
    Mesh m = mt::make_mesh({{0, 0, 0},
                            {s, 0, 0},
                            {s / 2, s * std::sqrt(3.0) / 2, 0},
                            {s / 2, s * std::sqrt(3.0) / 6, s * std::sqrt(2.0 / 3.0)}},
                           "tetra", {{0, 2, 1, 3}});
    QualityReport r = compute_quality(m);
    EXPECT_LT(q_val(r, "quality:volume", 0, 0), 0.0);
    EXPECT_EQ(q_val(r, "quality:inverted", 0, 0), 1.0);
    EXPECT_EQ(r.mNumInverted, 1);
}

TEST(Quality, DegenerateTriangleIsFlagged) {
    // Three collinear points -> zero area.
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}, "triangle", {{0, 1, 2}});
    QualityReport r = compute_quality(m);
    EXPECT_EQ(q_val(r, "quality:degenerate", 0, 0), 1.0);
    EXPECT_EQ(r.mNumDegenerate, 1);
    EXPECT_TRUE(std::isnan(q_val(r, "quality:aspect_ratio", 0, 0)));
}

TEST(Quality, AttachQualityAddsCellData) {
    Mesh m = mt::tet_mesh();
    Mesh out = meshioplusplus::attach_quality(m);
    EXPECT_TRUE(out.HasCellData("quality:volume"));
    EXPECT_TRUE(out.HasCellData("quality:inverted"));
    EXPECT_EQ(out.CellDataNumBlocks("quality:volume"), m.NumCellBlocks());
}

TEST(Quality, QuadraticTetraUsesCorners) {
    // tetra10: the metrics should match the underlying linear tetra corners.
    Mesh lin = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra", {{0, 1, 2, 3}});
    Mesh quad = mt::tet10_mesh();
    // Just assert tetra10 is scored (finite volume), exercising the corner path.
    EXPECT_TRUE(std::isfinite(q_val(compute_quality(quad), "quality:volume", 0, 0)));
    EXPECT_GT(q_val(compute_quality(lin), "quality:volume", 0, 0), 0.0);
}
