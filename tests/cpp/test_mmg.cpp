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
// The MMG operations (roadmap 3.2) called directly: with the backend the
// probes run, without it the refusal by name does. metric_from_hessian needs no
// MMG. tests/python/test_mmg.py covers the same through the binding.

// System includes
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/mmg.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

/// A 5x5-point unit square of 32 triangles with its four side regions.
Mesh mmg_test_square() {
    return meshioplusplus::rectangle({{4, 4}}, {{1.0, 1.0}}, {{0.0, 0.0, 0.0}}, "triangle");
}

double mmg_test_mean_edge(const Mesh& rMesh) {
    double total = 0.0;
    std::size_t count = 0;
    const NDArray& p = rMesh.Points();
    const std::size_t d = rMesh.PointDim();
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.Type() != "triangle")
            continue;
        const std::int64_t* c = cb.Conn().As<std::int64_t>();
        for (std::size_t i = 0; i < cb.NumCells(); ++i)
            for (int e = 0; e < 3; ++e) {
                const std::int64_t a = c[i * 3 + e], z = c[i * 3 + (e + 1) % 3];
                double s = 0.0;
                for (std::size_t k = 0; k < d; ++k) {
                    const double x = p.As<double>()[a * d + k] - p.As<double>()[z * d + k];
                    s += x * x;
                }
                total += std::sqrt(s);
                ++count;
            }
    }
    return total / static_cast<double>(count);
}

}  // namespace

TEST(Mmg, AdaptRefusesByNameWithoutTheBackend) {
    if (meshioplusplus::mmg_available())
        GTEST_SKIP() << "built with MMG";
    meshioplusplus::AdaptOptions o;
    o.mHsiz = 0.1;
    try {
        meshioplusplus::adapt(mmg_test_square(), o);
        FAIL() << "adapt ran without MMG";
    } catch (const meshioplusplus::Unsupported& e) {
        EXPECT_NE(std::string(e.what()).find("MESHIOPLUSPLUS_WITH_MMG"), std::string::npos);
    }
}

TEST(Mmg, AdaptToAConstantSize) {
    if (!meshioplusplus::mmg_available())
        GTEST_SKIP() << "built without MMG";
    meshioplusplus::AdaptOptions o;
    o.mHsiz = 0.1;
    const meshioplusplus::AdaptResult r = meshioplusplus::adapt(mmg_test_square(), o);
    EXPECT_EQ(r.mReport.mLibrary, "mmg2d");
    EXPECT_EQ(r.mReport.mStatus, meshioplusplus::MmgStatus::Success);
    EXPECT_NEAR(mmg_test_mean_edge(r.mMesh), 0.1, 0.02);
    EXPECT_EQ(r.mMesh.PointDim(), 3u);  // the input's z comes back
    // The four side regions, under their names.
    EXPECT_EQ(r.mMesh.NumRegions(), 4u);
    for (std::size_t i = 0; i < r.mMesh.NumRegions(); ++i)
        EXPECT_GE(r.mMesh.Region(i).NumEntries(), 8u) << r.mMesh.Region(i).mName;
}

TEST(Mmg, AdaptKeepsFrozenPointsExactly) {
    if (!meshioplusplus::mmg_available())
        GTEST_SKIP() << "built without MMG";
    const Mesh square = mmg_test_square();
    meshioplusplus::AdaptOptions o;
    o.mHsiz = 0.1;
    o.mParameters.mFrozen.assign(square.NumPoints(), 0);
    o.mParameters.mFrozen[12] = 1;  // the centre
    const Mesh out = meshioplusplus::adapt(square, o).mMesh;
    const double* p = out.Points().As<double>();
    const double* q = square.Points().As<double>() + 12 * 3;
    bool found = false;
    for (std::size_t i = 0; i < out.NumPoints() && !found; ++i)
        found = p[i * 3] == q[0] && p[i * 3 + 1] == q[1] && p[i * 3 + 2] == q[2];
    EXPECT_TRUE(found);
}

TEST(Mmg, AdaptRefusesWhatMmgDoesNotTake) {
    if (!meshioplusplus::mmg_available())
        GTEST_SKIP() << "built without MMG";
    meshioplusplus::AdaptOptions o;
    EXPECT_THROW(meshioplusplus::adapt(meshioplusplus::box(), o), std::invalid_argument);
    o.mHsiz = 0.1;
    o.mSize = "nope";
    EXPECT_THROW(meshioplusplus::adapt(mmg_test_square(), o), std::invalid_argument);
}

TEST(Mmg, LevelsetAndDomainRefuseByNameWithoutTheBackend) {
    if (meshioplusplus::mmg_available())
        GTEST_SKIP() << "built with MMG";
    meshioplusplus::DiscretizeLevelsetOptions lo;
    lo.mLevelSet = "phi";
    EXPECT_THROW(meshioplusplus::discretize_levelset(mmg_test_square(), lo),
                 meshioplusplus::Unsupported);
    EXPECT_THROW(meshioplusplus::mesh_domain(mmg_test_square(), {}), meshioplusplus::Unsupported);
}

TEST(Mmg, DiscretizeLevelsetCutsADisc) {
    if (!meshioplusplus::mmg_available())
        GTEST_SKIP() << "built without MMG";
    Mesh square =
        meshioplusplus::rectangle({{10, 10}}, {{2.0, 2.0}}, {{-1.0, -1.0, 0.0}}, "triangle");
    NDArray phi = NDArray::Uninit(DType::Float64, {square.NumPoints()});
    const double* p = square.Points().As<double>();
    for (std::size_t i = 0; i < square.NumPoints(); ++i)
        phi.As<double>()[i] = std::hypot(p[i * 3], p[i * 3 + 1]) - 0.5;
    square.AddPointData("phi", std::move(phi));
    meshioplusplus::DiscretizeLevelsetOptions o;
    o.mLevelSet = "phi";
    const meshioplusplus::DiscretizeLevelsetResult r =
        meshioplusplus::discretize_levelset(square, o);
    EXPECT_EQ(r.mReport.mLibrary, "mmg2d");
    EXPECT_NE(r.mMesh.FindRegion("inside", meshioplusplus::RegionKind::Cell), Mesh::npos);
    EXPECT_NE(r.mMesh.FindRegion("outside", meshioplusplus::RegionKind::Cell), Mesh::npos);
    const std::size_t iface = r.mMesh.FindRegion("interface", meshioplusplus::RegionKind::Side);
    ASSERT_NE(iface, Mesh::npos);
    EXPECT_GT(r.mMesh.Region(iface).NumEntries(), 8u);
    // The inside's area is the disc's, pi r^2, to the polygon's error.
    const std::size_t in = r.mMesh.FindRegion("inside", meshioplusplus::RegionKind::Cell);
    const meshioplusplus::Region& inside = r.mMesh.Region(in);
    const auto cb = r.mMesh.Cells(0);
    ASSERT_EQ(cb.Type(), "triangle");
    const std::int64_t* c = cb.Conn().As<std::int64_t>();
    const double* q = r.mMesh.Points().As<double>();
    const std::size_t d = r.mMesh.PointDim();
    double area = 0.0;
    for (std::size_t k = 0; k < inside.NumEntries(); ++k) {
        const std::int64_t t = inside.Entries()[k];
        const double* a = q + c[t * 3] * d;
        const double* b = q + c[t * 3 + 1] * d;
        const double* e = q + c[t * 3 + 2] * d;
        area += 0.5 * std::abs((b[0] - a[0]) * (e[1] - a[1]) - (b[1] - a[1]) * (e[0] - a[0]));
    }
    EXPECT_NEAR(area, M_PI * 0.25, 0.02);
}

TEST(Mmg, MeshDomainKeepsTheHoleOut) {
    if (!meshioplusplus::mmg_available())
        GTEST_SKIP() << "built without MMG";
    // Two squares, [0, 1]^2 around [0.4, 0.6]^2.
    Mesh boundary;
    const double xy[16] = {0, 0, 1, 0, 1, 1, 0, 1, 0.4, 0.4, 0.6, 0.4, 0.6, 0.6, 0.4, 0.6};
    NDArray points = NDArray::Uninit(DType::Float64, {8, 2});
    std::copy(xy, xy + 16, points.As<double>());
    boundary.AssignPoints(std::move(points));
    NDArray conn = NDArray::Uninit(DType::Int64, {8, 2});
    const std::int64_t lines[16] = {0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4};
    std::copy(lines, lines + 16, conn.As<std::int64_t>());
    boundary.AddCellBlock("line", std::move(conn));
    meshioplusplus::MeshDomainOptions o;
    o.mHsiz = 0.05;
    const meshioplusplus::MeshDomainResult both = meshioplusplus::mesh_domain(boundary, o);
    EXPECT_EQ(both.mNumSubdomains, 2);
    EXPECT_EQ(both.mNumHoles, 0);
    o.mHoles = {0.5, 0.5};
    const meshioplusplus::MeshDomainResult r = meshioplusplus::mesh_domain(boundary, o);
    EXPECT_EQ(r.mReport.mLibrary, "mmg2d");
    EXPECT_EQ(r.mNumSubdomains, 1);
    EXPECT_EQ(r.mNumHoles, 1);
    // No point strictly inside the hole.
    const double* q = r.mMesh.Points().As<double>();
    const std::size_t d = r.mMesh.PointDim();
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        const bool in_hole = q[i * d] > 0.4 + 1e-9 && q[i * d] < 0.6 - 1e-9 &&
                             q[i * d + 1] > 0.4 + 1e-9 && q[i * d + 1] < 0.6 - 1e-9;
        EXPECT_FALSE(in_hole) << i;
    }
    o.mHoles = {5.0, 5.0};
    EXPECT_THROW(meshioplusplus::mesh_domain(boundary, o), std::invalid_argument);
}

TEST(Mmg, MoveRefusesByNameWithoutTheBackendOrElas) {
    if (meshioplusplus::mmg_available() && meshioplusplus::mmg_has_elas())
        GTEST_SKIP() << "built with MMG and ELAS";
    meshioplusplus::MoveOptions o;
    o.mDisplacement = "d";
    o.mMoving = {"xmin"};
    try {
        meshioplusplus::move(mmg_test_square(), o);
        FAIL() << "move ran without MMG or ELAS";
    } catch (const meshioplusplus::Unsupported& e) {
        const std::string what = e.what();
        EXPECT_NE(what.find(meshioplusplus::mmg_available() ? "ELAS" : "MESHIOPLUSPLUS_WITH_MMG"),
                  std::string::npos)
            << what;
    }
}

TEST(Mmg, MoveTranslatesAHole) {
    if (!meshioplusplus::mmg_available() || !meshioplusplus::mmg_has_elas())
        GTEST_SKIP() << "built without MMG or ELAS";
    // [0, 1]^2 around the hole [0.4, 0.6]^2, the hole's lines a cell region.
    Mesh boundary;
    const double xy[16] = {0, 0, 1, 0, 1, 1, 0, 1, 0.4, 0.4, 0.6, 0.4, 0.6, 0.6, 0.4, 0.6};
    NDArray points = NDArray::Uninit(DType::Float64, {8, 2});
    std::copy(xy, xy + 16, points.As<double>());
    boundary.AssignPoints(std::move(points));
    NDArray conn = NDArray::Uninit(DType::Int64, {8, 2});
    const std::int64_t lines[16] = {0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4};
    std::copy(lines, lines + 16, conn.As<std::int64_t>());
    boundary.AddCellBlock("line", std::move(conn));
    NDArray hole = NDArray::Uninit(DType::Int64, {4});
    for (std::int64_t k = 0; k < 4; ++k)
        hole.As<std::int64_t>()[k] = 4 + k;
    boundary.AddRegion(
        meshioplusplus::Region("hole", meshioplusplus::RegionKind::Cell, -1, 0, std::move(hole)));
    meshioplusplus::MeshDomainOptions mo;
    mo.mHsiz = 0.05;
    mo.mHoles = {0.5, 0.5};
    Mesh domain = meshioplusplus::mesh_domain(boundary, mo).mMesh;

    // Translate the hole's points by +0.1 in x.
    const std::size_t np = domain.NumPoints();
    const std::size_t d = domain.PointDim();
    const double* p = domain.Points().As<double>();
    NDArray disp = NDArray::Uninit(DType::Float64, {np, 2});
    for (std::size_t i = 0; i < np; ++i) {
        const double x = p[i * d], y = p[i * d + 1];
        const bool on_hole = x > 0.39 && x < 0.61 && y > 0.39 && y < 0.61;
        disp.As<double>()[i * 2] = on_hole ? 0.1 : 0.0;
        disp.As<double>()[i * 2 + 1] = 0.0;
    }
    domain.AddPointData("d", std::move(disp));
    meshioplusplus::MoveOptions o;
    o.mDisplacement = "d";
    o.mMoving = {"hole"};
    const meshioplusplus::AdaptResult r = meshioplusplus::move(domain, o);
    EXPECT_EQ(r.mReport.mLibrary, "mmg2d");
    EXPECT_EQ(r.mReport.mStatus, meshioplusplus::MmgStatus::Success);
    // Every triangle keeps its orientation, and the hole region comes back.
    const Mesh& out = r.mMesh;
    const double* q = out.Points().As<double>();
    const std::size_t e = out.PointDim();
    for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
        const auto cb = out.Cells(b);
        const std::int64_t* c = cb.Conn().As<std::int64_t>();
        for (std::size_t i = 0; i < cb.NumCells(); ++i) {
            if (cb.Type() == "triangle") {
                const double* a = q + c[i * 3] * e;
                const double* bb = q + c[i * 3 + 1] * e;
                const double* cc = q + c[i * 3 + 2] * e;
                EXPECT_GT((bb[0] - a[0]) * (cc[1] - a[1]) - (bb[1] - a[1]) * (cc[0] - a[0]), 0.0);
            }
        }
    }
    const std::size_t hole_region = out.FindRegion("hole", meshioplusplus::RegionKind::Cell);
    ASSERT_NE(hole_region, Mesh::npos);
    EXPECT_GT(out.Region(hole_region).NumEntries(), 0u);

    o.mMoving = {"nope"};
    EXPECT_THROW(meshioplusplus::move(domain, o), std::invalid_argument);
}

TEST(Mmg, MetricFromHessianScalesAndClamps) {
    Mesh mesh = mt::tri_mesh();
    NDArray h(DType::Float64, {mesh.NumPoints(), 9});
    for (std::size_t i = 0; i < mesh.NumPoints(); ++i) {
        h.As<double>()[i * 9 + 0] = 2.0;
        h.As<double>()[i * 9 + 4] = 20.0;
    }
    mesh.AddPointData("h", std::move(h));
    meshioplusplus::MetricFromHessianOptions o;
    o.mHessian = "h";
    o.mError = 0.01;
    o.mHmax = 10.0;
    const Mesh out = meshioplusplus::metric_from_hessian(mesh, o);
    const double* m = out.PointData("metric").As<double>();
    const double c = 0.5 * 0.75 * 0.75;
    EXPECT_NEAR(m[0], c * 2.0 / 0.01, 1e-9);
    EXPECT_NEAR(m[1], 0.0, 1e-12);
    EXPECT_NEAR(m[2], c * 20.0 / 0.01, 1e-9);
    EXPECT_NEAR(m[5], 1.0 / 100.0, 1e-12);  // a zero curvature clamps to 1/hmax^2
    o.mError = 0.0;
    EXPECT_THROW(meshioplusplus::metric_from_hessian(mesh, o), std::invalid_argument);
}
