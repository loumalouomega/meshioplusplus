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
