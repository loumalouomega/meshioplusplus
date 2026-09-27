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

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/feature_edges.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/region.hpp"
#include "mesh_fixtures.hpp"

// Project includes (private, not installed)
#include "../../src/cpp/src/detail/crease_edges.hpp"

using namespace meshioplusplus;

namespace {

using mt::make_mesh;

Mesh cube_quads() {
    return make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "quad",
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 4, 7, 3}, {1, 2, 6, 5}});
}

// A flat grid of quads in z = 0 -- `Half` columns either side of x = Half,
// `Rows` rows -- with the right half folded up about x = Half by `FoldDeg`.
Mesh folded_grid(double FoldDeg, int Half = 1, int Rows = 2) {
    const double a = FoldDeg * 3.14159265358979323846 / 180.0;
    const int nx = 2 * Half + 1;
    std::vector<std::vector<double>> pts;
    for (int j = 0; j <= Rows; ++j)
        for (int i = 0; i < nx; ++i) {
            if (i <= Half)
                pts.push_back({double(i), double(j), 0.0});
            else
                pts.push_back(
                    {Half + (i - Half) * std::cos(a), double(j), (i - Half) * std::sin(a)});
        }
    std::vector<std::vector<std::int64_t>> quads;
    for (int j = 0; j < Rows; ++j)
        for (int i = 0; i + 1 < nx; ++i) {
            const std::int64_t p = j * nx + i;
            quads.push_back({p, p + 1, p + 1 + nx, p + nx});
        }
    return make_mesh(pts, "quad", quads);
}

std::vector<std::int32_t> kinds(const Mesh& rMesh) {
    const NDArray& a = rMesh.CellData(kFeatureKindName, 0);
    return std::vector<std::int32_t>(a.As<std::int32_t>(), a.As<std::int32_t>() + a.Size());
}

std::size_t count_kind(const Mesh& rMesh, FeatureEdgeKind Kind) {
    std::size_t c = 0;
    for (std::int32_t k : kinds(rMesh))
        c += k == static_cast<std::int32_t>(Kind) ? 1 : 0;
    return c;
}

}  // namespace

TEST(FeatureEdges, CubeHasTwelveSharpEdges) {
    const FeatureEdgeResult r = feature_edges(cube_quads());
    EXPECT_EQ(r.mNumFeature, 12);
    EXPECT_EQ(r.mNumBoundary, 0);
    EXPECT_EQ(r.mNumNonManifold, 0);
    EXPECT_EQ(r.mNumInconsistent, 0);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "line");
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 12u);
    EXPECT_EQ(r.mMesh.NumPoints(), 8u);
    const NDArray& angle = r.mMesh.CellData(kFeatureAngleName, 0);
    for (std::size_t i = 0; i < angle.Size(); ++i)
        EXPECT_NEAR(angle.As<double>()[i], 90.0, 1e-9);
    // Ascending (low, high) order.
    const NDArray& conn = r.mMesh.Cells(0).Conn();
    for (std::size_t i = 0; i < 12; ++i) {
        const std::int64_t* e = conn.As<std::int64_t>() + i * 2;
        EXPECT_LT(e[0], e[1]);
        if (i > 0)
            EXPECT_TRUE(e[-2] < e[0] || (e[-2] == e[0] && e[-1] < e[1]));
    }
}

TEST(FeatureEdges, AngleThresholdIsStrict) {
    FeatureEdgeOptions o;
    o.mBoundary = false;
    o.mFeatureAngleDeg = 30.0;
    EXPECT_EQ(feature_edges(folded_grid(20.0), o).mNumFeature, 0);
    EXPECT_EQ(feature_edges(folded_grid(40.0), o).mNumFeature, 2);
    const Mesh out = feature_edges(folded_grid(40.0), o).mMesh;
    EXPECT_EQ(out.Cells(0).NumCells(), 2u);
    EXPECT_NEAR(out.CellData(kFeatureAngleName, 0).As<double>()[0], 40.0, 1e-9);
}

TEST(FeatureEdges, SmoothSphereHasNoFeatures) {
    FeatureEdgeOptions o;
    o.mFeatureAngleDeg = 60.0;
    const FeatureEdgeResult r = feature_edges(mt::icosphere(2, 1.0), o);
    EXPECT_EQ(r.mNumFeature, 0);
    EXPECT_EQ(r.mNumBoundary, 0);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 0u);
}

TEST(FeatureEdges, OpenSquareReportsItsBoundaryLoop) {
    const FeatureEdgeResult r = feature_edges(folded_grid(0.0));
    EXPECT_EQ(r.mNumBoundary, 8);
    EXPECT_EQ(count_kind(r.mMesh, FeatureEdgeKind::Boundary), 8u);
    FeatureEdgeOptions o;
    o.mBoundary = false;
    EXPECT_EQ(feature_edges(folded_grid(0.0), o).mMesh.Cells(0).NumCells(), 0u);
}

TEST(FeatureEdges, NonManifoldFinIsReported) {
    // Three triangles sharing edge (0, 1).
    const Mesh fin = make_mesh({{0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0.5, -1, 0}, {0.5, 0, 1}},
                               "triangle", {{0, 1, 2}, {1, 0, 3}, {0, 1, 4}});
    const FeatureEdgeResult r = feature_edges(fin);
    EXPECT_EQ(r.mNumNonManifold, 1);
    EXPECT_EQ(count_kind(r.mMesh, FeatureEdgeKind::NonManifold), 1u);
    EXPECT_TRUE(std::isnan(r.mMesh.CellData(kFeatureAngleName, 0).As<double>()[0]));
}

TEST(FeatureEdges, InconsistentPairIsReorientedBeforeTheAngle) {
    // Two coplanar triangles wound the same way along (1, 2): after the
    // reorientation they are flat, so the pair is inconsistent but not sharp.
    const Mesh m =
        make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}, {0, 2, 3}});
    const Mesh bad =
        make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}, {2, 0, 3}});
    EXPECT_EQ(feature_edges(m).mNumInconsistent, 0);
    const FeatureEdgeResult r = feature_edges(bad);
    EXPECT_EQ(r.mNumInconsistent, 1);
    EXPECT_EQ(r.mNumFeature, 0);
    EXPECT_EQ(count_kind(r.mMesh, FeatureEdgeKind::Inconsistent), 1u);
}

TEST(FeatureEdges, VolumeMeshUsesItsSkin) {
    const FeatureEdgeResult r = feature_edges(mt::hex_mesh());
    EXPECT_EQ(r.mNumFeature, 12);
    EXPECT_EQ(r.mNumBoundary, 0);
}

TEST(FeatureEdges, RegionRestrictsTheCells) {
    Mesh m = folded_grid(0.0);
    Region reg;
    reg.mName = "left";
    reg.mKind = RegionKind::Cell;
    reg.mEntries = NDArray::Uninit(DType::Int64, {1});
    reg.mEntries.As<std::int64_t>()[0] = 0;
    m.AddRegion(reg);
    FeatureEdgeOptions o;
    o.mRegion = "left";
    EXPECT_EQ(feature_edges(m, o).mNumBoundary, 4);
    o.mRegion = "nope";
    EXPECT_THROW(feature_edges(m, o), std::invalid_argument);
}

TEST(FeatureEdges, RejectsAnAngleOutOfRange) {
    FeatureEdgeOptions o;
    o.mFeatureAngleDeg = 181.0;
    EXPECT_THROW(feature_edges(cube_quads(), o), std::invalid_argument);
    o.mFeatureAngleDeg = -1.0;
    EXPECT_THROW(feature_edges(cube_quads(), o), std::invalid_argument);
}

TEST(FeatureEdges, KernelSkipsCollapsedAndSelfRevisitedEdges) {
    // A degenerate ring (0, 0, 1) has one real edge; a ring walking (0,1) twice
    // is not a surface edge.
    const std::vector<std::int64_t> start = {0, 3};
    const std::vector<std::int64_t> nodes = {0, 0, 1};
    const std::vector<double> normals = {0.0, 0.0, 0.0};
    const auto edges = detail::crease_edges(start, nodes, normals, 30.0);
    // (0,1) used by 0->1 and 1->0 of the same face: skipped.
    EXPECT_TRUE(edges.empty());
}

TEST(FeatureEdges, DecimateNoLongerPinsASmoothCoarseSphere) {
    // Around every vertex of a once-subdivided icosphere, faces on opposite
    // sides are about 40 degrees apart, which the old per-vertex test pinned
    // (so nothing could collapse); adjacent faces are about 21 degrees apart,
    // so no edge is a crease at 30 degrees and decimation now proceeds.
    DecimateOptions o;
    o.mTargetRatio = 0.5;
    o.mFeatureAngleDeg = 30.0;
    const DecimateResult r = decimate(mt::icosphere(1, 1.0), o);
    EXPECT_LT(r.mMesh.Cells(0).NumCells(), 80u);
}

TEST(FeatureEdges, SmoothKeepsTheCreasesOfABentSurface) {
    // A 6x4 grid folded 90 degrees about x = 3: the interior crease nodes
    // (x = 3, rows 1..3) are not on the boundary, so only the crease test can
    // hold them -- before v16.23.0 smooth never looked at a surface's creases.
    Mesh m = folded_grid(90.0, 3, 4);
    // Knock an interior node of the flat half off the regular lattice, so
    // there is something to smooth.
    NDArray pts = detail::data_owned_copy(m.Points());
    pts.As<double>()[(2 * 7 + 1) * 3] += 0.3;
    m.AssignPoints(std::move(pts));
    SmoothOptions o;
    o.mIterations = 5;
    const SmoothResult r = smooth(m, o);
    const NDArray& a = m.Points();
    const NDArray& b = r.mMesh.Points();
    for (std::size_t j = 1; j < 4; ++j) {
        const std::size_t node = j * 7 + 3;
        for (std::size_t d = 0; d < 3; ++d)
            EXPECT_EQ(a.As<double>()[node * 3 + d], b.As<double>()[node * 3 + d]) << node;
    }
    EXPECT_GT(r.mNumNodesMoved, 0);
}
