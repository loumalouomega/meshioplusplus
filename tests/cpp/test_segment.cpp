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
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"
#include "analytic_surfaces.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

std::vector<std::int64_t> patches(const Mesh& rMesh, std::size_t Block = 0) {
    const NDArray& a = rMesh.CellData(kSegmentPatchName, Block);
    return std::vector<std::int64_t>(a.As<std::int64_t>(), a.As<std::int64_t>() + a.Size());
}

std::vector<std::int64_t> adjacency(const Mesh& rMesh) {
    const NDArray& a = rMesh.FieldData(kSegmentAdjacencyName);
    return std::vector<std::int64_t>(a.As<std::int64_t>(), a.As<std::int64_t>() + a.Size());
}

// A cylinder closed by two flat quad-and-triangle caps sharing its rims.
Mesh capped_cylinder() {
    Mesh m = analytic::cylinder(1.0, 2.0, 4);
    const std::size_t n = 32;
    std::vector<std::vector<double>> pts;
    const NDArray& p = m.Points();
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        pts.push_back(
            {p.As<double>()[i * 3], p.As<double>()[i * 3 + 1], p.As<double>()[i * 3 + 2]});
    const std::int64_t bottom = static_cast<std::int64_t>(pts.size());
    pts.push_back({0.0, 0.0, 0.0});
    const std::int64_t top = bottom + 1;
    pts.push_back({0.0, 0.0, 2.0});
    std::vector<std::vector<std::int64_t>> quads;
    const NDArray& c = m.Cells(0).Conn();
    for (std::size_t i = 0; i < m.Cells(0).NumCells(); ++i)
        quads.push_back({c.As<std::int64_t>()[i * 4], c.As<std::int64_t>()[i * 4 + 1],
                         c.As<std::int64_t>()[i * 4 + 2], c.As<std::int64_t>()[i * 4 + 3]});
    std::vector<std::vector<std::int64_t>> tris;
    const std::int64_t last = 4 * static_cast<std::int64_t>(n);
    for (std::int64_t k = 0; k < static_cast<std::int64_t>(n); ++k) {
        const std::int64_t k1 = (k + 1) % static_cast<std::int64_t>(n);
        tris.push_back({bottom, k1, k});
        tris.push_back({top, last + k, last + k1});
    }
    Mesh out;
    out.AssignPoints(mt::points_from(pts));
    out.AddCellBlock("quad", mt::conn_from(quads));
    out.AddCellBlock("triangle", mt::conn_from(tris));
    return out;
}

}  // namespace

TEST(Segment, BoxSkinGivesSixFacesAndTwelveAdjacencies) {
    const SegmentResult r = segment(analytic::box_skin(3));
    EXPECT_EQ(r.mNumPatches, 6);
    EXPECT_EQ(r.mNumMerged, 0);
    EXPECT_EQ(r.mNumBoundaryEdges, 0);
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    EXPECT_EQ(r.mNumFeatureEdges, 12 * 3);
    const std::vector<std::int64_t> adj = adjacency(r.mMesh);
    ASSERT_EQ(adj.size(), 12u * 3u);
    for (std::size_t row = 0; row < 12; ++row) {
        EXPECT_LT(adj[row * 3], adj[row * 3 + 1]);
        EXPECT_EQ(adj[row * 3 + 2], 3);  // three quad edges along each cube edge
    }
    // Each face is nine cells, and every region names one patch.
    std::vector<std::int64_t> count(6, 0);
    for (std::int64_t p : patches(r.mMesh))
        ++count[static_cast<std::size_t>(p)];
    for (std::int64_t c : count)
        EXPECT_EQ(c, 9);
    for (int k = 0; k < 6; ++k) {
        const std::size_t i = r.mMesh.FindRegion("patch:" + std::to_string(k), RegionKind::Cell);
        ASSERT_NE(i, Mesh::npos);
        EXPECT_EQ(r.mMesh.Region(i).NumEntries(), 9u);
    }
}

TEST(Segment, InputIsLeftUnchanged) {
    const Mesh in = analytic::box_skin(2);
    const SegmentResult r = segment(in);
    ASSERT_EQ(r.mMesh.NumPoints(), in.NumPoints());
    ASSERT_EQ(r.mMesh.NumCellBlocks(), in.NumCellBlocks());
    EXPECT_EQ(std::memcmp(r.mMesh.Points().Data(), in.Points().Data(), in.Points().Nbytes()), 0);
    EXPECT_EQ(std::memcmp(r.mMesh.Cells(0).Conn().Data(), in.Cells(0).Conn().Data(),
                          in.Cells(0).Conn().Nbytes()),
              0);
}

TEST(Segment, CappedCylinderGivesThreePatches) {
    const SegmentResult r = segment(capped_cylinder());
    EXPECT_EQ(r.mNumPatches, 3);
    EXPECT_EQ(r.mNumBoundaryEdges, 0);
    const std::vector<std::int64_t> side = patches(r.mMesh, 0), caps = patches(r.mMesh, 1);
    EXPECT_EQ(std::set<std::int64_t>(side.begin(), side.end()).size(), 1u);
    EXPECT_EQ(std::set<std::int64_t>(caps.begin(), caps.end()).size(), 2u);
}

TEST(Segment, NormalToleranceSplitsASmoothSurface) {
    SegmentOptions o;
    EXPECT_EQ(segment(analytic::cylinder(), o).mNumPatches, 1);
    o.mNormalToleranceDeg = 40.0;
    EXPECT_GT(segment(analytic::cylinder(), o).mNumPatches, 1);
}

TEST(Segment, CurvatureToleranceSeparatesCylinderFromPlane) {
    // A flat strip tangent to a half cylinder: no crease anywhere.
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> quads;
    const int arc = 16, flat = 6, rows = 4;
    const int cols = arc + flat + 1;
    for (int j = 0; j <= rows; ++j) {
        for (int i = 0; i <= arc; ++i) {
            const double a = analytic::kPi * i / arc;
            pts.push_back({std::cos(a), std::sin(a), double(j)});
        }
        for (int i = 1; i <= flat; ++i)
            pts.push_back({-1.0, -0.5 * i, double(j)});
    }
    for (int j = 0; j < rows; ++j)
        for (int i = 0; i + 1 < cols; ++i) {
            const std::int64_t p = j * cols + i;
            quads.push_back({p, p + 1, p + 1 + cols, p + cols});
        }
    const Mesh m = mt::make_mesh(pts, "quad", quads);
    SegmentOptions o;
    EXPECT_EQ(segment(m, o).mNumPatches, 1);
    o.mCurvatureTolerance = 0.5;
    EXPECT_EQ(segment(m, o).mNumPatches, 2);
}

TEST(Segment, SmallPatchesMergeIntoTheirNeighbour) {
    SegmentOptions o;
    o.mNormalToleranceDeg = 40.0;
    const std::int64_t before = segment(analytic::torus(), o).mNumPatches;
    o.mMinPatchCells = 40;
    const SegmentResult r = segment(analytic::torus(), o);
    EXPECT_GT(r.mNumMerged, 0);
    EXPECT_EQ(r.mNumPatches + r.mNumMerged, before);
    std::vector<std::int64_t> count(static_cast<std::size_t>(r.mNumPatches), 0);
    for (std::int64_t p : patches(r.mMesh))
        ++count[static_cast<std::size_t>(p)];
    for (std::int64_t c : count)
        EXPECT_GT(c, 0);
}

TEST(Segment, RegionRestrictsAndLinesAreUnlabelled) {
    Mesh m = analytic::box_skin(2);
    std::vector<std::int64_t> first(8);
    for (std::int64_t i = 0; i < 8; ++i)
        first[static_cast<std::size_t>(i)] = i;
    NDArray e = NDArray::Uninit(DType::Int64, {first.size()});
    std::copy(first.begin(), first.end(), e.As<std::int64_t>());
    m.AddRegion(Region("half", RegionKind::Cell, std::move(e)));
    m.AddCellBlock("line", mt::conn_from({{0, 1}}));
    SegmentOptions o;
    o.mRegion = "half";
    const SegmentResult r = segment(m, o);
    const std::vector<std::int64_t> p = patches(r.mMesh, 0);
    for (std::size_t c = 8; c < p.size(); ++c)
        EXPECT_EQ(p[c], -1);
    EXPECT_EQ(patches(r.mMesh, 1), std::vector<std::int64_t>{-1});
    EXPECT_EQ(r.mNumPatches, 2);  // the first two faces of the cube
}

TEST(Segment, RerunReplacesThePatchRegions) {
    SegmentOptions o;
    o.mNormalToleranceDeg = 40.0;
    const SegmentResult many = segment(analytic::cylinder(), o);
    const SegmentResult one = segment(many.mMesh);
    EXPECT_EQ(one.mNumPatches, 1);
    std::size_t n = 0;
    for (std::size_t i = 0; i < one.mMesh.NumRegions(); ++i)
        n += one.mMesh.Region(i).mName.rfind("patch:", 0) == 0 ? 1 : 0;
    EXPECT_EQ(n, 1u);
}

TEST(Segment, RefusesVolumesAndBadOptions) {
    EXPECT_THROW(segment(mt::tet_mesh()), std::invalid_argument);
    SegmentOptions o;
    o.mFeatureAngleDeg = 200.0;
    EXPECT_THROW(segment(analytic::box_skin(1), o), std::invalid_argument);
    o = {};
    o.mMinPatchCells = 0;
    EXPECT_THROW(segment(analytic::box_skin(1), o), std::invalid_argument);
    o = {};
    o.mRegion = "nope";
    EXPECT_THROW(segment(analytic::box_skin(1), o), std::invalid_argument);
    o = {};
    o.mPrefix = "";
    EXPECT_THROW(segment(analytic::box_skin(1), o), std::invalid_argument);
}

TEST(Segment, PolygonsAndInconsistentWinding) {
    // Two coplanar squares as polygons, the second wound the other way.
    const Mesh m = [] {
        Mesh out;
        out.AssignPoints(
            mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}}));
        out.AddPolygonBlock("polygon",
                            std::vector<std::vector<std::int64_t>>{{0, 1, 2, 3}, {1, 2, 5, 4}});
        return out;
    }();
    const SegmentResult r = segment(m);
    EXPECT_EQ(r.mNumPatches, 1);
    EXPECT_EQ(patches(r.mMesh), (std::vector<std::int64_t>{0, 0}));
}

namespace {

std::size_t distinct(const std::vector<std::int64_t>& rLabels) {
    return std::set<std::int64_t>(rLabels.begin(), rLabels.end()).size();
}

SegmentOptions primitive_method(double Tol = 1e-3) {
    SegmentOptions o;
    o.mMethod = SegmentMethod::Primitive;
    o.mTolerance = Tol;
    return o;
}

}  // namespace

TEST(Segment, PrimitiveMethodSplitsATangentFillet) {
    // The fillet meets its two sides at smooth edges: angles alone give one
    // patch for the three, the primitive method a plane, a cylinder, a plane.
    const Mesh plate = analytic::filleted_plate();
    const SegmentResult by_angle = segment(plate);
    const SegmentResult by_surface = segment(plate, primitive_method());
    EXPECT_EQ(by_surface.mNumPatches, 8);  // six planes, the bore and the fillet
    EXPECT_LT(by_angle.mNumPatches, by_surface.mNumPatches);
    EXPECT_EQ(distinct(patches(by_surface.mMesh)), 8u);
    // Every face of the input is labelled, and the output is reproducible.
    for (std::int64_t p : patches(by_surface.mMesh))
        EXPECT_GE(p, 0);
    EXPECT_EQ(patches(segment(plate, primitive_method()).mMesh), patches(by_surface.mMesh));
}

TEST(Segment, PrimitiveMethodLeavesWhatNoPrimitiveExplains) {
    // A freeform wall fits small spheres and tori anywhere; the primitive
    // method keeps it whole, as the angle method does.
    const Mesh wavy = analytic::wavy_capped();
    EXPECT_EQ(segment(wavy, primitive_method()).mNumPatches, segment(wavy).mNumPatches);
    // An analytic shape is one patch per surface either way.
    EXPECT_EQ(segment(analytic::box_skin(3), primitive_method()).mNumPatches, 6);
}

TEST(Segment, PrimitiveMethodRefusesBadOptions) {
    SegmentOptions o = primitive_method(0.0);
    EXPECT_THROW(segment(analytic::box_skin(1), o), std::invalid_argument);
    o = primitive_method(std::nan(""));
    EXPECT_THROW(segment(analytic::box_skin(1), o), std::invalid_argument);
    o = {};
    o.mMethod = static_cast<SegmentMethod>(7);
    EXPECT_THROW(segment(analytic::box_skin(1), o), std::invalid_argument);
}
