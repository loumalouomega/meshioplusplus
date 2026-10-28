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
// extrude and revolve (roadmap 3.1.2, 3.1.3): surfaces swept into volumes and
// curves into surfaces, checked against `box`, analytic volumes and the
// conformity check.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"
#include "meshioplusplus/operations/extrude.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/quality.hpp"

using meshioplusplus::ExtrudeOptions;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::RegionKind;
using meshioplusplus::RevolveOptions;
using meshioplusplus::SweepResult;

namespace {

using Vec3 = meshioplusplus::detail::Vec3;
constexpr double kPi = 3.14159265358979323846;

std::vector<double> ex_volumes(const Mesh& rM) {
    std::vector<double> out;
    meshioplusplus::detail::CellRings rings;
    std::vector<Vec3> xyz;
    for (const auto cb : rM.CellRange())
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            if (meshioplusplus::detail::cell_rings(cb, c, rM.Points(), rM.PointDim(), rings, xyz))
                out.push_back(meshioplusplus::detail::poly_measure(rings, xyz.data()).mVolume);
    return out;
}

double ex_total(const Mesh& rM) {
    double s = 0.0;
    for (const double v : ex_volumes(rM)) {
        EXPECT_GT(v, 0.0);
        s += v;
    }
    return s;
}

std::vector<std::int64_t> ex_entries(const Mesh& rM, const std::string& rName,
                                     RegionKind Kind = RegionKind::Side) {
    const std::size_t i = rM.FindRegion(rName, Kind);
    if (i == Mesh::npos)
        return {};
    const auto& r = rM.Region(i);
    const std::int64_t* e = r.Entries();
    return e ? std::vector<std::int64_t>(e, e + r.mEntries.Size()) : std::vector<std::int64_t>{};
}

std::vector<std::pair<std::int64_t, std::int64_t>> ex_pairs(const std::vector<std::int64_t>& rE) {
    std::vector<std::pair<std::int64_t, std::int64_t>> out;
    for (std::size_t k = 0; k + 1 < rE.size(); k += 2)
        out.emplace_back(rE[k], rE[k + 1]);
    std::sort(out.begin(), out.end());
    return out;
}

Vec3 ex_point(const Mesh& rM, std::int64_t P) {
    Vec3 x{0.0, 0.0, 0.0};
    for (std::size_t d = 0; d < rM.PointDim() && d < 3; ++d)
        x[d] = meshioplusplus::detail::read_double(rM.Points(),
                                                   static_cast<std::size_t>(P) * rM.PointDim() + d);
    return x;
}

// The signed volume a closed triangle and quad surface encloses.
double ex_enclosed(const Mesh& rM) {
    double v = 0.0;
    for (const auto cb : rM.CellRange()) {
        const std::size_t k = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::vector<Vec3> x;
            for (std::size_t i = 0; i < k; ++i)
                x.push_back(ex_point(rM, meshioplusplus::detail::read_int(cb.Conn(), c * k + i)));
            for (std::size_t i = 1; i + 1 < k; ++i) {
                const Vec3 &a = x[0], &b = x[i], &z = x[i + 1];
                v += (a[0] * (b[1] * z[2] - b[2] * z[1]) - a[1] * (b[0] * z[2] - b[2] * z[0]) +
                      a[2] * (b[0] * z[1] - b[1] * z[0])) /
                     6.0;
            }
        }
    }
    return v;
}

// Every directed edge of a surface once, and its reverse too: closed and
// consistently wound.
void ex_expect_closed_surface(const Mesh& rM) {
    std::set<std::pair<std::int64_t, std::int64_t>> directed;
    for (const auto cb : rM.CellRange()) {
        const std::size_t k = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            for (std::size_t i = 0; i < k; ++i) {
                const std::int64_t a = meshioplusplus::detail::read_int(cb.Conn(), c * k + i);
                const std::int64_t b =
                    meshioplusplus::detail::read_int(cb.Conn(), c * k + (i + 1) % k);
                EXPECT_TRUE(directed.insert({a, b}).second) << a << " " << b;
            }
    }
    for (const auto& [a, b] : directed)
        EXPECT_TRUE(directed.count({b, a})) << a << " " << b;
}

Mesh ex_polyline(const std::vector<std::vector<double>>& rPts, bool Closed) {
    std::vector<std::vector<std::int64_t>> lines;
    const auto n = static_cast<std::int64_t>(rPts.size());
    for (std::int64_t i = 0; i + 1 < n; ++i)
        lines.push_back({i, i + 1});
    if (Closed)
        lines.push_back({n - 1, 0});
    return mt::make_mesh(rPts, "line", lines);
}

}  // namespace

TEST(Extrude, ARectangleExtrudesToTheBox) {
    const Mesh rect = meshioplusplus::rectangle({{3, 2}}, {{1.0, 1.0}});
    ExtrudeOptions o;
    o.mDirection = {{0.0, 0.0, 1.0}};
    o.mNumLayers = 4;
    const SweepResult r = meshioplusplus::extrude(rect, o);
    const Mesh box = meshioplusplus::box({{3, 2, 4}});
    ASSERT_EQ(r.mMesh.NumPoints(), box.NumPoints());
    ASSERT_EQ(r.mMesh.PointDim(), 3u);
    for (std::size_t i = 0; i < 3 * box.NumPoints(); ++i)
        ASSERT_EQ(meshioplusplus::detail::read_double(r.mMesh.Points(), i),
                  meshioplusplus::detail::read_double(box.Points(), i))
            << i;
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "hexahedron");
    for (std::size_t i = 0; i < box.Cells(0).Conn().Size(); ++i)
        ASSERT_EQ(meshioplusplus::detail::read_int(r.mMesh.Cells(0).Conn(), i),
                  meshioplusplus::detail::read_int(box.Cells(0).Conn(), i));
    EXPECT_EQ(r.mNumHexahedra, 24);
    EXPECT_EQ(r.mNumFlipped, 0);
    EXPECT_DOUBLE_EQ(r.mMinScaledJacobian, 1.0);
    for (const char* side : {"xmin", "xmax", "ymin", "ymax"})
        EXPECT_EQ(ex_pairs(ex_entries(r.mMesh, side)), ex_pairs(ex_entries(box, side))) << side;
    EXPECT_EQ(ex_pairs(ex_entries(r.mMesh, "bottom")), ex_pairs(ex_entries(box, "zmin")));
    EXPECT_EQ(ex_pairs(ex_entries(r.mMesh, "top")), ex_pairs(ex_entries(box, "zmax")));
    std::vector<std::int64_t> sides;
    for (const char* side : {"xmin", "xmax", "ymin", "ymax"}) {
        const auto e = ex_entries(box, side);
        sides.insert(sides.end(), e.begin(), e.end());
    }
    EXPECT_EQ(ex_pairs(ex_entries(r.mMesh, "side")), ex_pairs(sides));
    EXPECT_EQ(r.mMesh.Region(r.mMesh.FindRegion("side", RegionKind::Side)).mDim, 2);
}

TEST(Extrude, TrianglesGiveWedgesAndGradingSpacesTheLevels) {
    const Mesh rect =
        meshioplusplus::rectangle({{2, 2}}, {{1.0, 1.0}}, {{0.0, 0.0, 0.0}}, "triangle");
    ExtrudeOptions o;
    o.mDirection = {{0.0, 0.0, 2.0}};
    o.mNumLayers = 3;
    o.mRatio = 2.0;
    o.mRecordLayers = true;
    const SweepResult r = meshioplusplus::extrude(rect, o);
    EXPECT_EQ(r.mNumWedges, 8 * 3);
    EXPECT_NEAR(ex_total(r.mMesh), 2.0, 1e-12);
    std::set<double> z;
    for (std::size_t p = 0; p < r.mMesh.NumPoints(); ++p)
        z.insert(ex_point(r.mMesh, static_cast<std::int64_t>(p))[2]);
    ASSERT_EQ(z.size(), 4u);
    const std::vector<double> want{0.0, 2.0 / 7.0, 6.0 / 7.0, 2.0};
    std::size_t i = 0;
    for (const double v : z)
        EXPECT_NEAR(v, want[i++], 1e-15);
    const NDArray& layer = r.mMesh.CellData(meshioplusplus::kExtrudeLayerName, 0);
    EXPECT_EQ(meshioplusplus::detail::read_int(layer, 0), 1);
    EXPECT_EQ(meshioplusplus::detail::read_int(layer, layer.Size() - 1), 3);
    // Explicit offsets give the same levels.
    ExtrudeOptions e;
    e.mDirection = o.mDirection;
    e.mOffsets = {1.0 / 7.0, 3.0 / 7.0, 1.0};
    const SweepResult re = meshioplusplus::extrude(rect, e);
    EXPECT_EQ(re.mNumWedges, 24);
}

TEST(Extrude, FacesTurnToTheDirection) {
    const Mesh rect = meshioplusplus::rectangle({{2, 3}});
    ExtrudeOptions o;
    o.mDirection = {{0.1, 0.0, -1.0}};
    o.mNumLayers = 2;
    const SweepResult r = meshioplusplus::extrude(rect, o);
    EXPECT_EQ(r.mNumFlipped, 6);
    EXPECT_EQ(meshioplusplus::compute_quality(r.mMesh).mNumInverted, 0);
    EXPECT_NEAR(ex_total(r.mMesh), 1.0, 1e-12);
    const auto conf = meshioplusplus::check_conformity(r.mMesh);
    EXPECT_TRUE(conf.mConforming);
    EXPECT_EQ(conf.mNumBoundaryFacets,
              static_cast<std::int64_t>(ex_entries(r.mMesh, "bottom").size() / 2 +
                                        ex_entries(r.mMesh, "top").size() / 2 +
                                        ex_entries(r.mMesh, "side").size() / 2));
}

TEST(Extrude, ACurveExtrudesToASurface) {
    std::vector<std::vector<double>> circle;
    for (int i = 0; i < 12; ++i)
        circle.push_back({std::cos(2 * kPi * i / 12), std::sin(2 * kPi * i / 12), 0.0});
    ExtrudeOptions o;
    o.mNumLayers = 3;
    const SweepResult r = meshioplusplus::extrude(ex_polyline(circle, true), o);
    EXPECT_EQ(r.mNumQuads, 36);
    // Counter-clockwise seen from +z, swept up: the quads face outward.
    const auto cb = r.mMesh.Cells(0);
    for (std::size_t c = 0; c < cb.NumCells(); ++c) {
        std::array<Vec3, 4> x;
        for (std::size_t k = 0; k < 4; ++k)
            x[k] = ex_point(r.mMesh, meshioplusplus::detail::read_int(cb.Conn(), c * 4 + k));
        const Vec3 e1{x[1][0] - x[0][0], x[1][1] - x[0][1], 0.0};
        const Vec3 e2{x[3][0] - x[0][0], x[3][1] - x[0][1], x[3][2] - x[0][2]};
        const Vec3 nrm{e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], 0.0};
        EXPECT_GT(nrm[0] * (x[0][0] + x[1][0]) + nrm[1] * (x[0][1] + x[1][1]), 0.0);
    }
    EXPECT_EQ(ex_entries(r.mMesh, "bottom").size(), 24u);
    EXPECT_EQ(ex_entries(r.mMesh, "top").size(), 24u);
    EXPECT_TRUE(ex_entries(r.mMesh, "side").empty());
    // An open curve has two ends, each sweeping to an edge per layer.
    circle.pop_back();
    const SweepResult ro = meshioplusplus::extrude(ex_polyline(circle, false), o);
    EXPECT_EQ(ex_entries(ro.mMesh, "side").size(), 2u * 2u * 3u);
    EXPECT_EQ(ro.mMesh.Region(ro.mMesh.FindRegion("side", RegionKind::Side)).mDim, 1);
}

TEST(Extrude, RegionsAndDataAreCarried) {
    Mesh rect = meshioplusplus::rectangle({{2, 1}});
    // A boundary line along y = 0 named as gmsh names it, a cell and a point.
    {
        std::vector<std::vector<std::int64_t>> rows{{0, 1}, {1, 2}};
        rect.AddCellBlock("line", mt::conn_from(rows));
        NDArray lines(meshioplusplus::DType::Int64, {2});
        lines.As<std::int64_t>()[0] = 2;
        lines.As<std::int64_t>()[1] = 3;
        rect.AddRegion(meshioplusplus::Region("wall", RegionKind::Cell, 1, 5, lines));
        NDArray cell(meshioplusplus::DType::Int64, {1});
        cell.As<std::int64_t>()[0] = 1;
        rect.AddRegion(meshioplusplus::Region("right", RegionKind::Cell, 2, -1, cell));
        NDArray pt(meshioplusplus::DType::Int64, {1});
        pt.As<std::int64_t>()[0] = 0;
        rect.AddRegion(meshioplusplus::Region("corner", RegionKind::Point, pt));
        NDArray id(meshioplusplus::DType::Int32, {2});
        id.As<std::int32_t>()[0] = 10;
        id.As<std::int32_t>()[1] = 20;
        NDArray lid(meshioplusplus::DType::Int32, {2});
        rect.AddCellData("id", {id, lid});
        NDArray t(meshioplusplus::DType::Float64, {rect.NumPoints()});
        for (std::size_t p = 0; p < rect.NumPoints(); ++p)
            t.As<double>()[p] = static_cast<double>(p);
        rect.AddPointData("t", t);
    }
    ExtrudeOptions o;
    o.mNumLayers = 2;
    const SweepResult r = meshioplusplus::extrude(rect, o);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    const NDArray& id = r.mMesh.CellData("id", 0);
    std::vector<std::int64_t> ids;
    for (std::size_t i = 0; i < 4; ++i)
        ids.push_back(meshioplusplus::detail::read_int(id, i));
    EXPECT_EQ(ids, (std::vector<std::int64_t>{10, 20, 10, 20}));
    EXPECT_EQ(ex_entries(r.mMesh, "right", RegionKind::Cell), (std::vector<std::int64_t>{1, 3}));
    EXPECT_EQ(ex_entries(r.mMesh, "corner", RegionKind::Point).size(), 3u);
    const std::size_t wall = r.mMesh.FindRegion("wall", RegionKind::Side);
    ASSERT_NE(wall, Mesh::npos);
    EXPECT_EQ(r.mMesh.Region(wall).mDim, 2);
    EXPECT_EQ(r.mMesh.Region(wall).mTag, 5);
    EXPECT_EQ(ex_pairs(ex_entries(r.mMesh, "wall")), ex_pairs(ex_entries(r.mMesh, "ymin")));
    const NDArray& t = r.mMesh.PointData("t");
    EXPECT_EQ(t.As<double>()[r.mMesh.NumPoints() - 1], static_cast<double>(rect.NumPoints() - 1));
    EXPECT_FALSE(r.mWarnings.empty());  // the line block's cell data is not carried
}

TEST(Extrude, Refusals) {
    const Mesh rect = meshioplusplus::rectangle({{2, 2}});
    ExtrudeOptions o;
    o.mDirection = {{0.0, 1.0, 0.0}};
    EXPECT_THROW(meshioplusplus::extrude(rect, o), std::invalid_argument);
    o.mDirection = {{0.0, 0.0, 0.0}};
    EXPECT_THROW(meshioplusplus::extrude(rect, o), std::invalid_argument);
    o.mDirection = {{0.0, 0.0, 1.0}};
    o.mNumLayers = 10;
    o.mMaxCells = 39;
    EXPECT_THROW(meshioplusplus::extrude(rect, o), std::invalid_argument);
    o.mMaxCells = 40;
    EXPECT_NO_THROW(meshioplusplus::extrude(rect, o));
    o.mNumLayers = 0;
    EXPECT_THROW(meshioplusplus::extrude(rect, o), std::invalid_argument);
    o.mNumLayers = 1;
    o.mOffsets = {0.5, 0.9};
    EXPECT_THROW(meshioplusplus::extrude(rect, o), std::invalid_argument);
    o.mOffsets = {0.5, 1.0};
    o.mRatio = 1.5;
    EXPECT_THROW(meshioplusplus::extrude(rect, o), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::extrude(meshioplusplus::box(), ExtrudeOptions{}),
                 meshioplusplus::Unsupported);
    EXPECT_THROW(meshioplusplus::extrude(mt::quad8_mesh(), ExtrudeOptions{}),
                 meshioplusplus::Unsupported);
}

TEST(Revolve, AFullTurnOfARectangleTouchingTheAxis) {
    const double R = 2.0, H = 1.5;
    const Mesh rect = meshioplusplus::rectangle({{4, 3}}, {{R, H}});
    RevolveOptions o;
    o.mNumLayers = 12;
    const SweepResult r = meshioplusplus::revolve(rect, o);
    EXPECT_TRUE(r.mFull);
    EXPECT_EQ(r.mNumOnAxis, 4);
    EXPECT_EQ(r.mNumWedges, 3 * 12);
    EXPECT_EQ(r.mNumHexahedra, 9 * 12);
    EXPECT_EQ(r.mNumPyramids + r.mNumTetrahedra, 0);
    EXPECT_EQ(meshioplusplus::compute_quality(r.mMesh).mNumInverted, 0);
    EXPECT_GT(r.mMinScaledJacobian, 0.0);
    const double L = 12.0;
    EXPECT_NEAR(ex_total(r.mMesh), L / 2.0 * std::sin(2.0 * kPi / L) * R * R * H, 1e-12);
    // Points: the profile's, then 11 new levels of the 16 off-axis points.
    EXPECT_EQ(r.mMesh.NumPoints(), 20u + 11u * 16u);
    const auto conf = meshioplusplus::check_conformity(r.mMesh);
    EXPECT_TRUE(conf.mConforming);
    EXPECT_EQ(conf.mNumBoundaryFacets, 132);
    EXPECT_EQ(ex_entries(r.mMesh, "side").size() / 2, 132u);
    EXPECT_EQ(r.mMesh.FindRegion("start", RegionKind::Side), Mesh::npos);
    EXPECT_EQ(r.mMesh.FindRegion("end", RegionKind::Side), Mesh::npos);
    // The axis edge of the rectangle sweeps to nothing; its xmin region is empty.
    EXPECT_TRUE(ex_entries(r.mMesh, "xmin").empty());
    EXPECT_EQ(ex_entries(r.mMesh, "xmax").size() / 2, 36u);
}

TEST(Revolve, APartialTurnHasAStartAndAnEnd) {
    const double R = 1.0, H = 2.0;
    const Mesh rect = meshioplusplus::rectangle({{2, 2}}, {{R, H}}, {{0.0, 0.0, 0.0}});
    RevolveOptions o;
    o.mAngle = -90.0;
    o.mNumLayers = 3;
    const SweepResult r = meshioplusplus::revolve(rect, o);
    EXPECT_FALSE(r.mFull);
    const double th = kPi / 2.0;
    EXPECT_NEAR(ex_total(r.mMesh), 3.0 / 2.0 * std::sin(th / 3.0) * R * R * H, 1e-12);
    EXPECT_EQ(ex_entries(r.mMesh, "start").size() / 2, 4u);
    EXPECT_EQ(ex_entries(r.mMesh, "end").size() / 2, 4u);
    const auto conf = meshioplusplus::check_conformity(r.mMesh);
    EXPECT_TRUE(conf.mConforming);
    EXPECT_EQ(conf.mNumBoundaryFacets,
              static_cast<std::int64_t>((ex_entries(r.mMesh, "start").size() +
                                         ex_entries(r.mMesh, "end").size() +
                                         ex_entries(r.mMesh, "side").size()) /
                                        2));
}

TEST(Revolve, ASemicircleTurnsIntoASphere) {
    // From the top pole down to the bottom one, so the surface faces out.
    std::vector<std::vector<double>> arc;
    const int N = 8;
    for (int i = 0; i <= N; ++i) {
        const double phi = kPi * (N - i) / N;
        arc.push_back({i == 0 || i == N ? 0.0 : std::sin(phi), -std::cos(phi), 0.0});
    }
    RevolveOptions o;
    o.mNumLayers = 16;
    const SweepResult r = meshioplusplus::revolve(ex_polyline(arc, false), o);
    EXPECT_EQ(r.mNumTriangles, 2 * 16);
    EXPECT_EQ(r.mNumQuads, (N - 2) * 16);
    EXPECT_EQ(r.mNumOnAxis, 2);
    ex_expect_closed_surface(r.mMesh);
    const double v = ex_enclosed(r.mMesh);
    EXPECT_GT(v, 0.0);
    // An inscribed polyhedron: a little less than the ball's volume.
    EXPECT_LT(v, 4.0 / 3.0 * kPi);
    EXPECT_GT(v, 0.9 * 4.0 / 3.0 * kPi);
    EXPECT_TRUE(ex_entries(r.mMesh, "side").empty());
}

TEST(Revolve, CellsAroundTheAxisCollapse) {
    // One triangle with a corner on the axis (y): pyramids.
    const Mesh tri = mt::make_mesh({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}},
                                   "triangle", {{0, 1, 2}});
    RevolveOptions o;
    o.mNumLayers = 6;
    const SweepResult rt = meshioplusplus::revolve(tri, o);
    EXPECT_EQ(rt.mNumPyramids, 6);
    EXPECT_EQ(meshioplusplus::compute_quality(rt.mMesh).mNumInverted, 0);
    // One quad with a corner on the axis: split into a pyramid and a wedge.
    const Mesh quad = mt::make_mesh(
        {{1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {0.0, 1.0, 0.0}}, "quad",
        {{0, 1, 2, 3}});
    const SweepResult rq = meshioplusplus::revolve(quad, o);
    EXPECT_EQ(rq.mNumPyramids + rq.mNumWedges, 12);
    EXPECT_EQ(rq.mNumPyramids, 6);
    for (const double v : ex_volumes(rq.mMesh))
        EXPECT_GT(v, 0.0);
    EXPECT_TRUE(meshioplusplus::check_conformity(rq.mMesh).mConforming);
}

TEST(Revolve, Refusals) {
    RevolveOptions o;
    // One quad across the axis, and two on either side of it.
    const Mesh across =
        meshioplusplus::rectangle({{1, 1}}, {{2.0, 1.0}}, {{-1.0, 0.0, 0.0}});
    EXPECT_THROW(meshioplusplus::revolve(across, o), std::invalid_argument);
    const Mesh both = meshioplusplus::rectangle({{2, 1}}, {{2.0, 1.0}}, {{-1.0, 0.0, 0.0}});
    EXPECT_THROW(meshioplusplus::revolve(both, o), std::invalid_argument);
    RevolveOptions half = o;
    half.mAngle = 90.0;
    half.mNumLayers = 2;
    const SweepResult rb = meshioplusplus::revolve(both, half);
    EXPECT_EQ(rb.mWarnings.size(), 1u);
    EXPECT_EQ(meshioplusplus::compute_quality(rb.mMesh).mNumInverted, 0);
    const Mesh rect = meshioplusplus::rectangle({{2, 1}}, {{1.0, 1.0}}, {{1.0, 0.0, 0.0}});
    o.mNumLayers = 2;
    EXPECT_THROW(meshioplusplus::revolve(rect, o), std::invalid_argument);
    o.mAngle = 200.0;
    o.mNumLayers = 1;
    EXPECT_THROW(meshioplusplus::revolve(rect, o), std::invalid_argument);
    o.mAngle = 0.0;
    EXPECT_THROW(meshioplusplus::revolve(rect, o), std::invalid_argument);
    o.mAngle = 90.0;
    o.mAxis = {{0.0, 0.0, 1.0}};  // the rectangle turns within its own plane
    EXPECT_THROW(meshioplusplus::revolve(rect, o), std::invalid_argument);
}
