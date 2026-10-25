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
// mesh_cad_surface: conforming surface meshes of B-reps. The solids of
// brep_solids.hpp (a box, a capped cylinder through a seam, a sphere through
// its poles) and the trimmed plate of CAD JSON are meshed and checked for:
// watertightness (every triangle edge used exactly twice, in opposite
// directions), area and enclosed volume against the analytic values with
// second-order convergence, every point on the exact surface, edge lengths
// and angles on planar faces, the Hausdorff distance to a dense reference,
// determinism across runs and thread counts, a face that fails being
// reported and never meshed, and, where OpenCASCADE is linked, the area of
// BRepMesh's triangulation of the same STEP solid.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

// Project includes
#include "analytic_surfaces.hpp"
#include "brep_solids.hpp"
#include "mesh_fixtures.hpp"
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::MeshCadFaceStatus;
using mio::MeshCadSurfaceOptions;
using mio::MeshCadSurfaceResult;

constexpr double kPi = 3.14159265358979323846;

MeshCadSurfaceOptions size(double H) {
    MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    return o;
}

struct Surface {
    double mArea = 0.0;
    double mVolume = 0.0;  // by the divergence theorem
    std::size_t mTriangles = 0;
};

// The triangle block's area and enclosed volume.
Surface measure(const mio::Mesh& rM) {
    Surface s;
    const double* x = rM.Points().As<double>();
    const auto cells = rM.Cells(0);
    EXPECT_EQ(cells.Type(), "triangle");
    const mio::NDArray& c = cells.Conn();
    s.mTriangles = cells.NumCells();
    for (std::size_t t = 0; t < cells.NumCells(); ++t) {
        const double* a = x + 3 * mio::detail::read_int(c, 3 * t);
        const double* b = x + 3 * mio::detail::read_int(c, 3 * t + 1);
        const double* d = x + 3 * mio::detail::read_int(c, 3 * t + 2);
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double v[3] = {d[0] - a[0], d[1] - a[1], d[2] - a[2]};
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                             u[0] * v[1] - u[1] * v[0]};
        s.mArea += 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        s.mVolume += (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
                      a[2] * (b[0] * d[1] - b[1] * d[0])) /
                     6.0;
    }
    return s;
}

// Every triangle edge used exactly twice, once each way; no degenerate
// triangle (repeated vertex).
void expect_watertight(const MeshCadSurfaceResult& rR) {
    EXPECT_EQ(rR.mNumFailed, 0);
    EXPECT_EQ(rR.mNumOpenEdges, 0);
    EXPECT_EQ(rR.mNumNonManifoldEdges, 0);
    EXPECT_EQ(rR.mNumMisorientedEdges, 0);
    std::map<std::pair<std::int64_t, std::int64_t>, int> directed;
    const auto cells = rR.mMesh.Cells(0);
    const mio::NDArray& c = cells.Conn();
    for (std::size_t t = 0; t < cells.NumCells(); ++t)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t a = mio::detail::read_int(c, 3 * t + k);
            const std::int64_t b = mio::detail::read_int(c, 3 * t + (k + 1) % 3);
            ASSERT_NE(a, b);
            ++directed[{a, b}];
        }
    for (const auto& [e, n] : directed) {
        EXPECT_EQ(n, 1) << e.first << "-" << e.second << " used twice the same way";
        const auto it = directed.find({e.second, e.first});
        EXPECT_TRUE(it != directed.end() && it->second == 1)
            << e.first << "-" << e.second << " has no opposite use";
    }
}

std::uint64_t digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

TEST(MeshCad, PlanarSquare) {
    const double h = 0.1;
    const auto r = mio::mesh_cad_surface(brep_solids::square(1.0), size(h));
    ASSERT_EQ(r.mNumFailed, 0) << r.mFaces[0].mReason;
    const Surface s = measure(r.mMesh);
    EXPECT_NEAR(s.mArea, 1.0, 1e-12);
    EXPECT_EQ(r.mNumOpenEdges, 40);  // the square's boundary, ten segments a side
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    EXPECT_EQ(r.mNumMisorientedEdges, 0);
    EXPECT_EQ(r.mNumEdgeSegments, 40);
    EXPECT_GE(r.mMinEdge, 0.3 * h);
    EXPECT_LE(r.mMaxEdge, 1.5 * h);
    EXPECT_GT(r.mMinAngleDeg, 10.0);
    // Counter-clockwise seen from +z: the face's normal.
    const double* x = r.mMesh.Points().As<double>();
    const auto cells = r.mMesh.Cells(0);
    const mio::NDArray& c = cells.Conn();
    for (std::size_t t = 0; t < cells.NumCells(); ++t) {
        const double* a = x + 3 * mio::detail::read_int(c, 3 * t);
        const double* b = x + 3 * mio::detail::read_int(c, 3 * t + 1);
        const double* d = x + 3 * mio::detail::read_int(c, 3 * t + 2);
        EXPECT_GT((b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0]), 0.0);
    }
}

TEST(MeshCad, OutputLayout) {
    const auto r = mio::mesh_cad_surface(brep_solids::box(), size(0.25));
    const mio::Mesh& m = r.mMesh;
    ASSERT_EQ(m.NumCellBlocks(), 2u);
    EXPECT_EQ(m.Cells(0).Type(), "triangle");
    EXPECT_EQ(m.Cells(1).Type(), "line");
    EXPECT_EQ(static_cast<std::int64_t>(m.Cells(1).NumCells()), r.mNumEdgeSegments);
    EXPECT_EQ(r.mNumEdgeSegments, 12 * 4);
    ASSERT_TRUE(m.HasCellData(mio::kMeshCadFaceName));
    ASSERT_TRUE(m.HasCellData(mio::kMeshCadEdgeName));
    const mio::NDArray& tri_face = m.CellData(mio::kMeshCadFaceName, 0);
    const mio::NDArray& line_face = m.CellData(mio::kMeshCadFaceName, 1);
    const mio::NDArray& tri_edge = m.CellData(mio::kMeshCadEdgeName, 0);
    const mio::NDArray& line_edge = m.CellData(mio::kMeshCadEdgeName, 1);
    for (std::size_t i = 0; i < tri_face.Size(); ++i) {
        const std::int64_t f = mio::detail::read_int(tri_face, i);
        EXPECT_GE(f, 0);
        EXPECT_LT(f, 6);
        EXPECT_EQ(mio::detail::read_int(tri_edge, i), -1);
    }
    for (std::size_t i = 0; i < line_face.Size(); ++i) {
        EXPECT_EQ(mio::detail::read_int(line_face, i), -1);
        EXPECT_GE(mio::detail::read_int(line_edge, i), 0);
    }
    std::size_t in_regions = 0;
    for (int f = 0; f < 6; ++f) {
        const std::size_t ri = m.FindRegion("face:" + std::to_string(f), mio::RegionKind::Cell);
        ASSERT_NE(ri, mio::Mesh::npos);
        const mio::Region& region = m.Region(ri);
        in_regions += region.NumEntries();
        for (std::size_t k = 0; k < region.NumEntries(); ++k)
            EXPECT_EQ(
                mio::detail::read_int(tri_face, static_cast<std::size_t>(region.Entries()[k])), f);
    }
    EXPECT_EQ(in_regions, m.Cells(0).NumCells());
    MeshCadSurfaceOptions bare = size(0.25);
    bare.mEdges = false;
    bare.mRecordFaceIds = false;
    const auto b = mio::mesh_cad_surface(brep_solids::box(), bare);
    EXPECT_EQ(b.mMesh.NumCellBlocks(), 1u);
    EXPECT_FALSE(b.mMesh.HasCellData(mio::kMeshCadFaceName));
    EXPECT_FALSE(b.mMesh.HasRegion("face:0"));
    EXPECT_EQ(digest(mio::mesh_cad_surface(brep_solids::box(), bare).mMesh), digest(b.mMesh));
}

TEST(MeshCad, UnitBoxIsWatertight) {
    for (const double h : {0.3, 0.1, 0.05}) {
        const auto r = mio::mesh_cad_surface(brep_solids::box(), size(h));
        expect_watertight(r);
        const Surface s = measure(r.mMesh);
        EXPECT_NEAR(s.mArea, 6.0, 1e-12);
        EXPECT_NEAR(s.mVolume, 1.0, 1e-12);
        EXPECT_GE(r.mMinEdge, 0.3 * h) << h;
        EXPECT_LE(r.mMaxEdge, 1.5 * h) << h;
        EXPECT_GT(r.mMinAngleDeg, 10.0) << h;
    }
    // A box of other proportions: still exact.
    const auto r = mio::mesh_cad_surface(brep_solids::box(2.0, 0.5, 1.3), size(0.1));
    expect_watertight(r);
    const Surface s = measure(r.mMesh);
    EXPECT_NEAR(s.mArea, 2.0 * (2.0 * 0.5 + 2.0 * 1.3 + 0.5 * 1.3), 1e-12);
    EXPECT_NEAR(s.mVolume, 2.0 * 0.5 * 1.3, 1e-12);
}

TEST(MeshCad, CappedCylinderConvergesAtSecondOrder) {
    const double R = 1.0, H = 2.0;
    const double area = 2.0 * kPi * R * H + 2.0 * kPi * R * R;
    const double volume = kPi * R * R * H;
    double prev_area = 0.0, prev_volume = 0.0;
    for (const double h : {0.2, 0.1, 0.05}) {
        const auto r = mio::mesh_cad_surface(brep_solids::capped_cylinder(R, H), size(h));
        expect_watertight(r);
        const Surface s = measure(r.mMesh);
        const double ea = std::fabs(s.mArea - area);
        const double ev = std::fabs(s.mVolume - volume);
        EXPECT_LT(ea, 0.5 * h * h * area) << h;
        EXPECT_LT(ev, 0.5 * h * h * volume) << h;
        if (prev_area > 0.0) {
            EXPECT_LT(ea, prev_area / 2.5) << h;  // about 4 for second order
            EXPECT_LT(ev, prev_volume / 2.5) << h;
        }
        prev_area = ea;
        prev_volume = ev;
        // Every point on the exact surface: the wall or a cap.
        const double* x = r.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
            const double rad = std::hypot(x[3 * i], x[3 * i + 1]);
            const bool wall = std::fabs(rad - R) < 1e-12;
            const bool cap =
                (std::fabs(x[3 * i + 2]) < 1e-12 || std::fabs(x[3 * i + 2] - H) < 1e-12) &&
                rad <= R + 1e-12;
            EXPECT_TRUE(wall || cap) << i;
        }
    }
}

TEST(MeshCad, CappedCylinderHausdorffToADenseReference) {
    const double R = 1.0, H = 2.0, h = 0.1;
    MeshCadSurfaceOptions o = size(h);
    o.mEdges = false;
    const auto r = mio::mesh_cad_surface(brep_solids::capped_cylinder(R, H), o);
    // The reference: a fine faceting whose own sagitta is about 1e-5.
    const mio::Mesh ref = analytic::capped_cylinder(R, H, 64, 720);
    mio::HausdorffOptions ho;
    ho.mFaceSamples = 3;
    const auto d = mio::hausdorff_distance(r.mMesh, ref, ho);
    // A chord of length l on a circle of radius R sags l^2 / (8 R): the
    // longest edge bounds the distance, plus the reference's own sagitta.
    EXPECT_LE(r.mMaxEdge, 1.5 * h);
    EXPECT_LT(d.mDistance, r.mMaxEdge * r.mMaxEdge / (8.0 * R) + 2e-5);
    EXPECT_GT(d.mDistance, 0.0);
}

TEST(MeshCad, SphereThroughThePoles) {
    const double R = 1.0;
    double prev = 0.0;
    for (const double h : {0.2, 0.1}) {
        const auto r = mio::mesh_cad_surface(brep_solids::sphere(R), size(h));
        expect_watertight(r);
        const Surface s = measure(r.mMesh);
        const double ea = std::fabs(s.mArea - 4.0 * kPi * R * R);
        EXPECT_LT(ea, 0.5 * h * h * 4.0 * kPi) << h;
        EXPECT_NEAR(s.mVolume, 4.0 / 3.0 * kPi, 0.5 * h * h * 4.0) << h;
        if (prev > 0.0)
            EXPECT_LT(ea, prev / 2.5);
        prev = ea;
        const double* x = r.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
            EXPECT_NEAR(std::sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] +
                                  x[3 * i + 2] * x[3 * i + 2]),
                        R, 1e-12);
    }
}

// tests/cpp/test_iga.cpp's plate: a 4 x 4 plane with a diamond hole, a
// single degree-1 trim whose corners are knots (kept as samples).
mio::BrepInfo trimmed_plate() {
    const char* plate = R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2,
        "swapped_surface_normal": false, "surface": {"is_trimmed": true, "is_rational": false,
        "degrees": [1, 1], "knot_vectors": [[0, 4], [0, 4]],
        "control_points": [[1, [0, 0, 0, 1]], [2, [4, 0, 0, 1]], [3, [0, 4, 0, 1]], [4, [4, 4, 0, 1]]]},
        "boundary_loops": [
          {"loop_type": "outer", "trimming_curves": [
            {"trim_index": 0, "curve_direction": true, "parameter_curve": {"degree": 1,
              "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 0, 0, 1]], [[4, 0, 0, 1]]]}},
            {"trim_index": 1, "curve_direction": true, "parameter_curve": {"degree": 1,
              "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[4, 0, 0, 1]], [[4, 4, 0, 1]]]}},
            {"trim_index": 2, "curve_direction": true, "parameter_curve": {"degree": 1,
              "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[4, 4, 0, 1]], [[0, 4, 0, 1]]]}},
            {"trim_index": 3, "curve_direction": true, "parameter_curve": {"degree": 1,
              "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 4, 0, 1]], [[0, 0, 0, 1]]]}}]},
          {"loop_type": "inner", "trimming_curves": [
            {"trim_index": 4, "curve_direction": true, "parameter_curve": {"degree": 1,
              "knot_vector": [0, 0, 0.25, 0.5, 0.75, 1, 1], "active_range": [0, 1],
              "control_points": [[[2, 0.8, 0, 1]], [[3.2, 2, 0, 1]], [[2, 3.2, 0, 1]], [[0.8, 2, 0, 1]], [[2, 0.8, 0, 1]]]}}]}],
        "embedded_edges": []}], "edges": [], "vertices": []}]})";
    const std::string path = mt::temp_path(".cad.json");
    std::ofstream(path) << plate;
    return mio::read_cad_json_info(path);
}

TEST(MeshCad, TrimmedPlateWithAHole) {
    const mio::BrepInfo info = trimmed_plate();
    for (const double h : {0.5, 0.2, 0.1}) {
        const auto r = mio::mesh_cad_surface(info, size(h));
        ASSERT_EQ(r.mNumFailed, 0) << r.mFaces[0].mReason;
        EXPECT_NEAR(measure(r.mMesh).mArea, 16.0 - 2.88, 1e-12) << h;
        EXPECT_GT(r.mMinAngleDeg, 10.0) << h;
        // No refinement yet: beside the hole's corners, which sit at 45
        // degrees to the lattice, an edge can reach about 1.56 h.
        EXPECT_LE(r.mMaxEdge, 1.6 * h) << h;
        EXPECT_GE(r.mMinEdge, 0.3 * h) << h;
    }
}

// A quarter of a cylinder of radius R and length L as one rational patch with
// no loops (the face is its whole parameter rectangle), its second parameter
// running over [0, V].
mio::BrepInfo quarter_cylinder(double R, double L, double V) {
    const double w = std::sqrt(0.5);
    mio::BrepInfo info;
    mio::NurbsPatch p;
    p.mDegrees = {2, 1};
    p.mKnots = {{0, 0, 0, 1, 1, 1}, {0, 0, V, V}};
    p.mControlPoints = mio::NDArray(mio::DType::Float64, {6, 3});
    p.mWeights = mio::NDArray(mio::DType::Float64, {6});
    const double ring[3][3] = {{R, 0, 1}, {R, R, w}, {0, R, 1}};
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 3; ++i) {
            double* c = p.mControlPoints.As<double>() + 3 * (j * 3 + i);
            c[0] = ring[i][0];
            c[1] = ring[i][1];
            c[2] = j * L;
            p.mWeights.As<double>()[j * 3 + i] = ring[i][2];
        }
    info.mPatches.push_back(p);
    mio::BrepFace f;
    f.mId = 1;
    f.mPatch = 0;
    f.mRational = true;
    info.mFaces.push_back(f);
    return info;
}

TEST(MeshCad, UntrimmedPatchUsesItsParameterRectangle) {
    // The second case has the Scordelis-Lo roof's proportions: |dS/du| about
    // 40 times |dS/dv|, which a Delaunay triangulation in raw parameters turns
    // into triangles across the whole roof.
    struct Case {
        double mR, mL, mV, mH;
    };
    for (const Case c : {Case{2.0, 3.0, 1.0, 0.2}, Case{25.0, 50.0, 50.0, 4.0}}) {
        const mio::BrepInfo info = quarter_cylinder(c.mR, c.mL, c.mV);
        const double exact = 0.5 * kPi * c.mR * c.mL;
        double prev = 0.0;
        for (const double h : {c.mH, 0.5 * c.mH}) {
            const auto r = mio::mesh_cad_surface(info, size(h));
            ASSERT_EQ(r.mNumFailed, 0) << r.mFaces[0].mReason;
            EXPECT_EQ(r.mNumMisorientedEdges, 0);
            EXPECT_GT(r.mMinAngleDeg, 10.0) << c.mR << " " << h;
            EXPECT_LE(r.mMaxEdge, 1.5 * h) << c.mR << " " << h;
            const double e = std::fabs(measure(r.mMesh).mArea - exact);
            EXPECT_LT(e, 0.5 * (h / c.mR) * (h / c.mR) * exact) << c.mR << " " << h;
            if (prev > 0.0)
                EXPECT_LT(e, prev / 2.5) << c.mR << " " << h;
            prev = e;
        }
    }
}

TEST(MeshCad, WithoutEdgeTopologyTheTrimsArePaired) {
    mio::BrepInfo info = brep_solids::box();
    info.mEdges.clear();
    const auto paired = mio::mesh_cad_surface(info, size(0.2));
    expect_watertight(paired);
    EXPECT_EQ(paired.mNumEdges, 12);
    // Not paired: every face has its own copies of its edges' points.
    MeshCadSurfaceOptions o = size(0.2);
    o.mMatchUnreferencedTrims = false;
    const auto apart = mio::mesh_cad_surface(info, o);
    EXPECT_EQ(apart.mNumFailed, 0);
    EXPECT_EQ(apart.mNumEdges, 24);
    EXPECT_GT(apart.mNumOpenEdges, 0);
}

TEST(MeshCad, UntrimmedPatchesWithoutTopologyAreSewnByTheirImages) {
    const auto r = mio::mesh_cad_surface(brep_solids::patch_box(), size(0.1));
    expect_watertight(r);
    EXPECT_EQ(r.mNumEdges, 12);
    const Surface s = measure(r.mMesh);
    EXPECT_NEAR(s.mArea, 6.0, 1e-12);
    EXPECT_NEAR(s.mVolume, 1.0, 1e-12);
}

TEST(MeshCad, AFailedFaceIsReportedNeverMeshed) {
    mio::BrepInfo info = brep_solids::box();
    info.mFaces[2].mPatch = 99;
    const auto r = mio::mesh_cad_surface(info, size(0.25));
    EXPECT_EQ(r.mNumFailed, 1);
    EXPECT_EQ(r.mFaces[2].mStatus, MeshCadFaceStatus::Unsupported);
    EXPECT_FALSE(r.mFaces[2].mReason.empty());
    EXPECT_EQ(r.mFaces[2].mNumTriangles, 0);
    ASSERT_FALSE(r.mWarnings.empty());
    // The other five are meshed, and the hole is visible.
    EXPECT_EQ(r.mNumOpenEdges, 4 * 4);
    const std::size_t ri = r.mMesh.FindRegion("face:2", mio::RegionKind::Cell);
    ASSERT_NE(ri, mio::Mesh::npos);
    EXPECT_EQ(r.mMesh.Region(ri).NumEntries(), 0u);
    EXPECT_NEAR(measure(r.mMesh).mArea, 5.0, 1e-12);
    MeshCadSurfaceOptions strict = size(0.25);
    strict.mStrict = true;
    EXPECT_THROW(mio::mesh_cad_surface(info, strict), std::runtime_error);
}

TEST(MeshCad, AClosedEdgeStartedElsewhereIsMeshed) {
    // Roadmap 6.1.2: the top cap's rim starts a quarter turn, a half turn and
    // three sixteenths of a turn away from the band's. The cap's circle is a loop of
    // its own, so it is rotated onto the band's corner and the solid closes,
    // converging at second order as the consistent model does.
    const double R = 1.0, H = 2.0;
    const double area = 2.0 * kPi * R * H + 2.0 * kPi * R * R;
    const double volume = kPi * R * R * H;
    // (Starts whose 16 builder samples include the circle's axis extremes:
    // brep_nurbs_face sizes a plane's patch from those samples.)
    for (const double start : {0.5 * kPi, kPi, 0.375 * kPi}) {
        double prev_area = 0.0;
        for (const double h : {0.2, 0.1}) {
            const auto r =
                mio::mesh_cad_surface(brep_solids::capped_cylinder(R, H, start), size(h));
            expect_watertight(r);
            EXPECT_LT(r.mMaxEdgeGap, 1e-12) << start;
            const Surface s = measure(r.mMesh);
            const double ea = std::fabs(s.mArea - area);
            EXPECT_LT(ea, 0.5 * h * h * area) << start << " " << h;
            EXPECT_LT(std::fabs(s.mVolume - volume), 0.5 * h * h * volume) << start << " " << h;
            if (prev_area > 0.0)
                EXPECT_LT(ea, prev_area / 2.5) << start;
            prev_area = ea;
        }
    }
}

TEST(MeshCad, AClosedEdgeIsAnchoredOnTheUseThatCannotMove) {
    // The discs are listed first, so each rim's first use is a disc's; the
    // band's use runs between its seams and cannot start elsewhere, so it
    // becomes the master and the discs are rotated onto it.
    const double R = 1.0, H = 2.0;
    const auto r =
        mio::mesh_cad_surface(brep_solids::capped_cylinder(R, H, 0.5 * kPi, true), size(0.2));
    expect_watertight(r);
    const Surface s = measure(r.mMesh);
    EXPECT_NEAR(s.mVolume, kPi * R * R * H, 0.5 * 0.04 * kPi * R * R * H);
    // A model whose uses start together is untouched by the reordering apart
    // from the face order.
    const auto a = mio::mesh_cad_surface(brep_solids::capped_cylinder(R, H, 0.0, true), size(0.2));
    const auto b = mio::mesh_cad_surface(brep_solids::capped_cylinder(R, H), size(0.2));
    expect_watertight(a);
    EXPECT_EQ(a.mNumPoints, b.mNumPoints);
    EXPECT_EQ(a.mNumTriangles, b.mNumTriangles);
}

TEST(MeshCad, AClosedEdgeNoRotationResolvesIsReported) {
    // Both uses of the shared rim run between their bands' seams, a quarter
    // turn apart: the faces that cannot be carried are reported, never
    // meshed wrong, and the reason says why.
    const auto r = mio::mesh_cad_surface(brep_solids::stacked_bands(1.0, 1.0), size(0.2));
    EXPECT_GE(r.mNumFailed, 1);
    bool named = false;
    for (const auto& f : r.mFaces)
        if (f.mStatus != MeshCadFaceStatus::Ok) {
            EXPECT_EQ(f.mStatus, MeshCadFaceStatus::ProjectionFailed) << f.mReason;
            EXPECT_EQ(f.mNumTriangles, 0);
            named = named || f.mReason.find("start at different points") != std::string::npos;
        }
    EXPECT_TRUE(named);
    EXPECT_GT(r.mNumOpenEdges, 0);
}

TEST(MeshCad, ASeamListedFromEitherEndIsTheSameSeam) {
    // The band's seam is used twice, once up and once down (that is what a
    // seam is); listing one of its trims from the other end of its curve,
    // with the direction flag flipped to match, describes the same loop.
    const double R = 1.0, H = 2.0;
    mio::BrepInfo info = brep_solids::capped_cylinder(R, H);
    bool flipped = false;
    for (mio::BrepTrim& t : info.mFaces[0].mLoops[0].mTrims) {
        const double* p = t.mCurve.mControlPoints.As<double>();
        const std::size_t n = t.mCurve.mControlPoints.Shape()[0];
        if (!flipped && n >= 2 && p[0] == p[3 * (n - 1)]) {  // constant u: a seam
            t.mCurve = mio::detail::brep_curve_reverse(t.mCurve);
            t.mDirection = !t.mDirection;
            flipped = true;
        }
    }
    ASSERT_TRUE(flipped);
    const auto r = mio::mesh_cad_surface(info, size(0.2));
    const auto ref = mio::mesh_cad_surface(brep_solids::capped_cylinder(R, H), size(0.2));
    expect_watertight(r);
    EXPECT_EQ(r.mNumPoints, ref.mNumPoints);
    EXPECT_EQ(r.mNumTriangles, ref.mNumTriangles);
    EXPECT_NEAR(measure(r.mMesh).mVolume, measure(ref.mMesh).mVolume, 1e-12);
}

// Twice the signed area of the line block's segments of one B-rep edge
// (the shoelace sum over the edge as its master runs).
double shoelace(const mio::Mesh& rM, std::int64_t Edge) {
    const double* x = rM.Points().As<double>();
    const auto lines = rM.Cells(1);
    const mio::NDArray& c = lines.Conn();
    const mio::NDArray& ids = rM.CellData("cad:edge", 1);
    double sum = 0.0;
    for (std::size_t s = 0; s < lines.NumCells(); ++s) {
        if (mio::detail::read_int(ids, s) != Edge)
            continue;
        const double* a = x + 3 * mio::detail::read_int(c, 2 * s);
        const double* b = x + 3 * mio::detail::read_int(c, 2 * s + 1);
        sum += a[0] * b[1] - b[0] * a[1];
    }
    return sum;
}

TEST(MeshCad, ABoundaryCrossingAtACoarseSizeIsSplit) {
    // Roadmap 6.1.3's probe: a plate with a hole 0.02 from its rim, meshed at
    // h = 1.1. The rim's coarse chord passes inside the hole's sample, so the
    // sampled loops cross; the pieces involved are halved until they do not,
    // and the plate is meshed, its area exactly that of its boundary polygon
    // (the rim counter-clockwise, the hole clockwise, as the loops run).
    const double R = 2.0, Rh = 1.0, Cx = 0.98;
    const auto r = mio::mesh_cad_surface(brep_solids::disc_with_hole(R, Rh, Cx), size(1.1));
    ASSERT_EQ(r.mNumFailed, 0) << r.mFaces[0].mReason;
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    EXPECT_EQ(r.mNumMisorientedEdges, 0);
    const Surface s = measure(r.mMesh);
    const double polygon = 0.5 * (shoelace(r.mMesh, 0) + shoelace(r.mMesh, 1));
    EXPECT_GT(shoelace(r.mMesh, 0), 0.0);
    EXPECT_LT(shoelace(r.mMesh, 1), 0.0);
    EXPECT_NEAR(s.mArea, polygon, 1e-12 * polygon);
    EXPECT_LT(s.mArea, kPi * (R * R - Rh * Rh));
    EXPECT_GT(s.mArea, 0.9 * kPi * (R * R - Rh * Rh));
    // Both curves gained points where they come close, and every boundary
    // point is on its circle.
    const double* x = r.mMesh.Points().As<double>();
    std::size_t near_rim = 0;  // rim points beside the hole: 345 and 15 degrees, and 0
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        EXPECT_EQ(x[3 * i + 2], 0.0);
        const double ro = std::hypot(x[3 * i], x[3 * i + 1]);
        const double rh = std::hypot(x[3 * i] - Cx, x[3 * i + 1]);
        EXPECT_LE(ro, R + 1e-12);
        EXPECT_GE(rh, Rh - 1e-12);
        if (std::fabs(ro - R) < 1e-12 && x[3 * i] > 1.9)
            ++near_rim;
    }
    EXPECT_GE(near_rim, 3u);
    // The same mesh on every run and thread count.
    const std::uint64_t ref = digest(r.mMesh);
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
        EXPECT_EQ(
            digest(mio::mesh_cad_surface(brep_solids::disc_with_hole(R, Rh, Cx), size(1.1)).mMesh),
            ref)
            << threads;
    }
    omp_set_num_threads(saved);
#endif
    EXPECT_EQ(
        digest(mio::mesh_cad_surface(brep_solids::disc_with_hole(R, Rh, Cx), size(1.1)).mMesh),
        ref);
}

TEST(MeshCad, ACrossingIsSplitInEveryFaceThatUsesTheEdge) {
    // The same plate with a thickness: the caps' loops cross, and the rims
    // their splits move are shared with the bands, which take the same new
    // points, so the solid stays closed.
    const double R = 2.0, Rh = 1.0, Cx = 0.98, T = 0.5;
    const auto r =
        mio::mesh_cad_surface(brep_solids::thick_disc_with_hole(R, Rh, Cx, T), size(1.1));
    for (const auto& f : r.mFaces)
        EXPECT_EQ(f.mStatus, MeshCadFaceStatus::Ok) << f.mFace << ": " << f.mReason;
    expect_watertight(r);
    const double volume = kPi * (R * R - Rh * Rh) * T;
    EXPECT_NEAR(measure(r.mMesh).mVolume, volume, 0.1 * volume);
}

TEST(MeshCad, CurvesThatTouchStaySelfIntersecting) {
    // The hole touches the rim: no split separates the loops, and the
    // splitting stops once the pieces at the touch are within the tolerance.
    const auto r = mio::mesh_cad_surface(brep_solids::disc_with_hole(2.0, 1.0, 1.0), size(1.1));
    ASSERT_EQ(r.mNumFailed, 1);
    EXPECT_EQ(r.mFaces[0].mStatus, MeshCadFaceStatus::SelfIntersecting);
    EXPECT_NE(r.mFaces[0].mReason.find("curves themselves meet"), std::string::npos)
        << r.mFaces[0].mReason;
    EXPECT_EQ(r.mNumTriangles, 0);
    EXPECT_LT(r.mNumPoints, 400);  // the halving stopped near the touch
}

// How many triangles use each undirected edge.
std::map<std::pair<std::int64_t, std::int64_t>, int> triangle_edge_uses(const mio::Mesh& rM) {
    std::map<std::pair<std::int64_t, std::int64_t>, int> n;
    const auto cells = rM.Cells(0);
    const mio::NDArray& c = cells.Conn();
    for (std::size_t t = 0; t < cells.NumCells(); ++t)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t a = mio::detail::read_int(c, 3 * t + k);
            const std::int64_t b = mio::detail::read_int(c, 3 * t + (k + 1) % 3);
            ++n[{std::min(a, b), std::max(a, b)}];
        }
    return n;
}

// The line block's segments of one B-rep edge (cad:edge), as point pairs.
std::vector<std::pair<std::int64_t, std::int64_t>> edge_segments(const mio::Mesh& rM,
                                                                 std::int64_t Edge) {
    std::vector<std::pair<std::int64_t, std::int64_t>> out;
    const auto lines = rM.Cells(1);
    const mio::NDArray& c = lines.Conn();
    const mio::NDArray& ids = rM.CellData("cad:edge", 1);
    for (std::size_t s = 0; s < lines.NumCells(); ++s)
        if (mio::detail::read_int(ids, s) == Edge) {
            const std::int64_t a = mio::detail::read_int(c, 2 * s);
            const std::int64_t b = mio::detail::read_int(c, 2 * s + 1);
            out.emplace_back(std::min(a, b), std::max(a, b));
        }
    return out;
}

bool any_warning_with(const MeshCadSurfaceResult& rR, const std::string& rText) {
    for (const std::string& w : rR.mWarnings)
        if (w.find(rText) != std::string::npos)
            return true;
    return false;
}

TEST(MeshCad, TwoPatchBeamMeshesItsEmbeddedCurve) {
    // Roadmap 6.1.1's probe: the two-patch beam's first face carries an
    // embedded edge (trim 4, a line across it that no BrepEdge names). It is
    // now constrained into the face's triangulation, so it is a run of mesh
    // edges, interior to the face, and the warning is gone.
    const std::filesystem::path input = std::filesystem::path(__FILE__).parent_path() / ".." /
                                        "python" / "meshes" / "cad_json" /
                                        "two_patch_beam.cad.json";
    ASSERT_TRUE(std::filesystem::exists(input)) << input;
    const mio::BrepInfo info = mio::read_cad_json_info(input.string());
    const auto r = mio::mesh_cad_surface(info, size(0.2));
    EXPECT_EQ(r.mNumFailed, 0);
    EXPECT_FALSE(any_warning_with(r, "embedded")) << r.mWarnings.front();
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    EXPECT_EQ(r.mNumMisorientedEdges, 0);
    // The same model without the embedded edge: the extra line cells are the
    // curve, and each is an edge of two triangles of the first face.
    mio::BrepInfo bare = info;
    for (mio::BrepFace& f : bare.mFaces)
        f.mEmbeddedEdges.clear();
    const auto b = mio::mesh_cad_surface(bare, size(0.2));
    const auto extra = edge_segments(r.mMesh, -1).size() - edge_segments(b.mMesh, -1).size();
    EXPECT_GE(extra, 5u);  // the line is about 3.0 long in a patch 5 long: h = 0.2
    const auto uses = triangle_edge_uses(r.mMesh);
    std::size_t interior = 0;
    for (const auto& seg : edge_segments(r.mMesh, -1)) {
        const auto it = uses.find(seg);
        ASSERT_TRUE(it != uses.end()) << seg.first << "-" << seg.second << " is not a mesh edge";
        interior += it->second == 2;
    }
    EXPECT_GE(interior, extra);
    EXPECT_NEAR(measure(r.mMesh).mArea, measure(b.mMesh).mArea, 1e-3 * measure(b.mMesh).mArea);
}

TEST(MeshCad, AnEmbeddedEdgeSharedWithAnotherFaceSharesItsPoints) {
    // A fin standing on a plate: the fin's bottom edge (BrepEdge 4) is an
    // embedded edge of the plate. Both faces use its points, so the fin is
    // attached to the plate along it: each of its segments is an edge of two
    // plate triangles and one fin triangle (a T, by design non-manifold).
    const auto r = mio::mesh_cad_surface(brep_solids::plate_with_fin(), size(0.1));
    ASSERT_EQ(r.mNumFailed, 0);
    EXPECT_TRUE(r.mWarnings.empty()) << r.mWarnings.front();
    const auto segs = edge_segments(r.mMesh, 4);
    ASSERT_GE(segs.size(), 6u);
    const auto uses = triangle_edge_uses(r.mMesh);
    for (const auto& seg : segs) {
        const auto it = uses.find(seg);
        ASSERT_TRUE(it != uses.end());
        EXPECT_EQ(it->second, 3) << seg.first << "-" << seg.second;
    }
    EXPECT_EQ(r.mNumNonManifoldEdges, static_cast<std::int64_t>(segs.size()));
    EXPECT_EQ(r.mNumMisorientedEdges, 0);
    EXPECT_NEAR(measure(r.mMesh).mArea, 1.0 + 0.6 * 0.5, 1e-12);
}

TEST(MeshCad, AClosedEmbeddedCurveIsMeshed) {
    // A circle inside a square plate, named by an edge: a closed curve with
    // no loop of its own, one corner, constrained into the plate.
    mio::BrepInfo info = brep_solids::square(4.0);
    mio::BrepTrim t;
    t.mTrimIndex = 100;
    t.mCurve = mio::detail::brep_curve_arc({2, 2, 0}, {1, 0, 0}, {0, 1, 0}, 1.0, 0.0, 2 * kPi);
    const std::int64_t e = brep_solids::embed(info, 0, t);
    const auto r = mio::mesh_cad_surface(info, size(0.25));
    ASSERT_EQ(r.mNumFailed, 0);
    EXPECT_TRUE(r.mWarnings.empty()) << r.mWarnings.front();
    const auto segs = edge_segments(r.mMesh, e);
    EXPECT_GE(segs.size(), 25u);
    const auto uses = triangle_edge_uses(r.mMesh);
    for (const auto& seg : segs) {
        const auto it = uses.find(seg);
        ASSERT_TRUE(it != uses.end());
        EXPECT_EQ(it->second, 2);
    }
    EXPECT_NEAR(measure(r.mMesh).mArea, 16.0, 1e-12);
    // Closed: as many segments as distinct points.
    std::set<std::int64_t> pts;
    for (const auto& seg : segs) {
        pts.insert(seg.first);
        pts.insert(seg.second);
    }
    EXPECT_EQ(pts.size(), segs.size());
}

TEST(MeshCad, AnEmbeddedEdgeThatCannotBeInsertedIsLeftOut) {
    // Outside the plate, ending on its boundary away from a corner, and
    // crossing it: each is left out with a warning naming it, and the plate
    // is still meshed (the face is not lost for a bad curve).
    for (const auto& ends : std::vector<std::array<double, 4>>{
             {5.0, 5.0, 6.0, 6.0}, {1.0, 1.0, 2.0, 0.0}, {1.0, 1.0, 2.0, -1.0}}) {
        mio::BrepInfo info = brep_solids::square(4.0);
        brep_solids::embed(info, 0,
                           brep_solids::param_line(ends[0], ends[1], ends[2], ends[3], 100));
        const auto r = mio::mesh_cad_surface(info, size(0.5));
        EXPECT_EQ(r.mNumFailed, 0);
        EXPECT_TRUE(any_warning_with(r, "embedded edge (trim 100) is not meshed"))
            << ends[2] << "," << ends[3];
        EXPECT_NEAR(measure(r.mMesh).mArea, 16.0, 1e-12);
        EXPECT_EQ(r.mNumNonManifoldEdges, 0);
    }
}

TEST(MeshCad, EmbeddedEdgesAreTheSameOnEveryThreadCount) {
    const mio::BrepInfo info = brep_solids::plate_with_fin();
    const std::uint64_t ref = digest(mio::mesh_cad_surface(info, size(0.05)).mMesh);
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
        EXPECT_EQ(digest(mio::mesh_cad_surface(info, size(0.05)).mMesh), ref) << threads;
    }
    omp_set_num_threads(saved);
#endif
    EXPECT_EQ(digest(mio::mesh_cad_surface(info, size(0.05)).mMesh), ref);
}

TEST(MeshCad, AFaceBeyondThePredicatesRangeIsReported) {
    // A box 1e40 on a side: every face's plane coordinates exceed 2^128, where
    // the exact predicates stop being exact, so every face is reported as
    // unsupported, by name, rather than triangulated wrong.
    const auto r = mio::mesh_cad_surface(brep_solids::box(1e40, 1e40, 1e40), size(1e40));
    EXPECT_EQ(r.mNumFailed, 6);
    for (const auto& f : r.mFaces) {
        EXPECT_EQ(f.mStatus, MeshCadFaceStatus::Unsupported);
        EXPECT_NE(f.mReason.find("exact predicates"), std::string::npos) << f.mReason;
    }
    EXPECT_EQ(r.mNumTriangles, 0);
}

TEST(MeshCad, ChordalDeviationCapsTheBoundary) {
    // With a target size larger than the cylinder, the rims are cut by the
    // sagitta alone: every boundary chord within the cap.
    const double dev = 1e-3;
    MeshCadSurfaceOptions o = size(10.0);
    o.mMaxChordalDeviation = dev;
    const auto r = mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), o);
    expect_watertight(r);
    const double* x = r.mMesh.Points().As<double>();
    const auto lines = r.mMesh.Cells(1);
    const mio::NDArray& c = lines.Conn();
    std::size_t rim = 0;
    for (std::size_t s = 0; s < lines.NumCells(); ++s) {
        const double* a = x + 3 * mio::detail::read_int(c, 2 * s);
        const double* b = x + 3 * mio::detail::read_int(c, 2 * s + 1);
        if (std::fabs(a[2] - b[2]) > 1e-12)
            continue;  // the seam
        ++rim;
        const double half = 0.5 * std::hypot(b[0] - a[0], b[1] - a[1]);
        EXPECT_LE(1.0 - std::sqrt(1.0 - half * half), dev * (1.0 + 1e-9));
    }
    EXPECT_GT(rim, 2u * 30u);
}

TEST(MeshCad, IsTheSameOnEveryThreadCount) {
    const mio::BrepInfo info = brep_solids::capped_cylinder(1.0, 2.0);
    const std::uint64_t ref = digest(mio::mesh_cad_surface(info, size(0.1)).mMesh);
    EXPECT_EQ(digest(mio::mesh_cad_surface(info, size(0.1)).mMesh), ref);
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
        EXPECT_EQ(digest(mio::mesh_cad_surface(info, size(0.1)).mMesh), ref) << threads;
    }
    omp_set_num_threads(saved);
#endif
}

TEST(MeshCad, RoundTripsThroughVtuAndGmsh) {
    const auto r = mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), size(0.25));
    const std::string vtu = mt::temp_path(".vtu");
    mio::write_vtu(vtu, r.mMesh, true, false);
    const mio::Mesh v = mio::read_vtu(vtu);
    const std::string msh = mt::temp_path(".msh");
    mio::write_gmsh41(msh, r.mMesh, false);
    const mio::Mesh g = mio::read_gmsh(msh);
    for (const mio::Mesh* p_m : {&v, &g}) {
        ASSERT_EQ(p_m->NumPoints(), r.mMesh.NumPoints());
        ASSERT_GE(p_m->NumCellBlocks(), 1u);
        const double* a = p_m->Points().As<double>();
        const double* b = r.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < 3 * r.mMesh.NumPoints(); ++i)
            ASSERT_EQ(a[i], b[i]);
        std::size_t tris = 0;
        for (std::size_t k = 0; k < p_m->NumCellBlocks(); ++k)
            if (p_m->Cells(k).Type() == "triangle")
                tris += p_m->Cells(k).NumCells();
        EXPECT_EQ(tris, r.mMesh.Cells(0).NumCells());
    }
    ASSERT_TRUE(v.HasCellData(mio::kMeshCadFaceName));
    const mio::NDArray& ca = v.CellData(mio::kMeshCadFaceName, 0);
    const mio::NDArray& cb = r.mMesh.CellData(mio::kMeshCadFaceName, 0);
    ASSERT_EQ(ca.Size(), cb.Size());
    for (std::size_t i = 0; i < cb.Size(); ++i)
        EXPECT_EQ(mio::detail::read_int(ca, i), mio::detail::read_int(cb, i));
}

// Roadmap 6.1.2: every face of the one triangle block is a `face:<k>` region,
// and Gmsh 4.1 places a physical group on an entity, so the writer splits the
// block into one per face; the box opens with six named physical surfaces.
TEST(MeshCad, GmshFourOneGivesEachFaceItsPhysicalGroup) {
    const auto r = mio::mesh_cad_surface(brep_solids::box(), size(0.25));
    ASSERT_EQ(r.mNumFailed, 0);
    const mio::NDArray& faces = r.mMesh.CellData(mio::kMeshCadFaceName, 0);
    std::map<std::int64_t, std::size_t> per_face;
    for (std::size_t i = 0; i < faces.Size(); ++i)
        ++per_face[mio::detail::read_int(faces, i)];
    ASSERT_EQ(per_face.size(), 6u);
    for (const bool binary : {false, true}) {
        const std::string msh = mt::temp_path(".msh");
        mio::write_gmsh41(msh, r.mMesh, binary);
        const mio::Mesh g = mio::read_gmsh(msh);
        std::size_t named = 0;
        for (std::size_t i = 0; i < g.NumRegions(); ++i) {
            const mio::Region& reg = g.Region(i);
            if (reg.mName.rfind("face:", 0) != 0)
                continue;
            ++named;
            EXPECT_EQ(reg.mKind, mio::RegionKind::Cell);
            EXPECT_EQ(reg.mDim, 2);
            const std::int64_t k = std::stoll(reg.mName.substr(5));
            EXPECT_EQ(reg.NumEntries(), per_face.at(k)) << reg.mName << " binary " << binary;
        }
        EXPECT_EQ(named, 6u) << "binary " << binary;
        // The triangles still carry their face, now one block per face.
        std::size_t tris = 0;
        for (std::size_t b = 0; b < g.NumCellBlocks(); ++b) {
            if (g.Cells(b).Type() != "triangle")
                continue;
            tris += g.Cells(b).NumCells();
            const mio::NDArray& f = g.CellData(mio::kMeshCadFaceName, b);
            for (std::size_t i = 1; i < f.Size(); ++i)
                ASSERT_EQ(mio::detail::read_int(f, i), mio::detail::read_int(f, 0));
        }
        EXPECT_EQ(tris, r.mMesh.Cells(0).NumCells());
    }
}

// Roadmap 6.1.1: the pipeline carries the B-rep a CAD read returns from Input
// to the steps that mesh it, and drops it at the first step that does not.
std::string pipeline_cube_path() {
    const std::filesystem::path input = std::filesystem::path(__FILE__).parent_path() / ".." /
                                        "python" / "meshes" / "cad_json" / "unit_cube.cad.json";
    return input.string();
}

mio::PipelineStep pipeline_step(std::string Op,
                                std::map<std::string, mio::PipelineValue> Params = {}) {
    mio::PipelineStep s;
    s.mOp = std::move(Op);
    s.mParams = std::move(Params);
    return s;
}

TEST(MeshCad, PipelineMeshesTheModelItsInputReads) {
    const std::string input = pipeline_cube_path();
    ASSERT_TRUE(std::filesystem::exists(input)) << input;
    mio::Pipeline pipe;
    pipe.mInput.mPath = input;
    pipe.mOutput.mPath = mt::temp_path(".vtu");
    pipe.mSteps = {pipeline_step("MeshCadSurface", {{"TargetSize", 0.25}, {"MinAngle", 25.0}})};
    for (const int version : {1, 2}) {
        pipe.mVersion = version;
        const mio::PipelineReport report = mio::run_pipeline(pipe);
        ASSERT_EQ(report.mSteps.size(), 1u);
        MeshCadSurfaceOptions o = size(0.25);
        o.mMinAngleDeg = 25.0;
        const auto direct = mio::mesh_cad_surface(mio::read_cad_json_info(input), o);
        const mio::Mesh out = mio::read_vtu(pipe.mOutput.mPath);
        ASSERT_EQ(out.NumPoints(), direct.mMesh.NumPoints());
        const double* a = out.Points().As<double>();
        const double* b = direct.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < 3 * out.NumPoints(); ++i)
            ASSERT_EQ(a[i], b[i]);
        const auto& counters = report.mSteps[0].mCounters;
        ASSERT_FALSE(counters.empty());
        EXPECT_EQ(counters[0].first, "NumPoints");
        EXPECT_EQ(counters[0].second, static_cast<double>(direct.mNumPoints));
    }
}

TEST(MeshCad, PipelineRefinesThenMeshesTheRefinedModel) {
    mio::Pipeline pipe;
    pipe.mVersion = 2;
    pipe.mInput.mPath = pipeline_cube_path();
    pipe.mOutput.mPath = mt::temp_path(".vtu");
    pipe.mSteps = {
        pipeline_step("IgaRefine", {{"Subdivisions", std::vector<std::int64_t>{2, 3}}}),
        pipeline_step("MeshCadSurface", {{"TargetSize", 0.3}, {"Order", std::int64_t(2)}})};
    const mio::PipelineReport report = mio::run_pipeline(pipe);
    ASSERT_EQ(report.mSteps.size(), 2u);
    EXPECT_EQ(report.mSteps[0].mOp, "IgaRefine");
    const mio::Mesh out = mio::read_vtu(pipe.mOutput.mPath);
    EXPECT_EQ(out.Cells(0).Type(), "triangle6");
}

TEST(MeshCad, PipelineModelStepsFailByNameWithoutAModel) {
    const auto message = [](const mio::Pipeline& rPipe) {
        try {
            mio::run_pipeline(rPipe);
        } catch (const std::invalid_argument& r_e) {
            return std::string(r_e.what());
        }
        return std::string();
    };
    mio::Pipeline pipe;
    pipe.mInput.mPath = pipeline_cube_path();
    pipe.mOutput.mPath = mt::temp_path(".vtu");
    // A step that changes the mesh drops the model.
    pipe.mSteps = {pipeline_step("Quality"),
                   pipeline_step("MeshCadSurface", {{"TargetSize", 0.25}})};
    EXPECT_NE(message(pipe).find("step 1 'Quality' changed the mesh"), std::string::npos)
        << message(pipe);
    // A format whose read carries none.
    const std::string vtu = mt::temp_path(".vtu");
    mio::write_vtu(vtu, mio::mesh_cad_surface(brep_solids::box(), size(0.5)).mMesh, true, false);
    pipe.mInput.mPath = vtu;
    pipe.mSteps = {pipeline_step("IgaSample")};
    EXPECT_NE(message(pipe).find("Input format 'vtu' carries none"), std::string::npos)
        << message(pipe);
    // TargetSize is required.
    pipe.mInput.mPath = pipeline_cube_path();
    pipe.mSteps = {pipeline_step("MeshCadSurface")};
    EXPECT_NE(message(pipe).find("'TargetSize' is required"), std::string::npos);
    // A step applied on its own has no run to carry a model.
    mio::PipelineReport report;
    EXPECT_THROW(mio::apply_pipeline_step(
                     mio::Mesh(), pipeline_step("MeshCadSurface", {{"TargetSize", 0.25}}), report),
                 std::invalid_argument);
}

// Roadmap 6.7.1: a closed `mesh_cad_surface` result is a skin `octree_mesh`
// fills with tetrahedra. Under uniform skin grading at the surface's own size
// the enclosed volume follows the skin's at second order, the boundary stays
// within a few hundredths of `h` of it and the dihedral angles keep inside
// Labelle and Shewchuk's 10.7-164.8 degree bound.
double tet_volume(const mio::Mesh& rM) {
    const double* x = rM.Points().As<double>();
    double v = 0.0;
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
        const auto cells = rM.Cells(b);
        if (cells.Type() != "tetra")
            continue;
        const mio::NDArray& c = cells.Conn();
        for (std::size_t t = 0; t < cells.NumCells(); ++t) {
            const double* p[4];
            for (int k = 0; k < 4; ++k)
                p[k] = x + 3 * mio::detail::read_int(c, 4 * t + k);
            const double u[3] = {p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]};
            const double w[3] = {p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]};
            const double z[3] = {p[3][0] - p[0][0], p[3][1] - p[0][1], p[3][2] - p[0][2]};
            v += std::abs(u[0] * (w[1] * z[2] - w[2] * z[1]) - u[1] * (w[0] * z[2] - w[2] * z[0]) +
                          u[2] * (w[0] * z[1] - w[1] * z[0])) /
                 6.0;
        }
    }
    return v;
}

mio::OctreeMeshResult fill(const mio::Mesh& rSkin, double H) {
    mio::OctreeMeshOptions o;
    o.mSkinGrading = mio::OctreeMeshSkinGrading::Uniform;
    o.mSurfaceSize = H;
    return mio::octree_mesh(rSkin, o);
}

TEST(MeshCad, OctreeMeshFillsTheSkinWithTetrahedra) {
    MeshCadSurfaceOptions o = size(0.2);
    o.mMinAngleDeg = 25.0;
    std::vector<double> errors;
    for (const double h : {0.2, 0.1}) {
        o.mTargetSize = h;
        const auto skin = mio::mesh_cad_surface(brep_solids::sphere(), o);
        ASSERT_EQ(skin.mNumOpenEdges, 0);
        const auto vol = fill(skin.mMesh, h);
        const double v = tet_volume(vol.mMesh);
        const double skin_volume = measure(skin.mMesh).mVolume;
        errors.push_back(std::abs(v - skin_volume));
        // Measured 0.30% and 0.075% of the sphere's volume.
        EXPECT_LT(errors.back(), 0.005 * skin_volume) << h;
        // Measured 0.036 h and 0.028 h.
        EXPECT_GE(vol.mHausdorffToSkin, 0.0);
        EXPECT_LT(vol.mHausdorffToSkin, 0.06 * h) << h;
        // Measured 25.6-132.4 and 21.6-144.2 degrees.
        EXPECT_GT(vol.mMinDihedral, 10.7) << h;
        EXPECT_LT(vol.mMaxDihedral, 164.8) << h;
        // And the sphere's own volume to the skin's chordal error.
        EXPECT_LT(std::abs(v - 4.0 * kPi / 3.0),
                  2.0 * std::abs(skin_volume - 4.0 * kPi / 3.0) + errors.back())
            << h;
    }
    EXPECT_GT(std::log2(errors[0] / errors[1]), 1.7);  // measured 2.0
    // A box is filled exactly, by BCC tetrahedra only.
    o.mTargetSize = 0.2;
    const auto box = mio::mesh_cad_surface(brep_solids::box(), o);
    const auto vol = fill(box.mMesh, 0.2);
    EXPECT_NEAR(tet_volume(vol.mMesh), 1.0, 1e-12);
    EXPECT_GT(vol.mMinDihedral, 35.0);
    EXPECT_LT(vol.mMaxDihedral, 126.0);
    // The curved elements are skins only once linear: the octree refuses them.
    o.mOrder = 2;
    EXPECT_THROW(fill(mio::mesh_cad_surface(brep_solids::box(), o).mMesh, 0.2), std::exception);
}

TEST(MeshCad, Refusals) {
    const mio::BrepInfo box = brep_solids::box();
    for (const double h : {0.0, -1.0, std::nan("")})
        EXPECT_THROW(mio::mesh_cad_surface(box, size(h)), std::invalid_argument) << h;
    MeshCadSurfaceOptions o = size(0.1);
    o.mMaxChordalDeviation = -1.0;
    EXPECT_THROW(mio::mesh_cad_surface(box, o), std::invalid_argument);
    o = size(0.1);
    o.mMinPointsPerEdge = -1;
    EXPECT_THROW(mio::mesh_cad_surface(box, o), std::invalid_argument);
    EXPECT_THROW(mio::mesh_cad_surface(mio::BrepInfo{}, size(0.1)), std::invalid_argument);
    o = size(1e-4);
    o.mMaxPoints = 1000;
    try {
        const auto r = mio::mesh_cad_surface(box, o);
        EXPECT_EQ(r.mNumFailed, 6);
    } catch (const std::length_error&) {
        SUCCEED();
    }
}

TEST(TriangulatePolygon2d, SquareWithAHole) {
    const std::vector<std::vector<std::array<double, 2>>> loops = {
        {{0, 0}, {4, 0}, {4, 4}, {0, 4}}, {{1, 1}, {1, 3}, {3, 3}, {3, 1}}};
    const mio::Mesh m = mio::triangulate_polygon_2d(loops);
    EXPECT_EQ(m.NumPoints(), 8u);
    EXPECT_EQ(m.Cells(0).NumCells(), 8u);
    EXPECT_NEAR(measure(m).mArea, 12.0, 1e-12);
    const mio::Mesh fine = mio::triangulate_polygon_2d(loops, 0.25);
    EXPECT_NEAR(measure(fine).mArea, 12.0, 1e-12);
    EXPECT_GT(fine.NumPoints(), 100u);
    EXPECT_THROW(mio::triangulate_polygon_2d({{{0, 0}, {1, 1}, {1, 0}, {0, 1}}}),
                 std::invalid_argument);
    EXPECT_THROW(mio::triangulate_polygon_2d(loops, -1.0), std::invalid_argument);
}

#ifdef MESHIOPLUSPLUS_HAS_OCCT
// The closed unit cube as six outward quads.
mio::Mesh unit_cube_skin() {
    return mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "quad",
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}});
}

TEST(MeshCadOcct, StepCubeIsWatertightAndMatchesBRepMesh) {
    const std::string path = mt::temp_path(".step");
    mio::write_step(path, unit_cube_skin());
    mio::BrepInfo info;
    const mio::Mesh occt = mio::read_step(path, info);
    ASSERT_EQ(info.mFaces.size(), 6u);
    for (const double h : {0.25, 0.1}) {
        const auto r = mio::mesh_cad_surface(info, size(h));
        expect_watertight(r);
        const Surface s = measure(r.mMesh);
        EXPECT_NEAR(s.mArea, 6.0, 1e-9);
        EXPECT_NEAR(std::fabs(s.mVolume), 1.0, 1e-9);
        EXPECT_GT(s.mVolume, 0.0) << "outward";
        EXPECT_NEAR(s.mArea, measure(occt).mArea, 1e-9);
        EXPECT_LE(r.mMaxEdge, 1.5 * h);
    }
}
#endif

}  // namespace

namespace {

// One untrimmed bilinear patch, the square [0, 1]^2 of the plane z = 0, its
// knots in u starting at U0.
mio::BrepInfo mcadt_square_from(double U0) {
    mio::BrepInfo b;
    mio::NurbsPatch p;
    p.mDegrees = {1, 1};
    p.mKnots = {{U0, U0, 1, 1}, {0, 0, 1, 1}};
    p.mControlPoints = mio::NDArray(mio::DType::Float64, {4, 3});
    p.mWeights = mio::NDArray(mio::DType::Float64, {4});
    const double pts[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 3; ++k)
            p.mControlPoints.As<double>()[3 * i + k] = pts[i][k];
        p.mWeights.As<double>()[i] = 1.0;
    }
    b.mPatches.push_back(p);
    mio::BrepFace f;
    f.mId = 1;
    f.mPatch = 0;
    b.mFaces.push_back(f);
    return b;
}

}  // namespace

// Roadmap 6.5.7: real STEP files write parameters a hair off zero (7e-44
// where 0 was meant), under the 2^-128 the exact predicates hold for; 95
// faces of the corpus failed on it. Such a face is triangulated with those
// coordinates flushed to zero; a face whose parameters are in range is
// meshed on exactly the coordinates it always was.
TEST(MeshCad, ParametersAHairOffZeroAreMeshed) {
    mio::MeshCadSurfaceOptions o;
    o.mTargetSize = 0.25;
    const auto clean = mio::mesh_cad_surface(mcadt_square_from(0.0), o);
    for (const double u0 : {7.3e-44, -2.1e-45, 3.2e-48}) {
        const auto r = mio::mesh_cad_surface(mcadt_square_from(u0), o);
        ASSERT_EQ(r.mNumFailed, 0) << u0 << ": " << r.mFaces[0].mReason;
        EXPECT_EQ(r.mNumTriangles, clean.mNumTriangles) << u0;
        EXPECT_EQ(r.mNumOpenEdges, clean.mNumOpenEdges) << u0;
    }
    // In range, untouched.
    const auto r = mio::mesh_cad_surface(mcadt_square_from(1e-30), o);
    EXPECT_EQ(r.mNumFailed, 0);
}
