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
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
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

TEST(MeshCad, TrimmedPlateWithAHole) {
    // tests/cpp/test_iga.cpp's plate: a 4 x 4 plane with a diamond hole, a
    // single degree-1 trim whose corners are knots (kept as samples).
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
    const mio::BrepInfo info = mio::read_cad_json_info(path);
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

TEST(MeshCad, AClosedEdgeStartedElsewhereIsReported) {
    // The top cap's rim starts a quarter turn away from the band's: the
    // topology joins two points that are not one corner.
    mio::BrepInfo info = brep_solids::capped_cylinder(1.0, 2.0, /*TopCapStart=*/0.5 * kPi);
    const auto r = mio::mesh_cad_surface(info, size(0.2));
    EXPECT_GE(r.mNumFailed, 1);
    bool reported = false;
    for (const auto& f : r.mFaces)
        if (f.mStatus != MeshCadFaceStatus::Ok) {
            reported = true;
            EXPECT_FALSE(f.mReason.empty());
        }
    EXPECT_TRUE(reported);
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
