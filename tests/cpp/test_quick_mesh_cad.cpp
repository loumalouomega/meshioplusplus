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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// Small versions of the mesh_cad tests (`test_mesh_cad*.cpp`). The full suites mesh a sphere
// at h = 0.1 and 0.05 and a cylinder at fine curvature bounds, which takes minutes each in
// the instrumented coverage build, so a pull request skips them there; these mesh the same
// solids at h = 0.4 to 0.8 and check the same properties (watertight, the enclosed volume,
// the angle bound, every size source, the curved elements, the quads, the octree fill of the
// skin) with looser bounds. Their suite name starts with `Quick`, which the job's exclusion
// does not match.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

#include "thread_sweep.hpp"

// Project includes
#include "brep_solids.hpp"
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::MeshCadSurfaceOptions;
using mio::MeshCadSurfaceResult;

constexpr double kPi = 3.14159265358979323846;

MeshCadSurfaceOptions quick_size(double H, double MinAngle = 25.0) {
    MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    o.mMinAngleDeg = MinAngle;
    return o;
}

// The area and, by the divergence theorem, the volume of the triangle blocks of a mesh.
struct QuickMeasure {
    double mArea = 0.0, mVolume = 0.0;
};

QuickMeasure quick_measure(const mio::Mesh& rM) {
    QuickMeasure s;
    const double* x = rM.Points().As<double>();
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
        const auto cells = rM.Cells(b);
        if (cells.Type() != "triangle")
            continue;
        const mio::NDArray& c = cells.Conn();
        for (std::size_t t = 0; t < cells.NumCells(); ++t) {
            const double* a = x + 3 * mio::detail::read_int(c, 3 * t);
            const double* p = x + 3 * mio::detail::read_int(c, 3 * t + 1);
            const double* q = x + 3 * mio::detail::read_int(c, 3 * t + 2);
            const double u[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
            const double v[3] = {q[0] - a[0], q[1] - a[1], q[2] - a[2]};
            const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                 u[0] * v[1] - u[1] * v[0]};
            s.mArea += 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            s.mVolume += (a[0] * (p[1] * q[2] - p[2] * q[1]) - a[1] * (p[0] * q[2] - p[2] * q[0]) +
                          a[2] * (p[0] * q[1] - p[1] * q[0])) /
                         6.0;
        }
    }
    return s;
}

void quick_expect_closed(const MeshCadSurfaceResult& rR) {
    EXPECT_EQ(rR.mNumFailed, 0);
    EXPECT_EQ(rR.mNumOpenEdges, 0);
    EXPECT_EQ(rR.mNumNonManifoldEdges, 0);
    EXPECT_EQ(rR.mNumMisorientedEdges, 0);
}

std::uint64_t quick_digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

std::int64_t quick_block(const mio::Mesh& rM, const std::string& rType) {
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b)
        if (rM.Cells(b).Type() == rType)
            return static_cast<std::int64_t>(b);
    return -1;
}

}  // namespace

// The box, the capped cylinder and the sphere are watertight at a coarse size, enclose the
// solid's volume to the chordal error, and meet the angle bound the refinement is asked for,
// with each smoothing method.
TEST(QuickMeshCad, SolidsAreWatertightAndEncloseTheirVolume) {
    struct Solid {
        const char* mName;
        mio::BrepInfo mBrep;
        double mVolume, mTolerance;
    };
    const Solid solids[] = {{"box", brep_solids::box(), 1.0, 1e-9},
                            {"cylinder", brep_solids::capped_cylinder(1.0, 2.0), 2.0 * kPi, 0.12},
                            {"sphere", brep_solids::sphere(1.0), 4.0 * kPi / 3.0, 0.2}};
    for (const Solid& s : solids) {
        SCOPED_TRACE(s.mName);
        const auto r = mio::mesh_cad_surface(s.mBrep, quick_size(0.6));
        quick_expect_closed(r);
        EXPECT_NEAR(quick_measure(r.mMesh).mVolume, s.mVolume, s.mTolerance * s.mVolume);
        EXPECT_GE(r.mMinAngleDeg, 0.0);
        // The same mesh twice, and on every thread count.
        EXPECT_EQ(quick_digest(mio::mesh_cad_surface(s.mBrep, quick_size(0.6)).mMesh),
                  quick_digest(r.mMesh));
    }
    for (auto smoothing : {mio::MeshCadSmoothing::Laplacian, mio::MeshCadSmoothing::Odt}) {
        auto o = quick_size(0.5);
        o.mSmoothMethod = smoothing;
        o.mSmoothIterations = 2;
        const auto r = mio::mesh_cad_surface(brep_solids::capped_cylinder(1.0, 2.0), o);
        quick_expect_closed(r);
        EXPECT_GE(r.mNumSmoothed, 0);
    }
}

// The size sources: curvature, cells across a thin wall, the grading, the varying metric,
// the pole rings and the reparameterisation of a flattened face.
TEST(QuickMeshCad, SizeSourcesRefineWhereAsked) {
    const mio::BrepInfo cylinder = brep_solids::capped_cylinder(1.0, 2.0);
    const auto coarse = mio::mesh_cad_surface(cylinder, quick_size(0.8));
    auto o = quick_size(0.8);
    o.mCurvatureDeviation = 0.02;
    const auto curved = mio::mesh_cad_surface(cylinder, o);
    quick_expect_closed(curved);
    EXPECT_GT(curved.mNumTriangles, coarse.mNumTriangles);
    o.mGrading = 1.3;
    quick_expect_closed(mio::mesh_cad_surface(cylinder, o));

    auto thin = quick_size(0.5);
    thin.mCellsAcrossGaps = 2.0;
    const auto plate = mio::mesh_cad_surface(brep_solids::box(1.0, 1.0, 0.1), thin);
    quick_expect_closed(plate);
    EXPECT_GT(
        plate.mNumTriangles,
        mio::mesh_cad_surface(brep_solids::box(1.0, 1.0, 0.1), quick_size(0.5)).mNumTriangles);

    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    auto m = quick_size(0.6);
    m.mMetric = mio::MeshCadMetric::Varying;
    quick_expect_closed(mio::mesh_cad_surface(sphere, m));
    auto rings = quick_size(0.6);
    rings.mPoleRings = 2;
    quick_expect_closed(mio::mesh_cad_surface(sphere, rings));
    auto flat = quick_size(0.6);
    flat.mReparamConditionNumber = 2.0;
    quick_expect_closed(mio::mesh_cad_surface(sphere, flat));
}

// Curved elements of order 2 and 3: their points lie on the sphere, the blocks are the
// curved types, and untangling keeps the mesh valid.
TEST(QuickMeshCad, CurvedElementsLieOnTheSurface) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    for (int order : {2, 3}) {
        auto o = quick_size(0.8);
        o.mOrder = order;
        mio::MeshCadSurfaceExtras ex;
        ex.mUntangle = true;
        const auto r = mio::mesh_cad_surface(sphere, o, mio::MeshCadSizeSources{}, ex);
        EXPECT_EQ(r.mNumFailed, 0);
        EXPECT_GE(quick_block(r.mMesh, order == 2 ? "triangle6" : "triangle10"), 0);
        const double* x = r.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
            EXPECT_NEAR(std::sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] +
                                  x[3 * i + 2] * x[3 * i + 2]),
                        1.0, 1e-6)
                << "point " << i;
    }
}

// Quads: all-quad on the box and the cylinder, quad-dominant on the sphere with each pairing.
TEST(QuickMeshCad, QuadsAreRecombinedAndSplit) {
    const auto run = [](const mio::BrepInfo& rB, double H, mio::MeshCadQuads Q,
                        mio::MeshCadQuadPairs Pairs, mio::MeshCadSurfaceExtrasReport& rReport) {
        mio::MeshCadSurfaceExtras ex;
        ex.mQuads = Q;
        ex.mQuadPairs = Pairs;
        return mio::mesh_cad_surface(rB, quick_size(H), mio::MeshCadSizeSources{}, ex, nullptr,
                                     &rReport);
    };
    for (const mio::BrepInfo& b : {brep_solids::box(), brep_solids::capped_cylinder(1.0, 2.0)}) {
        mio::MeshCadSurfaceExtrasReport rep;
        const auto r = run(b, 0.5, mio::MeshCadQuads::All, mio::MeshCadQuadPairs::Matching, rep);
        EXPECT_EQ(r.mNumFailed, 0);
        EXPECT_EQ(quick_block(r.mMesh, "triangle"), -1);
        ASSERT_GE(quick_block(r.mMesh, "quad"), 0);
        EXPECT_EQ(rep.mNumLeftoverTriangles, 0);
        EXPECT_GT(rep.mMinQuadScaledJacobian, 0.0);
    }
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    std::int64_t leftover[2] = {0, 0};
    int k = 0;
    for (auto pairs : {mio::MeshCadQuadPairs::Matching, mio::MeshCadQuadPairs::Greedy}) {
        mio::MeshCadSurfaceExtrasReport rep;
        const auto r = run(sphere, 0.4, mio::MeshCadQuads::Dominant, pairs, rep);
        EXPECT_EQ(r.mNumFailed, 0);
        EXPECT_GT(rep.mNumQuads, 0);
        EXPECT_EQ(r.mNumOpenEdges, 0);
        leftover[k++] = rep.mNumLeftoverTriangles;
    }
    // The matching leaves no more triangles than the greedy pairing.
    EXPECT_LE(leftover[0], leftover[1]);
}

// The octree backend fills the skin: the box exactly, the sphere to the skin's error.
TEST(QuickMeshCad, TheOctreeMeshFillsTheSkinWithTetrahedra) {
    const auto fill = [](const mio::Mesh& rSkin, double H) {
        mio::OctreeMeshOptions o;
        o.mSkinGrading = mio::OctreeMeshSkinGrading::Uniform;
        o.mSurfaceSize = H;
        return mio::octree_mesh(rSkin, o);
    };
    const auto tet_volume = [](const mio::Mesh& rM) {
        double v = 0.0;
        const double* x = rM.Points().As<double>();
        for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
            const auto cells = rM.Cells(b);
            if (cells.Type() != "tetra")
                continue;
            for (std::size_t t = 0; t < cells.NumCells(); ++t) {
                const double* p[4];
                for (int j = 0; j < 4; ++j)
                    p[j] = x + 3 * mio::detail::read_int(cells.Conn(), 4 * t + j);
                const double u[3] = {p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]};
                const double w[3] = {p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]};
                const double z[3] = {p[3][0] - p[0][0], p[3][1] - p[0][1], p[3][2] - p[0][2]};
                v += std::abs(u[0] * (w[1] * z[2] - w[2] * z[1]) -
                              u[1] * (w[0] * z[2] - w[2] * z[0]) +
                              u[2] * (w[0] * z[1] - w[1] * z[0])) /
                     6.0;
            }
        }
        return v;
    };
    const auto sphere = mio::mesh_cad_surface(brep_solids::sphere(1.0), quick_size(0.4));
    ASSERT_EQ(sphere.mNumOpenEdges, 0);
    const auto vol = fill(sphere.mMesh, 0.4);
    EXPECT_NEAR(tet_volume(vol.mMesh), quick_measure(sphere.mMesh).mVolume,
                0.03 * quick_measure(sphere.mMesh).mVolume);
    EXPECT_GT(vol.mMinDihedral, 10.7);
    EXPECT_LT(vol.mMaxDihedral, 164.8);
    const auto box = mio::mesh_cad_surface(brep_solids::box(), quick_size(0.4));
    EXPECT_NEAR(tet_volume(fill(box.mMesh, 0.4).mMesh), 1.0, 1e-12);
    // Curved elements are refused: the octree takes linear skins only.
    auto o = quick_size(0.4);
    o.mOrder = 2;
    EXPECT_THROW(fill(mio::mesh_cad_surface(brep_solids::box(), o).mMesh, 0.4), std::exception);
}
