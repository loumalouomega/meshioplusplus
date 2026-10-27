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
// project_to_cad (roadmap 6.10.1): a mesh's points moved back onto the CAD
// faces and edges its cells name.

// System includes
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

// Project includes
#include "brep_solids.hpp"
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/project_to_cad.hpp"

namespace {

namespace mio = meshioplusplus;

std::uint64_t digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

mio::Mesh surface(const mio::BrepInfo& rB, double H) {
    mio::MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    return mio::mesh_cad_surface(rB, o).mMesh;
}

/// The mesh with every point moved by `rMove`.
template <class TMove>
mio::Mesh moved(const mio::Mesh& rM, TMove&& rMove) {
    mio::Mesh out = mio::detail::clone_mesh(rM);
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {out.NumPoints(), std::size_t{3}});
    const double* src = rM.Points().As<double>();
    for (std::size_t i = 0; i < out.NumPoints(); ++i) {
        std::array<double, 3> x{src[3 * i], src[3 * i + 1], src[3 * i + 2]};
        rMove(i, x);
        for (std::size_t k = 0; k < 3; ++k)
            p.As<double>()[3 * i + k] = x[k];
    }
    out.AssignPoints(std::move(p));
    return out;
}

}  // namespace

// The sphere's mesh shrunk by 3% and jittered: every point goes back onto
// the unit sphere, none is left, and the result does not depend on threads.
TEST(ProjectToCad, TheShrunkSphereGoesBack) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    const mio::Mesh m = surface(sphere, 0.2);
    const mio::Mesh off = moved(m, [](std::size_t I, std::array<double, 3>& rX) {
        const double s = 0.97 + 0.01 * std::sin(double(I));
        for (double& r_c : rX)
            r_c *= s;
    });
    const auto r = mio::project_to_cad(off, sphere);
    std::printf("[6.10.1] sphere: %lld projected (%lld on edges), largest move %.3e\n",
                static_cast<long long>(r.mNumProjected), static_cast<long long>(r.mNumOnEdges),
                r.mMaxMove);
    EXPECT_EQ(r.mNumProjected, static_cast<std::int64_t>(m.NumPoints()));
    EXPECT_EQ(r.mNumAmbiguous + r.mNumFailed + r.mNumTooFar + r.mNumOutsideTrim, 0);
    // On the sphere; a point on an edge goes onto the edge's own 3-D curve
    // (6.10.1), which the fixture's seams fit within 1e-8 of it, as a model's
    // curves meet its surfaces within the model's tolerance.
    const double* p = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
        EXPECT_NEAR(std::hypot(p[3 * i], p[3 * i + 1], p[3 * i + 2]), 1.0, 1e-8) << i;
    EXPECT_GT(r.mMaxMove, 0.01);
    EXPECT_EQ(r.mMesh.NumCellBlocks(), m.NumCellBlocks());
    // Points already on it stay, to round-off, those on the seams to within
    // the fit of their curves.
    const auto again = mio::project_to_cad(m, sphere);
    EXPECT_LT(again.mMaxMove, 1e-8);
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const std::uint64_t ref = digest(r.mMesh);
    const int saved = omp_get_max_threads();
    for (int threads : {1, 3, 8}) {
        omp_set_num_threads(threads);
        EXPECT_EQ(digest(mio::project_to_cad(off, sphere).mMesh), ref) << threads;
    }
    omp_set_num_threads(saved);
#endif
}

// The capped cylinder moved off: the rim points (on the circles' lines) go
// back onto their circle, the side's onto the cylinder, the caps' onto
// their planes; a corner-free rim has no corner.
TEST(ProjectToCad, TheCylinderKeepsItsRims) {
    const mio::BrepInfo cyl = brep_solids::capped_cylinder(1.0, 2.0);
    const mio::Mesh m = surface(cyl, 0.25);
    const mio::Mesh off = moved(m, [](std::size_t I, std::array<double, 3>& rX) {
        rX[0] += 0.02 * std::cos(3.0 * double(I));
        rX[1] += 0.02 * std::sin(5.0 * double(I));
        rX[2] += 0.02 * std::sin(7.0 * double(I));
    });
    const auto r = mio::project_to_cad(off, cyl);
    EXPECT_EQ(r.mNumAmbiguous + r.mNumFailed, 0);
    EXPECT_GT(r.mNumOnEdges, 0);
    const double* p = r.mMesh.Points().As<double>();
    const double* p0 = m.Points().As<double>();
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        const double rr = std::hypot(p[3 * i], p[3 * i + 1]);
        const double z = p[3 * i + 2];
        const bool on_side = std::fabs(rr - 1.0) < 1e-10 && z > -1e-10 && z < 2.0 + 1e-10;
        const bool on_cap =
            (std::fabs(z) < 1e-10 || std::fabs(z - 2.0) < 1e-10) && rr < 1.0 + 1e-10;
        EXPECT_TRUE(on_side || on_cap) << i;
        // A rim point stays on its rim.
        const double z0 = p0[3 * i + 2];
        if (std::fabs(std::hypot(p0[3 * i], p0[3 * i + 1]) - 1.0) < 1e-12 &&
            (std::fabs(z0) < 1e-12 || std::fabs(z0 - 2.0) < 1e-12)) {
            EXPECT_NEAR(rr, 1.0, 1e-10) << i;
            EXPECT_NEAR(z, z0, 1e-10) << i;
        }
    }
}

// Roadmap 6.10.1: an edge's own 3-D curve, where the model gives one that
// its trace follows, takes the edge's points: every rim point goes onto it.
// A curve the trace leaves (the rim's moved up 0.1) is warned and the trace
// kept, so the points still land on the rim.
TEST(ProjectToCad, EdgePointsGoOntoTheEdgesCurves) {
    mio::BrepInfo cyl = brep_solids::capped_cylinder(1.0, 2.0);
    const mio::Mesh m = surface(cyl, 0.25);
    const mio::Mesh off = moved(m, [](std::size_t I, std::array<double, 3>& rX) {
        rX[0] += 0.01 * std::cos(3.0 * double(I));
        rX[2] += 0.01 * std::sin(7.0 * double(I));
    });
    const auto r = mio::project_to_cad(off, cyl);
    std::printf("[6.10.1] cylinder: %lld on edges, %lld onto their 3-D curves\n",
                static_cast<long long>(r.mNumOnEdges), static_cast<long long>(r.mNumOnCurves));
    EXPECT_GT(r.mNumOnEdges, 0);
    EXPECT_EQ(r.mNumOnCurves, r.mNumOnEdges);
    std::size_t with_curve = 0;
    for (mio::BrepEdge& r_e : cyl.mEdges)
        if (r_e.mHasCurve) {
            ++with_curve;
            double* pc = r_e.mCurve.mControlPoints.As<double>();
            for (std::size_t i = 0; i < r_e.mCurve.mControlPoints.Shape()[0]; ++i)
                pc[3 * i + 2] += 0.1;
        }
    ASSERT_GT(with_curve, 0u);
    const auto w = mio::project_to_cad(off, cyl);
    EXPECT_EQ(w.mNumOnCurves, 0);
    EXPECT_EQ(w.mNumOnEdges, r.mNumOnEdges);
    bool warned = false;
    for (const std::string& r_w : w.mWarnings)
        warned = warned || r_w.find("3-D curve") != std::string::npos;
    EXPECT_TRUE(warned);
    const double* p = w.mMesh.Points().As<double>();
    const double* q = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < w.mMesh.NumPoints(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            EXPECT_NEAR(p[3 * i + k], q[3 * i + k], 1e-9) << i;
}

// A point of the holed disc's cap moved into the hole lands on the plane
// there, outside the face's trims: it is moved, counted and named.
TEST(ProjectToCad, AProjectionOutsideTheTrimsIsReported) {
    const mio::BrepInfo disc = brep_solids::thick_disc_with_hole();
    const mio::Mesh m = surface(disc, 0.3);
    // The first point of the bottom cap that no edge runs through.
    const auto r0 = mio::project_to_cad(m, disc);
    std::int64_t pick = -1;
    {
        std::vector<char> on_edge(m.NumPoints(), 0);
        for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
            const auto cells = m.Cells(b);
            if (cells.Type() != "line")
                continue;
            for (std::size_t c = 0; c < 2 * cells.NumCells(); ++c)
                on_edge[static_cast<std::size_t>(mio::detail::read_int(cells.Conn(), c))] = 1;
        }
        const double* p = m.Points().As<double>();
        for (std::size_t i = 0; i < m.NumPoints() && pick < 0; ++i)
            if (!on_edge[i] && std::fabs(p[3 * i + 2]) < 1e-12)
                pick = static_cast<std::int64_t>(i);
    }
    ASSERT_GE(pick, 0);
    EXPECT_EQ(r0.mNumOutsideTrim, 0);
    const mio::Mesh off = moved(m, [&](std::size_t I, std::array<double, 3>& rX) {
        if (static_cast<std::int64_t>(I) == pick)
            rX = {0.98, 0.0, 0.1};
    });
    const auto r = mio::project_to_cad(off, disc);
    ASSERT_EQ(r.mNumOutsideTrim, 1);
    EXPECT_EQ(r.mOutsideTrim[0], pick);
    const double* p = r.mMesh.Points().As<double>();
    EXPECT_NEAR(p[3 * pick], 0.98, 1e-12);
    EXPECT_NEAR(p[3 * pick + 2], 0.0, 1e-12);
    EXPECT_FALSE(r.mWarnings.empty());
    // With a distance limit below the move it stays.
    mio::ProjectToCadOptions o;
    o.mMaxDistance = 0.05;
    const auto far = mio::project_to_cad(off, disc, o);
    EXPECT_EQ(far.mNumTooFar, 1);
    EXPECT_EQ(far.mNumOutsideTrim, 0);
}

TEST(ProjectToCad, RefusalsAreByName) {
    const mio::BrepInfo sphere = brep_solids::sphere(1.0);
    mio::Mesh bare;
    bare.AssignPoints(mio::NDArray::Uninit(mio::DType::Float64, {0, 3}));
    EXPECT_THROW((void)mio::project_to_cad(bare, sphere), std::invalid_argument);
    const mio::Mesh m = surface(sphere, 0.4);
    mio::ProjectToCadOptions o;
    o.mTrimSamples = 1;
    EXPECT_THROW((void)mio::project_to_cad(m, sphere, o), std::invalid_argument);
    o = {};
    o.mMaxDistance = -1.0;
    EXPECT_THROW((void)mio::project_to_cad(m, sphere, o), std::invalid_argument);
    // An index beyond the model's faces.
    EXPECT_THROW((void)mio::project_to_cad(m, brep_solids::square()), std::invalid_argument);
}
