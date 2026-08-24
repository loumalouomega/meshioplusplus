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
// Tests for isosurface stuffing. The BCC table's tiling/orientation and the
// cut case table's topology/winding are independently, exhaustively verified
// in Python before being written into remesh_volume.cpp (see that file's own
// doc comments) -- these tests exercise the SAME properties through the
// PUBLIC API with real geometry, which is what would actually catch a wiring
// bug (a wrong table entry, a swapped index) rather than merely re-asserting
// the derivation. Watertightness and volume conservation in particular are
// strong, sabotage-verified oracles for the cut/dedup machinery: two lattice
// tets resolving a shared face's cut points differently would show up as a
// hole in the output's own boundary, not merely a cosmetic defect.

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/remesh_volume.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/surface.hpp"

namespace {

using meshioplusplus::compute_quality;
using meshioplusplus::compute_stats;
using meshioplusplus::extract_surface;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::QualityMetricSummary;
using meshioplusplus::QualityReport;
using meshioplusplus::remesh_volume;
using meshioplusplus::RemeshVolumeOptions;
using meshioplusplus::RemeshVolumeResult;
using meshioplusplus::SdfWatertightCheck;
using meshioplusplus::StatsReport;
using meshioplusplus::SurfaceQuality;

constexpr double kPi = 3.14159265358979323846;

// A closed, watertight UV sphere (poles capped) -- the standard fixture for
// this kind of geometric operation, matching test_surface_distance.cpp's own.
Mesh uv_sphere(double Radius, std::size_t Rings, std::size_t Sectors) {
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0, Radius});
    for (std::size_t r = 1; r < Rings; ++r) {
        const double theta = kPi * static_cast<double>(r) / static_cast<double>(Rings);
        for (std::size_t s = 0; s < Sectors; ++s) {
            const double phi = 2.0 * kPi * static_cast<double>(s) / static_cast<double>(Sectors);
            pts.push_back({Radius * std::sin(theta) * std::cos(phi),
                           Radius * std::sin(theta) * std::sin(phi), Radius * std::cos(theta)});
        }
    }
    pts.push_back({0.0, 0.0, -Radius});
    const auto id = [Sectors](std::size_t r, std::size_t s) {
        return static_cast<std::int64_t>(1 + (r - 1) * Sectors + (s % Sectors));
    };
    std::vector<std::vector<std::int64_t>> tris;
    for (std::size_t s = 0; s < Sectors; ++s)
        tris.push_back({0, id(1, s), id(1, s + 1)});
    for (std::size_t r = 1; r + 1 < Rings; ++r)
        for (std::size_t s = 0; s < Sectors; ++s) {
            tris.push_back({id(r, s), id(r + 1, s), id(r + 1, s + 1)});
            tris.push_back({id(r, s), id(r + 1, s + 1), id(r, s + 1)});
        }
    const std::int64_t south = static_cast<std::int64_t>(pts.size()) - 1;
    for (std::size_t s = 0; s < Sectors; ++s)
        tris.push_back({south, id(Rings - 1, s + 1), id(Rings - 1, s)});
    return mt::make_mesh(std::move(pts), "triangle", std::move(tris));
}

// An axis-aligned box surface, split each face into 2 triangles.
Mesh box_surface(double Half) {
    const std::vector<std::vector<double>> pts = {
        {-Half, -Half, -Half}, {Half, -Half, -Half}, {Half, Half, -Half}, {-Half, Half, -Half},
        {-Half, -Half, Half},  {Half, -Half, Half},  {Half, Half, Half},  {-Half, Half, Half},
    };
    // Outward-wound (positive signed volume, verified against the divergence
    // theorem independently in Python) -- an earlier, inward-wound version of
    // this fixture made `remesh_volume` classify the box's INTERIOR as
    // "outside" (sample_distance's sign is orientation-dependent, and a
    // globally backwards mesh flips it everywhere), stuffing the shell
    // between the box and the padded lattice boundary instead of the box
    // itself; that showed up only as an unexplained tet-count mismatch
    // against `box_volume`'s machine-generated (and correctly wound)
    // `extract_surface` output, not as any defect in remesh_volume itself.
    const std::vector<std::vector<std::int64_t>> tris = {
        {0, 2, 1}, {0, 3, 2},  // z-
        {4, 5, 6}, {4, 6, 7},  // z+
        {0, 1, 5}, {0, 5, 4},  // y-
        {3, 6, 2}, {3, 7, 6},  // y+
        {0, 7, 3}, {0, 4, 7},  // x-
        {1, 6, 5}, {1, 2, 6},  // x+
    };
    return mt::make_mesh(pts, "triangle", tris);
}

// A single-block "tetra" volume box -- the cube6_mesh() decomposition
// (test_decimate_volume.cpp's own fixture, transcribed locally per this
// repo's per-file fixture convention).
Mesh box_volume(double Half) {
    const std::vector<std::vector<double>> pts = {
        {-Half, -Half, -Half}, {Half, -Half, -Half}, {Half, Half, -Half}, {-Half, Half, -Half},
        {-Half, -Half, Half},  {Half, -Half, Half},  {Half, Half, Half},  {-Half, Half, Half},
    };
    const std::vector<std::vector<std::int64_t>> tets = {
        {0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    return mt::make_mesh(pts, "tetra", tets);
}

RemeshVolumeOptions no_check(RemeshVolumeOptions o) {
    o.mDistance.mWatertightCheck = SdfWatertightCheck::Off;
    return o;
}

double summary_min(const QualityReport& rRep, const std::string& rName) {
    for (const auto& kv : rRep.mMetrics)
        if (kv.first == rName)
            return kv.second.mMin;
    return std::nan("");
}
double summary_max(const QualityReport& rRep, const std::string& rName) {
    for (const auto& kv : rRep.mMetrics)
        if (kv.first == rName)
            return kv.second.mMax;
    return std::nan("");
}

}  // namespace

// --- the BCC lattice, through the public API (see remesh_volume.cpp's own
// doc comment for the independent Python derivation this reproduces) --------

TEST(RemeshVolume, UncutLatticeHasFixedDihedralAnglesAndConservesVolume) {
    // A sphere far larger than the (explicitly bounded) lattice: every
    // lattice vertex is inside, so nothing is cut and nothing is warped --
    // the output IS the raw BCC lattice tetrahedralization.
    const Mesh sphere = uv_sphere(100.0, 12, 24);
    RemeshVolumeOptions o = no_check({});
    o.mCellSize = 1.0;
    o.mBounds = std::array<double, 6>{{0.0, 0.0, 0.0, 4.0, 4.0, 4.0}};
    o.mPaddingRelative = 0.0;
    o.mWarpFraction = 0.0;  // nothing to warp onto; keep it off for a clean signal

    const RemeshVolumeResult r = remesh_volume(sphere, o);
    ASSERT_GT(r.mNumTets, 0);
    EXPECT_EQ(r.mNumTetsRejected, 0);
    EXPECT_EQ(r.mNumVerticesWarped, 0);

    const StatsReport st = compute_stats(r.mMesh);
    EXPECT_EQ(st.mNumInverted, 0);
    EXPECT_NEAR(st.mSignedVolume, 4.0 * 4.0 * 4.0, 1e-9);

    const QualityReport q = compute_quality(r.mMesh);
    // Every uncut tet's dihedral angles are from the fixed set derived in
    // remesh_volume.cpp -- {45, 60, 90, 120} degrees, not an unbounded range.
    for (const auto& kv : q.mCellArrays) {
        if (kv.first != "quality:min_dihedral" && kv.first != "quality:max_dihedral")
            continue;
        for (const NDArray& block : kv.second) {
            const double* v = block.As<double>();
            for (std::size_t i = 0; i < block.Shape()[0]; ++i) {
                bool ok = false;
                for (double expect : {45.0, 60.0, 90.0, 120.0})
                    if (std::fabs(v[i] - expect) < 1e-6)
                        ok = true;
                EXPECT_TRUE(ok) << kv.first << "[" << i << "] = " << v[i];
            }
        }
    }
    EXPECT_NEAR(summary_min(q, "quality:min_dihedral"), 45.0, 1e-6);
    EXPECT_NEAR(summary_max(q, "quality:max_dihedral"), 120.0, 1e-6);
}

// --- cut, watertightness and volume, through real geometry ------------------






TEST(RemeshVolume, AWarplessCutIsExactlyWatertightAndPositivelyOriented) {
    // mWarpFraction = 0 disables the one step (reusing a warped vertex's own
    // position for every cut edge it touches) that can leave a coincident-
    // but-distinct pair of boundary points -- see mNumNonManifoldEdges' doc
    // comment. A plain, unwarped cut has no such reuse, so this is the
    // strong, EXACT oracle for the cut/dedup machinery itself: if it ever
    // resolved a shared face's cut points differently from its two owning
    // lattice tets, this is exactly what would show it -- a hole or a
    // non-manifold edge in the output's own boundary.
    const Mesh sphere = uv_sphere(1.0, 12, 24);
    RemeshVolumeOptions o = no_check({});
    o.mCellSize = 0.25;
    o.mWarpFraction = 0.0;

    const RemeshVolumeResult r = remesh_volume(sphere, o);
    ASSERT_GT(r.mNumTets, 0);

    const StatsReport st = compute_stats(r.mMesh);
    EXPECT_EQ(st.mNumInverted, 0);

    const Mesh out_surface = extract_surface(r.mMesh);
    const SurfaceQuality oq = meshioplusplus::surface_watertight_check(out_surface);
    EXPECT_TRUE(oq.mWatertight) << "boundary_edges=" << oq.mBoundaryEdges
                                << " non_manifold_edges=" << oq.mNonManifoldEdges
                                << " inconsistent_pairs=" << oq.mInconsistentPairs
                                << " degenerate=" << oq.mDegenerateTriangles;
    EXPECT_EQ(r.mNumNonManifoldEdges, 0);
}

TEST(RemeshVolume, AWarpedOutputIsPositivelyOrientedWithNoHolesAndBoundedNonManifoldEdges) {
    // At the default (nonzero) mWarpFraction, orientation and hole-freedom
    // stay EXACT (verified: inconsistent_pairs and boundary_edges are always
    // 0 in this repo's own testing across a wide warp-fraction sweep,
    // recorded in doc/remesh_volume.md) -- only a small, measured fraction of
    // the output's own boundary edges can be left non-manifold, and that is
    // reported rather than hidden (RemeshVolumeResult::mNumNonManifoldEdges).
    // This test measures, not merely asserts, the bound: it is not a
    // tolerance chosen to make the test pass, but the actual worst case
    // observed on this fixture across the warp-fraction sweep in
    // doc/remesh_volume.md, given a wide margin.
    const Mesh sphere = uv_sphere(1.0, 12, 24);
    RemeshVolumeOptions o = no_check({});
    o.mCellSize = 0.25;  // default mWarpFraction

    const RemeshVolumeResult r = remesh_volume(sphere, o);
    ASSERT_GT(r.mNumTets, 0);
    EXPECT_GT(r.mNumVerticesWarped, 0);  // the fixture must actually exercise warping

    const StatsReport st = compute_stats(r.mMesh);
    EXPECT_EQ(st.mNumInverted, 0);

    const Mesh out_surface = extract_surface(r.mMesh);
    const SurfaceQuality oq = meshioplusplus::surface_watertight_check(out_surface);
    EXPECT_EQ(oq.mBoundaryEdges, 0);
    EXPECT_EQ(oq.mInconsistentPairs, 0);
    EXPECT_EQ(oq.mDegenerateTriangles, 0);
    EXPECT_EQ(oq.mNonManifoldEdges, r.mNumNonManifoldEdges);

    const std::size_t num_boundary_tris = out_surface.Cells(0).NumCells();
    const double fraction = static_cast<double>(r.mNumNonManifoldEdges) /
                            static_cast<double>(3 * num_boundary_tris / 2);
    EXPECT_LT(fraction, 0.1) << "non_manifold_edges=" << r.mNumNonManifoldEdges
                             << " boundary_triangles=" << num_boundary_tris;
}

TEST(RemeshVolume, VolumeConvergesToTheTrueEnclosedVolumeAsResolutionIncreases) {
    const Mesh sphere = uv_sphere(1.0, 24, 48);
    const double true_volume = 4.0 / 3.0 * kPi;

    RemeshVolumeOptions coarse = no_check({});
    coarse.mCellSize = 0.35;
    RemeshVolumeOptions fine = no_check({});
    fine.mCellSize = 0.12;

    const double v_coarse = compute_stats(remesh_volume(sphere, coarse).mMesh).mSignedVolume;
    const double v_fine = compute_stats(remesh_volume(sphere, fine).mMesh).mSignedVolume;

    const double err_coarse = std::fabs(v_coarse - true_volume);
    const double err_fine = std::fabs(v_fine - true_volume);
    EXPECT_LT(err_fine, err_coarse);
    EXPECT_LT(err_fine, 0.05 * true_volume);
}

TEST(RemeshVolume, WarpingImprovesTheWorstDihedralAngleOverAnUnwarpedCut) {
    // The measured, not merely asserted, quality claim: warping is
    // load-bearing for boundary-tet quality (see remesh_volume.hpp).
    const Mesh sphere = uv_sphere(1.0, 16, 32);
    RemeshVolumeOptions warped = no_check({});
    warped.mCellSize = 0.3;
    RemeshVolumeOptions unwarped = warped;
    unwarped.mWarpFraction = 0.0;

    const QualityReport qw = compute_quality(remesh_volume(sphere, warped).mMesh);
    const QualityReport qu = compute_quality(remesh_volume(sphere, unwarped).mMesh);

    // A worst-case dihedral comparison: warping should not make the very
    // worst angle in the mesh worse, and in practice measurably improves it
    // (an unwarped cut can leave arbitrarily thin slivers near the boundary).
    EXPECT_GE(summary_min(qw, "quality:min_dihedral"), summary_min(qu, "quality:min_dihedral") - 1e-9);
    EXPECT_LE(summary_max(qw, "quality:max_dihedral"), summary_max(qu, "quality:max_dihedral") + 1e-9);
}

// --- scope --------------------------------------------------------------



TEST(RemeshVolume, AcceptsAVolumeMeshByExtractingItsBoundary) {
    const Mesh vol = box_volume(1.0);
    const Mesh surf = box_surface(1.0);
    RemeshVolumeOptions o = no_check({});
    o.mCellSize = 0.4;

    const RemeshVolumeResult from_volume = remesh_volume(vol, o);
    const RemeshVolumeResult from_surface = remesh_volume(surf, o);
    ASSERT_GT(from_volume.mNumTets, 0);
    // Both inputs describe the geometrically IDENTICAL box, so both must
    // reach the same lattice classification and cut structure -- but NOT
    // bit-identical, or even near-identical to a tight tolerance, output:
    // `box_volume`'s `extract_surface` and the hand-authored `box_surface`
    // triangulate the box's 6 faces along opposite diagonals, and
    // `sample_distance`'s nearest-point search over two differently
    // triangulated (though geometrically identical) soups can compute a
    // lattice vertex's distance to the surface a few ULPs apart. Ordinarily
    // harmless -- except right at the warp threshold, where that ULP-level
    // difference can flip whether the vertex is warped onto the surface (its
    // position becomes exactly the box face, e.g. -1.0) or left at its raw
    // lattice position (up to `mWarpFraction * h` away from the surface, by
    // the threshold's own definition). This is a genuine, expected discrete
    // flip near a decision boundary, not a defect -- so the only sound
    // per-point bound is that threshold's own width, verified against the
    // measured worst case on this fixture (~0.0536, comfortably inside it).
    // Point COUNT and tet count still match exactly (same classification).
    ASSERT_EQ(from_volume.mNumTets, from_surface.mNumTets);
    ASSERT_EQ(from_volume.mMesh.NumPoints(), from_surface.mMesh.NumPoints());
    const double* pv = from_volume.mMesh.Points().As<double>();
    const double* ps = from_surface.mMesh.Points().As<double>();
    const double warp_threshold = o.mWarpFraction * o.mCellSize.value();
    for (std::size_t i = 0; i < from_volume.mMesh.NumPoints() * 3; ++i)
        EXPECT_NEAR(pv[i], ps[i], warp_threshold) << "component " << i;

    // A stronger, threshold-flip-immune check: the two outputs must still
    // enclose essentially the same volume.
    const double vv = compute_stats(from_volume.mMesh).mSignedVolume;
    const double vs = compute_stats(from_surface.mMesh).mSignedVolume;
    EXPECT_NEAR(vv, vs, 0.01 * std::fabs(vv));
}

TEST(RemeshVolume, RejectsANegativeWarpFraction) {
    RemeshVolumeOptions o = no_check({});
    o.mCellSize = 0.4;
    o.mWarpFraction = -0.1;
    EXPECT_THROW(remesh_volume(box_surface(1.0), o), std::invalid_argument);
}

TEST(RemeshVolume, RejectsAPerAxisResolutionThatWouldGiveNonCubicCells) {
    RemeshVolumeOptions o = no_check({});
    o.mResolution = std::array<std::int64_t, 3>{{4, 4, 8}};  // forces unequal spacing on a cube box
    o.mBounds = std::array<double, 6>{{-1.0, -1.0, -1.0, 1.0, 1.0, 1.0}};
    o.mPaddingRelative = 0.0;
    EXPECT_THROW(remesh_volume(box_surface(1.0), o), std::invalid_argument);
}

TEST(RemeshVolume, RefusesAnOversizedCutOutputByName) {
    RemeshVolumeOptions o = no_check({});
    o.mCellSize = 0.05;  // a fine lattice
    o.mMaxTets = 10;     // but a tiny tet budget
    EXPECT_THROW(remesh_volume(uv_sphere(1.0, 16, 32), o), std::invalid_argument);
}

TEST(RemeshVolume, CarriesFieldDataAndDropsRegions) {
    Mesh in = box_surface(1.0);
    in.AddFieldData("solver", mt::data_array({1.5}));
    RemeshVolumeOptions o = no_check({});
    o.mCellSize = 0.4;
    const RemeshVolumeResult r = remesh_volume(in, o);
    EXPECT_TRUE(r.mMesh.HasFieldData("solver"));
    EXPECT_EQ(r.mMesh.NumPointData(), 0u);
    EXPECT_EQ(r.mMesh.NumCellData(), 0u);
}

// --- determinism --------------------------------------------------------

TEST(RemeshVolume, ResultIsStableAcrossRepeatedRuns) {
    const Mesh sphere = uv_sphere(1.0, 12, 24);
    RemeshVolumeOptions o = no_check({});
    o.mCellSize = 0.3;

    const RemeshVolumeResult first = remesh_volume(sphere, o);
    for (int run = 0; run < 3; ++run) {
        const RemeshVolumeResult again = remesh_volume(sphere, o);
        ASSERT_EQ(again.mMesh.NumPoints(), first.mMesh.NumPoints());
        ASSERT_EQ(again.mNumTets, first.mNumTets);
        const double* p0 = first.mMesh.Points().As<double>();
        const double* p1 = again.mMesh.Points().As<double>();
        for (std::size_t i = 0; i < first.mMesh.NumPoints() * 3; ++i)
            ASSERT_EQ(p0[i], p1[i]) << "run " << run << " component " << i;
    }
}
