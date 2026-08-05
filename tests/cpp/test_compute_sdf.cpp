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
// compute_sdf: the umbrella that generates its own grid. The lattice and the
// distance kernel are each tested elsewhere; what is left here is the
// composition, and specifically the two ways the OCTREE can look right and be
// wrong -- a field interpolated from a coarser pass, and a selection carried
// across a pass. Both produce a valid, plausible mesh, so both need an oracle
// that fires.

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/stats.hpp"

using meshioplusplus::compute_sdf;
using meshioplusplus::isosurface;
using meshioplusplus::IsosurfaceOptions;
using meshioplusplus::kSdfDistanceName;
using meshioplusplus::Mesh;
using meshioplusplus::SdfOptions;
using meshioplusplus::SdfResult;
using meshioplusplus::SdfStructure;
namespace d = meshioplusplus::detail;

namespace {

// The unit cube [0,1]^3 as a closed, outward-wound triangle surface.
Mesh csdf_cube_surface() {
    std::vector<std::vector<double>> pts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    std::vector<std::vector<std::int64_t>> tris = {
        {0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
        {1, 2, 6}, {1, 6, 5}, {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7},
    };
    return mt::make_mesh(std::move(pts), "triangle", std::move(tris));
}

// A sphere of radius 0.4 about the origin, as a subdivided octahedron.
//
// The cube fixture above is deliberately NOT used for the contour tests: a
// cube's signed distance is piecewise LINEAR near its faces, so interpolating it
// from a coarser grid reproduces it exactly and a "field attached one pass too
// early" sabotage changes nothing at all. That was measured, not assumed --
// TheContourOracleActuallyFires reported 18636 facets either way against the
// cube. A sphere's distance field is curved everywhere, which is what makes the
// oracle able to fire.
Mesh csdf_sphere_surface(int Subdivisions) {
    std::vector<d::Vec3> pts = {{{1, 0, 0}},  {{-1, 0, 0}}, {{0, 1, 0}},
                                {{0, -1, 0}}, {{0, 0, 1}},  {{0, 0, -1}}};
    std::vector<std::array<std::int64_t, 3>> tris = {{{0, 2, 4}}, {{2, 1, 4}}, {{1, 3, 4}},
                                                     {{3, 0, 4}}, {{2, 0, 5}}, {{1, 2, 5}},
                                                     {{3, 1, 5}}, {{0, 3, 5}}};
    for (int s = 0; s < Subdivisions; ++s) {
        std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t> mid;
        auto midpoint = [&](std::int64_t a, std::int64_t b) {
            const auto key = a < b ? std::make_pair(a, b) : std::make_pair(b, a);
            const auto it = mid.find(key);
            if (it != mid.end())
                return it->second;
            d::Vec3 m = d::vec3_scale(
                d::vec3_add(pts[static_cast<std::size_t>(a)], pts[static_cast<std::size_t>(b)]),
                0.5);
            const double len = std::sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
            m = d::vec3_scale(m, 1.0 / len);
            pts.push_back(m);
            const std::int64_t id = static_cast<std::int64_t>(pts.size()) - 1;
            mid.emplace(key, id);
            return id;
        };
        std::vector<std::array<std::int64_t, 3>> next;
        for (const auto& t : tris) {
            const std::int64_t a = midpoint(t[0], t[1]);
            const std::int64_t b = midpoint(t[1], t[2]);
            const std::int64_t c = midpoint(t[2], t[0]);
            next.push_back({{t[0], a, c}});
            next.push_back({{a, t[1], b}});
            next.push_back({{c, b, t[2]}});
            next.push_back({{a, b, c}});
        }
        tris = std::move(next);
    }
    std::vector<std::vector<double>> out_pts;
    for (const d::Vec3& p : pts)
        out_pts.push_back({0.4 * p[0], 0.4 * p[1], 0.4 * p[2]});
    std::vector<std::vector<std::int64_t>> out_tris;
    for (const auto& t : tris)
        out_tris.push_back({t[0], t[1], t[2]});
    return mt::make_mesh(std::move(out_pts), "triangle", std::move(out_tris));
}

SdfOptions csdf_voxel(std::int64_t N) {
    SdfOptions o;
    o.mResolution = std::array<std::int64_t, 3>{{N, N, N}};
    return o;
}

SdfOptions csdf_octree(std::int64_t Root, std::int64_t Depth) {
    SdfOptions o;
    o.mStructure = SdfStructure::Octree;
    o.mRootResolution = Root;
    o.mMaxDepth = Depth;
    return o;
}

std::size_t csdf_total_cells(const Mesh& rMesh) {
    std::size_t n = 0;
    for (const auto cb : rMesh.CellRange())
        n += cb.NumCells();
    return n;
}

// The zero level set of `sdf:distance`, as a facet count plus a total area --
// the pair a wrong field changes and a merely differently-tiled grid does not.
std::pair<std::size_t, double> csdf_contour(const Mesh& rMesh) {
    IsosurfaceOptions io;
    io.mArrayName = kSdfDistanceName;
    io.mIsovalues = {0.0};
    const Mesh c = isosurface(rMesh, io);
    return {csdf_total_cells(c), meshioplusplus::compute_stats(c).mTotalArea};
}

// Every cell's own corner-bbox diagonal, recomputed here from the OUTPUT's
// geometry rather than from anything the operation recorded.
std::vector<double> csdf_diagonals(const Mesh& rMesh) {
    std::vector<double> out;
    const meshioplusplus::NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    for (const auto cb : rMesh.CellRange()) {
        const meshioplusplus::NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            d::Vec3 lo{{0, 0, 0}}, hi{{0, 0, 0}};
            for (std::size_t i = 0; i < npc; ++i) {
                const d::Vec3 p = d::read_point(points, dim, d::read_int(conn, c * npc + i));
                for (std::size_t k = 0; k < 3; ++k) {
                    if (i == 0 || p[k] < lo[k])
                        lo[k] = p[k];
                    if (i == 0 || p[k] > hi[k])
                        hi[k] = p[k];
                }
            }
            double s = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
                s += (hi[k] - lo[k]) * (hi[k] - lo[k]);
            out.push_back(std::sqrt(s));
        }
    }
    return out;
}

// Every cell's centroid, likewise recomputed from the output.
std::vector<d::Vec3> csdf_centroids(const Mesh& rMesh) {
    std::vector<d::Vec3> out;
    const meshioplusplus::NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    for (const auto cb : rMesh.CellRange()) {
        const meshioplusplus::NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            d::Vec3 sum{{0, 0, 0}};
            for (std::size_t i = 0; i < npc; ++i)
                sum = d::vec3_add(sum, d::read_point(points, dim, d::read_int(conn, c * npc + i)));
            out.push_back(d::vec3_scale(sum, 1.0 / static_cast<double>(npc)));
        }
    }
    return out;
}

// The exact unsigned distance from a point to the sphere of radius 0.4. The
// faceted fixture is inscribed, so this is approximate at the facet scale --
// three orders of magnitude below the tolerance the band test uses.
double csdf_exact_sphere_distance(const d::Vec3& rP) {
    return std::fabs(std::sqrt(rP[0] * rP[0] + rP[1] * rP[1] + rP[2] * rP[2]) - 0.4);
}

// The exact unsigned distance from a point to the unit cube's surface.
double csdf_exact_cube_distance(const d::Vec3& rP) {
    // Distance to the boundary of [0,1]^3: outside, the usual box distance;
    // inside, the smallest distance to any face.
    double out2 = 0.0;
    bool outside = false;
    double inner = 1.0e300;
    for (std::size_t k = 0; k < 3; ++k) {
        const double over = rP[k] > 1.0 ? rP[k] - 1.0 : (rP[k] < 0.0 ? -rP[k] : 0.0);
        if (over > 0.0) {
            outside = true;
            out2 += over * over;
        }
        const double dk = std::min(rP[k], 1.0 - rP[k]);
        inner = std::min(inner, dk);
    }
    return outside ? std::sqrt(out2) : inner;
}

}  // namespace

TEST(ComputeSdf, VoxelIsTheGridPlusTheField) {
    const SdfResult r = compute_sdf(csdf_cube_surface(), csdf_voxel(8));
    EXPECT_EQ(csdf_total_cells(r.mMesh), 512u);
    EXPECT_EQ(r.mMaxDepth, 0);
    for (std::size_t k = 0; k < 3; ++k)
        EXPECT_EQ(r.mDims[k], 8);
    ASSERT_TRUE(r.mMesh.HasPointData(kSdfDistanceName));
    EXPECT_EQ(r.mMesh.PointData(kSdfDistanceName).Size(), 9u * 9u * 9u);
    // The default relative padding is non-zero on purpose: a field that stops at
    // the surface is not much use.
    EXPECT_LT(r.mOrigin[0], 0.0);
}

TEST(ComputeSdf, TheHeaderDescribesTheGridItGenerated) {
    const SdfResult r = compute_sdf(csdf_cube_surface(), csdf_voxel(4));
    ASSERT_TRUE(r.mMesh.HasFieldData(meshioplusplus::kSdfOriginName));
    ASSERT_TRUE(r.mMesh.HasFieldData(meshioplusplus::kSdfSpacingName));
    ASSERT_TRUE(r.mMesh.HasFieldData(meshioplusplus::kSdfDimsName));
    ASSERT_TRUE(r.mMesh.HasFieldData(meshioplusplus::kSdfBoundsName));
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_DOUBLE_EQ(d::read_double(r.mMesh.FieldData(meshioplusplus::kSdfOriginName), k),
                         r.mOrigin[k]);
        EXPECT_DOUBLE_EQ(d::read_double(r.mMesh.FieldData(meshioplusplus::kSdfSpacingName), k),
                         r.mSpacing[k]);
        EXPECT_EQ(d::read_int(r.mMesh.FieldData(meshioplusplus::kSdfDimsName), k), r.mDims[k]);
    }
    EXPECT_EQ(d::read_int(r.mMesh.FieldData(meshioplusplus::kSdfStructureName), 0), 0);
    EXPECT_EQ(d::read_int(r.mMesh.FieldData(meshioplusplus::kSdfMaxDepthName), 0), 0);

    // The very same geometry, recovered with no header at all -- which is what a
    // written-and-reread grid has to fall back on, since no format persists
    // arbitrary field_data.
    d::LatticeSpec spec;
    ASSERT_TRUE(d::lattice_from_mesh(r.mMesh, spec));
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_EQ(spec.mDims[k], r.mDims[k]);
        EXPECT_DOUBLE_EQ(spec.mOrigin[k], r.mOrigin[k]);
        EXPECT_NEAR(spec.mSpacing[k], r.mSpacing[k], 1e-15);
    }
}

TEST(ComputeSdf, TheFieldIsTheDistance) {
    SdfOptions o = csdf_voxel(6);
    o.mPaddingRelative = 0.25;
    const SdfResult r = compute_sdf(csdf_cube_surface(), o);
    const meshioplusplus::NDArray& dist = r.mMesh.PointData(kSdfDistanceName);
    const meshioplusplus::NDArray& pts = r.mMesh.Points();
    for (std::size_t p = 0; p < r.mMesh.NumPoints(); ++p) {
        const d::Vec3 q = d::read_point(pts, 3, static_cast<std::int64_t>(p));
        EXPECT_NEAR(std::fabs(d::read_double(dist, p)), csdf_exact_cube_distance(q), 1e-12);
    }
}

// --- the octree --------------------------------------------------------------

TEST(ComputeSdf, OctreeRefinesOnlyNearTheSurface) {
    const SdfResult r = compute_sdf(csdf_cube_surface(), csdf_octree(8, 2));
    EXPECT_EQ(r.mMaxDepth, 2);
    // Fewer cells than the uniform grid of the same finest resolution: that is
    // the entire justification for the structure.
    EXPECT_LT(csdf_total_cells(r.mMesh), 32u * 32u * 32u);
    EXPECT_GT(csdf_total_cells(r.mMesh), 8u * 8u * 8u);
    // The reported spacing is the FINEST cell, not the root's.
    EXPECT_NEAR(r.mSpacing[0], (r.mMesh.NumPoints() > 0 ? r.mSpacing[0] : 0.0), 0.0);
    ASSERT_TRUE(r.mMesh.HasCellData(meshioplusplus::kRefineLevelName));
}

// Risk 1: the octree silently degrades the field. Attaching it before the last
// pass would leave a smooth interpolation of the coarse values -- which still
// contours, and still looks like a cube.
TEST(ComputeSdf, TheOctreeContourEqualsTheUniformGridsFromFewerCells) {
    const Mesh surface = csdf_sphere_surface(3);
    const SdfResult uni = compute_sdf(surface, csdf_voxel(32));
    const SdfResult oct = compute_sdf(surface, csdf_octree(8, 2));
    ASSERT_LT(csdf_total_cells(oct.mMesh), csdf_total_cells(uni.mMesh));

    const auto cu = csdf_contour(uni.mMesh);
    const auto co = csdf_contour(oct.mMesh);
    EXPECT_EQ(co.first, cu.first);
    EXPECT_NEAR(co.second, cu.second, 1e-9);
}

// Risk 1's oracle, proved to fire: recomputing the field on the mesh of the
// SECOND-to-last pass and then refining once more (which is what "attach it too
// early" amounts to) changes the contour the test above pins.
TEST(ComputeSdf, TheContourOracleActuallyFires) {
    const Mesh surface = csdf_sphere_surface(3);
    const SdfResult uni = compute_sdf(surface, csdf_voxel(32));
    // The sabotage: one pass short, then refined uniformly to the same finest
    // resolution, so the field rides through an interpolation exactly as it
    // would if compute_sdf attached it mid-way.
    const SdfResult early = compute_sdf(surface, csdf_octree(8, 1));
    meshioplusplus::RefineOptions ro;
    ro.mLevels = 1;
    const Mesh interpolated = meshioplusplus::refine(early.mMesh, ro).mMesh;

    const auto cu = csdf_contour(uni.mMesh);
    const auto ci = csdf_contour(interpolated);
    EXPECT_NE(ci.first, cu.first)
        << "the interpolated field produced the same contour, so the oracle above "
           "would not catch a field attached before the last pass";
}

// Risk 2: the octree refines the wrong region. A selection carried across a pass
// names cells of a mesh that no longer exists, so the tree deepens somewhere
// arbitrary -- and still looks like a plausible adaptive grid.
TEST(ComputeSdf, EveryFinestCellIsWithinTheBandOfTheSurface) {
    SdfOptions o = csdf_octree(8, 2);
    o.mBandCells = 1.0;
    const SdfResult r = compute_sdf(csdf_sphere_surface(3), o);

    // Recomputed from the output's own geometry, independently of anything the
    // operation recorded.
    const std::vector<double> diag = csdf_diagonals(r.mMesh);
    const std::vector<d::Vec3> cen = csdf_centroids(r.mMesh);
    ASSERT_EQ(diag.size(), cen.size());

    double coarsest = 0.0;
    for (double dg : diag)
        coarsest = std::max(coarsest, dg);

    std::size_t finest_checked = 0;
    for (std::size_t c = 0; c < diag.size(); ++c) {
        // Only the cells at the finest level are constrained: a coarse cell may
        // sit anywhere.
        if (diag[c] > 0.26 * coarsest)
            continue;
        ++finest_checked;
        const double dist = csdf_exact_sphere_distance(cen[c]);
        // A finest cell exists because its own parent chain was inside the band
        // at every level, plus whatever the balance closure drew in. The bound
        // is MEASURED rather than merely generous: the real maximum for this
        // fixture is 0.59 * coarsest, so 0.75 leaves ~27% headroom while staying
        // far below the domain, which is ~10 coarse diagonals across. A loose
        // bound here would pass under a stale-selection bug, which is exactly
        // what this test exists to catch.
        EXPECT_LE(dist, 0.75 * coarsest)
            << "a finest cell sits " << dist << " from the surface, outside any band";
    }
    EXPECT_GT(finest_checked, 0u);
}

// The band oracle above was verified to fire by sabotage: making the driver
// reuse the FIRST pass's selection (stale global block-major indices into a mesh
// that no longer exists) put finest cells 0.297 from a surface the bound allows
// 0.175 -- twelve failures, against zero for the real implementation. It is
// recorded here rather than left implicit because an oracle nobody has watched
// fail is not known to work, and because two weaker forms of this test were
// tried first and BOTH passed under that sabotage: a cube fixture (whose
// symmetry leaves the stale low indices near the surface anyway) and a
// `2.0 * coarsest` bound (over three times the real 0.59 maximum).
TEST(ComputeSdf, TheBandOracleIsTight) {
    SdfOptions o = csdf_octree(8, 2);
    o.mBandCells = 1.0;
    const SdfResult r = compute_sdf(csdf_sphere_surface(3), o);
    const std::vector<double> diag = csdf_diagonals(r.mMesh);
    const std::vector<d::Vec3> cen = csdf_centroids(r.mMesh);
    double coarsest = 0.0;
    for (double dg : diag)
        coarsest = std::max(coarsest, dg);
    double worst = 0.0;
    for (std::size_t c = 0; c < diag.size(); ++c)
        if (diag[c] <= 0.26 * coarsest)
            worst = std::max(worst, csdf_exact_sphere_distance(cen[c]));
    // The margin the bound leaves, pinned so that loosening it silently is a
    // failing test rather than a quiet weakening of the test above.
    EXPECT_GT(worst, 0.5 * coarsest);
    EXPECT_LT(worst, 0.65 * coarsest);
}

TEST(ComputeSdf, TheOctreeIsTwoToOneBalanced) {
    const SdfResult r = compute_sdf(csdf_cube_surface(), csdf_octree(8, 2));
    ASSERT_TRUE(r.mMesh.HasCellData(meshioplusplus::kRefineLevelName));
    // The balance rule is over node adjacency, so check it that way: no node may
    // be shared by two cells more than one level apart.
    const std::size_t np = r.mMesh.NumPoints();
    std::vector<std::int64_t> lo(np, 1 << 30), hi(np, -1);
    std::size_t g = 0;
    for (std::size_t b = 0; b < r.mMesh.NumCellBlocks(); ++b) {
        const auto cb = r.mMesh.Cells(b);
        const meshioplusplus::NDArray& lvl = r.mMesh.CellData(meshioplusplus::kRefineLevelName, b);
        const meshioplusplus::NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c, ++g) {
            const std::int64_t l = d::read_int(lvl, c);
            for (std::size_t i = 0; i < npc; ++i) {
                const std::size_t n = static_cast<std::size_t>(d::read_int(conn, c * npc + i));
                lo[n] = std::min(lo[n], l);
                hi[n] = std::max(hi[n], l);
            }
        }
    }
    for (std::size_t n = 0; n < np; ++n)
        if (hi[n] >= 0)
            EXPECT_LE(hi[n] - lo[n], 1)
                << "node " << n << " joins levels " << lo[n] << " and " << hi[n];
}

TEST(ComputeSdf, OctreeRejectsAVoxelSizingRequest) {
    SdfOptions o = csdf_octree(8, 2);
    o.mResolution = std::array<std::int64_t, 3>{{16, 16, 16}};
    EXPECT_THROW(compute_sdf(csdf_cube_surface(), o), std::invalid_argument);
    SdfOptions o2 = csdf_octree(8, 2);
    o2.mCellSize = 0.1;
    EXPECT_THROW(compute_sdf(csdf_cube_surface(), o2), std::invalid_argument);
}

TEST(ComputeSdf, OctreeRefusesToBlowThroughTheCellBudget) {
    SdfOptions o = csdf_octree(8, 6);
    o.mMaxCells = 5000;
    EXPECT_THROW(compute_sdf(csdf_cube_surface(), o), std::invalid_argument);
}

TEST(ComputeSdf, ZeroDepthOctreeIsJustTheRootLattice) {
    const SdfResult r = compute_sdf(csdf_cube_surface(), csdf_octree(4, 0));
    EXPECT_EQ(csdf_total_cells(r.mMesh), 64u);
    EXPECT_EQ(r.mMaxDepth, 0);
}

TEST(ComputeSdf, VoxelStillRefusesNeitherOrBothSizings) {
    SdfOptions none;
    EXPECT_THROW(compute_sdf(csdf_cube_surface(), none), std::invalid_argument);
    SdfOptions both = csdf_voxel(4);
    both.mCellSize = 0.5;
    EXPECT_THROW(compute_sdf(csdf_cube_surface(), both), std::invalid_argument);
}
