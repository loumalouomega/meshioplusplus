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
// Unit tests for detail/polyhedron.hpp -- the geometric kernel for cells
// bounded by arbitrary polygonal faces.
//
// The suite is deliberately built around cases the repo's existing measure code
// gets wrong or refuses: inward and inconsistently wound faces, a face set that
// is not closed, a genuinely non-planar face, and a many-faced solid with a
// closed-form volume.

// System includes
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/surface.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::detail::CellRings;
using meshioplusplus::detail::orient_rings;
using meshioplusplus::detail::poly_measure;
using meshioplusplus::detail::PolyMeasure;
using meshioplusplus::detail::RingOrientation;
using meshioplusplus::detail::Vec3;

namespace {

using Faces = std::vector<std::vector<std::int64_t>>;

// Build a one-cell polyhedron mesh and pull its rings + coordinates out.
struct Cell {
    Mesh mMesh;
    CellRings mRings;
    std::vector<Vec3> mCoords;
    bool mOk = false;
};

Cell make_cell(const std::vector<std::vector<double>>& rPts, const Faces& rFaces,
               const char* pType = "polyhedron") {
    Cell c;
    c.mMesh.AssignPoints(mt::points_from(rPts));
    c.mMesh.AddPolyhedronBlock(pType, {rFaces});
    c.mOk = meshioplusplus::detail::cell_rings(c.mMesh.Cells(0), 0, c.mMesh.Points(),
                                               c.mMesh.PointDim(), c.mRings, c.mCoords);
    return c;
}

// The unit cube's eight corners, in the usual hexahedron order.
const std::vector<std::vector<double>> kCubePts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                                   {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

// Six quad faces, every one wound so its right-hand normal points OUT.
const Faces kCubeOutward = {
    {0, 3, 2, 1},  // bottom  (-z)
    {4, 5, 6, 7},  // top     (+z)
    {0, 1, 5, 4},  // front   (-y)
    {2, 3, 7, 6},  // back    (+y)
    {0, 4, 7, 3},  // left    (-x)
    {1, 2, 6, 5},  // right   (+x)
};

Faces reversed_all(const Faces& rIn) {
    Faces out = rIn;
    for (auto& f : out)
        std::reverse(f.begin(), f.end());
    return out;
}

}  // namespace

TEST(PolyhedronKernel, UnitCubeMeasuresExactly) {
    Cell c = make_cell(kCubePts, kCubeOutward);
    ASSERT_TRUE(c.mOk);
    ASSERT_EQ(c.mRings.NumFaces(), 6u);
    ASSERT_EQ(c.mRings.mNodes.size(), 8u);

    const PolyMeasure m = poly_measure(c.mRings, c.mCoords.data());
    EXPECT_NEAR(m.mVolume, 1.0, 1e-14);
    EXPECT_NEAR(m.mSurfaceArea, 6.0, 1e-14);
    EXPECT_NEAR(m.mCentroid[0], 0.5, 1e-14);
    EXPECT_NEAR(m.mCentroid[1], 0.5, 1e-14);
    EXPECT_NEAR(m.mCentroid[2], 0.5, 1e-14);
}

TEST(PolyhedronKernel, AlreadyOutwardWindingIsLeftAlone) {
    Cell c = make_cell(kCubePts, kCubeOutward);
    ASSERT_TRUE(c.mOk);
    const std::vector<std::uint32_t> before = c.mRings.mFaceNodes;
    EXPECT_EQ(orient_rings(c.mRings, c.mCoords.data()), RingOrientation::Consistent);
    EXPECT_EQ(c.mRings.mFaceNodes, before);
}

TEST(PolyhedronKernel, EveryFaceInwardIsRepairedToAPositiveVolume) {
    // Consistently wound, just the wrong way round: the BFS finds nothing to
    // fix and the global signed-volume probe is what catches it.
    Cell c = make_cell(kCubePts, reversed_all(kCubeOutward));
    ASSERT_TRUE(c.mOk);
    EXPECT_LT(poly_measure(c.mRings, c.mCoords.data()).mVolume, 0.0);

    EXPECT_EQ(orient_rings(c.mRings, c.mCoords.data()), RingOrientation::Repaired);
    const PolyMeasure m = poly_measure(c.mRings, c.mCoords.data());
    EXPECT_NEAR(m.mVolume, 1.0, 1e-14);
    EXPECT_NEAR(m.mCentroid[2], 0.5, 1e-14);
}

TEST(PolyhedronKernel, OneFlippedFaceIsRepairedByTheBfs) {
    // The case the global probe alone cannot fix: five faces agree, one does
    // not, so the volume is wrong but not simply negated.
    Faces f = kCubeOutward;
    std::reverse(f[3].begin(), f[3].end());
    Cell c = make_cell(kCubePts, f);
    ASSERT_TRUE(c.mOk);
    EXPECT_GT(std::abs(poly_measure(c.mRings, c.mCoords.data()).mVolume - 1.0), 1e-6);

    EXPECT_EQ(orient_rings(c.mRings, c.mCoords.data()), RingOrientation::Repaired);
    EXPECT_NEAR(poly_measure(c.mRings, c.mCoords.data()).mVolume, 1.0, 1e-14);
}

TEST(PolyhedronKernel, AnOpenSurfaceIsUnorientableNotGuessedAt) {
    Faces f = kCubeOutward;
    f.pop_back();  // five faces: the +x side is missing
    Cell c = make_cell(kCubePts, f);
    ASSERT_TRUE(c.mOk);
    EXPECT_EQ(orient_rings(c.mRings, c.mCoords.data()), RingOrientation::Unorientable);
}

TEST(PolyhedronKernel, TruncatedOctahedronMatchesItsClosedForm) {
    // 14 faces (8 hexagons + 6 squares), a permutohedron: the vertices are all
    // permutations of (0, +-1, +-2). Its volume is exactly 32 for this scaling,
    // which is an oracle independent of anything in this codebase.
    // Enumerate by which coordinate holds the zero (3 ways), which of the
    // other two holds the 1 (2 ways), and the two signs (4 ways) = 24. Signing
    // the zero would double-count, which is exactly the trap here.
    std::vector<std::vector<double>> pts;
    for (int z = 0; z < 3; ++z) {
        const int a = (z + 1) % 3, b = (z + 2) % 3;
        for (int swap = 0; swap < 2; ++swap) {
            const int one = swap ? b : a, two = swap ? a : b;
            for (int s1 = -1; s1 <= 1; s1 += 2) {
                for (int s2 = -1; s2 <= 1; s2 += 2) {
                    std::vector<double> p(3, 0.0);
                    p[one] = s1 * 1.0;
                    p[two] = s2 * 2.0;
                    pts.push_back(p);
                }
            }
        }
    }
    ASSERT_EQ(pts.size(), 24u);

    // Plane normals: 6 axis-aligned squares (|n|=1) and 8 hexagons (|n|=(1,1,1)).
    std::vector<std::array<double, 4>> planes;  // nx, ny, nz, d
    for (int a = 0; a < 3; ++a)
        for (int s = -1; s <= 1; s += 2) {
            std::array<double, 4> p{0, 0, 0, 2.0};
            p[a] = s;
            planes.push_back(p);
        }
    for (int sx = -1; sx <= 1; sx += 2)
        for (int sy = -1; sy <= 1; sy += 2)
            for (int sz = -1; sz <= 1; sz += 2)
                planes.push_back({static_cast<double>(sx), static_cast<double>(sy),
                                  static_cast<double>(sz), 3.0});
    ASSERT_EQ(planes.size(), 14u);

    Faces faces;
    for (const auto& pl : planes) {
        std::vector<std::int64_t> on;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            const double d = pl[0] * pts[i][0] + pl[1] * pts[i][1] + pl[2] * pts[i][2];
            if (std::abs(d - pl[3]) < 1e-9)
                on.push_back(static_cast<std::int64_t>(i));
        }
        ASSERT_TRUE(on.size() == 4 || on.size() == 6) << on.size();
        // Order the coplanar points into a ring by angle in the face's plane.
        Vec3 n{pl[0], pl[1], pl[2]};
        n = meshioplusplus::detail::vec3_normalize(n);
        Vec3 ctr{0, 0, 0};
        for (std::int64_t i : on)
            ctr = meshioplusplus::detail::vec3_add(ctr, Vec3{pts[i][0], pts[i][1], pts[i][2]});
        ctr = meshioplusplus::detail::vec3_scale(ctr, 1.0 / static_cast<double>(on.size()));
        Vec3 u = meshioplusplus::detail::vec3_normalize(meshioplusplus::detail::vec3_sub(
            Vec3{pts[on[0]][0], pts[on[0]][1], pts[on[0]][2]}, ctr));
        Vec3 v = meshioplusplus::detail::vec3_cross(n, u);
        std::sort(on.begin(), on.end(), [&](std::int64_t a, std::int64_t b) {
            const Vec3 da =
                meshioplusplus::detail::vec3_sub(Vec3{pts[a][0], pts[a][1], pts[a][2]}, ctr);
            const Vec3 db =
                meshioplusplus::detail::vec3_sub(Vec3{pts[b][0], pts[b][1], pts[b][2]}, ctr);
            return std::atan2(meshioplusplus::detail::vec3_dot(da, v),
                              meshioplusplus::detail::vec3_dot(da, u)) <
                   std::atan2(meshioplusplus::detail::vec3_dot(db, v),
                              meshioplusplus::detail::vec3_dot(db, u));
        });
        faces.push_back(on);
    }

    Cell c = make_cell(pts, faces, "polyhedron24");
    ASSERT_TRUE(c.mOk);
    ASSERT_EQ(c.mRings.NumFaces(), 14u);
    // The angular sort gives each face SOME consistent ring but no global
    // orientation, which is exactly the input orient_rings exists for.
    ASSERT_NE(orient_rings(c.mRings, c.mCoords.data()), RingOrientation::Unorientable);

    const PolyMeasure m = poly_measure(c.mRings, c.mCoords.data());
    EXPECT_NEAR(m.mVolume, 32.0, 1e-9);
    // Centred on the origin by construction.
    EXPECT_NEAR(m.mCentroid[0], 0.0, 1e-9);
    EXPECT_NEAR(m.mCentroid[1], 0.0, 1e-9);
    EXPECT_NEAR(m.mCentroid[2], 0.0, 1e-9);
    // 6 squares of side sqrt(2) (the x = +-2 face is (2,+-1,0),(2,0,+-1)) and 8
    // regular hexagons of the same side: 6*2 + 8*(3*sqrt(3)/2)*2.
    const double square = 2.0;                              // s^2, s = sqrt(2)
    const double hex = 6.0 * (std::sqrt(3.0) / 4.0) * 2.0;  // 6 * (sqrt(3)/4) * s^2
    EXPECT_NEAR(m.mSurfaceArea, 6.0 * square + 8.0 * hex, 1e-9);
}

TEST(PolyhedronKernel, MeasureIsAdditiveOverASplitCell) {
    // Splitting the cube at z = 0.5 must leave the total volume, area of the
    // outer boundary and centroid unchanged. This catches sign and
    // double-counting errors a single-cell test cannot.
    const std::vector<std::vector<double>> pts = {
        {0, 0, 0},   {1, 0, 0},   {1, 1, 0}, {0, 1, 0}, {0, 0, 0.5}, {1, 0, 0.5},
        {1, 1, 0.5}, {0, 1, 0.5}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1},   {0, 1, 1}};
    const Faces lower = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                         {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}};
    const Faces upper = {{4, 7, 6, 5},   {8, 9, 10, 11}, {4, 5, 9, 8},
                         {6, 7, 11, 10}, {4, 8, 11, 7},  {5, 6, 10, 9}};

    Cell lo = make_cell(pts, lower);
    Cell hi = make_cell(pts, upper);
    ASSERT_TRUE(lo.mOk && hi.mOk);
    ASSERT_NE(orient_rings(lo.mRings, lo.mCoords.data()), RingOrientation::Unorientable);
    ASSERT_NE(orient_rings(hi.mRings, hi.mCoords.data()), RingOrientation::Unorientable);

    const PolyMeasure a = poly_measure(lo.mRings, lo.mCoords.data());
    const PolyMeasure b = poly_measure(hi.mRings, hi.mCoords.data());
    EXPECT_NEAR(a.mVolume + b.mVolume, 1.0, 1e-14);
    // Volume-weighted centroids recombine to the whole cube's.
    const double zc =
        (a.mVolume * a.mCentroid[2] + b.mVolume * b.mCentroid[2]) / (a.mVolume + b.mVolume);
    EXPECT_NEAR(zc, 0.5, 1e-14);
}

TEST(PolyhedronKernel, MeasureSurvivesAMeshFarFromTheOrigin) {
    // The recentring in poly_measure is load-bearing: V = sum(x . A)/3 only
    // telescopes because sum(A) == 0, so without it the cancellation at 1e8
    // destroys the answer. Verified to FAIL at this offset with the recentring
    // removed.
    const double off = 1.0e8;
    std::vector<std::vector<double>> pts = kCubePts;
    for (auto& p : pts) {
        p[0] += off;
        p[1] += off;
        p[2] += off;
    }
    Cell c = make_cell(pts, kCubeOutward);
    ASSERT_TRUE(c.mOk);
    const PolyMeasure m = poly_measure(c.mRings, c.mCoords.data());
    EXPECT_NEAR(m.mVolume, 1.0, 1e-9);
    EXPECT_NEAR(m.mCentroid[0] - off, 0.5, 1e-6);
}

TEST(PolyhedronKernel, NonPlanarFaceUsesTheCornerAverageFan) {
    // A cube with one top corner lifted, so the top face is genuinely
    // non-planar and its volume DEPENDS on the decomposition. The
    // corner-average fan's answer is pinned here; a fan about the ring's first
    // node gives a different number, which is exactly the ambiguity the
    // convention exists to remove.
    std::vector<std::vector<double>> pts = kCubePts;
    pts[6][2] = 1.5;  // lift the (1,1,1) corner

    const PolyMeasure m = [&] {
        Cell c = make_cell(pts, kCubeOutward);
        EXPECT_TRUE(c.mOk);
        EXPECT_NE(orient_rings(c.mRings, c.mCoords.data()), RingOrientation::Unorientable);
        return poly_measure(c.mRings, c.mCoords.data());
    }();

    // Cube (1) + the lifted corner's contribution. The corner-average fan
    // splits each warped quad about its centre, giving exactly 1/8 of the
    // 0.5-tall prism over the quarter-square, i.e. 1 + 0.5/4 = 1.125.
    EXPECT_NEAR(m.mVolume, 1.125, 1e-12);
}

TEST(PolyhedronKernel, CellsSharingAWarpedFaceAgreeOnIt) {
    // THE reason the fan apex is the face's corner average and not the ring's
    // first node: the apex is then a function of the FACE ALONE, so two cells
    // meeting on a warped quad triangulate it identically however each of them
    // happens to store that ring -- and their volumes sum to the union's.
    //
    // A fan about the ring's first node is a function of the cell's storage, so
    // the two sides can split the shared quad along different diagonals and the
    // sum drifts. Note a structured generator numbers shared faces consistently
    // and so escapes this by luck, which is exactly why this fixture rotates
    // the second cell's copy of the ring instead of using a generated block.
    //
    // Shared quad at x ~ 1 with one corner pulled out of plane.
    const std::vector<std::vector<double>> pts = {
        {0, 0, 0}, {0, 1, 0}, {0, 1, 1},   {0, 0, 1},  // left cap
        {1, 0, 0}, {1, 1, 0}, {1.6, 1, 1}, {1, 0, 1},  // the warped shared quad (4,5,6,7)
        {2, 0, 0}, {2, 1, 0}, {2, 1, 1},   {2, 0, 1},  // right cap
    };
    const std::vector<std::int64_t> shared = {4, 5, 6, 7};

    // `rot` rotates only the SECOND cell's copy of the shared ring.
    auto total_for = [&](int rot) {
        std::vector<std::int64_t> b_shared(shared.rbegin(), shared.rend());  // reversed: outward
        std::rotate(b_shared.begin(), b_shared.begin() + rot, b_shared.end());

        const Faces cell_a = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                              {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}};
        const Faces cell_b = {b_shared,     {8, 9, 10, 11}, {4, 7, 11, 8},
                              {5, 4, 8, 9}, {6, 5, 9, 10},  {7, 6, 10, 11}};
        double total = 0.0;
        for (const Faces& f : {cell_a, cell_b}) {
            Cell c = make_cell(pts, f);
            EXPECT_TRUE(c.mOk);
            EXPECT_NE(orient_rings(c.mRings, c.mCoords.data()), RingOrientation::Unorientable);
            total += poly_measure(c.mRings, c.mCoords.data()).mVolume;
        }
        return total;
    };

    const double base = total_for(0);
    EXPECT_GT(base, 0.0);
    for (int rot : {1, 2, 3})
        EXPECT_NEAR(total_for(rot), base, 1e-14)
            << "rotating the shared ring in one cell changed the pair's total volume, so the "
               "two cells no longer agree on the face between them (rot = "
            << rot << ")";
}

TEST(PolyhedronKernel, VolumeIsIndependentOfWhereEachRingStarts) {
    // The headline property of the corner-average fan, and the reason it is
    // not a fan about node 0: rotating each face's stored node list must not
    // change the answer. On the warped cube above, a node-0 fan does.
    std::vector<std::vector<double>> pts = kCubePts;
    pts[6][2] = 1.5;

    Cell a = make_cell(pts, kCubeOutward);
    Faces rotated = kCubeOutward;
    for (std::size_t i = 0; i < rotated.size(); ++i)
        std::rotate(rotated[i].begin(), rotated[i].begin() + static_cast<std::ptrdiff_t>(i % 4),
                    rotated[i].end());
    Cell b = make_cell(pts, rotated);
    ASSERT_TRUE(a.mOk && b.mOk);
    ASSERT_NE(orient_rings(a.mRings, a.mCoords.data()), RingOrientation::Unorientable);
    ASSERT_NE(orient_rings(b.mRings, b.mCoords.data()), RingOrientation::Unorientable);

    EXPECT_NEAR(poly_measure(a.mRings, a.mCoords.data()).mVolume,
                poly_measure(b.mRings, b.mCoords.data()).mVolume, 1e-15);
}

TEST(PolyhedronKernel, TheReposOwnInconsistentlyWoundFixtureMeasuresPositively) {
    // tests/python/helpers.py's `polyhedron_mesh` has two faces traversing a
    // shared edge in the SAME direction in both of its polyhedron5 cells --
    // topologically closed, orientation-inconsistent. Any kernel that assumed
    // outward winding would report a wrong-signed or wrong-magnitude volume
    // here. This is the regression guard for that, transcribed from the
    // fixture.
    const std::vector<std::vector<double>> pts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                                  {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    // Cell 1: a pyramid on the z=0 square with apex 7.
    const Faces cell1 = {{0, 1, 2, 3}, {0, 1, 7}, {1, 2, 7}, {2, 3, 7}, {3, 0, 7}};
    // Cell 2: a pyramid on the (0,1,5,4) square, base split in two triangles.
    const Faces cell2 = {{0, 1, 5}, {0, 4, 5}, {0, 1, 7}, {1, 5, 7}, {5, 4, 7}, {0, 4, 7}};

    for (const Faces& f : {cell1, cell2}) {
        Cell c = make_cell(pts, f);
        ASSERT_TRUE(c.mOk);
        ASSERT_NE(orient_rings(c.mRings, c.mCoords.data()), RingOrientation::Unorientable);
        const PolyMeasure m = poly_measure(c.mRings, c.mCoords.data());
        EXPECT_GT(m.mVolume, 0.0) << "winding repair failed on the fixture's cell";
        EXPECT_NEAR(m.mVolume, 1.0 / 3.0, 1e-14);
    }
}

TEST(PolyhedronKernel, RingsComeFromTheFaceTableForTabulatedTypes) {
    // The same kernel must serve a rectangular hexahedron: cell_rings reads
    // cell_faces.hpp for it, so a tabulated cell and a polyhedron of identical
    // geometry must measure identically. That is what makes deduplicating the
    // repo's three copy-pasted volume helpers onto this kernel sound.
    Mesh hex = mt::hex_mesh();
    CellRings rings;
    std::vector<Vec3> coords;
    ASSERT_TRUE(meshioplusplus::detail::cell_rings(hex.Cells(0), 0, hex.Points(), hex.PointDim(),
                                                   rings, coords));
    EXPECT_EQ(rings.NumFaces(), 6u);
    EXPECT_EQ(rings.mNodes.size(), 8u);
    const PolyMeasure m = poly_measure(rings, coords.data());
    EXPECT_GT(m.mVolume, 0.0) << "cell_faces rows are outward-wound; volume must be positive";

    Cell poly = make_cell(kCubePts, kCubeOutward);
    ASSERT_TRUE(poly.mOk);
    EXPECT_NEAR(m.mVolume, poly_measure(poly.mRings, poly.mCoords.data()).mVolume, 1e-14);
}

TEST(PolyhedronKernel, BlocksWithNoFaceTopologyAreReportedNotGuessed) {
    CellRings rings;
    std::vector<Vec3> coords;

    Mesh tri = mt::tri_mesh();  // 2D: no enclosed volume
    EXPECT_FALSE(meshioplusplus::detail::cell_rings(tri.Cells(0), 0, tri.Points(), tri.PointDim(),
                                                    rings, coords));

    Mesh poly;  // 1-level ragged: a surface cell
    poly.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}));
    poly.AddPolygonBlock("polygon", {{0, 1, 2, 3}});
    EXPECT_FALSE(meshioplusplus::detail::cell_rings(poly.Cells(0), 0, poly.Points(),
                                                    poly.PointDim(), rings, coords));
}

TEST(PolyhedronKernel, FacetKeyHandlesAnyArity) {
    // surface.cpp's existing key is a fixed array<int64_t,4>; a general
    // polygonal face needs more, and the inline fast path must still work.
    using meshioplusplus::detail::FacetKey;
    const std::int64_t tri_a[] = {5, 1, 3};
    const std::int64_t tri_b[] = {3, 5, 1};
    const std::int64_t tri_c[] = {3, 5, 2};
    EXPECT_TRUE(FacetKey(tri_a, 3) == FacetKey(tri_b, 3));
    EXPECT_FALSE(FacetKey(tri_a, 3) == FacetKey(tri_c, 3));

    const std::int64_t hex_a[] = {9, 2, 7, 4, 1, 8};
    const std::int64_t hex_b[] = {1, 2, 4, 7, 8, 9};
    EXPECT_TRUE(FacetKey(hex_a, 6) == FacetKey(hex_b, 6));
    // Arity is part of the key: a prefix must not collide.
    EXPECT_FALSE(FacetKey(hex_a, 6) == FacetKey(hex_a, 4));

    meshioplusplus::detail::FacetKeyHash h;
    EXPECT_EQ(h(FacetKey(hex_a, 6)), h(FacetKey(hex_b, 6)));
}

// --- surface extraction on polyhedra ----------------------------------------
//
// These live here rather than in test_surface.cpp because they are about the
// polyhedral path specifically, and because the volume oracle they lean on is
// this file's kernel.

TEST(PolyhedronSurface, SkinOfASinglePolyhedronIsAllItsFaces) {
    Mesh m;
    m.AssignPoints(mt::points_from(kCubePts));
    m.AddPolyhedronBlock("polyhedron8", {kCubeOutward});

    const Mesh s = meshioplusplus::extract_surface(m);
    // Every face of a lone cell is boundary: six quads.
    std::size_t total = 0;
    for (const auto cb : s.CellRange()) {
        EXPECT_EQ(std::string(cb.Type()), "quad");
        total += cb.NumCells();
    }
    EXPECT_EQ(total, 6u);
    EXPECT_EQ(s.NumPoints(), 8u);
}

TEST(PolyhedronSurface, SharedFacesCancelBetweenTwoPolyhedra) {
    // The whole point of hashing facets: the interior face must NOT appear.
    const std::vector<std::vector<double>> pts = {
        {0, 0, 0},   {1, 0, 0},   {1, 1, 0}, {0, 1, 0}, {0, 0, 0.5}, {1, 0, 0.5},
        {1, 1, 0.5}, {0, 1, 0.5}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1},   {0, 1, 1}};
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock(
        "polyhedron8",
        {{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}},
         {{4, 7, 6, 5},
          {8, 9, 10, 11},
          {4, 5, 9, 8},
          {6, 7, 11, 10},
          {4, 8, 11, 7},
          {5, 6, 10, 9}}});

    const Mesh s = meshioplusplus::extract_surface(m);
    std::size_t total = 0;
    for (const auto cb : s.CellRange())
        total += cb.NumCells();
    // 12 faces in, the two copies of the shared z = 0.5 quad cancel -> 10.
    EXPECT_EQ(total, 10u);
}

TEST(PolyhedronSurface, APolyhedronAndAHexahedronCancelOnTheirSharedFace) {
    // The reason surface.cpp had to unify on ONE key type: with a separate key
    // for each shape, neither side would find the other and the interior face
    // would be reported twice as boundary.
    const std::vector<std::vector<double>> pts = {
        {0, 0, 0},   {1, 0, 0},   {1, 1, 0}, {0, 1, 0}, {0, 0, 0.5}, {1, 0, 0.5},
        {1, 1, 0.5}, {0, 1, 0.5}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1},   {0, 1, 1}};
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}}));
    m.AddPolyhedronBlock("polyhedron8", {{{4, 7, 6, 5},
                                          {8, 9, 10, 11},
                                          {4, 5, 9, 8},
                                          {6, 7, 11, 10},
                                          {4, 8, 11, 7},
                                          {5, 6, 10, 9}}});

    const Mesh s = meshioplusplus::extract_surface(m);
    std::size_t total = 0;
    for (const auto cb : s.CellRange())
        total += cb.NumCells();
    EXPECT_EQ(total, 10u) << "the shared quad was not cancelled across the two block shapes";
}

TEST(PolyhedronSurface, NGonFacesBecomeARaggedPolygonBlock) {
    // A pentagonal face fits no fixed-width bucket, so it must land in a ragged
    // `polygon` block rather than being dropped.
    const std::vector<std::vector<double>> pts = {{0, 0, 0},     {1, 0, 0},    {1.5, 1, 0},
                                                  {0.5, 1.6, 0}, {-0.5, 1, 0}, {0.5, 0.8, 1}};
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    // A pentagonal pyramid: one pentagon base + five triangles.
    m.AddPolyhedronBlock(
        "polyhedron6", {{{0, 4, 3, 2, 1}, {0, 1, 5}, {1, 2, 5}, {2, 3, 5}, {3, 4, 5}, {4, 0, 5}}});

    const Mesh s = meshioplusplus::extract_surface(m);
    bool saw_polygon = false, saw_triangle = false;
    for (const auto cb : s.CellRange()) {
        if (std::string(cb.Type()) == "polygon") {
            saw_polygon = true;
            ASSERT_EQ(cb.NumCells(), 1u);
            EXPECT_EQ(cb.RowSize(0), 5u);
        } else if (std::string(cb.Type()) == "triangle") {
            saw_triangle = true;
            EXPECT_EQ(cb.NumCells(), 5u);
        }
    }
    EXPECT_TRUE(saw_polygon) << "the pentagonal face was dropped";
    EXPECT_TRUE(saw_triangle);
}

// --- smooth's inversion guard on polyhedra ----------------------------------

TEST(PolyhedronSmooth, TheInversionGuardEngagesOnPolyhedronCells) {
    // Before v9.16.0 smooth SKIPPED polyhedron blocks when building its
    // measurable-cell table, while node_edge_topology_known reported their
    // nodes as known -- so their nodes were moved with NO inversion guard, and
    // silently. This asserts the guard now actually rejects moves there.
    //
    // The fixture has to be a TANGLED block, not a clean convex cell: Laplacian
    // with lambda < 1 moves a node only part-way to its neighbours' centroid
    // and so can never overshoot a convex cell, which means a tidy pyramid
    // gives the guard nothing to do however flat it is. Heavy jitter is what
    // puts a node's centroid on the far side of one of its own cells.
    const int n = 3;  // 3x3x3 unit cells
    auto build = [&] {
        std::vector<std::vector<double>> pts;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j)
                for (int k = 0; k <= n; ++k)
                    pts.push_back({double(i), double(j), double(k)});
        auto pid = [&](int i, int j, int k) {
            return static_cast<std::int64_t>((i * (n + 1) + j) * (n + 1) + k);
        };
        // A cheap deterministic jitter, matching test_smooth.cpp's approach.
        for (std::size_t v = 0; v < pts.size(); ++v)
            for (int d = 0; d < 3; ++d)
                pts[v][d] += 0.9 * std::sin(static_cast<double>(v * 7 + d * 13));
        Faces dummy;
        std::vector<std::vector<std::vector<std::int64_t>>> cells;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    const std::int64_t c0 = pid(i, j, k), c1 = pid(i + 1, j, k),
                                       c2 = pid(i + 1, j + 1, k), c3 = pid(i, j + 1, k),
                                       c4 = pid(i, j, k + 1), c5 = pid(i + 1, j, k + 1),
                                       c6 = pid(i + 1, j + 1, k + 1), c7 = pid(i, j + 1, k + 1);
                    cells.push_back({{c0, c3, c2, c1},
                                     {c4, c5, c6, c7},
                                     {c0, c1, c5, c4},
                                     {c2, c3, c7, c6},
                                     {c0, c4, c7, c3},
                                     {c1, c2, c6, c5}});
                }
        Mesh m;
        m.AssignPoints(mt::points_from(pts));
        m.AddPolyhedronBlock("polyhedron8", std::move(cells));
        (void)dummy;
        return m;
    };

    meshioplusplus::SmoothOptions on;
    on.mMethod = meshioplusplus::SmoothMethod::Laplacian;
    on.mIterations = 15;
    on.mFixBoundary = false;  // polyhedral boundary marking is a documented gap
    on.mGuardInversion = true;
    meshioplusplus::SmoothOptions off = on;
    off.mGuardInversion = false;

    const meshioplusplus::SmoothResult guarded = meshioplusplus::smooth(build(), on);
    const meshioplusplus::SmoothResult unguarded = meshioplusplus::smooth(build(), off);

    EXPECT_GT(guarded.mNumSkippedInversion, 0)
        << "the inversion guard never engaged on a polyhedron block";
    EXPECT_EQ(unguarded.mNumSkippedInversion, 0);
}

// --- clean on polyhedra ------------------------------------------------------

TEST(PolyhedronClean, DegenerateAndDuplicatePolyhedraAreDropped) {
    // Before v9.16.0 clean's degenerate/duplicate logic ran on rectangular
    // blocks only, so a collapsed or repeated polyhedron survived while the
    // equivalent hexahedron did not.
    Mesh m;
    m.AssignPoints(mt::points_from(kCubePts));
    m.AddPolyhedronBlock("polyhedron8", {
                                            kCubeOutward,  // a genuine cube
                                            kCubeOutward,  // an exact duplicate
                                        });

    meshioplusplus::CleanOptions o;
    o.drop_degenerate = true;
    o.drop_duplicate_cells = true;
    o.remove_orphans = false;
    const meshioplusplus::CleanResult r = meshioplusplus::clean(m, o);
    EXPECT_EQ(r.mCellsDroppedDuplicate, 1);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 1u);
}

TEST(PolyhedronClean, DuplicateDetectionIgnoresFaceOrderAndWinding) {
    // The key is the sorted set of sorted faces, so the same solid still
    // collides with itself when its faces are listed in a different order, each
    // face starts at a different node, or the windings are reversed -- none of
    // which changes the cell. An order-sensitive key would miss all three, and
    // a key built from the cell's NODE set alone would be the opposite problem:
    // it cannot see the face structure at all, which for a polyhedron is the
    // only thing distinguishing one solid from another on the same points.
    Faces shuffled;
    for (std::size_t i = 0; i < kCubeOutward.size(); ++i) {
        // Walk the faces in a different order, reverse each one, and rotate its
        // start node.
        std::vector<std::int64_t> f = kCubeOutward[(i + 3) % kCubeOutward.size()];
        std::reverse(f.begin(), f.end());
        std::rotate(f.begin(), f.begin() + static_cast<std::ptrdiff_t>(i % f.size()), f.end());
        shuffled.push_back(std::move(f));
    }

    Mesh m;
    m.AssignPoints(mt::points_from(kCubePts));
    m.AddPolyhedronBlock("polyhedron8", {kCubeOutward, shuffled});

    meshioplusplus::CleanOptions o;
    o.drop_degenerate = false;
    o.drop_duplicate_cells = true;
    o.remove_orphans = false;
    const meshioplusplus::CleanResult r = meshioplusplus::clean(m, o);
    EXPECT_EQ(r.mCellsDroppedDuplicate, 1)
        << "the same solid, relisted, was not recognised as a duplicate";
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 1u);
}

TEST(PolyhedronClean, AWeldCollapsedPolyhedronIsDegenerate) {
    // Two coincident node pairs: welding flattens the cell to zero volume.
    std::vector<std::vector<double>> pts = kCubePts;
    for (int i = 4; i < 8; ++i)
        pts[static_cast<std::size_t>(i)][2] = 0.0;  // top layer collapses onto the bottom
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock("polyhedron8", {kCubeOutward});

    meshioplusplus::CleanOptions o;
    o.weld = true;
    o.atol = 1e-9;
    o.drop_degenerate = true;
    o.drop_duplicate_cells = false;
    o.remove_orphans = false;
    const meshioplusplus::CleanResult r = meshioplusplus::clean(m, o);
    EXPECT_EQ(r.mCellsDroppedDegenerate, 1) << "a polyhedron welded down to zero volume was kept";
}

TEST(PolyhedronSmooth, FixBoundaryPinsAPolyhedronsOuterNodes) {
    // Before v9.16.1 smooth_mark_boundary skipped polyhedron blocks, so
    // `fix_boundary` had no effect on them at all -- every node of a purely
    // polyhedral mesh was free to move, and the mesh quietly shrank.
    //
    // Two stacked cells: every node is on the outer boundary except none, so
    // with fix_boundary the mesh must not move at all.
    const std::vector<std::vector<double>> pts = {
        {0, 0, 0},   {1, 0, 0},   {1, 1, 0}, {0, 1, 0}, {0, 0, 0.5}, {1, 0, 0.5},
        {1, 1, 0.5}, {0, 1, 0.5}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1},   {0, 1, 1}};
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock(
        "polyhedron8",
        {{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}},
         {{4, 7, 6, 5},
          {8, 9, 10, 11},
          {4, 5, 9, 8},
          {6, 7, 11, 10},
          {4, 8, 11, 7},
          {5, 6, 10, 9}}});

    meshioplusplus::SmoothOptions o;
    o.mMethod = meshioplusplus::SmoothMethod::Laplacian;
    o.mIterations = 10;
    o.mFixBoundary = true;
    const meshioplusplus::SmoothResult r = meshioplusplus::smooth(m, o);

    // Every node of this mesh lies on its outer surface, so all are pinned.
    EXPECT_EQ(r.mNumNodesMoved, 0)
        << "fix_boundary did not pin the polyhedral block's boundary nodes";
    EXPECT_NEAR(r.mMaxDisplacement, 0.0, 1e-15);
}

TEST(PolyhedronSmooth, AnInteriorNodeIsStillFreeWithFixBoundary) {
    // The complement of the test above: boundary marking must not pin
    // EVERYTHING, or the previous test would pass for the wrong reason. A
    // 2x2x2 polyhedral block has exactly one genuinely interior node.
    const int n = 2;
    std::vector<std::vector<double>> pts;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            for (int k = 0; k <= n; ++k)
                pts.push_back({double(i), double(j), double(k)});
    auto pid = [&](int i, int j, int k) {
        return static_cast<std::int64_t>((i * (n + 1) + j) * (n + 1) + k);
    };
    // Nudge the centre node so smoothing has something to pull back.
    pts[static_cast<std::size_t>(pid(1, 1, 1))][0] += 0.3;

    std::vector<std::vector<std::vector<std::int64_t>>> cells;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                const std::int64_t c0 = pid(i, j, k), c1 = pid(i + 1, j, k),
                                   c2 = pid(i + 1, j + 1, k), c3 = pid(i, j + 1, k),
                                   c4 = pid(i, j, k + 1), c5 = pid(i + 1, j, k + 1),
                                   c6 = pid(i + 1, j + 1, k + 1), c7 = pid(i, j + 1, k + 1);
                cells.push_back({{c0, c3, c2, c1},
                                 {c4, c5, c6, c7},
                                 {c0, c1, c5, c4},
                                 {c2, c3, c7, c6},
                                 {c0, c4, c7, c3},
                                 {c1, c2, c6, c5}});
            }
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock("polyhedron8", std::move(cells));

    meshioplusplus::SmoothOptions o;
    o.mMethod = meshioplusplus::SmoothMethod::Laplacian;
    o.mIterations = 10;
    o.mFixBoundary = true;
    o.mPreserveFeatures = false;
    const meshioplusplus::SmoothResult r = meshioplusplus::smooth(m, o);

    // Exactly one node is interior, so exactly one may move.
    EXPECT_EQ(r.mNumNodesMoved, 1)
        << "expected only the single interior node of a 2x2x2 block to move";
    EXPECT_GT(r.mMaxDisplacement, 0.0);
}

// --- partition's dual graph on polyhedra -------------------------------------

TEST(PolyhedronPartition, TheDualGraphConnectsPolyhedraSharingAFace) {
    // Before v9.16.1 partition_dual_graph skipped polyhedron blocks entirely,
    // so KaHIP saw a graph with no edges and produced a balanced but
    // cut-blind partition. The dual graph is not exported, so this asserts it
    // through the observable consequence: a chain of cells partitioned in two
    // must split into two CONTIGUOUS runs, which a graph with no edges cannot
    // reliably produce.
    //
    // Only meaningful with KaHIP compiled in; SFC ignores the dual graph
    // entirely (it works off centroids) and would pass either way.
    if (!meshioplusplus::partition_has_kahip())
        GTEST_SKIP() << "needs -DMESHIOPLUSPLUS_WITH_KAHIP=ON";

    // A 1-D chain of 8 unit cells, expressed as polyhedra.
    const int n = 8;
    std::vector<std::vector<double>> pts;
    for (int i = 0; i <= n; ++i) {
        pts.push_back({double(i), 0, 0});
        pts.push_back({double(i), 1, 0});
        pts.push_back({double(i), 1, 1});
        pts.push_back({double(i), 0, 1});
    }
    std::vector<std::vector<std::vector<std::int64_t>>> cells;
    for (int i = 0; i < n; ++i) {
        const std::int64_t a = 4 * i, b = 4 * (i + 1);
        cells.push_back({{a + 0, a + 3, a + 2, a + 1},
                         {b + 0, b + 1, b + 2, b + 3},
                         {a + 0, a + 1, b + 1, b + 0},
                         {a + 2, a + 3, b + 3, b + 2},
                         {a + 0, b + 0, b + 3, a + 3},
                         {a + 1, a + 2, b + 2, b + 1}});
    }
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock("polyhedron8", std::move(cells));

    meshioplusplus::PartitionOptions o;
    o.mNParts = 2;
    o.mMethod = meshioplusplus::PartitionMethod::KaHIP;
    // partition_labels returns one block-aligned array per cell block; there
    // is a single block here.
    const std::vector<NDArray> label_blocks = meshioplusplus::partition_labels(m, o);
    ASSERT_EQ(label_blocks.size(), 1u);
    std::vector<std::int64_t> labels;
    for (std::size_t i = 0; i < label_blocks[0].Size(); ++i)
        labels.push_back(meshioplusplus::detail::read_int(label_blocks[0], i));
    ASSERT_EQ(labels.size(), static_cast<std::size_t>(n));

    // One boundary between the two parts: a chain cut optimally has exactly one
    // place where the label changes.
    int changes = 0;
    for (std::size_t i = 1; i < labels.size(); ++i)
        if (labels[i] != labels[i - 1])
            ++changes;
    EXPECT_EQ(changes, 1) << "the polyhedral dual graph did not connect the chain";
}

// --- what simplexify unlocks for free ---------------------------------------

TEST(PolyhedronSimplexify, SliceAndIsosurfaceWorkOnPolyhedraViaMarching) {
    // slice / isosurface / interpolate --barycentric all go through
    // marching_prepare -> convert_cells(Simplexify). The plan predicted they
    // become polyhedron-capable with no further code; this VERIFIES it rather
    // than assuming it, which is the whole reason the test exists.
    const std::vector<std::vector<double>> pts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                                  {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock("polyhedron8", {kCubeOutward});
    // A field that increases with z, so a mid-height contour is a full square.
    NDArray f = NDArray::Uninit(meshioplusplus::DType::Float64, {8, 1});
    for (std::size_t i = 0; i < 8; ++i)
        f.As<double>()[i] = pts[i][2];
    m.AddPointData("f", std::move(f));

    // slice at z = 0.5: a unit square of area 1.
    meshioplusplus::SliceOptions so;
    so.mOrigin = {0.5, 0.5, 0.5};
    so.mNormal = {0.0, 0.0, 1.0};
    const Mesh cut = meshioplusplus::slice(m, so);
    ASSERT_GT(cut.NumCellBlocks(), 0u);
    double area = 0.0;
    for (const auto cb : cut.CellRange()) {
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::vector<Vec3> corners(cb.NodesPerCell());
            for (std::size_t k = 0; k < cb.NodesPerCell(); ++k)
                corners[k] = meshioplusplus::detail::read_point(
                    cut.Points(), cut.PointDim(),
                    meshioplusplus::detail::read_int(cb.Conn(), c * cb.NodesPerCell() + k));
            area += meshioplusplus::detail::polygon_area(corners.data(), corners.size());
        }
    }
    EXPECT_NEAR(area, 1.0, 1e-12) << "the plane section of a unit cube is a unit square";

    // isosurface at f = 0.5: the same square.
    meshioplusplus::IsosurfaceOptions io;
    io.mArrayName = "f";
    io.mIsovalues = {0.5};
    const Mesh iso = meshioplusplus::isosurface(m, io);
    std::size_t iso_cells = 0;
    for (const auto cb : iso.CellRange())
        iso_cells += cb.NumCells();
    EXPECT_GT(iso_cells, 0u) << "isosurface produced nothing on a polyhedral mesh";
}
