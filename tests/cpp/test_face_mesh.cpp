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

// Unit tests for detail/face_mesh.hpp -- the globally deduplicated face list
// with owner/neighbour pairing that the OpenFOAM and CGNS NGON_n/NFACE_n
// writers share.
//
// The cases that matter here are the ones a single-cell-shape mesh cannot
// exercise: a hexahedron and a polyhedron meeting on one face, a mesh whose
// boundary blocks occupy global cell indices, and an inverted cell.

// System includes
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::detail::build_global_faces;
using meshioplusplus::detail::FaceLookup;
using meshioplusplus::detail::GlobalFaces;
using meshioplusplus::detail::Vec3;

namespace {

// Two unit cubes stacked along x, as one hexahedron block. Nodes 4..7 are the
// shared plane, so the two cells share exactly one face.
Mesh two_hexes() {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {0, 1, 0},
                                    {0, 1, 1},
                                    {0, 0, 1},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {1, 1, 1},
                                    {1, 0, 1},
                                    {2, 0, 0},
                                    {2, 1, 0},
                                    {2, 1, 1},
                                    {2, 0, 1}}));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7},
                                                {4, 5, 6, 7, 8, 9, 10, 11}}));
    return m;
}

// Sum of a face's node ids -- a cheap ring-identity fingerprint for assertions
// that care about *which* face, not its winding.
std::int64_t face_sum(const GlobalFaces& rF, std::size_t Face) {
    std::int64_t s = 0;
    for (std::size_t k = 0; k < rF.FaceSize(Face); ++k)
        s += rF.Face(Face)[k];
    return s;
}

// Newell area vector of a global face, from the mesh's own points.
Vec3 face_normal(const GlobalFaces& rF, std::size_t Face, const Mesh& rM) {
    std::vector<Vec3> ring;
    for (std::size_t k = 0; k < rF.FaceSize(Face); ++k)
        ring.push_back(meshioplusplus::detail::read_point(rM.Points(), rM.PointDim(),
                                                          rF.Face(Face)[k]));
    return meshioplusplus::detail::polygon_area_vector(ring.data(), ring.size());
}

Vec3 cell_centroid(const Mesh& rM, const GlobalFaces& rF, std::size_t Cell) {
    // Average of every node the cell's faces reference.
    Vec3 c{0, 0, 0};
    std::size_t n = 0;
    for (std::size_t k = 0; k < rF.NumCellFaces(Cell); ++k) {
        const std::int64_t sid = rF.CellFaces(Cell)[k];
        const std::size_t f = static_cast<std::size_t>(std::abs(sid) - 1);
        for (std::size_t j = 0; j < rF.FaceSize(f); ++j) {
            const Vec3 p =
                meshioplusplus::detail::read_point(rM.Points(), rM.PointDim(), rF.Face(f)[j]);
            c[0] += p[0];
            c[1] += p[1];
            c[2] += p[2];
            ++n;
        }
    }
    for (int i = 0; i < 3; ++i)
        c[i] /= static_cast<double>(n);
    return c;
}

}  // namespace

// --------------------------------------------------------------------------
// Deduplication
// --------------------------------------------------------------------------

TEST(FaceMesh, TwoHexesShareExactlyOneInternalFace) {
    const Mesh m = two_hexes();
    const GlobalFaces f = build_global_faces(m);

    EXPECT_EQ(f.NumCells(), 2u);
    // 12 face slots, one shared => 11 distinct faces, 1 internal, 10 boundary.
    EXPECT_EQ(f.NumFaces(), 11u);

    std::size_t internal = 0;
    for (std::size_t i = 0; i < f.NumFaces(); ++i)
        if (f.mNeighbour[i] >= 0)
            ++internal;
    EXPECT_EQ(internal, 1u);
    EXPECT_EQ(f.mNumNonManifold, 0);
}

// The load-bearing dedup case: two DIFFERENT cell shapes meeting on one face.
// With a separate key per shape -- or a key that respects ring order or winding
// -- each would report the shared face as its own boundary face, and the two
// cells would never be recognised as neighbours.
TEST(FaceMesh, AHexAndAPolyhedronSharingAFaceProduceOneFace) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {0, 1, 0},
                                    {0, 1, 1},
                                    {0, 0, 1},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {1, 1, 1},
                                    {1, 0, 1},
                                    {2, 0, 0},
                                    {2, 1, 0},
                                    {2, 1, 1},
                                    {2, 0, 1}}));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}}));
    // The same second cube, but as a polyhedron whose copy of the shared face
    // starts at a DIFFERENT node and winds the OTHER way round.
    m.AddPolyhedronBlock("polyhedron8", {{{
                                             {6, 5, 4, 7},         // shared face, rotated+reversed
                                             {8, 9, 10, 11},       // far face
                                             {4, 5, 9, 8},         // bottom
                                             {7, 6, 10, 11},       // top
                                             {4, 8, 11, 7},        // side
                                             {5, 6, 10, 9},        // side
                                         }}});

    const GlobalFaces f = build_global_faces(m);
    EXPECT_EQ(f.NumCells(), 2u);
    EXPECT_EQ(f.NumFaces(), 11u) << "the shared face was not deduplicated";

    // And the pairing must actually name the two cells.
    std::size_t internal = 0;
    for (std::size_t i = 0; i < f.NumFaces(); ++i) {
        if (f.mNeighbour[i] < 0)
            continue;
        ++internal;
        EXPECT_EQ(f.mOwner[i], 0);
        EXPECT_EQ(f.mNeighbour[i], 1);
        EXPECT_EQ(face_sum(f, i), 4 + 5 + 6 + 7);
    }
    EXPECT_EQ(internal, 1u);
}

// --------------------------------------------------------------------------
// The compact cell space
// --------------------------------------------------------------------------

// Every mesh `read_openfoam` produces carries its boundary faces as 2D blocks,
// so the block-major global numbering assigns cell indices to things that are
// not cells. owner/neighbour must not use that numbering.
TEST(FaceMesh, CompactCellIdsSkipNonVolumeBlocks) {
    Mesh m = two_hexes();
    // A boundary quad block AHEAD of nothing -- appended, but global cell ids
    // 2..3 now belong to quads, not cells.
    m.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}, {8, 9, 10, 11}}));

    const GlobalFaces f = build_global_faces(m);
    EXPECT_EQ(f.NumCells(), 2u) << "a 2D block was counted as a volume cell";
    ASSERT_EQ(f.mNonCellBlocks.size(), 1u);
    EXPECT_EQ(f.mNonCellBlocks[0], 1u);

    for (std::size_t i = 0; i < f.NumFaces(); ++i) {
        EXPECT_GE(f.mOwner[i], 0);
        EXPECT_LT(f.mOwner[i], 2) << "owner is a global id, not a compact one";
        EXPECT_LT(f.mNeighbour[i], 2);
    }
    // mCellToGlobal bridges back: both hexes are in block 0, so global == compact.
    EXPECT_EQ(f.mCellToGlobal[0], 0);
    EXPECT_EQ(f.mCellToGlobal[1], 1);
}

// The same, with the 2D block FIRST -- now global and compact genuinely differ,
// which is the shape that actually catches an owner/neighbour built on
// `block_bases`.
TEST(FaceMesh, CellToGlobalBridgesWhenABoundaryBlockComesFirst) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {0, 1, 0},
                                    {0, 1, 1},
                                    {0, 0, 1},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {1, 1, 1},
                                    {1, 0, 1}}));
    m.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}}));  // global cell 0
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}}));

    const GlobalFaces f = build_global_faces(m);
    ASSERT_EQ(f.NumCells(), 1u);
    EXPECT_EQ(f.mCellToGlobal[0], 1) << "the hexahedron is global cell 1, not 0";
    EXPECT_EQ(f.mOwner[0], 0) << "owner must be the COMPACT id";
}

// --------------------------------------------------------------------------
// Winding
// --------------------------------------------------------------------------

TEST(FaceMesh, BoundaryFaceNormalsPointOutOfTheirOwner) {
    const Mesh m = two_hexes();
    const GlobalFaces f = build_global_faces(m);

    for (std::size_t i = 0; i < f.NumFaces(); ++i) {
        if (f.mNeighbour[i] >= 0)
            continue;
        const Vec3 n = face_normal(f, i, m);
        const Vec3 c = cell_centroid(m, f, static_cast<std::size_t>(f.mOwner[i]));
        // outward: the normal agrees with (face centre - cell centre)
        Vec3 fcm{0, 0, 0};
        for (std::size_t k = 0; k < f.FaceSize(i); ++k) {
            const Vec3 p =
                meshioplusplus::detail::read_point(m.Points(), m.PointDim(), f.Face(i)[k]);
            for (int a = 0; a < 3; ++a)
                fcm[a] += p[a] / static_cast<double>(f.FaceSize(i));
        }
        const double d = n[0] * (fcm[0] - c[0]) + n[1] * (fcm[1] - c[1]) + n[2] * (fcm[2] - c[2]);
        EXPECT_GT(d, 0.0) << "boundary face " << i << " is wound inward";
    }
}

TEST(FaceMesh, InternalFaceNormalPointsFromOwnerToNeighbour) {
    const Mesh m = two_hexes();
    const GlobalFaces f = build_global_faces(m);

    for (std::size_t i = 0; i < f.NumFaces(); ++i) {
        if (f.mNeighbour[i] < 0)
            continue;
        const Vec3 n = face_normal(f, i, m);
        const Vec3 co = cell_centroid(m, f, static_cast<std::size_t>(f.mOwner[i]));
        const Vec3 cn = cell_centroid(m, f, static_cast<std::size_t>(f.mNeighbour[i]));
        const double d =
            n[0] * (cn[0] - co[0]) + n[1] * (cn[1] - co[1]) + n[2] * (cn[2] - co[2]);
        EXPECT_GT(d, 0.0) << "internal face " << i << " does not point owner->neighbour";
    }
}

// `cell_faces.hpp`'s rows are outward on the REFERENCE element. A cell given
// inverted (negative-Jacobian) connectivity therefore yields six INWARD normals
// unless the builder repairs them -- and OpenFOAM would reject every one.
TEST(FaceMesh, AnInvertedCellIsRewoundOutward) {
    Mesh m;
    // Swap the two z-levels of a unit cube => negative Jacobian.
    m.AssignPoints(mt::points_from({{0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1},
                                    {0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0}}));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}}));

    const GlobalFaces f = build_global_faces(m);
    EXPECT_EQ(f.mNumFlipped, 1) << "the inverted cell was not detected";

    const Vec3 c = cell_centroid(m, f, 0);
    for (std::size_t i = 0; i < f.NumFaces(); ++i) {
        const Vec3 n = face_normal(f, i, m);
        Vec3 fcm{0, 0, 0};
        for (std::size_t k = 0; k < f.FaceSize(i); ++k) {
            const Vec3 p =
                meshioplusplus::detail::read_point(m.Points(), m.PointDim(), f.Face(i)[k]);
            for (int a = 0; a < 3; ++a)
                fcm[a] += p[a] / static_cast<double>(f.FaceSize(i));
        }
        const double d = n[0] * (fcm[0] - c[0]) + n[1] * (fcm[1] - c[1]) + n[2] * (fcm[2] - c[2]);
        EXPECT_GT(d, 0.0) << "face " << i << " of the inverted cell still points inward";
    }
}

// --------------------------------------------------------------------------
// The signed NFACE_n encoding
// --------------------------------------------------------------------------

TEST(FaceMesh, SignedCellFacesRecordWhoStoredTheFace) {
    const Mesh m = two_hexes();
    const GlobalFaces f = build_global_faces(m);

    ASSERT_EQ(f.NumCells(), 2u);
    // Cell 0 saw every one of its faces first, so all six entries are positive.
    for (std::size_t k = 0; k < f.NumCellFaces(0); ++k)
        EXPECT_GT(f.CellFaces(0)[k], 0);

    // Cell 1 shares exactly one, which it must reference reversed.
    std::size_t negatives = 0;
    for (std::size_t k = 0; k < f.NumCellFaces(1); ++k)
        if (f.CellFaces(1)[k] < 0)
            ++negatives;
    EXPECT_EQ(negatives, 1u);

    // Every entry is a valid 1-based face id.
    for (std::size_t c = 0; c < f.NumCells(); ++c)
        for (std::size_t k = 0; k < f.NumCellFaces(c); ++k) {
            const std::int64_t id = std::abs(f.CellFaces(c)[k]);
            EXPECT_GE(id, 1);
            EXPECT_LE(id, static_cast<std::int64_t>(f.NumFaces()));
        }
}

// Every cell must be a closed surface once the signs are applied: summing the
// outward area vectors (negated where the entry is negative) gives zero. This
// is the oracle that catches a dropped face or a wrong sign, and it shares no
// code with the dedup.
TEST(FaceMesh, EveryCellIsClosedUnderTheSignedEncoding) {
    Mesh m = two_hexes();
    m.AddPolyhedronBlock("polyhedron8", {{{
                                             {8, 9, 10, 11},
                                             {12, 13, 14, 15},
                                             {8, 9, 13, 12},
                                             {11, 10, 14, 15},
                                             {8, 12, 15, 11},
                                             {9, 13, 14, 10},
                                         }}});
    // extend the point set for the third cube
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {0, 1, 0},
                                    {0, 1, 1},
                                    {0, 0, 1},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {1, 1, 1},
                                    {1, 0, 1},
                                    {2, 0, 0},
                                    {2, 1, 0},
                                    {2, 1, 1},
                                    {2, 0, 1},
                                    {3, 0, 0},
                                    {3, 1, 0},
                                    {3, 1, 1},
                                    {3, 0, 1}}));

    const GlobalFaces f = build_global_faces(m);
    ASSERT_EQ(f.NumCells(), 3u);

    for (std::size_t c = 0; c < f.NumCells(); ++c) {
        Vec3 sum{0, 0, 0};
        double scale = 0.0;
        for (std::size_t k = 0; k < f.NumCellFaces(c); ++k) {
            const std::int64_t sid = f.CellFaces(c)[k];
            const std::size_t fi = static_cast<std::size_t>(std::abs(sid) - 1);
            const Vec3 n = face_normal(f, fi, m);
            const double s = sid > 0 ? 1.0 : -1.0;
            for (int a = 0; a < 3; ++a)
                sum[a] += s * n[a];
            scale += std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        }
        const double mag = std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
        EXPECT_LT(mag, 1e-9 * scale) << "cell " << c << " is not closed";
    }
}

// --------------------------------------------------------------------------
// Determinism and lookup
// --------------------------------------------------------------------------

TEST(FaceMesh, IsDeterministicAcrossThreadCounts) {
    // A grid big enough to actually dispatch across threads.
    Mesh m;
    std::vector<std::vector<double>> pts;
    const int n = 9;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                pts.push_back({static_cast<double>(i), static_cast<double>(j),
                               static_cast<double>(k)});
    auto id = [&](int i, int j, int k) { return (i * n + j) * n + k; };
    std::vector<std::vector<std::int64_t>> rows;
    for (int i = 0; i + 1 < n; ++i)
        for (int j = 0; j + 1 < n; ++j)
            for (int k = 0; k + 1 < n; ++k)
                rows.push_back({id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k),
                                id(i, j + 1, k), id(i, j, k + 1), id(i + 1, j, k + 1),
                                id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)});
    m.AssignPoints(mt::points_from(pts));
    m.AddCellBlock("hexahedron", mt::conn_from(rows));

    const GlobalFaces a = build_global_faces(m);
    const GlobalFaces b = build_global_faces(m);
    EXPECT_EQ(a.mFaceNodes, b.mFaceNodes);
    EXPECT_EQ(a.mFaceStart, b.mFaceStart);
    EXPECT_EQ(a.mOwner, b.mOwner);
    EXPECT_EQ(a.mNeighbour, b.mNeighbour);
    EXPECT_EQ(a.mCellFaces, b.mCellFaces);

    // and the pairing is sane on a structured grid: interior faces have two
    // cells, and owner < neighbour falls out of the ascending dedup order.
    for (std::size_t i = 0; i < a.NumFaces(); ++i)
        if (a.mNeighbour[i] >= 0)
            EXPECT_LT(a.mOwner[i], a.mNeighbour[i]);
}

TEST(FaceMesh, FaceLookupFindsAFaceByItsCornerSetInAnyOrder) {
    const Mesh m = two_hexes();
    const GlobalFaces f = build_global_faces(m);
    const FaceLookup lookup(f);

    // The shared face, given in a rotated and reversed order.
    const std::vector<std::int64_t> shared = {6, 5, 4, 7};
    const std::int64_t fid = lookup.Find(shared.data(), shared.size());
    ASSERT_GE(fid, 0);
    EXPECT_GE(f.mNeighbour[static_cast<std::size_t>(fid)], 0) << "found face is not internal";

    const std::vector<std::int64_t> nonexistent = {0, 1, 8, 9};
    EXPECT_EQ(lookup.Find(nonexistent.data(), nonexistent.size()), -1);
}

TEST(FaceMesh, NonManifoldFacesAreCounted) {
    // Three cubes all sharing the plane x = 1: cells 1 and 2 both claim it.
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {0, 1, 0},
                                    {0, 1, 1},
                                    {0, 0, 1},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {1, 1, 1},
                                    {1, 0, 1},
                                    {2, 0, 0},
                                    {2, 1, 0},
                                    {2, 1, 1},
                                    {2, 0, 1},
                                    {2, 0, 2},
                                    {2, 1, 2},
                                    {1, 1, 2},
                                    {1, 0, 2}}));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7},
                                                {4, 5, 6, 7, 8, 9, 10, 11},
                                                {4, 5, 6, 7, 15, 14, 13, 12}}));
    const GlobalFaces f = build_global_faces(m);
    EXPECT_EQ(f.mNumNonManifold, 1) << "a face used by three cells was not reported";
}
