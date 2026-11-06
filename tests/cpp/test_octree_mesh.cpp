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
// octree_mesh: the castellated output. Its inside region converges to the
// ball's volume and equals the volume the surrogate boundary encloses (the
// divergence theorem: a closed surrogate surface); every surrogate face has
// no inside cell across it; hanging nodes, the signed distance and the
// closest points equal brute forces; the output is the same on every thread
// count; a holed or inverted skin still meshes; and bad options are refused
// by name.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <set>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

#include "thread_sweep.hpp"

// Project includes
#include "analytic_surfaces.hpp"
#include "mesh_fixtures.hpp"
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "../../src/cpp/src/detail/bcc_octree_tets.hpp"
#include "../../src/cpp/src/detail/feature_snap.hpp"
#include "../../src/cpp/src/detail/hex_fit.hpp"
#include "../../src/cpp/src/detail/tri_bvh.hpp"
#include "../../src/cpp/src/operations/octree_mesh_internal.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::detail::Vec3;

constexpr double kPi = 3.14159265358979323846;
// Hexahedron facets (detail/cell_faces), outward.
constexpr int kFacet[6][4] = {{0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4},
                              {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};

struct OmBox {
    Vec3 mLo, mHi;
};

std::vector<OmBox> om_boxes(const mio::Mesh& rM) {
    const double* p = rM.Points().As<double>();
    const auto block = rM.Cells(0);
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    std::vector<OmBox> out(block.NumCells());
    for (std::size_t n = 0; n < out.size(); ++n) {
        out[n].mLo = {1e300, 1e300, 1e300};
        out[n].mHi = {-1e300, -1e300, -1e300};
        for (std::size_t v = 0; v < 8; ++v)
            for (std::size_t k = 0; k < 3; ++k) {
                const double x = p[static_cast<std::size_t>(c[n * 8 + v]) * 3 + k];
                out[n].mLo[k] = std::min(out[n].mLo[k], x);
                out[n].mHi[k] = std::max(out[n].mHi[k], x);
            }
    }
    return out;
}

double om_volume(const OmBox& rB) {
    return (rB.mHi[0] - rB.mLo[0]) * (rB.mHi[1] - rB.mLo[1]) * (rB.mHi[2] - rB.mLo[2]);
}

std::vector<std::int64_t> om_region(const mio::Mesh& rM, const char* pName, mio::RegionKind Kind) {
    const std::size_t r = rM.FindRegion(pName, Kind);
    EXPECT_NE(r, mio::Mesh::npos) << pName;
    const auto& e = rM.Region(r).mEntries;
    return std::vector<std::int64_t>(e.As<std::int64_t>(), e.As<std::int64_t>() + e.Size());
}

mio::OctreeMeshOptions om_castellated(double MaxSize, double MinSize) {
    mio::OctreeMeshOptions o;
    o.mOutput = mio::OctreeMeshOutput::Castellated;
    o.mMaxSize = MaxSize;
    o.mMinSize = MinSize;
    o.mPaddingRelative = 0.05;
    return o;
}

}  // namespace

TEST(OctreeMesh, InsideConvergesToTheBallAndTheSurrogateBoundaryEnclosesIt) {
    const mio::Mesh sphere = mio::icosa_sphere(1.0, 4);
    double prev = 1e300;
    for (double hmin : {0.2, 0.1, 0.05}) {
        auto o = om_castellated(0.4, hmin);
        o.mSkinGrading = mio::OctreeMeshSkinGrading::Uniform;
        o.mSurfaceSize = hmin;
        const auto r = mio::octree_mesh(sphere, o);
        EXPECT_EQ(r.mSurfaceLevel, r.mMaxLevel);
        const mio::Mesh& m = r.mMesh;
        const auto boxes = om_boxes(m);
        const auto inside = om_region(m, mio::kOctreeMeshInsideRegion, mio::RegionKind::Cell);
        EXPECT_EQ(static_cast<std::int64_t>(inside.size()), r.mNumInside);
        EXPECT_EQ(inside.size(), boxes.size());  // fill = inside
        double vol = 0.0;
        for (std::int64_t n : inside)
            vol += om_volume(boxes[static_cast<std::size_t>(n)]);
        // The ball discretised by the icosphere: a little below 4/3 pi.
        const double err = std::abs(vol - 4.0 / 3.0 * kPi);
        EXPECT_LT(err, prev);
        prev = err;
        // Divergence theorem over the surrogate faces: (1/3) sum x . n A.
        const auto sides = om_region(m, mio::kOctreeMeshSurrogateRegion, mio::RegionKind::Side);
        EXPECT_EQ(static_cast<std::int64_t>(sides.size() / 2), r.mNumSurrogateFaces);
        const double* p = m.Points().As<double>();
        const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
        double flux = 0.0;
        for (std::size_t s = 0; s < sides.size(); s += 2) {
            const auto cell = static_cast<std::size_t>(sides[s]);
            const auto f = static_cast<std::size_t>(sides[s + 1]);
            Vec3 q[4];
            for (std::size_t k = 0; k < 4; ++k)
                for (std::size_t d = 0; d < 3; ++d)
                    q[k][d] = p[static_cast<std::size_t>(c[cell * 8 + kFacet[f][k]]) * 3 + d];
            const Vec3 n = mio::detail::vec3_cross(mio::detail::vec3_sub(q[2], q[0]),
                                                   mio::detail::vec3_sub(q[3], q[1]));
            Vec3 mid{0, 0, 0};
            for (const Vec3& x : q)
                mid = mio::detail::vec3_add(mid, mio::detail::vec3_scale(x, 0.25));
            flux += mio::detail::vec3_dot(mid, n) * 0.5 / 3.0;
        }
        EXPECT_NEAR(flux, vol, 1e-9 * vol);
        EXPECT_EQ(r.mNumInverted, 0);
        EXPECT_NEAR(r.mMinScaledJacobian, 1.0, 1e-12);
    }
    EXPECT_LT(prev, 0.06 * 4.0 / 3.0 * kPi);
}

TEST(OctreeMesh, SurrogateFacesHaveNoInsideCellAcrossThem) {
    const auto r = mio::octree_mesh(analytic::torus(1.0, 0.4, 48, 24), [] {
        auto o = om_castellated(0.4, 0.05);
        o.mFill = mio::OctreeMeshFill::All;
        return o;
    }());
    const mio::Mesh& m = r.mMesh;
    const auto boxes = om_boxes(m);
    const auto inside = om_region(m, mio::kOctreeMeshInsideRegion, mio::RegionKind::Cell);
    std::vector<std::uint8_t> is_in(boxes.size(), 0);
    for (std::int64_t n : inside)
        is_in[static_cast<std::size_t>(n)] = 1;
    const auto sides = om_region(m, mio::kOctreeMeshSurrogateRegion, mio::RegionKind::Side);
    ASSERT_GT(sides.size(), 0u);
    // Direction of each facet: -x, +x, -y, +y, -z, +z.
    const int axis[6] = {0, 0, 1, 1, 2, 2};
    const double sign[6] = {-1, 1, -1, 1, -1, 1};
    for (std::size_t s = 0; s < sides.size(); s += 2) {
        const auto cell = static_cast<std::size_t>(sides[s]);
        const auto f = static_cast<std::size_t>(sides[s + 1]);
        EXPECT_TRUE(is_in[cell]);
        const OmBox& b = boxes[cell];
        Vec3 probe;
        for (std::size_t k = 0; k < 3; ++k)
            probe[k] = 0.5 * (b.mLo[k] + b.mHi[k]);
        const auto a = static_cast<std::size_t>(axis[f]);
        const double h = b.mHi[a] - b.mLo[a];
        probe[a] = (sign[f] > 0 ? b.mHi[a] : b.mLo[a]) + sign[f] * 1e-6 * h;
        // Every cell touching the face's centre from outside is not inside.
        for (std::size_t n = 0; n < boxes.size(); ++n) {
            bool contains = true;
            for (std::size_t k = 0; k < 3; ++k)
                contains = contains && probe[k] > boxes[n].mLo[k] && probe[k] < boxes[n].mHi[k];
            if (contains)
                EXPECT_FALSE(is_in[n]) << "cell " << cell << " facet " << f;
        }
    }
}

TEST(OctreeMesh, HangingNodesAndDistancesEqualBruteForces) {
    const mio::Mesh skin = analytic::torus(1.0, 0.4, 32, 16);
    const auto r = mio::octree_mesh(skin, om_castellated(0.5, 0.08));
    const mio::Mesh& m = r.mMesh;
    const auto boxes = om_boxes(m);
    const double* p = m.Points().As<double>();
    const std::int64_t* hang = m.PointData(mio::kOctreeHangingName).As<std::int64_t>();
    std::int64_t count = 0;
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        bool hanging = false;
        for (const OmBox& b : boxes) {
            bool on = true, corner = true;
            for (std::size_t k = 0; k < 3; ++k) {
                const double x = p[i * 3 + k];
                on = on && x >= b.mLo[k] && x <= b.mHi[k];
                corner = corner && (x == b.mLo[k] || x == b.mHi[k]);
            }
            hanging = hanging || (on && !corner);
        }
        EXPECT_EQ(hang[i], hanging ? 1 : 0) << "point " << i;
        count += hanging;
    }
    EXPECT_EQ(count, r.mNumHanging);
    EXPECT_GT(r.mNumHanging, 0);
    // The distances and closest points are the distance kernel's, signed by
    // the exact winding number.
    const auto soup = mio::detail::build_triangle_soup(skin, "");
    mio::SurfaceDistanceOptions dopt;
    dopt.mWatertightCheck = mio::SdfWatertightCheck::Off;
    const auto query = mio::detail::build_distance_query(soup, dopt);
    const auto tree = mio::detail::build_tri_bvh(soup);
    std::vector<Vec3> pts(m.NumPoints());
    for (std::size_t i = 0; i < pts.size(); ++i)
        pts[i] = {p[i * 3], p[i * 3 + 1], p[i * 3 + 2]};
    const auto hits = mio::detail::query_closest_points(query, pts);
    const double* d = m.PointData(mio::kOctreeMeshDistanceName).As<double>();
    const double* cp = m.PointData(mio::kOctreeMeshClosestPointName).As<double>();
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const double w = mio::detail::fast_winding(tree, pts[i], 0.0);
        EXPECT_DOUBLE_EQ(std::abs(d[i]), hits[i].mDistance);
        if (std::abs(w - 0.5) > 0.2)
            EXPECT_EQ(d[i] < 0.0, w > 0.5) << "point " << i;
        for (std::size_t k = 0; k < 3; ++k)
            EXPECT_EQ(cp[i * 3 + k], hits[i].mPoint[k]);
    }
}

TEST(OctreeMesh, RegionsAndFillsAgree) {
    const mio::Mesh skin = mio::icosa_sphere(1.0, 3);
    auto o = om_castellated(0.5, 0.1);
    o.mFill = mio::OctreeMeshFill::All;
    const auto all = mio::octree_mesh(skin, o);
    EXPECT_EQ(all.mNumCells, all.mNumLeaves);
    o.mFill = mio::OctreeMeshFill::Cover;
    const auto cover = mio::octree_mesh(skin, o);
    o.mFill = mio::OctreeMeshFill::Inside;
    const auto inside = mio::octree_mesh(skin, o);
    EXPECT_EQ(inside.mNumCells, inside.mNumInside);
    EXPECT_GT(cover.mNumCells, inside.mNumCells);
    EXPECT_LT(cover.mNumCells, all.mNumCells);
    EXPECT_EQ(cover.mNumInside, inside.mNumInside);
    EXPECT_GT(inside.mNumBoundaryCells, 0);
    // The skin's cut leaves are all at the finest level the sizes allow there.
    const auto lv = all.mMesh.CellData(mio::kOctreeLevelName, 0).As<std::int64_t>();
    const auto ct = all.mMesh.CellData(mio::kOctreeCutName, 0).As<std::int64_t>();
    for (std::size_t n = 0; n < static_cast<std::size_t>(all.mNumCells); ++n)
        if (ct[n])
            EXPECT_GE(lv[n], all.mMaxLevel - 1);
    EXPECT_EQ(all.mLeavesPerLevel.size(), static_cast<std::size_t>(all.mMaxLevel) + 1);
}

TEST(OctreeMesh, SizesFollowTheirCriteria) {
    // A sphere of radius 2: the curvature term asks for 2 sin(7.5 deg) * 2 ~ 0.52.
    const mio::Mesh sphere = mio::icosa_sphere(2.0, 4);
    auto o = om_castellated(1.0, 1.0 / 64.0);
    o.mFill = mio::OctreeMeshFill::Cover;
    o.mGapCells = 0.0;
    const auto r = mio::octree_mesh(sphere, o);
    const auto boxes = om_boxes(r.mMesh);
    const auto ct = r.mMesh.CellData(mio::kOctreeCutName, 0).As<std::int64_t>();
    const double want = 2.0 * std::sin(7.5 * kPi / 180.0) * 2.0;
    for (std::size_t n = 0; n < boxes.size(); ++n)
        if (ct[n]) {
            const double edge = boxes[n].mHi[0] - boxes[n].mLo[0];
            EXPECT_LE(edge, want * 1.05);
            EXPECT_GT(edge, want / 2.2);
        }
    // A thin plate: gap_cells leaves across its thickness.
    const mio::Mesh plate =
        mio::extract_surface(mio::box({{16, 16, 1}}, {{2.0, 2.0, 0.125}}, {{-1.0, -1.0, -0.0625}}));
    auto g = om_castellated(0.5, 1.0 / 128.0);
    g.mFill = mio::OctreeMeshFill::Cover;
    g.mCurvatureAngle = 0.0;
    g.mFeatureAngle = 0.0;
    g.mGapCells = 4.0;
    const auto pr = mio::octree_mesh(plate, g);
    const auto pb = om_boxes(pr.mMesh);
    const auto pc = pr.mMesh.CellData(mio::kOctreeCutName, 0).As<std::int64_t>();
    for (std::size_t n = 0; n < pb.size(); ++n)
        if (pc[n])
            EXPECT_LE(pb[n].mHi[0] - pb[n].mLo[0], 0.125 / 4.0 * 1.0001);
    EXPECT_GT(pr.mNumInside, 0);
}

TEST(OctreeMesh, DirtySkinsStillMesh) {
    const mio::Mesh sphere = mio::icosa_sphere(1.0, 3);
    auto o = om_castellated(0.5, 0.1);
    const auto clean = mio::octree_mesh(sphere, o);
    // Every 9th triangle removed: holes the winding number closes.
    std::vector<std::vector<double>> pts;
    const double* p = sphere.Points().As<double>();
    for (std::size_t i = 0; i < sphere.NumPoints(); ++i)
        pts.push_back({p[3 * i], p[3 * i + 1], p[3 * i + 2]});
    std::vector<std::vector<std::int64_t>> holed, flipped;
    const std::int64_t* c = sphere.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < sphere.Cells(0).NumCells(); ++t) {
        if (t % 9 != 0)
            holed.push_back({c[3 * t], c[3 * t + 1], c[3 * t + 2]});
        flipped.push_back({c[3 * t], c[3 * t + 2], c[3 * t + 1]});
    }
    const auto h = mio::octree_mesh(mt::make_mesh(pts, "triangle", holed), o);
    EXPECT_NEAR(static_cast<double>(h.mNumInside), static_cast<double>(clean.mNumInside),
                0.05 * static_cast<double>(clean.mNumInside));
    // Inside out: the same cells, and a warning.
    const auto f = mio::octree_mesh(mt::make_mesh(pts, "triangle", flipped), o);
    EXPECT_EQ(f.mNumInside, clean.mNumInside);
    ASSERT_FALSE(f.mWarnings.empty());
    EXPECT_NE(f.mWarnings[0].find("inside out"), std::string::npos);
}

TEST(OctreeMesh, IsTheSameOnEveryThreadCount) {
    const mio::Mesh skin = analytic::torus(1.0, 0.4, 32, 16);
    std::uint64_t ref = 0;
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    const int saved = mio_test::max_threads();
    for (int threads : {1, 2, 4}) {
        mio_test::set_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const auto r = mio::octree_mesh(skin, om_castellated(0.5, 0.05));
        mio::bench::MeshDigest d;
        d.Of(r.mMesh);
        if (ref == 0)
            ref = d.Value();
        EXPECT_EQ(d.Value(), ref);
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(saved);
#endif
}

namespace {

// Switch the flooded inside test off (the internal test hook) for a scope.
struct OmNoFlood {
    OmNoFlood() { om_setenv("0"); }
    ~OmNoFlood() { om_setenv(nullptr); }
    static void om_setenv(const char* pValue) {
#ifdef _WIN32
        _putenv_s("MESHIOPLUSPLUS_OCTREE_MESH_FLOOD", pValue != nullptr ? pValue : "");
#else
        if (pValue != nullptr)
            setenv("MESHIOPLUSPLUS_OCTREE_MESH_FLOOD", pValue, 1);
        else
            unsetenv("MESHIOPLUSPLUS_OCTREE_MESH_FLOOD");
#endif
    }
};

std::uint64_t om_digest(const mio::Mesh& rMesh) {
    mio::bench::MeshDigest d;
    d.Of(rMesh);
    return d.Value();
}

}  // namespace

// Roadmap 3.1.4.3: deep inside a closed skin the winding number is flooded in
// rather than evaluated, and every output is bit for bit what it was.
TEST(OctreeMesh, TheFloodedInsideTestChangesNoOutput) {
    const mio::Mesh sphere = mio::icosa_sphere(1.0, 3);
    for (auto output : {mio::OctreeMeshOutput::Tetra, mio::OctreeMeshOutput::Castellated,
                        mio::OctreeMeshOutput::Hexa}) {
        mio::OctreeMeshOptions o;
        o.mOutput = output;
        o.mMaxSize = 0.12;
        o.mMinSize = 0.06;
        const auto flooded = mio::octree_mesh(sphere, o);
        EXPECT_GT(flooded.mNumWindingSkipped, 0) << static_cast<int>(output);
        std::uint64_t plain_digest = 0;
        std::int64_t plain_rechecks = 0;
        {
            OmNoFlood off;
            const auto plain = mio::octree_mesh(sphere, o);
            EXPECT_EQ(plain.mNumWindingSkipped, 0);
            plain_digest = om_digest(plain.mMesh);
            plain_rechecks = plain.mNumWindingRechecks;
        }
        EXPECT_EQ(om_digest(flooded.mMesh), plain_digest) << static_cast<int>(output);
        EXPECT_EQ(flooded.mNumWindingRechecks, plain_rechecks) << static_cast<int>(output);
    }
}

// A cylinder whose caps lie exactly on lattice planes: the cut flags' closed
// test may miss a cap in the face between two leaves, and a flood through it
// would fill the outside. The flood's own test grows the boxes by a margin.
TEST(OctreeMesh, AFloodDoesNotCrossASkinOnALatticePlane) {
    const mio::Mesh skin = mio::convert_cells(mio::extract_surface(mio::cylinder(0.5, 2.0, 48, 16)),
                                              {mio::ConvertCellsMode::Simplexify, false})
                               .mMesh;
    for (auto output : {mio::OctreeMeshOutput::Tetra, mio::OctreeMeshOutput::Castellated}) {
        mio::OctreeMeshOptions o;
        o.mOutput = output;
        o.mMinSize = 0.1;
        o.mMaxSize = 0.5;
        const auto flooded = mio::octree_mesh(skin, o);
        EXPECT_GT(flooded.mNumWindingSkipped, 0);
        OmNoFlood off;
        const auto plain = mio::octree_mesh(skin, o);
        EXPECT_EQ(om_digest(flooded.mMesh), om_digest(plain.mMesh)) << static_cast<int>(output);
        EXPECT_EQ(flooded.mNumCells, plain.mNumCells);
    }
}

// Roadmap 3.1.4.3: a cylinder's rims on fine graded skins. The guard counted
// the slivers outside the skin, which the output drops, and held rim vertices
// off their chains; every chain is recovered now, and nothing kept inverts.
TEST(OctreeMesh, ACylindersRimsAreRecoveredOnFineSkins) {
    const mio::Mesh skin = mio::convert_cells(mio::extract_surface(mio::cylinder(0.5, 2.0, 48, 16)),
                                              {mio::ConvertCellsMode::Simplexify, false})
                               .mMesh;
    for (double hmin : {0.1, 0.05}) {
        mio::OctreeMeshOptions o;
        o.mMinSize = hmin;
        o.mMaxSize = 0.5;
        const auto r = mio::octree_mesh(skin, o);
        EXPECT_EQ(r.mNumFeatureChains, 4) << hmin;
        EXPECT_EQ(r.mNumFeatureChainsMissed, 0) << hmin;
        EXPECT_EQ(r.mNumCornersSnapped, r.mNumCorners) << hmin;
        EXPECT_EQ(r.mNumInverted, 0) << hmin;
        EXPECT_LT(r.mHausdorffToSkin, 0.25 * hmin) << hmin;
        EXPECT_TRUE(mio::check_conformity(r.mMesh).mConforming) << hmin;
    }
}

TEST(OctreeMesh, AnOpenSkinIsNotFlooded) {
    // The sphere less one triangle: its winding number is not an integer.
    const mio::Mesh closed = mio::icosa_sphere(1.0, 2);
    const auto block = closed.Cells(0);
    std::vector<std::vector<std::int64_t>> tris;
    for (std::size_t t = 1; t < block.NumCells(); ++t)
        tris.push_back({block.Conn().As<std::int64_t>()[t * 3],
                        block.Conn().As<std::int64_t>()[t * 3 + 1],
                        block.Conn().As<std::int64_t>()[t * 3 + 2]});
    std::vector<std::vector<double>> pts;
    const double* p = closed.Points().As<double>();
    for (std::size_t i = 0; i < closed.NumPoints(); ++i)
        pts.push_back({p[i * 3], p[i * 3 + 1], p[i * 3 + 2]});
    mio::OctreeMeshOptions o = om_castellated(0.2, 0.1);
    const auto r = mio::octree_mesh(mt::make_mesh(pts, "triangle", tris), o);
    EXPECT_EQ(r.mNumWindingSkipped, 0);
}

TEST(OctreeMesh, RefusesBadOptionsByName) {
    const mio::Mesh sphere = mio::icosa_sphere(1.0, 1);
    auto expect_throw = [&](mio::OctreeMeshOptions o, const char* pWhat) {
        try {
            mio::octree_mesh(sphere, o);
            ADD_FAILURE() << "no throw for " << pWhat;
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find(pWhat), std::string::npos) << e.what();
        }
    };
    auto o = om_castellated(0.5, 0.1);
    auto bad = o;
    bad.mMinSize = 1.0;
    expect_throw(bad, "min_size");
    bad = o;
    bad.mWindingThreshold = 1.5;
    expect_throw(bad, "winding_threshold");
    bad = o;
    bad.mMaxCells = 10;
    expect_throw(bad, "max_cells");
    bad = o;
    bad.mMinSize = 1e-7;
    expect_throw(bad, "levels");
    bad = o;
    bad.mSizeArray = "nope";
    expect_throw(bad, "size_array");
    bad = o;
    bad.mCrossingMargin = 0.5;
    expect_throw(bad, "crossing_margin");
    bad = o;
    bad.mZeroVolumeFactor = -1.0;
    expect_throw(bad, "zero_volume_factor");
    EXPECT_THROW(mio::octree_mesh(mt::line_mesh(), o), std::invalid_argument);
    EXPECT_THROW(mio::octree_mesh_output_from_name("quad"), std::invalid_argument);
    EXPECT_EQ(mio::octree_mesh_output_from_name("tet"), mio::OctreeMeshOutput::Tetra);
}

namespace {

// The checks every tetrahedral output must pass: positive volumes, every
// interior triangle shared by two tetrahedra with opposite orientations, and
// a boundary that is a closed 2-manifold (every boundary edge on exactly two
// boundary triangles).
struct OmTetCheck {
    double mVolume = 0.0;
    double mMinDihedral = 180.0, mMaxDihedral = 0.0;
    std::size_t mBoundaryFaces = 0;
};

OmTetCheck om_check_tets(const mio::Mesh& rM) {
    OmTetCheck out;
    const double* p = rM.Points().As<double>();
    const auto block = rM.Cells(0);
    EXPECT_EQ(block.Type(), "tetra");
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    auto at = [&](std::int64_t v) { return Vec3{p[3 * v], p[3 * v + 1], p[3 * v + 2]}; };
    std::map<std::array<std::int64_t, 3>, int> faces;
    constexpr int kFace[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    for (std::size_t t = 0; t < block.NumCells(); ++t) {
        const std::int64_t* v = c + 4 * t;
        const double vol = mio::detail::triple_product(mio::detail::vec3_sub(at(v[1]), at(v[0])),
                                                       mio::detail::vec3_sub(at(v[2]), at(v[0])),
                                                       mio::detail::vec3_sub(at(v[3]), at(v[0]))) /
                           6.0;
        EXPECT_GT(vol, 0.0) << "tet " << t;
        out.mVolume += vol;
        const auto a = mio::detail::tet_dihedral_angles(at(v[0]), at(v[1]), at(v[2]), at(v[3]));
        out.mMinDihedral = std::min(out.mMinDihedral, a[0]);
        out.mMaxDihedral = std::max(out.mMaxDihedral, a[5]);
        for (const auto& f : kFace) {
            std::array<std::int64_t, 3> k{v[f[0]], v[f[1]], v[f[2]]};
            int parity = 0;
            for (int x = 0; x < 3; ++x)
                for (int y = x + 1; y < 3; ++y)
                    parity ^= k[x] > k[y];
            std::sort(k.begin(), k.end());
            const int s = parity ? -1 : 1;
            auto [it, fresh] = faces.emplace(k, s);
            if (!fresh) {
                EXPECT_EQ(it->second, -s) << "a triangle seen twice the same way";
                it->second = 2;  // matched
            }
        }
    }
    std::map<std::array<std::int64_t, 2>, int> edges;
    for (const auto& [k, s] : faces) {
        if (s == 2)
            continue;
        ++out.mBoundaryFaces;
        for (int e = 0; e < 3; ++e)
            ++edges[{std::min(k[e], k[(e + 1) % 3]), std::max(k[e], k[(e + 1) % 3])}];
    }
    for (const auto& [e, n] : edges)
        EXPECT_EQ(n, 2) << "boundary edge " << e[0] << "-" << e[1] << " on " << n << " triangles";
    return out;
}

mio::OctreeMeshOptions om_tet(double MaxSize, double MinSize, bool Uniform) {
    mio::OctreeMeshOptions o;
    o.mOutput = mio::OctreeMeshOutput::Tetra;
    o.mMaxSize = MaxSize;
    o.mMinSize = MinSize;
    o.mPaddingRelative = 0.05;
    if (Uniform) {
        o.mSkinGrading = mio::OctreeMeshSkinGrading::Uniform;
        o.mSurfaceSize = 2.0 * MinSize;
    }
    return o;
}

// The volume a closed, outward surface encloses.
double om_enclosed(const mio::Mesh& rSkin) {
    const auto soup = mio::detail::build_triangle_soup(rSkin, "");
    double v = 0.0;
    for (std::size_t t = 0; t < soup.NumTriangles(); ++t)
        v += mio::detail::vec3_dot(
                 soup.mCorners[3 * t],
                 mio::detail::vec3_cross(soup.mCorners[3 * t + 1], soup.mCorners[3 * t + 2])) /
             6.0;
    return v;
}

}  // namespace

TEST(OctreeMeshTet, UniformSkinsKeepTheStuffingAngleBounds) {
    const mio::Mesh skins[] = {mio::icosa_sphere(1.0, 4), analytic::torus(1.0, 0.4, 64, 32)};
    for (const mio::Mesh& skin : skins) {
        auto o = om_tet(0.5, 0.05, true);
        o.mFitIterations = 0;
        const auto r = mio::octree_mesh(skin, o);
        const OmTetCheck c = om_check_tets(r.mMesh);
        EXPECT_GE(c.mMinDihedral, 10.7);
        EXPECT_LE(c.mMaxDihedral, 164.8);
        EXPECT_TRUE(r.mGuaranteeApplies);
        EXPECT_EQ(r.mNumInverted, 0);
        EXPECT_NEAR(c.mVolume, om_enclosed(skin), 0.02 * om_enclosed(skin));
        EXPECT_GT(r.mNumWarped, 0);
        EXPECT_NEAR(r.mMinDihedral, c.mMinDihedral, 1e-9);
    }
}

TEST(OctreeMeshTet, CubeCornersAndEdgesArePreserved) {
    const double lo[3] = {-0.61, -0.67, -0.7}, L = 1.37;
    const mio::Mesh cube =
        mio::convert_cells(
            mio::extract_surface(mio::box({{8, 8, 8}}, {{L, L, L}}, {{lo[0], lo[1], lo[2]}})),
            {mio::ConvertCellsMode::Simplexify, false})
            .mMesh;
    for (bool uniform : {true, false}) {
        auto o = om_tet(0.5, 0.05, uniform);
        const auto r = mio::octree_mesh(cube, o);
        om_check_tets(r.mMesh);
        EXPECT_EQ(r.mNumCorners, 8);
        EXPECT_EQ(r.mNumCornersSnapped, 8);
        EXPECT_EQ(r.mNumFeatureChains, 12);
        EXPECT_EQ(r.mNumFeatureChainsMissed, 0);
        EXPECT_FALSE(r.mGuaranteeApplies);
        const double* p = r.mMesh.Points().As<double>();
        for (int i = 0; i < 8; ++i) {
            const Vec3 corner{lo[0] + ((i & 1) ? L : 0.0), lo[1] + ((i & 2) ? L : 0.0),
                              lo[2] + ((i & 4) ? L : 0.0)};
            double best = 1e300;
            for (std::size_t v = 0; v < r.mMesh.NumPoints(); ++v)
                best = std::min(best, mio::detail::vec3_norm(mio::detail::vec3_sub(
                                          Vec3{p[3 * v], p[3 * v + 1], p[3 * v + 2]}, corner)));
            EXPECT_LE(best, 1e-12);
        }
        // The feature-edge nodes lie on the cube's edges.
        for (std::int64_t v :
             om_region(r.mMesh, mio::kOctreeMeshFeatureEdgesRegion, mio::RegionKind::Point)) {
            int on = 0;
            for (int k = 0; k < 3; ++k) {
                const double x = p[3 * v + k];
                on += std::abs(x - lo[k]) < 1e-12 || std::abs(x - lo[k] - L) < 1e-12;
            }
            EXPECT_GE(on, 2) << "node " << v;
        }
    }
}

TEST(OctreeMeshTet, CleavingNestedSpheresGivesTwoConformingMaterials) {
    const mio::Mesh a = mio::icosa_sphere(1.0, 4);
    const mio::Mesh b = mio::icosa_sphere(0.5, 3);
    const mio::Mesh both = mio::merge({&a, &b}, mio::MergeOptions{}).mMesh;
    auto o = om_tet(0.4, 0.05, false);
    o.mConformity = mio::OctreeMeshConformity::Cleaving;
    const auto r = mio::octree_mesh(both, o);
    const OmTetCheck c = om_check_tets(r.mMesh);
    EXPECT_NEAR(c.mVolume, om_enclosed(a), 0.02 * om_enclosed(a));
    const auto m1 = om_region(r.mMesh, "material_1", mio::RegionKind::Cell);
    const auto m2 = om_region(r.mMesh, "material_2", mio::RegionKind::Cell);
    EXPECT_GT(m1.size(), 0u);
    EXPECT_GT(m2.size(), 0u);
    EXPECT_EQ(m1.size() + m2.size(), r.mMesh.Cells(0).NumCells());
    EXPECT_EQ(r.mNumInverted, 0);
}

TEST(OctreeMeshTet, KeepOutsideFillsTheLattice) {
    const mio::Mesh skin = mio::icosa_sphere(1.0, 3);
    auto o = om_tet(0.4, 0.05, false);
    o.mKeepOutside = true;
    const auto r = mio::octree_mesh(skin, o);
    const OmTetCheck c = om_check_tets(r.mMesh);
    double box = 1.0;
    for (int k = 0; k < 3; ++k)
        box *= static_cast<double>(r.mRootCells[static_cast<std::size_t>(k)]) * r.mRootSize;
    EXPECT_NEAR(c.mVolume, box, 1e-9 * box);
    const auto in = om_region(r.mMesh, mio::kOctreeMeshInsideRegion, mio::RegionKind::Cell);
    const auto out = om_region(r.mMesh, mio::kOctreeMeshOutsideRegion, mio::RegionKind::Cell);
    EXPECT_EQ(in.size() + out.size(), r.mMesh.Cells(0).NumCells());
}

TEST(OctreeMeshTet, DirtyAndGradedSkinsStillGiveValidMeshes) {
    const mio::Mesh sphere = mio::icosa_sphere(1.0, 4);
    std::vector<std::vector<double>> pts;
    const double* p = sphere.Points().As<double>();
    for (std::size_t i = 0; i < sphere.NumPoints(); ++i)
        pts.push_back({p[3 * i], p[3 * i + 1], p[3 * i + 2]});
    std::vector<std::vector<std::int64_t>> holed;
    const std::int64_t* c = sphere.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < sphere.Cells(0).NumCells(); ++t)
        if (t % 13 != 0)
            holed.push_back({c[3 * t], c[3 * t + 1], c[3 * t + 2]});
    const auto h =
        mio::octree_mesh(mt::make_mesh(pts, "triangle", holed), om_tet(0.4, 0.05, false));
    const OmTetCheck hc = om_check_tets(h.mMesh);
    EXPECT_NEAR(hc.mVolume, om_enclosed(sphere), 0.05 * om_enclosed(sphere));
    // A graded skin: sizes from a point array varying along the skin.
    mio::Mesh graded = mio::icosa_sphere(1.0, 4);
    mio::NDArray size = mio::NDArray::Uninit(mio::DType::Float64, {graded.NumPoints()});
    const double* gp = graded.Points().As<double>();
    for (std::size_t i = 0; i < graded.NumPoints(); ++i)
        size.As<double>()[i] = 0.05 + 0.15 * (gp[3 * i] + 1.0);
    graded.AddPointData("h", std::move(size));
    auto o = om_tet(0.4, 0.04, false);
    o.mSizeArray = "h";
    const auto g = mio::octree_mesh(graded, o);
    const OmTetCheck gc = om_check_tets(g.mMesh);
    EXPECT_GT(gc.mMinDihedral, 0.0);
    EXPECT_GT(g.mFinestLevel, g.mLeavesPerLevel.size() > 0 ? 0 : -1);
    EXPECT_EQ(g.mNumInverted, 0);
}

TEST(OctreeMeshTet, IsTheSameOnEveryThreadCount) {
    const mio::Mesh skin = analytic::torus(1.0, 0.4, 32, 16);
    std::uint64_t ref = 0;
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    const int saved = mio_test::max_threads();
    for (int threads : {1, 2, 4}) {
        mio_test::set_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const auto r = mio::octree_mesh(skin, om_tet(0.5, 0.05, false));
        mio::bench::MeshDigest d;
        d.Of(r.mMesh);
        if (ref == 0)
            ref = d.Value();
        EXPECT_EQ(d.Value(), ref);
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(saved);
#endif
}

TEST(OctreeMeshTet, OptimisationNeverWorsensTheExtremeAngles) {
    const double lo[3] = {-0.61, -0.67, -0.7}, L = 1.37;
    const mio::Mesh skins[] = {
        mio::icosa_sphere(1.0, 4), analytic::torus(1.0, 0.4, 64, 32),
        mio::convert_cells(
            mio::extract_surface(mio::box({{8, 8, 8}}, {{L, L, L}}, {{lo[0], lo[1], lo[2]}})),
            {mio::ConvertCellsMode::Simplexify, false})
            .mMesh};
    for (const mio::Mesh& skin : skins)
        for (bool uniform : {true, false}) {
            auto o = om_tet(0.5, 0.05, uniform);
            o.mFitIterations = 0;
            const auto before = mio::octree_mesh(skin, o);
            o.mFitIterations = 5;
            const auto after = mio::octree_mesh(skin, o);
            const OmTetCheck c = om_check_tets(after.mMesh);
            EXPECT_EQ(after.mNumInverted, 0);
            EXPECT_GE(c.mMinDihedral, before.mMinDihedral - 1e-9);
            EXPECT_LE(c.mMaxDihedral, before.mMaxDihedral + 1e-9);
            EXPECT_GT(after.mNumOptimizedMoves + after.mNumFlips, 0);
            EXPECT_EQ(before.mNumOptimizedMoves, 0);
            EXPECT_NEAR(c.mVolume, om_enclosed(skin), 0.02 * std::abs(om_enclosed(skin)));
            // A uniform skin keeps the stuffing bounds through optimisation.
            if (uniform && before.mGuaranteeApplies) {
                EXPECT_GE(c.mMinDihedral, 10.7);
                EXPECT_LE(c.mMaxDihedral, 164.8);
                EXPECT_TRUE(after.mGuaranteeApplies);
            }
        }
}

TEST(OctreeMeshTet, OptimisedMeshesHaveNoInvertedOrOverlappingTetrahedra) {
    const mio::Mesh skins[] = {mio::icosa_sphere(1.0, 4), analytic::torus(1.0, 0.4, 64, 32)};
    for (const mio::Mesh& skin : skins) {
        // om_check_tets requires every volume to be positive.
        const auto r = mio::octree_mesh(skin, om_tet(0.5, 0.05, true));
        om_check_tets(r.mMesh);
        EXPECT_EQ(r.mNumInverted, 0);
        EXPECT_EQ(r.mNumOverlaps, 0);
    }
}

TEST(OctreeMeshTet, TheVolumeFixCanBeSwitchedOffAndTuned) {
    const mio::Mesh skin = analytic::torus(1.0, 0.4, 64, 32);
    auto options = om_tet(0.5, 0.05, true);
    options.mVolumeFix = false;
    const auto off = mio::octree_mesh(skin, options);
    EXPECT_EQ(off.mNumCollapsed + off.mNumCapsRemoved + off.mNumSplit + off.mNumSlivers +
                  off.mNumUntangled + off.mNumReordered + off.mNumOverlaps,
              0);
    // A zero-volume factor of 0 collapses nothing (no tetrahedron has a volume under 0 times
    // an edge cubed beyond the coordinate noise), and the result is still valid.
    options.mVolumeFix = true;
    options.mZeroVolumeFactor = 0.0;
    const auto tuned = mio::octree_mesh(skin, options);
    om_check_tets(tuned.mMesh);
    EXPECT_EQ(tuned.mNumInverted, 0);
}

TEST(OctreeMeshTet, ReportsTheBoundarysDistanceAndManifoldness) {
    const mio::Mesh skin = mio::icosa_sphere(1.0, 4);
    const auto r = mio::octree_mesh(skin, om_tet(0.5, 0.05, false));
    om_check_tets(r.mMesh);
    EXPECT_EQ(r.mNumBoundaryEdgesOpen, 0);
    EXPECT_EQ(r.mNumNonManifoldBoundaryEdges, 0);
    EXPECT_GT(r.mHausdorffToSkin, 0.0);
    EXPECT_LT(r.mHausdorffToSkin, 0.5 * r.mMinSize + 0.02);
    EXPECT_GT(r.mQuality.mNumCells, 0);
    for (const char* stage : {"skin", "size_field", "octree", "stuffing", "optimize", "output"}) {
        bool found = false;
        for (const auto& [name, seconds] : r.mTimings)
            found = found || name == stage;
        EXPECT_TRUE(found) << stage;
    }
}

namespace {

struct OmHexCheck {
    double mVolume = 0.0;
    double mMinScaledJacobian = 1.0;
};

// Every cell a positive hexahedron; the boundary closed and edge-manifold.
OmHexCheck om_check_hexes(const mio::Mesh& rM) {
    OmHexCheck out;
    EXPECT_EQ(rM.NumCellBlocks(), 1u);
    const auto block = rM.Cells(0);
    EXPECT_EQ(block.Type(), "hexahedron");
    mio::detail::HexMesh hm;
    const double* p = rM.Points().As<double>();
    for (std::size_t v = 0; v < rM.NumPoints(); ++v)
        hm.mPoints.push_back({p[3 * v], p[3 * v + 1], p[3 * v + 2]});
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    for (std::size_t h = 0; h < block.NumCells(); ++h) {
        std::array<std::int64_t, 8> hex{};
        for (std::size_t k = 0; k < 8; ++k)
            hex[k] = c[8 * h + k];
        hm.mHexes.push_back(hex);
        out.mMinScaledJacobian =
            std::min(out.mMinScaledJacobian, mio::detail::hex_min_scaled_jacobian(hm.mPoints, hex));
    }
    const mio::detail::HexBoundary b = mio::detail::hex_boundary(hm);
    EXPECT_EQ(b.mNumOpenEdges, 0);
    EXPECT_EQ(b.mNumNonManifoldEdges, 0);
    for (const auto& q : b.mQuads) {
        const Vec3& a = hm.mPoints[static_cast<std::size_t>(q[0])];
        out.mVolume += (mio::detail::triple_product(a, hm.mPoints[static_cast<std::size_t>(q[1])],
                                                    hm.mPoints[static_cast<std::size_t>(q[2])]) +
                        mio::detail::triple_product(a, hm.mPoints[static_cast<std::size_t>(q[2])],
                                                    hm.mPoints[static_cast<std::size_t>(q[3])])) /
                       6.0;
    }
    return out;
}

mio::OctreeMeshOptions om_hex(double SurfaceSize) {
    mio::OctreeMeshOptions o;
    o.mOutput = mio::OctreeMeshOutput::Hexa;
    o.mSurfaceSize = SurfaceSize;
    return o;
}

}  // namespace

TEST(OctreeMeshHex, SpheresAndToriAreAllHexahedraFittedToTheSkin) {
    const mio::Mesh skins[] = {mio::icosa_sphere(1.0, 4), analytic::torus(1.0, 0.4, 64, 32)};
    for (const mio::Mesh& skin : skins)
        for (bool pillow : {true, false}) {
            auto o = om_hex(0.15);
            o.mPillow = pillow;
            const auto r = mio::octree_mesh(skin, o);
            const OmHexCheck c = om_check_hexes(r.mMesh);
            if (pillow) {
                EXPECT_GE(c.mMinScaledJacobian, 0.4);
                EXPECT_TRUE(r.mWarnings.empty());
            } else if (c.mMinScaledJacobian < 0.2) {
                // Without a pillow, cells with several faces on the skin can
                // stay poor; that is reported, never silent.
                EXPECT_FALSE(r.mWarnings.empty());
            }
            EXPECT_EQ(r.mNumInverted, 0);
            EXPECT_EQ(r.mNumHexahedra, static_cast<std::int64_t>(r.mMesh.Cells(0).NumCells()));
            EXPECT_EQ(r.mNumBoundaryEdgesOpen, 0);
            EXPECT_EQ(r.mNumNonManifoldBoundaryEdges, 0);
            EXPECT_LT(r.mHausdorffToSkin, 0.05);
            EXPECT_NEAR(c.mVolume, om_enclosed(skin), 0.03 * om_enclosed(skin));
            EXPECT_NEAR(r.mMinScaledJacobian, c.mMinScaledJacobian, 0.05);
        }
}

TEST(OctreeMeshHex, CubeCornersAndEdgesArePreserved) {
    const double lo[3] = {-0.61, -0.67, -0.7}, L = 1.37;
    const mio::Mesh cube =
        mio::convert_cells(
            mio::extract_surface(mio::box({{8, 8, 8}}, {{L, L, L}}, {{lo[0], lo[1], lo[2]}})),
            {mio::ConvertCellsMode::Simplexify, false})
            .mMesh;
    const auto r = mio::octree_mesh(cube, om_hex(0.15));
    const OmHexCheck c = om_check_hexes(r.mMesh);
    EXPECT_GE(c.mMinScaledJacobian, 0.5);
    EXPECT_EQ(r.mNumCorners, 8);
    EXPECT_EQ(r.mNumCornersSnapped, 8);
    EXPECT_EQ(r.mNumFeatureChains, 12);
    EXPECT_EQ(r.mNumFeatureChainsMissed, 0);
    EXPECT_NEAR(c.mVolume, L * L * L, 1e-9);
    const double* p = r.mMesh.Points().As<double>();
    for (std::int64_t v :
         om_region(r.mMesh, mio::kOctreeMeshCornersRegion, mio::RegionKind::Point)) {
        int on = 0;
        for (int k = 0; k < 3; ++k) {
            const double x = p[3 * v + k];
            on += std::abs(x - lo[k]) < 1e-12 || std::abs(x - lo[k] - L) < 1e-12;
        }
        EXPECT_EQ(on, 3) << "corner node " << v;
    }
    // The pillow layer is tagged, and the boundary region is its outer faces.
    const auto boundary = om_region(r.mMesh, mio::kOctreeMeshBoundaryRegion, mio::RegionKind::Side);
    EXPECT_FALSE(boundary.empty());
    for (std::size_t s = 1; s < boundary.size(); s += 2)
        EXPECT_EQ(boundary[s], 5);
}

TEST(OctreeMeshHex, IsTheSameOnEveryThreadCount) {
    const mio::Mesh skin = analytic::torus(1.0, 0.4, 32, 16);
    std::uint64_t ref = 0;
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    const int saved = mio_test::max_threads();
    for (int threads : {1, 2, 4}) {
        mio_test::set_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const auto r = mio::octree_mesh(skin, om_hex(0.12));
        mio::bench::MeshDigest d;
        d.Of(r.mMesh);
        if (ref == 0)
            ref = d.Value();
        EXPECT_EQ(d.Value(), ref);
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(saved);
#endif
}

TEST(OctreeMeshHex, RefusesAGridAboveMaxCellsByName) {
    auto o = om_hex(0.01);
    o.mMaxCells = 1000;
    try {
        (void)mio::octree_mesh(mio::icosa_sphere(1.0, 3), o);
        FAIL() << "expected a refusal";
    } catch (const std::invalid_argument& rE) {
        EXPECT_NE(std::string(rE.what()).find("max_cells"), std::string::npos);
    }
}

TEST(OctreeMeshHex, GradedHexahedraFollowTheOctree) {
    // The sphere at a fine curvature angle and a capped cylinder whose skin
    // grades from its rims to its caps: the graded dual has fewer hexahedra
    // than one grid at the skin's size, none inverted, a closed manifold
    // boundary on the skin, and every corner of scaled Jacobian 0.2 or more.
    struct Case {
        mio::Mesh mSkin;
        mio::OctreeMeshOptions mOptions;
    };
    std::vector<Case> cases;
    {
        auto o = om_hex(0.0);
        o.mCurvatureAngle = 5.0;
        o.mMaxSize = 0.5;
        cases.push_back({mio::icosa_sphere(1.0, 5), o});
    }
    {
        auto o = om_hex(0.0);
        o.mMinSize = 0.05;
        cases.push_back({analytic::capped_cylinder(0.5, 2.0, 6, 48), o});
    }
    for (const Case& c : cases) {
        auto o = c.mOptions;
        o.mHexGrading = mio::OctreeMeshHexGrading::Graded;
        const auto r = mio::octree_mesh(c.mSkin, o);
        const auto u = mio::octree_mesh(c.mSkin, c.mOptions);
        const OmHexCheck h = om_check_hexes(r.mMesh);
        EXPECT_LT(r.mNumHexahedra, u.mNumHexahedra);
        EXPECT_EQ(r.mNumInverted, 0);
        EXPECT_GE(h.mMinScaledJacobian, 0.2);
        EXPECT_EQ(r.mNumBoundaryEdgesOpen, 0);
        EXPECT_EQ(r.mNumNonManifoldBoundaryEdges, 0);
        EXPECT_LT(r.mHausdorffToSkin, 0.05);
        EXPECT_NEAR(h.mVolume, om_enclosed(c.mSkin), 0.03 * om_enclosed(c.mSkin));
        const mio::ConformityResult cc = mio::check_conformity(r.mMesh);
        EXPECT_EQ(cc.mNumHangingNodes, 0);
        EXPECT_EQ(cc.mNumNonManifoldFacets, 0);
        // More than one level, and the levels recorded per cell.
        const mio::NDArray& lv = r.mMesh.CellData(mio::kOctreeLevelName, 0);
        const std::int64_t* l = lv.As<std::int64_t>();
        std::set<std::int64_t> levels(l, l + r.mMesh.Cells(0).NumCells());
        EXPECT_GT(levels.size(), 1u);
    }
}

TEST(OctreeMeshHex, GradedHexahedraAreTheSameOnEveryThreadCount) {
    auto o = om_hex(0.0);
    o.mMinSize = 0.05;
    o.mHexGrading = mio::OctreeMeshHexGrading::Graded;
    const mio::Mesh skin = analytic::capped_cylinder(0.5, 2.0, 6, 48);
    std::uint64_t ref = 0;
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    const int saved = mio_test::max_threads();
    for (int threads : {1, 2, 4}) {
        mio_test::set_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const std::uint64_t d = om_digest(mio::octree_mesh(skin, o).mMesh);
        if (ref == 0)
            ref = d;
        EXPECT_EQ(d, ref);
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(saved);
#endif
}

namespace {

// Set the hexahedral fit's internal test hook for a scope ("body_fitted": the
// 18.0.0 fit alone; "robust": the robust fit alone).
struct OmHexFit {
    explicit OmHexFit(const char* pValue) { om_setenv(pValue); }
    ~OmHexFit() { om_setenv(nullptr); }
    static void om_setenv(const char* pValue) {
#ifdef _WIN32
        _putenv_s("MESHIOPLUSPLUS_OCTREE_MESH_HEX_FIT", pValue != nullptr ? pValue : "");
#else
        if (pValue != nullptr)
            setenv("MESHIOPLUSPLUS_OCTREE_MESH_HEX_FIT", pValue, 1);
        else
            unsetenv("MESHIOPLUSPLUS_OCTREE_MESH_HEX_FIT");
#endif
    }
};

// The hexahedra's corners with a non-positive scaled Jacobian.
std::int64_t om_flat_corners(const mio::Mesh& rM) {
    std::vector<Vec3> pts;
    const double* p = rM.Points().As<double>();
    for (std::size_t v = 0; v < rM.NumPoints(); ++v)
        pts.push_back({p[3 * v], p[3 * v + 1], p[3 * v + 2]});
    const auto block = rM.Cells(0);
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    std::int64_t flat = 0;
    for (std::size_t h = 0; h < block.NumCells(); ++h) {
        std::array<std::int64_t, 8> hex{};
        for (std::size_t k = 0; k < 8; ++k)
            hex[k] = c[8 * h + k];
        for (double j : mio::detail::hex_corner_scaled_jacobians(pts, hex))
            flat += j <= 0.0;
    }
    return flat;
}

// A unit box on the lattice's planes, its faces split into triangles.
mio::Mesh om_unit_box() {
    return mio::convert_cells(mio::extract_surface(mio::box({{4, 4, 4}}, {{1.0, 1.0, 1.0}})),
                              {mio::ConvertCellsMode::Simplexify, false})
        .mMesh;
}

}  // namespace

// The capped cylinder of the hexahedral table (radius 0.5, height 2,
// min_size 0.05): its rims are chains a staircase follows by turning inside
// quadrilaterals, which the robust fit releases from the chain rather than
// flatten, and every corner keeps a scaled Jacobian of at least 0.2.
TEST(OctreeMeshHex, ACappedCylindersHexahedraKeepAFifth) {
    const mio::Mesh skin = analytic::capped_cylinder(0.5, 2.0, 6, 48);
    for (auto grading : {mio::OctreeMeshHexGrading::Uniform, mio::OctreeMeshHexGrading::Graded}) {
        auto o = om_hex(0.0);
        o.mMinSize = 0.05;
        o.mHexGrading = grading;
        const auto r = mio::octree_mesh(skin, o);
        const OmHexCheck c = om_check_hexes(r.mMesh);
        EXPECT_GE(c.mMinScaledJacobian, 0.2) << static_cast<int>(grading);
        EXPECT_EQ(om_flat_corners(r.mMesh), 0);
        EXPECT_EQ(r.mNumInverted, 0);
        EXPECT_EQ(r.mNumCornersSnapped, r.mNumCorners);
        EXPECT_EQ(r.mNumFeatureChainsMissed, 0);
        EXPECT_LT(r.mHausdorffToSkin, 0.01);
        EXPECT_NEAR(c.mVolume, om_enclosed(skin), 0.03 * om_enclosed(skin));
        EXPECT_TRUE(r.mWarnings.empty());
    }
}

// Features off (feature_angle 0, or preserve_features false): near a box's
// edges and corners a graded boundary node is about as near two or three
// faces, and its closest point put a face's pillow copies on the side face,
// flattening the pillow cells (the 18.0.0 fit, kept here as the hook, has
// hundreds of such corners). Placed along the boundary's own normal, none is
// flat.
TEST(OctreeMeshHex, ABoxWithoutFeaturesHasNoFlatCorner) {
    const mio::Mesh box = om_unit_box();
    for (auto grading : {mio::OctreeMeshHexGrading::Uniform, mio::OctreeMeshHexGrading::Graded})
        for (bool angle_off : {true, false}) {
            auto o = om_hex(0.0);
            o.mMinSize = 0.05;
            o.mHexGrading = grading;
            if (angle_off)
                o.mFeatureAngle = 0.0;
            else
                o.mPreserveFeatures = false;
            const auto r = mio::octree_mesh(box, o);
            EXPECT_EQ(r.mNumCorners, 0);
            EXPECT_EQ(om_flat_corners(r.mMesh), 0) << static_cast<int>(grading) << angle_off;
            const OmHexCheck c = om_check_hexes(r.mMesh);
            EXPECT_GE(c.mMinScaledJacobian, 0.2);
            EXPECT_NEAR(c.mVolume, 1.0, 0.03);
            EXPECT_TRUE(r.mWarnings.empty());
            if (grading == mio::OctreeMeshHexGrading::Graded) {
                const OmHexFit hook("body_fitted");
                EXPECT_GT(om_flat_corners(mio::octree_mesh(box, o).mMesh), 100);
            }
        }
}

// With features on, a chain node beside a corner went to the corner itself
// (two nodes on one point), and the same box kept a few flat corners.
TEST(OctreeMeshHex, ABoxWithFeaturesHasNoFlatCorner) {
    auto o = om_hex(0.0);
    o.mMinSize = 0.05;
    o.mHexGrading = mio::OctreeMeshHexGrading::Graded;
    const auto r = mio::octree_mesh(om_unit_box(), o);
    EXPECT_EQ(om_flat_corners(r.mMesh), 0);
    EXPECT_EQ(r.mNumCornersSnapped, 8);
    EXPECT_EQ(r.mNumFeatureChainsMissed, 0);
    EXPECT_GE(om_check_hexes(r.mMesh).mMinScaledJacobian, 0.2);
}

// The robust fit is compared with the body-fitted fit of 18.0.0 when it
// misses a fifth, and the better is kept: the hook that forces the old fit
// gives its output bit for bit (digests recorded from 18.0.0).
TEST(OctreeMeshHex, TheBodyFittedHookIsTheFitOfVersion18) {
    const OmHexFit hook("body_fitted");
    const auto torus = mio::octree_mesh(analytic::torus(1.0, 0.4, 32, 16), om_hex(0.12));
    auto o = om_hex(0.0);
    o.mMinSize = 0.05;
    o.mHexGrading = mio::OctreeMeshHexGrading::Graded;
    const auto cylinder = mio::octree_mesh(analytic::capped_cylinder(0.5, 2.0, 6, 48), o);
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    EXPECT_EQ(om_digest(torus.mMesh), 0x7d1e198a5fb3cc33ULL) << std::hex << om_digest(torus.mMesh);
    EXPECT_EQ(om_digest(cylinder.mMesh), 0x1c72bfbe7f7003f1ULL)
        << std::hex << om_digest(cylinder.mMesh);
#else
    EXPECT_GT(torus.mNumHexahedra, 0);
    EXPECT_GT(cylinder.mNumHexahedra, 0);
#endif
}

namespace {

// The volume of each material of a hexahedral output (each hexahedron as six
// tetrahedra around its 0-6 diagonal), indexed by material.
std::vector<double> om_material_volumes(const mio::Mesh& rM) {
    const double* p = rM.Points().As<double>();
    const std::int64_t* c = rM.Cells(0).Conn().As<std::int64_t>();
    const std::int64_t* mat = rM.CellData(mio::kOctreeMeshMaterialName, 0).As<std::int64_t>();
    std::vector<double> out;
    static constexpr int kTets[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
                                        {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    for (std::size_t h = 0; h < rM.Cells(0).NumCells(); ++h) {
        const auto m = static_cast<std::size_t>(mat[h]);
        if (out.size() <= m)
            out.resize(m + 1, 0.0);
        auto at = [&](int K) {
            const std::int64_t v = c[8 * h + static_cast<std::size_t>(K)];
            return Vec3{p[3 * v], p[3 * v + 1], p[3 * v + 2]};
        };
        for (const auto& t : kTets) {
            const Vec3 o = at(t[0]);
            out[m] += mio::detail::triple_product(mio::detail::vec3_sub(at(t[1]), o),
                                                  mio::detail::vec3_sub(at(t[2]), o),
                                                  mio::detail::vec3_sub(at(t[3]), o)) /
                      6.0;
        }
    }
    return out;
}

mio::Mesh om_nested_spheres(const std::vector<double>& rRadii) {
    std::vector<mio::Mesh> parts;
    for (double r : rRadii)
        parts.push_back(mio::icosa_sphere(r, 4));
    std::vector<const mio::Mesh*> ptrs;
    for (const auto& m : parts)
        ptrs.push_back(&m);
    return mio::merge(ptrs, mio::MergeOptions{}).mMesh;
}

}  // namespace

// Roadmap 7.9.3: `conformity = Cleaving` gives the hexahedral output one
// material per winding level, each interface pillowed on both sides: the
// materials conform (no hanging node, no facet in three cells), each takes
// its shell's volume, and every corner keeps a fifth.
TEST(OctreeMeshHex, CleavingNestedShellsGivesConformingMaterials) {
    const double pi = kPi;
    for (const std::vector<double>& radii :
         {std::vector<double>{1.0, 0.5}, std::vector<double>{1.0, 0.6, 0.3}}) {
        const mio::Mesh skin = om_nested_spheres(radii);
        for (auto grading :
             {mio::OctreeMeshHexGrading::Uniform, mio::OctreeMeshHexGrading::Graded}) {
            auto o = om_hex(grading == mio::OctreeMeshHexGrading::Uniform ? 0.1 : 0.0);
            o.mMinSize = radii.size() > 2 ? 0.05 : 0.1;
            if (grading == mio::OctreeMeshHexGrading::Uniform)
                o.mSurfaceSize = o.mMinSize;
            o.mHexGrading = grading;
            o.mConformity = mio::OctreeMeshConformity::Cleaving;
            const auto r = mio::octree_mesh(skin, o);
            const OmHexCheck c = om_check_hexes(r.mMesh);
            EXPECT_GE(c.mMinScaledJacobian, 0.2) << radii.size() << static_cast<int>(grading);
            EXPECT_EQ(om_flat_corners(r.mMesh), 0);
            EXPECT_EQ(r.mNumInverted, 0);
            EXPECT_TRUE(r.mWarnings.empty());
            EXPECT_LT(r.mHausdorffToSkin, 0.05);
            const mio::ConformityResult cc = mio::check_conformity(r.mMesh);
            EXPECT_EQ(cc.mNumHangingNodes, 0);
            EXPECT_EQ(cc.mNumNonManifoldFacets, 0);
            const std::vector<double> vol = om_material_volumes(r.mMesh);
            ASSERT_EQ(vol.size(), radii.size() + 1);
            for (std::size_t k = 0; k < radii.size(); ++k) {
                const double ball = 4.0 / 3.0 * pi * std::pow(radii[k], 3);
                const double hole =
                    k + 1 < radii.size() ? 4.0 / 3.0 * pi * std::pow(radii[k + 1], 3) : 0.0;
                EXPECT_NEAR(vol[k + 1], ball - hole, 0.03 * (ball - hole)) << k;
                EXPECT_FALSE(om_region(r.mMesh, ("material_" + std::to_string(k + 1)).c_str(),
                                       mio::RegionKind::Cell)
                                 .empty());
            }
        }
    }
}

// A skin with one shell has one material: cleaving gives stuffing's bytes.
TEST(OctreeMeshHex, CleavingOneShellIsStuffing) {
    for (auto grading : {mio::OctreeMeshHexGrading::Uniform, mio::OctreeMeshHexGrading::Graded}) {
        auto o = om_hex(0.15);
        o.mHexGrading = grading;
        const std::uint64_t plain = om_digest(mio::octree_mesh(mio::icosa_sphere(1.0, 3), o).mMesh);
        o.mConformity = mio::OctreeMeshConformity::Cleaving;
        EXPECT_EQ(om_digest(mio::octree_mesh(mio::icosa_sphere(1.0, 3), o).mMesh), plain);
    }
}

// Graded hexahedra used to drop the dual centres within a tenth of a cell of
// any skin, which opened a cavity at a nested shell: the inner sphere became
// an inner boundary, pillowed and folded. Only the boundary's centres drop now.
TEST(OctreeMeshHex, GradedNestedShellsLeaveNoCavity) {
    auto o = om_hex(0.0);
    o.mMinSize = 0.1;
    o.mHexGrading = mio::OctreeMeshHexGrading::Graded;
    const mio::Mesh skin = om_nested_spheres({1.0, 0.5});
    const auto r = mio::octree_mesh(skin, o);
    const OmHexCheck c = om_check_hexes(r.mMesh);
    EXPECT_GE(c.mMinScaledJacobian, 0.2);
    EXPECT_EQ(r.mNumInverted, 0);
    // The whole ball is meshed, and the boundary is the outer sphere only.
    EXPECT_NEAR(c.mVolume, 4.0 / 3.0 * kPi, 0.03 * 4.0 / 3.0 * kPi);
    const double* p = r.mMesh.Points().As<double>();
    const std::int64_t* conn = r.mMesh.Cells(0).Conn().As<std::int64_t>();
    const auto sides = om_region(r.mMesh, mio::kOctreeMeshBoundaryRegion, mio::RegionKind::Side);
    for (std::size_t s = 0; s < sides.size(); s += 2)
        for (int k = 0; k < 4; ++k) {
            const std::int64_t v =
                conn[8 * sides[s] + static_cast<std::size_t>(kFacet[sides[s + 1]][k])];
            EXPECT_NEAR(std::sqrt(p[3 * v] * p[3 * v] + p[3 * v + 1] * p[3 * v + 1] +
                                  p[3 * v + 2] * p[3 * v + 2]),
                        1.0, 0.05);
        }
}

TEST(OctreeMeshHex, CleavingIsTheSameOnEveryThreadCount) {
    auto o = om_hex(0.0);
    o.mMinSize = 0.1;
    o.mHexGrading = mio::OctreeMeshHexGrading::Graded;
    o.mConformity = mio::OctreeMeshConformity::Cleaving;
    const mio::Mesh skin = om_nested_spheres({1.0, 0.5});
    std::uint64_t ref = 0;
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    const int saved = mio_test::max_threads();
    for (int threads : {1, 2, 4}) {
        mio_test::set_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const std::uint64_t d = om_digest(mio::octree_mesh(skin, o).mMesh);
        if (ref == 0)
            ref = d;
        EXPECT_EQ(d, ref);
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(saved);
#endif
}

// Roadmap 7.9.4: a crease shorter than a cell (a scan's noise) is not
// snapped, and says so; the cube's edges, longer than a cell, still are.
TEST(OctreeMeshHex, CreasesShorterThanACellAreNotSnapped) {
    // A unit box whose face centre is pushed out: four short creases there.
    mio::Mesh box =
        mio::convert_cells(
            mio::extract_surface(mio::box({{8, 8, 8}}, {{1.0, 1.0, 1.0}}, {{-0.51, -0.53, -0.55}})),
            {mio::ConvertCellsMode::Simplexify, false})
            .mMesh;
    {
        const std::size_t np = box.NumPoints();
        mio::NDArray pts = mio::NDArray::Uninit(mio::DType::Float64, {np, std::size_t{3}});
        const double* src = box.Points().As<double>();
        double* p = pts.As<double>();
        std::copy(src, src + 3 * np, p);
        for (std::size_t v = 0; v < np; ++v)
            if (std::abs(p[3 * v] - 0.49) < 1e-9 && std::abs(p[3 * v + 1] + 0.03) < 1e-9 &&
                std::abs(p[3 * v + 2] + 0.05) < 1e-9)
                p[3 * v] += 0.1;
        box.AssignPoints(std::move(pts));
    }
    const auto r = mio::octree_mesh(box, om_hex(0.3));
    EXPECT_EQ(om_flat_corners(r.mMesh), 0);
    EXPECT_GE(r.mNumCornersSnapped, 8);
    bool said = false;
    for (const std::string& w : r.mWarnings)
        said = said || w.find("shorter than a cell") != std::string::npos;
    EXPECT_TRUE(said);
    EXPECT_GT(r.mNumFeatureChainsMissed, 0);
}

// ---- volume fixing and feature snapping on hand-made tetrahedra ----

namespace {

using meshioplusplus::detail::OctreeMeshTets;

double fix_test_volume(const OctreeMeshTets& rT, const std::array<std::int64_t, 4>& rTet) {
    const auto& a = rT.mPoints[static_cast<std::size_t>(rTet[0])];
    const auto sub = [&](const int I) {
        const auto& p = rT.mPoints[static_cast<std::size_t>(rTet[static_cast<std::size_t>(I)])];
        return std::array<double, 3>{p[0] - a[0], p[1] - a[1], p[2] - a[2]};
    };
    const auto u = sub(1), v = sub(2), w = sub(3);
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) +
            u[2] * (v[0] * w[1] - v[1] * w[0])) /
           6.0;
}

OctreeMeshTets fix_test_tets(const std::vector<std::array<double, 3>>& rPoints,
                             const std::vector<std::array<std::int64_t, 4>>& rTets) {
    OctreeMeshTets t;
    t.mPoints = rPoints;
    t.mTets = rTets;
    t.mMaterial.assign(rTets.size(), 1);
    t.mCut.assign(rTets.size(), 0);
    t.mLevel.assign(rTets.size(), 0);
    t.mShape.assign(rTets.size(), 0);
    t.mOnSkin.assign(rPoints.size(), 0);
    t.mFeature.assign(rPoints.size(), 0);
    t.mNumMaterials = 2;
    return t;
}

}  // namespace

TEST(OctreeMeshFixVolumes, UntanglesAFreeNode) {
    // Node 3 lies below the triangle 0-1-2 although its tetrahedron stands above it; the centroid of its neighbours (0, 1, 2 and 4) puts it back above
    OctreeMeshTets tets = fix_test_tets(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.33, 0.33, -0.1}, {0.33, 0.33, 0.4}},
        {{0, 1, 2, 3}, {0, 1, 2, 4}, {1, 2, 3, 4}});
    meshioplusplus::detail::OctreeMeshSetup setup;
    meshioplusplus::OctreeMeshOptions options;
    meshioplusplus::OctreeMeshResult result;
    meshioplusplus::detail::octree_mesh_fix_volumes(setup, options, result, tets);
    EXPECT_EQ(result.mNumUntangled, 1);
    EXPECT_EQ(result.mNumReordered, 0);
    EXPECT_EQ(tets.mTets.size(), 3u);
    for (const auto& t : tets.mTets)
        EXPECT_GT(fix_test_volume(tets, t), 0.0);
}

TEST(OctreeMeshFixVolumes, SplitsAFlatSliver) {
    // A flat tetrahedron 0-1-2-3 whose opposite edges 0-1 and 2-3 cross at (0.5, 0, 0), with all its nodes on features so that no collapse can remove it
    OctreeMeshTets tets = fix_test_tets(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, -0.5, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.2, 1.0}, {0.8, -0.2, 2.0}},
        {{0, 1, 2, 3}, {0, 1, 4, 5}, {2, 3, 4, 5}});
    // The neighbours around the two edges are positive
    for (std::size_t t = 1; t < 3; ++t)
        if (fix_test_volume(tets, tets.mTets[t]) < 0.0)
            std::swap(tets.mTets[t][0], tets.mTets[t][1]);
    tets.mFeature = {1, 1, 1, 1, 0, 0};
    meshioplusplus::detail::OctreeMeshSetup setup;
    meshioplusplus::OctreeMeshOptions options;
    meshioplusplus::OctreeMeshResult result;
    meshioplusplus::detail::octree_mesh_fix_volumes(setup, options, result, tets);
    EXPECT_EQ(result.mNumSlivers, 1);
    EXPECT_EQ(result.mNumReordered, 0);
    // The sliver is gone and each neighbour is replaced by two
    EXPECT_EQ(tets.mTets.size(), 4u);
    EXPECT_EQ(tets.mPoints.size(), 7u);
    for (const auto& t : tets.mTets)
        EXPECT_GT(fix_test_volume(tets, t), 0.0);
}

TEST(OctreeMeshFeatureSnap, KeepsGuardedTetrahedraPositive) {
    // A cube of 6x6x6 cells, each split in six tetrahedra around its diagonal, and a ridge that runs along the edge x = y = 0 slightly inside it
    constexpr int n = 6;
    const double h = 1.0 / n;
    const auto id = [&](const int I, const int J, const int K) {
        return static_cast<std::int64_t>((I * (n + 1) + J) * (n + 1) + K);
    };
    OctreeMeshTets m;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            for (int k = 0; k <= n; ++k)
                m.mPoints.push_back({i * h, j * h, k * h});
    const int path[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                for (const auto& path_i : path) {
                    int c[3] = {i, j, k};
                    std::array<std::int64_t, 4> tet{};
                    tet[0] = id(c[0], c[1], c[2]);
                    for (int s = 0; s < 3; ++s) {
                        ++c[path_i[s]];
                        tet[static_cast<std::size_t>(s) + 1] = id(c[0], c[1], c[2]);
                    }
                    if (fix_test_volume(m, tet) < 0.0)
                        std::swap(tet[0], tet[1]);
                    m.mTets.push_back(tet);
                }
    std::vector<std::array<double, 3>> ridge;
    std::vector<std::array<std::int64_t, 2>> edges;
    for (int k = 0; k <= n; ++k) {
        ridge.push_back({0.3 * h, 0.2 * h, k * h});
        if (k > 0)
            edges.push_back({k - 1, k});
    }
    const auto graph = meshioplusplus::detail::feature_graph(ridge, edges, 45.0);
    const std::vector<std::uint8_t> inside(m.mTets.size(), 1);
    std::vector<std::uint8_t> pinned(m.mPoints.size(), 0);
    std::vector<std::uint8_t> feature(m.mPoints.size(), 0);
    const double min_volume = 1.0e-12 * h * h * h;
    meshioplusplus::detail::feature_snap(
        m.mPoints, m.mTets, inside, inside, graph, ridge, [&](const std::array<double, 3>&) { return h; },
        min_volume, pinned, feature);
    EXPECT_GT(std::count(feature.begin(), feature.end(), 2) + std::count(feature.begin(), feature.end(), 1), 0);
    // However the moves went, no tetrahedron the snap was told to guard may end up flat or inverted
    for (const auto& t : m.mTets)
        EXPECT_GT(fix_test_volume(m, t), min_volume);
}
