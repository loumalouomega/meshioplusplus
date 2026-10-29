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

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

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
#include "../../src/cpp/src/detail/hex_fit.hpp"
#include "../../src/cpp/src/detail/tri_bvh.hpp"

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
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
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
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
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
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
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
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
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
    EXPECT_GE(c.mMinScaledJacobian, 0.2);
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
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
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
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
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
    // boundary on the skin, and every corner of positive scaled Jacobian.
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
        EXPECT_GT(h.mMinScaledJacobian, 0.0);
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
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const std::uint64_t d = om_digest(mio::octree_mesh(skin, o).mMesh);
        if (ref == 0)
            ref = d;
        EXPECT_EQ(d, ref);
    }
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
#endif
}
