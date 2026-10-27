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
// Tests for tetrahedralize (roadmap 6.7): tetrahedra filling a closed
// triangle surface through the Netgen backend, its triangles kept.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "brep_solids.hpp"
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/boundary_layer.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/tetrahedralize.hpp"

namespace {

namespace mio = meshioplusplus;

constexpr double kPi = 3.14159265358979323846;

std::uint64_t digest(const mio::Mesh& rM) {
    mio::bench::MeshDigest d;
    d.Of(rM);
    return d.Value();
}

/// The triangles of a mesh_cad_surface result, alone.
mio::Mesh skin(const mio::BrepInfo& rB, double H) {
    mio::MeshCadSurfaceOptions o;
    o.mTargetSize = H;
    o.mEdges = false;
    return mio::mesh_cad_surface(rB, o).mMesh;
}

/// The axis-aligned box [Lo, Hi]^3 as 12 triangles, facing out (or in).
void add_box(std::vector<double>& rX, std::vector<std::int64_t>& rT, double Lo, double Hi,
             bool Inward = false) {
    const auto base = static_cast<std::int64_t>(rX.size() / 3);
    const double p[8][3] = {{Lo, Lo, Lo}, {Hi, Lo, Lo}, {Hi, Hi, Lo}, {Lo, Hi, Lo},
                            {Lo, Lo, Hi}, {Hi, Lo, Hi}, {Hi, Hi, Hi}, {Lo, Hi, Hi}};
    for (const auto& r_p : p)
        rX.insert(rX.end(), r_p, r_p + 3);
    const std::int64_t t[12][3] = {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
                                   {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                                   {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    for (const auto& r_t : t) {
        rT.push_back(base + r_t[0]);
        rT.push_back(base + (Inward ? r_t[2] : r_t[1]));
        rT.push_back(base + (Inward ? r_t[1] : r_t[2]));
    }
}

mio::Mesh triangles(const std::vector<double>& rX, const std::vector<std::int64_t>& rT) {
    mio::Mesh m;
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {rX.size() / 3, std::size_t{3}});
    std::copy(rX.begin(), rX.end(), p.As<double>());
    m.AssignPoints(std::move(p));
    mio::NDArray c = mio::NDArray::Uninit(mio::DType::Int64, {rT.size() / 3, std::size_t{3}});
    std::copy(rT.begin(), rT.end(), c.As<std::int64_t>());
    m.AddCellBlock("triangle", std::move(c));
    return m;
}

/// The facets of the tetra block, sorted, with how many tetrahedra use each.
std::map<std::array<std::int64_t, 3>, int> facets(const mio::Mesh& rM) {
    std::map<std::array<std::int64_t, 3>, int> out;
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
        const auto cells = rM.Cells(b);
        if (cells.Type() != "tetra")
            continue;
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            std::array<std::int64_t, 4> t{};
            for (std::size_t k = 0; k < 4; ++k)
                t[k] = mio::detail::read_int(cells.Conn(), 4 * c + k);
            for (std::size_t skip = 0; skip < 4; ++skip) {
                std::array<std::int64_t, 3> f{};
                std::size_t j = 0;
                for (std::size_t k = 0; k < 4; ++k)
                    if (k != skip)
                        f[j++] = t[k];
                std::sort(f.begin(), f.end());
                ++out[f];
            }
        }
    }
    return out;
}

}  // namespace

// Netgen asked for where it is not built in is refused by name; Auto falls
// back to the native octree backend with a warning.
TEST(Tetrahedralize, AutoTakesNetgenWhenBuiltElseTheOctree) {
    std::vector<double> x;
    std::vector<std::int64_t> t;
    add_box(x, t, 0.0, 1.0);
    const auto r = mio::tetrahedralize(triangles(x, t));
    if (mio::tetrahedralize_has_netgen()) {
        EXPECT_EQ(r.mBackend, "netgen");
        EXPECT_TRUE(r.mWarnings.empty());
        return;
    }
    EXPECT_EQ(r.mBackend, "octree");
    ASSERT_FALSE(r.mWarnings.empty());
    EXPECT_NE(r.mWarnings[0].find("MESHIOPLUSPLUS_WITH_NETGEN"), std::string::npos);
    mio::TetrahedralizeOptions o;
    o.mBackend = mio::TetrahedralizeBackend::Netgen;
    try {
        (void)mio::tetrahedralize(triangles(x, t), o);
        FAIL() << "no Netgen, no throw";
    } catch (const mio::Unsupported& rE) {
        EXPECT_NE(std::string(rE.what()).find("MESHIOPLUSPLUS_WITH_NETGEN"), std::string::npos);
    }
}

// The native backend (roadmap 6.7): the sphere's skin filled by isosurface
// stuffing. Every tetrahedron is positive, the boundary is closed and within
// a small fraction of h of the skin (not kept: boundary_preserved is false),
// the volume is the skin's within that, and two runs agree.
TEST(Tetrahedralize, TheOctreeBackendFillsTheSphere) {
    const double h = 0.2;
    const mio::Mesh s = skin(brep_solids::sphere(1.0), h);
    mio::TetrahedralizeOptions o;
    o.mBackend = mio::TetrahedralizeBackend::Octree;
    const auto r = mio::tetrahedralize(s, o);
    std::printf(
        "[6.7] sphere, octree: %lld tetrahedra, boundary within %.4f, dihedral %.1f to %.1f, "
        "volume %.6f\n",
        static_cast<long long>(r.mNumTetrahedra), r.mHausdorffToSurface, r.mMinDihedralDeg,
        r.mMaxDihedralDeg, r.mVolume);
    EXPECT_EQ(r.mBackend, "octree");
    EXPECT_TRUE(r.mWarnings.empty());  // asked for by name
    EXPECT_FALSE(r.mBoundaryPreserved);
    EXPECT_GT(r.mHausdorffToSurface, 0.0);
    EXPECT_LT(r.mHausdorffToSurface, 0.1 * h);
    EXPECT_EQ(r.mNumRegions, 1);
    EXPECT_EQ(r.mNumNewPoints, r.mNumPoints);
    EXPECT_GT(r.mMinDihedralDeg, 0.0);
    EXPECT_LT(r.mMaxDihedralDeg, 180.0);
    // The skin's own enclosed volume, by the divergence theorem.
    double skin_volume = 0.0;
    {
        const double* p = s.Points().As<double>();
        const auto tri = s.Cells(0);
        for (std::size_t c = 0; c < tri.NumCells(); ++c) {
            const double* a = p + 3 * mio::detail::read_int(tri.Conn(), 3 * c);
            const double* b = p + 3 * mio::detail::read_int(tri.Conn(), 3 * c + 1);
            const double* d = p + 3 * mio::detail::read_int(tri.Conn(), 3 * c + 2);
            skin_volume +=
                (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
                 a[2] * (b[0] * d[1] - b[1] * d[0])) /
                6.0;
        }
    }
    EXPECT_NEAR(r.mVolume, skin_volume, 0.01 * skin_volume);
    // The triangle block is the tetrahedra's boundary: each a facet of one
    // tetrahedron, every other facet of two.
    const auto f = facets(r.mMesh);
    std::size_t once = 0;
    for (const auto& [k, n] : f) {
        EXPECT_LE(n, 2);
        once += n == 1 ? 1 : 0;
    }
    EXPECT_EQ(static_cast<std::int64_t>(once), r.mNumBoundaryTriangles);
    EXPECT_EQ(r.mMesh.Cells(1).NumCells(), once);
    EXPECT_EQ(digest(mio::tetrahedralize(s, o).mMesh), digest(r.mMesh));
}

// The octree backend keeps the regions and the CAD faces: a box in a box
// gives two regions (not sharing points, warned), and the capped cylinder's
// faces become Side regions on the fitted boundary.
TEST(Tetrahedralize, TheOctreeBackendKeepsRegionsAndFaces) {
    std::vector<double> x;
    std::vector<std::int64_t> t;
    add_box(x, t, 0.0, 1.0);
    add_box(x, t, 0.3, 0.6, true);
    mio::TetrahedralizeOptions o;
    o.mBackend = mio::TetrahedralizeBackend::Octree;
    o.mMaxSize = 0.1;
    const auto r = mio::tetrahedralize(triangles(x, t), o);
    ASSERT_EQ(r.mNumRegions, 2);
    ASSERT_EQ(r.mWarnings.size(), 1u);
    EXPECT_NE(r.mWarnings[0].find("regions"), std::string::npos);
    double vol[2] = {0.0, 0.0};
    const double* p = r.mMesh.Points().As<double>();
    const auto cells = r.mMesh.Cells(0);
    const mio::NDArray& region = r.mMesh.CellData(mio::kTetrahedralizeRegion, 0);
    for (std::size_t c = 0; c < cells.NumCells(); ++c) {
        std::array<const double*, 4> q{};
        for (std::size_t k = 0; k < 4; ++k)
            q[k] = p + 3 * mio::detail::read_int(cells.Conn(), 4 * c + k);
        const double a[3] = {q[1][0] - q[0][0], q[1][1] - q[0][1], q[1][2] - q[0][2]};
        const double b[3] = {q[2][0] - q[0][0], q[2][1] - q[0][1], q[2][2] - q[0][2]};
        const double d[3] = {q[3][0] - q[0][0], q[3][1] - q[0][1], q[3][2] - q[0][2]};
        const double v = (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
                          a[2] * (b[0] * d[1] - b[1] * d[0])) /
                         6.0;
        EXPECT_GT(v, 0.0) << c;
        vol[mio::detail::read_int(region, c)] += v;
    }
    EXPECT_NEAR(vol[0], 1.0 - 0.027, 0.02);
    EXPECT_NEAR(vol[1], 0.027, 0.005);

    const mio::Mesh s = skin(brep_solids::capped_cylinder(1.0, 2.0), 0.4);
    const auto c = mio::tetrahedralize(s, o);
    std::size_t sides = 0, entries = 0;
    for (std::size_t i = 0; i < c.mMesh.NumRegions(); ++i) {
        const mio::Region& r_r = c.mMesh.Region(i);
        if (r_r.mKind == mio::RegionKind::Side) {
            ++sides;
            entries += r_r.NumEntries();
        }
    }
    EXPECT_EQ(sides, 3u);  // the side and the two caps
    EXPECT_EQ(static_cast<std::int64_t>(entries), c.mNumBoundaryTriangles);
    const mio::NDArray& face = c.mMesh.CellData("cad:face", 1);
    for (std::size_t i = 0; i < c.mMesh.Cells(1).NumCells(); ++i)
        EXPECT_GE(mio::detail::read_int(face, i), 0) << i;
    // max_cells bounds it too.
    o = {};
    o.mBackend = mio::TetrahedralizeBackend::Octree;
    o.mMaxCells = 2;
    EXPECT_THROW((void)mio::tetrahedralize(s, o), std::length_error);
}

// Roadmap 6.7.6: the sphere's skin filled. Every tetrahedron is positive,
// the surface's triangles are facets of exactly one, the volume is the
// surface's own enclosed volume (so within the skin's chordal error of
// 4/3 pi), and two runs agree.
TEST(Tetrahedralize, TheSphereIsFilledKeepingItsSkin) {
    if (!mio::tetrahedralize_has_netgen())
        GTEST_SKIP() << "built without Netgen";
    const mio::Mesh s = skin(brep_solids::sphere(1.0), 0.2);
    mio::TetrahedralizeOptions o;
    o.mMaxSize = 0.2;
    const auto r = mio::tetrahedralize(s, o);
    std::printf(
        "[6.7] sphere: %lld tetrahedra, %lld new points, dihedral %.1f to %.1f, volume "
        "%.6f (4/3 pi %.6f)\n",
        static_cast<long long>(r.mNumTetrahedra), static_cast<long long>(r.mNumNewPoints),
        r.mMinDihedralDeg, r.mMaxDihedralDeg, r.mVolume, 4.0 / 3.0 * kPi);
    EXPECT_TRUE(r.mBoundaryPreserved);
    EXPECT_EQ(r.mHausdorffToSurface, 0.0);
    EXPECT_EQ(r.mNumRegions, 1);
    EXPECT_GT(r.mNumNewPoints, 0);
    EXPECT_GT(r.mMinDihedralDeg, 0.0);
    EXPECT_LT(r.mMaxDihedralDeg, 180.0);
    EXPECT_NEAR(r.mVolume, 4.0 / 3.0 * kPi, 0.05);
    // The input's points first, unchanged.
    const double* p0 = s.Points().As<double>();
    const double* p1 = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < 3 * s.NumPoints(); ++i)
        ASSERT_EQ(p0[i], p1[i]) << i;
    // Each surface triangle a facet of one tetrahedron, every other facet of
    // two.
    const auto f = facets(r.mMesh);
    std::set<std::array<std::int64_t, 3>> surface;
    const auto tri = s.Cells(0);
    for (std::size_t c = 0; c < tri.NumCells(); ++c) {
        std::array<std::int64_t, 3> k{};
        for (std::size_t j = 0; j < 3; ++j)
            k[j] = mio::detail::read_int(tri.Conn(), 3 * c + j);
        std::sort(k.begin(), k.end());
        surface.insert(k);
        ASSERT_TRUE(f.count(k)) << c;
        EXPECT_EQ(f.at(k), 1) << c;
    }
    for (const auto& [k, n] : f)
        if (!surface.count(k))
            EXPECT_EQ(n, 2);
    EXPECT_EQ(digest(mio::tetrahedralize(s, o).mMesh), digest(r.mMesh));
    // A smaller size gives more tetrahedra.
    o.mMaxSize = 0.1;
    EXPECT_GT(mio::tetrahedralize(s, o).mNumTetrahedra, r.mNumTetrahedra);
}

// Roadmap 6.7.4 and 6.7.6: a box inside a box gives two cell regions that
// share the inner box's triangles; facing in, a shell is turned.
TEST(Tetrahedralize, ABoxInABoxGivesTwoRegions) {
    if (!mio::tetrahedralize_has_netgen())
        GTEST_SKIP() << "built without Netgen";
    std::vector<double> x;
    std::vector<std::int64_t> t;
    add_box(x, t, 0.0, 1.0);
    add_box(x, t, 0.3, 0.6, true);  // facing in: turned
    mio::TetrahedralizeOptions o;
    o.mMaxSize = 0.25;
    const auto r = mio::tetrahedralize(triangles(x, t), o);
    ASSERT_EQ(r.mNumRegions, 2);
    EXPECT_NEAR(r.mVolume, 1.0, 1e-12);
    double vol[2] = {0.0, 0.0};
    const double* p = r.mMesh.Points().As<double>();
    const auto cells = r.mMesh.Cells(0);
    const mio::NDArray& region = r.mMesh.CellData(mio::kTetrahedralizeRegion, 0);
    for (std::size_t c = 0; c < cells.NumCells(); ++c) {
        std::array<const double*, 4> q{};
        for (std::size_t k = 0; k < 4; ++k)
            q[k] = p + 3 * mio::detail::read_int(cells.Conn(), 4 * c + k);
        const double a[3] = {q[1][0] - q[0][0], q[1][1] - q[0][1], q[1][2] - q[0][2]};
        const double b[3] = {q[2][0] - q[0][0], q[2][1] - q[0][1], q[2][2] - q[0][2]};
        const double d[3] = {q[3][0] - q[0][0], q[3][1] - q[0][1], q[3][2] - q[0][2]};
        const double v = (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
                          a[2] * (b[0] * d[1] - b[1] * d[0])) /
                         6.0;
        EXPECT_GT(v, 0.0) << c;
        vol[mio::detail::read_int(region, c)] += v;
    }
    EXPECT_NEAR(vol[0], 1.0 - 0.027, 1e-12);
    EXPECT_NEAR(vol[1], 0.027, 1e-12);
    EXPECT_TRUE(r.mMesh.HasRegion("region:0", mio::RegionKind::Cell));
    EXPECT_TRUE(r.mMesh.HasRegion("region:1", mio::RegionKind::Cell));
    // The inner box's triangles are a facet of a tetrahedron on each side.
    const auto f = facets(r.mMesh);
    for (std::size_t c = 12; c < 24; ++c) {
        std::array<std::int64_t, 3> k{t[3 * c], t[3 * c + 1], t[3 * c + 2]};
        std::sort(k.begin(), k.end());
        EXPECT_EQ(f.at(k), 2) << c;
    }
}

// The CAD faces become Side regions on the tetrahedra's facets, one per
// face of the model.
TEST(Tetrahedralize, CadFacesBecomeSideRegions) {
    if (!mio::tetrahedralize_has_netgen())
        GTEST_SKIP() << "built without Netgen";
    const mio::Mesh s = skin(brep_solids::capped_cylinder(1.0, 2.0), 0.4);
    const auto r = mio::tetrahedralize(s);
    std::size_t sides = 0, entries = 0;
    for (std::size_t i = 0; i < r.mMesh.NumRegions(); ++i) {
        const mio::Region& r_r = r.mMesh.Region(i);
        if (r_r.mKind == mio::RegionKind::Side) {
            ++sides;
            entries += r_r.NumEntries();
        }
    }
    EXPECT_EQ(sides, 3u);  // the side and the two caps
    EXPECT_EQ(entries, s.Cells(0).NumCells());
    EXPECT_TRUE(r.mMesh.HasCellData("cad:face"));
}

// Roadmap 6.9.1 with 6.7: the volume behind a boundary layer. The sphere's
// layers grown inward, their cap filled: the cap's triangles are facets of
// the tetrahedra, so the two meshes conform along it.
TEST(Tetrahedralize, TheVolumeBehindABoundaryLayerConforms) {
    if (!mio::tetrahedralize_has_netgen())
        GTEST_SKIP() << "built without Netgen";
    const mio::Mesh s = skin(brep_solids::sphere(1.0), 0.25);
    mio::BoundaryLayerOptions bo;
    bo.mFirstHeight = 0.02;
    bo.mNumLayers = 3;
    bo.mReverse = true;
    const auto layers = mio::boundary_layer(s, bo);
    ASSERT_TRUE(layers.mCapClosed);
    // The cap alone: the triangle cells of the cap region.
    const mio::Mesh& lm = layers.mMesh;
    const std::size_t cap_at = lm.FindRegion(mio::kBoundaryLayerCap, mio::RegionKind::Cell);
    ASSERT_NE(cap_at, mio::Mesh::npos);
    const mio::Region& cap = lm.Region(cap_at);
    std::vector<std::int64_t> global(cap.NumEntries());
    for (std::size_t i = 0; i < global.size(); ++i)
        global[i] = mio::detail::read_int(cap.mEntries, i);
    std::vector<std::int64_t> conn;
    std::size_t first = 0;
    for (std::size_t b = 0; b < lm.NumCellBlocks(); ++b) {
        const auto cells = lm.Cells(b);
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            const auto g = static_cast<std::int64_t>(first + c);
            if (cells.Type() == "triangle" && std::binary_search(global.begin(), global.end(), g))
                for (std::size_t k = 0; k < 3; ++k)
                    conn.push_back(mio::detail::read_int(cells.Conn(), 3 * c + k));
        }
        first += cells.NumCells();
    }
    ASSERT_FALSE(conn.empty());
    std::vector<double> x(lm.Points().As<double>(), lm.Points().As<double>() + 3 * lm.NumPoints());
    // The cap faces away from the wall, into the sphere's middle: turned
    // outward by the operation, as every shell is.
    const auto r = mio::tetrahedralize(triangles(x, conn));
    EXPECT_TRUE(r.mBoundaryPreserved);
    const auto f = facets(r.mMesh);
    for (std::size_t c = 0; c < conn.size() / 3; ++c) {
        std::array<std::int64_t, 3> k{conn[3 * c], conn[3 * c + 1], conn[3 * c + 2]};
        std::sort(k.begin(), k.end());
        EXPECT_EQ(f.count(k), 1u) << c;
    }
}

TEST(Tetrahedralize, RefusalsAreByName) {
    std::vector<double> x;
    std::vector<std::int64_t> t;
    add_box(x, t, 0.0, 1.0);
    // Open: a triangle left out.
    std::vector<std::int64_t> open(t.begin(), t.end() - 3);
    EXPECT_THROW((void)mio::tetrahedralize(triangles(x, open)), std::invalid_argument);
    mio::TetrahedralizeOptions o;
    o.mGrading = 2.0;
    EXPECT_THROW((void)mio::tetrahedralize(triangles(x, t), o), std::invalid_argument);
    o = {};
    o.mMaxSize = -1.0;
    EXPECT_THROW((void)mio::tetrahedralize(triangles(x, t), o), std::invalid_argument);
    o = {};
    o.mBackend = static_cast<mio::TetrahedralizeBackend>(9);
    EXPECT_THROW((void)mio::tetrahedralize(triangles(x, t), o), std::invalid_argument);
    // A quad block.
    mio::Mesh q;
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {4, 3});
    const double c4[12] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    std::copy(c4, c4 + 12, p.As<double>());
    q.AssignPoints(std::move(p));
    mio::NDArray c = mio::NDArray::Uninit(mio::DType::Int64, {1, 4});
    for (std::int64_t i = 0; i < 4; ++i)
        c.As<std::int64_t>()[i] = i;
    q.AddCellBlock("quad", std::move(c));
    EXPECT_THROW((void)mio::tetrahedralize(q), std::invalid_argument);
    if (mio::tetrahedralize_has_netgen()) {
        o = {};
        o.mMaxCells = 2;
        o.mMaxSize = 0.1;
        EXPECT_THROW((void)mio::tetrahedralize(triangles(x, t), o), std::length_error);
    }
}
