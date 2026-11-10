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

// System includes
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/adapt_volume.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/primitives.hpp"

// Project includes (private, not installed)
#include "../../src/cpp/src/detail/fem_kernels.hpp"
#include "../../src/cpp/src/operations/adapt_volume_internal.hpp"

namespace {

namespace mio = meshioplusplus;
using meshioplusplus::adapt_volume;
using meshioplusplus::AdaptVolumeOptions;
using meshioplusplus::AdaptVolumeResult;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::RegionKind;

Mesh tet_box(std::int64_t N, double Size = 1.0) {
    return mio::box({{N, N, N}}, {{Size, Size, Size}}, {{0.0, 0.0, 0.0}}, "tetra");
}

// A cylinder of radius 1 and height 1: a Kuhn-split box on [-1, 1]^2 x [0, 1]
// pushed onto the disk by the squircle map (its square's sides become arcs of
// the unit circle), so the tetrahedra conform and the rims are smooth loops.
Mesh tet_cylinder(std::int64_t N = 6) {
    Mesh b = mio::box({{N, N, 2}}, {{2.0, 2.0, 1.0}}, {{-1.0, -1.0, 0.0}}, "tetra");
    const NDArray& p = b.Points();
    NDArray q = NDArray::Uninit(mio::DType::Float64, {b.NumPoints(), 3});
    for (std::size_t i = 0; i < b.NumPoints(); ++i) {
        const double x = mio::detail::read_double(p, i * 3),
                     y = mio::detail::read_double(p, i * 3 + 1);
        mio::detail::write_double(q, i * 3, x * std::sqrt(1.0 - 0.5 * y * y));
        mio::detail::write_double(q, i * 3 + 1, y * std::sqrt(1.0 - 0.5 * x * x));
        mio::detail::write_double(q, i * 3 + 2, mio::detail::read_double(p, i * 3 + 2));
    }
    Mesh m;
    m.AssignPoints(std::move(q));
    NDArray conn = b.Cells(0).Conn();
    m.AddCellBlock("tetra", std::move(conn));
    return m;
}

AdaptVolumeOptions passes_off() {
    AdaptVolumeOptions o;
    o.mEnableCollapse = false;
    o.mEnableSplit = false;
    o.mEnableSwap = false;
    o.mEnableSmooth = false;
    return o;
}

std::vector<std::int64_t> entries(const mio::Region& rR) {
    return std::vector<std::int64_t>(rR.Entries(), rR.Entries() + rR.mEntries.Size());
}

TEST(AdaptVolume, EveryPassOffIsTheIdentity) {
    const Mesh in = tet_box(3);
    const AdaptVolumeResult r = adapt_volume(in, passes_off());
    const Mesh& out = r.mMesh;
    ASSERT_EQ(out.NumPoints(), in.NumPoints());
    ASSERT_EQ(out.NumCellBlocks(), in.NumCellBlocks());
    EXPECT_EQ(r.mNumTetsBefore, r.mNumTetsAfter);
    EXPECT_EQ(r.mNumPointsBefore, r.mNumPointsAfter);
    for (std::size_t i = 0; i < in.NumPoints() * 3; ++i)
        EXPECT_EQ(mio::detail::read_double(out.Points(), i),
                  mio::detail::read_double(in.Points(), i));
    // same tetrahedra (as vertex sets) in the same order
    const NDArray& ci = in.Cells(0).Conn();
    const NDArray& co = out.Cells(0).Conn();
    ASSERT_EQ(ci.Size(), co.Size());
    for (std::size_t c = 0; c < ci.Size() / 4; ++c) {
        std::array<std::int64_t, 4> a{}, b{};
        for (std::size_t k = 0; k < 4; ++k) {
            a[k] = mio::detail::read_int(ci, c * 4 + k);
            b[k] = mio::detail::read_int(co, c * 4 + k);
        }
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        EXPECT_EQ(a, b);
    }
    // every side region comes back with the same faces
    ASSERT_EQ(out.NumRegions(), in.NumRegions());
    for (std::size_t i = 0; i < in.NumRegions(); ++i) {
        const mio::Region& ri = in.Region(i);
        ASSERT_TRUE(out.HasRegion(ri.mName, ri.mKind)) << ri.mName;
        const mio::Region& ro = out.Region(out.FindRegion(ri.mName, ri.mKind));
        if (ri.mKind != RegionKind::Side) {
            EXPECT_EQ(entries(ri), entries(ro)) << ri.mName;
            continue;
        }
        // same faces as vertex sets (the facet index follows the stored orientation)
        auto faces = [](const Mesh& rM, const mio::Region& rR) {
            std::vector<std::array<std::int64_t, 3>> f;
            constexpr int kFacet[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};
            const NDArray& c = rM.Cells(0).Conn();
            for (std::size_t e = 0; e < rR.NumEntries(); ++e) {
                const std::size_t cell = static_cast<std::size_t>(rR.Entries()[2 * e]);
                const int facet = static_cast<int>(rR.Entries()[2 * e + 1]);
                std::array<std::int64_t, 3> k{};
                for (int q = 0; q < 3; ++q)
                    k[static_cast<std::size_t>(q)] = mio::detail::read_int(
                        c, cell * 4 + static_cast<std::size_t>(kFacet[facet][q]));
                std::sort(k.begin(), k.end());
                f.push_back(k);
            }
            std::sort(f.begin(), f.end());
            return f;
        };
        EXPECT_EQ(faces(in, ri), faces(out, ro)) << ri.mName;
    }
    EXPECT_DOUBLE_EQ(r.mVolumeBefore, r.mVolumeAfter);
    EXPECT_DOUBLE_EQ(r.mMinDihedralBefore, r.mMinDihedralAfter);
}

TEST(AdaptVolume, ClassifiesTheCornersAndEdgesOfACube) {
    const AdaptVolumeResult r = adapt_volume(tet_box(4), passes_off());
    EXPECT_EQ(r.mNumCorners, 8);
    EXPECT_EQ(r.mNumChains, 12);
    EXPECT_EQ(r.mNumRidgeVertices, 12 * 3);
    EXPECT_EQ(r.mNumPatches, 6);
}

TEST(AdaptVolume, ClassifiesTheRimsOfACylinder) {
    const AdaptVolumeResult r = adapt_volume(tet_cylinder(), passes_off());
    // each rim is a closed loop with no corner: cut in two, so two corners each
    EXPECT_EQ(r.mNumChains, 4);
    EXPECT_EQ(r.mNumCorners, 4);
    EXPECT_EQ(r.mNumPatches, 3);
    EXPECT_EQ(r.mNumRidgeVertices, 2 * 4 * 6 - 4);  // 24 rim points per rim
}

TEST(AdaptVolume, SizeBoundsAndHausdorffDefaults) {
    const Mesh in = tet_box(4, 2.0);
    const AdaptVolumeResult r = adapt_volume(in, passes_off());
    const double diag = 2.0 * std::sqrt(3.0);
    EXPECT_NEAR(r.mHmax, diag / 8.0, 1e-12);
    EXPECT_NEAR(r.mHausd, 1e-3 * diag, 1e-12);
    EXPECT_GT(r.mHmin, 0.0);
    EXPECT_LE(r.mHmin, r.mHmax);
    AdaptVolumeOptions o = passes_off();
    o.mHmin = 0.1;
    o.mHmax = 0.4;
    o.mHausd = 0.01;
    o.mHausdRelative = true;
    const AdaptVolumeResult r2 = adapt_volume(in, o);
    EXPECT_DOUBLE_EQ(r2.mHmin, 0.1);
    EXPECT_DOUBLE_EQ(r2.mHmax, 0.4);
    EXPECT_NEAR(r2.mHausd, 0.01 * diag, 1e-12);
}

TEST(AdaptVolume, KeepsMaterialArraysAndDropsTheRest) {
    Mesh in = tet_box(2);
    const std::size_t nc = in.Cells(0).NumCells();
    NDArray mat = NDArray::Uninit(mio::DType::Int32, {nc});
    NDArray other = NDArray::Uninit(mio::DType::Float64, {nc});
    for (std::size_t c = 0; c < nc; ++c) {
        mio::detail::write_int(mat, c, c < nc / 2 ? 1 : 2);
        mio::detail::write_double(other, c, 0.5);
    }
    in.AddCellData("octree_mesh:material", {mat});
    in.AddCellData("quality", {other});
    const AdaptVolumeResult r = adapt_volume(in, passes_off());
    ASSERT_TRUE(r.mMesh.HasCellData("octree_mesh:material"));
    EXPECT_FALSE(r.mMesh.HasCellData("quality"));
    const NDArray& m = r.mMesh.CellData("octree_mesh:material", 0);
    for (std::size_t c = 0; c < nc; ++c)
        EXPECT_EQ(mio::detail::read_int(m, c), mio::detail::read_int(mat, c));
    EXPECT_FALSE(r.mWarnings.empty());
}

TEST(AdaptVolume, RejectsBadInputByName) {
    EXPECT_THROW(adapt_volume(mt::tri_mesh()), std::invalid_argument);
    AdaptVolumeOptions o;
    o.mHgrad = 0.9;
    EXPECT_THROW(adapt_volume(tet_box(1), o), std::invalid_argument);
    o = AdaptVolumeOptions{};
    o.mSizeArray = "nope";
    EXPECT_THROW(adapt_volume(tet_box(1), o), std::invalid_argument);
    o = AdaptVolumeOptions{};
    o.mHmin = 1.0;
    o.mHmax = 0.5;
    EXPECT_THROW(adapt_volume(tet_box(1), o), std::invalid_argument);
    EXPECT_EQ(mio::adapt_volume_quality_from_name("amips"), mio::AdaptVolumeQuality::Amips);
    EXPECT_STREQ(mio::adapt_volume_quality_name(mio::AdaptVolumeQuality::MeanRatio), "mean_ratio");
    EXPECT_THROW(mio::adapt_volume_quality_from_name("x"), std::invalid_argument);
}

// --- M2: the collapse guards -----------------------------------------------------

// The 20 tetrahedra from the centre (point 12) of an icosahedron to its faces.
Mesh icosa_star() {
    const Mesh ico = mio::icosa_sphere(1.0, 0);
    std::vector<std::vector<double>> pts;
    for (std::size_t i = 0; i < ico.NumPoints(); ++i)
        pts.push_back({mio::detail::read_double(ico.Points(), i * 3),
                       mio::detail::read_double(ico.Points(), i * 3 + 1),
                       mio::detail::read_double(ico.Points(), i * 3 + 2)});
    pts.push_back({0.0, 0.0, 0.0});
    std::vector<std::vector<std::int64_t>> tets;
    const NDArray& c = ico.Cells(0).Conn();
    for (std::size_t f = 0; f < ico.Cells(0).NumCells(); ++f)
        tets.push_back({mio::detail::read_int(c, f * 3), mio::detail::read_int(c, f * 3 + 2),
                        mio::detail::read_int(c, f * 3 + 1), 12});
    return mt::make_mesh(pts, "tetra", tets);
}

Mesh with_sizes(Mesh m, double H) {
    NDArray h = NDArray::Uninit(mio::DType::Float64, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        mio::detail::write_double(h, i, H);
    m.AddPointData("h", std::move(h));
    return m;
}

std::int64_t find_point(const Mesh& rM, double X, double Y, double Z) {
    for (std::size_t i = 0; i < rM.NumPoints(); ++i)
        if (std::abs(mio::detail::read_double(rM.Points(), i * 3) - X) < 1e-12 &&
            std::abs(mio::detail::read_double(rM.Points(), i * 3 + 1) - Y) < 1e-12 &&
            std::abs(mio::detail::read_double(rM.Points(), i * 3 + 2) - Z) < 1e-12)
            return static_cast<std::int64_t>(i);
    return -1;
}

TEST(AdaptVolume, CollapseGuardsOnAnIcosahedralStar) {
    // hmin 10: every size is clamped up, so no edge is metric-long
    AdaptVolumeOptions o;
    o.mSizeArray = "h";
    o.mHmin = 10.0;
    o.mHmax = 100.0;
    // a convex star: the centre collapses onto any corner, every cone stays positive
    EXPECT_EQ(mio::detail::adapt_volume_probe_collapse(with_sizes(icosa_star(), 10.0), o, 12, 0),
              "none");
    // a boundary point may not leave the boundary (M2: the boundary stays)
    EXPECT_EQ(mio::detail::adapt_volume_probe_collapse(with_sizes(icosa_star(), 10.0), o, 0, 12),
              "class");
    // the quality floor: a cone over 15 faces is worse than the regular-ish star
    AdaptVolumeOptions q = o;
    q.mQualityMin = 89.0;
    EXPECT_EQ(mio::detail::adapt_volume_probe_collapse(with_sizes(icosa_star(), 10.0), q, 12, 0),
              "quality");
    // tiny sizes: the edges the collapse would create are far above 4/3
    AdaptVolumeOptions l = o;
    l.mHmin = 1e-3;
    EXPECT_EQ(mio::detail::adapt_volume_probe_collapse(with_sizes(icosa_star(), 0.01), l, 12, 0),
              "length");
    // frozen
    AdaptVolumeOptions f = o;
    f.mFrozen.assign(13, 0);
    f.mFrozen[12] = 1;
    EXPECT_EQ(mio::detail::adapt_volume_probe_collapse(with_sizes(icosa_star(), 10.0), f, 12, 0),
              "class");
}

TEST(AdaptVolume, CollapseIntoADentIsAnInversion) {
    // the icosahedral star with corner 0 pushed below its neighbours: its link is not
    // convex, and moving the centre onto some neighbours of the dent flips a cone
    Mesh m = icosa_star();
    NDArray p = m.Points();
    p.MakeOwned();
    for (std::size_t d = 0; d < 3; ++d)
        mio::detail::write_double(p, d, 0.15 * mio::detail::read_double(p, d));
    Mesh dented;
    dented.AssignPoints(std::move(p));
    dented.AddCellBlock("tetra", NDArray(m.Cells(0).Conn()));
    AdaptVolumeOptions o;
    o.mHmin = 10.0;
    o.mHmax = 100.0;
    int inversions = 0;
    for (std::int64_t k = 0; k < 12; ++k) {
        const std::string why = mio::detail::adapt_volume_probe_collapse(dented, o, 12, k);
        EXPECT_TRUE(why == "none" || why == "inversion" || why == "quality") << k << ": " << why;
        inversions += why == "inversion";
        // on the convex star every one is valid
        EXPECT_EQ(mio::detail::adapt_volume_probe_collapse(m, o, 12, k), "none") << k;
    }
    EXPECT_GT(inversions, 0);
}

// The sorted boundary triangles of a mesh, as coordinate triples.
std::vector<std::array<std::array<double, 3>, 3>> boundary_triangles(const Mesh& rM) {
    std::map<std::array<std::int64_t, 3>, int> count;
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
        const NDArray& c = rM.Cells(b).Conn();
        for (std::size_t t = 0; t < rM.Cells(b).NumCells(); ++t) {
            std::array<std::int64_t, 4> v{};
            for (std::size_t k = 0; k < 4; ++k)
                v[k] = mio::detail::read_int(c, t * 4 + k);
            constexpr int kF[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};
            for (const auto& f : kF) {
                std::array<std::int64_t, 3> k{v[f[0]], v[f[1]], v[f[2]]};
                std::sort(k.begin(), k.end());
                ++count[k];
            }
        }
    }
    std::vector<std::array<std::array<double, 3>, 3>> out;
    for (const auto& [k, n] : count) {
        if (n != 1)
            continue;
        std::array<std::array<double, 3>, 3> tri{};
        for (int i = 0; i < 3; ++i)
            for (int d = 0; d < 3; ++d)
                tri[i][d] =
                    mio::detail::read_double(rM.Points(), static_cast<std::size_t>(k[i]) * 3 + d);
        std::sort(tri.begin(), tri.end());
        out.push_back(tri);
    }
    std::sort(out.begin(), out.end());
    return out;
}

double min_dihedral(const Mesh& rM) {
    double m = 180.0;
    for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
        const NDArray& c = rM.Cells(b).Conn();
        for (std::size_t t = 0; t < rM.Cells(b).NumCells(); ++t) {
            mio::detail::Vec3 p[4];
            for (std::size_t k = 0; k < 4; ++k) {
                const std::size_t v = static_cast<std::size_t>(mio::detail::read_int(c, t * 4 + k));
                p[k] = {mio::detail::read_double(rM.Points(), v * 3),
                        mio::detail::read_double(rM.Points(), v * 3 + 1),
                        mio::detail::read_double(rM.Points(), v * 3 + 2)};
            }
            EXPECT_GT(mio::detail::tet_volume6(p[0], p[1], p[2], p[3]), 0.0);
            m = std::min(m, mio::detail::tet_min_dihedral(p[0], p[1], p[2], p[3]));
        }
    }
    return m;
}

struct CheckInvariants {
    CheckInvariants() { setenv("MESHIOPLUSPLUS_ADAPT_VOLUME_CHECK", "1", 1); }
    ~CheckInvariants() { unsetenv("MESHIOPLUSPLUS_ADAPT_VOLUME_CHECK"); }
};

TEST(AdaptVolume, InteriorCollapseCoarsensAFrozenBoundaryBox) {
    const CheckInvariants check;
    const Mesh in = tet_box(8);
    AdaptVolumeOptions o;
    o.mPreserveBoundary = true;
    o.mHmax = 0.5;
    o.mEnableSplit = false;
    o.mEnableSwap = false;
    o.mEnableSmooth = false;
    const AdaptVolumeResult r = adapt_volume(in, o);
    EXPECT_GT(r.mNumCollapsed, 0);
    EXPECT_LT(r.mNumTetsAfter, r.mNumTetsBefore);
    EXPECT_EQ(boundary_triangles(r.mMesh), boundary_triangles(in));
    EXPECT_NEAR(r.mVolumeAfter, r.mVolumeBefore, 1e-12);
    EXPECT_GE(min_dihedral(r.mMesh), std::min(r.mMinDihedralBefore, 10.0) - 1e-9);
    EXPECT_NEAR(r.mMinDihedralAfter, min_dihedral(r.mMesh), 1e-9);
}

TEST(AdaptVolume, CollapseIsDeterministic) {
    const Mesh in = tet_box(6);
    AdaptVolumeOptions o;
    o.mPreserveBoundary = true;
    o.mHmax = 0.5;
    auto digest = [](const AdaptVolumeResult& rR) {
        meshioplusplus::bench::MeshDigest d;
        d.Of(rR.mMesh);
        d.U64(static_cast<std::uint64_t>(rR.mNumCollapsed));
        return d.Value();
    };
    const std::uint64_t first = digest(adapt_volume(in, o));
    for (int run = 0; run < 3; ++run)
        EXPECT_EQ(digest(adapt_volume(in, o)), first);
}

// The unit sphere through octree_mesh's default graded tetrahedra.
Mesh octree_sphere(double MinSize = 0.0) {
    mio::OctreeMeshOptions o;
    if (MinSize > 0.0)
        o.mMinSize = MinSize;
    return mio::octree_mesh(mio::icosa_sphere(1.0, 4), o).mMesh;
}

void print_report(const char* pName, const AdaptVolumeResult& rR) {
    std::printf(
        "[ adapt    ] %s: tets %lld -> %lld, points %lld -> %lld, dihedral [%.2f, %.2f] -> "
        "[%.2f, %.2f], collapsed %lld split %lld swapped %lld smoothed %lld, iterations %d, "
        "rejected link %lld inv %lld qual %lld len %lld class %lld hausd %lld normal %lld\n",
        pName, static_cast<long long>(rR.mNumTetsBefore), static_cast<long long>(rR.mNumTetsAfter),
        static_cast<long long>(rR.mNumPointsBefore), static_cast<long long>(rR.mNumPointsAfter),
        rR.mMinDihedralBefore, rR.mMaxDihedralBefore, rR.mMinDihedralAfter, rR.mMaxDihedralAfter,
        static_cast<long long>(rR.mNumCollapsed), static_cast<long long>(rR.mNumSplit),
        static_cast<long long>(rR.mNumSwapped), static_cast<long long>(rR.mNumSmoothed),
        rR.mIterations, static_cast<long long>(rR.mNumRejectedByLink),
        static_cast<long long>(rR.mNumRejectedByInversion),
        static_cast<long long>(rR.mNumRejectedByQuality),
        static_cast<long long>(rR.mNumRejectedByLength),
        static_cast<long long>(rR.mNumRejectedByClass),
        static_cast<long long>(rR.mNumRejectedByHausdorff),
        static_cast<long long>(rR.mNumRejectedByNormal));
}

TEST(AdaptVolume, OctreeSphereWithAFrozenBoundary) {
    const CheckInvariants check;
    const Mesh in = octree_sphere();
    AdaptVolumeOptions o;
    o.mPreserveBoundary = true;
    const AdaptVolumeResult r = adapt_volume(in, o);
    print_report("octree sphere, boundary frozen", r);
    EXPECT_LT(r.mNumTetsAfter, r.mNumTetsBefore);
    EXPECT_EQ(boundary_triangles(r.mMesh), boundary_triangles(in));
    EXPECT_NEAR(r.mVolumeAfter, r.mVolumeBefore, 1e-12 * r.mVolumeBefore);
    EXPECT_GE(min_dihedral(r.mMesh), std::min(r.mMinDihedralBefore, 10.0) - 1e-9);
}

// --- M3: splits ---------------------------------------------------------------------

TEST(AdaptVolume, SplitRefinesToTheSizeAndKeepsTheBoxAndItsSides) {
    const CheckInvariants check;
    const Mesh in = tet_box(2);
    AdaptVolumeOptions o;
    o.mHmin = 0.1;
    o.mHmax = 0.15;  // well below the input's 0.5 edges
    o.mSizeArray = "h";
    o.mEnableCollapse = false;
    o.mEnableSwap = false;
    o.mEnableSmooth = false;
    const AdaptVolumeResult r = adapt_volume(with_sizes(in, 0.12), o);
    print_report("box split", r);
    EXPECT_GT(r.mNumSplit, 0);
    EXPECT_GT(r.mNumTetsAfter, 4 * r.mNumTetsBefore);
    EXPECT_NEAR(r.mVolumeAfter, 1.0, 1e-12);
    // every face of a side region still lies on its side of the box
    const Mesh& m = r.mMesh;
    constexpr int kF[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};
    const std::map<std::string, std::pair<int, double>> planes{
        {"xmin", {0, 0.0}}, {"xmax", {0, 1.0}}, {"ymin", {1, 0.0}},
        {"ymax", {1, 1.0}}, {"zmin", {2, 0.0}}, {"zmax", {2, 1.0}}};
    for (const auto& [name, plane] : planes) {
        ASSERT_TRUE(m.HasRegion(name, RegionKind::Side)) << name;
        const mio::Region& reg = m.Region(m.FindRegion(name, RegionKind::Side));
        EXPECT_GT(reg.NumEntries(), 8u) << name;
        double area = 0.0;
        for (std::size_t e = 0; e < reg.NumEntries(); ++e) {
            const std::size_t cell = static_cast<std::size_t>(reg.Entries()[2 * e]);
            const int facet = static_cast<int>(reg.Entries()[2 * e + 1]);
            mio::detail::Vec3 p[3];
            for (int q = 0; q < 3; ++q) {
                const std::size_t v = static_cast<std::size_t>(mio::detail::read_int(
                    m.Cells(0).Conn(), cell * 4 + static_cast<std::size_t>(kF[facet][q])));
                for (int d = 0; d < 3; ++d)
                    p[q][static_cast<std::size_t>(d)] =
                        mio::detail::read_double(m.Points(), v * 3 + static_cast<std::size_t>(d));
                EXPECT_NEAR(p[q][static_cast<std::size_t>(plane.first)], plane.second, 1e-12)
                    << name;
            }
            area +=
                0.5 * mio::detail::vec3_norm(mio::detail::vec3_cross(
                          mio::detail::vec3_sub(p[1], p[0]), mio::detail::vec3_sub(p[2], p[0])));
        }
        EXPECT_NEAR(area, 1.0, 1e-12) << name;
    }
    EXPECT_GE(min_dihedral(m), std::min(r.mMinDihedralBefore, 10.0) - 1e-9);
}

TEST(AdaptVolume, CollapseAndSplitDoNotOscillate) {
    const CheckInvariants check;
    // a uniform size between the box's two edge lengths: both passes have work
    AdaptVolumeOptions o;
    o.mHmin = 0.18;
    o.mHmax = 0.18;
    o.mEnableSwap = false;
    o.mEnableSmooth = false;
    const AdaptVolumeResult r = adapt_volume(tet_box(4), o);
    print_report("box collapse+split", r);
    const AdaptVolumeResult again = adapt_volume(r.mMesh, o);
    print_report("box collapse+split, again", again);
    EXPECT_LE(again.mNumCollapsed + again.mNumSplit, r.mNumTetsAfter / 10 + 1);
}

}  // namespace
