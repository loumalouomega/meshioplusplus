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
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/adapt_volume.hpp"
#include "meshioplusplus/operations/merge.hpp"
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
Mesh tet_cylinder(std::int64_t N = 6, std::int64_t Nz = 2) {
    Mesh b = mio::box({{N, N, Nz}}, {{2.0, 2.0, 1.0}}, {{-1.0, -1.0, 0.0}}, "tetra");
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

Mesh with_sizes(Mesh&& rM, double H) {
    Mesh m = std::move(rM);
    NDArray h = NDArray::Uninit(mio::DType::Float64, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        mio::detail::write_double(h, i, H);
    m.AddPointData("h", std::move(h));
    return m;
}

// A KRATOS mesh is not copy-constructible: an lvalue goes through the uniform clone.
Mesh with_sizes(const Mesh& rM, double H) {
    return with_sizes(mio::detail::clone_mesh(rM), H);
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

std::size_t count_below(const Mesh& rM, double Deg) {
    std::size_t n = 0;
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
            n += mio::detail::tet_min_dihedral(p[0], p[1], p[2], p[3]) < Deg;
        }
    }
    return n;
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
    // the default sphere's interior is near the size field already: the
    // passes trade a few vertices at the 4/5 and 4/3 thresholds, and the least
    // angle rises
    EXPECT_LE(r.mNumTetsAfter, r.mNumTetsBefore + r.mNumTetsBefore / 50);
    EXPECT_GE(r.mMinDihedralAfter, r.mMinDihedralBefore);
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

// --- M4: swaps ----------------------------------------------------------------------

TEST(AdaptVolume, SwapsRaiseTheLeastAngleOfSplitTetrahedra) {
    const CheckInvariants check;
    AdaptVolumeOptions o;
    o.mHmin = 0.1;
    o.mHmax = 0.15;
    o.mSizeArray = "h";
    o.mEnableCollapse = false;
    o.mEnableSmooth = false;
    o.mEnableSwap = false;
    const AdaptVolumeResult split = adapt_volume(with_sizes(tet_box(2), 0.12), o);
    // swaps alone on the split mesh: the boundary stays, the least angle rises
    AdaptVolumeOptions w;
    w.mEnableCollapse = false;
    w.mEnableSplit = false;
    w.mEnableSmooth = false;
    const AdaptVolumeResult swapped = adapt_volume(split.mMesh, w);
    print_report("box split, then swaps", swapped);
    EXPECT_GT(swapped.mNumSwapped, 0);
    // the worst tetrahedra here touch the boundary, which no swap changes, so
    // the least angle holds while the number of bad ones drops
    EXPECT_GE(swapped.mMinDihedralAfter, split.mMinDihedralAfter);
    EXPECT_LT(count_below(swapped.mMesh, 20.0), count_below(split.mMesh, 20.0) / 2);
    EXPECT_EQ(swapped.mNumPointsAfter, split.mNumPointsAfter);
    EXPECT_NEAR(swapped.mVolumeAfter, 1.0, 1e-12);
    EXPECT_EQ(boundary_triangles(swapped.mMesh), boundary_triangles(split.mMesh));
}

TEST(AdaptVolume, SwapsNeverCrossAMaterialInterface) {
    const CheckInvariants check;
    Mesh in = tet_box(4);
    const std::size_t nc = in.Cells(0).NumCells();
    NDArray mat = NDArray::Uninit(mio::DType::Int64, {nc});
    for (std::size_t c = 0; c < nc; ++c) {
        double x = 0.0;
        for (std::size_t k = 0; k < 4; ++k)
            x += mio::detail::read_double(
                in.Points(),
                static_cast<std::size_t>(mio::detail::read_int(in.Cells(0).Conn(), c * 4 + k)) * 3);
        mio::detail::write_int(mat, c, x < 2.0 ? 1 : 2);  // centroid x < 0.5
    }
    in.AddCellData("octree_mesh:material", {mat});
    AdaptVolumeOptions o;
    o.mHmin = 0.08;
    o.mHmax = 0.12;
    o.mEnableCollapse = false;
    o.mEnableSmooth = false;
    const AdaptVolumeResult r = adapt_volume(in, o);
    print_report("two materials split+swap", r);
    EXPECT_GT(r.mNumSwapped, 0);
    // every tetrahedron of material 1 lies left of x = 0.5, of material 2 right of it
    const NDArray& m = r.mMesh.CellData("octree_mesh:material", 0);
    const NDArray& c = r.mMesh.Cells(0).Conn();
    for (std::size_t t = 0; t < r.mMesh.Cells(0).NumCells(); ++t) {
        const std::int64_t mt = mio::detail::read_int(m, t);
        for (std::size_t k = 0; k < 4; ++k) {
            const double x = mio::detail::read_double(
                r.mMesh.Points(),
                static_cast<std::size_t>(mio::detail::read_int(c, t * 4 + k)) * 3);
            if (mt == 1)
                EXPECT_LE(x, 0.5 + 1e-12);
            else
                EXPECT_GE(x, 0.5 - 1e-12);
        }
    }
    ASSERT_EQ(r.mRegionVolumeChange.size(), 2u);
    EXPECT_NEAR(r.mRegionVolumeChange[0].second, 0.0, 1e-12);
    EXPECT_NEAR(r.mRegionVolumeChange[1].second, 0.0, 1e-12);
}

// --- M5: smoothing ------------------------------------------------------------------

TEST(AdaptVolume, SmoothingRaisesTheInteriorWithoutTouchingTheBoundary) {
    const CheckInvariants check;
    // a split box (bad interior tetrahedra), then smoothing alone
    AdaptVolumeOptions o;
    o.mHmin = 0.1;
    o.mHmax = 0.15;
    o.mSizeArray = "h";
    o.mEnableCollapse = false;
    o.mEnableSmooth = false;
    o.mEnableSwap = false;
    const Mesh rough = adapt_volume(with_sizes(tet_box(2), 0.12), o).mMesh;
    for (const auto metric :
         {mio::AdaptVolumeQuality::MinDihedral, mio::AdaptVolumeQuality::Amips}) {
        AdaptVolumeOptions s;
        s.mPreserveBoundary = true;
        s.mQualityMetric = metric;
        s.mEnableCollapse = false;
        s.mEnableSplit = false;
        s.mEnableSwap = false;
        const AdaptVolumeResult r = adapt_volume(rough, s);
        print_report(metric == mio::AdaptVolumeQuality::Amips ? "smoothing alone, amips"
                                                              : "smoothing alone, odt",
                     r);
        EXPECT_GT(r.mNumSmoothed, 0);
        EXPECT_GE(r.mMinDihedralAfter, r.mMinDihedralBefore);
        EXPECT_GT(r.mMeanQualityAfter, r.mMeanQualityBefore);
        EXPECT_EQ(boundary_triangles(r.mMesh), boundary_triangles(rough));
        EXPECT_NEAR(r.mVolumeAfter, r.mVolumeBefore, 1e-12);
    }
}

// --- M6: the boundary ---------------------------------------------------------------

TEST(AdaptVolume, BoundarySmoothingSlidesWithinTheSidesOfABox) {
    const CheckInvariants check;
    AdaptVolumeOptions o;
    o.mHmin = 0.1;
    o.mHmax = 0.15;
    o.mSizeArray = "h";
    o.mEnableCollapse = false;
    o.mEnableSmooth = false;
    o.mEnableSwap = false;
    const Mesh rough = adapt_volume(with_sizes(tet_box(2), 0.12), o).mMesh;
    AdaptVolumeOptions s;
    s.mEnableCollapse = false;
    s.mEnableSplit = false;
    s.mEnableSwap = false;
    const AdaptVolumeResult r = adapt_volume(rough, s);
    print_report("box smoothing with the boundary", r);
    EXPECT_NE(boundary_triangles(r.mMesh), boundary_triangles(rough));  // sides smoothed
    EXPECT_NEAR(r.mVolumeAfter, 1.0, 1e-12);  // within the planes of the box
    EXPECT_LE(r.mMaxBoundaryDeviation, 1e-12);
    EXPECT_GE(r.mMinDihedralAfter, r.mMinDihedralBefore);
}

TEST(AdaptVolume, CoarsensTheBoundaryOfAnOverRefinedBoxKeepingItsCornersAndEdges) {
    const CheckInvariants check;
    const Mesh in = tet_box(8);
    AdaptVolumeOptions o;
    o.mHmin = 0.3;
    o.mHmax = 0.4;
    const AdaptVolumeResult r = adapt_volume(in, o);
    print_report("box 8^3, h 0.3", r);
    EXPECT_LT(r.mNumTetsAfter, r.mNumTetsBefore / 4);
    EXPECT_NEAR(r.mVolumeAfter, 1.0, 1e-12);
    EXPECT_LE(r.mMaxBoundaryDeviation, 1e-12);
    // the eight corners are still points of the mesh
    for (const double x : {0.0, 1.0})
        for (const double y : {0.0, 1.0})
            for (const double z : {0.0, 1.0})
                EXPECT_GE(find_point(r.mMesh, x, y, z), 0);
    EXPECT_EQ(r.mNumCorners, 8);
    EXPECT_GE(r.mMinDihedralAfter, std::min(r.mMinDihedralBefore, 10.0) - 1e-9);
}

TEST(AdaptVolume, SphereStaysWithinHausd) {
    const CheckInvariants check;
    const Mesh in = octree_sphere(0.04);
    AdaptVolumeOptions o;
    const AdaptVolumeResult r = adapt_volume(in, o);
    print_report("octree sphere min_size 0.04", r);
    std::printf(
        "[ adapt    ] hmin %.3g hmax %.3g hausd %.3g, boundary deviation %.3g, volume %.6f -> "
        "%.6f\n",
        r.mHmin, r.mHmax, r.mHausd, r.mMaxBoundaryDeviation, r.mVolumeBefore, r.mVolumeAfter);
    EXPECT_LE(r.mMaxBoundaryDeviation, r.mHausd);
    EXPECT_NEAR(r.mVolumeAfter, r.mVolumeBefore, 0.01 * r.mVolumeBefore);
    EXPECT_GE(r.mMinDihedralAfter, std::min(r.mMinDihedralBefore, 10.0) - 1e-9);
}

TEST(AdaptVolume, CylinderKeepsItsRimsAndCaps) {
    const CheckInvariants check;
    const Mesh in = tet_cylinder(16, 8);  // edges 0.125
    AdaptVolumeOptions o;
    o.mHmin = 0.25;
    o.mHmax = 0.4;
    const AdaptVolumeResult r = adapt_volume(in, o);
    print_report("cylinder, h 0.25", r);
    EXPECT_LT(r.mNumTetsAfter, r.mNumTetsBefore);
    EXPECT_EQ(r.mNumCorners, 4);
    EXPECT_LE(r.mMaxBoundaryDeviation, r.mHausd);
    // the caps stay flat and every point stays inside the cylinder
    const Mesh& m = r.mMesh;
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        const double x = mio::detail::read_double(m.Points(), i * 3),
                     y = mio::detail::read_double(m.Points(), i * 3 + 1),
                     z = mio::detail::read_double(m.Points(), i * 3 + 2);
        EXPECT_LE(std::hypot(x, y), 1.0 + 1e-12);
        EXPECT_GE(z, -1e-12);
        EXPECT_LE(z, 1.0 + 1e-12);
    }
}

TEST(AdaptVolume, NestedShellsKeepTheirMaterials) {
    const CheckInvariants check;
    const Mesh a = mio::icosa_sphere(1.0, 4);
    const Mesh b = mio::icosa_sphere(0.5, 3);
    const Mesh both = mio::merge({&a, &b}, mio::MergeOptions{}).mMesh;
    mio::OctreeMeshOptions oo;
    oo.mMaxSize = 0.4;
    oo.mMinSize = 0.05;
    oo.mPaddingRelative = 0.05;
    oo.mConformity = mio::OctreeMeshConformity::Cleaving;
    const Mesh in = mio::octree_mesh(both, oo).mMesh;
    const AdaptVolumeResult r = adapt_volume(in, AdaptVolumeOptions{});
    print_report("nested shells (cleaving)", r);
    const Mesh& m = r.mMesh;
    ASSERT_TRUE(m.HasRegion("material_1", RegionKind::Cell));
    ASSERT_TRUE(m.HasRegion("material_2", RegionKind::Cell));
    const auto& m1 = m.Region(m.FindRegion("material_1", RegionKind::Cell));
    const auto& m2 = m.Region(m.FindRegion("material_2", RegionKind::Cell));
    EXPECT_GT(m1.NumEntries(), 0u);
    EXPECT_GT(m2.NumEntries(), 0u);
    EXPECT_EQ(m1.NumEntries() + m2.NumEntries(), m.Cells(0).NumCells());
    ASSERT_TRUE(m.HasCellData("octree_mesh:material"));
    // material 2 (the inner ball) stays inside the inner shell, within hausd
    // and the input's own distance from the skin (0.05 of slack)
    const NDArray& mat = m.CellData("octree_mesh:material", 0);
    const NDArray& c = m.Cells(0).Conn();
    for (std::size_t t = 0; t < m.Cells(0).NumCells(); ++t) {
        double rmax = 0.0, rmin = 1e9;
        for (std::size_t k = 0; k < 4; ++k) {
            const std::size_t v = static_cast<std::size_t>(mio::detail::read_int(c, t * 4 + k));
            const double rr = std::hypot(mio::detail::read_double(m.Points(), v * 3),
                                         mio::detail::read_double(m.Points(), v * 3 + 1),
                                         mio::detail::read_double(m.Points(), v * 3 + 2));
            rmax = std::max(rmax, rr);
            rmin = std::min(rmin, rr);
        }
        if (mio::detail::read_int(mat, t) == 2)
            EXPECT_LE(rmax, 0.55);
        else
            EXPECT_GE(rmin, 0.45);
    }
    EXPECT_LE(r.mMaxBoundaryDeviation, r.mHausd);
    for (const auto& [label, dv] : r.mRegionVolumeChange)
        EXPECT_LE(std::abs(dv), 6.0 * r.mHausd) << label;  // hausd times the shells' area
}

TEST(AdaptVolume, ASecondRunChangesLittle) {
    AdaptVolumeOptions o;
    o.mHmin = 0.3;
    o.mHmax = 0.4;
    const AdaptVolumeResult first = adapt_volume(tet_box(8), o);
    const AdaptVolumeResult second = adapt_volume(first.mMesh, o);
    print_report("box 8^3 again", second);
    // the passes keep trading edges at the 4/5 and 4/3 thresholds (measured:
    // about a tenth of them), but the mesh they leave is about the same size
    EXPECT_LE(std::abs(second.mNumTetsAfter - first.mNumTetsAfter), first.mNumTetsAfter / 10);
    EXPECT_GE(second.mMinDihedralAfter, std::min(first.mMinDihedralAfter, 10.0));
}

TEST(AdaptVolume, OutputIsDeterministic) {
    const Mesh in = octree_sphere();
    auto digest = [](const AdaptVolumeResult& rR) {
        meshioplusplus::bench::MeshDigest d;
        d.Of(rR.mMesh);
        d.U64(static_cast<std::uint64_t>(rR.mNumCollapsed));
        d.U64(static_cast<std::uint64_t>(rR.mNumSplit));
        d.U64(static_cast<std::uint64_t>(rR.mNumSwapped));
        d.U64(static_cast<std::uint64_t>(rR.mNumSmoothed));
        return d.Value();
    };
    const std::uint64_t first = digest(adapt_volume(in, AdaptVolumeOptions{}));
    for (int run = 0; run < 2; ++run)
        EXPECT_EQ(digest(adapt_volume(in, AdaptVolumeOptions{})), first);
}

TEST(AdaptVolume, RefusesATangledMesh) {
    // two tetrahedra on the same side of their shared face
    const Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.2, 0.2, 1}, {0.3, 0.3, 0.5}},
                                 "tetra", {{0, 1, 2, 3}, {0, 2, 1, 4}});
    EXPECT_THROW(adapt_volume(m), std::invalid_argument);
}

// --- M7: octree_mesh's post-stage ---------------------------------------------------

TEST(AdaptVolume, OctreeMeshPostAdaptIsOptIn) {
    const Mesh skin = mio::icosa_sphere(1.0, 4);
    mio::OctreeMeshOptions o;
    o.mMinSize = 0.04;
    const auto plain = mio::octree_mesh(skin, o);
    EXPECT_FALSE(plain.mAdapted);
    o.mPostAdapt = true;
    o.mAdapt.mHgrad = 2.0;
    const auto adapted = mio::octree_mesh(skin, o);
    ASSERT_TRUE(adapted.mAdapted);
    EXPECT_EQ(adapted.mAdapt.mNumTetsBefore, plain.mNumTetrahedra);
    EXPECT_EQ(adapted.mAdapt.mNumTetsAfter, adapted.mNumTetrahedra);
    EXPECT_LT(adapted.mNumTetrahedra, plain.mNumTetrahedra);
    EXPECT_EQ(adapted.mNumInverted, 0);
    EXPECT_FALSE(adapted.mGuaranteeApplies);
    EXPECT_EQ(adapted.mNumBoundaryEdgesOpen, 0);
    EXPECT_EQ(adapted.mNumNonManifoldBoundaryEdges, 0);
    EXPECT_GE(adapted.mMinDihedral, std::min(plain.mMinDihedral, 10.0) - 1e-9);
    EXPECT_LE(adapted.mHausdorffToSkin, plain.mHausdorffToSkin + adapted.mAdapt.mHausd);
    bool lap = false;
    for (const auto& t : adapted.mTimings)
        lap = lap || t.first == "adapt";
    EXPECT_TRUE(lap);
    std::printf(
        "[ adapt    ] octree_mesh post_adapt (hgrad 2): tets %lld -> %lld, dihedral %.2f -> %.2f, "
        "hausdorff to skin %.4g -> %.4g\n",
        static_cast<long long>(plain.mNumTetrahedra),
        static_cast<long long>(adapted.mNumTetrahedra), plain.mMinDihedral, adapted.mMinDihedral,
        plain.mHausdorffToSkin, adapted.mHausdorffToSkin);
}

// --- geometric boundary sizing ------------------------------------------------------

TEST(AdaptVolume, GeometricBoundarySizingCoarsensFlatFaces) {
    const CheckInvariants check;
    const Mesh in = tet_box(8);
    AdaptVolumeOptions o;
    o.mBoundarySizing = mio::AdaptVolumeBoundarySizing::Geometry;
    o.mHmax = 0.5;
    const AdaptVolumeResult r = adapt_volume(in, o);
    print_report("box 8^3, boundary sized by geometry", r);
    AdaptVolumeOptions input = o;
    input.mBoundarySizing = mio::AdaptVolumeBoundarySizing::Input;
    const AdaptVolumeResult ri = adapt_volume(in, input);
    print_report("box 8^3, boundary sized by the input", ri);
    EXPECT_LT(r.mNumTetsAfter, ri.mNumTetsAfter);
    EXPECT_NEAR(r.mVolumeAfter, 1.0, 1e-12);
    EXPECT_LE(r.mMaxBoundaryDeviation, r.mHausd);
    for (const double x : {0.0, 1.0})
        for (const double y : {0.0, 1.0})
            for (const double z : {0.0, 1.0})
                EXPECT_GE(find_point(r.mMesh, x, y, z), 0);
    // on a curved boundary the curvature bound and the guards hold the shape
    AdaptVolumeOptions g;
    g.mBoundarySizing = mio::AdaptVolumeBoundarySizing::Geometry;
    const AdaptVolumeResult rs = adapt_volume(octree_sphere(0.04), g);
    print_report("octree sphere, boundary sized by geometry", rs);
    EXPECT_LE(rs.mMaxBoundaryDeviation, rs.mHausd);
    EXPECT_GE(rs.mMinDihedralAfter, std::min(rs.mMinDihedralBefore, 10.0) - 1e-9);
    EXPECT_EQ(mio::adapt_volume_boundary_sizing_from_name("geometry"),
              mio::AdaptVolumeBoundarySizing::Geometry);
    EXPECT_THROW(mio::adapt_volume_boundary_sizing_from_name("x"), std::invalid_argument);
}

}  // namespace
