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

// Unit tests of adapt_volume's private building blocks: the quality kernels
// of detail/fem_kernels.hpp and the store and
// classification of detail/tet_remesh.hpp.

// System includes
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/fem_kernels.hpp"
#include "../../src/cpp/src/detail/tet_remesh.hpp"

namespace {

using meshioplusplus::detail::TetCorners;
using meshioplusplus::detail::TetRemeshStore;
using meshioplusplus::detail::TetVertexKind;
using meshioplusplus::detail::TriKey;
using meshioplusplus::detail::TriKeyHash;
using meshioplusplus::detail::Vec3;

const Vec3 kRegular[4] = {{1.0, 1.0, 1.0}, {1.0, -1.0, -1.0}, {-1.0, 1.0, -1.0}, {-1.0, -1.0, 1.0}};

TEST(TetQualityKernels, RegularTetrahedronIsIdeal) {
    namespace d = meshioplusplus::detail;
    // (1,1,1) (1,-1,-1) (-1,1,-1) (-1,-1,1) is negatively oriented; swap two.
    const Vec3 a = kRegular[0], b = kRegular[1], c = kRegular[3], e = kRegular[2];
    ASSERT_GT(d::tet_volume6(a, b, c, e), 0.0);
    EXPECT_NEAR(d::tet_mean_ratio(a, b, c, e), 1.0, 1e-12);
    EXPECT_NEAR(d::tet_amips(a, b, c, e), 1.0, 1e-12);
    double mn = 0.0, mx = 0.0;
    d::tet_dihedral_range(a, b, c, e, mn, mx);
    const double regular = std::acos(1.0 / 3.0) * 180.0 / 3.14159265358979323846;
    EXPECT_NEAR(mn, regular, 1e-9);
    EXPECT_NEAR(mx, regular, 1e-9);
    EXPECT_NEAR(d::tet_min_dihedral(a, b, c, e), regular, 1e-9);
    // inverted: no mean ratio, infinite energy; the angles do not see the sign
    EXPECT_EQ(d::tet_mean_ratio(a, b, e, c), 0.0);
    EXPECT_TRUE(std::isinf(d::tet_amips(a, b, e, c)));
    EXPECT_NEAR(d::tet_min_dihedral(a, b, e, c), regular, 1e-9);
}

TEST(TetQualityKernels, FlatteningLowersMeanRatioAndAngles) {
    namespace d = meshioplusplus::detail;
    const Vec3 a{0, 0, 0}, b{1, 0, 0}, c{0, 1, 0};
    double last = 2.0, last_angle = 200.0;
    for (const double h : {1.0, 0.3, 0.1, 0.01}) {
        const Vec3 e{0.3, 0.3, h};
        const double q = d::tet_mean_ratio(a, b, c, e);
        const double angle = d::tet_min_dihedral(a, b, c, e);
        EXPECT_LT(q, last);
        EXPECT_LT(angle, last_angle);
        last = q;
        last_angle = angle;
    }
    EXPECT_EQ(d::tet_min_dihedral(a, b, c, Vec3{0.2, 0.2, 0.0}), 0.0);
}

TEST(TetRemeshStore, OrientationFollowsMeshioConvention) {
    const Vec3 a{0, 0, 0}, b{1, 0, 0}, c{0, 1, 0}, d{0, 0, 1};
    EXPECT_TRUE(TetRemeshStore::Positive(a, b, c, d));
    EXPECT_EQ(TetRemeshStore::OrientSign(a, b, c, d), 1);
    EXPECT_EQ(TetRemeshStore::OrientSign(a, c, b, d), -1);
    EXPECT_EQ(TetRemeshStore::OrientSign(a, b, c, Vec3{0.3, 0.3, 0.0}), 0);
}

// The 6-tet Kuhn split of the unit cube.
std::vector<TetCorners> kuhn_cube(std::vector<double>& rXyz) {
    rXyz.clear();
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i)
                rXyz.insert(rXyz.end(), {double(i), double(j), double(k)});
    // vertex id = i + 2 j + 4 k; the six paths from 0 to 7
    return {{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
}

TEST(TetRemeshStore, BuildFindsBoundaryAndClassifiesACube) {
    std::vector<double> xyz;
    std::vector<TetCorners> tets = kuhn_cube(xyz);
    TetRemeshStore store;
    store.Build(xyz, tets, std::vector<std::int32_t>(tets.size(), 0), {});
    EXPECT_EQ(store.mBoundary.size(), 12u);  // two triangles per side
    std::string why;
    EXPECT_TRUE(store.CheckInvariants(why)) << why;
    const auto report = store.Classify(30.0, 30.0, true);
    EXPECT_EQ(report.mNumPatches, 6);
    // every vertex of a cube is a corner; the diagonals of the sides are not ridges
    EXPECT_EQ(report.mNumCorners, 8);
    EXPECT_EQ(report.mNumChains, 12);
    EXPECT_EQ(store.mRidgeEdges.size(), 12u);
    EXPECT_EQ(store.RidgeChain(0, 3), -1);  // a side's diagonal
    EXPECT_GE(store.RidgeChain(0, 1), 0);
    // without features only the closed outer boundary is left: one patch
    TetRemeshStore smooth;
    smooth.Build(xyz, tets, std::vector<std::int32_t>(tets.size(), 0), {});
    const auto r2 = smooth.Classify(30.0, 30.0, false);
    EXPECT_EQ(r2.mNumPatches, 1);
    EXPECT_EQ(r2.mNumCorners, 0);
    for (std::int64_t v = 0; v < 8; ++v)
        EXPECT_EQ(smooth.Kind(v), TetVertexKind::Surface);
}

TEST(TetRemeshStore, MaterialInterfaceIsAnInternalBoundary) {
    // two cubes side by side, x in [0,1] class 0 and [1,2] class 1
    std::vector<double> xyz;
    std::vector<TetCorners> a = kuhn_cube(xyz);
    const std::size_t n = xyz.size() / 3;
    std::vector<double> all = xyz;
    for (std::size_t i = 0; i < n; ++i)
        all.insert(all.end(), {xyz[i * 3] + 1.0, xyz[i * 3 + 1], xyz[i * 3 + 2]});
    // weld the shared side x = 1: vertex i + 1 of the first cube (i even) is i of the second
    std::vector<std::int64_t> map(n);
    for (std::size_t i = 0; i < n; ++i)
        map[i] = (i % 2 == 0) ? static_cast<std::int64_t>(i + 1) : static_cast<std::int64_t>(n + i);
    std::vector<TetCorners> tets = a;
    std::vector<std::int32_t> cls(a.size(), 0);
    for (const TetCorners& t : a) {
        tets.push_back({map[static_cast<std::size_t>(t[0])], map[static_cast<std::size_t>(t[1])],
                        map[static_cast<std::size_t>(t[2])], map[static_cast<std::size_t>(t[3])]});
        cls.push_back(1);
    }
    TetRemeshStore store;
    store.Build(all, tets, cls, {});
    std::string why;
    ASSERT_TRUE(store.CheckInvariants(why)) << why;
    std::size_t interface = 0;
    for (const auto& kv : store.mBoundary)
        interface += kv.second.mClassB >= 0;
    EXPECT_EQ(interface, 2u);
    EXPECT_EQ(store.mBoundary.size(), 22u);  // 10 outer sides x 2 + the interface's 2
    // with features off, the rim of the interface is still a ridge (three faces meet there)
    const auto r = store.Classify(30.0, 30.0, false);
    EXPECT_EQ(r.mNumPatches, 3);  // left skin, right skin, interface
    for (const std::int64_t v : {1, 3, 5, 7})
        EXPECT_NE(store.Kind(v), TetVertexKind::Surface) << v;
}

TEST(TetRemeshStore, EditsKeepIncidenceConsistent) {
    std::vector<double> xyz;
    std::vector<TetCorners> tets = kuhn_cube(xyz);
    TetRemeshStore store;
    store.Build(xyz, tets, std::vector<std::int32_t>(tets.size(), 0), {});
    store.Classify(30.0, 30.0, true);
    std::vector<std::int64_t> shared;
    store.SharedTets(0, 7, shared);
    EXPECT_EQ(shared.size(), 6u);  // the main diagonal is in every tet
    // split the diagonal at its midpoint: each tet becomes two
    const std::int64_t m =
        store.AddVertex({0.5, 0.5, 0.5}, 1.0, TetVertexKind::Interior, -1, 0.0, 0);
    for (const std::int64_t t : shared) {
        TetCorners c = store.mTets[static_cast<std::size_t>(t)];
        TetCorners lo = c, hi = c;
        for (int k = 0; k < 4; ++k) {
            if (c[k] == 7)
                lo[k] = m;
            if (c[k] == 0)
                hi[k] = m;
        }
        store.KillTet(t);
        store.AddTet(lo, 0);
        store.AddTet(hi, 0);
    }
    std::string why;
    EXPECT_TRUE(store.CheckInvariants(why)) << why;
    EXPECT_EQ(store.NumAliveTets(), 12);
    std::vector<std::int64_t> ring;
    store.VertexRing(m, ring);
    EXPECT_EQ(ring.size(), 8u);
    EXPECT_FALSE(store.OnBoundary(m));
    EXPECT_TRUE(store.OnBoundary(0));
}

TEST(TetRemeshStore, ChainEvaluatesByArcLength) {
    meshioplusplus::detail::TetRidgeChain c;
    c.mPoints = {{0, 0, 0}, {1, 0, 0}, {1, 2, 0}};
    c.mArc = {0.0, 1.0, 3.0};
    const Vec3 p = c.Eval(2.0);
    EXPECT_DOUBLE_EQ(p[0], 1.0);
    EXPECT_DOUBLE_EQ(p[1], 1.0);
    EXPECT_DOUBLE_EQ(c.Eval(-1.0)[0], 0.0);
    EXPECT_DOUBLE_EQ(c.Eval(10.0)[1], 2.0);
}

}  // namespace
