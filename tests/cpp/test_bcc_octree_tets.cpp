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
// detail/bcc_octree_tets on seeded random balanced trees: the tetrahedra's
// exact volumes add up to the forest's; every interior triangle is shared by
// exactly two tetrahedra that see it with opposite orientations, and every
// other one lies on the forest's boundary; every tetrahedron is positive;
// each shape class has the dihedral angles of the table, and a region of
// equal leaves is T1 only; no vertex lies inside an edge; the result is the
// same on every thread count; and an unbalanced tree is refused.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <random>
#include <stdexcept>
#include <vector>

// External includes
#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

// Project includes
#include "../../src/cpp/src/detail/bcc_octree_tets.hpp"
#include "../../src/cpp/src/detail/linear_octree.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::detail::BccLattice;
using mio::detail::LinearOctree;
using mio::detail::Vec3;

LinearOctree bt_tree(const std::array<std::int64_t, 3>& rRoots, std::int32_t MaxLevel,
                     std::uint32_t Seed, int Rounds, std::uint32_t Odds, bool Balance = true) {
    LinearOctree t = mio::detail::linear_octree_roots(rRoots, MaxLevel);
    std::mt19937 rng(Seed);
    for (int r = 0; r < Rounds; ++r) {
        std::vector<std::uint8_t> mark(t.Size());
        for (auto& m : mark)
            m = (rng() % Odds) == 0;
        mio::detail::linear_octree_split(t, mark, 0);
    }
    if (Balance)
        mio::detail::linear_octree_balance(t, mio::detail::OctreeAdjacency::Full, 0);
    return t;
}

Vec3 bt_pos(const BccLattice& rL, std::int64_t V) {
    const auto& p = rL.mIjk2[static_cast<std::size_t>(V)];
    return {static_cast<double>(p[0]), static_cast<double>(p[1]), static_cast<double>(p[2])};
}

// Six times the signed volume, exactly (small doubled coordinates).
std::int64_t bt_vol6(const BccLattice& rL, const std::array<std::int64_t, 4>& rT) {
    std::int64_t d[3][3];
    for (std::size_t r = 0; r < 3; ++r)
        for (std::size_t k = 0; k < 3; ++k)
            d[r][k] = rL.mIjk2[static_cast<std::size_t>(rT[r + 1])][k] -
                      rL.mIjk2[static_cast<std::size_t>(rT[0])][k];
    return d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) -
           d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0]) +
           d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]);
}

void bt_check(const LinearOctree& rTree) {
    const BccLattice L = mio::detail::bcc_octree_lattice(rTree, 0);
    ASSERT_FALSE(L.mTets.empty());
    // Keys ascending, coordinates consistent.
    for (std::size_t p = 1; p < L.mKey.size(); ++p)
        ASSERT_LT(L.mKey[p - 1], L.mKey[p]);
    // Positive volumes adding up to the forest (doubled: 8x the volume).
    std::int64_t vol = 0;
    for (const auto& t : L.mTets) {
        const std::int64_t v = bt_vol6(L, t);
        ASSERT_GT(v, 0);
        vol += v;
    }
    const std::int64_t box = 8 * rTree.Extent(0) * rTree.Extent(1) * rTree.Extent(2);
    EXPECT_EQ(vol, 6 * box);
    // Every triangle: twice with opposite orientation, or once on the boundary.
    std::map<std::array<std::int64_t, 3>, std::vector<int>> faces;
    constexpr int kFace[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    for (const auto& t : L.mTets)
        for (const auto& f : kFace) {
            std::array<std::int64_t, 3> v{t[f[0]], t[f[1]], t[f[2]]};
            // Orientation parity of the sort.
            int parity = 0;
            for (int a = 0; a < 3; ++a)
                for (int b = a + 1; b < 3; ++b)
                    parity ^= v[a] > v[b];
            std::sort(v.begin(), v.end());
            faces[v].push_back(parity ? -1 : 1);
        }
    for (const auto& [key, uses] : faces) {
        if (uses.size() == 2) {
            EXPECT_EQ(uses[0] + uses[1], 0) << "a face seen with the same orientation twice";
            continue;
        }
        ASSERT_EQ(uses.size(), 1u) << "a face shared by " << uses.size() << " tetrahedra";
        // On the forest's boundary: one coordinate the same at 0 or the extent.
        bool on = false;
        for (std::size_t k = 0; k < 3 && !on; ++k) {
            const std::int64_t c = L.mIjk2[static_cast<std::size_t>(key[0])][k];
            on = (c == 0 || c == 2 * rTree.Extent(k)) &&
                 L.mIjk2[static_cast<std::size_t>(key[1])][k] == c &&
                 L.mIjk2[static_cast<std::size_t>(key[2])][k] == c;
        }
        EXPECT_TRUE(on) << "an open face inside the forest";
    }
    // No vertex strictly inside an edge.
    std::map<std::array<std::int64_t, 2>, int> edges;
    for (const auto& t : L.mTets)
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b)
                edges[{std::min(t[a], t[b]), std::max(t[a], t[b])}] = 1;
    for (const auto& [e, unused] : edges) {
        const auto& a = L.mIjk2[static_cast<std::size_t>(e[0])];
        const auto& b = L.mIjk2[static_cast<std::size_t>(e[1])];
        for (std::size_t p = 0; p < L.mIjk2.size(); ++p) {
            if (static_cast<std::int64_t>(p) == e[0] || static_cast<std::int64_t>(p) == e[1])
                continue;
            const auto& q = L.mIjk2[p];
            bool inside = true;
            std::int64_t cross = 0;
            for (std::size_t k = 0; k < 3 && inside; ++k)
                inside = q[k] >= std::min(a[k], b[k]) && q[k] <= std::max(a[k], b[k]);
            if (!inside)
                continue;
            const std::int64_t u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const std::int64_t w[3] = {q[0] - a[0], q[1] - a[1], q[2] - a[2]};
            cross = std::abs(u[1] * w[2] - u[2] * w[1]) + std::abs(u[2] * w[0] - u[0] * w[2]) +
                    std::abs(u[0] * w[1] - u[1] * w[0]);
            EXPECT_NE(cross, 0) << "vertex " << p << " inside an edge";
        }
    }
}

// The dihedral angles of each shape class.
const std::array<double, 6> kBtAngles[7] = {
    {},
    {60, 60, 60, 60, 90, 90},
    {45, 60, 60, 90, 90, 90},
    {45, 60, 60, 90, 90, 90},
    {45, 45, 60, 90, 90, 90},
    {35.2643896827546, 45, 60, 90, 90, 125.264389682754},
    {45, 54.7356103172453, 54.7356103172453, 90, 90, 90},
};

}  // namespace

TEST(BccOctreeTets, RandomBalancedTreesAreConforming) {
    const std::array<std::int64_t, 3> forests[] = {{{1, 1, 1}}, {{2, 1, 2}}};
    std::uint32_t seed = 3;
    for (const auto& roots : forests)
        for (int trial = 0; trial < 4; ++trial, ++seed)
            bt_check(bt_tree(roots, 4, seed, 4, 3));
}

TEST(BccOctreeTets, ShapesHaveTheirTabulatedAngles) {
    const LinearOctree t = bt_tree({{2, 2, 2}}, 4, 21u, 4, 3);
    const BccLattice L = mio::detail::bcc_octree_lattice(t, 0);
    std::array<int, 7> seen{};
    for (std::size_t i = 0; i < L.mTets.size(); ++i) {
        const auto& v = L.mTets[i];
        const auto a = mio::detail::tet_dihedral_angles(bt_pos(L, v[0]), bt_pos(L, v[1]),
                                                        bt_pos(L, v[2]), bt_pos(L, v[3]));
        const int s = L.mShape[i];
        ASSERT_GE(s, 1);
        ASSERT_LE(s, 6);
        ++seen[static_cast<std::size_t>(s)];
        for (std::size_t k = 0; k < 6; ++k)
            EXPECT_NEAR(a[k], kBtAngles[s][k], 1e-9) << "class T" << s;
    }
    for (int s = 1; s <= 6; ++s)
        EXPECT_GT(seen[static_cast<std::size_t>(s)], 0) << "class T" << s << " never occurs";
}

TEST(BccOctreeTets, EqualLeavesAreT1AwayFromTheBoundary) {
    const LinearOctree t = mio::detail::linear_octree_roots({{3, 3, 3}}, 0);
    const BccLattice L = mio::detail::bcc_octree_lattice(t, 0);
    // 54 interior faces of 4 BCC tetrahedra, 54 boundary faces of 4 half-BCC ones.
    EXPECT_EQ(L.mTets.size(), 54u * 4u + 54u * 4u);
    for (std::size_t i = 0; i < L.mTets.size(); ++i) {
        bool boundary = false;
        for (std::int64_t v : L.mTets[i])
            boundary = boundary || L.mType[static_cast<std::size_t>(v)] ==
                                       mio::detail::BccVertexType::FaceCentre;
        EXPECT_EQ(L.mShape[i], boundary ? 3 : 1);
    }
}

TEST(BccOctreeTets, IsTheSameOnEveryThreadCount) {
    const LinearOctree t = bt_tree({{2, 2, 1}}, 6, 9u, 6, 3);
    BccLattice ref;
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const BccLattice L = mio::detail::bcc_octree_lattice(t, 0);
        if (ref.mTets.empty())
            ref = L;
        EXPECT_EQ(L.mKey, ref.mKey);
        EXPECT_EQ(L.mTets, ref.mTets);
        EXPECT_EQ(L.mShape, ref.mShape);
    }
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
#endif
}

TEST(BccOctreeTets, RefusesAnUnbalancedTreeAndABudget) {
    const LinearOctree t = bt_tree({{1, 1, 1}}, 5, 4u, 5, 2, /*Balance=*/false);
    std::string why;
    const auto full = mio::detail::OctreeAdjacency::Full;
    ASSERT_FALSE(mio::detail::linear_octree_check(t, &full, &why));
    EXPECT_THROW(mio::detail::bcc_octree_lattice(t, 0), std::logic_error);
    const LinearOctree b = bt_tree({{1, 1, 1}}, 3, 4u, 3, 2);
    EXPECT_THROW(mio::detail::bcc_octree_lattice(b, 10), std::invalid_argument);
}
