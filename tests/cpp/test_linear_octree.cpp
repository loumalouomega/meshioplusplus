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
// detail/linear_octree: property tests on seeded random trees. The invariants
// (keys ascending, leaves aligned and inside the forest, no overlap, full
// coverage) hold after every operation; balance equals the fixed point of the
// naive "split any leaf two levels coarser than a neighbour" rule on every
// adjacency, is idempotent and only refines; neighbour queries equal a brute
// force over all leaves; and balancing the unbalanced output of `octree`
// reproduces its own balanced output leaf for leaf.

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

// Project includes
#include "analytic_surfaces.hpp"
#include "meshioplusplus/operations/octree.hpp"
#include "../../src/cpp/src/detail/linear_octree.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::detail::LinearOctree;
using mio::detail::OctreeAdjacency;

constexpr std::int64_t kNoBudget = 0;

// A random tree: `Rounds` passes splitting each leaf with probability 1 / `Odds`,
// from raw mt19937 output (its sequence is fixed by the standard).
LinearOctree lo_random_tree(const std::array<std::int64_t, 3>& rRoots, std::int32_t MaxLevel,
                            std::uint32_t Seed, int Rounds, std::uint32_t Odds) {
    LinearOctree t = mio::detail::linear_octree_roots(rRoots, MaxLevel);
    std::mt19937 rng(Seed);
    for (int r = 0; r < Rounds; ++r) {
        std::vector<std::uint8_t> mark(t.Size());
        for (auto& m : mark)
            m = (rng() % Odds) == 0;
        mio::detail::linear_octree_split(t, mark, kNoBudget);
    }
    return t;
}

// Whether two leaves' closures meet in a way `Kind` counts: across a face (a
// shared square of positive area), an edge (a shared segment) or a corner.
bool lo_adjacent(const LinearOctree& rT, std::size_t A, std::size_t B, OctreeAdjacency Kind) {
    const auto a = rT.Lo(A), b = rT.Lo(B);
    const std::int64_t sa = rT.Span(A), sb = rT.Span(B);
    int touching = 0;  // axes where the boxes only touch
    for (std::size_t k = 0; k < 3; ++k) {
        const std::int64_t lo = std::max(a[k], b[k]), hi = std::min(a[k] + sa, b[k] + sb);
        if (hi < lo)
            return false;
        if (hi == lo)
            ++touching;
    }
    if (touching == 0)
        return false;  // overlapping: the same leaf
    return touching <= (Kind == OctreeAdjacency::Face ? 1 : Kind == OctreeAdjacency::Edge ? 2 : 3);
}

// The naive fixed point: split any leaf adjacent to one two or more levels
// finer, until none is.
LinearOctree lo_brute_balance(LinearOctree t, OctreeAdjacency Kind) {
    for (;;) {
        std::vector<std::uint8_t> mark(t.Size(), 0);
        bool any = false;
        for (std::size_t i = 0; i < t.Size(); ++i)
            for (std::size_t j = 0; j < t.Size(); ++j)
                if (t.mLevel[j] > t.mLevel[i] + 1 && lo_adjacent(t, i, j, Kind)) {
                    mark[i] = 1;
                    any = true;
                    break;
                }
        if (!any)
            return t;
        mio::detail::linear_octree_split(t, mark, kNoBudget);
    }
}

void lo_expect_valid(const LinearOctree& rT, const OctreeAdjacency* pKind) {
    std::string why;
    EXPECT_TRUE(mio::detail::linear_octree_check(rT, pKind, &why)) << why;
}

}  // namespace

TEST(LinearOctree, KeysRoundTrip) {
    for (std::int64_t x : {0, 1, 5, 1023, (1 << 21) - 1})
        for (std::int64_t y : {0, 3, 77, (1 << 20) + 9})
            for (std::int64_t z : {0, 2, 999})
                EXPECT_EQ(mio::detail::linoct_decode(mio::detail::linoct_key(x, y, z)),
                          (std::array<std::int64_t, 3>{{x, y, z}}));
}

TEST(LinearOctree, RootsTileTheForestAndLocateFindsEveryCell) {
    const LinearOctree t = mio::detail::linear_octree_roots({{3, 2, 1}}, 2);
    EXPECT_EQ(t.Size(), 6u);
    lo_expect_valid(t, nullptr);
    for (std::int64_t z = 0; z < 4; ++z)
        for (std::int64_t y = 0; y < 8; ++y)
            for (std::int64_t x = 0; x < 12; ++x) {
                const std::int64_t j = t.Locate(x, y, z);
                ASSERT_GE(j, 0);
                const auto lo = t.Lo(static_cast<std::size_t>(j));
                EXPECT_EQ(lo[0], x / 4 * 4);
                EXPECT_EQ(lo[1], y / 4 * 4);
                EXPECT_EQ(lo[2], 0);
            }
    EXPECT_EQ(t.Locate(12, 0, 0), -1);
    EXPECT_EQ(t.Locate(-1, 0, 0), -1);
}

TEST(LinearOctree, RejectsAForestTooFineForTheKeys) {
    EXPECT_THROW(mio::detail::linear_octree_roots({{3, 1, 1}}, 20), std::invalid_argument);
    EXPECT_THROW(mio::detail::linear_octree_roots({{0, 1, 1}}, 2), std::invalid_argument);
    EXPECT_NO_THROW(mio::detail::linear_octree_roots({{2, 1, 1}}, 20));
}

TEST(LinearOctree, SplitKeepsTheInvariantsAndTheBudget) {
    LinearOctree t = lo_random_tree({{2, 1, 3}}, 6, 11u, 5, 3);
    lo_expect_valid(t, nullptr);
    std::vector<std::uint8_t> all(t.Size(), 1);
    EXPECT_THROW(mio::detail::linear_octree_split(t, all, static_cast<std::int64_t>(t.Size())),
                 std::invalid_argument);
}

TEST(LinearOctree, BalanceIsTheNaiveFixedPointOnRandomTrees) {
    const OctreeAdjacency kinds[] = {OctreeAdjacency::Face, OctreeAdjacency::Edge,
                                     OctreeAdjacency::Full};
    const std::array<std::int64_t, 3> forests[] = {{{1, 1, 1}}, {{3, 2, 1}}};
    std::uint32_t seed = 1;
    for (const auto& roots : forests)
        for (int trial = 0; trial < 6; ++trial, ++seed) {
            const LinearOctree t = lo_random_tree(roots, 5, seed, 5, 4);
            for (OctreeAdjacency kind : kinds) {
                LinearOctree b = t;
                std::vector<std::int64_t> origin;
                mio::detail::linear_octree_balance(b, kind, kNoBudget, &origin);
                lo_expect_valid(b, &kind);
                const LinearOctree naive = lo_brute_balance(t, kind);
                EXPECT_EQ(b.mKey, naive.mKey) << "seed " << seed;
                EXPECT_EQ(b.mLevel, naive.mLevel) << "seed " << seed;
                // Refinement only: every new leaf lies in the leaf it came from,
                // and the pieces of a leaf add up to it.
                ASSERT_EQ(origin.size(), b.Size());
                std::vector<std::uint64_t> volume(t.Size(), 0);
                for (std::size_t i = 0; i < b.Size(); ++i) {
                    const auto o = static_cast<std::size_t>(origin[i]);
                    const auto lo = b.Lo(i), olo = t.Lo(o);
                    for (std::size_t k = 0; k < 3; ++k) {
                        EXPECT_GE(lo[k], olo[k]);
                        EXPECT_LE(lo[k] + b.Span(i), olo[k] + t.Span(o));
                    }
                    const auto s = static_cast<std::uint64_t>(b.Span(i));
                    volume[o] += s * s * s;
                }
                for (std::size_t o = 0; o < t.Size(); ++o) {
                    const auto s = static_cast<std::uint64_t>(t.Span(o));
                    EXPECT_EQ(volume[o], s * s * s);
                }
                // Idempotent.
                LinearOctree again = b;
                mio::detail::linear_octree_balance(again, kind, kNoBudget);
                EXPECT_EQ(again.mKey, b.mKey);
                EXPECT_EQ(again.mLevel, b.mLevel);
            }
        }
}

TEST(LinearOctree, NeighboursEqualABruteForce) {
    const LinearOctree t = lo_random_tree({{2, 2, 1}}, 4, 99u, 4, 3);
    const auto dirs = mio::detail::linear_octree_directions(OctreeAdjacency::Full);
    ASSERT_EQ(dirs.size(), 26u);
    ASSERT_EQ(mio::detail::linear_octree_directions(OctreeAdjacency::Face).size(), 6u);
    ASSERT_EQ(mio::detail::linear_octree_directions(OctreeAdjacency::Edge).size(), 18u);
    for (std::size_t i = 0; i < t.Size(); ++i) {
        const auto a = t.Lo(i);
        const std::int64_t sa = t.Span(i);
        for (const auto& d : dirs) {
            std::vector<std::int64_t> want;
            for (std::size_t j = 0; j < t.Size(); ++j) {
                const auto b = t.Lo(j);
                const std::int64_t sb = t.Span(j);
                // Overlaps the same-size cube on that side, and meets the leaf.
                bool ok = j != i;
                for (std::size_t k = 0; k < 3 && ok; ++k) {
                    const std::int64_t c = a[k] + d[k] * sa;
                    ok = std::max(c, b[k]) < std::min(c + sa, b[k] + sb) &&
                         std::max(a[k], b[k]) <= std::min(a[k] + sa, b[k] + sb);
                }
                if (ok)
                    want.push_back(static_cast<std::int64_t>(j));
            }
            EXPECT_EQ(mio::detail::linear_octree_neighbours(t, i, d[0], d[1], d[2]), want)
                << "leaf " << i << " dir " << d[0] << d[1] << d[2];
        }
    }
}

TEST(LinearOctree, DilateRequiresTheRingsAtTheSeedLevel) {
    LinearOctree t = mio::detail::linear_octree_roots({{1, 1, 1}}, 4);
    // Refine the leaf at the origin down to level 4.
    for (int l = 0; l < 4; ++l) {
        std::vector<std::uint8_t> mark(t.Size(), 0);
        mark[0] = 1;
        mio::detail::linear_octree_split(t, mark, kNoBudget);
    }
    std::vector<std::uint8_t> seed(t.Size(), 0);
    seed[0] = 1;
    mio::detail::linear_octree_dilate(t, seed, 2, kNoBudget);
    for (std::int64_t z = 0; z <= 2; ++z)
        for (std::int64_t y = 0; y <= 2; ++y)
            for (std::int64_t x = 0; x <= 2; ++x)
                EXPECT_EQ(t.mLevel[static_cast<std::size_t>(t.Locate(x, y, z))], 4);
    lo_expect_valid(t, nullptr);
}

TEST(LinearOctree, BalanceIsTheSameOnEveryThreadCount) {
    const LinearOctree t = lo_random_tree({{3, 3, 2}}, 7, 5u, 7, 3);
    LinearOctree ref;
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        LinearOctree b = t;
        mio::detail::linear_octree_balance(b, OctreeAdjacency::Full, kNoBudget);
        if (ref.Size() == 0)
            ref = b;
        EXPECT_EQ(b.mKey, ref.mKey);
        EXPECT_EQ(b.mLevel, ref.mLevel);
    }
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
#endif
}

namespace {

// The leaves of an `octree` result, read back from the hexahedra's lo corners.
LinearOctree lo_from_octree(const mio::OctreeResult& rR, const std::array<std::int64_t, 3>& rRoots,
                            std::int32_t MaxLevel) {
    const mio::Mesh& m = rR.mMesh;
    const double fine = rR.mRootSize / static_cast<double>(std::int64_t{1} << MaxLevel);
    const double* p = m.Points().As<double>();
    const auto block = m.Cells(0);
    const std::int64_t* conn = block.Conn().As<std::int64_t>();
    const std::int64_t* level = m.CellData(mio::kOctreeLevelName, 0).As<std::int64_t>();
    std::vector<std::array<std::int64_t, 4>> leaves;
    for (std::size_t c = 0; c < block.NumCells(); ++c) {
        std::array<std::int64_t, 4> leaf{{-1, -1, -1, level[c]}};
        for (std::size_t k = 0; k < 3; ++k) {
            double lo = p[static_cast<std::size_t>(conn[c * 8]) * 3 + k];
            for (std::size_t v = 1; v < 8; ++v)
                lo = std::min(lo, p[static_cast<std::size_t>(conn[c * 8 + v]) * 3 + k]);
            leaf[k] = std::llround((lo - rR.mOrigin[k]) / fine);
        }
        leaves.push_back(leaf);
    }
    return mio::detail::linear_octree_from_leaves(rRoots, MaxLevel, leaves);
}

}  // namespace

TEST(LinearOctree, BalancingOctreesUnbalancedTreeGivesItsBalancedOne) {
    const mio::Mesh skin = analytic::torus(1.0, 0.35, 24, 12);
    const std::array<std::int64_t, 3> roots{{2, 2, 1}};
    for (auto [kind, balance] : {std::pair{OctreeAdjacency::Face, mio::OctreeBalance::Face},
                                 std::pair{OctreeAdjacency::Full, mio::OctreeBalance::Full}}) {
        mio::OctreeOptions o;
        o.mRootCells = roots;
        o.mMaxLevel = 5;
        o.mPaddingRelative = 0.05;
        o.mBalance = mio::OctreeBalance::None;
        const mio::OctreeResult raw = mio::octree(skin, o);
        o.mBalance = balance;
        const mio::OctreeResult balanced = mio::octree(skin, o);
        LinearOctree t = lo_from_octree(raw, roots, 5);
        lo_expect_valid(t, nullptr);
        mio::detail::linear_octree_balance(t, kind, kNoBudget);
        const LinearOctree want = lo_from_octree(balanced, roots, 5);
        EXPECT_GT(want.Size(), t.Size() / 2);
        EXPECT_EQ(t.mKey, want.mKey);
        EXPECT_EQ(t.mLevel, want.mLevel);
    }
}

TEST(LinearOctree, PairingSplitsTheRootsAndLonelyLeaves) {
    // One root split down two levels in one corner: its other children are
    // leaves beside a split sibling, and the other roots have no parent.
    for (const std::uint32_t seed : {1u, 2u, 3u, 4u}) {
        LinearOctree t = lo_random_tree({{2, 1, 2}}, 4, seed, 3, 5);
        EXPECT_FALSE(mio::detail::linear_octree_is_paired(t));
        EXPECT_TRUE(mio::detail::linear_octree_pair(t, OctreeAdjacency::Full, kNoBudget));
        const OctreeAdjacency full = OctreeAdjacency::Full;
        lo_expect_valid(t, &full);
        std::string why;
        EXPECT_TRUE(mio::detail::linear_octree_is_paired(t, &why)) << why;
        // A paired tree is a fixed point.
        const LinearOctree again = t;
        EXPECT_FALSE(mio::detail::linear_octree_pair(t, OctreeAdjacency::Full, kNoBudget));
        EXPECT_EQ(t.mKey, again.mKey);
        EXPECT_EQ(t.mLevel, again.mLevel);
    }
}

TEST(LinearOctree, PairingKeepsTheBudget) {
    LinearOctree t = lo_random_tree({{1, 1, 1}}, 5, 7, 4, 3);
    EXPECT_THROW(mio::detail::linear_octree_pair(t, OctreeAdjacency::Full, 10),
                 std::invalid_argument);
}
