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
#pragma once

/**
 * @file detail/linear_octree.hpp
 * @brief A linear octree: the leaves of a forest of cubic roots as sorted
 * Morton keys, with neighbour queries and 2:1 balance, for `octree_mesh`.
 *
 * Every leaf is an aligned cube on the finest lattice (`2^L` cells per root
 * edge, `L = mMaxLevel`) and is stored as the Morton key of its lo corner
 * (`sfc_morton_key`, 21 bits per axis) and its level. An aligned cube of span
 * `s` covers the contiguous key range `[key, key + s^3)`, so the sorted keys
 * are a preorder of the forest and locating the leaf that holds a finest cell
 * is one `upper_bound`: no pointers, no hash, nothing to rebuild.
 *
 * Balance is the 2:1 rule across faces, faces and edges, or faces, edges and
 * corners, enforced finest level first (Sundar, Sampath and Biros, "Bottom-up
 * construction and 2:1 balance refinement of linear octrees in parallel",
 * SIAM J. Sci. Comput. 30(5), 2008; Burstedde, Wilcox and Ghattas, "p4est",
 * SIAM J. Sci. Comput. 33(3), 2011): for every level `l` the parents of the
 * level-`l` leaves require their neighbour octants to exist at level `l - 1`,
 * and the leaves that hold those octants are split down to them. A split at
 * step `l` only creates leaves of level `l - 1` or coarser, which the later
 * steps see, so one pass suffices. The minimal balanced refinement of a tree
 * is unique, so the result does not depend on the order of the work, and every
 * sort here is a total order: the output is identical across parallel
 * backends and thread counts.
 *
 * A **core-private** header (the `cut_graph.hpp` precedent): nothing installed
 * names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// Which neighbours the 2:1 rule covers.
enum class OctreeAdjacency {
    Face = 6,   ///< across faces
    Edge = 18,  ///< across faces and edges
    Full = 26,  ///< across faces, edges and corners
};

/// A Morton key of finest-lattice coordinates (each below `2^21`).
std::uint64_t linoct_key(std::int64_t X, std::int64_t Y, std::int64_t Z);

/// The finest-lattice coordinates of a Morton key.
std::array<std::int64_t, 3> linoct_decode(std::uint64_t Key);

/// The leaves of a forest of cubic roots.
struct LinearOctree {
    /// The finest level `L`: a root edge is `2^L` finest cells.
    std::int32_t mMaxLevel = 0;
    /// Root cells per axis.
    std::array<std::int64_t, 3> mRoots{{1, 1, 1}};
    /// The Morton key of each leaf's lo corner, strictly ascending.
    std::vector<std::uint64_t> mKey;
    /// Each leaf's level (0: a root).
    std::vector<std::uint8_t> mLevel;

    std::size_t Size() const { return mKey.size(); }
    /// Finest cells per axis of the whole forest.
    std::int64_t Extent(std::size_t Axis) const { return mRoots[Axis] << mMaxLevel; }
    /// Lo corner of leaf `I` on the finest lattice.
    std::array<std::int64_t, 3> Lo(std::size_t I) const { return linoct_decode(mKey[I]); }
    /// Edge of leaf `I` in finest cells.
    std::int64_t Span(std::size_t I) const {
        return std::int64_t{1} << (mMaxLevel - static_cast<std::int32_t>(mLevel[I]));
    }
    /// Whether the finest cell `(X, Y, Z)` lies in the forest.
    bool Contains(std::int64_t X, std::int64_t Y, std::int64_t Z) const {
        return X >= 0 && Y >= 0 && Z >= 0 && X < Extent(0) && Y < Extent(1) && Z < Extent(2);
    }
    /// The leaf that holds the finest cell `(X, Y, Z)`; -1 outside the forest.
    std::int64_t Locate(std::int64_t X, std::int64_t Y, std::int64_t Z) const;
};

/// The forest's roots, every one a leaf at level 0. Throws when an axis has no
/// root or `Roots[k] << MaxLevel` reaches `2^21`.
LinearOctree linear_octree_roots(const std::array<std::int64_t, 3>& rRoots, std::int32_t MaxLevel);

/// A tree from explicit leaves `{x, y, z, level}` (lo corners on the finest
/// lattice), sorted into key order. The leaves must tile the forest; check with
/// `linear_octree_check`.
LinearOctree linear_octree_from_leaves(const std::array<std::int64_t, 3>& rRoots,
                                       std::int32_t MaxLevel,
                                       const std::vector<std::array<std::int64_t, 4>>& rLeaves);

/// Split every leaf with `rMark[i] != 0` (and a level below `mMaxLevel`) into
/// its eight children. When @p pOrigin is given it receives, for every new
/// leaf, the index of the leaf it came from. Throws when the tree would exceed
/// @p MaxLeaves leaves.
void linear_octree_split(LinearOctree& rTree, const std::vector<std::uint8_t>& rMark,
                         std::int64_t MaxLeaves, std::vector<std::int64_t>* pOrigin = nullptr);

/// Refine the tree so that every octant `{key, level}` of @p rRequired is
/// covered by leaves at its level or finer: a leaf coarser than an octant it
/// holds is split down the path to it (its other children stay leaves).
/// @p pOrigin as for `linear_octree_split`.
void linear_octree_require(LinearOctree& rTree,
                           std::vector<std::pair<std::uint64_t, std::int32_t>> Required,
                           std::int64_t MaxLeaves, std::vector<std::int64_t>* pOrigin = nullptr);

/// The minimal 2:1-balanced refinement of the tree for @p Kind.
/// @p pOrigin as for `linear_octree_split`.
void linear_octree_balance(LinearOctree& rTree, OctreeAdjacency Kind, std::int64_t MaxLeaves,
                           std::vector<std::int64_t>* pOrigin = nullptr);

/// Require, around every leaf with `rSeed[i] != 0`, every octant of the
/// seed's own level within Chebyshev distance @p Rings (in units of the seed's
/// span) to exist at that level. Balance afterwards. @p pOrigin as above.
void linear_octree_dilate(LinearOctree& rTree, const std::vector<std::uint8_t>& rSeed,
                          std::int32_t Rings, std::int64_t MaxLeaves,
                          std::vector<std::int64_t>* pOrigin = nullptr);

/// The leaves adjacent to leaf @p Leaf across the face, edge or corner in the
/// direction `{Dx, Dy, Dz}` (each -1, 0 or 1, not all 0), in key order: the
/// leaves that overlap the cube of the leaf's size on that side and meet the
/// leaf's closure. A coarser leaf holding that cube is returned whatever else
/// it shares with the leaf (a corner query can return a face neighbour), so a
/// leaf can appear in several directions. Empty on the boundary of the forest.
/**
 * @brief Pair the tree for the dual hexahedral transition schemes (Pitzalis,
 * Livesu, Cherchi, Gobbetti and Scateni, "Generalized adaptive refinement for
 * grid-based hexahedral meshing", ACM TOG 40(6), 2021): every leaf has a parent
 * (the roots are split once) and the eight children of every parent are all
 * leaves or all split. Splits and rebalances for @p Kind until both hold.
 * @return whether anything was split.
 */
bool linear_octree_pair(LinearOctree& rTree, OctreeAdjacency Kind, std::int64_t MaxLeaves);

/// Whether every leaf has a parent and every octet of siblings is all leaves
/// or all split; @p pWhy receives the first violation.
bool linear_octree_is_paired(const LinearOctree& rTree, std::string* pWhy = nullptr);

std::vector<std::int64_t> linear_octree_neighbours(const LinearOctree& rTree, std::size_t Leaf,
                                                   int Dx, int Dy, int Dz);

/// The neighbour directions @p Kind covers: the 6 faces first, then the 12
/// edges, then the 8 corners.
std::vector<std::array<int, 3>> linear_octree_directions(OctreeAdjacency Kind);

/**
 * @brief Check the tree's invariants: keys strictly ascending, every leaf
 * aligned to its span and inside the forest, no two leaves overlapping, the
 * leaves covering every root, and (unless @p Kind is null) the 2:1 rule.
 * @param pWhy receives the first violation.
 */
bool linear_octree_check(const LinearOctree& rTree, const OctreeAdjacency* pKind,
                         std::string* pWhy = nullptr);

}  // namespace detail
}  // namespace meshioplusplus
