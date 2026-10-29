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
// The linear octree of detail/linear_octree.hpp. Anonymous-namespace helpers
// are prefixed `linoct_` (the amalgamation concatenates every translation unit
// into one).

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "linear_octree.hpp"
#include "meshioplusplus/detail/space_filling.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr const char* kLinoctPrefix = "meshio++: linear_octree: ";

using LinoctOctant = std::pair<std::uint64_t, std::int32_t>;

std::uint64_t linoct_compact(std::uint64_t V) {
    V &= 0x1249249249249249ULL;
    V = (V | (V >> 2)) & 0x10c30c30c30c30c3ULL;
    V = (V | (V >> 4)) & 0x100f00f00f00f00fULL;
    V = (V | (V >> 8)) & 0x1f0000ff0000ffULL;
    V = (V | (V >> 16)) & 0x1f00000000ffffULL;
    V = (V | (V >> 32)) & 0x1fffffULL;
    return V;
}

// Finest cells in an aligned cube of edge `Span` (a power of two).
std::uint64_t linoct_volume(std::int64_t Span) {
    const auto s = static_cast<std::uint64_t>(Span);
    return s * s * s;
}

std::int64_t linoct_span(std::int32_t MaxLevel, std::int32_t Level) {
    return std::int64_t{1} << (MaxLevel - Level);
}

void linoct_check_budget(std::size_t Count, std::int64_t MaxLeaves) {
    if (MaxLeaves > 0 && Count > static_cast<std::size_t>(MaxLeaves))
        throw std::invalid_argument(std::string(kLinoctPrefix) + "the tree would have " +
                                    std::to_string(Count) + " leaves, more than max_cells (" +
                                    std::to_string(MaxLeaves) + ")");
}

// The leaves an aligned octant `{Key, Level}` of a leaf `{LeafKey, LeafLevel}`
// is refined into so that every octant of `[First, Last)` (sorted, all inside
// the leaf) is covered at its level or finer: the leaf is split while a
// required octant inside the current cube is finer than the cube.
void linoct_refine(std::int32_t MaxLevel, std::uint64_t Key, std::int32_t Level,
                   const LinoctOctant* pFirst, const LinoctOctant* pLast,
                   std::vector<LinoctOctant>& rOut) {
    bool split = false;
    for (const LinoctOctant* p = pFirst; p != pLast; ++p)
        if (p->second > Level) {
            split = true;
            break;
        }
    if (!split) {
        rOut.emplace_back(Key, Level);
        return;
    }
    const std::uint64_t child = linoct_volume(linoct_span(MaxLevel, Level + 1));
    const LinoctOctant* p = pFirst;
    for (std::uint64_t c = 0; c < 8; ++c) {
        const std::uint64_t lo = Key + c * child, hi = lo + child;
        const LinoctOctant* q = p;
        while (q != pLast && q->first < hi)
            ++q;
        linoct_refine(MaxLevel, lo, Level + 1, p, q, rOut);
        p = q;
    }
}

// Replace leaf `i` by `rRepl[g]` for every group `g` holding leaf `rGroupLeaf[g]`
// (ascending), keeping every other leaf.
void linoct_splice(LinearOctree& rTree, const std::vector<std::size_t>& rGroupLeaf,
                   const std::vector<std::vector<LinoctOctant>>& rRepl, std::int64_t MaxLeaves,
                   std::vector<std::int64_t>* pOrigin) {
    const std::size_t n = rTree.Size();
    std::vector<std::int64_t> count(n, 1);
    for (std::size_t g = 0; g < rGroupLeaf.size(); ++g)
        count[rGroupLeaf[g]] = static_cast<std::int64_t>(rRepl[g].size());
    std::vector<std::int64_t> offset(n);
    const std::int64_t total =
        parallel_exclusive_scan(count.data(), n, offset.data(), std::int64_t{0});
    linoct_check_budget(static_cast<std::size_t>(total), MaxLeaves);
    std::vector<std::int64_t> group(n, -1);
    for (std::size_t g = 0; g < rGroupLeaf.size(); ++g)
        group[rGroupLeaf[g]] = static_cast<std::int64_t>(g);
    std::vector<std::uint64_t> key(static_cast<std::size_t>(total));
    std::vector<std::uint8_t> level(static_cast<std::size_t>(total));
    std::vector<std::int64_t> origin(pOrigin ? static_cast<std::size_t>(total) : 0);
    parallel_for(n, [&](std::size_t i) {
        const auto at = static_cast<std::size_t>(offset[i]);
        if (group[i] < 0) {
            key[at] = rTree.mKey[i];
            level[at] = rTree.mLevel[i];
            if (pOrigin)
                origin[at] = static_cast<std::int64_t>(i);
            return;
        }
        const auto& repl = rRepl[static_cast<std::size_t>(group[i])];
        for (std::size_t k = 0; k < repl.size(); ++k) {
            key[at + k] = repl[k].first;
            level[at + k] = static_cast<std::uint8_t>(repl[k].second);
            if (pOrigin)
                origin[at + k] = static_cast<std::int64_t>(i);
        }
    });
    rTree.mKey.swap(key);
    rTree.mLevel.swap(level);
    if (pOrigin)
        pOrigin->swap(origin);
}

// `rOuter[i] = rOuter[rInner[i]]`: compose a new-to-old map with an earlier one.
void linoct_compose(std::vector<std::int64_t>& rOuter, const std::vector<std::int64_t>& rInner) {
    std::vector<std::int64_t> out(rInner.size());
    for (std::size_t i = 0; i < rInner.size(); ++i)
        out[i] = rOuter[static_cast<std::size_t>(rInner[i])];
    rOuter.swap(out);
}

}  // namespace

std::uint64_t linoct_key(std::int64_t X, std::int64_t Y, std::int64_t Z) {
    const std::uint32_t q[3] = {static_cast<std::uint32_t>(X), static_cast<std::uint32_t>(Y),
                                static_cast<std::uint32_t>(Z)};
    return sfc_morton_key(q);
}

std::array<std::int64_t, 3> linoct_decode(std::uint64_t Key) {
    return {{static_cast<std::int64_t>(linoct_compact(Key)),
             static_cast<std::int64_t>(linoct_compact(Key >> 1)),
             static_cast<std::int64_t>(linoct_compact(Key >> 2))}};
}

std::int64_t LinearOctree::Locate(std::int64_t X, std::int64_t Y, std::int64_t Z) const {
    if (!Contains(X, Y, Z) || mKey.empty())
        return -1;
    const std::uint64_t key = linoct_key(X, Y, Z);
    auto it = std::upper_bound(mKey.begin(), mKey.end(), key);
    if (it == mKey.begin())
        return -1;
    const auto i = static_cast<std::size_t>(it - mKey.begin()) - 1;
    if (key - mKey[i] >= linoct_volume(Span(i)))
        return -1;  // a gap: the leaves do not tile the forest
    return static_cast<std::int64_t>(i);
}

LinearOctree linear_octree_roots(const std::array<std::int64_t, 3>& rRoots, std::int32_t MaxLevel) {
    if (MaxLevel < 0 || MaxLevel > sfc_bits - 1)
        throw std::invalid_argument(std::string(kLinoctPrefix) + "max_level must be in [0, " +
                                    std::to_string(sfc_bits - 1) + "]");
    for (std::size_t k = 0; k < 3; ++k)
        if (rRoots[k] < 1 || (rRoots[k] << MaxLevel) > (std::int64_t{1} << sfc_bits))
            throw std::invalid_argument(std::string(kLinoctPrefix) +
                                        "root_cells << max_level must be in [1, 2^21] per axis");
    std::vector<std::array<std::int64_t, 4>> leaves;
    const std::int64_t root = std::int64_t{1} << MaxLevel;
    for (std::int64_t z = 0; z < rRoots[2]; ++z)
        for (std::int64_t y = 0; y < rRoots[1]; ++y)
            for (std::int64_t x = 0; x < rRoots[0]; ++x)
                leaves.push_back({{x * root, y * root, z * root, 0}});
    return linear_octree_from_leaves(rRoots, MaxLevel, leaves);
}

LinearOctree linear_octree_from_leaves(const std::array<std::int64_t, 3>& rRoots,
                                       std::int32_t MaxLevel,
                                       const std::vector<std::array<std::int64_t, 4>>& rLeaves) {
    LinearOctree tree;
    tree.mMaxLevel = MaxLevel;
    tree.mRoots = rRoots;
    std::vector<std::pair<std::uint64_t, std::uint8_t>> leaves(rLeaves.size());
    for (std::size_t i = 0; i < rLeaves.size(); ++i) {
        const auto& l = rLeaves[i];
        if (l[3] < 0 || l[3] > MaxLevel)
            throw std::invalid_argument(std::string(kLinoctPrefix) +
                                        "a leaf level is out of range");
        leaves[i] = {linoct_key(l[0], l[1], l[2]), static_cast<std::uint8_t>(l[3])};
    }
    parallel_sort(leaves.begin(), leaves.end());
    tree.mKey.resize(leaves.size());
    tree.mLevel.resize(leaves.size());
    for (std::size_t i = 0; i < leaves.size(); ++i) {
        tree.mKey[i] = leaves[i].first;
        tree.mLevel[i] = leaves[i].second;
    }
    return tree;
}

void linear_octree_split(LinearOctree& rTree, const std::vector<std::uint8_t>& rMark,
                         std::int64_t MaxLeaves, std::vector<std::int64_t>* pOrigin) {
    const std::size_t n = rTree.Size();
    std::vector<std::int64_t> count(n);
    for (std::size_t i = 0; i < n; ++i)
        count[i] = (rMark[i] && rTree.mLevel[i] < rTree.mMaxLevel) ? 8 : 1;
    std::vector<std::int64_t> offset(n);
    const std::int64_t total =
        parallel_exclusive_scan(count.data(), n, offset.data(), std::int64_t{0});
    linoct_check_budget(static_cast<std::size_t>(total), MaxLeaves);
    std::vector<std::uint64_t> key(static_cast<std::size_t>(total));
    std::vector<std::uint8_t> level(static_cast<std::size_t>(total));
    std::vector<std::int64_t> origin(pOrigin ? static_cast<std::size_t>(total) : 0);
    parallel_for(n, [&](std::size_t i) {
        const auto at = static_cast<std::size_t>(offset[i]);
        if (count[i] == 1) {
            key[at] = rTree.mKey[i];
            level[at] = rTree.mLevel[i];
            if (pOrigin)
                origin[at] = static_cast<std::int64_t>(i);
            return;
        }
        const auto lvl = static_cast<std::uint8_t>(rTree.mLevel[i] + 1);
        const std::uint64_t child = linoct_volume(linoct_span(rTree.mMaxLevel, lvl));
        for (std::size_t c = 0; c < 8; ++c) {
            key[at + c] = rTree.mKey[i] + c * child;
            level[at + c] = lvl;
            if (pOrigin)
                origin[at + c] = static_cast<std::int64_t>(i);
        }
    });
    rTree.mKey.swap(key);
    rTree.mLevel.swap(level);
    if (pOrigin)
        pOrigin->swap(origin);
}

void linear_octree_require(LinearOctree& rTree, std::vector<LinoctOctant> Required,
                           std::int64_t MaxLeaves, std::vector<std::int64_t>* pOrigin) {
    parallel_sort(Required.begin(), Required.end());
    Required.erase(std::unique(Required.begin(), Required.end()), Required.end());
    // Keep the octants a coarser leaf holds, with that leaf.
    std::vector<std::int64_t> holder(Required.size());
    parallel_for(Required.size(), [&](std::size_t r) {
        const auto lo = linoct_decode(Required[r].first);
        const std::int64_t j = rTree.Locate(lo[0], lo[1], lo[2]);
        holder[r] =
            (j >= 0 && rTree.mLevel[static_cast<std::size_t>(j)] < Required[r].second) ? j : -1;
    });
    std::vector<LinoctOctant> need;
    std::vector<std::size_t> need_leaf;
    for (std::size_t r = 0; r < Required.size(); ++r)
        if (holder[r] >= 0) {
            need.push_back(Required[r]);
            need_leaf.push_back(static_cast<std::size_t>(holder[r]));
        }
    if (need.empty()) {
        if (pOrigin) {
            pOrigin->resize(rTree.Size());
            for (std::size_t i = 0; i < rTree.Size(); ++i)
                (*pOrigin)[i] = static_cast<std::int64_t>(i);
        }
        return;
    }
    // Octants sorted by key are grouped by their holding leaf, in leaf order.
    std::vector<std::size_t> group_leaf, group_begin;
    for (std::size_t r = 0; r < need.size(); ++r)
        if (r == 0 || need_leaf[r] != need_leaf[r - 1]) {
            group_leaf.push_back(need_leaf[r]);
            group_begin.push_back(r);
        }
    group_begin.push_back(need.size());
    std::vector<std::vector<LinoctOctant>> repl(group_leaf.size());
    parallel_for(
        group_leaf.size(),
        [&](std::size_t g) {
            const std::size_t leaf = group_leaf[g];
            linoct_refine(rTree.mMaxLevel, rTree.mKey[leaf], rTree.mLevel[leaf],
                          need.data() + group_begin[g], need.data() + group_begin[g + 1], repl[g]);
        },
        64);
    linoct_splice(rTree, group_leaf, repl, MaxLeaves, pOrigin);
}

std::vector<std::array<int, 3>> linear_octree_directions(OctreeAdjacency Kind) {
    std::vector<std::array<int, 3>> dirs;
    const int want = static_cast<int>(Kind);
    for (int nz = 1; nz <= 3; ++nz) {  // faces, then edges, then corners
        if ((nz == 2 && want < 18) || (nz == 3 && want < 26))
            break;
        for (int dz = -1; dz <= 1; ++dz)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx)
                    if ((dx != 0) + (dy != 0) + (dz != 0) == nz)
                        dirs.push_back({{dx, dy, dz}});
    }
    return dirs;
}

void linear_octree_balance(LinearOctree& rTree, OctreeAdjacency Kind, std::int64_t MaxLeaves,
                           std::vector<std::int64_t>* pOrigin) {
    const auto dirs = linear_octree_directions(Kind);
    std::vector<std::int64_t> total;
    if (pOrigin) {
        total.resize(rTree.Size());
        for (std::size_t i = 0; i < rTree.Size(); ++i)
            total[i] = static_cast<std::int64_t>(i);
    }
    for (std::int32_t l = rTree.mMaxLevel; l >= 2; --l) {
        // The parents of the level-l leaves, ascending and unique.
        const std::uint64_t parent_vol = linoct_volume(linoct_span(rTree.mMaxLevel, l - 1));
        std::vector<std::uint64_t> parents;
        for (std::size_t i = 0; i < rTree.Size(); ++i)
            if (rTree.mLevel[i] == l) {
                const std::uint64_t p = rTree.mKey[i] - rTree.mKey[i] % parent_vol;
                if (parents.empty() || parents.back() != p)
                    parents.push_back(p);
            }
        if (parents.empty())
            continue;
        const std::int64_t ps = linoct_span(rTree.mMaxLevel, l - 1);
        const std::size_t nd = dirs.size();
        std::vector<LinoctOctant> req(parents.size() * nd, {0, -1});
        parallel_for(parents.size(), [&](std::size_t p) {
            const auto lo = linoct_decode(parents[p]);
            for (std::size_t d = 0; d < nd; ++d) {
                const std::int64_t x = lo[0] + dirs[d][0] * ps, y = lo[1] + dirs[d][1] * ps,
                                   z = lo[2] + dirs[d][2] * ps;
                if (rTree.Contains(x, y, z))
                    req[p * nd + d] = {linoct_key(x, y, z), l - 1};
            }
        });
        req.erase(std::remove_if(req.begin(), req.end(),
                                 [](const LinoctOctant& r) { return r.second < 0; }),
                  req.end());
        std::vector<std::int64_t> step;
        linear_octree_require(rTree, std::move(req), MaxLeaves, pOrigin ? &step : nullptr);
        if (pOrigin)
            linoct_compose(total, step);
    }
    if (pOrigin)
        pOrigin->swap(total);
}

void linear_octree_dilate(LinearOctree& rTree, const std::vector<std::uint8_t>& rSeed,
                          std::int32_t Rings, std::int64_t MaxLeaves,
                          std::vector<std::int64_t>* pOrigin) {
    std::vector<LinoctOctant> req;
    for (std::size_t i = 0; i < rTree.Size(); ++i) {
        if (!rSeed[i])
            continue;
        const auto lo = rTree.Lo(i);
        const std::int64_t s = rTree.Span(i);
        for (std::int32_t c = -Rings; c <= Rings; ++c)
            for (std::int32_t b = -Rings; b <= Rings; ++b)
                for (std::int32_t a = -Rings; a <= Rings; ++a) {
                    const std::int64_t x = lo[0] + a * s, y = lo[1] + b * s, z = lo[2] + c * s;
                    if ((a || b || c) && rTree.Contains(x, y, z))
                        req.emplace_back(linoct_key(x, y, z), rTree.mLevel[i]);
                }
    }
    linear_octree_require(rTree, std::move(req), MaxLeaves, pOrigin);
}

namespace {

// The leaves that break the pairing: roots, and leaves with a split sibling.
std::vector<std::uint8_t> linoct_unpaired(const LinearOctree& rTree) {
    std::vector<std::uint8_t> mark(rTree.Size(), 0);
    parallel_for(rTree.Size(), [&](std::size_t i) {
        if (rTree.mLevel[i] == 0) {
            mark[i] = 1;
            return;
        }
        const auto lo = rTree.Lo(i);
        const std::int64_t s = rTree.Span(i);
        std::array<std::int64_t, 3> parent{};
        for (std::size_t k = 0; k < 3; ++k)
            parent[k] = lo[k] - lo[k] % (2 * s);
        for (std::int64_t c = 0; c < 8 && !mark[i]; ++c) {
            const std::int64_t j =
                rTree.Locate(parent[0] + (c & 1) * s, parent[1] + ((c >> 1) & 1) * s,
                             parent[2] + ((c >> 2) & 1) * s);
            mark[i] = j >= 0 && rTree.mLevel[static_cast<std::size_t>(j)] > rTree.mLevel[i];
        }
    });
    return mark;
}

}  // namespace

bool linear_octree_pair(LinearOctree& rTree, OctreeAdjacency Kind, std::int64_t MaxLeaves) {
    bool changed = false;
    for (;;) {
        const std::size_t before = rTree.Size();
        linear_octree_balance(rTree, Kind, MaxLeaves);
        changed = changed || rTree.Size() != before;
        std::vector<std::uint8_t> mark = linoct_unpaired(rTree);
        bool any = false;
        for (std::size_t i = 0; i < mark.size(); ++i) {
            mark[i] = mark[i] && rTree.mLevel[i] < rTree.mMaxLevel;
            any = any || mark[i];
        }
        if (!any)
            return changed;
        linear_octree_split(rTree, mark, MaxLeaves);
        changed = true;
    }
}

bool linear_octree_is_paired(const LinearOctree& rTree, std::string* pWhy) {
    const std::vector<std::uint8_t> mark = linoct_unpaired(rTree);
    for (std::size_t i = 0; i < mark.size(); ++i)
        if (mark[i]) {
            if (pWhy) {
                const auto lo = rTree.Lo(i);
                *pWhy = "leaf " + std::to_string(i) + " at (" + std::to_string(lo[0]) + ", " +
                        std::to_string(lo[1]) + ", " + std::to_string(lo[2]) + "), level " +
                        std::to_string(rTree.mLevel[i]) +
                        (rTree.mLevel[i] == 0 ? " has no parent" : " has a split sibling");
            }
            return false;
        }
    return true;
}

std::vector<std::int64_t> linear_octree_neighbours(const LinearOctree& rTree, std::size_t Leaf,
                                                   int Dx, int Dy, int Dz) {
    std::vector<std::int64_t> out;
    const auto lo = rTree.Lo(Leaf);
    const std::int64_t s = rTree.Span(Leaf);
    const int d[3] = {Dx, Dy, Dz};
    std::int64_t b[3];
    for (std::size_t k = 0; k < 3; ++k)
        b[k] = lo[k] + d[k] * s;
    if (!rTree.Contains(b[0], b[1], b[2]))
        return out;
    const std::int64_t j = rTree.Locate(b[0], b[1], b[2]);
    if (j < 0)
        return out;
    if (rTree.mLevel[static_cast<std::size_t>(j)] <= rTree.mLevel[Leaf]) {
        out.push_back(j);
        return out;
    }
    // The same-size box on that side is subdivided: keep the leaves inside it
    // that touch the leaf's face, edge or corner.
    const std::uint64_t first = linoct_key(b[0], b[1], b[2]);
    const std::uint64_t last = first + linoct_volume(s);
    auto it = std::lower_bound(rTree.mKey.begin(), rTree.mKey.end(), first);
    for (; it != rTree.mKey.end() && *it < last; ++it) {
        const auto i = static_cast<std::size_t>(it - rTree.mKey.begin());
        const auto q = rTree.Lo(i);
        const std::int64_t qs = rTree.Span(i);
        bool touches = true;
        for (std::size_t k = 0; k < 3 && touches; ++k) {
            if (d[k] > 0)
                touches = q[k] == b[k];
            else if (d[k] < 0)
                touches = q[k] + qs == b[k] + s;
        }
        if (touches)
            out.push_back(static_cast<std::int64_t>(i));
    }
    return out;
}

bool linear_octree_check(const LinearOctree& rTree, const OctreeAdjacency* pKind,
                         std::string* pWhy) {
    auto fail = [&](const std::string& rWhat) {
        if (pWhy)
            *pWhy = rWhat;
        return false;
    };
    if (rTree.mKey.size() != rTree.mLevel.size())
        return fail("keys and levels differ in length");
    std::uint64_t covered = 0;
    for (std::size_t i = 0; i < rTree.Size(); ++i) {
        if (rTree.mLevel[i] > rTree.mMaxLevel)
            return fail("leaf " + std::to_string(i) + " is below the finest level");
        const auto lo = rTree.Lo(i);
        const std::int64_t s = rTree.Span(i);
        for (std::size_t k = 0; k < 3; ++k) {
            if (lo[k] % s != 0)
                return fail("leaf " + std::to_string(i) + " is not aligned to its span");
            if (lo[k] + s > rTree.Extent(k))
                return fail("leaf " + std::to_string(i) + " leaves the forest");
        }
        if (i + 1 < rTree.Size()) {
            if (rTree.mKey[i] >= rTree.mKey[i + 1])
                return fail("keys are not strictly ascending at leaf " + std::to_string(i));
            if (rTree.mKey[i] + linoct_volume(s) > rTree.mKey[i + 1])
                return fail("leaf " + std::to_string(i) + " overlaps the next one");
        }
        covered += linoct_volume(s);
    }
    const std::uint64_t forest = static_cast<std::uint64_t>(rTree.Extent(0)) *
                                 static_cast<std::uint64_t>(rTree.Extent(1)) *
                                 static_cast<std::uint64_t>(rTree.Extent(2));
    if (covered != forest)
        return fail("the leaves cover " + std::to_string(covered) + " of " +
                    std::to_string(forest) + " finest cells");
    if (!pKind)
        return true;
    const auto dirs = linear_octree_directions(*pKind);
    for (std::size_t i = 0; i < rTree.Size(); ++i)
        for (const auto& d : dirs)
            for (std::int64_t j : linear_octree_neighbours(rTree, i, d[0], d[1], d[2]))
                if (rTree.mLevel[static_cast<std::size_t>(j)] > rTree.mLevel[i] + 1)
                    return fail("leaves " + std::to_string(i) + " and " + std::to_string(j) +
                                " are more than one level apart");
    return true;
}

}  // namespace detail
}  // namespace meshioplusplus
