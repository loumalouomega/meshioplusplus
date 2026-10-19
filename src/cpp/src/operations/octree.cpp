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
// Octree / quadtree leaves around a skin. See operations/octree.hpp for the
// contract; the numpy twin in python/meshioplusplus/_octree.py mirrors every
// step below, in the same order, so the two agree byte for byte.
//
// Anonymous-namespace helpers are prefixed `oct_` (the amalgamation concatenates
// every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/planar_skin.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/tri_box.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kOctPrefix = "meshio++: octree: ";
constexpr std::int64_t kOctMaxLevel = 20;

// A skin primitive: a triangle in 3-D, a segment (a, b) in 2-D.
struct OctPrim {
    detail::Vec3 mA, mB, mC;
};

// One leaf: its lo corner on the finest lattice, its level, and the primitives
// that reach its (possibly inflated) box. Children inherit and filter the list.
struct OctLeaf {
    std::int64_t mX = 0, mY = 0, mZ = 0;
    std::int32_t mLevel = 0;
    bool mCut = false;
    std::vector<std::int32_t> mCand;
};

// The lattice every step shares.
struct OctGrid {
    bool mTwoD = false;
    std::int64_t mFinest = 0;                   // max level
    std::array<std::int64_t, 3> mRoots{{1, 1, 1}};
    std::array<double, 3> mOrigin{{0, 0, 0}};
    double mRootSize = 0.0;
    double mFine = 0.0;      // edge of a finest cell
    double mHalfFine = 0.0;  // half of it, exact
    double mBand = 0.0;
    std::int64_t NumAxes() const { return mTwoD ? 2 : 3; }
    std::int64_t Extent(std::size_t k) const { return mRoots[k] << mFinest; }
};

std::uint64_t oct_spread(std::uint64_t v) {
    v &= 0x1fffffULL;
    v = (v | (v << 32)) & 0x1f00000000ffffULL;
    v = (v | (v << 16)) & 0x1f0000ff0000ffULL;
    v = (v | (v << 8)) & 0x100f00f00f00f00fULL;
    v = (v | (v << 4)) & 0x10c30c30c30c30c3ULL;
    v = (v | (v << 2)) & 0x1249249249249249ULL;
    return v;
}

std::uint64_t oct_morton(std::int64_t x, std::int64_t y, std::int64_t z) {
    return oct_spread(static_cast<std::uint64_t>(x)) |
           (oct_spread(static_cast<std::uint64_t>(y)) << 1) |
           (oct_spread(static_cast<std::uint64_t>(z)) << 2);
}

bool oct_prim_hits(const OctGrid& rGrid, const OctPrim& rPrim, const detail::Vec3& rCentre,
                   const detail::Vec3& rHalf) {
    if (rGrid.mTwoD)
        return detail::segment_box_overlap_2d(rCentre, rHalf, rPrim.mA, rPrim.mB);
    return detail::tri_box_overlap(rCentre, rHalf, rPrim.mA, rPrim.mB, rPrim.mC);
}

// Centre and half-extent of a leaf's box, from integer coordinates only.
void oct_box(const OctGrid& rGrid, std::int64_t X, std::int64_t Y, std::int64_t Z, std::int32_t Level,
             detail::Vec3& rCentre, detail::Vec3& rHalf) {
    const std::int64_t s = std::int64_t{1} << (rGrid.mFinest - Level);
    const std::int64_t idx[3] = {X, Y, Z};
    for (std::size_t k = 0; k < 3; ++k) {
        if (k == 2 && rGrid.mTwoD) {
            rCentre[2] = rGrid.mOrigin[2];
            rHalf[2] = 0.0;
            continue;
        }
        rCentre[k] = rGrid.mOrigin[k] + static_cast<double>(2 * idx[k] + s) * rGrid.mHalfFine;
        rHalf[k] = static_cast<double>(s) * rGrid.mHalfFine;
    }
}

OctLeaf oct_make_leaf(const OctGrid& rGrid, const std::vector<OctPrim>& rPrims,
                      const std::vector<std::int32_t>& rParentCand, std::int64_t X, std::int64_t Y,
                      std::int64_t Z, std::int32_t Level) {
    OctLeaf leaf;
    leaf.mX = X;
    leaf.mY = Y;
    leaf.mZ = Z;
    leaf.mLevel = Level;
    detail::Vec3 centre{}, half{};
    oct_box(rGrid, X, Y, Z, Level, centre, half);
    detail::Vec3 inflated = half;
    for (std::size_t k = 0; k < 3; ++k)
        inflated[k] += rGrid.mBand;
    for (std::int32_t p : rParentCand) {
        const OctPrim& prim = rPrims[static_cast<std::size_t>(p)];
        if (!oct_prim_hits(rGrid, prim, centre, inflated))
            continue;
        leaf.mCand.push_back(p);
        if (!leaf.mCut && oct_prim_hits(rGrid, prim, centre, half))
            leaf.mCut = true;
    }
    return leaf;
}

// Children of a leaf, x fastest. Order is irrelevant: the output is Morton sorted.
void oct_split(const OctGrid& rGrid, const std::vector<OctPrim>& rPrims, const OctLeaf& rLeaf,
               std::vector<OctLeaf>& rOut) {
    const std::int64_t half = std::int64_t{1} << (rGrid.mFinest - rLeaf.mLevel - 1);
    const std::int64_t nz = rGrid.mTwoD ? 1 : 2;
    for (std::int64_t cz = 0; cz < nz; ++cz)
        for (std::int64_t cy = 0; cy < 2; ++cy)
            for (std::int64_t cx = 0; cx < 2; ++cx)
                rOut.push_back(oct_make_leaf(rGrid, rPrims, rLeaf.mCand, rLeaf.mX + cx * half,
                                             rLeaf.mY + cy * half, rLeaf.mZ + cz * half,
                                             rLeaf.mLevel + 1));
}

void oct_check_budget(std::size_t Leaves, std::int64_t MaxCells) {
    if (MaxCells > 0 && static_cast<std::int64_t>(Leaves) > MaxCells)
        throw std::invalid_argument(std::string(kOctPrefix) + "the tree has " +
                                    std::to_string(Leaves) + " leaves, above the limit of " +
                                    std::to_string(MaxCells) +
                                    " (raise max_cells, lower max_level, or narrow the band)");
}

}  // namespace

OctreeBalance octree_balance_from_name(const std::string& rName) {
    if (rName == "none")
        return OctreeBalance::None;
    if (rName == "face")
        return OctreeBalance::Face;
    if (rName == "full")
        return OctreeBalance::Full;
    throw std::invalid_argument(std::string(kOctPrefix) + "unknown balance '" + rName +
                                "' (expected one of: none, face, full)");
}

OctreeResult octree(const Mesh& rMesh, const OctreeOptions& rOptions) {
    // ---- options --------------------------------------------------------------
    if (rOptions.mMinLevel < 0 || rOptions.mMaxLevel < rOptions.mMinLevel ||
        rOptions.mMaxLevel > kOctMaxLevel)
        throw std::invalid_argument(std::string(kOctPrefix) + "need 0 <= min_level <= max_level <= " +
                                    std::to_string(kOctMaxLevel) + ", got min_level " +
                                    std::to_string(rOptions.mMinLevel) + " and max_level " +
                                    std::to_string(rOptions.mMaxLevel));
    if (rOptions.mSampleDistance < 0.0)
        throw std::invalid_argument(std::string(kOctPrefix) + "sample_distance must not be negative");
    if (rOptions.mDimension != 0 && rOptions.mDimension != 2 && rOptions.mDimension != 3)
        throw std::invalid_argument(std::string(kOctPrefix) + "dimension must be 2, 3, or 0 (auto)");

    OctGrid grid;
    grid.mFinest = rOptions.mMaxLevel;
    grid.mBand = rOptions.mSampleDistance;
    grid.mRoots = rOptions.mRootCells;

    bool two_d = rOptions.mDimension == 2;
    if (rOptions.mDimension == 0)
        two_d = detail::mesh_is_polyline_only(rMesh);
    grid.mTwoD = two_d;
    if (two_d)
        grid.mRoots[2] = 1;
    for (std::size_t k = 0; k < 3; ++k)
        if (grid.mRoots[k] <= 0)
            throw std::invalid_argument(std::string(kOctPrefix) + "root_cells must be positive, got " +
                                        std::to_string(grid.mRoots[k]) + " on axis " +
                                        std::to_string(k));

    // Keys: every finest-lattice coordinate must fit 21 bits for the Morton code,
    // and the linear point key must fit 63 bits.
    for (std::size_t k = 0; k < 3; ++k)
        if ((grid.mRoots[k] << grid.mFinest) >= (std::int64_t{1} << 21))
            throw std::invalid_argument(std::string(kOctPrefix) +
                                        "the finest lattice is too fine: root_cells << max_level "
                                        "must stay below 2^21 on every axis");

    // ---- the skin -------------------------------------------------------------
    std::vector<OctPrim> prims;
    detail::PlanarSkin skin;
    detail::TriangleSoup soup;
    if (two_d) {
        skin = detail::collect_planar_skin(rMesh, kOctPrefix);
        if (rOptions.mFill == VoxelFill::Inside || rOptions.mAttachDistance)
            detail::require_closed_polyline(skin, kOctPrefix);
        prims.reserve(skin.mSegments.size());
        for (const auto& seg : skin.mSegments)
            prims.push_back(OctPrim{seg[0], seg[1], seg[1]});
    } else {
        soup = detail::build_triangle_soup(rMesh, rOptions.mDistance.mSurfaceRegion);
        prims.reserve(soup.NumTriangles());
        for (std::size_t t = 0; t < soup.NumTriangles(); ++t)
            prims.push_back(
                OctPrim{soup.mCorners[t * 3], soup.mCorners[t * 3 + 1], soup.mCorners[t * 3 + 2]});
    }

    // ---- the root lattice -----------------------------------------------------
    std::array<double, 3> lo{}, hi{};
    if (rOptions.mBounds) {
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = (*rOptions.mBounds)[k];
            hi[k] = (*rOptions.mBounds)[3 + k];
            if (!(hi[k] >= lo[k]))
                throw std::invalid_argument(std::string(kOctPrefix) + "bounds are inverted on axis " +
                                            std::to_string(k));
        }
    } else if (!detail::point_bbox(rMesh, lo, hi)) {
        throw std::invalid_argument(std::string(kOctPrefix) +
                                    "the mesh has no points, so it has no bounding box to cover "
                                    "(pass explicit bounds)");
    }
    {
        double diag2 = 0.0;
        for (std::size_t k = 0; k < 3; ++k)
            diag2 += (hi[k] - lo[k]) * (hi[k] - lo[k]);
        const double pad = rOptions.mPadding + rOptions.mPaddingRelative * std::sqrt(diag2);
        if (pad < 0.0)
            throw std::invalid_argument(std::string(kOctPrefix) + "padding is negative");
        // A planar skin keeps its plane: padding applies to the two in-plane axes.
        for (std::size_t k = 0; k < (two_d ? 2u : 3u); ++k) {
            lo[k] -= pad;
            hi[k] += pad;
        }
    }
    double root = 0.0;
    for (std::int64_t k = 0; k < grid.NumAxes(); ++k) {
        const double h = (hi[static_cast<std::size_t>(k)] - lo[static_cast<std::size_t>(k)]) /
                         static_cast<double>(grid.mRoots[static_cast<std::size_t>(k)]);
        root = h > root ? h : root;
    }
    if (!(root > 0.0))
        throw std::invalid_argument(std::string(kOctPrefix) +
                                    "the bounding box is degenerate, so no root cell can cover it "
                                    "(pass explicit bounds or padding)");
    grid.mOrigin = lo;
    grid.mRootSize = root;
    grid.mFine = root / static_cast<double>(std::int64_t{1} << grid.mFinest);
    grid.mHalfFine = grid.mFine * 0.5;

    detail::warn_regions_dropped(rMesh, "octree");

    // ---- refinement, level by level --------------------------------------------
    std::vector<std::int32_t> all_cand(prims.size());
    for (std::size_t i = 0; i < prims.size(); ++i)
        all_cand[i] = static_cast<std::int32_t>(i);

    std::vector<OctLeaf> current;
    for (std::int64_t k = 0; k < grid.mRoots[2]; ++k)
        for (std::int64_t j = 0; j < grid.mRoots[1]; ++j)
            for (std::int64_t i = 0; i < grid.mRoots[0]; ++i)
                current.push_back(oct_make_leaf(grid, prims, all_cand, i << grid.mFinest,
                                                j << grid.mFinest, k << grid.mFinest, 0));
    oct_check_budget(current.size(), rOptions.mMaxCells);

    std::vector<OctLeaf> done;
    for (std::int64_t level = 0; level < grid.mFinest; ++level) {
        std::vector<std::vector<OctLeaf>> kids(current.size());
        parallel_for(current.size(), [&](std::size_t n) {
            const OctLeaf& leaf = current[n];
            if (level < rOptions.mMinLevel || !leaf.mCand.empty())
                oct_split(grid, prims, leaf, kids[n]);
        });
        std::vector<OctLeaf> next;
        for (std::size_t n = 0; n < current.size(); ++n) {
            if (kids[n].empty()) {
                done.push_back(std::move(current[n]));
            } else {
                for (OctLeaf& kid : kids[n])
                    next.push_back(std::move(kid));
            }
            oct_check_budget(done.size() + next.size(), rOptions.mMaxCells);
        }
        current = std::move(next);
    }
    for (OctLeaf& leaf : current)
        done.push_back(std::move(leaf));
    std::vector<OctLeaf> leaves = std::move(done);

    // ---- 2:1 balance, rippled to a fixed point ---------------------------------
    if (rOptions.mBalance != OctreeBalance::None && grid.mFinest >= 2) {
        const bool full = rOptions.mBalance == OctreeBalance::Full;
        const std::int64_t nz = grid.mTwoD ? 0 : 1;
        for (;;) {
            std::vector<std::unordered_map<std::uint64_t, std::size_t>> maps(
                static_cast<std::size_t>(grid.mFinest + 1));
            for (std::size_t n = 0; n < leaves.size(); ++n)
                maps[static_cast<std::size_t>(leaves[n].mLevel)]
                    [oct_morton(leaves[n].mX, leaves[n].mY, leaves[n].mZ)] = n;
            std::vector<char> flag(leaves.size(), 0);
            parallel_for(leaves.size(), [&](std::size_t n) {
                const OctLeaf& leaf = leaves[n];
                if (leaf.mLevel < 2)
                    return;
                const std::int64_t s = std::int64_t{1} << (grid.mFinest - leaf.mLevel);
                for (std::int64_t dz = -nz; dz <= nz; ++dz)
                    for (std::int64_t dy = -1; dy <= 1; ++dy)
                        for (std::int64_t dx = -1; dx <= 1; ++dx) {
                            const int nonzero = (dx != 0) + (dy != 0) + (dz != 0);
                            if (nonzero == 0 || (!full && nonzero > 1))
                                continue;
                            const std::int64_t nx = leaf.mX + dx * s;
                            const std::int64_t ny = leaf.mY + dy * s;
                            const std::int64_t nzc = leaf.mZ + dz * s;
                            if (nx < 0 || ny < 0 || nzc < 0 || nx >= grid.Extent(0) ||
                                ny >= grid.Extent(1) || nzc >= grid.Extent(2))
                                continue;
                            for (std::int64_t lev = 0; lev <= leaf.mLevel - 2; ++lev) {
                                const std::int64_t sh = grid.mFinest - lev;
                                const auto it = maps[static_cast<std::size_t>(lev)].find(oct_morton(
                                    (nx >> sh) << sh, (ny >> sh) << sh, (nzc >> sh) << sh));
                                if (it != maps[static_cast<std::size_t>(lev)].end()) {
                                    std::atomic_ref<char>(flag[it->second])
                                        .store(1, std::memory_order_relaxed);
                                    break;
                                }
                            }
                        }
            });
            bool any = false;
            std::vector<OctLeaf> next;
            for (std::size_t n = 0; n < leaves.size(); ++n) {
                if (flag[n]) {
                    any = true;
                    oct_split(grid, prims, leaves[n], next);
                } else {
                    next.push_back(std::move(leaves[n]));
                }
                oct_check_budget(next.size(), rOptions.mMaxCells);
            }
            leaves = std::move(next);
            if (!any)
                break;
        }
    }

    // ---- Morton order ---------------------------------------------------------
    std::vector<std::pair<std::uint64_t, std::size_t>> sorted(leaves.size());
    for (std::size_t n = 0; n < leaves.size(); ++n)
        sorted[n] = {oct_morton(leaves[n].mX, leaves[n].mY, leaves[n].mZ), n};
    std::sort(sorted.begin(), sorted.end());
    {
        std::vector<OctLeaf> ordered;
        ordered.reserve(leaves.size());
        for (const auto& kv : sorted)
            ordered.push_back(std::move(leaves[kv.second]));
        leaves = std::move(ordered);
    }
    const std::size_t nleaves = leaves.size();

    // ---- fill rule and distance -----------------------------------------------
    std::vector<detail::Vec3> centres(nleaves);
    parallel_for(nleaves, [&](std::size_t n) {
        detail::Vec3 half{};
        oct_box(grid, leaves[n].mX, leaves[n].mY, leaves[n].mZ, leaves[n].mLevel, centres[n], half);
    });
    std::vector<double> distance;
    const bool need_distance = rOptions.mFill == VoxelFill::Inside || rOptions.mAttachDistance;
    if (need_distance) {
        distance.assign(nleaves, 0.0);
        if (two_d) {
            parallel_for(nleaves, [&](std::size_t n) {
                distance[n] = detail::polyline_signed_distance_2d(skin, centres[n], rOptions.mEvenOdd);
            });
        } else {
            SurfaceDistanceOptions dopts = rOptions.mDistance;
            if (rOptions.mFill == VoxelFill::Inside && dopts.mSign == SdfSign::Unsigned)
                throw std::invalid_argument(std::string(kOctPrefix) +
                                            "fill 'inside' needs a sign, but sign='unsigned' was given");
            dopts.mBand = 0.0;
            const detail::DistanceQuery query = detail::build_distance_query(soup, dopts);
            const std::vector<detail::DistanceHit> hits =
                detail::query_distances(query, centres, dopts);
            for (std::size_t n = 0; n < nleaves; ++n)
                distance[n] = hits[n].mSignedDistance;
        }
    }

    std::vector<std::size_t> kept;
    for (std::size_t n = 0; n < nleaves; ++n) {
        bool keep = true;
        if (rOptions.mFill == VoxelFill::Surface)
            keep = leaves[n].mCut;
        else if (rOptions.mFill == VoxelFill::Inside)
            keep = distance[n] < 0.0;
        if (keep)
            kept.push_back(n);
    }

    OctreeResult result;
    result.mDimension = two_d ? 2 : 3;
    result.mOrigin = grid.mOrigin;
    result.mRootSize = grid.mRootSize;
    result.mNumLeaves = static_cast<std::int64_t>(kept.size());
    for (std::size_t n : kept)
        result.mFinestLevel = std::max<std::int64_t>(result.mFinestLevel, leaves[n].mLevel);
    if (kept.empty())
        return result;

    // ---- points: integer corners, deduplicated by sorting ----------------------
    const std::size_t nc = two_d ? 4 : 8;
    const std::int64_t px = grid.Extent(0) + 1;
    const std::int64_t py = grid.Extent(1) + 1;
    if (!two_d && px > (std::int64_t{1} << 62) / py / (grid.Extent(2) + 1))
        throw std::invalid_argument(std::string(kOctPrefix) +
                                    "the finest lattice is too large for 63-bit point keys");
    static const int kCorner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                      {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    std::vector<std::int64_t> keys(kept.size() * nc);
    parallel_for(kept.size(), [&](std::size_t n) {
        const OctLeaf& leaf = leaves[kept[n]];
        const std::int64_t s = std::int64_t{1} << (grid.mFinest - leaf.mLevel);
        for (std::size_t c = 0; c < nc; ++c) {
            const std::int64_t x = leaf.mX + kCorner[c][0] * s;
            const std::int64_t y = leaf.mY + kCorner[c][1] * s;
            const std::int64_t z = leaf.mZ + kCorner[c][2] * s;
            keys[n * nc + c] = (z * py + y) * px + x;
        }
    });
    std::vector<std::int64_t> uniq = keys;
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    const std::size_t npoints = uniq.size();

    NDArray points = NDArray::Uninit(DType::Float64, {npoints, std::size_t{3}});
    double* pdst = points.As<double>();
    parallel_for(npoints, [&](std::size_t p) {
        const std::int64_t key = uniq[p];
        const std::int64_t x = key % px;
        const std::int64_t y = (key / px) % py;
        const std::int64_t z = key / (px * py);
        pdst[p * 3 + 0] = grid.mOrigin[0] + static_cast<double>(x) * grid.mFine;
        pdst[p * 3 + 1] = grid.mOrigin[1] + static_cast<double>(y) * grid.mFine;
        pdst[p * 3 + 2] = two_d ? grid.mOrigin[2] : grid.mOrigin[2] + static_cast<double>(z) * grid.mFine;
    });

    NDArray conn = NDArray::Uninit(DType::Int64, {kept.size(), nc});
    std::int64_t* cdst = conn.As<std::int64_t>();
    parallel_for(kept.size() * nc, [&](std::size_t i) {
        cdst[i] = static_cast<std::int64_t>(
            std::lower_bound(uniq.begin(), uniq.end(), keys[i]) - uniq.begin());
    });

    // ---- hanging nodes: relative to the kept leaves ----------------------------
    std::vector<std::unordered_map<std::uint64_t, std::size_t>> kept_maps(
        static_cast<std::size_t>(grid.mFinest + 1));
    for (std::size_t n = 0; n < kept.size(); ++n) {
        const OctLeaf& leaf = leaves[kept[n]];
        kept_maps[static_cast<std::size_t>(leaf.mLevel)][oct_morton(leaf.mX, leaf.mY, leaf.mZ)] = n;
    }
    NDArray hanging = NDArray::Uninit(DType::Int64, {npoints});
    std::int64_t* hdst = hanging.As<std::int64_t>();
    parallel_for(npoints, [&](std::size_t p) {
        const std::int64_t key = uniq[p];
        const std::int64_t pt[3] = {key % px, (key / px) % py, key / (px * py)};
        std::int64_t flag = 0;
        const int nzo = two_d ? 1 : 2;
        for (int oz = 0; oz < nzo && !flag; ++oz)
            for (int oy = 0; oy < 2 && !flag; ++oy)
                for (int ox = 0; ox < 2 && !flag; ++ox) {
                    const std::int64_t c[3] = {pt[0] - 1 + ox, pt[1] - 1 + oy,
                                               two_d ? 0 : pt[2] - 1 + oz};
                    bool inside = true;
                    for (std::size_t k = 0; k < 3; ++k)
                        if (c[k] < 0 || c[k] >= grid.Extent(k))
                            inside = false;
                    if (!inside)
                        continue;
                    for (std::int64_t lev = grid.mFinest; lev >= 0; --lev) {
                        const std::int64_t sh = grid.mFinest - lev;
                        const std::int64_t ax = (c[0] >> sh) << sh, ay = (c[1] >> sh) << sh,
                                           az = (c[2] >> sh) << sh;
                        const auto it = kept_maps[static_cast<std::size_t>(lev)].find(
                            oct_morton(ax, ay, az));
                        if (it == kept_maps[static_cast<std::size_t>(lev)].end())
                            continue;
                        const std::int64_t s = std::int64_t{1} << sh;
                        const std::int64_t a[3] = {ax, ay, az};
                        for (std::size_t k = 0; k < static_cast<std::size_t>(grid.NumAxes()); ++k)
                            if (pt[k] != a[k] && pt[k] != a[k] + s)
                                flag = 1;
                        break;
                    }
                }
        hdst[p] = flag;
    });
    for (std::size_t p = 0; p < npoints; ++p)
        result.mNumHanging += hdst[p];

    // ---- assemble ------------------------------------------------------------
    result.mMesh.AssignPoints(std::move(points));
    result.mMesh.AddCellBlock(two_d ? "quad" : "hexahedron", std::move(conn));
    NDArray level = NDArray::Uninit(DType::Int64, {kept.size()});
    NDArray cut = NDArray::Uninit(DType::Int64, {kept.size()});
    for (std::size_t n = 0; n < kept.size(); ++n) {
        level.As<std::int64_t>()[n] = leaves[kept[n]].mLevel;
        cut.As<std::int64_t>()[n] = leaves[kept[n]].mCut ? 1 : 0;
    }
    {
        std::vector<NDArray> b;
        b.push_back(std::move(level));
        result.mMesh.AddCellData(kOctreeLevelName, std::move(b));
    }
    {
        std::vector<NDArray> b;
        b.push_back(std::move(cut));
        result.mMesh.AddCellData(kOctreeCutName, std::move(b));
    }
    if (rOptions.mAttachDistance) {
        NDArray d = NDArray::Uninit(DType::Float64, {kept.size()});
        for (std::size_t n = 0; n < kept.size(); ++n)
            d.As<double>()[n] = distance[kept[n]];
        std::vector<NDArray> b;
        b.push_back(std::move(d));
        result.mMesh.AddCellData(kSdfDistanceName, std::move(b));
    }
    result.mMesh.AddPointData(kOctreeHangingName, std::move(hanging));
    return result;
}

}  // namespace meshioplusplus
