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
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/operations/curvature.hpp"
#include "meshioplusplus/operations/feature_edges.hpp"
#include "meshioplusplus/operations/shrinkwrap.hpp"
#include "meshioplusplus/detail/planar_skin.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/tri_box.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kOctPrefix = "meshio++: octree: ";
constexpr std::int64_t kOctMaxLevel = 20;

// A skin primitive: a triangle in 3-D, a segment (a, b) in 2-D.
struct OctPrim {
    detail::Vec3 mA, mB, mC;
    /// A refinement criterion (curvature, feature edge) asks for `max_level` here.
    bool mRefine = false;
    /// The smallest target edge length of the primitive's vertices (infinity: none).
    double mSize = std::numeric_limits<double>::infinity();
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

OctreeTransition octree_transition_from_name(const std::string& rName) {
    if (rName == "hanging")
        return OctreeTransition::Hanging;
    if (rName == "polyhedra")
        return OctreeTransition::Polyhedra;
    if (rName == "templates")
        return OctreeTransition::Templates;
    throw std::invalid_argument(std::string(kOctPrefix) + "unknown transition '" + rName +
                                "' (expected one of: hanging, polyhedra, templates)");
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
    const std::int64_t surface_level =
        rOptions.mSurfaceLevel < 0 ? rOptions.mMaxLevel : rOptions.mSurfaceLevel;
    if (surface_level > rOptions.mMaxLevel)
        throw std::invalid_argument(std::string(kOctPrefix) + "surface_level " +
                                    std::to_string(surface_level) + " is above max_level " +
                                    std::to_string(rOptions.mMaxLevel));
    if (rOptions.mCurvatureThreshold < 0.0)
        throw std::invalid_argument(std::string(kOctPrefix) +
                                    "curvature_threshold must not be negative");
    if (rOptions.mFeatureAngle < 0.0 || rOptions.mFeatureAngle > 180.0)
        throw std::invalid_argument(std::string(kOctPrefix) +
                                    "feature_angle must lie in [0, 180]");
    if (rOptions.mSnap && rOptions.mFill == VoxelFill::All)
        throw std::invalid_argument(
            std::string(kOctPrefix) +
            "snap needs fill 'surface' or 'inside': with fill 'all' there is no outer boundary "
            "to move onto the skin");
    if (rOptions.mTransition != OctreeTransition::Hanging && rOptions.mBalance == OctreeBalance::None)
        throw std::invalid_argument(std::string(kOctPrefix) +
                                    "a conforming transition needs a balanced tree "
                                    "(balance 'face' or 'full')");

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
    // 2-D runs in the frame where the skin's plane is xy; 3-D has the identity.
    detail::PlanarFrame frame;
    if (two_d)
        frame = detail::planar_frame(rMesh, rOptions.mPlane);
    std::vector<OctPrim> prims;
    detail::PlanarSkin skin;
    detail::TriangleSoup soup;
    if (two_d) {
        if (rOptions.mCurvatureThreshold > 0.0 || rOptions.mFeatureAngle > 0.0)
            throw std::invalid_argument(std::string(kOctPrefix) +
                                        "curvature_threshold and feature_angle need a surface "
                                        "skin, but the skin is a polyline");
        skin = detail::collect_planar_skin(rMesh, kOctPrefix, frame);
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

    // Refinement criteria beyond "the skin cuts the leaf".
    std::vector<std::array<std::int64_t, 2>> feature_segments;
    if (rOptions.mCurvatureThreshold > 0.0) {
        CurvatureOptions copts;
        copts.mGaussian = false;
        copts.mRegion = rOptions.mDistance.mSurfaceRegion;
        const CurvatureResult cres = compute_curvature(rMesh, copts);
        const NDArray& h = cres.mMesh.PointData("curvature:mean");
        for (std::size_t t = 0; t < soup.NumTriangles(); ++t)
            for (std::size_t i = 0; i < 3; ++i) {
                const double v = std::fabs(
                    detail::read_double(h, static_cast<std::size_t>(soup.mVertices[t][i])));
                if (v > rOptions.mCurvatureThreshold)  // NaN compares false
                    prims[t].mRefine = true;
            }
    }
    if (rOptions.mFeatureAngle > 0.0) {
        FeatureEdgeOptions fopts;
        fopts.mFeatureAngleDeg = rOptions.mFeatureAngle;
        fopts.mInconsistent = false;
        fopts.mRegion = rOptions.mDistance.mSurfaceRegion;
        const FeatureEdgeResult fres = feature_edges(rMesh, fopts);
        // The edges themselves become refining primitives: flagging whole
        // triangles would refine the middle of a big flat face whose boundary is
        // sharp. A degenerate triangle (a, b, b) is the segment a-b for the SAT.
        std::set<std::array<std::int64_t, 2>> sharp;
        for (const auto fb : fres.mMesh.CellRange()) {
            if (fb.NodesPerCell() != 2)
                continue;
            for (std::size_t c = 0; c < fb.NumCells(); ++c) {
                const std::int64_t a = detail::read_int(fb.Conn(), c * 2);
                const std::int64_t b = detail::read_int(fb.Conn(), c * 2 + 1);
                sharp.insert({std::min(a, b), std::max(a, b)});
            }
        }
        feature_segments.assign(sharp.begin(), sharp.end());
    }
    if (!rOptions.mSizeArray.empty()) {
        if (!rMesh.HasPointData(rOptions.mSizeArray))
            throw std::invalid_argument(std::string(kOctPrefix) + "size_array '" +
                                        rOptions.mSizeArray + "' is not a point_data array");
        const NDArray& sz = rMesh.PointData(rOptions.mSizeArray);
        if (sz.Size() != rMesh.NumPoints())
            throw std::invalid_argument(std::string(kOctPrefix) + "size_array '" +
                                        rOptions.mSizeArray + "' must have one value per point");
        for (std::size_t t = 0; t < prims.size(); ++t) {
            double smallest = std::numeric_limits<double>::infinity();
            const std::size_t nv = two_d ? 2 : 3;
            for (std::size_t i = 0; i < nv; ++i) {
                const std::int64_t id = two_d ? skin.mIds[t][i] : soup.mVertices[t][i];
                const double v = detail::read_double(sz, static_cast<std::size_t>(id));
                if (v > 0.0 && v < smallest)
                    smallest = v;
            }
            prims[t].mSize = smallest;
        }
    }

    for (const auto& seg : feature_segments) {
        OctPrim e{soup.mPoints[static_cast<std::size_t>(seg[0])],
                  soup.mPoints[static_cast<std::size_t>(seg[1])],
                  soup.mPoints[static_cast<std::size_t>(seg[1])]};
        e.mRefine = true;
        prims.push_back(e);
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
    if (two_d) {  // world -> the skin's frame
        const detail::Vec3 flo = frame.Permute(detail::Vec3{{lo[0], lo[1], lo[2]}});
        const detail::Vec3 fhi = frame.Permute(detail::Vec3{{hi[0], hi[1], hi[2]}});
        lo = {flo[0], flo[1], flo[2]};
        hi = {fhi[0], fhi[1], fhi[2]};
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
            bool split = level < rOptions.mMinLevel;
            if (!split && !leaf.mCand.empty()) {
                if (level < surface_level) {
                    split = true;
                } else {
                    const double edge = grid.mRootSize / static_cast<double>(std::int64_t{1} << level);
                    for (std::int32_t p : leaf.mCand) {
                        const OctPrim& prim = prims[static_cast<std::size_t>(p)];
                        if (prim.mRefine || edge > prim.mSize) {
                            split = true;
                            break;
                        }
                    }
                }
            }
            if (split)
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

    // Coordinates are held in the skin's frame until the end.
    std::vector<double> coords(npoints * 3);
    parallel_for(npoints, [&](std::size_t p) {
        const std::int64_t key = uniq[p];
        const std::int64_t x = key % px;
        const std::int64_t y = (key / px) % py;
        const std::int64_t z = key / (px * py);
        coords[p * 3 + 0] = grid.mOrigin[0] + static_cast<double>(x) * grid.mFine;
        coords[p * 3 + 1] = grid.mOrigin[1] + static_cast<double>(y) * grid.mFine;
        coords[p * 3 + 2] =
            two_d ? grid.mOrigin[2] : grid.mOrigin[2] + static_cast<double>(z) * grid.mFine;
    });

    std::vector<std::int64_t> connv(kept.size() * nc);
    parallel_for(kept.size() * nc, [&](std::size_t i) {
        connv[i] = static_cast<std::int64_t>(
            std::lower_bound(uniq.begin(), uniq.end(), keys[i]) - uniq.begin());
    });

    // ---- leaf lookup over the kept leaves ----------------------------------------
    std::vector<std::unordered_map<std::uint64_t, std::size_t>> kept_maps(
        static_cast<std::size_t>(grid.mFinest + 1));
    for (std::size_t n = 0; n < kept.size(); ++n) {
        const OctLeaf& leaf = leaves[kept[n]];
        kept_maps[static_cast<std::size_t>(leaf.mLevel)][oct_morton(leaf.mX, leaf.mY, leaf.mZ)] = n;
    }
    // The kept leaf containing the finest cell @p rC, or -1.
    const auto leaf_containing = [&](const std::array<std::int64_t, 3>& rC) -> std::int64_t {
        for (std::size_t k = 0; k < 3; ++k)
            if (rC[k] < 0 || rC[k] >= grid.Extent(k))
                return -1;
        for (std::int64_t lev = grid.mFinest; lev >= 0; --lev) {
            const std::int64_t sh = grid.mFinest - lev;
            const auto it = kept_maps[static_cast<std::size_t>(lev)].find(
                oct_morton((rC[0] >> sh) << sh, (rC[1] >> sh) << sh, (rC[2] >> sh) << sh));
            if (it != kept_maps[static_cast<std::size_t>(lev)].end())
                return static_cast<std::int64_t>(it->second);
        }
        return -1;
    };

    // ---- hanging nodes: relative to the kept leaves ----------------------------
    std::vector<std::int64_t> hang(npoints, 0);
    parallel_for(npoints, [&](std::size_t p) {
        const std::int64_t key = uniq[p];
        const std::int64_t pt[3] = {key % px, (key / px) % py, key / (px * py)};
        std::int64_t flag = 0;
        const int nzo = two_d ? 1 : 2;
        for (int oz = 0; oz < nzo && !flag; ++oz)
            for (int oy = 0; oy < 2 && !flag; ++oy)
                for (int ox = 0; ox < 2 && !flag; ++ox) {
                    const std::array<std::int64_t, 3> c{pt[0] - 1 + ox, pt[1] - 1 + oy,
                                                        two_d ? 0 : pt[2] - 1 + oz};
                    for (std::int64_t lev = grid.mFinest; lev >= 0; --lev) {
                        bool inside = true;
                        for (std::size_t k = 0; k < 3; ++k)
                            if (c[k] < 0 || c[k] >= grid.Extent(k))
                                inside = false;
                        if (!inside)
                            break;
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
        hang[p] = flag;
    });

    // ---- snap: move the outer nodes onto the skin -------------------------------
    if (rOptions.mSnap) {
        // Corner indices of the face (edge, in 2-D) on [axis][side] of a leaf.
        static const int kFace3[3][2][4] = {{{0, 4, 7, 3}, {1, 2, 6, 5}},
                                            {{0, 1, 5, 4}, {3, 7, 6, 2}},
                                            {{0, 3, 2, 1}, {4, 5, 6, 7}}};
        static const int kEdge2[2][2][2] = {{{3, 0}, {1, 2}}, {{0, 1}, {2, 3}}};
        std::vector<char> outer(npoints, 0);
        const std::size_t naxes = static_cast<std::size_t>(grid.NumAxes());
        for (std::size_t n = 0; n < kept.size(); ++n) {
            const OctLeaf& leaf = leaves[kept[n]];
            const std::int64_t s = std::int64_t{1} << (grid.mFinest - leaf.mLevel);
            const std::int64_t lo3[3] = {leaf.mX, leaf.mY, leaf.mZ};
            for (std::size_t a = 0; a < naxes; ++a)
                for (std::int64_t sd = 0; sd < 2; ++sd) {
                    // Probe the finest cells on the outside of the face around its
                    // centre: any kept leaf there means the face is shared.
                    bool covered = false;
                    const std::size_t u = two_d ? 1 - a : (a + 1) % 3, v = (a + 2) % 3;
                    const std::int64_t offs[2] = {s / 2 - 1 < 0 ? 0 : s / 2 - 1, s / 2};
                    for (std::int64_t ou : offs)
                        for (std::int64_t ov : offs) {
                            std::array<std::int64_t, 3> c{lo3[0], lo3[1], lo3[2]};
                            c[a] = sd ? lo3[a] + s : lo3[a] - 1;
                            c[u] = lo3[u] + ou;
                            if (!two_d)
                                c[v] = lo3[v] + ov;
                            if (two_d)
                                c[2] = 0;
                            if (leaf_containing(c) >= 0)
                                covered = true;
                        }
                    if (covered)
                        continue;
                    if (two_d) {
                        for (int e = 0; e < 2; ++e)
                            outer[static_cast<std::size_t>(
                                connv[n * 4 + static_cast<std::size_t>(kEdge2[a][sd][e])])] = 1;
                    } else {
                        for (int e = 0; e < 4; ++e)
                            outer[static_cast<std::size_t>(
                                connv[n * 8 + static_cast<std::size_t>(kFace3[a][sd][e])])] = 1;
                    }
                }
        }
        for (std::size_t p = 0; p < npoints; ++p)
            if (hang[p])
                outer[p] = 0;  // a hanging node must stay on its parent edge/face

        std::vector<double> target = coords;
        if (two_d) {
            parallel_for(npoints, [&](std::size_t p) {
                if (!outer[p])
                    return;
                const detail::Vec3 q = detail::polyline_closest_point_2d(
                    skin, detail::Vec3{{coords[p * 3], coords[p * 3 + 1], coords[p * 3 + 2]}});
                for (std::size_t k = 0; k < 3; ++k)
                    target[p * 3 + k] = q[k];
            });
        } else {
            Mesh work;
            NDArray wpts = NDArray::Uninit(DType::Float64, {npoints, std::size_t{3}});
            std::memcpy(wpts.Data(), coords.data(), coords.size() * sizeof(double));
            work.AssignPoints(std::move(wpts));
            NDArray wflag = NDArray::Uninit(DType::Int64, {npoints});
            for (std::size_t p = 0; p < npoints; ++p)
                wflag.As<std::int64_t>()[p] = outer[p];
            work.AddPointData("octree:snap", std::move(wflag));
            ShrinkwrapOptions sopts;
            sopts.mWeights = "octree:snap";
            sopts.mTargetRegion = rOptions.mDistance.mSurfaceRegion;
            const ShrinkwrapResult sres = shrinkwrap(work, rMesh, sopts);
            const double* moved = sres.mMesh.Points().As<double>();
            for (std::size_t p = 0; p < npoints; ++p)
                if (outer[p])
                    for (std::size_t k = 0; k < 3; ++k)
                        target[p * 3 + k] = moved[p * 3 + k];
        }

        // Guard: a leaf with a non-positive corner Jacobian gets every node it
        // owns put back; repeat until none is left. Reverting only ever restores
        // the lattice, whose leaves are valid, so this terminates.
        static const int kAdj3[8][3] = {{1, 3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2, 7},
                                        {7, 5, 0}, {4, 6, 1}, {5, 7, 2}, {6, 4, 3}};
        std::vector<double> cur = target;
        std::vector<char> moved_flag(npoints, 0), refused(npoints, 0);
        for (std::size_t p = 0; p < npoints; ++p)
            for (std::size_t k = 0; k < 3; ++k)
                if (cur[p * 3 + k] != coords[p * 3 + k])
                    moved_flag[p] = 1;
        const auto bad_cell = [&](std::size_t n) {
            for (std::size_t c = 0; c < nc; ++c) {
                const double* p0 = &cur[static_cast<std::size_t>(connv[n * nc + c]) * 3];
                if (two_d) {
                    const double* pn = &cur[static_cast<std::size_t>(connv[n * 4 + (c + 1) % 4]) * 3];
                    const double* pp = &cur[static_cast<std::size_t>(connv[n * 4 + (c + 3) % 4]) * 3];
                    const double det = (pn[0] - p0[0]) * (pp[1] - p0[1]) -
                                       (pn[1] - p0[1]) * (pp[0] - p0[0]);
                    if (!(det > 0.0))
                        return true;
                } else {
                    double e[3][3];
                    for (int j = 0; j < 3; ++j) {
                        const double* pj =
                            &cur[static_cast<std::size_t>(connv[n * 8 + static_cast<std::size_t>(kAdj3[c][j])]) * 3];
                        for (std::size_t k = 0; k < 3; ++k)
                            e[j][k] = pj[k] - p0[k];
                    }
                    const double det = e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1]) -
                                       e[0][1] * (e[1][0] * e[2][2] - e[1][2] * e[2][0]) +
                                       e[0][2] * (e[1][0] * e[2][1] - e[1][1] * e[2][0]);
                    if (!(det > 0.0))
                        return true;
                }
            }
            return false;
        };
        for (std::size_t round = 0; round <= npoints; ++round) {
            bool any = false;
            for (std::size_t n = 0; n < kept.size(); ++n) {
                if (!bad_cell(n))
                    continue;
                for (std::size_t c = 0; c < nc; ++c) {
                    const std::size_t p = static_cast<std::size_t>(connv[n * nc + c]);
                    if (!moved_flag[p])
                        continue;
                    for (std::size_t k = 0; k < 3; ++k)
                        cur[p * 3 + k] = coords[p * 3 + k];
                    moved_flag[p] = 0;
                    refused[p] = 1;
                    any = true;
                }
            }
            if (!any)
                break;
        }
        coords = cur;
        for (std::size_t p = 0; p < npoints; ++p) {
            result.mNumSnapped += moved_flag[p] ? 1 : 0;
            result.mNumRefused += refused[p] ? 1 : 0;
        }
    }

    // ---- cell data every output carries ----------------------------------------
    // One value per kept leaf; a conforming transition repeats a leaf's value on
    // each cell cut from it.
    const auto leaf_values = [&](const std::vector<std::size_t>& rParents) {
        struct Out {
            NDArray mLevel, mCut, mDist, mLeaf;
        } out;
        out.mLevel = NDArray::Uninit(DType::Int64, {rParents.size()});
        out.mCut = NDArray::Uninit(DType::Int64, {rParents.size()});
        out.mLeaf = NDArray::Uninit(DType::Int64, {rParents.size()});
        if (rOptions.mAttachDistance)
            out.mDist = NDArray::Uninit(DType::Float64, {rParents.size()});
        for (std::size_t i = 0; i < rParents.size(); ++i) {
            const std::size_t n = rParents[i];
            out.mLevel.As<std::int64_t>()[i] = leaves[kept[n]].mLevel;
            out.mCut.As<std::int64_t>()[i] = leaves[kept[n]].mCut ? 1 : 0;
            out.mLeaf.As<std::int64_t>()[i] = static_cast<std::int64_t>(n);
            if (rOptions.mAttachDistance)
                out.mDist.As<double>()[i] = distance[kept[n]];
        }
        return out;
    };
    std::vector<std::vector<std::size_t>> block_parents;  // per cell block, in order
    const auto finish_cell_data = [&]() {
        std::vector<NDArray> lv, ct, ds, lf;
        for (const auto& parents : block_parents) {
            auto vals = leaf_values(parents);
            lv.push_back(std::move(vals.mLevel));
            ct.push_back(std::move(vals.mCut));
            lf.push_back(std::move(vals.mLeaf));
            if (rOptions.mAttachDistance)
                ds.push_back(std::move(vals.mDist));
        }
        result.mMesh.AddCellData(kOctreeLevelName, std::move(lv));
        result.mMesh.AddCellData(kOctreeCutName, std::move(ct));
        if (rOptions.mAttachDistance)
            result.mMesh.AddCellData(kSdfDistanceName, std::move(ds));
        if (rOptions.mTransition != OctreeTransition::Hanging)
            result.mMesh.AddCellData(kOctreeLeafName, std::move(lf));
    };
    // Frame -> world, applied to the final point array.
    const auto emit_points = [&](std::vector<double>& rC) {
        const std::size_t np = rC.size() / 3;
        NDArray out = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
        double* dst = out.As<double>();
        for (std::size_t p = 0; p < np; ++p) {
            const detail::Vec3 w =
                frame.Unpermute(detail::Vec3{{rC[p * 3], rC[p * 3 + 1], rC[p * 3 + 2]}});
            dst[p * 3 + 0] = w[0];
            dst[p * 3 + 1] = w[1];
            dst[p * 3 + 2] = w[2];
        }
        result.mMesh.AssignPoints(std::move(out));
    };
    const auto matrix = [&](DType Type, std::size_t Rows, std::size_t Cols,
                            const std::vector<std::int64_t>& rValues) {
        NDArray out = NDArray::Uninit(Type, {Rows, Cols});
        if (!rValues.empty())
            std::memcpy(out.Data(), rValues.data(), rValues.size() * sizeof(std::int64_t));
        return out;
    };

    if (rOptions.mTransition == OctreeTransition::Hanging) {
        for (std::size_t p = 0; p < npoints; ++p)
            result.mNumHanging += hang[p];
        emit_points(coords);
        result.mMesh.AddCellBlock(two_d ? "quad" : "hexahedron",
                                  matrix(DType::Int64, kept.size(), nc, connv));
        std::vector<std::size_t> all(kept.size());
        for (std::size_t n = 0; n < all.size(); ++n)
            all[n] = n;
        block_parents.push_back(std::move(all));
        finish_cell_data();
        NDArray hanging = NDArray::Uninit(DType::Int64, {npoints});
        for (std::size_t p = 0; p < npoints; ++p)
            hanging.As<std::int64_t>()[p] = hang[p];
        result.mMesh.AddPointData(kOctreeHangingName, std::move(hanging));
        return result;
    }

    // ---- conforming transitions ---------------------------------------------------
    // Everything below works from the integer lattice: which nodes exist on a leaf's
    // boundary decides its faces, and both cells sharing a boundary look at the same
    // global node set, so they agree by construction.
    std::unordered_map<std::int64_t, std::int64_t> index_of;
    index_of.reserve(npoints * 2);
    for (std::size_t p = 0; p < npoints; ++p)
        index_of[uniq[p]] = static_cast<std::int64_t>(p);
    using Lat = std::array<std::int64_t, 3>;
    const auto node_at = [&](const Lat& rL) -> std::int64_t {
        const auto it = index_of.find((rL[2] * py + rL[1]) * px + rL[0]);
        return it == index_of.end() ? -1 : it->second;
    };
    // Interior nodes of the lattice segment a -> b, in order (bisecting while the
    // midpoint exists).
    const std::function<void(const Lat&, const Lat&, std::vector<std::int64_t>&)> edge_nodes =
        [&](const Lat& rA, const Lat& rB, std::vector<std::int64_t>& rOut) {
            Lat mid{};
            std::int64_t len = 0;
            for (std::size_t k = 0; k < 3; ++k) {
                len = std::max<std::int64_t>(len, std::llabs(rB[k] - rA[k]));
                mid[k] = (rA[k] + rB[k]) / 2;
            }
            if (len < 2)
                return;
            const std::int64_t id = node_at(mid);
            if (id < 0)
                return;
            edge_nodes(rA, mid, rOut);
            rOut.push_back(id);
            edge_nodes(mid, rB, rOut);
        };
    // One face polygon: the nodes around a lattice square, plus its four corners.
    struct FacePoly {
        std::vector<std::int64_t> mNodes;
        std::array<std::int64_t, 4> mCorner{};
        std::array<std::int64_t, 3> mCentre2{};  // doubled lattice coordinates
    };
    const auto square_poly = [&](const std::array<Lat, 4>& rCorners) {
        FacePoly poly;
        for (std::size_t i = 0; i < 4; ++i) {
            const std::int64_t id = node_at(rCorners[i]);
            poly.mCorner[i] = id;
            poly.mNodes.push_back(id);
            edge_nodes(rCorners[i], rCorners[(i + 1) % 4], poly.mNodes);
        }
        for (std::size_t k = 0; k < 3; ++k)
            poly.mCentre2[k] = rCorners[0][k] + rCorners[2][k];
        return poly;
    };
    // Faces of the square of side @p Size at (u0, v0) on the plane of leaf face
    // (@p A, @p Sd): the square itself, or its four quarters if its centre exists.
    const std::function<void(const Lat&, std::size_t, std::int64_t, std::int64_t, std::int64_t,
                             std::int64_t, std::vector<FacePoly>&)>
        face_squares = [&](const Lat& rLo, std::size_t A, std::int64_t Sd, std::int64_t U0,
                           std::int64_t V0, std::int64_t Size, std::vector<FacePoly>& rOut) {
            const std::size_t u = (A + 1) % 3, v = (A + 2) % 3;
            const auto at = [&](std::int64_t du, std::int64_t dv) {
                Lat l = rLo;
                l[u] = U0 + du;
                l[v] = V0 + dv;
                return l;
            };
            // `rLo[A]` already holds the plane's coordinate.
            if (Size >= 2) {
                const Lat centre = at(Size / 2, Size / 2);
                if (node_at(centre) >= 0) {
                    const std::int64_t h = Size / 2;
                    for (std::int64_t qv = 0; qv < 2; ++qv)
                        for (std::int64_t qu = 0; qu < 2; ++qu)
                            face_squares(rLo, A, Sd, U0 + qu * h, V0 + qv * h, h, rOut);
                    return;
                }
            }
            std::array<Lat, 4> c = {at(0, 0), at(Size, 0), at(Size, Size), at(0, Size)};
            if (Sd == 0)
                std::swap(c[1], c[3]);  // outward is -axis: reverse the winding
            rOut.push_back(square_poly(c));
        };

    std::vector<std::int64_t> hex_conn;                  // untouched leaves
    std::vector<std::size_t> hex_parent;
    // polyhedra
    std::vector<std::int64_t> poly_flat, poly_rows{0}, poly_faces{0};
    std::vector<std::size_t> poly_parent;
    // templates
    std::vector<std::int64_t> pyr_conn, tet_conn, tri_conn;
    std::vector<std::size_t> pyr_parent, tet_parent, tri_parent;
    std::vector<std::array<std::int64_t, 4>> fan_corners;     // per face-centre point
    std::map<std::array<std::int64_t, 3>, std::size_t> fan_slot;  // centre2 -> fan_corners
    struct Pending {
        std::size_t mLeaf;
        std::vector<FacePoly> mFaces;
    };
    std::vector<Pending> pending;  // transition leaves of the templates path
    const bool templates = rOptions.mTransition == OctreeTransition::Templates;

    for (std::size_t n = 0; n < kept.size(); ++n) {
        const OctLeaf& leaf = leaves[kept[n]];
        const std::int64_t s = std::int64_t{1} << (grid.mFinest - leaf.mLevel);
        const Lat lo{leaf.mX, leaf.mY, leaf.mZ};
        std::vector<FacePoly> faces;
        if (two_d) {
            const std::array<Lat, 4> c = {Lat{lo[0], lo[1], 0}, Lat{lo[0] + s, lo[1], 0},
                                          Lat{lo[0] + s, lo[1] + s, 0}, Lat{lo[0], lo[1] + s, 0}};
            faces.push_back(square_poly(c));
        } else {
            for (std::size_t a = 0; a < 3; ++a)
                for (std::int64_t sd = 0; sd < 2; ++sd) {
                    Lat plane = lo;
                    plane[a] = lo[a] + sd * s;
                    face_squares(plane, a, sd, lo[(a + 1) % 3], lo[(a + 2) % 3], s, faces);
                }
        }
        bool plain = faces.size() == (two_d ? 1u : 6u);
        for (const FacePoly& f : faces)
            if (f.mNodes.size() != 4)
                plain = false;
        if (plain) {
            for (std::size_t c = 0; c < nc; ++c)
                hex_conn.push_back(connv[n * nc + c]);
            hex_parent.push_back(n);
            continue;
        }
        if (!templates) {
            if (two_d) {
                poly_flat.insert(poly_flat.end(), faces[0].mNodes.begin(), faces[0].mNodes.end());
                poly_rows.push_back(static_cast<std::int64_t>(poly_flat.size()));
            } else {
                for (const FacePoly& f : faces) {
                    poly_flat.insert(poly_flat.end(), f.mNodes.begin(), f.mNodes.end());
                    poly_rows.push_back(static_cast<std::int64_t>(poly_flat.size()));
                }
                poly_faces.push_back(static_cast<std::int64_t>(poly_rows.size()) - 1);
            }
            poly_parent.push_back(n);
        } else {
            for (const FacePoly& f : faces)
                if (f.mNodes.size() != 4 && !two_d) {
                    const auto ins = fan_slot.emplace(f.mCentre2, fan_corners.size());
                    if (ins.second)
                        fan_corners.push_back(f.mCorner);
                }
            pending.push_back(Pending{n, std::move(faces)});
        }
    }

    if (templates) {
        // New points: face centres (in key order), then leaf centroids (in leaf order).
        const std::size_t nfan = fan_slot.size();
        std::size_t next_id = npoints;
        std::map<std::array<std::int64_t, 3>, std::int64_t> fan_id;
        coords.resize((npoints + nfan + pending.size()) * 3);
        const auto mean = [&](const std::int64_t* pIds, std::size_t N, std::size_t Dest) {
            for (std::size_t k = 0; k < 3; ++k) {
                double sum = 0.0;
                for (std::size_t i = 0; i < N; ++i)
                    sum += coords[static_cast<std::size_t>(pIds[i]) * 3 + k];
                coords[Dest * 3 + k] = sum / static_cast<double>(N);
            }
        };
        for (const auto& kv : fan_slot) {
            fan_id[kv.first] = static_cast<std::int64_t>(next_id);
            mean(fan_corners[kv.second].data(), 4, next_id++);
        }
        for (const Pending& pl : pending) {
            mean(&connv[pl.mLeaf * nc], nc, next_id);
            const std::int64_t apex = static_cast<std::int64_t>(next_id++);
            if (two_d) {
                const auto& ring = pl.mFaces[0].mNodes;
                for (std::size_t i = 0; i < ring.size(); ++i) {
                    tri_conn.push_back(ring[i]);
                    tri_conn.push_back(ring[(i + 1) % ring.size()]);
                    tri_conn.push_back(apex);
                    tri_parent.push_back(pl.mLeaf);
                }
                continue;
            }
            for (const FacePoly& f : pl.mFaces) {
                if (f.mNodes.size() == 4) {
                    pyr_conn.insert(pyr_conn.end(), {f.mNodes[0], f.mNodes[3], f.mNodes[2],
                                                     f.mNodes[1], apex});
                    pyr_parent.push_back(pl.mLeaf);
                    continue;
                }
                const std::int64_t centre = fan_id.at(f.mCentre2);
                for (std::size_t i = 0; i < f.mNodes.size(); ++i) {
                    tet_conn.insert(tet_conn.end(), {centre, f.mNodes[(i + 1) % f.mNodes.size()],
                                                     f.mNodes[i], apex});
                    tet_parent.push_back(pl.mLeaf);
                }
            }
        }
    }

    emit_points(coords);
    const auto add_block = [&](const char* pType, std::size_t Width, std::vector<std::int64_t>& rConn,
                               std::vector<std::size_t>& rParents) {
        if (rParents.empty())
            return;
        result.mMesh.AddCellBlock(pType, matrix(DType::Int64, rParents.size(), Width, rConn));
        block_parents.push_back(std::move(rParents));
    };
    add_block(two_d ? "quad" : "hexahedron", nc, hex_conn, hex_parent);
    if (!templates) {
        if (!poly_parent.empty()) {
            if (two_d)
                result.mMesh.AddPolygonBlock("polygon", std::move(poly_flat), std::move(poly_rows));
            else
                result.mMesh.AddPolyhedronBlock("polyhedron", std::move(poly_flat),
                                                std::move(poly_rows), std::move(poly_faces));
            block_parents.push_back(std::move(poly_parent));
        }
    } else if (two_d) {
        add_block("triangle", 3, tri_conn, tri_parent);
    } else {
        add_block("pyramid", 5, pyr_conn, pyr_parent);
        add_block("tetra", 4, tet_conn, tet_parent);
    }
    finish_cell_data();
    return result;
}

}  // namespace meshioplusplus
