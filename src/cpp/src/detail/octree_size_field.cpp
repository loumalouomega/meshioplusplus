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
// The skin, size field and adaptive octree of detail/octree_size_field.hpp.
// Anonymous-namespace helpers are prefixed `szf_` (the amalgamation
// concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "crease_edges.hpp"
#include "octree_size_field.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/tri_box.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/curvature.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kSzfPi = 3.14159265358979323846;

bool szf_has_volume_cells(const Mesh& rMesh) {
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        if (cell_type_dimension(cell_type_from_name(rMesh.Cells(b).Type())) == 3)
            return true;
    return false;
}

// A leaf under construction.
struct SzfLeaf {
    std::int64_t mX = 0, mY = 0, mZ = 0;
    std::int32_t mLevel = 0;
    std::vector<std::int32_t> mCand;
};

}  // namespace

OctreeSkin octree_skin(const Mesh& rMesh, const std::string& rRegion, const char* pPrefix) {
    OctreeSkin skin;
    // A volume mesh contributes its boundary.
    if (szf_has_volume_cells(rMesh)) {
        skin.mExtracted = extract_surface(rMesh);
        skin.mFromVolume = true;
        skin.mSoup = build_triangle_soup(skin.mExtracted, rRegion);
    } else
        skin.mSoup = build_triangle_soup(rMesh, rRegion);
    const TriangleSoup& soup = skin.mSoup;
    const std::size_t nt = soup.NumTriangles();
    if (nt == 0)
        throw std::invalid_argument(std::string(pPrefix) +
                                    "the skin has no surface triangles (a mesh of `line` cells is "
                                    "a planar skin: use `octree`)");
    // Weld exactly coincident points: sort the referenced points by their
    // coordinates, then their id; a group's representative is its lowest id,
    // and welded ids follow the representatives' input order.
    std::vector<std::int64_t> used;
    used.reserve(nt * 3);
    for (const auto& tri : soup.mVertices)
        for (std::int64_t v : tri)
            used.push_back(v);
    std::sort(used.begin(), used.end());
    used.erase(std::unique(used.begin(), used.end()), used.end());
    std::vector<std::int64_t> by_coord = used;
    auto coord = [&](std::int64_t v) -> const Vec3& {
        return soup.mPoints[static_cast<std::size_t>(v)];
    };
    parallel_sort(by_coord.begin(), by_coord.end(), [&](std::int64_t a, std::int64_t b) {
        const Vec3 &pa = coord(a), &pb = coord(b);
        if (pa != pb)
            return pa < pb;
        return a < b;
    });
    std::vector<std::int64_t> rep(soup.mPoints.size(), -1);
    for (std::size_t i = 0; i < by_coord.size(); ++i) {
        const std::int64_t v = by_coord[i];
        if (i > 0 && coord(by_coord[i - 1]) == coord(v))
            rep[static_cast<std::size_t>(v)] = rep[static_cast<std::size_t>(by_coord[i - 1])];
        else
            rep[static_cast<std::size_t>(v)] = v;
    }
    std::vector<std::int64_t> welded(soup.mPoints.size(), -1);
    for (std::int64_t v : used)
        if (rep[static_cast<std::size_t>(v)] == v) {
            welded[static_cast<std::size_t>(v)] = static_cast<std::int64_t>(skin.mPoints.size());
            skin.mPoints.push_back(coord(v));
            skin.mSourcePoint.push_back(v);
        }
    skin.mTris.resize(nt);
    for (std::size_t t = 0; t < nt; ++t)
        for (std::size_t k = 0; k < 3; ++k)
            skin.mTris[t][k] = welded[static_cast<std::size_t>(
                rep[static_cast<std::size_t>(soup.mVertices[t][k])])];
    // Area-weighted vertex normals, accumulated in triangle order.
    skin.mNormal.assign(skin.mPoints.size(), {0.0, 0.0, 0.0});
    for (std::size_t t = 0; t < nt; ++t) {
        const Vec3* c = &soup.mCorners[3 * t];
        const Vec3 n = vec3_cross(vec3_sub(c[1], c[0]), vec3_sub(c[2], c[0]));
        for (std::int64_t v : skin.mTris[t])
            skin.mNormal[static_cast<std::size_t>(v)] =
                vec3_add(skin.mNormal[static_cast<std::size_t>(v)], n);
    }
    for (auto& n : skin.mNormal)
        n = vec3_normalize(n);
    skin.mLo = soup.mCorners[0];
    skin.mHi = soup.mCorners[0];
    for (const Vec3& c : soup.mCorners)
        for (std::size_t k = 0; k < 3; ++k) {
            skin.mLo[k] = std::min(skin.mLo[k], c[k]);
            skin.mHi[k] = std::max(skin.mHi[k], c[k]);
        }
    return skin;
}

std::vector<OctreeFeatureEdge> octree_skin_features(const OctreeSkin& rSkin,
                                                    double FeatureAngleDeg) {
    std::vector<OctreeFeatureEdge> out;
    if (!(FeatureAngleDeg > 0.0))
        return out;
    const std::size_t nt = rSkin.mTris.size();
    std::vector<std::int64_t> corners(nt * 3);
    std::vector<double> normal(nt * 3);
    for (std::size_t t = 0; t < nt; ++t) {
        const Vec3* c = &rSkin.mSoup.mCorners[3 * t];
        const Vec3 n = vec3_normalize(vec3_cross(vec3_sub(c[1], c[0]), vec3_sub(c[2], c[0])));
        for (std::size_t k = 0; k < 3; ++k) {
            corners[3 * t + k] = rSkin.mTris[t][k];
            normal[3 * t + k] = n[k];
        }
    }
    for (const CreaseEdge& e : crease_edges_triangles(corners, normal, FeatureAngleDeg))
        if (e.IsCrease())
            out.push_back({e.mLo, e.mHi});
    return out;
}

std::int32_t OctreeFrame::LevelFor(double Size) const {
    if (!(Size < mRootSize))
        return 0;
    const double l = std::ceil(std::log2(mRootSize / Size) - 1e-9);
    return static_cast<std::int32_t>(std::clamp(l, 0.0, static_cast<double>(mMaxLevel)));
}

Vec3 OctreeFrame::Position(std::int64_t X, std::int64_t Y, std::int64_t Z) const {
    const double fine = Fine();
    return {mOrigin[0] + static_cast<double>(X) * fine, mOrigin[1] + static_cast<double>(Y) * fine,
            mOrigin[2] + static_cast<double>(Z) * fine};
}

OctreeFrame octree_frame(const Vec3& rLo, const Vec3& rHi, double MinSize, double MaxSize,
                         std::int64_t MaxLevelCap, const char* pPrefix) {
    if (!(MinSize > 0.0) || !(MaxSize >= MinSize) || !std::isfinite(MaxSize))
        throw std::invalid_argument(
            std::string(pPrefix) +
            "min_size and max_size must be positive with min_size <= max_size");
    OctreeFrame f;
    const double ratio = MaxSize / MinSize;
    std::int32_t levels = static_cast<std::int32_t>(std::ceil(std::log2(ratio) - 1e-9));
    levels = std::max(levels, 0);
    if (levels > 19)
        throw std::invalid_argument(std::string(pPrefix) + "max_size / min_size needs " +
                                    std::to_string(levels) + " octree levels, more than 19");
    f.mRootSize = MinSize * static_cast<double>(std::int64_t{1} << levels);
    f.mMaxLevel = levels;
    if (MaxLevelCap >= 0 && MaxLevelCap < levels)
        f.mMaxLevel = static_cast<std::int32_t>(MaxLevelCap);
    f.mMinSize = f.mRootSize / static_cast<double>(std::int64_t{1} << f.mMaxLevel);
    f.mMaxSize = MaxSize;
    for (std::size_t k = 0; k < 3; ++k) {
        const double ext = std::max(rHi[k] - rLo[k], 0.0);
        const double roots = std::max(1.0, std::ceil(ext / f.mRootSize - 1e-12));
        if (roots * static_cast<double>(std::int64_t{1} << f.mMaxLevel) >= 1048576.0)
            throw std::invalid_argument(std::string(pPrefix) +
                                        "the lattice needs 2^20 or more finest cells on an axis "
                                        "(raise min_size or lower max_level)");
        f.mRoots[k] = static_cast<std::int64_t>(roots);
        f.mOrigin[k] = 0.5 * (rLo[k] + rHi[k]) - 0.5 * roots * f.mRootSize;
    }
    return f;
}

SizeField octree_size_field(const OctreeSkin& rSkin, const TriBvh& rTree,
                            const std::vector<OctreeFeatureEdge>& rFeatures,
                            const SizeFieldSpec& rSpec) {
    SizeField out;
    const std::size_t np = rSkin.mPoints.size();
    const double inf = std::numeric_limits<double>::infinity();
    out.mPointSize.assign(np, inf);
    // Curvature: the chord subtending the angle on the largest principal curvature.
    if (rSpec.mCurvatureAngleDeg > 0.0) {
        NDArray pts = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t k = 0; k < 3; ++k)
                pts.As<double>()[i * 3 + k] = rSkin.mPoints[i][k];
        NDArray conn = NDArray::Uninit(DType::Int64, {rSkin.mTris.size(), std::size_t{3}});
        for (std::size_t t = 0; t < rSkin.mTris.size(); ++t)
            for (std::size_t k = 0; k < 3; ++k)
                conn.As<std::int64_t>()[t * 3 + k] = rSkin.mTris[t][k];
        Mesh tri;
        tri.AssignPoints(std::move(pts));
        tri.AddCellBlock("triangle", std::move(conn));
        CurvatureOptions co;
        co.mRecordPrincipal = true;
        const CurvatureResult cr = compute_curvature(tri, co);
        const double* k12 = cr.mMesh.PointData(kCurvaturePrincipalName).As<double>();
        const double chord = 2.0 * std::sin(0.5 * rSpec.mCurvatureAngleDeg * kSzfPi / 180.0);
        for (std::size_t i = 0; i < np; ++i) {
            const double kappa = std::max(std::abs(k12[2 * i]), std::abs(k12[2 * i + 1]));
            if (std::isfinite(kappa) && kappa > 0.0)
                out.mPointSize[i] = std::min(out.mPointSize[i], chord / kappa);
        }
    }
    // Thickness: an inward ray that leaves through the far wall.
    if (rSpec.mGapCells > 0.0) {
        const double tmin = 1e-3 * rSpec.mMinSize;
        const double tmax = std::numeric_limits<double>::infinity();
        std::vector<double> gap(np, inf);
        parallel_for(
            np,
            [&](std::size_t i) {
                const Vec3& n = rSkin.mNormal[i];
                if (vec3_norm_sq(n) == 0.0)
                    return;
                const Vec3 in = vec3_scale(n, -1.0);
                const TriRayHit h = tri_bvh_ray(rTree, rSkin.mPoints[i], in, tmin, tmax);
                if (h.mTriangle >= 0 && vec3_dot(h.mNormal, in) > 0.0)
                    gap[i] = h.mT;
                if (rSpec.mGapOutside) {
                    const TriRayHit o = tri_bvh_ray(rTree, rSkin.mPoints[i], n, tmin, tmax);
                    if (o.mTriangle >= 0 && vec3_dot(o.mNormal, n) < 0.0)
                        gap[i] = std::min(gap[i], o.mT);
                }
            },
            256);
        for (std::size_t i = 0; i < np; ++i)
            if (std::isfinite(gap[i]))
                out.mPointSize[i] = std::min(out.mPointSize[i], gap[i] / rSpec.mGapCells);
    }
    if (rSpec.pUserSize)
        for (std::size_t i = 0; i < np; ++i) {
            const double h = rSpec.pUserSize[static_cast<std::size_t>(rSkin.mSourcePoint[i])];
            if (std::isfinite(h) && h > 0.0)
                out.mPointSize[i] = std::min(out.mPointSize[i], h);
        }
    auto clamp = [&](double h) {
        return std::clamp(std::min(h, rSpec.mMaxSize), rSpec.mMinSize, rSpec.mMaxSize);
    };
    const auto& soup = rSkin.mSoup;
    out.mPrims.resize(rSkin.mTris.size());
    out.mSkinMinSize = rSpec.mMaxSize;
    for (std::size_t t = 0; t < rSkin.mTris.size(); ++t) {
        double h = inf;
        for (std::int64_t v : rSkin.mTris[t])
            h = std::min(h, out.mPointSize[static_cast<std::size_t>(v)]);
        SizedPrim& p = out.mPrims[t];
        p.mA = soup.mCorners[3 * t];
        p.mB = soup.mCorners[3 * t + 1];
        p.mC = soup.mCorners[3 * t + 2];
        p.mSize = clamp(h);
        out.mSkinMinSize = std::min(out.mSkinMinSize, p.mSize);
    }
    if (rSpec.mFeatureSize > 0.0)
        for (const OctreeFeatureEdge& e : rFeatures) {
            SizedPrim p;
            p.mA = rSkin.mPoints[static_cast<std::size_t>(e.mA)];
            p.mB = rSkin.mPoints[static_cast<std::size_t>(e.mB)];
            p.mC = p.mB;
            p.mSize = clamp(rSpec.mFeatureSize);
            out.mSkinMinSize = std::min(out.mSkinMinSize, p.mSize);
            out.mPrims.push_back(p);
        }
    return out;
}

bool octree_prim_hits(const OctreeFrame& rFrame, const SizedPrim& rPrim,
                      const std::array<std::int64_t, 3>& rLo, std::int64_t Span) {
    const double half_fine = 0.5 * rFrame.Fine();
    Vec3 centre, half;
    for (std::size_t k = 0; k < 3; ++k) {
        centre[k] = rFrame.mOrigin[k] + static_cast<double>(2 * rLo[k] + Span) * half_fine;
        half[k] = static_cast<double>(Span) * half_fine;
    }
    return tri_box_overlap(centre, half, rPrim.mA, rPrim.mB, rPrim.mC);
}

AdaptiveOctree octree_adaptive_build(const OctreeFrame& rFrame,
                                     const std::vector<SizedPrim>& rPrims,
                                     std::int32_t SurfaceLevel, std::int64_t MaxCells) {
    const std::int32_t L = rFrame.mMaxLevel;
    std::vector<std::int32_t> prim_level(rPrims.size());
    for (std::size_t p = 0; p < rPrims.size(); ++p)
        prim_level[p] = rFrame.LevelFor(rPrims[p].mSize);
    const std::int32_t volume_level = rFrame.LevelFor(rFrame.mMaxSize);
    // Roots, each with the primitives that reach it: each primitive is tested
    // against the roots its bounding box spans only, then scattered in
    // ascending primitive order, so every list is ascending.
    std::vector<SzfLeaf> active;
    const std::int64_t root = std::int64_t{1} << L;
    const std::int64_t nx = rFrame.mRoots[0], ny = rFrame.mRoots[1], nz = rFrame.mRoots[2];
    active.resize(static_cast<std::size_t>(nx * ny * nz));
    for (std::int64_t z = 0; z < nz; ++z)
        for (std::int64_t y = 0; y < ny; ++y)
            for (std::int64_t x = 0; x < nx; ++x) {
                SzfLeaf& leaf = active[static_cast<std::size_t>((z * ny + y) * nx + x)];
                leaf.mX = x * root;
                leaf.mY = y * root;
                leaf.mZ = z * root;
            }
    std::vector<std::vector<std::int64_t>> reach(rPrims.size());
    parallel_for(
        rPrims.size(),
        [&](std::size_t p) {
            const SizedPrim& prim = rPrims[p];
            std::array<std::int64_t, 3> r0{}, r1{};
            for (std::size_t k = 0; k < 3; ++k) {
                const double a = std::min({prim.mA[k], prim.mB[k], prim.mC[k]});
                const double b = std::max({prim.mA[k], prim.mB[k], prim.mC[k]});
                // One root of slack each way: the overlap test decides.
                r0[k] = std::max<std::int64_t>(0, static_cast<std::int64_t>(std::floor(
                                                      (a - rFrame.mOrigin[k]) / rFrame.mRootSize)) -
                                                      1);
                r1[k] = std::min<std::int64_t>(rFrame.mRoots[k] - 1,
                                               static_cast<std::int64_t>(std::floor(
                                                   (b - rFrame.mOrigin[k]) / rFrame.mRootSize)) +
                                                   1);
            }
            for (std::int64_t z = r0[2]; z <= r1[2]; ++z)
                for (std::int64_t y = r0[1]; y <= r1[1]; ++y)
                    for (std::int64_t x = r0[0]; x <= r1[0]; ++x)
                        if (octree_prim_hits(rFrame, prim, {{x * root, y * root, z * root}}, root))
                            reach[p].push_back((z * ny + y) * nx + x);
        },
        64);
    for (std::size_t p = 0; p < rPrims.size(); ++p)
        for (std::int64_t r : reach[p])
            active[static_cast<std::size_t>(r)].mCand.push_back(static_cast<std::int32_t>(p));
    std::vector<SzfLeaf> done;
    for (std::int32_t level = 0; !active.empty(); ++level) {
        std::vector<std::uint8_t> split(active.size(), 0);
        parallel_for(active.size(), [&](std::size_t i) {
            const SzfLeaf& leaf = active[i];
            if (level >= L)
                return;
            std::int32_t target = volume_level;
            if (!leaf.mCand.empty()) {
                if (SurfaceLevel >= 0)
                    target = std::max(target, SurfaceLevel);
                else
                    for (std::int32_t p : leaf.mCand)
                        target = std::max(target, prim_level[static_cast<std::size_t>(p)]);
            }
            split[i] = level < target;
        });
        std::size_t nsplit = 0;
        for (std::uint8_t s : split)
            nsplit += s;
        if (MaxCells > 0 &&
            done.size() + active.size() - nsplit + 8 * nsplit > static_cast<std::size_t>(MaxCells))
            throw std::invalid_argument(
                "meshio++: octree_mesh: the octree would have more than "
                "max_cells (" +
                std::to_string(MaxCells) + ") leaves; raise min_size or max_cells");
        std::vector<std::array<SzfLeaf, 8>> kids(active.size());
        const std::int64_t half = std::int64_t{1} << (L - level - 1 >= 0 ? L - level - 1 : 0);
        parallel_for(
            active.size(),
            [&](std::size_t i) {
                if (!split[i])
                    return;
                const SzfLeaf& leaf = active[i];
                for (std::int64_t c = 0; c < 8; ++c) {
                    SzfLeaf& kid = kids[i][static_cast<std::size_t>(c)];
                    kid.mX = leaf.mX + (c & 1) * half;
                    kid.mY = leaf.mY + ((c >> 1) & 1) * half;
                    kid.mZ = leaf.mZ + ((c >> 2) & 1) * half;
                    kid.mLevel = level + 1;
                    const std::array<std::int64_t, 3> lo{{kid.mX, kid.mY, kid.mZ}};
                    for (std::int32_t p : leaf.mCand)
                        if (octree_prim_hits(rFrame, rPrims[static_cast<std::size_t>(p)], lo, half))
                            kid.mCand.push_back(p);
                }
            },
            64);
        std::vector<SzfLeaf> next;
        next.reserve(8 * nsplit);
        for (std::size_t i = 0; i < active.size(); ++i) {
            if (split[i])
                for (auto& kid : kids[i])
                    next.push_back(std::move(kid));
            else
                done.push_back(std::move(active[i]));
        }
        active.swap(next);
    }
    // Key order.
    std::vector<std::pair<std::uint64_t, std::size_t>> order(done.size());
    for (std::size_t i = 0; i < done.size(); ++i)
        order[i] = {linoct_key(done[i].mX, done[i].mY, done[i].mZ), i};
    parallel_sort(order.begin(), order.end());
    AdaptiveOctree out;
    out.mTree.mMaxLevel = L;
    out.mTree.mRoots = rFrame.mRoots;
    out.mTree.mKey.resize(done.size());
    out.mTree.mLevel.resize(done.size());
    out.mCut.resize(done.size());
    out.mCand.resize(done.size());
    for (std::size_t i = 0; i < done.size(); ++i) {
        SzfLeaf& leaf = done[order[i].second];
        out.mTree.mKey[i] = order[i].first;
        out.mTree.mLevel[i] = static_cast<std::uint8_t>(leaf.mLevel);
        out.mCut[i] = !leaf.mCand.empty();
        out.mCand[i] = std::move(leaf.mCand);
    }
    return out;
}

void octree_adaptive_refresh(AdaptiveOctree& rOctree, const OctreeFrame& rFrame,
                             const std::vector<SizedPrim>& rPrims,
                             const std::vector<std::int64_t>& rOrigin,
                             std::vector<std::vector<std::int32_t>> OldCand) {
    const LinearOctree& tree = rOctree.mTree;
    std::vector<std::vector<std::int32_t>> cand(tree.Size());
    parallel_for(
        tree.Size(),
        [&](std::size_t i) {
            const auto& from = OldCand[static_cast<std::size_t>(rOrigin[i])];
            const auto lo = tree.Lo(i);
            const std::int64_t s = tree.Span(i);
            for (std::int32_t p : from)
                if (octree_prim_hits(rFrame, rPrims[static_cast<std::size_t>(p)], lo, s))
                    cand[i].push_back(p);
        },
        64);
    rOctree.mCut.assign(tree.Size(), 0);
    for (std::size_t i = 0; i < tree.Size(); ++i)
        rOctree.mCut[i] = !cand[i].empty();
    rOctree.mCand.swap(cand);
}

}  // namespace detail
}  // namespace meshioplusplus
