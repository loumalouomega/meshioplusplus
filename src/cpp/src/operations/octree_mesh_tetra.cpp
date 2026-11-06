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
// octree_mesh, tetrahedral output: the graded BCC lattice of the balanced
// octree, stuffed (or cleaved, one winding-number level at a time) to the
// skin. Anonymous-namespace helpers are prefixed `octmt_` (the amalgamation
// concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "../detail/bcc_octree_tets.hpp"
#include "../detail/feature_snap.hpp"
#include "../detail/isosurface_stuffing.hpp"
#include "../detail/mesh_from_soup.hpp"
#include "octree_mesh_internal.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The tetrahedron's facets (detail/cell_faces), each outward.
constexpr int kOctmtFacet[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};

NDArray octmt_ids(const std::vector<std::int64_t>& rIds, std::size_t Stride) {
    NDArray a = Stride == 1 ? NDArray::Uninit(DType::Int64, {rIds.size()})
                            : NDArray::Uninit(DType::Int64, {rIds.size() / Stride, Stride});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

}  // namespace

void octree_mesh_stuffed_lattice(OctreeMeshSetup& rS, const OctreeMeshOptions& rO,
                                 OctreeMeshResult& rR, OctreeMeshTets& rOut) {
    const LinearOctree& tree = rS.mOctree.mTree;
    octree_mesh_report_tree(tree, rR);
    const BccLattice lattice = bcc_octree_lattice(tree, rO.mMaxTets);
    const double half_fine = 0.5 * rS.mFrame.Fine();
    StuffMesh mesh;
    const std::size_t nv = lattice.mKey.size();
    mesh.mPoints.resize(nv);
    mesh.mCentre.resize(nv);
    mesh.mPinned.assign(nv, 0);
    for (std::size_t p = 0; p < nv; ++p) {
        for (std::size_t k = 0; k < 3; ++k)
            mesh.mPoints[p][k] =
                rS.mFrame.mOrigin[k] + static_cast<double>(lattice.mIjk2[p][k]) * half_fine;
        mesh.mCentre[p] = lattice.mType[p] == BccVertexType::Centre;
        // The forest's boundary stays where it is: a warp would move it off
        // its box.
        for (std::size_t k = 0; k < 3; ++k)
            if (lattice.mIjk2[p][k] == 0 || lattice.mIjk2[p][k] == 2 * tree.Extent(k))
                mesh.mPinned[p] = 1;
    }
    mesh.mTets = lattice.mTets;
    mesh.mMaterial.assign(mesh.mTets.size(), 0);
    mesh.mParent.resize(mesh.mTets.size());
    for (std::size_t t = 0; t < mesh.mTets.size(); ++t)
        mesh.mParent[t] = static_cast<std::int64_t>(t);
    mesh.mCut.assign(mesh.mTets.size(), 0);

    // The winding number at every lattice vertex; a vertex on the skin is 0.
    std::vector<double> w = octree_mesh_winding(rS, mesh.mPoints, rO, rR);
    std::vector<std::uint8_t> on_skin(nv, 0);
    {
        std::vector<std::size_t> near;
        for (std::size_t p = 0; p < nv; ++p)
            if (std::abs(w[p] - std::round(w[p])) > 0.25)
                near.push_back(p);
        std::vector<Vec3> pts(near.size());
        for (std::size_t k = 0; k < near.size(); ++k)
            pts[k] = mesh.mPoints[near[k]];
        const auto hits = query_closest_points(rS.mQuery, pts);
        for (std::size_t k = 0; k < near.size(); ++k)
            on_skin[near[k]] = hits[k].mDistance <= rS.mOnSkin;
    }
    // The levels to stuff: the threshold, and with cleaving every level above
    // it some vertex reaches (nested shells).
    std::int32_t levels = 1;
    if (rO.mConformity == OctreeMeshConformity::Cleaving) {
        double wmax = 0.0;
        for (double v : w)
            wmax = std::max(wmax, v);
        levels = std::max<std::int32_t>(
            1, static_cast<std::int32_t>(std::floor(wmax - rO.mWindingThreshold)) + 1);
        levels = std::min<std::int32_t>(levels, 16);
    }
    const double step = 1e-12 * rS.mDiagonal;
    std::vector<std::uint8_t> on_interface(nv, 0);
    for (std::int32_t level = 0; level < levels; ++level) {
        const double thr = rO.mWindingThreshold + static_cast<double>(level);
        std::vector<std::int8_t> sign(mesh.mPoints.size(), 1);
        for (std::size_t p = 0; p < nv; ++p) {
            if (on_interface[p])
                continue;  // on an earlier interface: below this one
            sign[p] = (on_skin[p] && std::abs(w[p] - thr) < 0.25) ? 0 : (w[p] > thr ? -1 : 1);
            if (sign[p] == 0) {
                mesh.mPinned[p] = 1;
                on_interface[p] = 1;
            }
        }
        // The interface on an edge: the first skin crossing from the inside
        // end that takes the winding number below the threshold, else (a hole)
        // bisection on the winding number.
        const std::vector<double>* pW = &w;
        StuffCut cut = [&, pW, thr](std::int64_t A, std::int64_t B) {
            const Vec3& a = mesh.mPoints[static_cast<std::size_t>(A)];
            const Vec3& b = mesh.mPoints[static_cast<std::size_t>(B)];
            const Vec3 ab = vec3_sub(b, a);
            double wa = static_cast<std::size_t>(A) < pW->size()
                            ? (*pW)[static_cast<std::size_t>(A)]
                            : thr + 0.5;
            const auto hits = tri_bvh_segment_hits(rS.mBvh, a, b);
            for (const TriSegmentHit& h : hits) {
                wa -= h.mSide * rS.mWindingSign;
                if (wa <= thr)
                    return vec3_add(a, vec3_scale(ab, h.mT));
            }
            double lo = 0.0, hi = 1.0;
            for (int it = 0; it < 40 && (hi - lo) * vec3_norm(ab) > step; ++it) {
                const double mid = 0.5 * (lo + hi);
                const double wm =
                    rS.mWindingSign * fast_winding(rS.mBvh, vec3_add(a, vec3_scale(ab, mid)), 0.0);
                (wm > thr ? lo : hi) = mid;
            }
            return vec3_add(a, vec3_scale(ab, 0.5 * (lo + hi)));
        };
        StuffOptions so;
        so.mAlphaLong = rO.mAlphaLong;
        so.mAlphaShort = rO.mAlphaShort;
        so.mMaterial = level;
        so.mMinVolume = 1e-12 * std::pow(rS.mFrame.mMinSize, 3.0);
        const StuffReport rep = isosurface_stuff(mesh, sign, cut, so);
        rR.mNumCutPoints += rep.mNumCutPoints;
        rR.mNumWarped += rep.mNumWarped;
        rR.mNumWarpRefused += rep.mNumWarpRefused;
        if (rep.mNumDegenerate > 0)
            octree_mesh_warn(
                rR, std::string(kOctreeMeshPrefix) + std::to_string(rep.mNumDegenerate) +
                        " stuffed tetrahedra are flat at interface " + std::to_string(level));
        // Warped vertices and cut points now lie on this interface.
        on_interface.resize(mesh.mPoints.size(), 1);
        for (std::size_t p = 0; p < nv; ++p)
            if (sign[p] == 0)
                on_interface[p] = 1;
        w.resize(mesh.mPoints.size(), rO.mWindingThreshold + static_cast<double>(level));
        on_skin.resize(mesh.mPoints.size(), 1);
    }
    if (rR.mNumWarpRefused > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(rR.mNumWarpRefused) +
                                 " warps were refused to keep the tetrahedra positive");
    rOut.mPoints = std::move(mesh.mPoints);
    rOut.mTets = std::move(mesh.mTets);
    rOut.mMaterial = std::move(mesh.mMaterial);
    rOut.mCut = std::move(mesh.mCut);
    rOut.mLevel.resize(rOut.mTets.size());
    rOut.mShape.resize(rOut.mTets.size());
    for (std::size_t t = 0; t < rOut.mTets.size(); ++t) {
        const auto parent = static_cast<std::size_t>(mesh.mParent[t]);
        rOut.mLevel[t] = tree.mLevel[static_cast<std::size_t>(lattice.mLeaf[parent])];
        rOut.mShape[t] = rOut.mCut[t] ? 0 : lattice.mShape[parent];
    }
    rOut.mOnSkin = std::move(on_interface);
    rOut.mNumMaterials = levels + 1;
}

double octree_mesh_local_size(const OctreeMeshSetup& rS, const Vec3& rP) {
    const LinearOctree& tree = rS.mOctree.mTree;
    const double fine = rS.mFrame.Fine();
    std::int64_t x[3];
    for (std::size_t k = 0; k < 3; ++k)
        x[k] = std::clamp<std::int64_t>(
            static_cast<std::int64_t>(std::floor((rP[k] - rS.mFrame.mOrigin[k]) / fine)), 0,
            tree.Extent(k) - 1);
    const std::int64_t j = tree.Locate(x[0], x[1], x[2]);
    return j < 0 ? rS.mFrame.mMinSize
                 : static_cast<double>(tree.Span(static_cast<std::size_t>(j))) * fine;
}

void octree_mesh_snap_features(OctreeMeshSetup& rS, const OctreeMeshOptions& rO,
                               OctreeMeshResult& rR, OctreeMeshTets& rTets) {
    rTets.mFeature.assign(rTets.mPoints.size(), 0);
    if (!rO.mPreserveFeatures || rS.mFeatures.empty())
        return;
    std::vector<std::array<std::int64_t, 2>> edges;
    for (const OctreeFeatureEdge& e : rS.mFeatures)
        edges.push_back({e.mA, e.mB});
    const FeatureGraph graph = feature_graph(rS.mSkin.mPoints, edges, rO.mCornerAngle);
    std::vector<std::uint8_t> inside(rTets.mTets.size()), guarded(rTets.mTets.size());
    for (std::size_t t = 0; t < inside.size(); ++t) {
        inside[t] = rTets.mMaterial[t] > 0;
        // The tetrahedra the output keeps; the others are about to be dropped.
        guarded[t] = inside[t] || rO.mKeepOutside;
    }
    std::vector<std::uint8_t> claimed(rTets.mPoints.size(), 0);
    const FeatureSnapReport rep = feature_snap(
        rTets.mPoints, rTets.mTets, inside, guarded, graph, rS.mSkin.mPoints,
        [&](const Vec3& rP) { return octree_mesh_local_size(rS, rP); },
        1e-12 * std::pow(rS.mFrame.mMinSize, 3.0), claimed, rTets.mFeature);
    rR.mNumCorners = rep.mNumCorners;
    rR.mNumCornersSnapped = rep.mNumCornersSnapped;
    rR.mNumFeatureChains = rep.mNumChains;
    rR.mNumFeatureChainsMissed = rep.mNumChainsMissed;
    rR.mNumFeatureMoves = rep.mNumMoved;
    rTets.mChainGaps = rep.mGaps;
    for (std::size_t p = 0; p < claimed.size(); ++p)
        if (claimed[p])
            rTets.mOnSkin[p] = 1;
    if (rep.mNumCornersSnapped < rep.mNumCorners)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) +
                                 std::to_string(rep.mNumCorners - rep.mNumCornersSnapped) + " of " +
                                 std::to_string(rep.mNumCorners) +
                                 " feature corners not recovered");
    if (rep.mNumChainsMissed > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(rep.mNumChainsMissed) +
                                 " of " + std::to_string(rep.mNumChains) +
                                 " feature chains not recovered");
}

void octree_mesh_fill_chain_gaps(OctreeMeshSetup& rS, const OctreeMeshOptions& rO,
                                 OctreeMeshResult& rR, OctreeMeshTets& rTets) {
    // The pocket of a notch holds outside tetrahedra, which would overlap the
    // filling if they were kept.
    if (rTets.mChainGaps.empty() || rO.mKeepOutside)
        return;
    using Index = std::int64_t;
    // The boundary triangles of the inside: the owner tetrahedron and the
    // point opposite the face.
    struct BoundaryTriangle {
        std::array<Index, 3> mNodes;
        std::size_t mTet;
        Index mOpposite;
    };
    std::vector<BoundaryTriangle> triangles;
    std::unordered_map<Index, std::vector<std::size_t>> incident;
    const auto build_boundary = [&]() {
        triangles.clear();
        incident.clear();
        std::vector<std::pair<std::array<Index, 3>, std::size_t>> faces;
        for (std::size_t t = 0; t < rTets.mTets.size(); ++t) {
            if (rTets.mMaterial[t] <= 0)
                continue;
            for (const auto& f : kOctmtFacet) {
                std::array<Index, 3> key{rTets.mTets[t][f[0]], rTets.mTets[t][f[1]],
                                         rTets.mTets[t][f[2]]};
                std::sort(key.begin(), key.end());
                faces.emplace_back(key, t);
            }
        }
        std::sort(faces.begin(), faces.end());
        for (std::size_t i = 0; i < faces.size();) {
            std::size_t j = i;
            while (j < faces.size() && faces[j].first == faces[i].first)
                ++j;
            if (j - i == 1) {
                const std::size_t t = faces[i].second;
                // The point of the tetrahedron which is not on the face.
                Index opposite = -1;
                for (const Index v : rTets.mTets[t])
                    if (v != faces[i].first[0] && v != faces[i].first[1] &&
                        v != faces[i].first[2])
                        opposite = v;
                for (const Index v : faces[i].first)
                    incident[v].push_back(triangles.size());
                triangles.push_back({faces[i].first, t, opposite});
            }
            i = j;
        }
    };
    const double volume_floor = 1e-12 * std::pow(rS.mFrame.mMinSize, 3.0);
    const auto point = [&](const Index I) -> const Vec3& {
        return rTets.mPoints[static_cast<std::size_t>(I)];
    };
    const auto volume = [&](const Index A, const Index B, const Index C, const Index D) {
        return vec3_dot(vec3_cross(vec3_sub(point(B), point(A)), vec3_sub(point(C), point(A))),
                        vec3_sub(point(D), point(A))) /
               6.0;
    };
    const auto is_chain = [&](const Index I) {
        return static_cast<std::size_t>(I) < rTets.mFeature.size() &&
               rTets.mFeature[static_cast<std::size_t>(I)] != 0;
    };
    // The tetrahedra filling the pyramid over the path from the apex, or none
    // when the fan does not make a valid pocket.
    const auto fill = [&](const Index C, const std::vector<Index>& rPath,
                          const std::vector<std::size_t>& rFan) {
        std::vector<std::array<Index, 4>> tets;
        const std::size_t k = rPath.size() - 1;
        const Index a = rPath.front();
        const Index b = rPath.back();
        // The polygon is planar: its points lie on the plane through the chain
        // line and the point furthest from it.
        const Vec3 line = vec3_sub(point(b), point(a));
        Vec3 normal{0.0, 0.0, 0.0};
        for (const Index p : rPath) {
            const Vec3 n = vec3_cross(line, vec3_sub(point(p), point(a)));
            if (vec3_norm(n) > vec3_norm(normal))
                normal = n;
        }
        const double normal_norm = vec3_norm(normal);
        const double tolerance = 1e-6 * octree_mesh_local_size(rS, point(a));
        if (normal_norm <= 0.0)
            return tets;
        for (const Index p : rPath)
            if (std::abs(vec3_dot(normal, vec3_sub(point(p), point(a)))) / normal_norm > tolerance)
                return tets;
        // One material over the patch.
        for (const std::size_t tri : rFan)
            if (rTets.mMaterial[triangles[tri].mTet] !=
                rTets.mMaterial[triangles[rFan.front()].mTet])
                return tets;
        // Fan of the cap from the chain node a: all tetrahedra on one side,
        // none flat.
        double sign = 0.0;
        for (std::size_t i = 1; i < k; ++i) {
            const double v = volume(C, a, rPath[i], rPath[i + 1]);
            if (std::abs(v) <= volume_floor || (sign != 0.0 && v * sign < 0.0))
                return std::vector<std::array<Index, 4>>{};
            sign = v;
            tets.push_back(v > 0.0 ? std::array<Index, 4>{C, a, rPath[i], rPath[i + 1]}
                                   : std::array<Index, 4>{C, a, rPath[i + 1], rPath[i]});
        }
        // Each patch triangle must border a new tetrahedron on the opposite
        // side from its owner.
        for (std::size_t j = 0; j < k; ++j) {
            const std::size_t i = std::min(std::max<std::size_t>(j, 1), k - 1);
            Index off = -1;
            for (const Index v : {a, rPath[i], rPath[i + 1]})
                if (v != rPath[j] && v != rPath[j + 1])
                    off = v;
            const double owner = volume(C, rPath[j], rPath[j + 1], triangles[rFan[j]].mOpposite);
            const double added = volume(C, rPath[j], rPath[j + 1], off);
            if (owner * added >= 0.0)
                return std::vector<std::array<Index, 4>>{};
        }
        return tets;
    };
    build_boundary();
    for (const auto& gap : rTets.mChainGaps) {
        const Index a = gap[0];
        const Index b = gap[1];
        // Nodes joined to a node by a boundary edge.
        const auto neighbours = [&](const Index I) {
            std::vector<Index> result;
            const auto it = incident.find(I);
            if (it != incident.end())
                for (const std::size_t tri : it->second)
                    for (const Index v : triangles[tri].mNodes)
                        if (v != I)
                            result.push_back(v);
            std::sort(result.begin(), result.end());
            result.erase(std::unique(result.begin(), result.end()), result.end());
            return result;
        };
        const std::vector<Index> near_a = neighbours(a);
        const std::vector<Index> near_b = neighbours(b);
        // A gap closed by an earlier fill has a boundary edge now.
        bool filled = std::binary_search(near_a.begin(), near_a.end(), b);
        std::vector<Index> apexes;
        std::set_intersection(near_a.begin(), near_a.end(), near_b.begin(), near_b.end(),
                              std::back_inserter(apexes));
        for (std::size_t n = 0; n < apexes.size() && !filled; ++n) {
            const Index c = apexes[n];
            // The link of the apex: per node, the nodes beside it on the fan,
            // with the triangle of each step.
            std::unordered_map<Index, std::vector<std::pair<Index, std::size_t>>> link;
            for (const std::size_t tri : incident[c]) {
                std::array<Index, 2> others{};
                std::size_t m = 0;
                for (const Index v : triangles[tri].mNodes)
                    if (v != c)
                        others[m++] = v;
                link[others[0]].push_back({others[1], tri});
                link[others[1]].push_back({others[0], tri});
            }
            // The two arcs of the fan from a to b, the shorter first.
            std::vector<std::pair<std::vector<Index>, std::vector<std::size_t>>> arcs;
            const auto la = link.find(a);
            if (la == link.end() || la->second.size() != 2)
                continue;
            for (const auto& [first, first_tri] : la->second) {
                std::vector<Index> path{a};
                std::vector<std::size_t> fan;
                Index previous = a;
                Index current = first;
                std::size_t tri = first_tri;
                bool closed = false;
                for (std::size_t step = 0; step < 64; ++step) {
                    path.push_back(current);
                    fan.push_back(tri);
                    if (current == b) {
                        closed = true;
                        break;
                    }
                    const auto lc = link.find(current);
                    if (current == a || is_chain(current) || lc == link.end() ||
                        lc->second.size() != 2)
                        break;
                    const auto& next = lc->second[0].first == previous ? lc->second[1]
                                                                       : lc->second[0];
                    previous = current;
                    current = next.first;
                    tri = next.second;
                }
                if (closed)
                    arcs.push_back({path, fan});
            }
            std::sort(arcs.begin(), arcs.end(), [](const auto& rA, const auto& rB) {
                return rA.first.size() < rB.first.size();
            });
            for (const auto& [path, fan] : arcs) {
                std::vector<std::array<Index, 4>> tets = fill(c, path, fan);
                if (tets.empty()) {
                    // The fan from the other chain node may avoid a flat tetrahedron.
                    const std::vector<Index> reversed(path.rbegin(), path.rend());
                    const std::vector<std::size_t> reversed_fan(fan.rbegin(), fan.rend());
                    tets = fill(c, reversed, reversed_fan);
                }
                if (tets.empty())
                    continue;
                const std::size_t owner = triangles[fan.front()].mTet;
                for (const auto& tet : tets) {
                    rTets.mTets.push_back(tet);
                    rTets.mMaterial.push_back(rTets.mMaterial[owner]);
                    rTets.mCut.push_back(1);
                    rTets.mLevel.push_back(rTets.mLevel[owner]);
                    rTets.mShape.push_back(0);
                }
                filled = true;
                build_boundary();
                break;
            }
        }
        if (filled)
            ++rR.mNumGapsFilled;
        else
            ++rR.mNumGapsSkipped;
    }
    if (rR.mNumGapsSkipped > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(rR.mNumGapsSkipped) +
                                 " gaps along feature chains not filled");
}

std::pair<double, double> octree_mesh_volume_and_tolerance(const std::vector<Vec3>& rPoints,
                                                           const std::array<std::int64_t, 4>& rTet,
                                                           const double ZeroVolumeFactor) {
    // Relative precision of the written coordinates (six significant digits): a volume below its
    // effect on the longest edge squared is noise.
    constexpr double kCoordinatePrecision = 1.0e-6;
    const auto point = [&](const std::int64_t I) -> const Vec3& {
        return rPoints[static_cast<std::size_t>(I)];
    };
    double longest = 0.0;
    double largest_coordinate = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (const double c : point(rTet[i]))
            largest_coordinate = std::max(largest_coordinate, std::abs(c));
        for (int j = i + 1; j < 4; ++j)
            longest = std::max(longest, vec3_norm(vec3_sub(point(rTet[j]), point(rTet[i]))));
    }
    const double volume = vec3_dot(vec3_cross(vec3_sub(point(rTet[1]), point(rTet[0])),
                                              vec3_sub(point(rTet[2]), point(rTet[0]))),
                                   vec3_sub(point(rTet[3]), point(rTet[0]))) /
                          6.0;
    return {volume, ZeroVolumeFactor * longest * longest * longest +
                        kCoordinatePrecision * largest_coordinate * longest * longest};
}

void octree_mesh_fix_volumes(OctreeMeshSetup& /*rS*/, const OctreeMeshOptions& rO,
                             OctreeMeshResult& rR, OctreeMeshTets& rTets) {
    using Index = std::int64_t;
    constexpr int kMaxCollapsePasses = 10;
    // Rank from which a node is a feature corner or chain node, which a collapse never moves.
    constexpr int kFeatureRank = 4;
    const auto point = [&](const Index I) -> const Vec3& {
        return rTets.mPoints[static_cast<std::size_t>(I)];
    };
    const auto volume_and_tolerance = [&](const std::array<Index, 4>& rT) {
        return octree_mesh_volume_and_tolerance(rTets.mPoints, rT, rO.mZeroVolumeFactor);
    };
    // Higher is preferred as the node a collapse keeps.
    const auto rank = [&](const Index I) {
        const std::size_t i = static_cast<std::size_t>(I);
        const int feature = i < rTets.mFeature.size() ? rTets.mFeature[i] : 0;
        const int on_skin = i < rTets.mOnSkin.size() ? rTets.mOnSkin[i] : 0;
        return (feature == 1 ? 3 : (feature == 2 ? 2 : 0)) * 2 + on_skin;
    };
    const auto considered = [&](const std::size_t T) {
        return rTets.mMaterial[T] > 0 || rO.mKeepOutside;
    };

    // Node-to-tetrahedron incidence in compressed form, rebuilt when needed: the tetrahedra of
    // node i are star[start[i]..start[i + 1]).
    std::vector<std::size_t> start;
    std::vector<std::size_t> star;
    const auto build_star = [&]() {
        start.assign(rTets.mPoints.size() + 1, 0);
        for (const auto& tet : rTets.mTets)
            for (const Index v : tet)
                ++start[static_cast<std::size_t>(v) + 1];
        for (std::size_t i = 0; i < rTets.mPoints.size(); ++i)
            start[i + 1] += start[i];
        star.resize(start.back());
        std::vector<std::size_t> fill(start.begin(), start.end() - 1);
        for (std::size_t t = 0; t < rTets.mTets.size(); ++t)
            for (const Index v : rTets.mTets[t])
                star[fill[static_cast<std::size_t>(v)]++] = t;
    };
    std::vector<Index> merge(rTets.mPoints.size());
    std::vector<std::uint8_t> touched(rTets.mPoints.size());
    for (int pass = 0; pass < kMaxCollapsePasses; ++pass) {
        std::fill(touched.begin(), touched.end(), std::uint8_t{0});
        bool any = false;
        for (std::size_t i = 0; i < merge.size(); ++i)
            merge[i] = static_cast<Index>(i);
        build_star();
        // The collapse of node A onto node B is valid when no tetrahedron around A, other than
        // those that contain B, ends up inverted or degenerate.
        const auto valid = [&](const Index A, const Index B) {
            // A skin node only moves along a skin edge, otherwise the collapse pinches the skin:
            // AB is a skin edge when a face ABC belongs to a single tetrahedron.
            if (rTets.mOnSkin[static_cast<std::size_t>(A)]) {
                std::map<Index, int> around_edge;
                for (std::size_t k = start[static_cast<std::size_t>(A)];
                     k < start[static_cast<std::size_t>(A) + 1]; ++k) {
                    const auto& tet = rTets.mTets[star[k]];
                    if (tet[0] == B || tet[1] == B || tet[2] == B || tet[3] == B)
                        for (const Index v : tet)
                            around_edge[v] += v != A && v != B;
                }
                const bool skin_edge =
                    std::any_of(around_edge.begin(), around_edge.end(),
                                [](const auto& rPair) { return rPair.second == 1; });
                if (!skin_edge)
                    return false;
            }
            for (std::size_t k = start[static_cast<std::size_t>(A)];
                 k < start[static_cast<std::size_t>(A) + 1]; ++k) {
                const auto& tet = rTets.mTets[star[k]];
                if (!considered(star[k]) || tet[0] == B || tet[1] == B || tet[2] == B ||
                    tet[3] == B)
                    continue;
                auto moved = tet;
                for (Index& v : moved)
                    v = v == A ? B : v;
                const auto [volume_before, tolerance_before] = volume_and_tolerance(tet);
                const auto [volume_after, tolerance_after] = volume_and_tolerance(moved);
                // A tetrahedron that was already bad only has to not get worse.
                if (volume_after <= tolerance_after && volume_after < volume_before)
                    return false;
            }
            return true;
        };
        for (std::size_t t = 0; t < rTets.mTets.size(); ++t) {
            const auto& tet = rTets.mTets[t];
            if (!considered(t) || touched[static_cast<std::size_t>(tet[0])] ||
                touched[static_cast<std::size_t>(tet[1])] ||
                touched[static_cast<std::size_t>(tet[2])] ||
                touched[static_cast<std::size_t>(tet[3])])
                continue;
            const auto [volume, tolerance] = volume_and_tolerance(tet);
            if (std::abs(volume) > tolerance)
                continue;
            // The edges as (length, A, B), with A the node that moves onto B: A never outranks
            // B and is never a feature corner or chain node.
            std::vector<std::tuple<double, Index, Index>> edges;
            for (int i = 0; i < 4; ++i)
                for (int j = i + 1; j < 4; ++j) {
                    const double length = vec3_norm(vec3_sub(point(tet[j]), point(tet[i])));
                    if (rank(tet[i]) < kFeatureRank && rank(tet[i]) <= rank(tet[j]))
                        edges.emplace_back(length, tet[i], tet[j]);
                    if (rank(tet[j]) < kFeatureRank && rank(tet[j]) <= rank(tet[i]))
                        edges.emplace_back(length, tet[j], tet[i]);
                }
            std::sort(edges.begin(), edges.end());
            for (const auto& [length, a, b] : edges) {
                if (!valid(a, b))
                    continue;
                // The nodes of every tetrahedron around A and B are frozen for this pass, as
                // their volumes were checked as they are.
                for (const Index n : {a, b})
                    for (std::size_t k = start[static_cast<std::size_t>(n)];
                         k < start[static_cast<std::size_t>(n) + 1]; ++k)
                        for (const Index v : rTets.mTets[star[k]])
                            touched[static_cast<std::size_t>(v)] = 1;
                merge[static_cast<std::size_t>(a)] = b;
                any = true;
                break;
            }
        }
        if (!any)
            break;
        // Remap every tetrahedron and drop those with a repeated node.
        std::size_t out = 0;
        for (std::size_t t = 0; t < rTets.mTets.size(); ++t) {
            std::array<Index, 4> tet;
            for (int k = 0; k < 4; ++k)
                tet[k] = merge[static_cast<std::size_t>(rTets.mTets[t][k])];
            if (tet[0] == tet[1] || tet[0] == tet[2] || tet[0] == tet[3] || tet[1] == tet[2] ||
                tet[1] == tet[3] || tet[2] == tet[3]) {
                ++rR.mNumCollapsed;
                continue;
            }
            rTets.mTets[out] = tet;
            rTets.mMaterial[out] = rTets.mMaterial[t];
            rTets.mCut[out] = rTets.mCut[t];
            rTets.mLevel[out] = rTets.mLevel[t];
            rTets.mShape[out] = rTets.mShape[t];
            ++out;
        }
        rTets.mTets.resize(out);
        rTets.mMaterial.resize(out);
        rTets.mCut.resize(out);
        rTets.mLevel.resize(out);
        rTets.mShape.resize(out);
    }

    // A cap is a flat tetrahedron (degenerate or slightly inverted) whose apex A lies over the
    // inside of the opposite face BCD. Swapping two of its nodes would turn it into a positive
    // tetrahedron outside the domain, so it is repaired instead: dropped when BCD is on the skin
    // or an interface, or else its neighbour across BCD is split at A.
    const auto compact = [&](const std::vector<std::uint8_t>& rDead) {
        std::size_t out = 0;
        for (std::size_t t = 0; t < rTets.mTets.size(); ++t) {
            if (rDead[t])
                continue;
            rTets.mTets[out] = rTets.mTets[t];
            rTets.mMaterial[out] = rTets.mMaterial[t];
            rTets.mCut[out] = rTets.mCut[t];
            rTets.mLevel[out] = rTets.mLevel[t];
            rTets.mShape[out] = rTets.mShape[t];
            ++out;
        }
        rTets.mTets.resize(out);
        rTets.mMaterial.resize(out);
        rTets.mCut.resize(out);
        rTets.mLevel.resize(out);
        rTets.mShape.resize(out);
    };
    constexpr int kMaxCapPasses = 3;
    for (int pass = 0; pass < kMaxCapPasses; ++pass) {
        // The tetrahedra of every face, by sorted node triplet.
        std::map<std::array<Index, 3>, std::vector<std::size_t>> faces_of;
        for (std::size_t t = 0; t < rTets.mTets.size(); ++t)
            for (const auto& facet : kOctmtFacet) {
                std::array<Index, 3> key{rTets.mTets[t][facet[0]], rTets.mTets[t][facet[1]],
                                         rTets.mTets[t][facet[2]]};
                std::sort(key.begin(), key.end());
                faces_of[key].push_back(t);
            }
        std::vector<std::uint8_t> dead(rTets.mTets.size(), 0);
        const std::size_t num_tets = rTets.mTets.size();
        std::size_t num_repaired = 0;
        for (std::size_t t = 0; t < num_tets; ++t) {
            if (dead[t] || !considered(t))
                continue;
            const auto [volume, tolerance] = volume_and_tolerance(rTets.mTets[t]);
            if (volume > tolerance)
                continue;
            const auto tet = rTets.mTets[t];
            for (int i = 0; i < 4; ++i) {
                const Index apex = tet[i];
                std::array<Index, 3> key{tet[(i + 1) % 4], tet[(i + 2) % 4], tet[(i + 3) % 4]};
                const Vec3& b = point(key[0]);
                const Vec3 edge_c = vec3_sub(point(key[1]), b);
                const Vec3 edge_d = vec3_sub(point(key[2]), b);
                const Vec3 normal = vec3_cross(edge_c, edge_d);
                const double area_2 = vec3_dot(normal, normal);
                if (area_2 <= 0.0)
                    continue;
                // Only flat tetrahedra: the height of the apex over the face is small against
                // the face size, otherwise it is a fold, not a cap.
                const Vec3 to_apex = vec3_sub(point(apex), b);
                const double height = vec3_dot(normal, to_apex);
                if (height * height > rO.mCapFlatness * 0.5 * area_2 * std::sqrt(area_2))
                    continue;
                // Barycentric coordinates of the projection of the apex on the face.
                const double lambda_c = vec3_dot(vec3_cross(to_apex, edge_d), normal) / area_2;
                const double lambda_d = vec3_dot(vec3_cross(edge_c, to_apex), normal) / area_2;
                if (lambda_c <= 0.0 || lambda_d <= 0.0 || lambda_c + lambda_d >= 1.0)
                    continue;
                std::sort(key.begin(), key.end());
                std::size_t neighbour = t;
                for (const std::size_t other : faces_of[key])
                    if (other != t)
                        neighbour = other;
                if (neighbour == t || rTets.mMaterial[neighbour] != rTets.mMaterial[t]) {
                    // Skin or interface face.
                    dead[t] = 1;
                    rTets.mOnSkin[static_cast<std::size_t>(apex)] = 1;
                    ++rR.mNumCapsRemoved;
                    ++num_repaired;
                    break;
                }
                if (dead[neighbour])
                    continue;
                // Replacing each node of the shared face by the apex splits the neighbour in
                // three, with the orientation kept.
                std::array<std::array<Index, 4>, 3> children;
                bool positive = true;
                int child = 0;
                for (const Index face_node : key) {
                    children[child] = rTets.mTets[neighbour];
                    for (Index& v : children[child])
                        v = v == face_node ? apex : v;
                    positive = positive && volume_and_tolerance(children[child]).first > 0.0;
                    ++child;
                }
                if (!positive)
                    continue;
                dead[t] = 1;
                dead[neighbour] = 1;
                for (const auto& r_child : children) {
                    rTets.mTets.push_back(r_child);
                    rTets.mMaterial.push_back(rTets.mMaterial[neighbour]);
                    rTets.mCut.push_back(1);
                    rTets.mLevel.push_back(rTets.mLevel[neighbour]);
                    rTets.mShape.push_back(0);
                }
                dead.resize(rTets.mTets.size(), 0);
                ++rR.mNumSplit;
                ++num_repaired;
                break;
            }
        }
        if (num_repaired == 0)
            break;
        compact(dead);
    }

    // A sliver is a flat tetrahedron whose opposite edges AB and CD cross. A node P is inserted
    // at the crossing and the shells of both edges are split at it, so the neighbours stay
    // conforming and the sliver, whose four children have no volume, disappears.
    constexpr int kMaxSliverPasses = 3;
    // The crossing has to lie inside both edges, and the edges have to be this close, relative
    // to the shorter one, to count as crossing.
    for (int pass = 0; pass < kMaxSliverPasses; ++pass) {
        build_star();
        std::fill(touched.begin(), touched.end(), std::uint8_t{0});
        touched.resize(rTets.mPoints.size(), 0);
        std::vector<std::uint8_t> dead(rTets.mTets.size(), 0);
        const std::size_t num_tets = rTets.mTets.size();
        std::size_t num_repaired = 0;
        const auto shell = [&](const Index A, const Index B) {
            std::vector<std::size_t> tets;
            for (std::size_t k = start[static_cast<std::size_t>(A)];
                 k < start[static_cast<std::size_t>(A) + 1]; ++k) {
                const auto& tet = rTets.mTets[star[k]];
                if (tet[0] == B || tet[1] == B || tet[2] == B || tet[3] == B)
                    tets.push_back(star[k]);
            }
            return tets;
        };
        for (std::size_t t = 0; t < num_tets; ++t) {
            const auto tet = rTets.mTets[t];
            if (dead[t] || !considered(t) || touched[static_cast<std::size_t>(tet[0])] ||
                touched[static_cast<std::size_t>(tet[1])] ||
                touched[static_cast<std::size_t>(tet[2])] ||
                touched[static_cast<std::size_t>(tet[3])])
                continue;
            const auto [volume, tolerance] = volume_and_tolerance(tet);
            if (volume > tolerance)
                continue;
            // The three pairs of opposite edges.
            for (const auto& [i, j, k, l] :
                 {std::array<int, 4>{0, 1, 2, 3}, std::array<int, 4>{0, 2, 1, 3},
                  std::array<int, 4>{0, 3, 1, 2}}) {
                const Vec3 u = vec3_sub(point(tet[j]), point(tet[i]));
                const Vec3 v = vec3_sub(point(tet[l]), point(tet[k]));
                const Vec3 w = vec3_sub(point(tet[i]), point(tet[k]));
                const double uu = vec3_dot(u, u);
                const double uv = vec3_dot(u, v);
                const double vv = vec3_dot(v, v);
                const double denominator = uu * vv - uv * uv;
                if (denominator <= 1.0e-12 * uu * vv)
                    continue;
                const double s = (uv * vec3_dot(v, w) - vv * vec3_dot(u, w)) / denominator;
                const double r = (uu * vec3_dot(v, w) - uv * vec3_dot(u, w)) / denominator;
                if (s <= rO.mCrossingMargin || s >= 1.0 - rO.mCrossingMargin ||
                    r <= rO.mCrossingMargin || r >= 1.0 - rO.mCrossingMargin)
                    continue;
                const Vec3 on_ab = vec3_add(point(tet[i]), vec3_scale(u, s));
                const Vec3 on_cd = vec3_add(point(tet[k]), vec3_scale(v, r));
                if (vec3_norm(vec3_sub(on_ab, on_cd)) >
                    rO.mCrossingGap * std::sqrt(std::min(uu, vv)))
                    continue;
                // Split each other tetrahedron of the two shells: replacing one end of the edge
                // by P keeps the orientation.
                const Index p = static_cast<Index>(rTets.mPoints.size());
                std::vector<std::pair<std::size_t, std::array<std::array<Index, 4>, 2>>> splits;
                for (const auto& [e0, e1] :
                     {std::pair<Index, Index>{tet[i], tet[j]},
                      std::pair<Index, Index>{tet[k], tet[l]}})
                    for (const std::size_t other : shell(e0, e1)) {
                        if (other == t)
                            continue;
                        std::array<std::array<Index, 4>, 2> children{rTets.mTets[other],
                                                                     rTets.mTets[other]};
                        for (Index& n : children[0])
                            n = n == e1 ? p : n;
                        for (Index& n : children[1])
                            n = n == e0 ? p : n;
                        splits.push_back({other, children});
                    }
                rTets.mPoints.push_back(vec3_scale(vec3_add(on_ab, on_cd), 0.5));
                bool splittable = true;
                for (const auto& [other, children] : splits)
                    for (const auto& r_child : children)
                        splittable = splittable &&
                                     (!considered(other) || volume_and_tolerance(r_child).first > 0.0);
                if (!splittable) {
                    rTets.mPoints.pop_back();
                    continue;
                }
                const auto on_skin = [&](const Index I) {
                    return rTets.mOnSkin[static_cast<std::size_t>(I)] != 0;
                };
                rTets.mOnSkin.push_back((on_skin(tet[i]) && on_skin(tet[j])) ||
                                        (on_skin(tet[k]) && on_skin(tet[l])));
                rTets.mFeature.push_back(0);
                touched.push_back(1);
                dead[t] = 1;
                for (const auto& [other, children] : splits) {
                    dead[other] = 1;
                    for (const Index n : rTets.mTets[other])
                        touched[static_cast<std::size_t>(n)] = 1;
                    for (const auto& r_child : children) {
                        rTets.mTets.push_back(r_child);
                        rTets.mMaterial.push_back(rTets.mMaterial[other]);
                        rTets.mCut.push_back(1);
                        rTets.mLevel.push_back(rTets.mLevel[other]);
                        rTets.mShape.push_back(0);
                    }
                }
                dead.resize(rTets.mTets.size(), 0);
                for (const Index n : tet)
                    touched[static_cast<std::size_t>(n)] = 1;
                ++rR.mNumSlivers;
                ++num_repaired;
                break;
            }
        }
        if (num_repaired == 0)
            break;
        compact(dead);
    }

    // An inverted tetrahedron is untangled by moving a free node (neither on the skin nor on a
    // feature) to the centroid of the nodes around it, when this leaves its whole star positive.
    build_star();
    for (std::size_t t = 0; t < rTets.mTets.size(); ++t) {
        if (!considered(t) || volume_and_tolerance(rTets.mTets[t]).first >= 0.0)
            continue;
        for (const Index node : rTets.mTets[t]) {
            if (rank(node) != 0)
                continue;
            const std::size_t n = static_cast<std::size_t>(node);
            Vec3 centroid{0.0, 0.0, 0.0};
            double count = 0.0;
            for (std::size_t k = start[n]; k < start[n + 1]; ++k)
                for (const Index v : rTets.mTets[star[k]])
                    if (v != node) {
                        centroid = vec3_add(centroid, point(v));
                        count += 1.0;
                    }
            const Vec3 previous = rTets.mPoints[n];
            rTets.mPoints[n] = vec3_scale(centroid, 1.0 / count);
            bool positive = true;
            for (std::size_t k = start[n]; k < start[n + 1] && positive; ++k)
                if (considered(star[k])) {
                    const auto [volume, tolerance] = volume_and_tolerance(rTets.mTets[star[k]]);
                    positive = volume > tolerance;
                }
            if (positive) {
                ++rR.mNumUntangled;
                break;
            }
            rTets.mPoints[n] = previous;
        }
    }

    // The remaining negative tetrahedra (including those the collapse left within the
    // tolerance) get two nodes swapped.
    for (std::size_t t = 0; t < rTets.mTets.size(); ++t) {
        if (!considered(t))
            continue;
        if (volume_and_tolerance(rTets.mTets[t]).first < 0.0) {
            std::swap(rTets.mTets[t][1], rTets.mTets[t][2]);
            rTets.mCut[t] = 1;
            ++rR.mNumReordered;
        }
    }

    // Two tetrahedra that share a face overlap when their opposite nodes lie on the same side
    // of it.
    std::vector<std::pair<std::array<Index, 3>, std::pair<std::size_t, Index>>> faces;
    for (std::size_t t = 0; t < rTets.mTets.size(); ++t) {
        if (!considered(t))
            continue;
        for (const auto& facet : kOctmtFacet) {
            std::array<Index, 3> key{rTets.mTets[t][facet[0]], rTets.mTets[t][facet[1]],
                                     rTets.mTets[t][facet[2]]};
            std::sort(key.begin(), key.end());
            // The node of the tetrahedron which is not on the face.
            Index opposite = -1;
            for (const Index v : rTets.mTets[t])
                if (v != key[0] && v != key[1] && v != key[2])
                    opposite = v;
            faces.push_back({key, {t, opposite}});
        }
    }
    std::sort(faces.begin(), faces.end());
    for (std::size_t i = 0; i + 1 < faces.size(); ++i) {
        if (faces[i].first != faces[i + 1].first ||
            (i + 2 < faces.size() && faces[i + 2].first == faces[i].first))
            continue;
        const Vec3& a = point(faces[i].first[0]);
        const Vec3 normal = vec3_cross(vec3_sub(point(faces[i].first[1]), a),
                                       vec3_sub(point(faces[i].first[2]), a));
        const double side_0 = vec3_dot(normal, vec3_sub(point(faces[i].second.second), a));
        const double side_1 = vec3_dot(normal, vec3_sub(point(faces[i + 1].second.second), a));
        rR.mNumOverlaps += side_0 * side_1 > 0.0;
    }

    if (rR.mNumCollapsed > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(rR.mNumCollapsed) +
                                 " degenerate tetrahedra collapsed");
    if (rR.mNumCapsRemoved > 0 || rR.mNumSplit > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(rR.mNumCapsRemoved) +
                                 " flat caps dropped and " + std::to_string(rR.mNumSplit) +
                                 " tetrahedra split to repair flat caps");
    if (rR.mNumSlivers > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(rR.mNumSlivers) +
                                 " flat slivers removed by inserting a node at the crossing of "
                                 "their opposite edges");
    if (rR.mNumUntangled > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(rR.mNumUntangled) +
                                 " inverted tetrahedra untangled by moving a free node");
    if (rR.mNumReordered > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(rR.mNumReordered) +
                                 " inverted tetrahedra reordered");
    if (rR.mNumOverlaps > 0)
        octree_mesh_warn(rR, std::string(kOctreeMeshPrefix) + std::to_string(rR.mNumOverlaps) +
                                 " overlapping tetrahedra pairs");
}

double octree_mesh_hausdorff_to_skin(const OctreeMeshSetup& rS, const std::vector<Vec3>& rPoints,
                                     const std::vector<std::array<std::int64_t, 3>>& rFaces) {
    const Mesh out_skin = mesh_from_soup(rPoints, rFaces, "triangle");
    const Mesh in_skin = mesh_from_soup(rS.mSkin.mPoints, rS.mSkin.mTris, "triangle");
    HausdorffOptions ho;
    ho.mFaceSamples = 2;
    return hausdorff_distance(out_skin, in_skin, ho).mDistance;
}

void octree_mesh_emit_tets(OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR,
                           const OctreeMeshTets& rIn) {
    // Keep the inside materials (and the outside on request), compact the points.
    std::vector<std::size_t> kept;
    for (std::size_t t = 0; t < rIn.mTets.size(); ++t)
        if (rIn.mMaterial[t] > 0 || rO.mKeepOutside)
            kept.push_back(t);
    std::vector<std::int64_t> remap(rIn.mPoints.size(), -1);
    for (std::size_t t : kept)
        for (std::int64_t v : rIn.mTets[t])
            remap[static_cast<std::size_t>(v)] = 0;
    std::int64_t np = 0;
    for (std::int64_t& r : remap)
        if (r == 0)
            r = np++;
    NDArray points =
        NDArray::Uninit(DType::Float64, {static_cast<std::size_t>(np), std::size_t{3}});
    std::vector<Vec3> pts(static_cast<std::size_t>(np));
    for (std::size_t p = 0; p < rIn.mPoints.size(); ++p)
        if (remap[p] >= 0) {
            pts[static_cast<std::size_t>(remap[p])] = rIn.mPoints[p];
            for (std::size_t k = 0; k < 3; ++k)
                points.As<double>()[static_cast<std::size_t>(remap[p]) * 3 + k] = rIn.mPoints[p][k];
        }
    NDArray conn = NDArray::Uninit(DType::Int64, {kept.size(), std::size_t{4}});
    NDArray level = NDArray::Uninit(DType::Int64, {kept.size()});
    NDArray shape = NDArray::Uninit(DType::Int64, {kept.size()});
    NDArray material = NDArray::Uninit(DType::Int64, {kept.size()});
    std::vector<std::vector<std::int64_t>> by_material(static_cast<std::size_t>(rIn.mNumMaterials));
    for (std::size_t n = 0; n < kept.size(); ++n) {
        const std::size_t t = kept[n];
        for (std::size_t c = 0; c < 4; ++c)
            conn.As<std::int64_t>()[n * 4 + c] = remap[static_cast<std::size_t>(rIn.mTets[t][c])];
        level.As<std::int64_t>()[n] = rIn.mLevel[t];
        shape.As<std::int64_t>()[n] = rIn.mShape[t];
        material.As<std::int64_t>()[n] = rIn.mMaterial[t];
        by_material[static_cast<std::size_t>(rIn.mMaterial[t])].push_back(
            static_cast<std::int64_t>(n));
    }
    Mesh& m = rR.mMesh;
    m.AssignPoints(std::move(points));
    m.AddCellBlock("tetra", std::move(conn));
    {
        std::vector<NDArray> lv, sh, mt;
        lv.push_back(std::move(level));
        sh.push_back(std::move(shape));
        mt.push_back(std::move(material));
        m.AddCellData(kOctreeLevelName, std::move(lv));
        m.AddCellData(kOctreeMeshShapeName, std::move(sh));
        if (rIn.mNumMaterials > 2 || rO.mKeepOutside)
            m.AddCellData(kOctreeMeshMaterialName, std::move(mt));
    }
    // Regions: inside (every material above 0), outside, each material, and
    // the boundary faces.
    std::vector<std::int64_t> inside;
    for (std::size_t mat = 1; mat < by_material.size(); ++mat)
        inside.insert(inside.end(), by_material[mat].begin(), by_material[mat].end());
    std::sort(inside.begin(), inside.end());
    rR.mNumInside = static_cast<std::int64_t>(inside.size());
    m.AddRegion(Region(kOctreeMeshInsideRegion, RegionKind::Cell, 3, -1, octmt_ids(inside, 1)));
    if (rO.mKeepOutside)
        m.AddRegion(Region(kOctreeMeshOutsideRegion, RegionKind::Cell, 3, -1,
                           octmt_ids(by_material[0], 1)));
    if (rIn.mNumMaterials > 2)
        for (std::size_t mat = 0; mat < by_material.size(); ++mat)
            if (mat > 0 || rO.mKeepOutside)
                m.AddRegion(Region("material_" + std::to_string(mat), RegionKind::Cell, 3,
                                   static_cast<std::int64_t>(mat), octmt_ids(by_material[mat], 1)));
    // Boundary faces: used once.
    std::vector<std::pair<std::array<std::int64_t, 3>, std::int64_t>> faces;
    const std::int64_t* cp = m.Cells(0).Conn().As<std::int64_t>();
    faces.reserve(kept.size() * 4);
    for (std::size_t n = 0; n < kept.size(); ++n)
        for (std::int64_t f = 0; f < 4; ++f) {
            std::array<std::int64_t, 3> key{cp[n * 4 + kOctmtFacet[f][0]],
                                            cp[n * 4 + kOctmtFacet[f][1]],
                                            cp[n * 4 + kOctmtFacet[f][2]]};
            std::sort(key.begin(), key.end());
            faces.emplace_back(key, static_cast<std::int64_t>(n) * 4 + f);
        }
    parallel_sort(faces.begin(), faces.end());
    std::vector<std::int64_t> sides;
    for (std::size_t i = 0; i < faces.size();) {
        std::size_t j = i;
        while (j < faces.size() && faces[j].first == faces[i].first)
            ++j;
        if (j - i == 1) {
            sides.push_back(faces[i].second / 4);
            sides.push_back(faces[i].second % 4);
        }
        i = j;
    }
    {
        // Ascending (cell, facet).
        std::vector<std::pair<std::int64_t, std::int64_t>> pairs;
        for (std::size_t s = 0; s < sides.size(); s += 2)
            pairs.emplace_back(sides[s], sides[s + 1]);
        std::sort(pairs.begin(), pairs.end());
        sides.clear();
        for (const auto& [c, f] : pairs) {
            sides.push_back(c);
            sides.push_back(f);
        }
    }
    m.AddRegion(Region(kOctreeMeshBoundaryRegion, RegionKind::Side, 2, -1, octmt_ids(sides, 2)));
    // The skin of the inside: its faces that no other inside tetrahedron
    // shares. Its edges must each lie on two of its faces, and its distance
    // to the input skin is the report's Hausdorff distance.
    {
        std::vector<std::array<std::int64_t, 3>> skin_faces;
        for (std::size_t i = 0; i < faces.size();) {
            std::size_t j = i;
            std::int64_t inside_uses = 0, last = -1;
            while (j < faces.size() && faces[j].first == faces[i].first) {
                const auto n = static_cast<std::size_t>(faces[j].second / 4);
                if (rIn.mMaterial[kept[n]] > 0) {
                    ++inside_uses;
                    last = faces[j].second;
                }
                ++j;
            }
            if (inside_uses == 1) {
                const auto n = static_cast<std::size_t>(last / 4);
                const auto f = static_cast<std::size_t>(last % 4);
                skin_faces.push_back({cp[n * 4 + kOctmtFacet[f][0]], cp[n * 4 + kOctmtFacet[f][1]],
                                      cp[n * 4 + kOctmtFacet[f][2]]});
            }
            i = j;
        }
        std::vector<std::array<std::int64_t, 2>> edges;
        for (const auto& f : skin_faces)
            for (int e = 0; e < 3; ++e)
                edges.push_back({std::min(f[e], f[(e + 1) % 3]), std::max(f[e], f[(e + 1) % 3])});
        std::sort(edges.begin(), edges.end());
        for (std::size_t i = 0; i < edges.size();) {
            std::size_t j = i;
            while (j < edges.size() && edges[j] == edges[i])
                ++j;
            if (j - i == 1)
                ++rR.mNumBoundaryEdgesOpen;
            else if (j - i > 2)
                ++rR.mNumNonManifoldBoundaryEdges;
            i = j;
        }
        if (!skin_faces.empty())
            rR.mHausdorffToSkin = octree_mesh_hausdorff_to_skin(rS, pts, skin_faces);
    }
    if (std::any_of(rIn.mFeature.begin(), rIn.mFeature.end(),
                    [](std::uint8_t f) { return f != 0; })) {
        std::vector<std::int64_t> corners, chains;
        for (std::size_t p = 0; p < rIn.mPoints.size(); ++p)
            if (remap[p] >= 0 && rIn.mFeature[p])
                (rIn.mFeature[p] == 1 ? corners : chains).push_back(remap[p]);
        std::sort(corners.begin(), corners.end());
        std::sort(chains.begin(), chains.end());
        m.AddRegion(
            Region(kOctreeMeshCornersRegion, RegionKind::Point, 0, -1, octmt_ids(corners, 1)));
        chains.insert(chains.end(), corners.begin(), corners.end());
        std::sort(chains.begin(), chains.end());
        m.AddRegion(
            Region(kOctreeMeshFeatureEdgesRegion, RegionKind::Point, 1, -1, octmt_ids(chains, 1)));
    }
    if (rO.mAttachDistance) {
        std::vector<double> dist;
        std::vector<Vec3> closest;
        octree_mesh_signed_distance(rS, pts, rO, rR, dist, closest);
        NDArray d = NDArray::Uninit(DType::Float64, {pts.size()});
        std::copy(dist.begin(), dist.end(), d.As<double>());
        m.AddPointData(kOctreeMeshDistanceName, std::move(d));
    }
}

}  // namespace detail
}  // namespace meshioplusplus
