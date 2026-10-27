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
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "../detail/bcc_octree_tets.hpp"
#include "../detail/feature_snap.hpp"
#include "../detail/isosurface_stuffing.hpp"
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

double octree_mesh_hausdorff_to_skin(const OctreeMeshSetup& rS, const std::vector<Vec3>& rPoints,
                                     const std::vector<std::array<std::int64_t, 3>>& rFaces) {
    Mesh out_skin;
    NDArray sp = NDArray::Uninit(DType::Float64, {rPoints.size(), std::size_t{3}});
    for (std::size_t q = 0; q < rPoints.size(); ++q)
        for (std::size_t k = 0; k < 3; ++k)
            sp.As<double>()[q * 3 + k] = rPoints[q][k];
    NDArray sc = NDArray::Uninit(DType::Int64, {rFaces.size(), std::size_t{3}});
    for (std::size_t q = 0; q < rFaces.size(); ++q)
        for (std::size_t k = 0; k < 3; ++k)
            sc.As<std::int64_t>()[q * 3 + k] = rFaces[q][k];
    out_skin.AssignPoints(std::move(sp));
    out_skin.AddCellBlock("triangle", std::move(sc));
    Mesh in_skin;
    NDArray ip = NDArray::Uninit(DType::Float64, {rS.mSkin.mPoints.size(), std::size_t{3}});
    for (std::size_t q = 0; q < rS.mSkin.mPoints.size(); ++q)
        for (std::size_t k = 0; k < 3; ++k)
            ip.As<double>()[q * 3 + k] = rS.mSkin.mPoints[q][k];
    NDArray ic = NDArray::Uninit(DType::Int64, {rS.mSkin.mTris.size(), std::size_t{3}});
    for (std::size_t q = 0; q < rS.mSkin.mTris.size(); ++q)
        for (std::size_t k = 0; k < 3; ++k)
            ic.As<std::int64_t>()[q * 3 + k] = rS.mSkin.mTris[q][k];
    in_skin.AssignPoints(std::move(ip));
    in_skin.AddCellBlock("triangle", std::move(ic));
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
